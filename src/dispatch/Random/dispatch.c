#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <emmintrin.h>
#include "pfOutput.h"
#include "pfCompute.h"
#include "rd_threads.h"

static const char DNA_LUT[4] = { 'A', 'C', 'G', 'T' };

static int thpool_add_work(threadpool_t * const restrict thpool_p, const RandomData_t * const restrict RD)
{
	if (!RD) return -1;
	
	struct timespec polling_interval;
	polling_interval.tv_sec  = 0;
	polling_interval.tv_nsec = 1000000;

	/* Allocating memory */
	rd_job_t * jobBatch1 = (rd_job_t*) malloc(2*JOBS_BATCH_SIZE*sizeof(rd_job_t));
	if (jobBatch1 == NULL) return -2;	
	rd_job_t * jobBatch2 = &jobBatch1[JOBS_BATCH_SIZE];
	const size_t SeqLength = RD->Length;
	const size_t TotalLength = RD->Length + 3*RD->nCAGSuffix;
	const size_t RD32bitLength = 1 + (SeqLength/32);
	char * const ReadsMemory = (char*) malloc( 2*JOBS_BATCH_SIZE*TotalLength*sizeof(char));
	if (ReadsMemory == NULL) return -3;
	
	{
		char * ptr = (char*) &ReadsMemory[0];
		for (size_t m=0;m<2*JOBS_BATCH_SIZE;m++) {
			jobBatch1[m].Reads = ptr; 
			jobBatch1[m].Length = TotalLength;
			ptr += TotalLength;
		}
	}
	size_t iJob = 0;
	
	srandom(RD->seed);

	const size_t N = RD->N;
	
	while (iJob < N) {
		size_t i = iJob;
		iJob += JOBS_BATCH_SIZE;
		const size_t lN = iJob < N ? iJob : N;
		int k = 0;
		while (i<lN) {
			
			/* Do LUT of that to generate the DNA sequence */
			{
				char * ptr = jobBatch1[k].Reads;
				uint32_t buffer = random();
				size_t l;
				for (l=0; l<SeqLength; l++) {
					*ptr++ = DNA_LUT[buffer & 0x3];
					buffer >>= 2;
					if ((l & (size_t) 7) == (size_t) 7) buffer = random();
				}
				while (l<TotalLength) {
					ptr[0] = 'C';
					ptr[1] = 'A';
					ptr[2] = 'G';
					ptr += 3;
					l += 3;
				}
			}
// 			jobqueue_push(thpool_p, &jobBatch1[k]);
			++i;
			++k;
		}
		
 		jobqueue_push(thpool_p, jobBatch1, k);
	
		/* Wait for jobs to be done before reloading */
		if (iJob >= 2*JOBS_BATCH_SIZE && iJob < N) {
			while ( (volatile int) (thpool_p->jobqueue_p->len) >= JOBS_BATCH_SIZE) {
				nanosleep(&polling_interval, NULL);
			}
			/*break*/;
		}
		
		/* swap pointers */
		register rd_job_t * tmp = jobBatch1;
		jobBatch1 = jobBatch2;
		jobBatch2 = tmp;
	}
	
	return SUCCESS;
}


int dispatchRandomData(const struct Profile * const restrict prf,
											 const RandomData_t * const restrict RD,
#ifdef PRF_CORE_PCRE
											 struct RegEx * const restrict regex,
#endif
											 const size_t nCPUs)
{
	int res = SUCCESS;
	
	/*************************************************************************/
  /*                    ALLOCATE HISTOGRAM MEMORY                          */
  /*************************************************************************/
	const size_t PaddedSize = (HistogramSize + 3UL) & ~(3UL);
	
  unsigned int * const restrict histograms = (unsigned int*) _mm_malloc(nCPUs*PaddedSize*sizeof(unsigned int),16); 
	if( histograms == NULL ) {
		fputs("Unable to allocate memory for the histograms\n", stderr);
		return -1;
	}
	memset(histograms, 0, nCPUs*PaddedSize*sizeof(unsigned int));
	
	/*************************************************************************/
  /*                    SET THE THREADS COMMON VARIABLES                   */
  /*************************************************************************/
	common_t common = { 	.profile = prf, .RD = RD
#ifdef PRF_CORE_PCRE
								, .regex = regex
#endif
	};
	
  /*************************************************************************/
  /*                         CHOOSE THREAD FUNCTION                        */
  /*************************************************************************/
	void* (*ThreadFct)(void*);
#ifdef PRF_CORE_PCRE
	if (regex != NULL && prf != NULL) {
		fputs("Combined histogram not implemented yet!\n", stderr);
		return -1;
	}
	else if (regex != NULL && prf == NULL) {
		ThreadFct = threadpool_random_regex_histogram;
	}
	else
#endif
	{
		if (prf) {
			ThreadFct = threadpool_random_repeat_histogram;
		}
		else {
			fprintf(stderr,"Please provide a profile and/or a regex string!\n");
			return -1;
		}
	}
	
 /*************************************************************************/
  /*                          PREPARE THREAD POOL                          */
  /*************************************************************************/  
  
  threads_on_hold   = 0;
  threads_keepalive = 1;
  
  // Make new thread pool
  threadpool_t thpool;
  
  pthread_mutex_init(&(thpool.thcount_lock), NULL);
  thpool.num_threads_alive   = 0;
  thpool.num_threads_working = 0;
  
  // Allocate memory for threads
  ThreadPool_Arg_t * const restrict threads_arg = (ThreadPool_Arg_t*) malloc(nCPUs*sizeof(ThreadPool_Arg_t));
  if (threads_arg == NULL) return -1;
  thpool.threads = (ThreadPool_Arg_t *) &threads_arg;
  
  // Initialize job queue
  if (jobqueue_init(&thpool) == -1 ) {
		res = -2;
		goto FIN2;
  }
      
  /*************************************************************************/
  /*                         START THE THREADS POOL                        */
  /*************************************************************************/

  // Initialize and start threads
  for (size_t i=0; i<nCPUs; i++) {
      threads_arg[i].ID     = i;
			threads_arg[i].common = &common;
			threads_arg[i].thpool = &thpool;
			threads_arg[i].histogram = histograms + i*PaddedSize;
			threads_arg[i].returnState = SUCCESS;

#ifdef PRF_USE_AFFINITY
			if (pthread_create (&(threads_arg[i].pthread), &threads_attr[i], ThreadFct,  (void*) &threads_arg[i]) != 0)
#else
			if (pthread_create (&(threads_arg[i].pthread),  NULL, ThreadFct, (void*) &threads_arg[i]) != 0)
#endif
			{
				fprintf(stderr, " Unable to start some threads\n"); 
			  int k = i;
			  while (--k>=0) pthread_kill(threads_arg[k].pthread, -9);
        res = -4;
        goto FIN0;
      }
      else {
			pthread_detach(threads_arg[i].pthread);
	  }
  }
  
  
  /*************************************************************************/
  /*                         WAIT FOR POOL READY                           */
  /*************************************************************************/
  while (thpool.num_threads_alive != nCPUs) {}
  
  /*************************************************************************/
  /*                        ADD TASKS TO JOB QUEUE                         */
  /*************************************************************************/
 
	struct timeval _t0, _t1;

	gettimeofday(&_t0,0);
	res = thpool_add_work(&thpool, RD);

  if (res != 0 ) {
		fprintf(stderr, "Master failed to send jobs, error code %i\n", res);
		goto KILL_THREADS;
	}
	
 /*************************************************************************/
  /*                       WAIT FOR TASK TO FINISH                         */
  /*************************************************************************/
  thpool_wait(&thpool);
	gettimeofday(&_t1,0);
	
	if (OutputVerbose) {
		const double t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
		fprintf(stderr,"This took %lf seconds to crunch on %lu cores.\n", t, nCPUs);

	}
	
	/*************************************************************************/
  /*                          FREE SOME MEMORY                             */
  /*************************************************************************/
  
  /*************************************************************************/
  /*                  TELL THREADS TO FLUSH AND TERMINATE                  */
  /*************************************************************************/
  thpool_destroy(&thpool);
  pthread_mutex_destroy(&(thpool.thcount_lock));
  
  /*************************************************************************/
  /*                         PROCESS DATA                                  */
  /*************************************************************************/ 
  for (size_t j=1; j<nCPUs; j++) {
		const unsigned int * const restrict ptr = &histograms[j*PaddedSize]; 
		for (size_t i=0; i<PaddedSize; i+=4) {
			__m128i _A = _mm_load_si128((__m128i*) &histograms[i]);
			_A = _mm_add_epi32(_A, *((__m128i*) &ptr[i]));
			_mm_store_si128((__m128i*)&histograms[i], _A);
		}
	}
	
  /*************************************************************************/
  /*             WAIT FOR THREADS TO TERMINATE BEFORE CLOSING              */
  /*************************************************************************/
#ifdef _AIO_
#define MAX_NANOSEC 999999999
#define CEIL(X) ((X-(int)(X)) > 0 ? (int)(X+1) : (int)(X))
  /* Exponential polling */
  {
    long init_nano =  1; /* MUST be above 0 */
    long new_nano;
    double multiplier = 1.01;
    int  max_secs = 20;
    
    struct timespec polling_interval;
    polling_interval.tv_sec  = 0;
    polling_interval.tv_nsec = init_nano;
    
    while (thpool.num_threads_alive > 0) {
			nanosleep(&polling_interval, NULL);
			if ( polling_interval.tv_sec < max_secs ){
			  new_nano = CEIL(polling_interval.tv_nsec * multiplier);
			  polling_interval.tv_nsec = new_nano % MAX_NANOSEC;
			  if ( new_nano > MAX_NANOSEC ) {
			    polling_interval.tv_sec ++;
			  }
			}
			else break;
    }
  }
#undef MAX_NANOSEC
#undef CEIL
#endif
  /*************************************************************************/
  /*                    CHECK FOR THREADS ERROR                            */
  /*************************************************************************/
  for (size_t i=0; i<nCPUs; i++) {
    if (threads_arg[i].returnState != SUCCESS) {
			res = threads_arg[i].returnState;
			goto FIN0;
    }
  }
  
  char FName[256];
	snprintf(FName, 256, "%s.histogram", OutputType->Specific.Histogram.BaseFileName);
	FILE * out = fopen(FName, "w");
	if (out != NULL) {
		for (unsigned int j=0; j<HistogramSize; j++) {
			fprintf(out, "%u\t%u\n", j, histograms[j]);
		}
		fclose(out);
	}
	else {
		fprintf(stderr, "Unable to create output histogram %s\n", FName);
	}
  
  /*************************************************************************/
  /*                         GATHER EXTRA DATA                             */
  /*************************************************************************/ 
	size_t MissedCycles = threads_arg[1].counter;
	for (size_t i=1; i<nCPUs; i++) {
		MissedCycles += threads_arg[i].counter;
	}
	if (MissedCycles) {
		fprintf(stderr, "Some sequences (%lu) bear alignment that are longer than the histogram bin number\n", MissedCycles);
	}
     
  /*************************************************************************/
  /*                      DESTROY THE THREAD POOL                          */
  /*************************************************************************/
  goto FIN0;
  
  KILL_THREADS:
  thpool_destroy(&thpool);
  pthread_mutex_destroy(&(thpool.thcount_lock));
  
  FIN0:
  free(threads_arg);
  FIN2:
  _mm_free(histograms);
	return res;
}

