#define _GNU_SOURCE
#include "prf_config.h"
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <sys/time.h>
#include <pthread.h>
#include "pfInput.h"
#include "pfCompute.h"
#include "pfOutput.h"
#include "FASTA_threads.h"
#ifdef PRF_CORE_PCRE
_Bool RegexCycleCount = false;      /* triggers to increase regex size up to histogram limit*/
#endif
static const struct timespec polling_interval = { .tv_sec = 0, .tv_nsec = 400};

int dispatchFASTAFile(const struct Profile * restrict prf,
                      const FASTAStructure * const restrict FASTA,
#ifdef PRF_CORE_PCRE
                      struct RegEx * const restrict regex,
#endif
                      const OutputType_t * const restrict OutputType,
                      const size_t nCPUs)
{
	struct timeval _t0, _t1;
	if (!prf) {
#ifdef PRF_CORE_PCRE
		if (regex) {
			if (OutputType->Type != HISTOGRAM) {
				fprintf(stderr,"Regex on FASTA file only for histogram yet!\n");
				return -1;
			}
		}
		else 
#endif
		{
			fprintf(stderr,"dispatchFASTAFile needs at least a profile and/or a regex!\n");
			return -1;
		}
	}
	
  /* Prepare structure common to filter and alignment */
  size_t * const restrict shares = alloca((nCPUs+1)*sizeof(size_t));

  /* Allocate stack memory for posix thread structures */
  pthread_t * const restrict threads = (pthread_t*) alloca(nCPUs*sizeof(pthread_t));
	
	
  /* Allocate stack memory for posix thread structures */
  struct ThreadData * const threads_arg = alloca(nCPUs*sizeof(struct ThreadData));
	memset(threads_arg, 0, nCPUs*sizeof(struct ThreadData));
  
  /* Dispatch work to threads share according to file size */
  size_t FileShare = (size_t) FASTA->FileSize / nCPUs;
  FileShare += ((size_t) FASTA->FileSize % nCPUs) > (nCPUs-1) ? 1 : 0;
  const s_Data * DataPtr = FASTA->DataPtr;
  register size_t counter = 0;
  shares[0] = 0;
  for (size_t i=1; i<nCPUs; ++i) {
    register size_t tmp = i*FileShare;
    while ( (size_t) DataPtr->Offset < tmp) { ++DataPtr; ++counter; }
    shares[i] = counter;
  }
  shares[nCPUs] = FASTA->SequenceCount;

  pthread_mutex_t PrintLock = PTHREAD_MUTEX_INITIALIZER;

#if (defined(PRF_CORE_STD) || defined(PRF_CORE_REPEAT))
#if (defined(PRF_CORE_STD) && defined(PRF_CORE_REPEAT))
	const Compute_t * const Model = (prf->isCircular) ? &Repeat_sse41 : &Standard_sse41;
#elif defined(PRF_CORE_STD)
	const Compute_t * const Model = &Standard_sse41;
#else
	const Compute_t * const Model = &Repeat_sse41;
#endif
#else
#error "Dispatching procedures requires at lest STD or REPEAT methods"
#endif
	
  /* Dispatch common information to threads */
  for (size_t i=0; i<nCPUs; ++i) {
		threads_arg[i].prf          = prf;
		threads_arg[i].FASTA        = FASTA;
		threads_arg[i].SequenceID   = NULL;
		threads_arg[i].threadId     = i;
		threads_arg[i].PrintLock    = &PrintLock;
		threads_arg[i].Compute      = Model;
		threads_arg[i].OutputType   = OutputType;
		threads_arg[i].start        = shares[i];          
		threads_arg[i].stop         = shares[i+1];
#ifdef PRF_CORE_PCRE
		threads_arg[i].regex        = regex;
#endif
  }
  
  size_t * restrict Histograms;
	size_t HistogramSize = 0UL;
	size_t DensitySize;

	/* Get the thread function */
	void * (*ThreadFct)(void*) = NULL;
	{
		const int index = GetDispatchThreadIndex(OutputType, prf);
		if (index >= 0) {
		ThreadFct = t_fap[index];
		}
		else {
			fputs("Error in the choice of thread function\n", stderr);
			return -1;
		}
	}
	
	if (OutputType->Type == HISTOGRAM) {
		if (OutputType->Specific.Histogram.CycleRatherThanScore) {
			HistogramSize = OutputType->CycleRange[1] - OutputType->CycleRange[0] + 1; 
		}
		else {
			HistogramSize = OutputType->ScoreRange[1] - OutputType->ScoreRange[0] + 1; 
		}
		Histograms = (size_t *) calloc(HistogramSize*nCPUs,sizeof(size_t));
		if (Histograms == NULL) {
				fprintf(stderr, "Unable to allocate memory for histograms, requested size was %lu bytes.\n",  HistogramSize*nCPUs*sizeof(size_t));
				return -1;
		}
		size_t * restrict ptr = Histograms;
		for (size_t i=0; i<nCPUs; ++i) {
			threads_arg[i].Histogram = ptr;
			ptr += HistogramSize;
		}
	}
	else if(OutputType->Type == DENSITY) {
		DensitySize = (OutputType->CycleRange[1] - OutputType->CycleRange[0] + 1)\
		             *(OutputType->ScoreRange[1] - OutputType->ScoreRange[0] + 1);
		Histograms = (size_t *) calloc(DensitySize*nCPUs, sizeof(size_t));
		if (Histograms == NULL) {
				fprintf(stderr, "Unable to allocate memory for histograms, requested size was %lu bytes.\n",  HistogramSize*nCPUs*sizeof(size_t));
				return -1;
		}
		size_t * restrict ptr = Histograms;
		for (size_t i=0; i<nCPUs; ++i) {
			threads_arg[i].Histogram = ptr;
			ptr += DensitySize;
		}
	}
  
  /* Dispatch to threads */
  gettimeofday(&_t0,0);
#ifdef PRF_USE_AFFINITY
	if (threads_attr) {
		for (size_t i=0; i<nCPUs; ++i) {
			if (pthread_create (&threads[i],  &threads_attr[i], ThreadFct,  (void*) &threads_arg[i]) != 0) {
				return -2;
			}
		}
	}
	else
#endif
	{
		for (size_t i=0; i<nCPUs; ++i) {
			if (pthread_create (&threads[i],  NULL, ThreadFct,  (void*) &threads_arg[i]) != 0) {
				return -2;
			}
		}
	}

	for (size_t i=0; i<nCPUs; i++) {
		pthread_join(threads[i], NULL);
	}
	gettimeofday(&_t1,0);

	pthread_mutex_destroy(&PrintLock);

	if (OutputType->Type == HISTOGRAM) {
		/* Gather thread histograms */
		size_t MissedCycles = threads_arg[0].counter;
		const size_t * restrict hptr = &Histograms[HistogramSize];
		for (size_t i=1; i<nCPUs; i++) {
			MissedCycles += threads_arg[i].counter;
			for (unsigned int j=0; j<HistogramSize; j++) Histograms[j] += hptr[j];
			hptr += HistogramSize;
		}
		if (MissedCycles) {
			fprintf(stderr, "Some sequences (%lu) bear alignment that are longer than the histogram bin number\n", MissedCycles);
		}
		
		char FName[256];
		snprintf(FName, 256, "%s.histogram", OutputType->Specific.Histogram.BaseFileName);
		FILE * out = fopen(FName, "w");
		if (out != NULL) {
			const int RangeStart = (OutputType->Specific.Histogram.CycleRatherThanScore) ? OutputType->CycleRange[0] :\
			                       OutputType->ScoreRange[0];
			for (unsigned int j=0; j<HistogramSize; j++) {
				fprintf(out, "%i\t%lu\n", RangeStart+j, Histograms[j]);
			}
			fclose(out);
		}
		else {
			fprintf(stderr, "Unable to create output histogram %s\n", FName);
			return -3;
		}
	}
	else if (OutputType->Type == DENSITY) {
		/* Gather thread histograms */
		size_t MissedCycles = threads_arg[0].counter;
		const size_t * restrict hptr = &Histograms[DensitySize];
		for (size_t i=1; i<nCPUs; i++) {
			MissedCycles += threads_arg[i].counter;
			for (unsigned int j=0; j<DensitySize; j++) Histograms[j] += hptr[j];
			hptr += DensitySize;
		}
		if (MissedCycles) {
			fprintf(stderr, "Some sequences (%lu) bear alignment that are longer than the histogram bin number\n", MissedCycles);
		}
		
		char FName[256];
		snprintf(FName, 256, "%s.density", OutputType->Specific.Histogram.BaseFileName);
		FILE * out = fopen(FName, "w");
		if (out != NULL) {
			hptr = Histograms;
			for (int i=OutputType->ScoreRange[0]; i<OutputType->ScoreRange[1]; i++) {
				for (unsigned int j=OutputType->CycleRange[0]; j<OutputType->CycleRange[1]; j++) {
					fprintf(out, "%i\t%u\t%lu\n", i, j,* hptr++);
				}
				fprintf(out, "\n");
			}
			fclose(out);
		}
		else {
			fprintf(stderr, "Unable to create output histogram %s\n", FName);
			return -3;
		}
	}
	else {
		if (OutputVerbose) {
			unsigned int AlignedSequencesCounter = threads_arg[0].counter;
			for (size_t i=1; i<nCPUs; i++) AlignedSequencesCounter += threads_arg[i].counter;
			const double t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr,"Overall there are %u aligned sequences found. These took %lf seconds to align on %li cores.\n", AlignedSequencesCounter, t, nCPUs);
		}
	}
	return SUCCESS;
}

int dispatchFASTAFileExt(const struct Profile * restrict prf,
                         const FASTAStructure * const restrict FASTA,
												 const Compute_t * const Model,
#ifdef PRF_CORE_PCRE
                         struct RegEx * const restrict regex,
#endif
                         const OutputType_t * const restrict OutputType,
                         const size_t nCPUs)
{
	struct timeval _t0, _t1;
	if (!prf) {
#ifdef PRF_CORE_PCRE
		if (regex) {
			if (OutputType->Type != HISTOGRAM) {
				fprintf(stderr,"Regex on FASTA file only for histogram yet!\n");
				return -1;
			}
		}
		else 
#endif
		{
			fprintf(stderr,"dispatchFASTAFile needs at least a profile and/or a regex!\n");
			return -1;
		}
	}
	
  /* Prepare structure common to filter and alignment */
  size_t * const restrict shares = alloca((nCPUs+1)*sizeof(size_t));

  /* Allocate stack memory for posix thread structures */
  pthread_t * const restrict threads = (pthread_t*) alloca(nCPUs*sizeof(pthread_t));
	
	
  /* Allocate stack memory for posix thread structures */
  struct ThreadData * const threads_arg = alloca(nCPUs*sizeof(struct ThreadData));
	memset(threads_arg, 0, nCPUs*sizeof(struct ThreadData));
  
  /* Dispatch work to threads share according to file size */
  size_t FileShare = (size_t) FASTA->FileSize / nCPUs;
  FileShare += ((size_t) FASTA->FileSize % nCPUs) > (nCPUs-1) ? 1 : 0;
  const s_Data * DataPtr = FASTA->DataPtr;
  register size_t counter = 0;
  shares[0] = 0;
  for (size_t i=1; i<nCPUs; ++i) {
    register size_t tmp = i*FileShare;
    while ( (size_t) DataPtr->Offset < tmp) { ++DataPtr; ++counter; }
    shares[i] = counter;
  }
  shares[nCPUs] = FASTA->SequenceCount;

  pthread_mutex_t PrintLock = PTHREAD_MUTEX_INITIALIZER;
	
  /* Dispatch common information to threads */
  for (size_t i=0; i<nCPUs; ++i) {
		threads_arg[i].prf          = prf;
		threads_arg[i].FASTA        = FASTA;
		threads_arg[i].SequenceID   = NULL;
		threads_arg[i].threadId     = i;
		threads_arg[i].PrintLock    = &PrintLock;
		threads_arg[i].Compute      = Model;
		threads_arg[i].OutputType   = OutputType;
		threads_arg[i].start        = shares[i];          
		threads_arg[i].stop         = shares[i+1];
#ifdef PRF_CORE_PCRE
		threads_arg[i].regex        = regex;
#endif
  }
  
  size_t * restrict Histograms;
	size_t HistogramSize = 0UL;
	size_t DensitySize;

	/* Get the thread function */
	void * (*ThreadFct)(void*) = NULL;
	{
		const int index = GetDispatchThreadIndex(OutputType, prf);
		if (index >= 0) {
		ThreadFct = t_fap[index];
		}
		else {
			fputs("Error in the choice of thread function\n", stderr);
			return -1;
		}
	}
	
	if (OutputType->Type == HISTOGRAM) {
		if (OutputType->Specific.Histogram.CycleRatherThanScore) {
			HistogramSize = OutputType->CycleRange[1] - OutputType->CycleRange[0] + 1; 
		}
		else {
			HistogramSize = OutputType->ScoreRange[1] - OutputType->ScoreRange[0] + 1; 
		}
		Histograms = (size_t *) calloc(HistogramSize*nCPUs,sizeof(size_t));
		if (Histograms == NULL) {
				fprintf(stderr, "Unable to allocate memory for histograms, requested size was %lu bytes.\n",  HistogramSize*nCPUs*sizeof(size_t));
				return -1;
		}
		size_t * restrict ptr = Histograms;
		for (size_t i=0; i<nCPUs; ++i) {
			threads_arg[i].Histogram = ptr;
			ptr += HistogramSize;
		}
	}
	else if(OutputType->Type == DENSITY) {
		DensitySize = (OutputType->CycleRange[1] - OutputType->CycleRange[0] + 1)\
		             *(OutputType->ScoreRange[1] - OutputType->ScoreRange[0] + 1);
		Histograms = (size_t *) calloc(DensitySize*nCPUs, sizeof(size_t));
		if (Histograms == NULL) {
				fprintf(stderr, "Unable to allocate memory for histograms, requested size was %lu bytes.\n",  HistogramSize*nCPUs*sizeof(size_t));
				return -1;
		}
		size_t * restrict ptr = Histograms;
		for (size_t i=0; i<nCPUs; ++i) {
			threads_arg[i].Histogram = ptr;
			ptr += DensitySize;
		}
	}
  
  /* Dispatch to threads */
  gettimeofday(&_t0,0);
#ifdef PRF_USE_AFFINITY
	if (threads_attr) {
		for (size_t i=0; i<nCPUs; ++i) {
			if (pthread_create (&threads[i],  &threads_attr[i], ThreadFct,  (void*) &threads_arg[i]) != 0) {
				return -2;
			}
		}
	}
	else
#endif
	{
		for (size_t i=0; i<nCPUs; ++i) {
			if (pthread_create (&threads[i],  NULL, ThreadFct,  (void*) &threads_arg[i]) != 0) {
				return -2;
			}
		}
	}

	for (size_t i=0; i<nCPUs; i++) {
		pthread_join(threads[i], NULL);
	}
	gettimeofday(&_t1,0);

	pthread_mutex_destroy(&PrintLock);

	if (OutputType->Type == HISTOGRAM) {
		/* Gather thread histograms */
		size_t MissedCycles = threads_arg[0].counter;
		const size_t * restrict hptr = &Histograms[HistogramSize];
		for (size_t i=1; i<nCPUs; i++) {
			MissedCycles += threads_arg[i].counter;
			for (unsigned int j=0; j<HistogramSize; j++) Histograms[j] += hptr[j];
			hptr += HistogramSize;
		}
		if (MissedCycles) {
			fprintf(stderr, "Some sequences (%lu) bear alignment that are longer than the histogram bin number\n", MissedCycles);
		}
		
		char FName[256];
		snprintf(FName, 256, "%s.histogram", OutputType->Specific.Histogram.BaseFileName);
		FILE * out = fopen(FName, "w");
		if (out != NULL) {
			const int RangeStart = (OutputType->Specific.Histogram.CycleRatherThanScore) ? OutputType->CycleRange[0] :\
			                       OutputType->ScoreRange[0];
			for (unsigned int j=0; j<HistogramSize; j++) {
				fprintf(out, "%i\t%lu\n", RangeStart+j, Histograms[j]);
			}
			fclose(out);
		}
		else {
			fprintf(stderr, "Unable to create output histogram %s\n", FName);
			return -3;
		}
	}
	else if (OutputType->Type == DENSITY) {
		/* Gather thread histograms */
		size_t MissedCycles = threads_arg[0].counter;
		const size_t * restrict hptr = &Histograms[DensitySize];
		for (size_t i=1; i<nCPUs; i++) {
			MissedCycles += threads_arg[i].counter;
			for (unsigned int j=0; j<DensitySize; j++) Histograms[j] += hptr[j];
			hptr += DensitySize;
		}
		if (MissedCycles) {
			fprintf(stderr, "Some sequences (%lu) bear alignment that are longer than the histogram bin number\n", MissedCycles);
		}
		
		char FName[256];
		snprintf(FName, 256, "%s.density", OutputType->Specific.Histogram.BaseFileName);
		FILE * out = fopen(FName, "w");
		if (out != NULL) {
			hptr = Histograms;
			for (int i=OutputType->ScoreRange[0]; i<OutputType->ScoreRange[1]; i++) {
				for (unsigned int j=OutputType->CycleRange[0]; j<OutputType->CycleRange[1]; j++) {
					fprintf(out, "%i\t%u\t%lu\n", i, j,* hptr++);
				}
				fprintf(out, "\n");
			}
			fclose(out);
		}
		else {
			fprintf(stderr, "Unable to create output histogram %s\n", FName);
			return -3;
		}
	}
	else {
		if (OutputVerbose) {
			unsigned int AlignedSequencesCounter = threads_arg[0].counter;
			for (size_t i=1; i<nCPUs; i++) AlignedSequencesCounter += threads_arg[i].counter;
			const double t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr,"Overall there are %u aligned sequences found. These took %lf seconds to align on %li cores.\n", AlignedSequencesCounter, t, nCPUs);
		}
	}
	return SUCCESS;
}


static inline int __ALWAYS_INLINE add_task(threadpool_t * const restrict thpool, FILE * const restrict Stream)
{
	int res = 0;
	char * LineBuffer = NULL;
	size_t LineBufferSize = 0UL;
	fasta_job_t * restrict newJob = NULL;
	size_t ReadSoFarInSequence = 0UL;
	size_t HeaderLength;
	
	while (!feof(Stream)) {
		const ssize_t tmpRead = getline(&LineBuffer, &LineBufferSize, Stream);
		if (tmpRead < (ssize_t) 0) {
			if (!feof(Stream)) res = -1;
			break;
		}
		const size_t nRead = (size_t) tmpRead;
		
		/* Is this a new sequence */
		if (LineBuffer[0] == '>') {
			/* Do we end an existing one ? */
			if (newJob) {
				newJob->Fasta_Sequence.ProfileData.Length = ReadSoFarInSequence - HeaderLength;
				newJob->Fasta_Sequence.Data.Header[HeaderLength-1] = '\0';
				newJob->Fasta_Sequence.ProfileData.ProfileIndex = (unsigned char*) newJob->Fasta_Sequence.Data.Header + HeaderLength;
				pthread_mutex_lock(&thpool->jobqueue_p.rwmutex);
				jobqueue_push(&thpool->jobqueue_p, (job_t*) newJob);
				pthread_mutex_unlock(&thpool->jobqueue_p.rwmutex);
			}
			
			/* Get new job slot */
			again:;
			pthread_mutex_lock(&thpool->donequeue_p.rwmutex);
			newJob = (fasta_job_t*) jobqueue_pull(&thpool->donequeue_p);
			pthread_mutex_unlock(&thpool->donequeue_p.rwmutex);
			
			if (!newJob) {
				if (thpool->num_threads_alive > 0) {
					/* Wait a bit and try again */
					nanosleep(&polling_interval, NULL);
					goto again;
				}
				else {
					res = -2;
					break;
				}
			}
			
			/* Check memory, increase if needed */
			if (newJob->Fasta_Sequence.Size < nRead) {
				newJob->Fasta_Sequence.Data.Header = realloc(newJob->Fasta_Sequence.Data.Header, nRead*sizeof(char));
				if (newJob->Fasta_Sequence.Data.Header == NULL) { res = -3; break; }
				newJob->Fasta_Sequence.Size = nRead;
			}
			memcpy(newJob->Fasta_Sequence.Data.Header, LineBuffer, nRead);
			
			HeaderLength = nRead;
			ReadSoFarInSequence = nRead;
		}
		else if (newJob) {
			const size_t ultmp = ReadSoFarInSequence;
			/* Check memory, increase if needed */
			ReadSoFarInSequence += nRead;
			if (newJob->Fasta_Sequence.Size < ReadSoFarInSequence) {
				newJob->Fasta_Sequence.Data.Header = realloc(newJob->Fasta_Sequence.Data.Header, ReadSoFarInSequence*sizeof(char));
				if (newJob->Fasta_Sequence.Data.Header == NULL) { res = -3; break; }
				newJob->Fasta_Sequence.Size = ReadSoFarInSequence;
			}
			memcpy(&(newJob->Fasta_Sequence.Data.Header[ultmp]), LineBuffer, nRead);
		}
	}
	if (newJob) {
		newJob->Fasta_Sequence.ProfileData.Length = ReadSoFarInSequence - HeaderLength;
		newJob->Fasta_Sequence.Data.Header[HeaderLength-1] = '\0';
		newJob->Fasta_Sequence.ProfileData.ProfileIndex = (unsigned char*) newJob->Fasta_Sequence.Data.Header + HeaderLength;
		pthread_mutex_lock(&thpool->jobqueue_p.rwmutex);
		jobqueue_push(&thpool->jobqueue_p, (job_t*) newJob);
		pthread_mutex_unlock(&thpool->jobqueue_p.rwmutex);
	}
	
	free(LineBuffer);
	return res;	
}

int dispatchStreamFASTA(const struct Profile * const restrict prf,
                        FILE * const restrict Stream,
#ifdef PRF_CORE_PCRE
                        struct RegEx * const restrict regex,
#endif
                        const OutputType_t * const restrict OutputType,
                        const size_t nCPUs)
{
	int res;	
	/*************************************************************************/
	/*                         CONSISTENCY CHECK                             */
	/*************************************************************************/

	/*************************************************************************/
	/*                         CHOOSE THREAD FUNCTION                        */
	/*************************************************************************/
#if (defined(PRF_CORE_STD) || defined(PRF_CORE_REPEAT))
#if (defined(PRF_CORE_STD) && defined(PRF_CORE_REPEAT))
	const Compute_t * CoreCompute = (prf->isCircular) ? &Repeat_sse41 : &Standard_sse41;
#elif defined(PRF_CORE_STD)
	const Compute_t * CoreCompute = &Standard_sse41;
#else
	const Compute_t * CoreCompute = &Repeat_sse41;
#endif
#else
#error "Dispatching procedures requires at lest STD or REPEAT methods"
#endif

 	void * (*ThreadFct)(threadarg_t* const restrict) = NULL;
	{
		const int index = GetDispatchThreadIndex(OutputType, prf);
		if (index >= 0) {
		ThreadFct = tp_fap[index];
		}
		else {
			fputs("Error in the choice of thread function\n", stderr);
			return -1;
		}
	}

	/*************************************************************************/
	/*                       EXTRA VARIABLES FOR THREADS                     */
	/*************************************************************************/
	size_t HistogramsMemSize = 0UL;
	if (OutputType->Type == HISTOGRAM) {
		if (OutputType->Specific.Histogram.CycleRatherThanScore) {
			HistogramsMemSize = OutputType->CycleRange[1] - OutputType->CycleRange[0] + 1; 
		}
		else {
			HistogramsMemSize = OutputType->ScoreRange[1] - OutputType->ScoreRange[0] + 1; 
		}
	}
	else if(OutputType->Type == DENSITY) {
		HistogramsMemSize = (OutputType->CycleRange[1] - OutputType->CycleRange[0] + 1)\
		                   *(OutputType->ScoreRange[1] - OutputType->ScoreRange[0] + 1);
	}
	
	size_t * const Histograms = (size_t *) calloc(HistogramsMemSize*nCPUs, sizeof(size_t));
	if (Histograms == NULL) {
		fprintf(stderr, "Unable to allocate memory for histograms, requested size was %lu bytes.\n",
						HistogramsMemSize*nCPUs*sizeof(size_t));
		return -1;
	}
	size_t * const Missed = (size_t *) calloc(nCPUs, sizeof(size_t));
	if (Missed == NULL) {
		fprintf(stderr, "Unable to allocate memory for missed values, requested size was %lu bytes.\n",
						nCPUs*sizeof(size_t));
		if (Histograms) free(Histograms);
		return -2;
	}

	/*************************************************************************/
	/*                    SET THE THREADS COMMON VARIABLES                   */
	/*************************************************************************/
	common_t common = {
		.profile = prf,
		.PrintLock = PTHREAD_MUTEX_INITIALIZER,
		.Histograms = Histograms,
		.Counters = Missed,
		.Compute = CoreCompute,
		.OutputType = OutputType
#ifdef PRF_CORE_PCRE
		, .regex = regex
#endif
	};
	
	/*************************************************************************/
	/*                          PREPARE THREAD POOL                          */
	/*************************************************************************/  
	threadpool_t * const restrict thpool = createThreadPool(ThreadFct,
#ifdef PRF_USE_AFFINITY
																													NULL,
#endif
                                                          sizeof(fasta_job_t),
                                                          (void*) &common,
                                                          nCPUs,0);
	if ( !thpool ) {
		if (Histograms) free(Histograms);
		if (Missed) free(Missed);
		return -3;
	}
	
	/*************************************************************************/
	/*                GETTING AND ANALYSING INPUT BAM FILE                   */
	/*************************************************************************/

	
	/*************************************************************************/
	/*                        ADD TASKS TO JOB QUEUE                         */
	/*************************************************************************/
	struct timeval _t0, _t1;
	gettimeofday(&_t0,0);
	
	res = add_task(thpool, Stream); 
	
	if (res != 0 ) {
		fprintf(stderr, "Master failed to send jobs, error code %i\n", res);
		destroyThreadPool(thpool);
		pthread_mutex_destroy(&(common.PrintLock));
		if (Histograms) free(Histograms);
		if (Missed) free(Missed);
		return -4;
	}
	
	
	/*************************************************************************/
	/*                       WAIT FOR TASK TO FINISH                         */
	/*************************************************************************/
	thpool_wait(thpool);
	gettimeofday(&_t1,0);
	
	if (OutputVerbose) {
		const double t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
		fprintf(stderr,"This took %lf seconds to crunch on %lu cores.\n", t, nCPUs);
	}
			
	/*************************************************************************/
	/*                  TELL THREADS TO FLUSH AND TERMINATE                  */
	/*************************************************************************/


	/*************************************************************************/
	/*                          GATHER_EXTRA_DATA                            */
	/*************************************************************************/ 
	if (HistogramsMemSize != 0) {
		size_t MissedValues = Missed[0];
		const size_t * restrict hptr = &Histograms[HistogramsMemSize];
		for (size_t i=1;i<nCPUs;i++) {
			MissedValues += Missed[i];
			for (size_t j=0; j<HistogramsMemSize; j++) Histograms[j] += hptr[j];
			hptr += HistogramsMemSize;
		}
		if (MissedValues) {
			fprintf(stderr,
							"Some sequences (%lu) bear alignment that are longer than the histogram bin number\n",
							MissedValues);
		}
		char FName[256];
		hptr = Histograms;
		if (OutputType->Type == HISTOGRAM) {
			snprintf(FName, 256, "%s.histogram", OutputType->Specific.Histogram.BaseFileName);
			FILE * const out = fopen(FName, "w");
			const int RangeStart = (OutputType->Specific.Histogram.CycleRatherThanScore) ? OutputType->CycleRange[0] :\
			                        OutputType->ScoreRange[0];
			if (out != NULL) {
				for (size_t j=0; j<HistogramsMemSize; j++) {
					fprintf(out, "%i\t%lu\n", RangeStart+(int)j, Histograms[j]);
				}
				fclose(out);
			}
			else {
				fprintf(stderr, "Unable to create output histogram %s\n", FName);
				res = -5;
			}
		}
		else if (OutputType->Type == DENSITY) {
			snprintf(FName, 256, "%s.density", OutputType->Specific.Histogram.BaseFileName);
			FILE * const out = fopen(FName, "w");
			if (out != NULL) {
				hptr = Histograms;
				for (int i=OutputType->ScoreRange[0]; i<OutputType->ScoreRange[1]; i++) {
					for (unsigned int j=OutputType->CycleRange[0]; j<OutputType->CycleRange[1]; j++) {
						fprintf(out, "%i\t%u\t%lu\n", i, j,* hptr++);
					}
					fprintf(out, "\n");
				}
				fclose(out);
			}
			else {
				fprintf(stderr, "Unable to create output histogram %s\n", FName);
				res = -5;
			}
		}
	}
	

	/*************************************************************************/
	/*                    CHECK FOR THREADS ERROR                            */
	/*************************************************************************/

			
	/*************************************************************************/
	/*                      DESTROY THE THREAD POOL                          */
	/*************************************************************************/

	destroyThreadPool(thpool);
	pthread_mutex_destroy(&(common.PrintLock));
	
	if (Histograms) free(Histograms);
	if (Missed) free(Missed);
	return res;
}

int dispatchStreamFASTAExt(FILE * const restrict Stream, threadpool_t * const restrict thpool )
{
	int res;
	/*************************************************************************/
	/*                        ADD TASKS TO JOB QUEUE                         */
	/*************************************************************************/
	struct timeval _t0, _t1;
	gettimeofday(&_t0,0);
	while (thpool->num_threads_alive == 0) nanosleep(&polling_interval, NULL);
	res = add_task(thpool, Stream);
	
	if (res != 0 ) {
		fprintf(stderr, "Master failed to send jobs, error code %i\n", res);
		return -4;
	}
	
	
	/*************************************************************************/
	/*                       WAIT FOR TASK TO FINISH                         */
	/*************************************************************************/
	thpool_wait(thpool);
	gettimeofday(&_t1,0);
	
	if (OutputVerbose) {
		const double t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
		fprintf(stderr,"This took %lf seconds to crunch on %i cores.\n", t, thpool->num_threads);
	}
			
	return res;
}
