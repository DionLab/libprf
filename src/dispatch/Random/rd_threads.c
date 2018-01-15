#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#define __USE_INLINE_FUNCTIONS__
#include "pfCompute.h"
#include "pfSequence.h"
#include "pfOutput.h"
#include "rd_threads.h"
#include "matrix_functions.h"

/*
 ************************************************************************************************
 *                                         GLOBALS                                              *
 ************************************************************************************************
 */

/*
 ************************************************************************************************
 *                                        FUNCTIONS                                             *
 ************************************************************************************************
 */

#if defined(PRF_CORE_PCRE)
THREAD_FUNCTION(threadpool_random_regex_histogram)
{
		int res = SUCCESS;
// 	  printf("Thread started fct THREAD_FCT_NAME\n");
		
	  /*************************************************************************/
    /*                          GET COMMON DATA                              */
    /*************************************************************************/
		const common_t * const common  = ((ThreadPool_Arg_t*)_Data)->common;
		const unsigned int SeqLength   = common->RD->Length + 3*common->RD->nCAGSuffix;
		const pcre * const restrict rg = common->regex->regexCompiled[0];
		
		/*************************************************************************/
    /*                         GET PRIVATE DATA                              */
    /*************************************************************************/
		unsigned int * const restrict Histogram = ((ThreadPool_Arg_t*) _Data)->histogram;
		
	  /*************************************************************************/
    /*                          ALLOCATE MEMORY                              */
    /*************************************************************************/
			
		/* Builds RevComp index */
		char * const restrict RevComp = (char*) malloc(SeqLength*sizeof(char));
		if (RevComp == NULL) {
			res = -1;
			goto FIN;
		}
		
    /*************************************************************************/
    /*                          CONFIGURE WORKER                             */
    /*************************************************************************/
    
    /* Assure all threads have been created before starting serving */
    threadpool_t* thpool = ((ThreadPool_Arg_t*)_Data)->thpool;
    
    /* Register signal handler */
//     struct sigaction act;
//     act.sa_handler = thread_hold;
//     if (sigaction(SIGUSR1, &act, NULL) == -1) {
// 	    fprintf(stderr, "thread_do(): cannot handle SIGUSR1");
//     }
    
    /* Mark thread as alive (initialized) */
    pthread_mutex_lock(&thpool->thcount_lock);
    thpool->num_threads_alive += 1;
    pthread_mutex_unlock(&thpool->thcount_lock);
    
    
    // Initialize some other data
    rd_job_t * task = NULL;
		size_t MissedCycles=0;
    
    while(threads_keepalive) {

			bsem_wait(thpool->jobqueue_p->has_jobs);
			
			if (threads_keepalive){
				pthread_mutex_lock(&thpool->thcount_lock);
				thpool->num_threads_working++;
				pthread_mutex_unlock(&thpool->thcount_lock);

				/* Read job from queue and execute it */
				pthread_mutex_lock(&thpool->jobqueue_p->rwmutex);
				task = (rd_job_t*) jobqueue_pull(thpool);
				pthread_mutex_unlock(&thpool->jobqueue_p->rwmutex);

				/***************************************************************/
				/*                     START THE TASK                          */
				/***************************************************************/
				if (task) {
					char * const restrict Sequence = task->Reads;
					
					/* Run the regex engine */
					int ovector[2*8];
					unsigned int offset = 0;
					int rc;
					size_t count = 0;
					while (offset < SeqLength && (rc = pcre_exec(rg, 0, Sequence, SeqLength, offset, 0, ovector, 16)) >= 0)
					{
						count += rc;
						offset = ovector[1];
					}
					
					{
						int k = (int) SeqLength;
						char * restrict cptr = RevComp;
						while (--k>=0) {
							switch(Sequence[k]) {
								case 'A': *cptr = 'T'; break;
								case 'T': *cptr = 'A'; break;
								case 'G': *cptr = 'C'; break;
								case 'C': *cptr = 'G'; break;
								default:
									*cptr = '?';
							}
							cptr++;
						}
					}
					
					offset = 0;
					size_t rccount = 0;
					
					while (offset < SeqLength && (rc = pcre_exec(rg, 0, RevComp, SeqLength, offset, 0, ovector, 16)) >= 0)
					{
						rccount += rc;
						offset = ovector[1];
					}
					
					count = (rccount > count) ? rccount : count;
					
					if (count < HistogramSize)
						Histogram[count]++;
					else
						MissedCycles++;
				}

				pthread_mutex_lock(&thpool->thcount_lock);
				thpool->num_threads_working--;
				pthread_mutex_unlock(&thpool->thcount_lock);
			}
    }
     
    ((ThreadPool_Arg_t*) _Data)->counter = (unsigned long) MissedCycles;
     
    FIN:
    pthread_mutex_lock(&thpool->thcount_lock);
    thpool->num_threads_alive--;
    pthread_mutex_unlock(&thpool->thcount_lock);
    
		free(RevComp);
		
    ((ThreadPool_Arg_t*)_Data)->returnState = res;
    if( res != SUCCESS) {
      fprintf(stderr, "Thread %lu error code %i\n", ((ThreadPool_Arg_t*)_Data)->ID, res);
      if (task) {
				pthread_mutex_lock(&thpool->thcount_lock);
				thpool->num_threads_working--;
				pthread_mutex_unlock(&thpool->thcount_lock);
				jobqueue_clear(thpool);
      }
      threads_keepalive = 0;
    }
    
//     printf("Thread done\n");
		return 0;
//     pthread_exit( (void*) ((intptr_t)res));
}

#endif // PRF_CORE_PCRE

THREAD_FUNCTION(threadpool_random_repeat_histogram)
{
		int res = SUCCESS;
		unsigned char RevComp_Alphabet_Mapping[ALPHABET_SIZE+2] __attribute__((aligned(16)));
// 	  printf("Thread started fct THREAD_FCT_NAME\n");
		
	  /*************************************************************************/
    /*                          GET COMMON DATA                              */
    /*************************************************************************/
		const common_t * const common  = ((ThreadPool_Arg_t*)_Data)->common;
		const struct Profile * const restrict prf = common->profile;
		const unsigned int SeqLength   = common->RD->Length + 3*common->RD->nCAGSuffix;

		
		/*************************************************************************/
    /*                         GET PRIVATE DATA                              */
    /*************************************************************************/
		unsigned int * const restrict Histogram = ((ThreadPool_Arg_t*) _Data)->histogram;
		
	  /*************************************************************************/
    /*                          ALLOCATE MEMORY                              */
    /*************************************************************************/
		const size_t LargestSequence = SeqLength+1;
		const size_t prfLength = prf->Length;
		const size_t WorkSize  = prf->Length + 1;
		
		union lScores * const restrict matrix   = _mm_malloc(WorkSize*LargestSequence*sizeof(union lScores), 64);
		union lScores * const restrict rvmatrix = _mm_malloc(WorkSize*LargestSequence*sizeof(union lScores), 64);
		int * const restrict WORK               = _mm_malloc(2*WorkSize*sizeof(union lScores)+63,64);
		char * restrict SequenceIndex     = (char*) malloc(LargestSequence*sizeof(char));
		if ( rvmatrix == NULL || matrix == NULL || WORK == NULL || SequenceIndex == NULL) {
			res = -1; 
			goto FIN;
		}
	
		PFSequence PFSeq;
		PFSeq.ProfileIndex = SequenceIndex;
	
		/* Builds RevComp index */
		memcpy(RevComp_Alphabet_Mapping, prf->Alphabet_Mapping, (ALPHABET_SIZE+2)*sizeof(unsigned char));
		RevComp_Alphabet_Mapping[(int) 'A' - (int) 'A'] =  prf->Alphabet_Mapping[(int) 'T' - (int) 'A' ];
		RevComp_Alphabet_Mapping[(int) 'C' - (int) 'A'] =  prf->Alphabet_Mapping[(int) 'G' - (int) 'A' ];
		RevComp_Alphabet_Mapping[(int) 'G' - (int) 'A'] =  prf->Alphabet_Mapping[(int) 'C' - (int) 'A' ];
		RevComp_Alphabet_Mapping[(int) 'T' - (int) 'A'] =  prf->Alphabet_Mapping[(int) 'A' - (int) 'A' ];
		
    /*************************************************************************/
    /*                          CONFIGURE WORKER                             */
    /*************************************************************************/
    
    /* Assure all threads have been created before starting serving */
    threadpool_t* thpool = ((ThreadPool_Arg_t*)_Data)->thpool;
    
    /* Register signal handler */
//     struct sigaction act;
//     act.sa_handler = thread_hold;
//     if (sigaction(SIGUSR1, &act, NULL) == -1) {
// 	    fprintf(stderr, "thread_do(): cannot handle SIGUSR1");
//     }
    
    /* Mark thread as alive (initialized) */
    pthread_mutex_lock(&thpool->thcount_lock);
    thpool->num_threads_alive += 1;
    pthread_mutex_unlock(&thpool->thcount_lock);
    
    
    // Initialize some other data
    rd_job_t * task = NULL;
		size_t MissedCycles=0;
    
    while(threads_keepalive) {

			bsem_wait(thpool->jobqueue_p->has_jobs);
			
			if (threads_keepalive){
				pthread_mutex_lock(&thpool->thcount_lock);
				thpool->num_threads_working++;
				pthread_mutex_unlock(&thpool->thcount_lock);

				/* Read job from queue and execute it */
				pthread_mutex_lock(&thpool->jobqueue_p->rwmutex);
				task = (rd_job_t*) jobqueue_pull(thpool);
				pthread_mutex_unlock(&thpool->jobqueue_p->rwmutex);

				/***************************************************************/
				/*                     START THE TASK                          */
				/***************************************************************/
				if (task) {
					char * const restrict seq = task->Reads;
					
					/* Copy sequence to local space index */ 
					memcpy(PFSeq.ProfileIndex, seq, SeqLength);
					PFSeq.Length = SeqLength;
					
					/* Translate into indices */
					TranslateSequenceToIndex(&PFSeq, prf->Alphabet_Mapping);
					
					/* Build matrix */
					Repeat_sse41.BuildMatrix(prf, PFSeq.ProfileIndex, matrix, WORK, NULL, 0, SeqLength);
					
					/* Copy sequence to local space index */ 
					memcpy(PFSeq.ProfileIndex, seq, SeqLength);
					
					/* Translate into indices */
					TranslateSequenceToIndex(&PFSeq, RevComp_Alphabet_Mapping);
					ReverseTranslatedSequence(&PFSeq);
					
					/* Build reverse matrix */
					Repeat_sse41.BuildMatrix(prf, PFSeq.ProfileIndex, rvmatrix, WORK, NULL, 0, SeqLength);
					
					/* Seek alignments */
					ScoreCycle_t SC, rvSC;
					GetBestScoreAndCycles(matrix, SeqLength, prfLength, &SC);
					GetBestScoreAndCycles(rvmatrix, SeqLength, prfLength, &rvSC);

					const int lmaxCycle = (SC.Cycles > rvSC.Cycles) ? SC.Cycles : rvSC.Cycles;
					
					if (lmaxCycle < HistogramSize)
						Histogram[lmaxCycle]++;
					else
						MissedCycles++;
				}

				pthread_mutex_lock(&thpool->thcount_lock);
				thpool->num_threads_working--;
				pthread_mutex_unlock(&thpool->thcount_lock);
			}
    }
     
    ((ThreadPool_Arg_t*) _Data)->counter = (unsigned long) MissedCycles;
     
    FIN:
    pthread_mutex_lock(&thpool->thcount_lock);
    thpool->num_threads_alive--;
    pthread_mutex_unlock(&thpool->thcount_lock);
    
		if (WORK) _mm_free(WORK);
		if (matrix) _mm_free(matrix);
		if (rvmatrix) _mm_free(rvmatrix);
		if (SequenceIndex) free(SequenceIndex);
		
    ((ThreadPool_Arg_t*)_Data)->returnState = res;
    if( res != SUCCESS) {
      fprintf(stderr, "Thread %lu error code %i\n", ((ThreadPool_Arg_t*)_Data)->ID, res);
      if (task) {
				pthread_mutex_lock(&thpool->thcount_lock);
				thpool->num_threads_working--;
				pthread_mutex_unlock(&thpool->thcount_lock);
				jobqueue_clear(thpool);
      }
      threads_keepalive = 0;
    }
    
//     printf("Thread done\n");
		return 0;
//     pthread_exit( (void*) ((intptr_t)res));
}
