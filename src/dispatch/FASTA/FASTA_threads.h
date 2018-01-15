#ifndef _FASTA_THREADS_H
#define _FASTA_THREADS_H
#include <stdlib.h>
#include <pthread.h>
#include "pfCompute.h"
#include "pfSequence.h"
#include "threadpool.h"
#include "pfOutput.h"

#if defined(__USE_WINAPI__)
#define pthread_exit(x) return (LPDWORD) x
#define THREAD_FUNCTION(F) LPDWORD WINAPI F(_In_ LPVOID _Data)
#else
#define THREAD_FUNCTION(Fct) void* Fct(void * _Data)
#define THREADPOOL_FUNCTION(Fct) void* Fct(threadarg_t * const restrict _Data)
#endif
/*
 ************************************************************************************************
 *                                       STRUCTURES                                             *
 ************************************************************************************************
 */

struct ThreadData {
	const struct Profile * prf;
	const FASTAStructure * FASTA;
#if defined(PRF_USE_MMAP) && defined(MMAP_DEBUG)
	size_t * maplength;
#endif
#if defined(PRF_CORE_PCRE)
	const struct RegEx * regex;
#endif
	unsigned int * SequenceID;
	const Compute_t * Compute;           /* pointer to core functions */
	pthread_mutex_t * PrintLock;         /* screen locker for printing */
	const OutputType_t * OutputType;
	size_t start;
	size_t stop;
	unsigned long counter; /* WARNING: on input in filter phase counter == 0 means stop when cutoff reached, 1 -> compute real filter value */
	size_t * restrict Histogram;
	size_t threadId;
// 	enum Version version;
};

typedef struct COMMON {
	const struct Profile * profile;
	const OutputType_t * OutputType;
	const Compute_t * Compute;           /* pointer to core functions */
	pthread_mutex_t PrintLock;
	size_t * Histograms;
	size_t * Counters;
#ifdef PRF_CORE_PCRE
  const struct RegEx * regex;
#endif
} common_t;

/* WARNING: This must be in agreement with job_t defined in threadpool.h for inheritance to function */
struct fasta_job_s {
	struct fasta_job_s * prev;          /* pointer to previous job   */
	Sequence_t Fasta_Sequence;
};

typedef struct fasta_job_s fasta_job_t;

// struct ThreadArrayData {
// 	const struct Profile * * prf;
// 	const PFSequence * restrict PFSeq;
// 	const char * restrict OutputFileName;
// 	float CutOff; 
// 	size_t MaxProfileLength;
// 	size_t start;
// 	size_t stop;
// 	size_t threadId;
// 	enum Version version;
// };

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

extern THREAD_FUNCTION( (* * t_fap) );
extern THREADPOOL_FUNCTION( (* * tp_fap) );

#ifdef PRF_CORE_PCRE
extern THREAD_FUNCTION( (* * t_far) );
#endif

#endif /*_FASTA_THREADS_H*/
