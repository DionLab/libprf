#ifndef _PB_THREADS_H
#define _PB_THREADS_H
#include "pfProfile.h"
#include "pfCompute.h"
#ifdef PRF_CORE_PCRE
#include "pfRegexp.h"
#endif
#include "pfOutput.h"
#include "threadpool.h"
#include "pb_common.h"

#if defined(__USE_WINAPI__)
#define pthread_exit(x) return (LPDWORD) x
#define THREADPOOL_FUNCTION(Fct) LPDWORD WINAPI Fct(_In_ LPVOID _Data)
#else
#define THREADPOOL_FUNCTION(Fct) void* Fct(threadarg_t * const restrict _Data)
#endif
/*
 ************************************************************************************************
 *                                       STRUCTURES                                             *
 ************************************************************************************************
 */

typedef struct COMMON {
	const struct Profile * profile;
	pthread_mutex_t PrintLock;
	const Compute_t * restrict Compute;           /* pointer to core functions */
	size_t * restrict Histograms;
	size_t * restrict Counters;
	const OutputType_t * restrict OutputType;
#ifdef PRF_CORE_PCRE
  struct RegEx * restrict regex;
#endif
} common_t;

/* WARNING: This must be in agreement with job_t defined in threadpool.h for inheritance to function */
struct pb_job_s {
	struct pb_job_s * prev;          /* pointer to previous job   */
	HoleReads_t Reads;
	HoleData_t HReg;
}; 
typedef struct pb_job_s pb_job_t;

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
extern THREADPOOL_FUNCTION( (* * tp_pbsp) );
extern THREADPOOL_FUNCTION( (* * tp_pbzp) );
THREADPOOL_FUNCTION(tp_pb_filtertest);
/*
 ************************************************************************************************
 *                                     INLINE FUNCTIONS                                         *
 ************************************************************************************************
 */

#endif /* _PB_THREADS */
