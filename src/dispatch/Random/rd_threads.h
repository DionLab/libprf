#ifndef _RD_THREADS_H
#define _RD_THREADS_H
#include "pfProfile.h"
#ifdef PRF_CORE_PCRE
#include "pfRegexp.h"
#endif
#include "threadpool.h"

#if defined(__USE_WINAPI__)
#define pthread_exit(x) return (LPDWORD) x
#define THREAD_FUNCTION(F) LPDWORD WINAPI F(_In_ LPVOID _Data)
#else
#define THREAD_FUNCTION(F) void* F(void * _Data)
#endif
/*
 ************************************************************************************************
 *                                       STRUCTURES                                             *
 ************************************************************************************************
 */

typedef struct COMMON {
	const struct Profile * profile;
#ifdef PRF_CORE_PCRE
  const struct RegEx * regex;
#endif
	const RandomData_t * RD;
} common_t;

/* WARNING: This must be in agreement with job_t defined in threadpool.h for inheritance to function */
struct rd_job_s {
	struct rd_job_s * prev;          /* pointer to previous job   */
	char * restrict Reads;
	size_t Length;
}; 
typedef struct rd_job_s rd_job_t;

typedef struct ThreadPool_Arg {
    size_t ID;                           /* friendly id               */
    pthread_t pthread;                   /* pointer to actual thread  */
    struct threadpool* thpool;           /* access to thpool          */
    void * common;                       /* common data to threadpool */
    unsigned int * restrict histogram;   /* individual histogram      */
    unsigned long counter;               /* number of counter of histogram size */
    int returnState;                     /* error return              */  
} ThreadPool_Arg_t;

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
THREAD_FUNCTION(threadpool_random_regex_histogram);
// THREAD_FUNCTION(threadpool_random_regex);
#endif // PRF_CORE_PCRE
THREAD_FUNCTION(threadpool_random_repeat_histogram);

/*
 ************************************************************************************************
 *                                     INLINE FUNCTIONS                                         *
 ************************************************************************************************
 */

#endif /* _RD_THREADS */
