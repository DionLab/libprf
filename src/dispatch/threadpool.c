#define _GNU_SOURCE
#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>
 
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <pthread.h>

#include "threadpool.h" 

/* *****************************************************************************************************
 * Author:       Johan Hanssen Seferidis
 * License:          MIT
 * Description:  Library providing a threading pool where you can add
 *               work. For usage, check the thpool.h file or README.md
 *
 *//** @file threadpool.h *//*
 * 
 ********************************/

#ifdef THREADPOOL_DEBUG
#define THREADPOOL_DEBUG 1
#else
#define THREADPOOL_DEBUG 0
#endif

#define MAX_NANOSEC 999999999
#define CEIL(X) ((X-(int)(X)) > 0 ? (int)(X+1) : (int)(X))

/* ========================== THREADPOOL ============================ */
threadpool_t * createThreadPool(void* (*Fct)(threadarg_t * const restrict),
#ifdef PRF_USE_AFFINITY
                                const cpu_set_t * const restrict affinities,
#endif
                                const size_t jobvarsize,
                                void * restrict common,
                                const int nThreads, 
                                const int onHold)
{
	threadpool_t * const restrict thpool = (threadpool_t*) malloc(sizeof(threadpool_t));
	if (thpool == NULL) return NULL;
	
	thpool->threads = (pthread_t*) malloc(nThreads*sizeof(pthread_t));
	if (thpool->threads == NULL) goto err1;
	
	thpool->threads_on_hold   = onHold;
	thpool->threads_keepalive = 1;
	
	thpool->num_threads         = nThreads;
	thpool->num_threads_alive   = 0;
	thpool->num_threads_working = 0;
	
	pthread_mutex_init(&(thpool->thcount_lock), NULL);
	
	if (jobqueue_init(&thpool->jobqueue_p) == -1) goto err2;

	if (jobvarsize > 0 ) {
		if (jobqueue_init(&thpool->donequeue_p) == -1) goto err3;
		
		/* Create free memory jobs and place them into donequeue, EXTRA could be NICE ?*/
		void * restrict Jobs = (void*) calloc((2*nThreads),jobvarsize);
		if (Jobs == NULL) goto err4;
		
		thpool->jobs = (job_t*) Jobs;
		for(int iJob=0; iJob<(2*nThreads); iJob++) {
			jobqueue_push(&thpool->donequeue_p, (job_t*) Jobs);
			Jobs += jobvarsize;
		}
	}
	else {
		thpool->jobs = NULL;
	}
	
	threadarg_t* args = (threadarg_t*) malloc(nThreads*sizeof(threadarg_t));
	if (args == NULL) {
		goto err6;
	}
	
	thpool->args = args;
	
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, (size_t) (2*1024*1024));
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
// 	pthread_attr_setguardsize(&attr, (size_t) (2*4096));
#ifdef PRF_USE_AFFINITY
	if (affinities) {
		for (int iThread = 0; iThread<nThreads; iThread++) {
			pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &affinities[iThread]);
			args[iThread].thpool = thpool;
			args[iThread].common = common;
			args[iThread].threadID = iThread;
			if (pthread_create(&(thpool->threads[iThread]), NULL, (void* (*)(void*)) Fct, &args[iThread]) != 0) {
				while (--iThread > 0) {
					pthread_kill(thpool->threads[iThread], SIGKILL);
				}
				goto err6;
			}
		}
	}
	else 
#endif
	{
		for (int iThread = 0; iThread<nThreads; iThread++) {
			args[iThread].thpool = thpool;
			args[iThread].common = common;
			args[iThread].threadID = iThread;
			if (pthread_create(&(thpool->threads[iThread]), &attr, (void* (*)(void*)) Fct, &args[iThread]) != 0) {
				while (--iThread > 0) {
					pthread_kill(thpool->threads[iThread], SIGKILL);
				}
				goto err6;
			}
		}
	}
	
	return thpool;
	err6:
		free(thpool->args);
	err5:
		jobqueue_destroy(&thpool->jobqueue_p);
	err4:
		jobqueue_destroy(&thpool->donequeue_p);
	err3:
		pthread_mutex_destroy(&(thpool->thcount_lock));
	err2:
		free(thpool->threads);
	err1:
		free(thpool);
	return NULL;
}

void destroyThreadPool(threadpool_t * restrict thpool)
{
	/* No need to destory if it's NULL */
	if (thpool == NULL) return ;

	/* End each thread 's infinite loop */
	thpool->threads_keepalive = 0;

	/* Give one second to kill idle threads */
	const double TIMEOUT = 1.0;
	time_t start, end;
	double tpassed = 0.0;
	time (&start);
	while (tpassed < TIMEOUT && thpool->num_threads_alive){
		bsem_post_all(thpool->jobqueue_p.has_items);
		time (&end);
		tpassed = difftime(end,start);
	}

	/* Poll remaining threads */
	while (thpool->num_threads_alive){
		bsem_post_all(thpool->jobqueue_p.has_items);
		sleep(1);
	}

	/* Job queue cleanup */
	jobqueue_destroy(&thpool->jobqueue_p);
	if (thpool->jobs) jobqueue_destroy(&thpool->donequeue_p);
	
	if (thpool->jobs) {
		free(thpool->jobs);
	}
	
	free(thpool->threads);
	free(thpool->args);
	free(thpool);
	thpool = NULL;
}


/* Wait until all jobs have finished */
void thpool_wait(threadpool_t * const thpool_p){

	/* Continuous polling */
	double timeout = 1.0;
	time_t start, end;
	double tpassed = 0.0;
	time (&start);
	while (tpassed < timeout && (thpool_p->jobqueue_p.len || thpool_p->num_threads_working))
	{
		time (&end);
		tpassed = difftime(end,start);
	}

	/* Exponential polling */
	long init_nano =  1; /* MUST be above 0 */
	long new_nano;
	double multiplier = 1.01;
	int  max_secs = 20;

	struct timespec polling_interval;
	polling_interval.tv_sec  = 0;
	polling_interval.tv_nsec = init_nano;

	while (thpool_p->jobqueue_p.len || thpool_p->num_threads_working)
	{
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

	/* Fall back to max polling */
	while (thpool_p->jobqueue_p.len || thpool_p->num_threads_working) sleep(max_secs);
}


/* Destroy the threadpool */
// void thpool_destroy(threadpool_t * const thpool_p){
// 
// 	volatile int threads_total = thpool_p->num_threads_alive;
// 
// 	/* End each thread 's infinite loop */
// 	thpool_p->threads_keepalive = 0;
// 
// 	/* Give one second to kill idle threads */
// 	double TIMEOUT = 1.0;
// 	time_t start, end;
// 	double tpassed = 0.0;
// 	time (&start);
// 	while (tpassed < TIMEOUT && thpool_p->num_threads_alive){
// 		bsem_post_all(thpool_p->jobqueue_p.has_items);
// 		time (&end);
// 		tpassed = difftime(end,start);
// 	}
// 
// 	/* Poll remaining threads */
// 	while (thpool_p->num_threads_alive){
// 		bsem_post_all(thpool_p->jobqueue_p.has_items);
// 		sleep(1);
// 	}
// 
// 	/* Job queue cleanup */
// 	jobqueue_destroy(&(thpool_p->jobqueue_p));
// 	jobqueue_destroy(&(thpool_p->donequeue_p));
// }

/* Pause threads */
// void thpool_pause(threadpool_t * const thpool_p) {
//         int n;
//         ThreadPool_Arg_t * const restrict threads = thpool_p->threads;
//         for (n=0; n < thpool_p->num_threads_alive; n++){
//                 pthread_kill(threads[n].pthread, SIGUSR1);
//         }
// }

/* ============================ THREAD ============================== */

/* Sets the calling thread on hold */
void thread_hold (threadpool_t * const thpool_p) {
	thpool_p->threads_on_hold = 1;
	while (thpool_p->threads_on_hold) sleep(1);
}

/* ============================ JOBS ============================== */
