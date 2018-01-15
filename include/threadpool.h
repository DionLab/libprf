#ifndef _THREAD_POOL_H
#define _THREAD_POOL_H
#include <stdlib.h> 
#include <stdio.h>
#include <pthread.h>
#include <signal.h>
#include "pfConfig.h"
 
/********************************** 
 * @author      Johan Hanssen Seferidis
 * License:     MIT
 * 
 **********************************/

/* ========================== STRUCTURES ============================ */

/* Binary semaphore */
typedef struct bsem {
	pthread_mutex_t mutex;
	pthread_cond_t   cond;
	int v;
} bsem;

/* Generic job */
typedef struct job {
	struct job * prev;
} job_t;

/* Job queue */
typedef struct jobqueue{
	pthread_mutex_t rwmutex;            /* used for queue r/w access */
	job_t  *front;                      /* pointer to front of queue */
	job_t  *rear;                       /* pointer to rear  of queue */
	bsem *has_items;                    /* flag as binary semaphore  */
	int   len;                          /* number of jobs in queue   */
} jobqueue_t;

/* Threadpool */
typedef struct threadpool {
	pthread_t * threads;                   /* pointer to threads                           */
	job_t * jobs;                          /* pointer to all jobs memory, only for freeing */
	struct threadArgs * args;              /* allocated spot for thread arguments          */
	pthread_mutex_t  thcount_lock;         /* used for thread count etc                    */
	jobqueue_t jobqueue_p;                 /* pointer to the job queue                     */
	jobqueue_t donequeue_p;                /* pointer to the free jobs queue               */
	int num_threads;                       /* total number of threads in the pool          */
	volatile int num_threads_alive;        /* threads currently alive                      */
	volatile int num_threads_working;      /* threads currently working                    */
	volatile int threads_keepalive;        /* Used to stop threads                         */
	volatile int threads_on_hold;
} threadpool_t;

typedef struct threadArgs {
	threadpool_t * thpool;    /* thread pool pointer                 */
	unsigned int threadID;    /* thread identification number [0,N]  */
	void * common;            /* place holder for common thread data */
	void * returnState;       /* returning values                    */
} threadarg_t;

/* =================================== API ======================================= */
threadpool_t * createThreadPool(void* (*Fct)(threadarg_t * const restrict),
#ifdef PRF_USE_AFFINITY
                                const cpu_set_t * const restrict affinities,
#endif
                                const size_t jobvarsize, 
                                void * const restrict common,
                                const int nThreads,
                                const int onHold);
void destroyThreadPool(threadpool_t * restrict thpool);
void thpool_wait(threadpool_t * const);
// void thpool_destroy(threadpool_t * const);
// void thpool_pause(threadpool_t * const thpool_p);

static __inline__ void __ALWAYS_INLINE thpool_resume(threadpool_t * const  thpool_p) {
	thpool_p->threads_on_hold = 0;
}

/* ======================== SYNCHRONISATION ========================= */

/* Init semaphore to 1 or 0 */
static __inline__  void __ALWAYS_INLINE bsem_init(bsem *bsem_p, int value) {
	if (value < 0 || value > 1) {
		fprintf(stderr, "bsem_init(): Binary semaphore can take only values 1 or 0");
		exit(1);
	}
	pthread_mutex_init(&(bsem_p->mutex), NULL);
	pthread_cond_init(&(bsem_p->cond), NULL);
	bsem_p->v = value;
}


/* Reset semaphore to 0 */
static __inline__  void __ALWAYS_INLINE bsem_reset(bsem *bsem_p) {
	bsem_init(bsem_p, 0);
}


/* Post to at least one thread */
static __inline__  void __ALWAYS_INLINE bsem_post(bsem *bsem_p) {
	pthread_mutex_lock(&bsem_p->mutex);
	bsem_p->v = 1;
	pthread_cond_signal(&bsem_p->cond);
	pthread_mutex_unlock(&bsem_p->mutex);
}


/* Post to all threads */
static __inline__  void __ALWAYS_INLINE bsem_post_all(bsem *bsem_p) {
	pthread_mutex_lock(&bsem_p->mutex);
	bsem_p->v = 1;
	pthread_cond_broadcast(&bsem_p->cond);
	pthread_mutex_unlock(&bsem_p->mutex);
}


/* Wait on semaphore until semaphore has value 0 */
static __inline__  void __ALWAYS_INLINE bsem_wait(bsem* bsem_p) {
	pthread_mutex_lock(&bsem_p->mutex);
	while (bsem_p->v != 1) {
		pthread_cond_wait(&bsem_p->cond, &bsem_p->mutex);
	}
	bsem_p->v = 0;
	pthread_mutex_unlock(&bsem_p->mutex);
}


/* ============================ JOB QUEUE =========================== */

/* Add (allocated) job to queue
 *
 * Notice: Caller MUST hold a mutex
 */
static __inline__  void __ALWAYS_INLINE
jobqueue_push(jobqueue_t * const restrict jobqueue_p, job_t * const restrict newjob)
{
	newjob->prev = NULL;
	switch(jobqueue_p->len){

		case 0:  /* if no jobs in queue */
					jobqueue_p->front = newjob;
					jobqueue_p->rear  = newjob;
					break;

		default: /* if jobs in queue */
					jobqueue_p->rear->prev = newjob;
					jobqueue_p->rear = newjob;
	}
	jobqueue_p->len++;
	bsem_post(jobqueue_p->has_items);
}

/*
 * Notice: No NEED for Caller to hold a mutex
 */
static __inline__  void __ALWAYS_INLINE
jobqueue_push_batch(jobqueue_t * const restrict jobqueue_p, job_t * const restrict newjobs, const size_t N)
{
	job_t * ljob         = &newjobs[N-1];
	job_t * const newjob = &newjobs[N-1];
// 	while ( (uintptr_t) ljob > (uintptr_t) newjobs) {
// 		ljob->prev = --ljob;
// 	}
// 	ljob->prev = NULL;
	for (int i=N-1; i>0; i--) newjobs[i].prev = &newjobs[i-1];
	newjobs[0].prev = NULL;
	
	pthread_mutex_lock(&jobqueue_p->rwmutex);
	switch(jobqueue_p->len){

		case 0:  /* if no jobs in queue */
					jobqueue_p->front = newjob;
					jobqueue_p->rear  = &newjobs[0];
					break;

		default: /* if jobs in queue */
					jobqueue_p->rear->prev = newjob;
					jobqueue_p->rear = &newjobs[0];
	}
	jobqueue_p->len += (int) N;
	
	bsem_post(jobqueue_p->has_items);
	pthread_mutex_unlock(&jobqueue_p->rwmutex);
}

/* Get first job from queue(removes it from queue)
 * 
 * Notice: Caller MUST hold a mutex
 */
static __inline__ job_t* __ALWAYS_INLINE
jobqueue_pull(jobqueue_t * const restrict jobqueue_p)
{
	job_t* job_p = jobqueue_p->front;

	switch(jobqueue_p->len){
		
		case 0:  /* if no jobs in queue */
					return NULL;
		
		case 1:  /* if one job in queue */
					jobqueue_p->front = NULL;
					jobqueue_p->rear  = NULL;
					break;
		
		default: /* if >1 jobs in queue */
					jobqueue_p->front = job_p->prev;
					
	}
	jobqueue_p->len--;
	
	/* Make sure has_items has right value */
	if (jobqueue_p->len > 0) {
		bsem_post(jobqueue_p->has_items);
	}

	return job_p;
}

/* Clear the queue */
static __inline__  void __ALWAYS_INLINE
jobqueue_clear(jobqueue_t * const restrict jobqueue_p)
{
// 	while(jobqueue_p->len){
// 		free(jobqueue_pull(jobqueue_p));
// 	}

	jobqueue_p->front = NULL;
	jobqueue_p->rear  = NULL;
	bsem_reset(jobqueue_p->has_items);
	jobqueue_p->len = 0;
}

/* Initialize queue */
static __inline__  int __ALWAYS_INLINE
jobqueue_init(jobqueue_t * const restrict jobqueue_p)
{
	pthread_mutex_init(&(jobqueue_p->rwmutex), NULL);
	if (jobqueue_p == NULL){
		return -1;
	}
	
	jobqueue_p->has_items = (struct bsem*)malloc(sizeof(struct bsem));
	if (jobqueue_p->has_items == NULL){
		return -1;
	}
	bsem_init(jobqueue_p->has_items, 0);
	jobqueue_p->len = 0;
	jobqueue_clear(jobqueue_p);
	return 0;
}


/* Free all queue resources back to the system */
static __inline__  void __ALWAYS_INLINE jobqueue_destroy(jobqueue_t * const restrict jobqueue_p){
	jobqueue_clear(jobqueue_p);
	free(jobqueue_p->has_items);
}

/* ============================ THREAD ============================== */


/* Sets the calling thread on hold */
void thread_hold (threadpool_t * const thpool_p); 

#endif /* _THREAD_POOL_H */
