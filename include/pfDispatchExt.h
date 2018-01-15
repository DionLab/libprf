#ifndef _DISPATCH_EXT_H
#define _DISPATCH_EXT_H
#include "pfConfig.h"
#include "pfInput.h"
#include "threadpool.h"

#if defined(__USE_WINAPI__)
#define pthread_exit(x) return (LPDWORD) x
#define THREADPOOL_FUNCTION(Fct) LPDWORD WINAPI Fct(_In_ LPVOID _Data)
#else
#define THREADPOOL_FUNCTION(Fct) void* Fct(threadarg_t * const restrict _Data)
#endif

#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
extern int PacBioHQThreshold;
enum PacBioSelection {
	PB_NONE                =    0,
	PB_DISPATCH_SUBREADS   =    1,
	PB_DISPATCH_ZMW        =    2,
	PB_DISPATCH_HQREGION   =    4,
	PB_DISPATCH_CONSENSUS  =    8,
	PB_DISPATCH_INVALID    =   16,
	PB_KEEP_INVALID        =   32,
	PB_DISCARD_FILTER      =   64,
	PB_DISCARD_HQREGION    =  128,
	PB_BEST_SUBREAD_OF_ZMW =  256,
	PB_NOT_ONLY_SEQUENCING =  512,
	PB_HAS_FILTER_SCORE    = 1024,
	PB_TEST_OUTPUT         = 2048
};

typedef struct PacBioDispatchOptions {
	unsigned int * ZMW; 
	unsigned int nZMW;
	unsigned int minReadAccuracy;
	unsigned int maxReadAccuracy;
	enum PacBioSelection Selection;
} PacBioDispatchOptions_t;

#define PB_DEFAULT_OPTIONS {\
	.ZMW = NULL,\
	.nZMW = 0U,\
	.minReadAccuracy = 0U,\
	.maxReadAccuracy = 1000U,\
	.Selection = PB_NONE\
}
#endif

/* ---------------------------------- PACBIO ---------------------------------- */
#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
#include "pb_common.h"
// typedef struct COMMON_PB {
// 	const struct Profile * profile;
// 	pthread_mutex_t PrintLock;
// 	const Compute_t * restrict Compute;           /* pointer to core functions */
// 	size_t * restrict Histograms;
// 	size_t * restrict Counters;
// 	const OutputType_t * restrict OutputType;
// #ifdef PRF_CORE_PCRE
//   struct RegEx * restrict regex;
// #endif
// } pb_common_t;

/* WARNING: This must be in agreement with job_t defined in threadpool.h for inheritance to function */
struct pb_job_s {
	struct pb_job_s * prev;          /* pointer to previous job   */
	HoleReads_t Reads;
	HoleData_t HReg;
}; 
typedef struct pb_job_s pb_job_t;


int dispatchPacBioBAMExt(const PacBioBAM_t * const restrict PBBAM,
												 threadpool_t * const restrict thpool,
                         const PacBioDispatchOptions_t * Options);

int dispatchPacBioExt(const PacBio_t * const restrict PB,
											threadpool_t * const restrict thpool,
                      const PacBioDispatchOptions_t * Options);

extern THREADPOOL_FUNCTION( (* * tp_pbsp) );
extern THREADPOOL_FUNCTION( (* * tp_pbzp) );
THREADPOOL_FUNCTION(tp_pb_filtertest);
#endif

#ifdef PRF_INPUT_FASTA
#include "pfSequence.h"
// typedef struct COMMON_FA {
// 	const struct Profile * profile;
// 	const OutputType_t * OutputType;
// 	const Compute_t * Compute;           /* pointer to core functions */
// 	pthread_mutex_t PrintLock;
// 	size_t * Histograms;
// 	size_t * Counters;
// #ifdef PRF_CORE_PCRE
//   const struct RegEx * regex;
// #endif
// } fa_common_t;

/* WARNING: This must be in agreement with job_t defined in threadpool.h for inheritance to function */
struct fasta_job_s {
	struct fasta_job_s * prev;          /* pointer to previous job   */
	Sequence_t Fasta_Sequence;
};

typedef struct fasta_job_s fasta_job_t;
int dispatchFASTAFileExt(const struct Profile * restrict prf,
                         const FASTAStructure * const restrict FASTA,
												 const Compute_t * const Model,
#ifdef PRF_CORE_PCRE
                         struct RegEx * const restrict regex,
#endif
                         const OutputType_t * const restrict OutputType,
                         const size_t nCPUs);
int dispatchStreamFASTAExt(FILE * const restrict Stream, threadpool_t * const restrict thpool);

extern THREADPOOL_FUNCTION( (* * tp_fap) );

#endif

#endif /* _DISPATCH_EXT_H */
