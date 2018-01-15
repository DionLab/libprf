#ifndef _DISPATCH_H
#define _DISPATCH_H
#include "pfConfig.h"
#include "pfProfile.h"
#include "pfDispatchExt.h"
#include "pfOutput.h"


int dispatchFASTAFile(const struct Profile * restrict prf,
                      const FASTAStructure * const restrict FASTA,
#ifdef PRF_CORE_PCRE
                      struct RegEx * const restrict regex,
#endif
                      const OutputType_t * const restrict OutputType,
                      const size_t nCPUs);

int dispatchStreamFASTA(const struct Profile * const restrict prf,
                        FILE * const restrict Stream,
#ifdef PRF_CORE_PCRE
                        struct RegEx * const restrict regex,
#endif
                        const OutputType_t * const restrict OutputType,
                        const size_t nCPUs);

#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
// #include "pb_common.h"
typedef struct COMMON_PB {
	const struct Profile * profile;
	pthread_mutex_t PrintLock;
	const Compute_t * restrict Compute;           /* pointer to core functions */
	size_t * restrict Histograms;
	size_t * restrict Counters;
	const OutputType_t * restrict OutputType;
#ifdef PRF_CORE_PCRE
  struct RegEx * restrict regex;
#endif
} pb_common_t;
#endif

#ifdef PRF_INPUT_HDF5
int dispatchPacBio(const struct Profile * const restrict prf,
                   const PacBio_t * const restrict PB,
#ifdef PRF_CORE_PCRE
                   struct RegEx * const restrict regex,
#endif
                   const OutputType_t * const restrict OutputType,
                   const PacBioDispatchOptions_t * Options,
                   const size_t nCPUs);

//int alignHoleReads(const PacBio_t * const restrict PB, const struct Profile * const restrict prf,
//                   const Compute_t * const restrict CoreCompute,
//                   const OutputMethod Output,
//                   const void * const restrict OutputMethodOptions,
//                   const char * const restrict OutputFileName,
//                   const enum PacBioDispatchOptions_t Options,
//                   const unsigned int HoleNumber);
#endif

#ifdef PRF_INPUT_PBBAM
int dispatchPacBioBAM(const struct Profile * const restrict prf,
                      const PacBioBAM_t * const restrict PBBAM,
#ifdef PRF_CORE_PCRE
                      struct RegEx * const restrict regex,
#endif
                      const OutputType_t * const restrict OutputType,
                      const PacBioDispatchOptions_t * Options,
                      const size_t nCPUs);
#endif

#ifdef PRF_INPUT_RANDOM
int dispatchRandomData(const struct Profile * const restrict prf,
                       const RandomData_t * const restrict RD,
#ifdef PRF_CORE_PCRE
                       struct RegEx * const restrict regex,
#endif
                       const size_t nCPUs);
#endif

#endif
