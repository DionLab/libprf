#ifndef _COMPUTE_EXT_H
#define _COMPUTE_EXT_H 
#include "pfOutput.h"
#include "pfCompute.h"

#ifdef PRF_INPUT_HDF5
int alignHoleReads(const PacBio_t * const restrict PB, const struct Profile * const restrict prf,
									 const Compute_t * const restrict CoreCompute,
									 const OutputMethod Output,
									 const char * const restrict OutputFileName,
									 void * const restrict OutputMethodOptions,
									 const PacBioDispatchOptions_t * const restrict Options);
#endif

#endif /* _COMPUTE_EXT_H */
