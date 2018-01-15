#ifndef PB_BAM_H_
#define PB_BAM_H_
#include <stdio.h>
#include "pb_common.h"
#ifndef __cplusplus
typedef void PacBioBAM_t;
#else
struct PacBioBAM_t {
	std::string subreads;
	std::string scraps;
};
extern "C" {
#endif

/* Common PacBio File */
int isPacBioBAM(const char * const restrict File);
PacBioBAM_t* OpenPacBioBAM(const char * const restrict File);
void ClosePacBioBAM(const PacBioBAM_t * const restrict PBBAM);

int getBAMSummary(const PacBioBAM_t * const restrict PBBAM, const unsigned int HQThreshold,
                  ZMWSummary_t * const restrict summary, FILE* restrict out);
#ifdef __cplusplus
}
#endif
#endif /* PB_BAM_H_ */
