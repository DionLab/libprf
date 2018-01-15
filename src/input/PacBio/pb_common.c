#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "pb_common.h"

int restrictReadstoHQ(const HoleData_t * const restrict HReg)
{
	/* Get the HQRegion from HReg starting for the end to speed things up */
	const unsigned int count = HReg->nRegions;
	HoleRegion_t * restrict RegPtr = &(HReg->Regions[count]);
	while ( (uintptr_t) --RegPtr >= (uintptr_t) &(HReg->Regions[0])) if (RegPtr->type == HQRegion) break;
	if (RegPtr->type != HQRegion) {
		fprintf(stderr,"Hole %u has no HQ region!!!\n", HReg->HoleNumber);
		return -1;
	}

	const int HQStart = RegPtr->start;
	const int HQStop  = RegPtr->stop;
	
	HoleRegion_t * restrict Regions = HReg->Regions;
	for (unsigned int r=0; r<count; r++) {
		if (Regions[r].type == Insert) {
			if (Regions[r].start < HQStart) Regions[r].start = HQStart;
			if (Regions[r].stop > HQStop) Regions[r].stop = HQStop;
			
			if (Regions[r].start >= Regions[r].stop) Regions[r].type = Useless;
		}
	}
	return HQStop - HQStart;
}
