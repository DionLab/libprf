#ifndef PB_COMMON_H_
#define PB_COMMON_H_

#define SUCCESS 0
#define ERROR   -15

enum RegionType { Adapter=0, Insert=1, HQRegion=2, Useless=3 };

typedef struct HoleRegions {
		int type;
		int start;
		int stop;
		int quality;
} HoleRegion_t;

typedef struct HoleData {
	unsigned int HoleNumber;
	unsigned int nRegions;
	unsigned int nAllocatedRegions;
	short int Coordinates[2];
	HoleRegion_t * Regions;
} HoleData_t;

typedef struct HoleRead {
	unsigned char * Stream;
	size_t Size;	
} HoleReads_t;

typedef struct ZMWSummary {
	unsigned int nSequencing;
	unsigned int nHQaboveThreshold;
	unsigned int nBelowHQThreshold;
	unsigned int nInvalidHQRegion;
	unsigned int nWithinHQSubreads;
	unsigned int nOutofHQSubreads;
	unsigned int nBelowHQThresholdWithin;
	unsigned int nBelowHQThresholdOutof;
	unsigned int nInvalidHQRegionSubreads;
} ZMWSummary_t;

int restrictReadstoHQ(const HoleData_t * const restrict HReg);

#endif /* PB_COMMON_H_ */
