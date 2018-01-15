#ifndef _PB_HDF5_H
#define _PB_HDF5_H
/*
 * Bases functions to address PacBio files
 */
#include "pfConfig.h"
#include <hdf5.h>
#include <string.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <xmmintrin.h>

#include "pb_common.h"

#ifdef _TEST_
# include "matrix_functions.h"
# include "pfProfile.h"
#endif

#define FILENAME_SIZE 256
#define DIRECTORY_SIZE 256
#define JOBS_BATCH_SIZE 1024

enum PBContent { BASECALLING=1, CONSENSUS_BASECALLING=2, PULSE=4, TRACE=8 }; 
enum HoleStatus { SEQUENCING=0,ANTIHOLE=1,FIDUCIAL=2,SUSPECT=3,
                  ANTIMIRROR=4,FDZMW=5,FBZMW=6,ANTIBEAMLET=7,
                  OUTSIDEFOV=8 };

typedef struct RegionIndex {
	unsigned int start;
	unsigned int stop;
	hsize_t BaseCallsOffset;
	unsigned int BaseCallsLength;
} RegionIndex_t;

typedef struct PacBioBas {
	RegionIndex_t * Regions;
} PacBioBas_t;

typedef struct PacBioTrc {
	float (*DecodeTable)[256];
	float *Bias;
	size_t Length;
} PacBioTrc_t;

typedef struct PacBio_s {
	char Directory[DIRECTORY_SIZE];
	char BaseFileName[FILENAME_SIZE];
	char * restrict * PartsFileName;
	hid_t Base;
	hid_t * restrict Parts;
	unsigned char * restrict HoleFileIndex;
	unsigned int * restrict HoleFilePositionIndex;
	unsigned char * Status;
	short int (*Coordinates)[2];
	unsigned int nParts;
	unsigned int nHoles;
	enum PBContent Content;
	PacBioBas_t BaseCalling;
	PacBioBas_t ConsensusBaseCalling;
	// PacBioPls_t Pulses;
	PacBioTrc_t Traces;
} PacBio_t;

typedef struct Read {
	char * start;
	char * stop;
} Read_t;

typedef HoleReads_t HoleQVs_t;

typedef struct BaseCallsPulses {
	unsigned int * restrict FramesStart;
	unsigned int * restrict FramesEnd;
	char * restrict Bases; /* this is aligned on a 16 byte boundary */
	size_t nBases;
} BaseCallsPulses_t;

typedef struct Traces_s {
	unsigned char * restrict T;
	unsigned char * restrict G;
	unsigned char * restrict A;
	unsigned char * restrict C;
	unsigned char * memory;
	unsigned int stride;
} Traces_t;

extern const char *const HoleStatusStrings[];

/* Common PacBio File */
int isPacBioH5(const char * const restrict File);
PacBio_t * PacBioOpen(const char * const restrict  BaseFileName);
int PacBioClose(PacBio_t * PB);
int IndexRegions(PacBio_t * const restrict PB);
int PopulateHoleStatus(PacBio_t * const restrict PB);
int PopulateHoleCoordinates(PacBio_t * const restrict PB);
int getSequencingHoles(const PacBio_t * const PB, unsigned int ** SequencingHoles);
int getHDFSummary(PacBio_t * const restrict PB, const unsigned int HQThreshold,
                  ZMWSummary_t * const restrict summary, FILE* restrict out);


/* Base Calling PacBio File */
int getHoleRegions(const PacBio_t * const restrict PB, HoleData_t * HReg, const unsigned int HoleNumber);
int getHoleReads_HN(const PacBio_t * const restrict PB, HoleReads_t *HRead, const unsigned int HoleNumber);
int getHoleQVs_HN(const PacBio_t * const restrict PB, HoleQVs_t *HRead, const unsigned int HoleNumber);
int getHoleReads_HR(const PacBio_t * const restrict PB, HoleReads_t *HRead, const HoleData_t * const restrict HReg);
int getHoleQVs_HR(const PacBio_t * const restrict PB, HoleQVs_t *HRead, const HoleData_t * const restrict HReg);
int getBasecallsPulses(const PacBio_t * const restrict PB, BaseCallsPulses_t * const restrict BP,
                       const unsigned int HoleNumber);


/* Consensus Base Calling PacBio File */
int IndexConsensus(PacBio_t * const restrict PB);

/* Trace PacBio File */
int PopulateDecodeTable(PacBio_t * const restrict PB);
Traces_t * getTraceIndex(const PacBio_t * const restrict PB, const unsigned int HoleNumber);
__m128* DecodeTracePositionIndex(const PacBio_t * const restrict PB, Traces_t * const restrict Indices,
                                 const unsigned int HoleNumber);
__m128* getTrace(const PacBio_t * const restrict PB, const unsigned int HoleNumber); 


/* Inline Functions */

extern inline void __ALWAYS_INLINE freeHoleReads(HoleReads_t* HRead)
{
	if (HRead->Stream) free(HRead->Stream);
	HRead->Size = 0;
}

extern inline void __ALWAYS_INLINE freeHoleRegions(HoleData_t * const HReg) {
	if (HReg->Regions) free(HReg->Regions);
	HReg->nRegions = 0U;
	HReg->nAllocatedRegions = 0U;
}

static inline void __ALWAYS_INLINE freeBaseCallsPulses(BaseCallsPulses_t * const BP) {
	if (BP->FramesStart) free(BP->FramesStart);
	if (BP->FramesEnd) free(BP->FramesEnd);
	if (BP->Bases) _mm_free(BP->Bases);
	BP->nBases = 0;
}

static inline void __ALWAYS_INLINE freeTraces(Traces_t * TR) {
	if (TR->memory) _mm_free(TR->memory);
	free(TR);
}

#define freeHoleQVs freeHoleReads
#define HOLE_DATA_INIT {.HoleNumber = 0, .nRegions = 0, .Regions = NULL }
#define HOLE_READS_INIT { .Stream = NULL, .Size = 0UL }
#define HOLE_QVS_INIT HOLE_READS_INIT

#endif /*_PB_HDF5_H*/
