/*******************************************************
                        PFTOOLS
 *******************************************************
  Apr 5, 2016 pfCompute.h
 *******************************************************
 (C) 2011-2016 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#ifndef _COMPUTE_H
#define _COMPUTE_H
#include <stdbool.h>
#include <pthread.h>
#include "pfConfig.h"
#include "pfProfile.h"

/*
 ************************************************************************************************
 *                                    DEFINITIONS                                               *
 ************************************************************************************************
 */

union lScores {
  /* SSE 4.1 can work on integer */
  int Element[4];
  __m128i xmm;
  /* Other have to rely upon float */
  float Elementf[4];
  __m128 xmmf;
};

union URegion {
  struct {
    int Begin;
    int End;
  };
  __m64 mm;
};

typedef struct Alignment {
  union {
    struct {
      union URegion column;
      union URegion row;
    };
    __m128i xmm;
		union URegion array[2];
  } Matrix __attribute__((aligned(16)));
	union {
    struct {
      union URegion profile;
      union URegion sequence;
    };
    __m128i xmm;
		union URegion array[2];
  } Region __attribute__((aligned(16)));
	
   /* The following order has to be kept due to sse zeroing, 128 bit alignment is required too*/
//    int ProtectedRegionBegin;  // JAL1
//    int ProtectedRegionEnd;    // JAL2;
//    int SequenceBegin;         // JALB;
//    int SequenceEnd;           // JALE;
   /* ---- */
   int Score;                    // JALS
   int IPMB;			 // Filled during xalit, used in coverage as length of found sequence
   int IPME;			 // Filled during xalit
   short int Orientation;
	 unsigned short int Cycles;
} Alignment_t;

typedef struct ScoreCycle {
	int Score;
	unsigned int Cycles;
	unsigned int Begin;
	unsigned int End;
} ScoreCycle_t;

/* This is the main structure holding the information to call computing procedures and exposes 
 * the function to output functions
 */
// pointer to the matrix building function
typedef void (*CoreComputeFct)(const struct Profile * const restrict prf,
                               const unsigned char * const restrict Sequence,
                               union lScores * const restrict matrix,
                               int * const restrict WORK,
                               unsigned short int * const restrict NumberOfInsertions,
                               const size_t BSEQ,
                               const size_t LSEQ);

// pointer to the sequence alignement search
typedef size_t (*GetBestAlignmentFct)(const union lScores * const restrict matrix,
                                      Alignment_t * const restrict Alignment,
                                      const size_t SeqLength,
                                      const size_t prfLength);
typedef size_t (*GetNextBestAlignmentFct)(const union lScores * const restrict matrix,
                                          const _Bool * const restrict HasAPath,
                                          Alignment_t * const restrict Alignment,
                                          const size_t SeqLength,
                                          const size_t prfLength);
// pointer to the sequence best score
typedef int (*GetBestScoreFct)(const union lScores * const restrict matrix,
                               const size_t SeqLength,
                               const size_t prfLength);

// pointer to function that return text sequence 
typedef int (*GetStateSequenceFct)(const union lScores * const restrict matrix,
                                   unsigned char * const restrict Alignment_Sequence,
                                   const union URegion * const restrict Alignment,
                                   const size_t matrix_ld,
                                   const size_t SeqLength);

typedef void (*GetAlignmentSequenceFct)(const union lScores * const restrict matrix,
                                        const unsigned char * const restrict Sequence,
                                        unsigned char * const restrict AlignmentSequence,
                                        const Alignment_t * const restrict Alignment,
                                        const size_t SeqLength,
                                        const size_t prfLength);
typedef int (*GetAlignmentsFct)(union lScores * const restrict matrix,
                                const struct Profile * const prf,
                                Alignment_t ** Alignments,
                                const size_t SeqLength);

// Global structure that holds it all
struct Compute {
	/* pointer to the matrix building function */
	const CoreComputeFct BuildMatrix;
	/* pointer to the sequence best score */
	const GetBestScoreFct GetBestScore;
	/* pointer to best alignement within matrix */
	const GetBestAlignmentFct GetBestAlignment;
	/* pointer to loop on best alignments score decreasing order */
	const GetNextBestAlignmentFct GetNextBestAlignment;
	/* pointer to the aligned state sequene */
	const GetStateSequenceFct GetStateSequence;
	/* pointer to get aligned sequence*/
	const GetAlignmentSequenceFct GetAlignmentSequence;
	/* pointer to get all alignments */
	const GetAlignmentsFct GetAlignments;
	/* Score data resolution */
	const unsigned int ScoreMask;
	const unsigned int LeftNegativeScoreMask;
	const unsigned char ScoreShift;
	/* State data resolution */
	const unsigned int StateMask;
	const unsigned char StateShift; 
	/* Matrix informations */
	const size_t MatrixColumnMultiplier;
	const size_t MatrixExtraColumn;
	const size_t MatrixExtraRow;
	/* Work array dimension */
	const size_t WorkCellSize;
};

typedef struct Compute Compute_t;


typedef struct RamdomData_s {
	unsigned int seed;
	unsigned int Length;
	unsigned int nCAGSuffix;
	size_t N;
} RandomData_t;

/*
 ************************************************************************************************
 *                                  VARIABLES REQUIRED                                          *
 ************************************************************************************************
 */
#ifdef PRF_CORE_STD
extern const Compute_t Standard_sse41;
#endif
#ifdef PRF_CORE_REPEAT
extern const Compute_t Repeat_sse41;
#endif
#ifdef PRF_CORE_FPGA
extern const Compute_t FPGA_sse41;
#endif
#ifdef PRF_CORE_EXT_PROFILE
extern const Compute_t ExtendedProfile_sse41;
#endif

#ifdef PRF_USE_AFFINITY
extern pthread_attr_t * restrict threads_attr;
#endif

/*
 ************************************************************************************************
 *                             SEQUENCE FUNCTION DECLARATIONS                                   *
 ************************************************************************************************
 */
#ifdef PRF_CORE_REPEAT
int GetBestScoreAndCycles(const union lScores * const matrix, const size_t SeqLength,
                          const size_t prfLength, ScoreCycle_t * const restrict SC);
#endif
/*
 ************************************************************************************************
 *                                   INLINE FUNCTIONS                                           *
 ************************************************************************************************
 */

static inline int __ALWAYS_INLINE GetElementScore (const union lScores * const restrict MatrixCell,
                                                   const Compute_t * const restrict compute,
                                                   const enum VectorPosition pos) {
	const int x = MatrixCell->Element[pos]; 
	return (int) (( x < 0) ? (((unsigned int) x) >> compute->ScoreShift) | compute->LeftNegativeScoreMask :\
                              ((unsigned int) x) >> compute->ScoreShift);
}

static inline int __ALWAYS_INLINE GetScore (const int MatrixCell,
                                            const Compute_t * const restrict compute) {
	const int x = MatrixCell; 
	return (int) (( x < 0) ? (((unsigned int) x) >> compute->ScoreShift) | compute->LeftNegativeScoreMask :\
	                          ((unsigned int) x) >> compute->ScoreShift);
}

static inline size_t __ALWAYS_INLINE GetMatrixLD(const struct Profile * const restrict prf,
                                                 const Compute_t * const restrict compute)
{
	return (prf->Length + compute->MatrixExtraColumn)*compute->MatrixColumnMultiplier;
}

static inline size_t __ALWAYS_INLINE GetWorkArraySize(const struct Profile * const restrict prf,
                                                      const Compute_t * const restrict compute)
{
	return (2UL*compute->WorkCellSize*(prf->Length)) + 63UL;
}

#endif /* _COMPUTE_H */
