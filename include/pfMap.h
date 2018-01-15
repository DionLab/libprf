#ifndef _PFMAP_H
#define _PFMAP_H
#include "pfConfig.h"
#include <smmintrin.h>

/*
 ************************************************************************************************
 *                                    DEFINITIONS                                               *
 ************************************************************************************************
 */

/*
  * INTEGER FORMAT  -------------------------------------
  *
  * |xxxxxxxxxxxxxxxxxxxxxxxxxxxxx|xx|x
  *                              |  | |
  *                              |  | + Mismatch or true match
  *                              |  +-- State
  *                              +----- Score
  */
#define SCORE_SHIFT               3
#define STATE_SHIFT               1
#define MATCH_SHIFT               0

/* 
 * CONSTANTS USED IN ALIGNMENT
 */
// Division by 8 is required to prevent tagging to underflow
// #  define NLOW             (((unsigned int) -536870912) >> SCORE_SHIFT )
#ifndef USE_SHORT_INT
#  define NLOW            -536870912/32
#  define STORED_INT_MIN  INT_MIN
#  define STORED_INT_MAX	INT_MAX
#  define SCORE_MASK                0xFFFFFFF8
#  define STATE_MASK                0x00000006
#  define MATCH_MASK                0x00000001
#  define LEFT_NEGATIVE_SCORE_MASK  0xE0000000
#  define CLEAR_MASK                0xFFFFFFF8
#else
#  define NLOW            (-((1<<12)>>3))
// #  define STORED_INT_MIN  INT_MIN
// #  define STORED_INT_MAX	INT_MAX
#  define SCORE_MASK                0xFFF8
#  define STATE_MASK                0x0006
#  define MATCH_MASK                0x0001
#  define LEFT_NEGATIVE_SCORE_MASK  0xE000
#  define CLEAR_MASK                0xFFF8
#endif

#define _I		 0
#define _IM		 0 
#define _II		 0 
#define _ID		 NLOW
#define _IX		 0

#define _D		-1
#define _DM		 0
#define _DI		 NLOW
#define _DD		 0
#define _DX		 0

#define _m		-1
#define _M		 3
#define _MM		 1
#define _MI		-15
#define _MD		-3
#define _MX		 0 

#define _BM		 0
#define _BD		-2


# define TO_SCORE(x) ((int) ((x < 0) ? (((unsigned int) x)>>SCORE_SHIFT) | LEFT_NEGATIVE_SCORE_MASK : ((unsigned int) x)>>SCORE_SHIFT))

//---------------------------------------------------------------
// ENUMERATION DEFINITIONS
//---------------------------------------------------------------
// This is the position within the SSE register tuple
// Pay attention that Storing functions and packing need to be aligned on that.
enum VectorPosition {
  /* Positions within both 4-tuple and array of 4-tuple */
  MATCH=2,
  INSERTION=3,
  DELETION=0,
  EXTRA=1
};

// This is the mask used to define the priority when performing comparisons
enum StatePriority {
  PRIORITY_MATCH     = 2,
  PRIORITY_INSERTION = 0,
  PRIORITY_DELETION  = 1,
  PRIORITY_EXTRA     = 3
};

#ifdef PRF_CORE_HEURISTIC
/********************* Heuristic structures *********************/
typedef union TransposeMatrix { const int * i; float * f;} TransposeMatrix;
#endif


typedef union {
	int elem[4];
	__m128i xmm;
} __my128i;

struct Map {
	__my128i RowInit;
	__my128i MatchScores;
	__my128i MismatchScores;
	__my128i SequenceNMatchScores;
	__my128i InsertionScores;
	__my128i DeletionScores;
	const char * Alphabet;
	const unsigned char * Alphabet_Mapping;
	int MatchScore;
	int MismatchScore;
	int NScore;
};

typedef struct Zone {
	int Score;
	unsigned int Begin;
	unsigned int End;
} Zone_t;

/*
 ************************************************************************************************
 *                                        VARIABLES                                             *
 ************************************************************************************************
 */

/*
 ************************************************************************************************
 *                                 FUNCTION DECLARATIONS                                        *
 ************************************************************************************************
 */
const struct Map * getDefaultMap();
int dumpDefaultMap(const char * FileName);
int loadMap(const char * const restrict, struct Map * const restrict);
void GetMapping(const struct Map * const restrict map, void * const restrict WorkSpace,
                const char * const restrict Genome, const char * const restrict Tag,
                Zone_t * const restrict zone, const size_t ReadLength, const size_t GenomeLength);

#endif /* _PFMAP_H */
