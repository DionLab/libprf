/*******************************************************
 *                        PFTOOLS
 *******************************************************
 *  Nov 18, 2014 coverage2_sse41.c
 *******************************************************
 * (C) 2012-14 Swiss Institute of Bioinformatics
 *     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h> 
#include <smmintrin.h>
#include "pfCompute.h"
#include "pfOutput.h"
#include "sse41_inline_fcts.h"


/* Priority */
enum StatePriority {
  PRIORITY_MATCH     = 1,
  PRIORITY_INSERTION = 2,
  PRIORITY_DELETION  = 0,
  PRIORITY_EXTRA     = 3
};

/*
  * INTEGER FORMAT  -------------------------------------
  *
  * |xxxxxxxxxxxxxxxxxxxxxxxxxxxxx|xx|x|
  *                              |  | |
	*                              |  | + Cycle
  *                              |  +-- State
  *                              +----- Score
  */
# define SCORE_SHIFT               3
# define SCORE_MASK                0xFFFFFFF8
# define STATE_SHIFT               1
# define STATE_MASK                0x00000006
# define CYCLE_SHIFT               0
# define CYCLE_MASK                0x00000001
# define LEFT_NEGATIVE_SCORE_MASK  0xE0000000
# define CLEAN_STATE               0xFFFFFFF8
// # define CLEAN_PHASE               0xFFFFFFFC
// # define CLEAN_STATE_AND_PHASE     0xFFFFFFF0


#define TO_SCORE(x) (int) ((x < 0) ? (((unsigned int) x)>>SCORE_SHIFT) | LEFT_NEGATIVE_SCORE_MASK : ((unsigned int) x)>>SCORE_SHIFT)
#define GET_STATE(x) ((x & STATE_MASK) >> STATE_SHIFT)
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MAX_INSERT_LENGTH 5000


void repeat_sse41(const struct Profile * const restrict prf, const unsigned char * const restrict Sequence,
                  union lScores * const restrict matrix, int * const restrict WORK,
                  unsigned short int * const restrict NumberOfInsertions, const size_t BSEQ, const size_t LSEQ)
/*
 * WARNING: for SSE version, WORK should be aligned on cache size (64b) and 2 times the 
 *            (profile size + 1)*sizeof(lscores) + 63 to align second part to cache line
 *          matrix should be of size at least (sequence size+1)*(profile size + 1)*sizeof(lscores) 
 *            and aligned on 16b.
 */
{
	int KOPD;
	const union lScores * restrict IOP_R;
	union lScores * restrict IOP_W = (union lScores *) WORK;

	const __m128i __MatchMask     = _mm_set1_epi32(PRIORITY_MATCH     << STATE_SHIFT);
	const __m128i __InsertionMask = _mm_set1_epi32(PRIORITY_INSERTION << STATE_SHIFT);
	const __m128i __DeletionMask  = _mm_set1_epi32(PRIORITY_DELETION  << STATE_SHIFT);
	const __m128i __ExtraMask     = _mm_set1_epi32(PRIORITY_EXTRA     << STATE_SHIFT);

	register const TransitionScores * const restrict Transitions = prf->Scores.Insertion.Transitions;
	const StoredIntegerFormat * const restrict Match             = prf->Scores.Match.Alphabet;
	const StoredIntegerFormat * const restrict Insertion         = prf->Scores.Insertion.Alphabet;
	const size_t AlignStep                                       = prf->Scores.Match.AlignStep;
	const size_t prfLength = prf->Length;
	const size_t matrixLD  = prf->Length+1; 
	
	const __m128i __NLOW = _mm_set1_epi32((unsigned int)NLOW << SCORE_SHIFT);

	/* NOTE: The following part could be replaced and performed only once for a profile as it
	 *       is profile dependent. Nevertheless it does a good job loading Match and Transition
	 *       matrices into the cache hierarchy.
	 */

	/*
	 * Initialize Insertion and Match Entrance Line using FirstSequenceProtein
	 */
	{
		register const StoredIntegerFormat * restrict lMatch = (const StoredIntegerFormat *) &Match[_D];
		register const ScoreTuple * restrict FirstSequenceProtein = prf->Scores.Insertion.FirstSequenceProtein;

		/*
		 * Set all insertion lengths to zero
		 */
// 		memset(NumberOfInsertions, 0, (1+prfLength)*sizeof(unsigned short int));

		/*
		 * PROFILE COLUMN 0 entrance
		 */
		IOP_W[0].Element[MATCH]     = (((int) FirstSequenceProtein[0].To[MATCH]) << SCORE_SHIFT)     | (PRIORITY_EXTRA << STATE_SHIFT);
		IOP_W[0].Element[INSERTION] = (((int) FirstSequenceProtein[0].To[INSERTION]) << SCORE_SHIFT) | (PRIORITY_EXTRA << STATE_SHIFT);
		IOP_W[0].Element[DELETION]  = (((int) FirstSequenceProtein[0].To[DELETION]) << SCORE_SHIFT)  | (PRIORITY_EXTRA << STATE_SHIFT);
		IOP_W[0].Element[EXTRA]     = (int) (((unsigned int) NLOW << SCORE_SHIFT) | (PRIORITY_EXTRA << STATE_SHIFT));
		KOPD                        = (int) FirstSequenceProtein[0].To[DELETION];

// 		{
// 			__m128i __FirstSequenceProtein = LoadStoredIntegerVector(&(FirstSequenceProtein[0]));
// 			__FirstSequenceProtein = _mm_slli_epi32(__FirstSequenceProtein, SCORE_SHIFT);
// 			__FirstSequenceProtein = _mm_or_si128(__FirstSequenceProtein, __ExtraMask);
// 		}

		register const TransitionScores (* restrict pTransitions) = &Transitions[1];
		register union lScores * restrict pIOP = &IOP_W[1];

		/*
		 * LOOP THROUGH THE REST OF THE PROFILE
		 */
		for (size_t iprf=1; iprf<=prfLength; ++iprf ) {
			register const int KD = KOPD + (int) *lMatch;
			lMatch += AlignStep;

			// Transform KD into a vector
			__m128i __KD = _mm_set1_epi32(KD);
			// Load Transitions and  Convert signed WORD into signed DWORD
			__m128i __TransitionsD = LoadStoredIntegerVector(&(pTransitions->From[DELETION]));

			// Add KD to Transitions
			__TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

			// Move to next profile transitions
			pTransitions++;

			// Load FirstSequenceProtein
			__m128i __FirstSequenceProtein = LoadStoredIntegerVector(&(FirstSequenceProtein[iprf]));

			// Paste index in lowest 2 bits
			__TransitionsD         = _mm_slli_epi32(__TransitionsD, SCORE_SHIFT);
			__FirstSequenceProtein = _mm_slli_epi32(__FirstSequenceProtein, SCORE_SHIFT);
			__TransitionsD         = _mm_or_si128(__TransitionsD, __DeletionMask);
			__FirstSequenceProtein = _mm_or_si128(__FirstSequenceProtein, __ExtraMask);

			// Get minimum ( this is SSE 4.1 )
			__m128i __max = _mm_max_epi32(__FirstSequenceProtein, __TransitionsD);

			// Store into matrix
//       _mm_store_si128(&matrix[iprf].xmm, __max);      
			_mm_store_si128( &(IOP_W[iprf].xmm), __max);
			
			// Clean extra bits
			__max = _mm_srai_epi32(__max, SCORE_SHIFT);

			// Store IOPI and IOPM and IOPD
// 			_mm_store_si128( &(IOP_W[iprf].xmm), __max);

			// Set KOPD ( this is SSE 4.1 )
			KOPD = _mm_extract_epi32(__max, DELETION);
		}
		
		/*
		 * HANDLES CIRCULAR POTENTIAL
		 */
		{
			const int itmp1 = (int) (((unsigned int) IOP_W[prfLength].Element[MATCH] )     | CYCLE_MASK);
			const int itmp2 = (int) (((unsigned int) IOP_W[prfLength].Element[INSERTION] ) | CYCLE_MASK);
			const int itmp3 = (int) (((unsigned int) IOP_W[prfLength].Element[DELETION] )  | CYCLE_MASK);
		
			IOP_W[0].Element[MATCH]     = MAX(itmp1, IOP_W[0].Element[MATCH]);
			IOP_W[0].Element[INSERTION] = MAX(itmp2, IOP_W[0].Element[INSERTION]);
			IOP_W[0].Element[DELETION]  = MAX(itmp3, IOP_W[0].Element[DELETION]);
		}
		
		_mm_store_si128(&matrix[0].xmm, IOP_W[0].xmm); 
		
		IOP_W[0].Element[MATCH]     = TO_SCORE(IOP_W[0].Element[MATCH]);
		IOP_W[0].Element[INSERTION] = TO_SCORE(IOP_W[0].Element[INSERTION]);
		KOPD                        = TO_SCORE(IOP_W[0].Element[DELETION]);
		
		lMatch = (const StoredIntegerFormat *) &Match[_D];
		for (size_t iprf=1; iprf<=prfLength; ++iprf) {
			const int KD = KOPD + (int) *lMatch;
			lMatch += AlignStep;

			__m128i __KD = _mm_set1_epi32(KD);
			__m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
			__TransitionsD = _mm_add_epi32(__TransitionsD, __KD);
			__TransitionsD = _mm_insert_epi32(__TransitionsD, NLOW, EXTRA);
			
			// Paste index in lowest 2 bits
			__TransitionsD = _mm_slli_epi32(__TransitionsD, SCORE_SHIFT);
			__TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);
			
			// Get minimum ( this is SSE 4.1 )
			__TransitionsD = _mm_max_epi32( IOP_W[iprf].xmm, __TransitionsD);
			
			// Store into matrix
			_mm_store_si128(&matrix[iprf].xmm, __TransitionsD); 
			
			// Clean extra bits
			__TransitionsD = _mm_srai_epi32(__TransitionsD, SCORE_SHIFT);
			
			// Store IOPI and IOPM and IOPD
			IOP_W[iprf].xmm = __TransitionsD;
			
			// Set KOPD ( this is SSE 4.1 )
			KOPD = _mm_extract_epi32(IOP_W[iprf].xmm, DELETION);
		}
	}

	// Swap and assign Read and write pointers
	IOP_R = IOP_W;
	IOP_W = (union lScores*) (((uintptr_t) &IOP_R[matrixLD] + 63) & ~63);

	/*
	 * LOOP THROUGH THE SEQUENCE STRING
	 */
	union lScores * restrict MatrixPtr = matrix + matrixLD;
	for ( size_t iseq=BSEQ; iseq < LSEQ-1; ++iseq) {
		register const size_t j1 = (size_t) Sequence[iseq];
 		int KOPM;
		register const StoredIntegerFormat * restrict lInsertion = Insertion;
		/*
		 * PROFILE COLUMN 0 entrance
		 */
		{
			register int KI = IOP_R[0].Element[INSERTION] + (int) lInsertion[j1];

			// Transform KI into a vector
			__m128i __KI = _mm_set1_epi32(KI);
			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[0].From[INSERTION]));
			// Add KI to Transition
			__TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[0].From[EXTRA]));

			// Paste index in lowest 2 bits
			__TransitionsI = _mm_slli_epi32(__TransitionsI, SCORE_SHIFT);
			__TransitionsX = _mm_slli_epi32(__TransitionsX, SCORE_SHIFT);
			__TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);
			__TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);

			// Get maximum( this is SSE 4.1 )
			__m128i __max = _mm_max_epi32(__TransitionsI, __TransitionsX);
			
			// Store IOPI and IOPM and IOPD
			IOP_W[0].xmm = __max;
			
			// Clean extra bits
			__max = _mm_srai_epi32(__max, SCORE_SHIFT);

			// Store KOPD
			KOPD = _mm_extract_epi32(__max, DELETION);
			
			// Store KOPM 
			KOPM = MAX(_mm_extract_epi32(__max, MATCH), IOP_R[0].Element[MATCH]);
		}

		lInsertion += AlignStep;
		register const StoredIntegerFormat * restrict lMatch = Match;

		/*
		 * LOOP THROUGH THE REST OF THE PROFILE
		 */
		for (size_t iprf=1; iprf<=prfLength; ++iprf ) {
			int KM = KOPM                           + (int) lMatch[j1];
			int KI = IOP_R[iprf].Element[INSERTION] + (int) lInsertion[j1];
			int KD = KOPD                           + (int) lMatch[_D];

			lMatch     += AlignStep;
			lInsertion += AlignStep;

			KOPM = IOP_R[iprf].Element[MATCH];

			// Transform KM into a vector
			__m128i __KM = _mm_set1_epi32(KM);
			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsM = LoadStoredIntegerVector(&(Transitions[iprf].From[MATCH]));
			// Add KM to Transition
			__TransitionsM = _mm_add_epi32(__TransitionsM, __KM);

			// Transform KI into a vector
			__m128i __KI = _mm_set1_epi32(KI);
			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[iprf].From[INSERTION]));
			// Add KI to Transition
			__TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

			// Paste index in lowest 2 bits
			__TransitionsM = _mm_slli_epi32(__TransitionsM, SCORE_SHIFT);
			__TransitionsI = _mm_slli_epi32(__TransitionsI, SCORE_SHIFT);
			__TransitionsM = _mm_or_si128(__TransitionsM, __MatchMask);
			__TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);

			// Get maximum ( this is SSE 4.1 )
			__m128i __max1 = _mm_max_epi32(__TransitionsM, __TransitionsI);

			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));

			// Transform KD into a vector
			__m128i __KD = _mm_set1_epi32(KD);
			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
			// Add KD to Transition
			__TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

			// Paste index in lowest 2 bits
			__TransitionsX = _mm_slli_epi32(__TransitionsX, SCORE_SHIFT);
			__TransitionsD = _mm_slli_epi32(__TransitionsD, SCORE_SHIFT);
			__TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);
			__TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);

			// Get maximum ( this is SSE 4.1 )
			__TransitionsD = _mm_max_epi32(__TransitionsD, __TransitionsX);
			__TransitionsD = _mm_max_epi32(__TransitionsD, __max1);
			

			// Store IOPI and IOPM and IOPD
			IOP_W[iprf].xmm = __TransitionsD;
			
			// Clean extra bits
			__TransitionsD = _mm_srai_epi32(__TransitionsD, SCORE_SHIFT);

			// Set KOPD ( this is SSE 4.1 )
			KOPD = _mm_extract_epi32(__TransitionsD, DELETION);
		}

		/*
		 * HANDLES CIRCULAR POTENTIAL
		 */
// 		if (iseq == 8) {
// 			printf("\nIOP_R\n%lu MATCH: ", iseq);
// 			for (int i=0; i<=prfLength; i++) {
// 				printf("%+5i  \t", (IOP_R[i].Element[MATCH]));  
// 			}
// 			printf("\n");
// 			
// 			printf("%lu DELET: ", iseq);
// 			for (int i=0; i<=prfLength; i++) {
// 				printf("%+5i  \t", (IOP_R[i].Element[DELETION]));  
// 			}
// 			printf("\n");
// 			
// 			printf("%lu INSER: ", iseq);
// 			for (int i=0; i<=prfLength; i++) {
// 				printf("%+5i  \t", (IOP_R[i].Element[INSERTION]));  
// 			}
// 			printf("\n");
// 			
// 			printf("\nIOP_W\n%lu MATCH: ", iseq);
// 			for (int i=0; i<=prfLength; i++) {
// 				printf("%+5i %i\t", TO_SCORE(IOP_W[i].Element[MATCH]), GET_STATE(IOP_W[i].Element[MATCH]));  
// 			}
// 			printf("\n");
// 			
// 			printf("%lu DELET: ", iseq);
// 			for (int i=0; i<=prfLength; i++) {
// 				printf("%+5i %i\t", TO_SCORE(IOP_W[i].Element[DELETION]), GET_STATE(IOP_W[i].Element[DELETION]));  
// 			}
// 			printf("\n");
// 			
// 			printf("%lu INSER: ", iseq);
// 			for (int i=0; i<=prfLength; i++) {
// 				printf("%+5i %i\t", TO_SCORE(IOP_W[i].Element[INSERTION]), GET_STATE(IOP_W[i].Element[INSERTION]));  
// 			}
// 			printf("\n");
// 		}
		{
			const int itmp1 = (int) (((unsigned int) IOP_W[prfLength].Element[MATCH])     | CYCLE_MASK);
			const int itmp2 = (int) (((unsigned int) IOP_W[prfLength].Element[INSERTION]) | CYCLE_MASK);
			const int itmp3 = (int) (((unsigned int) IOP_W[prfLength].Element[DELETION])  | CYCLE_MASK);
		
			IOP_W[0].Element[MATCH]     = MAX(itmp1, IOP_W[0].Element[MATCH]);
			IOP_W[0].Element[INSERTION] = MAX(itmp2, IOP_W[0].Element[INSERTION]);
			IOP_W[0].Element[DELETION]  = MAX(itmp3, IOP_W[0].Element[DELETION]);
		}
		assert(&MatrixPtr[0] >= matrix && &MatrixPtr[0] < &matrix[matrixLD*(LSEQ+1)]);
		_mm_store_si128(&MatrixPtr[0].xmm, IOP_W[0].xmm);
		
		IOP_W[0].Element[MATCH]     = TO_SCORE(IOP_W[0].Element[MATCH]);
		IOP_W[0].Element[INSERTION] = TO_SCORE(IOP_W[0].Element[INSERTION]);
		IOP_W[0].Element[DELETION]  = TO_SCORE(IOP_W[0].Element[DELETION]);
		KOPD = IOP_W[0].Element[DELETION];
		
		lMatch = (const StoredIntegerFormat *) &Match[_D];
		for (size_t iprf=1; iprf<=prfLength; ++iprf) {
			const int KD = KOPD + (int) *lMatch;
			lMatch += AlignStep;

			__m128i __KD = _mm_set1_epi32(KD);
			__m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
			__TransitionsD = _mm_add_epi32(__TransitionsD, __KD);
			__TransitionsD = _mm_insert_epi32(__TransitionsD, NLOW, EXTRA);
			
			// Paste index in lowest 2 bits
			__TransitionsD = _mm_slli_epi32(__TransitionsD, SCORE_SHIFT);
			__TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);
			
			// Get minimum ( this is SSE 4.1 )
			__TransitionsD = _mm_max_epi32( IOP_W[iprf].xmm, __TransitionsD);
			
			// Store into matrix
			assert(&MatrixPtr[iprf] >= matrix && &MatrixPtr[iprf] < &matrix[matrixLD*(LSEQ+1)]);
			_mm_store_si128(&MatrixPtr[iprf].xmm, __TransitionsD); 
			
			// Clean extra bits
			__TransitionsD = _mm_srai_epi32(__TransitionsD, SCORE_SHIFT);
			
			// Store IOPI and IOPM and IOPD
			IOP_W[iprf].xmm = __TransitionsD;
			
			// Set KOPD ( this is SSE 4.1 )
			KOPD = _mm_extract_epi32(__TransitionsD, DELETION);
		}
		
// 		if (iseq == 8) {
// 			printf("\nIOP_W CIRCULAR\n%lu MATCH: ", iseq);
// 			for (int i=0; i<=prfLength; i++) {
// 				printf("%+5i  \t", (IOP_W[i].Element[MATCH]));  
// 			}
// 			printf("\n");
// 			
// 			printf("%lu DELET: ", iseq);
// 			for (int i=0; i<=prfLength; i++) {
// 				printf("%+5i  \t", (IOP_W[i].Element[DELETION]));  
// 			}
// 			printf("\n");
// 			
// 			printf("%lu INSER: ", iseq);
// 			for (int i=0; i<=prfLength; i++) {
// 				printf("%+5i  \t", (IOP_W[i].Element[INSERTION]));  
// 			}
// 			printf("\n");
// 		}
		
		MatrixPtr += matrixLD;
		
		/*
		 * Swap Read and Write pointers
		 */
		{
			const register union lScores * const Iptr = IOP_W;
			IOP_W = (union lScores *) IOP_R;
			IOP_R = Iptr;
		}
	}

	/*
	 * Last position on the Sequence using LastSequenceProtein
	 */
	{
		register const StoredIntegerFormat * restrict lInsertion = Insertion;
		const int j1 = (int) Sequence[LSEQ-1];

		/*
		 * PROFILE COLUMN 0 entrance
		 */
		int KOPM     = IOP_R[0].Element[MATCH];
		int KI       = IOP_R[0].Element[INSERTION] + (int) lInsertion[j1];
		KI          += (int) Transitions[0].Element[_ID];

		const int tKOPD = MAX( KI, (int) Transitions[0].Element[_XD] );
		KOPD = MAX(KOPD, tKOPD);
		register const ScoreTuple * const restrict LastSequenceProtein = prf->Scores.Insertion.LastSequenceProtein;

		register const StoredIntegerFormat * restrict lMatch = Match;
		lInsertion += AlignStep;
		__m128i __Scores = __NLOW;

		/*
		 * LOOP THROUGH THE REST OF THE PROFILE
		 */
		__m128i __TransitionsD;
		for (size_t iprf=1; iprf<=prfLength; ++iprf) {
			const int KM = KOPM                           + lMatch[j1];
			KI           = IOP_R[iprf].Element[INSERTION] + lInsertion[j1];
			const int KD = KOPD                           + lMatch[_D];

			lMatch     += AlignStep;
			lInsertion += AlignStep;

			KOPM = IOP_R[iprf].Element[MATCH];

			/*----------------------------------------------------------------------------------------------------------*/
			// Transform KM into a vector
			__m128i __KM = _mm_set1_epi32(KM);
			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsM = LoadStoredIntegerVector(&(Transitions[iprf].From[MATCH]));
			// Insert LastProteinSequence
			__TransitionsM = _mm_insert_epi32(__TransitionsM, (int) LastSequenceProtein[iprf].From[MATCH], EXTRA);
			// Add KM to Transition
			__TransitionsM = _mm_add_epi32(__TransitionsM, __KM);

			// Transform KI into a vector
			__m128i __KI = _mm_set1_epi32(KI);
			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[iprf].From[INSERTION]));
			// Insert LastProteinSequence
			__TransitionsI = _mm_insert_epi32(__TransitionsI, (int) LastSequenceProtein[iprf].From[INSERTION], EXTRA);
			// Add KI to Transition
			__TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

			// Paste index in lowest 2 bits
			__TransitionsM = _mm_slli_epi32(__TransitionsM, SCORE_SHIFT);
			__TransitionsI = _mm_slli_epi32(__TransitionsI, SCORE_SHIFT);
			__TransitionsM = _mm_or_si128(__TransitionsM, __MatchMask);
			__TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);

			// Get maximum ( this is SSE 4.1 )
			__m128i __max1 = _mm_max_epi32(__TransitionsM, __TransitionsI);

			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));

			// Transform KD into a vector
			__m128i __KD = _mm_set1_epi32(KD);
			// Load Transitions and Convert signed WORD into signed DWORD
			__TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
			// Insert LastProteinSequence
			__TransitionsD = _mm_insert_epi32(__TransitionsD, (int) LastSequenceProtein[iprf].From[DELETION], EXTRA);
			// Add KD to Transition
			__TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

			// Paste index in lowest 2 bits
			__TransitionsX = _mm_slli_epi32(__TransitionsX, SCORE_SHIFT);
			__TransitionsD = _mm_slli_epi32(__TransitionsD, SCORE_SHIFT);
			__TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);
			__TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);

			// Get minimum ( this is SSE 4.1 )
			__TransitionsD = _mm_max_epi32(__TransitionsD, __TransitionsX);
			__TransitionsD = _mm_max_epi32(__TransitionsD, __max1);

			// Store into matrix
			IOP_W[iprf].xmm = __TransitionsD;
			
			// Clean extra bits
			__TransitionsD = _mm_srai_epi32(__TransitionsD, SCORE_SHIFT);

			// Set KOPD ( this is SSE 4.1 )
			KOPD = _mm_extract_epi32(__TransitionsD, DELETION);
		}
		
		/*
		 * HANDLES CIRCULAR POTENTIAL
		 */
		{
			const int itmp1 = (int) (((unsigned int) IOP_W[prfLength].Element[MATCH])     | CYCLE_MASK);
			const int itmp2 = (int) (((unsigned int) IOP_W[prfLength].Element[INSERTION]) | CYCLE_MASK);
			const int itmp3 = (int) (((unsigned int) IOP_W[prfLength].Element[DELETION])  | CYCLE_MASK);
		
			IOP_W[0].Element[MATCH]     = MAX(IOP_W[0].Element[MATCH],     itmp1);
			IOP_W[0].Element[INSERTION] = MAX(IOP_W[0].Element[INSERTION], itmp2);
			IOP_W[0].Element[DELETION]  = MAX(IOP_W[0].Element[DELETION],  itmp3);
		}
		assert(&MatrixPtr[0] >= matrix && &MatrixPtr[0] < &matrix[matrixLD*(LSEQ+1)]);
		_mm_store_si128(&MatrixPtr[0].xmm, IOP_W[0].xmm);
		
		IOP_W[0].Element[MATCH]     = TO_SCORE(IOP_W[0].Element[MATCH]);
		IOP_W[0].Element[INSERTION] = TO_SCORE(IOP_W[0].Element[INSERTION]);
		IOP_W[0].Element[DELETION]  = TO_SCORE(IOP_W[0].Element[DELETION]);
		KOPD = IOP_W[0].Element[DELETION];
		
		lMatch = (const StoredIntegerFormat *) &Match[_D];
		for (size_t iprf=1; iprf<=prfLength; ++iprf) {
			const int KD = KOPD + (int) *lMatch;
			lMatch += AlignStep;

			__m128i __KD = _mm_set1_epi32(KD);
			__m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
			__TransitionsD = _mm_add_epi32(__TransitionsD, __KD);
			__TransitionsD = _mm_insert_epi32(__TransitionsD, NLOW, EXTRA);
			
			// Paste index in lowest 2 bits
			__TransitionsD = _mm_slli_epi32(__TransitionsD, SCORE_SHIFT);
			__TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);
			
			// Get minimum ( this is SSE 4.1 )
			__TransitionsD = _mm_max_epi32( IOP_W[iprf].xmm, __TransitionsD);
			
			// Store into matrix
			assert(&MatrixPtr[iprf] >= matrix && &MatrixPtr[iprf] < &matrix[matrixLD*(LSEQ+1)]);
			_mm_store_si128(&MatrixPtr[iprf].xmm, __TransitionsD); 
			
			// Clean extra bits
			__TransitionsD = _mm_srai_epi32(__TransitionsD, SCORE_SHIFT);
			
			// Store IOPI and IOPM and IOPD
			IOP_W[iprf].xmm = __TransitionsD;
			
			// Set KOPD ( this is SSE 4.1 )
			KOPD = _mm_extract_epi32(__TransitionsD, DELETION);
		}
	}
}

static int GetBestScore(const union lScores * const restrict matrix,
                        const size_t SeqLength, const size_t prfLength)
{
	const size_t matrix_ld = prfLength + 1;
	const int (* restrict Scores)[4]  = (const int (* restrict)[4]) &matrix[matrix_ld+1];
	int extrema = NLOW;
	for (int i=0; i<SeqLength; ++i) {
		for (int j=0; j<prfLength; ++j) {
			const int Value = TO_SCORE(Scores[j][EXTRA]);
			if (Value >= extrema) {
				extrema = Value;
			}
		}
		Scores += matrix_ld;
	}
	return extrema;
}

static size_t GetBestAlignment(const union lScores * const restrict matrix, Alignment_t * const restrict Alignment,
                               const size_t SeqLength, const size_t prfLength)
{
	const size_t matrix_ld = prfLength + 1;
	const int (* restrict Scores)[4]  = (const int (* restrict)[4]) &matrix[matrix_ld+1];
	int index_r = -1, index_c = -1;
	int extrema = NLOW;
	for (int i=0; i<SeqLength; ++i) {
		for (int j=0; j<prfLength; ++j) {
			const int Value = TO_SCORE(Scores[j][EXTRA]);
			if (Value >= extrema) {
				extrema = Value;
				index_r = i;
				index_c = j;
			}
		}
		Scores += matrix_ld;
	}
	Alignment->Score = extrema;
	if (extrema == NLOW) return 0;
	/* Correct for initial entrance column and row in Scores array */
	index_c += 1;
	index_r += 1;
	
	Alignment->Matrix.row.End    = index_r;
	Alignment->Matrix.column.End = index_c;
	/* Correcting indexing to account extra top row and left column in Scores array */
	Alignment->Region.profile.End  = index_c - 1;
	Alignment->Region.sequence.End = index_r - 1;
	
	size_t counter     = 1; /* take into account initial X */
	int iprf           = index_c;
	size_t previousState, State = EXTRA, nCycles = 1UL;
	Scores  = (const int (* restrict)[4]) matrix;
	
	while ( iprf >= 0 && index_r >= 0) {
		previousState = State;
		const unsigned int MoveToState = (Scores[index_r*matrix_ld+iprf][State] & STATE_MASK);
		const unsigned int Cycle = (Scores[index_r*matrix_ld+iprf][State] & CYCLE_MASK);
		switch (MoveToState)
		{
			case PRIORITY_MATCH << STATE_SHIFT:
				if (Cycle) {
					iprf = prfLength;
					nCycles++;
				}
				else {
					iprf    -= 1;
					index_r -= 1;
				}
				State = MATCH;
				break;
			case PRIORITY_INSERTION << STATE_SHIFT:
				if (Cycle) { iprf = prfLength; nCycles++; }
				index_r -= 1;
				State = INSERTION;
				break;
			case PRIORITY_DELETION << STATE_SHIFT:
				if(Cycle) {
					iprf = prfLength;
					nCycles++;
				}
				else {
					--iprf;
				}
				State = DELETION;
				break;
			case PRIORITY_EXTRA << STATE_SHIFT:
				State = EXTRA;
				++counter;
				goto OUT;
				break;
		}
		++counter;
	}
	fprintf(stderr, "Potential issue as alignment does not start with an ENTRY point in GetBestAlignment!\n"
	                "Last position on matrix is iprf=%i, iseq=%i\n", iprf, index_r);
	exit(1);
	
	OUT:  
	Alignment->Matrix.column.Begin = iprf;
	Alignment->Matrix.row.Begin    = index_r;
	assert(Alignment->Matrix.column.Begin >= 0);
	assert(Alignment->Matrix.row.Begin  >= 0);
	
	/* Correcting indexing to account extra top row and left column */
	Alignment->Region.profile.Begin  = (previousState == INSERTION) ? iprf - 1: iprf /*-1+1*/;
	// 	Alignment->Region.sequence.Begin = (State == DELETION) ? index_r - 1 : index_r /*-1+1*/;
	Alignment->Region.sequence.Begin = index_r;
	assert(Alignment->Region.sequence.Begin >= 0);
	assert(Alignment->Region.profile.Begin >= 0);
	
	/* Number of cycles */
	Alignment->Cycles = (unsigned short int) (nCycles);
	
	return ++counter; // Account for terminal C code character \0
}

static size_t GetNextBestAlignment(const union lScores * const matrix, const _Bool * const HasAPath,
                                   Alignment_t * const Alignment,
                                   const size_t SeqLength, const size_t prfLength)
{
	const size_t matrix_ld = prfLength + 1;
	const int (* restrict Scores)[4]  = (const int (* restrict)[4]) &matrix[matrix_ld+1];
	int index_r = -1, index_c = -1;
	int extrema = NLOW;
	for (int i=0; i<SeqLength; ++i) {
		if (!HasAPath[i+1]) {
			for (int j=0; j<prfLength; ++j) {
				const int Value = TO_SCORE(Scores[j][EXTRA]);
				if (Value >= extrema) {
					extrema = Value;
					index_r = i;
					index_c = j;
				}
			}
		}
		Scores += matrix_ld;
	}
	Alignment->Score = extrema;
	if (extrema == NLOW) return 0;
	
	/* Correct for initial entrance column and row bypassed above */
	index_c += 1;
	index_r += 1;
	
	Alignment->Matrix.row.End      = index_r;
	Alignment->Matrix.column.End   = index_c;
	/* Correcting indexing to account extra top row and left column */
	Alignment->Region.profile.End  = index_c - 1;
	Alignment->Region.sequence.End = index_r - 1;

	size_t counter     = 1; /* take into account initial X */
	int iprf           = index_c;
	size_t State = EXTRA, nCycles = 1UL;
	Scores  = (const int (* restrict)[4]) matrix;
	
	while ( iprf >= 0 && index_r >= 0) {
		const unsigned int MoveToState = (Scores[index_r*matrix_ld+iprf][State] & STATE_MASK);
		const unsigned int Cycle = (Scores[index_r*matrix_ld+iprf][State] & CYCLE_MASK);
		switch (MoveToState)
		{
			case PRIORITY_MATCH << STATE_SHIFT:
				if (Cycle) {
					iprf = prfLength;
					nCycles++;
				}
				else {
					iprf    -= 1;
					index_r -= 1;
				}
				State = MATCH;
				break;
			case PRIORITY_INSERTION << STATE_SHIFT:
				// TODO: fix that !!! assert(Cycle == 0);
				if (Cycle) { iprf = prfLength; nCycles++; }
				index_r -= 1;
				State = INSERTION;
				break;
			case PRIORITY_DELETION << STATE_SHIFT:
				if(Cycle) {
					iprf = prfLength;
					nCycles++;
				}
				else {
					iprf -= 1;
				}
				State = DELETION;
				break;
			case PRIORITY_EXTRA << STATE_SHIFT:
				++counter;
				goto OUT;
				break;
		}
		++counter;
	}
	return 0UL;
	
	OUT:  
	Alignment->Matrix.column.Begin = iprf;
	Alignment->Matrix.row.Begin    = index_r;
	assert(Alignment->Matrix.column.Begin >= 0);
	assert(Alignment->Matrix.row.Begin  >= 0);
	
	/* Correcting indexing to account extra top row and left column */
	Alignment->Region.profile.Begin  = (State == INSERTION) ? iprf - 1: iprf /*-1+1*/;
	// 	Alignment->Region.sequence.Begin = (State == DELETION) ? index_r - 1 : index_r /*-1+1*/;
	Alignment->Region.sequence.Begin = index_r;
	assert(Alignment->Region.sequence.Begin >= 0);
	assert(Alignment->Region.profile.Begin >= 0);	

	/* Number of cycles */
	Alignment->Cycles = (unsigned short int) (nCycles);
	
	return ++counter; // Account for ternimal C code character \0
}

static void GetAlignmentSequence(const union lScores * const restrict matrix,
                                 const unsigned char * const restrict Sequence, unsigned char * const restrict AlignmentSequence,
                                 const Alignment_t * const restrict Alignment, const size_t OutputMemorySize, const size_t prfLength)
{
	const int (*Scores)[4] = (const int (*)[4]) matrix;
	
	int index          = Alignment->Matrix.row.End;
	size_t counter     = 0;
	int iprf           = Alignment->Matrix.column.End;
	size_t State       = EXTRA;
	const size_t ld    = prfLength+1;
	
	assert(iprf<=prfLength);
	while ( iprf >= 0 && index >= 0)
	{
		assert(counter < OutputMemorySize);
		const unsigned int Move = Scores[index*ld+iprf][State] & STATE_MASK;
		const unsigned int Cycle = (Scores[index*ld+iprf][State] & CYCLE_MASK);
		const unsigned char C = (Sequence[index-1] > 'Z') ? Sequence[index-1] - ('a' - 'A') : Sequence[index-1];
		switch (Move)
		{
			case PRIORITY_MATCH << STATE_SHIFT:
				if (Cycle) {
					iprf = prfLength;
				}
				else {
					--index;
					--iprf;
					AlignmentSequence[counter++] = C;
				}
				State = MATCH;
				break;
			case PRIORITY_INSERTION << STATE_SHIFT:
				if (Cycle) iprf = prfLength;
				--index;
				State = INSERTION;
				AlignmentSequence[counter++] = (char) (  C + ( 'a' - 'A'));
				break;
			case PRIORITY_DELETION << STATE_SHIFT:
				if(Cycle) {
					iprf = prfLength;
				}
				else {
					--iprf;
					AlignmentSequence[counter++] = '-';
				}
				State = DELETION;
				break;
			case PRIORITY_EXTRA << STATE_SHIFT:
				goto OUT;
				break;
		}
	}
	fprintf(stderr, "Potential issue as alignment does not start with an ENTRY point in GetBestSequence!\n"
	                "Last position on matrix is iprf=%i, iseq=%i\nSTA: %.*s\n", iprf, index, (int) counter, AlignmentSequence);
	exit(1);
	OUT:
	;
	
	/* Reverse the string */
	unsigned char * BackPtr = &AlignmentSequence[counter-1];
	
	for (size_t i=0; i<counter/2; ++i) {
		assert (BackPtr >= AlignmentSequence);
		const unsigned char c = AlignmentSequence[i];
		AlignmentSequence[i] = *BackPtr;
		*BackPtr-- = c;
	}
	if (counter & 0x1) {
		const unsigned char c = AlignmentSequence[counter/2];
		AlignmentSequence[counter/2] = *BackPtr;
		assert (BackPtr >= AlignmentSequence);
		*BackPtr = c;
	}
	assert(counter < OutputMemorySize);
	AlignmentSequence[counter] = '\0';
}

static int GetStateSequence(const union lScores * const matrix,
                            unsigned char * Alignment_Sequence, const union URegion * const Alignment, const size_t matrix_lda,
                            const size_t SeqLength)
{
	const size_t prfLength = matrix_lda - 1;
	unsigned char * SeqPtr = (Alignment_Sequence + SeqLength - 1);
	*SeqPtr-- = '\0';
	*SeqPtr-- = 'X';
	//unsigned char * const LastCharacter = SeqPtr;
	int counter = 1;
	int index = Alignment[1].End; // Removing first row score matrix offset
	int State = EXTRA;
	int iprf  = Alignment[0].End;
	while (1)
	{
		const unsigned int MoveToState = (matrix[index*matrix_lda+iprf].Element[State] & STATE_MASK);
		const unsigned int Cycle       = (matrix[index*matrix_lda+iprf].Element[State] & CYCLE_MASK);
		switch (MoveToState)
		{
			case PRIORITY_MATCH << STATE_SHIFT:
				if (Cycle) {
					iprf = prfLength;
					*SeqPtr-- = 'R';
				} else {
					*SeqPtr-- = 'M';
					iprf  -= 1;
					index -= 1;
				}
				State  = MATCH;
				break;
			case PRIORITY_INSERTION << STATE_SHIFT:
				if (Cycle) iprf = prfLength;
				*SeqPtr-- = 'I';
				index -= 1;
				State = INSERTION;
				break;
			case PRIORITY_DELETION << STATE_SHIFT:
				if(Cycle) {
					iprf = prfLength;
					*SeqPtr-- = 'R';
				}
				else {
					--iprf;
					*SeqPtr-- = 'D';
				}
				State = DELETION;
				break;
			case PRIORITY_EXTRA << STATE_SHIFT:
				*SeqPtr = 'X';
				counter++;
				State = EXTRA;
				goto OUT2;
				break;
		}
		counter++;
	}
	fputs("Potential issue as alignment does not start with an ENTRY point!\n", stderr);
	exit(1);
	
	OUT2: ;

	return counter;
}

static int GetAlignments(union lScores * const restrict matrix,
                         const struct Profile * const prf,
                         Alignment_t ** Alignments,
                         const size_t SeqLength)
{
	const char * const restrict ProfileSequence = prf->Sequence;
	const size_t prfLength = prf->Length;
	register const int (* restrict Scores)[4] = (const int (*)[4]) matrix;
	const size_t ld = prfLength+1;
	
	int Score;
	
	_Bool * const restrict HasAPath = (_Bool*) calloc(1+SeqLength,sizeof(_Bool));
	int AllocatedCount = 10;
	*Alignments = (struct Alignment*) malloc(AllocatedCount*sizeof(struct Alignment));
	if (!HasAPath || !(*Alignments)) {
		fputs("Unable to allocate memory\n", stderr);
		if (HasAPath) free(HasAPath);
		if (*Alignments) free(*Alignments);
		return -1;
	}
	unsigned char * StateSequence = NULL;
	size_t StateLength = 0;
	int count = 0;
	do {
		/* Get Best Alignment */
		struct Alignment lAlignment;
		const size_t Alignment_SeqLength = GetNextBestAlignment(matrix, HasAPath, &lAlignment, SeqLength, prfLength);
		if (Alignment_SeqLength == 0UL) break;
		lAlignment.Score = TO_SCORE(Scores[lAlignment.Matrix.row.End*ld+lAlignment.Matrix.column.End][EXTRA]);
		
		
		/* Prevent further alignment to use this region */
		_Bool AlreadyUsed = false;
		for(int k=lAlignment.Matrix.row.Begin; k<=lAlignment.Matrix.row.End; ++k) {
			AlreadyUsed |= HasAPath[k]; 
			HasAPath[k] = true;
		}
			
		if (!AlreadyUsed && (lAlignment.Score >= prf->CutOff) ) {
			
			/* Get Aligned Sequence */
			if (Alignment_SeqLength > StateLength) {
				StateLength = Alignment_SeqLength; 
				StateSequence = (unsigned char*) realloc(StateSequence, StateLength*sizeof(unsigned char));
				if (StateSequence == NULL) {
					fputs("Unable to allocate memory\n", stderr);
					if (*Alignments) free(*Alignments);
					return -2;
				}
			}
 			int counter = GetStateSequence(matrix, StateSequence, lAlignment.Matrix.array, ld, Alignment_SeqLength);

//  			printf("STATE: %s\n", StateSequence);
			
			// Correct beginning and ending
			switch(StateSequence[0]) {
				case 'M':
				case 'I':
					if (lAlignment.Matrix.row.Begin > 0) 
						lAlignment.Matrix.row.Begin += 1;
					else
						lAlignment.Matrix.row.Begin = 1;
					break;
				default:
					;
			}
			
			//lAlignment.Matrix.column.End -= 1-1; // substract left 3 columns, add 1 to get 1 based indexing
			lAlignment.Matrix.column.Begin = (lAlignment.Matrix.column.Begin > 1) ? lAlignment.Matrix.column.Begin-1 : 1;
			/* Compute motif start and end */
			lAlignment.IPMB = lAlignment.Matrix.column.Begin;
			lAlignment.IPME = lAlignment.Matrix.column.End - (int) prfLength - 1;
			
			if (OutputVerbose) {
				fprintf(stderr, "Alignment starts at (%i,%i) and ends at (%i,%i) with value %i and length %lu\n",
								lAlignment.Matrix.column.Begin, lAlignment.Matrix.row.Begin,
								lAlignment.Matrix.column.End, lAlignment.Matrix.row.End,
								lAlignment.Score, Alignment_SeqLength);
			}
			
			if (count == AllocatedCount) {
				AllocatedCount += 10;
				*Alignments = realloc(*Alignments, AllocatedCount*sizeof(struct Alignment));
				if (!(*Alignments)) {
					fputs("Unable to allocate memory\n", stderr);
					if (HasAPath) free(HasAPath);
					if (*Alignments) free(*Alignments);
					if (StateSequence) free(StateSequence);
					return -1;
				}
			}
			(*Alignments)[count++] = lAlignment;
		}
	} while (Score >= prf->CutOff);
	
	if (StateSequence) free(StateSequence);
	
	return count;
}


const Compute_t Repeat_sse41 = {
	.BuildMatrix = repeat_sse41,
	.GetBestScore = GetBestScore,
	.GetBestAlignment = GetBestAlignment,
	.GetNextBestAlignment = GetNextBestAlignment,
	.GetStateSequence = GetStateSequence,
	.GetAlignmentSequence = GetAlignmentSequence,
	.GetAlignments = GetAlignments,
	.ScoreMask = SCORE_MASK,
	.LeftNegativeScoreMask = LEFT_NEGATIVE_SCORE_MASK,
	.ScoreShift = SCORE_SHIFT,
	.StateShift = STATE_SHIFT,
	.StateMask = STATE_MASK,
	.MatrixColumnMultiplier = 1U,
	.MatrixExtraColumn = 1U,
	.MatrixExtraRow = 1U,
	.WorkCellSize = 4*sizeof(int)
};


int GetBestScoreAndCycles(const union lScores * const matrix, const size_t SeqLength, const size_t prfLength,
                          ScoreCycle_t * const restrict SC)
{
	/* Get the best alignment */
	int index_r = 0, index_c = 0;
	int extrema = NLOW;
	const size_t ld = prfLength+1;
	for (int i=1; i<=SeqLength; ++i) {
		for(int j=1; j<=prfLength; ++j) {
			const int val = TO_SCORE(matrix[i*ld+j].Element[EXTRA]);
			if (val >= extrema) {
				extrema = val;
				index_r = i;
				index_c = j;
			}
		}
	}
	SC->Score     = extrema;
	SC->End       = index_r;
	size_t Cycles = 1UL;
	int iprf      = index_c;
	size_t State  = EXTRA;

	int counter = 0;
	while ( iprf >= 0 && index_r >= 0)
	{
		const unsigned int Move = (matrix[index_r*ld+iprf].Element[State] & STATE_MASK);
		const unsigned int Cycle = (matrix[index_r*ld+iprf].Element[State] & CYCLE_MASK);
		switch (Move)
		{
			case PRIORITY_MATCH << STATE_SHIFT:
				if (Cycle) {
					iprf = prfLength;
					Cycles++;
				}
				else {
					--index_r;
					--iprf;
				}
				State = MATCH;
				counter++;
				break;
			case PRIORITY_INSERTION << STATE_SHIFT:
				if (Cycle) iprf = prfLength;
				--index_r;
				State = INSERTION;
				counter++;
				break;
			case PRIORITY_DELETION << STATE_SHIFT:
				if(Cycle) {
					iprf = prfLength;
					Cycles++;
				}
				else {
					--iprf;
				}
				State = DELETION;
				break;
			case PRIORITY_EXTRA << STATE_SHIFT:
				goto OUT;
				break;
			default:
				fprintf(stderr,"\nUnknown move (%u) encountered from cell (%i,%i)\n", Move, index_r, iprf);
				counter = -1;
				goto OUT;
				break;
		}
	}
	fputs("Potential issue as alignment does not start with an ENTRY point!\n", stderr);
	exit(1);
	
	OUT:
	SC->Cycles = Cycles;
	SC->Begin  = index_r;
	
	return counter;
}


#undef MAX
