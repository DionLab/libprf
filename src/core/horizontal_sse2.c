
/*******************************************************
                        PFTOOLS
 *******************************************************
  Oct 12, 2012 xali1_print_sse41.c
 *******************************************************
 (C) 2012 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include "pfCompute.h"
#include "pfOutput.h"
#include "sse2_inline_fcts.h"

#define TAG

/* Priority */
#if !defined(_BEST_IS_NEGATIVE_)
enum StatePriority {
  PRIORITY_MATCH     = 1,
  PRIORITY_INSERTION = 2,
  PRIORITY_DELETION  = 0,
  PRIORITY_EXTRA     = 3
};
#else
enum StatePriority {
  PRIORITY_MATCH     = 2,
  PRIORITY_INSERTION = 1,
  PRIORITY_DELETION  = 3,
  PRIORITY_EXTRA     = 0
};
#endif

/*
  * INTEGER FORMAT  -------------------------------------
  *
  * |xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|xx|
  *                               |  |
  *                               |  +-- State
  *                               +----- Score
  */
# define SCORE_SHIFT               2
# define SCORE_MASK                0xFFFFFFFC
# define STATE_SHIFT               0
# define STATE_MASK                0x00000003
# define LEFT_NEGATIVE_SCORE_MASK  0xC0000000
# define CLEAN_STATE               0xFFFFFFF3
// # define CLEAN_PHASE               0xFFFFFFFC
// # define CLEAN_STATE_AND_PHASE     0xFFFFFFF0


#define TO_SCORE(x) (int)(((x)<0) ? ((((unsigned int)(x)) >> SCORE_SHIFT) | LEFT_NEGATIVE_SCORE_MASK ) : (((unsigned int)x) >> SCORE_SHIFT))
#define MAX(a,b) ((a>b) ? a : b)

union sIOP {
    struct {
      int Match;
      int Insertion;
    } Element;
    __m64 mm;
};

#define operation(x, inout, row, column, lda) _mm_store_si128(&((inout)[(row)*(lda)+(column)]), x)

static void xali1_print_sse2(const struct Profile * const restrict prf,
                       const unsigned char * const restrict Sequence,
                       union lScores * const restrict matrix, int * const restrict WORK,
                       unsigned short int * const restrict NumberOfInsertions,
                       const size_t BSEQ, const size_t LSEQ)
/*
 * WARNING: for SSE version, WORK should be aligned on cache size (64b) and 4 times the (profile size + 1)*sizeof(int)
 *          + 63 to align to cache line
 *          matrix should be of size at least (profile size + 1)*(sequence size+1) and aligned on 16b.
 */
{
  int KOPD;
  const union sIOP * restrict IOP_R;
  union sIOP * restrict IOP_W = (union sIOP*) WORK;
  const __m128i __MatchMask     = _mm_set1_epi32(PRIORITY_MATCH);
  const __m128i __InsertionMask = _mm_set1_epi32(PRIORITY_INSERTION);
  const __m128i __DeletionMask  = _mm_set1_epi32(PRIORITY_DELETION);
  const __m128i __ExtraMask     = _mm_set1_epi32(PRIORITY_EXTRA);
  
  register const TransitionScores * const restrict Transitions = prf->Scores.Insertion.Transitions;
  const StoredIntegerFormat * const restrict Match             = prf->Scores.Match.Alphabet;
  const StoredIntegerFormat * const restrict Insertion         = prf->Scores.Insertion.Alphabet;
  const size_t AlignStep                                       = prf->Scores.Match.AlignStep;
  const size_t prfLength = prf->Length;
  
  /* Set matrix ptr according to BSEQ */
  __m128i * const restrict MatrixPtr = (BSEQ == 0) ? &matrix[1+prfLength].xmm : &matrix[BSEQ].xmm;
  
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
     * PROFILE COLUMN 0 entrance
     */
    IOP_W[0].Element.Match     = (int) FirstSequenceProtein[0].To[MATCH];
    IOP_W[0].Element.Insertion = (int) FirstSequenceProtein[0].To[INSERTION];
    KOPD                       = (int) FirstSequenceProtein[0].To[DELETION];


    {
			__m128i __FirstSequenceProtein = LoadStoredIntegerVector(&(FirstSequenceProtein[0]));
			__FirstSequenceProtein = _mm_slli_epi32(__FirstSequenceProtein, 2);
			__FirstSequenceProtein = _mm_or_si128(__FirstSequenceProtein, __ExtraMask);
			// Store all scores to matrix
			operation(__FirstSequenceProtein, MatrixPtr - (1+prfLength), 0, 0, 1+prfLength);
    }
    
    FirstSequenceProtein++;
    register const TransitionScores (* restrict pTransitions) = &Transitions[1];
    register union sIOP * restrict pIOP = &IOP_W[1];
    
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
      __m128i __FirstSequenceProtein = LoadStoredIntegerVector(&(FirstSequenceProtein[0]));

      // Move to next profile First Sequence
      FirstSequenceProtein++;
      
      // Paste index in lowest 2 bits
      __TransitionsD = _mm_slli_epi32(__TransitionsD, 2);
      __FirstSequenceProtein = _mm_slli_epi32(__FirstSequenceProtein, 2);
      __TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);
      __FirstSequenceProtein = _mm_or_si128(__FirstSequenceProtein, __ExtraMask);
      
      // Get maximum ( this is SSE 4.1 )
      __m128i __max = _my_max_epi32(__TransitionsD, __FirstSequenceProtein);

      // Store all scores to matrix
      operation(__max, MatrixPtr - (1+prfLength), 0, iprf, 1+prfLength);

      // Clean extra bits
      __max = _mm_srai_epi32(__max, 2);
      
      // Store IOPI and IOPM
      StoreMatchInsertion( &(pIOP->mm), (__m128) __max);
      pIOP++;

      // Set KOPD ( this is SSE 4.1 )
      KOPD = _mm_extract_epi32(__max, DELETION);
    }
  }

  // Swap and assign Read and write pointers
  IOP_R = IOP_W;
  IOP_W = (union sIOP*) (((uintptr_t) &WORK[2*(prf->Length+1)] + 63) & ~63);

  /*
   * LOOP THROUGH THE SEQUENCE STRING
   */
  for ( size_t iseq=BSEQ; iseq < LSEQ-1; ++iseq) {
    register const size_t j1 = (size_t) Sequence[iseq];
    int KOPM = IOP_R[0].Element.Match;
    register const StoredIntegerFormat * restrict lInsertion = Insertion;
    /*
     * PROFILE COLUMN 0 entrance
     */
    {
      register const int KI = IOP_R[0].Element.Insertion + (int) lInsertion[j1];

      // Transform KI into a vector
      __m128i __KI = _mm_set1_epi32(KI);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[0].From[INSERTION]));
     // Add KI to Transition
      __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

       // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[0].From[EXTRA]));

      // Insert lScore into __TransitionsX not necessary as profile loading store NLOW there automatically
      //__TransitionsX = _mm_insert_epi32(__TransitionsX, NLOW, DUMMY);
      
      // Paste index in lowest 2 bits
      __TransitionsI = _mm_slli_epi32(__TransitionsI, 2);
      __TransitionsX = _mm_slli_epi32(__TransitionsX, 2);
      __TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);
      __TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);

      // Get maximum ( this is SSE 4.1 )
      __m128i __max = _my_max_epi32(__TransitionsI, __TransitionsX);
      
      // Store all scores to matrix
      operation(__max, MatrixPtr, iseq, 0, 1+prfLength); //_mm_store_si128(&pmatrix[iprf-1], __max1);

      // Clean extra bits
      __max = _mm_srai_epi32(__max, 2);

			// Store IOPI and IOPM
      StoreMatchInsertion( &(IOP_W[0].mm), (__m128) __max);

      // Store KOPD
      KOPD = _mm_extract_epi32(__max, DELETION);

      // Backup new score to xmm register
      //lScore = _mm_extract_epi32(__max, DUMMY);
    }

    lInsertion += AlignStep;
    register const StoredIntegerFormat * restrict lMatch = Match;

    /*
     * LOOP THROUGH THE REST OF THE PROFILE
     */
    for (size_t iprf=1; iprf<=prfLength; ++iprf ) {
      const int KM = KOPM                          + (int) lMatch[j1];
      const int KI = IOP_R[iprf].Element.Insertion + (int) lInsertion[j1];
      const int KD = KOPD                          + (int) lMatch[_D];

      lMatch     += AlignStep;
      lInsertion += AlignStep;

      KOPM = IOP_R[iprf].Element.Match;

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
      __TransitionsM = _mm_slli_epi32(__TransitionsM, 2);
      __TransitionsI = _mm_slli_epi32(__TransitionsI, 2);
      __TransitionsM = _mm_or_si128(__TransitionsM, __MatchMask);
      __TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);

      // Get maximum ( this is SSE 4.1 )
      __m128i __max1 = _my_max_epi32(__TransitionsM, __TransitionsI);      

      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));
       // Insert lscore into TransitionX not necessary as profile loading should have already done it
      // __TransitionsX = _mm_insert_epi32(__TransitionsX, lScore, DUMMY);

      // Transform KD into a vector
      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
      // Add KD to Transition
      __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

      // Paste index in lowest 2 bits
      __TransitionsX = _mm_slli_epi32(__TransitionsX, 2);
      __TransitionsD = _mm_slli_epi32(__TransitionsD, 2);
      __TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);
      __TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);

      // Get maximum ( this is SSE 4.1 )
      __m128i __max2 = _my_max_epi32(__TransitionsD, __TransitionsX);
      __max1 = _my_max_epi32(__max1, __max2);

      // Store all scores to matrix
      operation(__max1, MatrixPtr, iseq, iprf, 1+prfLength); //_mm_store_si128(&pmatrix[iprf-1], __max1);

      // Clean extra bits
      __max1 = _mm_srai_epi32(__max1, 2);

      // Store IOPI and IOPM
      StoreMatchInsertion( &IOP_W[iprf].mm, (__m128) __max1);

      // Set KOPD ( this is SSE 4.1 )
      KOPD = _mm_extract_epi32(__max1, DELETION);

//       lScore = _mm_extract_epi32(__max1, DUMMY);

    }

    // Swap Read and Write pointers
    const register union sIOP * const ptr = IOP_W;
    IOP_W = (union sIOP *) IOP_R;
    IOP_R = ptr;
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
    int KOPM     = IOP_R[0].Element.Match;
    int KI       = IOP_R[0].Element.Insertion + (int) lInsertion[j1];

    KOPD = MAX( KI + (int) Transitions[0].Element[_ID], (int) Transitions[0].Element[_XD] );
    register const ScoreTuple * const restrict LastSequenceProtein = prf->Scores.Insertion.LastSequenceProtein;

    register const StoredIntegerFormat * restrict lMatch = Match;
    lInsertion += AlignStep;
    __m128i __Scores = _mm_set1_epi32((unsigned int) NLOW<<2);

		/*
     * LOOP THROUGH THE REST OF THE PROFILE
     */
    for (size_t iprf=1; iprf<=prfLength; ++iprf) {
      const int KM = KOPM                          + lMatch[j1];
      KI           = IOP_R[iprf].Element.Insertion + lInsertion[j1];
      const int KD = KOPD                          + lMatch[_D];

      lMatch     += AlignStep;
      lInsertion += AlignStep;

      KOPM = IOP_R[iprf].Element.Match;
#ifdef ONLY_NECESSARY
      int tM = ((KM + (int) Transitions[iprf].Element[_MD]) << 2 ) | PRIORITY_MATCH;
      int tX = ((     (int) Transitions[iprf].Element[_XD]) << 2 ) | PRIORITY_EXTRA;
      int tI = ((KI + (int) Transitions[iprf].Element[_ID]) << 2 ) | PRIORITY_INSERTION;
      int tD = ((KD + (int) Transitions[iprf].Element[_DD]) << 2 ) | PRIORITY_DELETION;
      int tIOPD1 = MAX(tM, tX);
      int tIOPD2 = MAX(tI, tD);
      KOPD     = MAX( tIOPD1, tIOPD2);
      __Scores = _mm_insert_epi32(__Scores, KOPD, DELETION);

      tM = ((KM + (int) LastSequenceProtein[iprf].From[MATCH])     << 2 ) | PRIORITY_MATCH;
      tI = ((KI + (int) LastSequenceProtein[iprf].From[INSERTION]) << 2 ) | PRIORITY_INSERTION;
      tD = ((KD + (int) LastSequenceProtein[iprf].From[DELETION])  << 2 ) | PRIORITY_DELETION;
      const int tIOPT1 = MAX(tM, tI);
      const int lScore = MAX(tIOPT1, tD);
      __Scores = _mm_insert_epi32(__Scores, lScore, DUMMY);

      operation(__Scores, MatrixPtr, LSEQ-1, iprf, prfLength+1); //_mm_store_si128(&pmatrix[iprf-1], __Scores);
      // Clean KOPD extra 2 bits
      KOPD = (KOPD < 0) ? (KOPD >> 2) | 0xC0000000 : (KOPD >> 2);
#else
      // Transform KM into a vector
      __m128i __KM = _mm_set1_epi32(KM);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsM = LoadStoredIntegerVector(&(Transitions[iprf].From[MATCH]));
      // Insert LastProteinSequence
      __TransitionsM = _mm_insert_epi32(__TransitionsM, (int) LastSequenceProtein[iprf].From[MATCH], DUMMY);
      // Add KM to Transition
      __TransitionsM = _mm_add_epi32(__TransitionsM, __KM);


      // Transform KI into a vector
      __m128i __KI = _mm_set1_epi32(KI);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[iprf].From[INSERTION]));
      // Insert LastProteinSequence
      __TransitionsI = _mm_insert_epi32(__TransitionsI, (int) LastSequenceProtein[iprf].From[INSERTION], DUMMY);
      // Add KI to Transition
      __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

      // Paste index in lowest 2 bits
      __TransitionsM = _mm_slli_epi32(__TransitionsM, 2);
      __TransitionsI = _mm_slli_epi32(__TransitionsI, 2);
      __TransitionsM = _mm_or_si128(__TransitionsM, __MatchMask);
      __TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);

      // Get maximum ( this is SSE 4.1 )
      __m128i __max1 = _my_max_epi32(__TransitionsM, __TransitionsI);

      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));
      // Insert lscore into TransitionX not necessary as profile loading should have already done it
      // __TransitionsX = _mm_insert_epi32(__TransitionsX, lScore, DUMMY);

      // Transform KD into a vector
      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
      // Insert LastProteinSequence
      __TransitionsD = _mm_insert_epi32(__TransitionsD, (int) LastSequenceProtein[iprf].From[DELETION], DUMMY);
      // Add KD to Transition
      __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

      // Paste index in lowest 2 bits
      __TransitionsX = _mm_slli_epi32(__TransitionsX, 2);
      __TransitionsD = _mm_slli_epi32(__TransitionsD, 2);
      __TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);
      __TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);

      // Get maximum ( this is SSE 4.1 )
      __m128i __max2 = _my_max_epi32(__TransitionsD, __TransitionsX);
      __max1 = _my_max_epi32(__max1, __max2);

      // Store all scores to matrix
      operation(__max1, MatrixPtr, LSEQ-1, iprf, prfLength+1); //_mm_store_si128(&pmatrix[iprf-1], __max1);

      // Clean extra bits
      __max1 = _mm_srai_epi32(__max1, 2);

      // Set KOPD ( this is SSE 4.1 )
      KOPD = _mm_extract_epi32(__max1, DELETION);
#endif
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
	
	/* Correct for initial entrance column and row */
	index_c += 1;
	index_r += 1;
	
	Alignment->Matrix.row.End      = index_r;
	Alignment->Matrix.column.End   = index_c;
	/* Correcting indexing to account extra top row and left column */
	Alignment->Region.profile.End  = index_c - 1;
	Alignment->Region.sequence.End = index_r - 1;
	
	size_t counter     = 1; /* take into account initial X */
	int iprf           = index_c;
	size_t State       = EXTRA;
	Scores  = (const int (* restrict)[4]) matrix;
	while ( iprf >= 0 && index_r >= 0) {
		const unsigned int MoveToState = (Scores[index_r*matrix_ld+iprf][State] & STATE_MASK);
		switch (MoveToState)
		{
			case PRIORITY_MATCH << STATE_SHIFT:
				iprf    -= 1;
				index_r -= 1;
				State = MATCH;
				break;
			case PRIORITY_INSERTION << STATE_SHIFT:
				index_r -= 1;
				State = INSERTION;
				break;
			case PRIORITY_DELETION << STATE_SHIFT:
					--iprf;
				State = DELETION;
				break;
			case PRIORITY_EXTRA << STATE_SHIFT:
				++counter;
				goto OUT;
				break;
		}
		++counter;
// 		index_r = (index_r < 0 ) ? 0 : index_r;
	}
	return 0;
	
	OUT:;
	Alignment->Matrix.column.Begin = iprf;
	Alignment->Matrix.row.Begin    = index_r;
	
	/* Correcting indexing to account extra top row and left column */
	Alignment->Region.profile.Begin  = (State == INSERTION) ? iprf - 1: iprf /*-1+1*/;
	Alignment->Region.sequence.Begin = (State == DELETION) ? index_r - 1 : index_r /*-1+1*/;
	
	return ++counter; // Account for ternimal C code character \0
}

static size_t GetNextBestAlignment(const union lScores * const matrix, const _Bool * const HasAPath, Alignment_t * const Alignment,
                                   const size_t SeqLength, const size_t prfLength)
{
	const size_t matrix_ld = prfLength + 1;
	const int (* restrict Scores)[4]  = (const int (* restrict)[4]) &matrix[matrix_ld+1];
	int index_r = -1, index_c = -1;
	int extrema = NLOW;
	for (int i=0; i<SeqLength; ++i) {
		if (!HasAPath[1+i]) {
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
	
	/* Correct for initial entrance column and row */
	index_c += 1;
	index_r += 1;
	
	Alignment->Matrix.row.End      = index_r;
	Alignment->Matrix.column.End   = index_c;
	/* Correcting indexing to account extra top row and left column */
	Alignment->Region.profile.End  = index_c - 1;
	Alignment->Region.sequence.End = index_r - 1;
	
	size_t counter     = 1; /* take into account initial X */
	int iprf           = index_c;
	size_t State       = EXTRA;
	Scores  = (const int (* restrict)[4]) matrix;
	while ( iprf >= 0 && index_r >= 0) {
		const unsigned int MoveToState = (Scores[index_r*matrix_ld+iprf][State] & STATE_MASK);
		switch (MoveToState)
		{
			case PRIORITY_MATCH << STATE_SHIFT:
				iprf    -= 1;
				index_r -= 1;
				State = MATCH;
				break;
			case PRIORITY_INSERTION << STATE_SHIFT:
				index_r -= 1;
				State = INSERTION;
				break;
			case PRIORITY_DELETION << STATE_SHIFT:
					--iprf;
				State = DELETION;
				break;
			case PRIORITY_EXTRA << STATE_SHIFT:
				++counter;
				goto OUT;
				break;
		}
		++counter;
// 		index_r = (index_r < 0 ) ? 0 : index_r;
	}
	return 0;
	
	OUT:;
	Alignment->Matrix.column.Begin = iprf;
	Alignment->Matrix.row.Begin    = index_r;
	
	/* Correcting indexing to account extra top row and left column */
	Alignment->Region.profile.Begin  = (State == INSERTION) ? iprf - 1: iprf /*-1+1*/;
	//Alignment->Region.sequence.Begin = (State == DELETION) ? index_r - 1 : index_r /*-1+1*/;
	
	return ++counter; // Account for ternimal C code character \0
}

static void GetAlignmentSequence(const union lScores * const restrict matrix,
			  const unsigned char * const restrict Sequence, unsigned char * const restrict AlignmentSequence,
			  const Alignment_t * const restrict Alignment, const size_t SeqLength, const size_t prfLength)
{
    const int (*Scores)[4]  = (const int (*)[4]) matrix;

    int index       = Alignment->Matrix.row.End;
    size_t counter  = 0;
    int iprf        = Alignment->Matrix.column.End;
    size_t State    = EXTRA;
    const size_t ld = prfLength+1;
    
    while ( iprf >= 0 && index >= 0)
    {
			const unsigned int Move = Scores[index*ld+iprf][State] & 0x3;
			const unsigned char C = (Sequence[index-1] > 'Z') ? Sequence[index-1] - ('a' - 'A') : Sequence[index-1];
			switch (Move)
			{
			  case PRIORITY_MATCH:
			    --index;
			    --iprf;
			    State = MATCH;
			    AlignmentSequence[counter] = C;
			    break;
			  case PRIORITY_INSERTION:
			    --index;
			    State = INSERTION;
			    AlignmentSequence[counter] = (char) (  C + ( 'a' - 'A'));
			    break;
			  case PRIORITY_DELETION:
			    --iprf;
			    State = DELETION;
			    AlignmentSequence[counter] = '-';
			    break;
			  case PRIORITY_EXTRA:
			    goto OUT;
			    break;
			  default:
			    fprintf(stderr,"\nUnknown move (%u) encountered from cell (%i,%i)\n", Move, index, iprf);
			    --iprf;
			    goto OUT;
			    break;
			}
			++counter;
    }
    fputs("Potential issue as alignment does not start with an ENTRY point in std GetAlignmentSequence!\n", stderr);

    OUT:;

    /* Reverse the string */
    unsigned char * BackPtr = &AlignmentSequence[counter-1];

    for (size_t i=0; i<counter/2; ++i) {
			const unsigned char c = AlignmentSequence[i];
			AlignmentSequence[i] = *BackPtr;
			*BackPtr-- = c;
    }
    if (counter & 0x1) {
			const unsigned char c = AlignmentSequence[counter/2];
			AlignmentSequence[counter/2] = *BackPtr;
			*BackPtr = c;
    }
    AlignmentSequence[counter] = '\0';
}

// This is probably BUGGY
static int GetStateSequence(const union lScores * const matrix,
                            unsigned char * Alignment_Sequence, const union URegion * const Alignment, const size_t matrix_lda,
                            const size_t SeqLength)
{
	const size_t prfLength = matrix_lda - 1;
	unsigned char * SeqPtr = (Alignment_Sequence + SeqLength - 1);
	*SeqPtr-- = '\0';
	*SeqPtr-- = 'X';
// 	unsigned char * const LastCharacter = SeqPtr;
	int counter = 2;
	int index = Alignment[1].End; // Removing first row score matrix offset
	int State = EXTRA;

	int iprf  = Alignment[0].End;
	while (1)
	{
		const unsigned int MoveToState = (matrix[index*matrix_lda+iprf].Element[State] & STATE_MASK);
		switch (MoveToState)
		{
			case PRIORITY_MATCH << STATE_SHIFT:
				*SeqPtr-- = 'M';
				iprf  -= 1;
				index -= 1;
				State  = MATCH;
				break;
			case PRIORITY_INSERTION << STATE_SHIFT:
				*SeqPtr-- = 'I';
				index -= 1;
				State = INSERTION;
				break;
			case PRIORITY_DELETION << STATE_SHIFT:
				--iprf;
				*SeqPtr-- = '-';
				State = DELETION;
				break;
			case PRIORITY_EXTRA << STATE_SHIFT:
				*SeqPtr = 'X';
				State = EXTRA;
				goto OUT2;
				break;
		}
		counter++;
		index = (index < 0 ) ? 0 : index;
		// 	printf("%s\n", SeqPtr+1);
	}
	fputs("Potential issue as alignment does not start with an ENTRY point in GetStateSequence!\n", stderr);
	exit(1);
	OUT2: ;
// 	const int Mod3End = Alignment[0].End % 3;
// 	switch (Mod3End) {
// 		case 2:
// 			if ((unsigned char) LastCharacter[-2] < 'a') {
// 				LastCharacter[-1] = '\0';
// 				counter -= 2;
// 			}
// 			break;
// 		case 1:
// 			if ((unsigned char) LastCharacter[-1] < 'a') {
// 				LastCharacter[0] = '\0';
// 				counter -= 1;
// 			}
// 			break;
// 		default:;
// 	}
	return counter;
}

static int GetAlignments(union lScores * const restrict matrix,
                         const struct Profile * const prf,
                         struct Alignment ** Alignments,
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
		Alignment_t lAlignment;
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


const Compute_t Standard_sse2 = {
	.BuildMatrix = xali1_print_sse2,
	.GetBestScore = GetBestScore,
	.GetBestAlignment = GetBestAlignment,
	.GetNextBestAlignment = GetNextBestAlignment,
	.GetStateSequence = GetStateSequence,
	.GetAlignmentSequence = GetAlignmentSequence,
	.GetAlignments = GetAlignments,
	.ScoreMask = SCORE_MASK,
	.LeftNegativeScoreMask = LEFT_NEGATIVE_SCORE_MASK,
	.ScoreShift = SCORE_SHIFT,
	.MatrixColumnMultiplier = 1U,
	.MatrixExtraColumn = 1U,
	.MatrixExtraRow = 1U,
	.WorkCellSize = 2*sizeof(int)
};

#undef MAX
