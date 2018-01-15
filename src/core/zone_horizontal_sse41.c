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
#include <smmintrin.h>
#include <assert.h>
#include "pfCompute.h"
#include "pfOutput.h"
#include "sse41_inline_fcts.h"

#define MAX(a,b) ((a>b) ? a : b)
#if ( CPP_MATCH != 2 ) || ( CPP_INSERTION != 3 ) || ( CPP_DELETION != 0 ) || ( CPP_EXTRA != 1 )
#error "Zone method was hardcoded with different State ordering"
#endif

#define _my_cmp(a, b) _mm_cmpgt_epi32(a, b)
/* Priority:
 *	Extra then Insertion then match then deletion
 */

void GetZoneAlignment(const struct Profile * const restrict prf,
                      const unsigned char * const restrict Sequence,
                      Zone_t * const restrict Alignment, int * const restrict WORK,
                      const size_t BSEQ, const size_t LSEQ)
/*
 * WARNING: for SSE version, WORK should be aligned on cache size (64b) and
 *            4 times the (profile size + 1)*sizeof(int)
 *          + 63 to align to cache line
 *          + 4 times the (profile size + 1)*sizeof(int)
 * ---------------------------------------------------------
 *            8*(profile size + 1)*sizeof(int) + 63
 */
{
  int KOPD;
  const __m128i * restrict IOP_R;
  __m128i * restrict IOP_W = (__m128i*) WORK;
  
  register const TransitionScores * const restrict Transitions = prf->Scores.Insertion.Transitions;
  const StoredIntegerFormat * const restrict Match             = prf->Scores.Match.Alphabet;
  const StoredIntegerFormat * const restrict Insertion         = prf->Scores.Insertion.Alphabet;
  const size_t AlignStep                                       = prf->Scores.Match.AlignStep;
  static const union { int elem[4]; __m128i xmm; } One   = { .elem = { 1,1,1,1 } };
	static const union { int elem[4]; __m128i xmm; } Adder = { .elem = { 0,NLOW,1,1 } }; 
  const size_t prfLength = prf->Length;
  
  /*
   * Initialize Insertion and Match Entrance Line using FirstSequenceProtein
   */  
	{
		register const StoredIntegerFormat * restrict lMatch = (const StoredIntegerFormat *) &Match[_D];
		register const ScoreTuple * restrict FirstSequenceProtein = prf->Scores.Insertion.FirstSequenceProtein;
		const __m128i __zero = _mm_setzero_si128();
		/*
		* PROFILE COLUMN 0 entrance
		*/
		{
			__m128i __dummy = __zero;
			_mm_insert_epi32(__dummy, (int) FirstSequenceProtein[0].To[MATCH], CPP_MATCH);
			_mm_insert_epi32(__dummy, (int) FirstSequenceProtein[0].To[INSERTION], CPP_INSERTION);
			_mm_store_si128(&IOP_W[0], __dummy);
		}
		KOPD                           = (int) FirstSequenceProtein[0].To[DELETION];
		__m128i __FirstSequenceProtein = LoadStoredIntegerVector(&(FirstSequenceProtein[0]));

		FirstSequenceProtein++;
		register const TransitionScores (* restrict pTransitions) = &Transitions[1];

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
			// Get maximum
			__m128i __max = _mm_max_epi32(__TransitionsD, __FirstSequenceProtein);
			// Set KOPD
			KOPD = _mm_extract_epi32(__max, CPP_DELETION);
			// Set starting position
			__max = _mm_unpackhi_epi64(__zero, __max);
			// Store IOPI and IOPM
			_mm_store_si128( &IOP_W[iprf], __max);
		}
	}

  // Swap and assign Read and write pointers
	IOP_R = IOP_W;
	IOP_W = (__m128i*) (((uintptr_t) &WORK[4*(prf->Length+1)] + 63) & ~63);
	__m128i __Scores   = _mm_set1_epi32(NLOW);
	__m128i __OutIndex = _mm_set1_epi32(-1);
	__m128i __InIndex  = __OutIndex;
  /*
   * LOOP THROUGH THE SEQUENCE STRING
   */
	__m128i __CurrentSequenceIndex = _mm_set1_epi32((int) BSEQ);
	for ( size_t iseq=BSEQ; iseq < LSEQ-1; ++iseq) {
		register const size_t j1 = (size_t) Sequence[iseq];
		register const StoredIntegerFormat * restrict lInsertion = Insertion;
		__m128i __Begin = _mm_load_si128(&IOP_R[0]);
		
		int KOPM = _mm_extract_epi32(__Begin, CPP_MATCH);
		/*
		 * PROFILE COLUMN 0 entrance
		 */
		{
			register const int KI = _mm_extract_epi32(__Begin, CPP_INSERTION) + (int) lInsertion[j1];
			// Transform KI into a vector
			__m128i __KI = _mm_set1_epi32(KI);
			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[0].From[INSERTION]));
			// Add KI to Transition
			__TransitionsI = _mm_add_epi32(__TransitionsI, __KI);
				// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[0].From[EXTRA]));
			// Set Insertion indices
			const __m128i __Iindices = _mm_shuffle_epi32(__Begin, 0b01010101);
			// Set Extra indices
			const __m128i __Xindices = _mm_add_epi32(__CurrentSequenceIndex, Adder.xmm);
			// Insert lScore into __TransitionsX not necessary as profile loading store NLOW there automatically
			//__TransitionsX = _mm_insert_epi32(__TransitionsX, NLOW, EXTRA);
			const __m128i __isBigger = _my_cmp(__TransitionsI, __TransitionsX);
			__TransitionsX = _mm_blendv_epi8(__TransitionsX, __TransitionsI, __isBigger); 
			__m128i __B    = _mm_blendv_epi8(__Xindices, __Iindices, __isBigger);
			// Store KOPD
			KOPD = _mm_extract_epi32(__TransitionsX, CPP_DELETION);
			// Store IOPI and IOPM
			__TransitionsX = _mm_unpackhi_epi64(__B, __TransitionsX); 
			_mm_store_si128( &IOP_W[0], __TransitionsX);
		}

		lInsertion += AlignStep;
		register const StoredIntegerFormat * restrict lMatch = Match;

		/*
			* LOOP THROUGH THE REST OF THE PROFILE
			*/
		size_t iprf = 1UL;
		while (iprf<prfLength) {
			const int KM = KOPM + (int) lMatch[j1];
			const int KD = KOPD + (int) lMatch[_D];
			__m128i __IOP = _mm_load_si128(&IOP_R[iprf]);
			
			const int KI = _mm_extract_epi32(__IOP, CPP_INSERTION) + (int) lInsertion[j1];
			KOPM         = _mm_extract_epi32(__IOP, CPP_MATCH);

			lMatch     += AlignStep;
			lInsertion += AlignStep;

			// Transform KM into a vector
			__m128i __KM = _mm_set1_epi32(KM);
			// Set Match indices
			const __m128i __Mindices = _mm_shuffle_epi32(__Begin, 0b00000000);
			__Begin = __IOP;
			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsM = LoadStoredIntegerVector(&(Transitions[iprf].From[MATCH]));
			// Add KM to Transition
			__TransitionsM = _mm_add_epi32(__TransitionsM, __KM);
			// Transform KI into a vector
			__m128i __KI = _mm_set1_epi32(KI);
			// Set Insertion indices
			const __m128i __Iindices = _mm_shuffle_epi32(__IOP, 0b01010101);
			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[iprf].From[INSERTION]));
			// Add KI to Transition
			__TransitionsI = _mm_add_epi32(__TransitionsI, __KI);
			// Get maximum
			const __m128i __isBiggerMI = _my_cmp(__TransitionsM, __TransitionsI); 
			const __m128i __maxMI = _mm_blendv_epi8(__TransitionsI, __TransitionsM, __isBiggerMI);
			const __m128i __MI_indices = _mm_blendv_epi8(__Iindices, __Mindices, __isBiggerMI);
			
			// Load Transitions and Convert signed WORD into signed DWORD
			const __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));
			// Insert lscore into TransitionX not necessary as profile loading should have already done it
			// __TransitionsX = _mm_insert_epi32(__TransitionsX, lScore, EXTRA);
			// Set Extra indices
			const __m128i __Xindices = _mm_add_epi32(__CurrentSequenceIndex, Adder.xmm);
			// Transform KD into a vector
			const __m128i __KD = _mm_set1_epi32(KD);
			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
			// Add KD to Transition
			__TransitionsD = _mm_add_epi32(__TransitionsD, __KD);
			// Get maximum
			const __m128i __isBiggerDX = _my_cmp(__TransitionsX, __TransitionsD);
			const __m128i __maxDX = _mm_blendv_epi8(__TransitionsD, __TransitionsX, __isBiggerDX);
			const __m128i __DXindices = _mm_blendv_epi8(__CurrentSequenceIndex, __Xindices, __isBiggerDX);

			const __m128i __isBigger = _my_cmp(__maxDX, __maxMI);
			const __m128i __max = _mm_blendv_epi8(__maxMI, __maxDX, __isBigger);
			const __m128i __B = _mm_blendv_epi8(__MI_indices, __DXindices, __isBigger);
			// Set KOPD
			KOPD = _mm_extract_epi32(__max, CPP_DELETION);
			// Store IOPI and IOPM
			_mm_store_si128( &IOP_W[iprf], _mm_unpackhi_epi64(__B, __max));
			++iprf;
		}
		/* Last profile position */
    {
			const int KM = KOPM + (int) lMatch[j1];
			const int KD = KOPD + (int) lMatch[_D];
			__m128i __IOP = _mm_load_si128(&IOP_R[iprf]);
			const int KI = _mm_extract_epi32(__IOP, CPP_INSERTION) + (int) lInsertion[j1];
			
			// Transform KM into a vector
			__m128i __KM = _mm_set1_epi32(KM);
			// Set Match indices
			const __m128i __Mindices = _mm_shuffle_epi32(__Begin, 0b00000000);
			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsM = LoadStoredIntegerVector(&(Transitions[iprf].From[MATCH]));
			// Add KM to Transition
			__TransitionsM = _mm_add_epi32(__TransitionsM, __KM);
			// Transform KI into a vector
			__m128i __KI = _mm_set1_epi32(KI);
			// Set Insertion indices
			__m128i __Iindices = _mm_shuffle_epi32(__IOP, 0b01010101);
			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[iprf].From[INSERTION]));
			// Add KI to Transition
			__TransitionsI = _mm_add_epi32(__TransitionsI, __KI);
			// Get maximum
			const __m128i __isBiggerMI = _my_cmp(__TransitionsM, __TransitionsI);
			__m128i __maxMI = _mm_blendv_epi8(__TransitionsI, __TransitionsM, __isBiggerMI);
			__m128i __MI_indices = _mm_blendv_epi8(__Iindices, __Mindices, __isBiggerMI);
			
			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));
			// Insert lscore into TransitionX not necessary as profile loading should have already done it
			// __TransitionsX = _mm_insert_epi32(__TransitionsX, lScore, EXTRA);
			const __m128i __Xindices = _mm_add_epi32(__CurrentSequenceIndex, Adder.xmm);
			// Transform KD into a vector
			__m128i __KD = _mm_set1_epi32(KD);
			// Load Transitions and Convert signed WORD into signed DWORD
			__m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
			// Add KD to Transition
			__TransitionsD = _mm_add_epi32(__TransitionsD, __KD);
			// Get maximum
			__m128i __isBiggerDX = _my_cmp(__TransitionsX, __TransitionsD);
			__m128i __maxDX = _mm_blendv_epi8(__TransitionsD, __TransitionsX, __isBiggerDX);
			__m128i __DXindices = _mm_blendv_epi8(__CurrentSequenceIndex, __Xindices, __isBiggerDX);

			__m128i __isBigger = _my_cmp(__maxDX, __maxMI);
			__m128i __max = _mm_blendv_epi8(__maxMI, __maxDX, __isBigger);
			const __m128i __B = _mm_blendv_epi8(__MI_indices, __DXindices, __isBigger);
			// Update scores 
			const __m128i __isBetter = _my_cmp(__max, __Scores);
			__Scores   = _mm_blendv_epi8(__Scores, __max, __isBetter);
			__OutIndex = _mm_blendv_epi8(__OutIndex, __CurrentSequenceIndex, __isBetter);
			__InIndex  = _mm_blendv_epi8(__InIndex, __B, __isBetter);
			// Store IOPI and IOPM
			_mm_store_si128(&IOP_W[iprf], _mm_unpackhi_epi64(__B, __max));
    }
		__CurrentSequenceIndex = _mm_add_epi32(__CurrentSequenceIndex, One.xmm);
    
		// Swap Read and Write pointers
    const register __m128i * const ptr = IOP_W;
    IOP_W = (__m128i *) IOP_R;
    IOP_R = ptr;
  }

  /*
   * Last position on the Sequence using LastSequenceProtein
   */
#if 0
  {
    register const StoredIntegerFormat * restrict lInsertion = Insertion;
    const int j1 = (int) Sequence[LSEQ-1];
    /*
     * PROFILE COLUMN 0 entrance
     */
    int KOPM = IOP_R[0].Element.ScoreMatch;
    int KI   = IOP_R[0].Element.ScoreInsertion + (int) lInsertion[j1];

    KOPD = MAX( KI + (int) Transitions[0].Element[_ID], (int) Transitions[0].Element[_XD] );
    register const ScoreTuple * const restrict LastSequenceProtein = prf->Scores.Insertion.LastSequenceProtein;

    register const StoredIntegerFormat * restrict lMatch = Match;
    lInsertion += AlignStep;

		/*
     * LOOP THROUGH THE REST OF THE PROFILE
     */
		size_t iprf = 1UL;
    while(iprf<prfLength) {
      const int KM = KOPM                          + lMatch[j1];
      KI           = IOP_R[iprf].Element.ScoreInsertion + lInsertion[j1];
      const int KD = KOPD                          + lMatch[_D];

      lMatch     += AlignStep;
      lInsertion += AlignStep;

      KOPM = IOP_R[iprf].Element.ScoreMatch;
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
      // Get maximum ( this is SSE 4.1 )
      __m128i __max1 = _mm_max_epi32(__TransitionsM, __TransitionsI);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));
      // Insert lscore into TransitionX not necessary as profile loading should have already done it
      // __TransitionsX = _mm_insert_epi32(__TransitionsX, lScore, EXTRA);
      // Transform KD into a vector
      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
      // Insert LastProteinSequence
      __TransitionsD = _mm_insert_epi32(__TransitionsD, (int) LastSequenceProtein[iprf].From[DELETION], DUMMY);
      // Add KD to Transition
      __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);
      // Get maximum ( this is SSE 4.1 )
      __m128i __max2 = _mm_max_epi32(__TransitionsD, __TransitionsX);
      __max1 = _mm_max_epi32(__max2, __max1);
      // Set KOPD ( this is SSE 4.1 )
      KOPD = _mm_extract_epi32(__max1, DELETION);
			++iprf;
    }
    /* Last profile position */
		{
      const int KM = KOPM                          + lMatch[j1];
      KI           = IOP_R[iprf].Element.ScoreInsertion + lInsertion[j1];
      const int KD = KOPD                          + lMatch[_D];

      KOPM = IOP_R[iprf].Element.ScoreMatch;
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
      // Get maximum ( this is SSE 4.1 )
      __m128i __max1 = _mm_max_epi32(__TransitionsM, __TransitionsI);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));
      // Insert lscore into TransitionX not necessary as profile loading should have already done it
      // __TransitionsX = _mm_insert_epi32(__TransitionsX, lScore, EXTRA);
      // Transform KD into a vector
      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
      // Insert LastProteinSequence
      __TransitionsD = _mm_insert_epi32(__TransitionsD, (int) LastSequenceProtein[iprf].From[DELETION], DUMMY);
      // Add KD to Transition
      __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);
      // Get maximum ( this is SSE 4.1 )
      __m128i __max2 = _mm_max_epi32(__TransitionsD, __TransitionsX);
      __max1 = _mm_max_epi32(__max2, __max1);
 			// Update scores 
			const __m128i __isBetter = _my_cmp(__max1, __Scores);
			__Scores   = _mm_blendv_epi8(__Scores, __max1, __isBetter);
			__OutIndex = _mm_blendv_epi8(__OutIndex, __CurrentSequenceIndex, __isBetter);
    }		
  }
#endif
	/*
	 * Get best score
	 */
	Alignment->Score = _mm_extract_epi32(__Scores, CPP_EXTRA);
	Alignment->Begin = _mm_extract_epi32(__InIndex, CPP_EXTRA);
	Alignment->End   = _mm_extract_epi32(__OutIndex, CPP_EXTRA);
}

#undef MAX

