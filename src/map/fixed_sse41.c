#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <smmintrin.h>
#include <assert.h> 
#include "pfMap.h"

/* IOP structure is 
 * BEGIN OF INSERTION / BEGIN OF MATCH / INSERTION SCORE / MATCH SCORE
 */

#define IOP_MATCH_BEGIN         0x2
#define IOP_INSERTION_BEGIN     0x3
#define IOP_ALL_MATCH_SCORE     0b00000000
#define IOP_ALL_INSERTION_SCORE 0b01010101
#define IOP_ALL_MATCH_BEGIN     0b10101010
#define IOP_ALL_INSERTION_BEGIN 0b11111111


#define ALL_DELETION ((DELETION<<6)|(DELETION<<4)|(DELETION<<2)|(DELETION))

static const __my128i MatchMask = { 
	.elem = {
		PRIORITY_MATCH << STATE_SHIFT,
		PRIORITY_MATCH << STATE_SHIFT,
		PRIORITY_MATCH << STATE_SHIFT,
		PRIORITY_MATCH << STATE_SHIFT
	}
};
static const __my128i InsertionMask = {
	.elem = {
		PRIORITY_INSERTION << STATE_SHIFT,
		PRIORITY_INSERTION << STATE_SHIFT,
		PRIORITY_INSERTION << STATE_SHIFT,
		PRIORITY_INSERTION << STATE_SHIFT
	}
};
static const __my128i DeletionMask = { 
	.elem = {
		PRIORITY_DELETION << STATE_SHIFT,
		PRIORITY_DELETION << STATE_SHIFT,
		PRIORITY_DELETION << STATE_SHIFT,
		PRIORITY_DELETION << STATE_SHIFT
	}
};
static const __my128i ClearMask = { 
	.elem = {CLEAR_MASK, CLEAR_MASK, CLEAR_MASK, CLEAR_MASK}
};
static const __my128i nlow = { 
	.elem = { (unsigned int)NLOW<<SCORE_SHIFT, (unsigned int)NLOW<<SCORE_SHIFT,
	          (unsigned int)NLOW<<SCORE_SHIFT, (unsigned int)NLOW<<SCORE_SHIFT}
};
static const union { unsigned int elem[4]; __m128i xmm; } BitOnes = {
	.elem = { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF }
};
static const __my128i Ones = { .elem = { 1, 1, 1, 1 } };

/*
 * WARNING: for SSE version, WORK should be aligned on cache size (64b) and
 *            2 times the (Read Length)*sizeof(__m128i)
 */
void GetMapping(const struct Map * const restrict map, void * const restrict WorkSpace,
                const char * const restrict Genome, const char * const restrict Tag,
                Zone_t * const restrict zone, const size_t TagLength, const size_t GenomeLength)
{
	const __m128i * restrict IOP_R;
	__m128i * restrict IOP_W = (__m128i *) WorkSpace;
	
	
	/*
	 * Initialize Insertion and Match Entrance Line using FirstSequenceProtein
	 */  
	{
		/*
		 * PROFILE COLUMN 0 entrance
		 */
		const register __m128i __Zeros = _mm_setzero_si128();
		__m128i __KOPD  = _mm_shuffle_epi32(map->RowInit.xmm, ALL_DELETION);
		/*
		 * LOOP THROUGH THE REST OF THE PROFILE
		 */
		
		const register __m128i __Deletions = map->DeletionScores.xmm;
		for (size_t iprf=0; iprf<TagLength; ++iprf ) {
			__KOPD = _mm_and_si128(__KOPD, ClearMask.xmm);	
			__m128i __TransitionsD = _mm_add_epi32(__Deletions, __KOPD);
			__KOPD = _mm_shuffle_epi32(__TransitionsD, ALL_DELETION);
			_mm_store_si128(&IOP_W[iprf], _mm_unpackhi_epi64(__TransitionsD, __Zeros));
		}
	}
		
	// Swap and assign Read and write pointers
	IOP_R = IOP_W;
	IOP_W = &((__m128i*) WorkSpace)[TagLength];
	
	/*
	 * LOOP THROUGH THE SEQUENCE STRING
	 */
		const register __m128i __Zeros = _mm_setzero_si128();
	__m128i __BestScore = nlow.xmm;
	__m128i __BestLocationEnd = BitOnes.xmm;
	__m128i __BestLocationStart = __Zeros;
	__m128i __location = __Zeros;
	const size_t MinTagLength = TagLength - 1UL;
	for ( size_t iseq=0; iseq < GenomeLength; iseq++) {
		const register char SequenceDNA = Genome[iseq];
		__m128i __KOPM = _mm_unpackhi_epi64(map->RowInit.xmm, __location);
		/*
		 * PROFILE COLUMN 0 entrance
		 */
		__m128i __KOPD = _mm_shuffle_epi32(map->RowInit.xmm, ALL_DELETION);
		__m128i __BeginD = __location;		
		/*
		 * LOOP THROUGH THE REST OF THE PROFILE
		 */
		size_t iprf;
		if (SequenceDNA != 'N') {
			for (iprf=0; iprf<MinTagLength; iprf++ ) {
				const int ReadIsN        = (Tag[iprf] == 'N') ? 1 : 0;
				const int DoWeMatch      = (Tag[iprf] == SequenceDNA) ? 1 : 0;
				const __m128i __movemask = (DoWeMatch || ReadIsN) ? BitOnes.xmm : __Zeros;
				const __m128i __AddtoM   = _mm_blendv_epi8(map->MismatchScores.xmm, map->MatchScores.xmm, __movemask);
				
				__m128i __KI = _mm_shuffle_epi32(IOP_R[iprf], IOP_ALL_INSERTION_SCORE);
				__m128i __KM = _mm_shuffle_epi32(__KOPM, IOP_ALL_MATCH_SCORE);
				const __m128i __BeginI = _mm_shuffle_epi32(IOP_R[iprf], IOP_ALL_INSERTION_BEGIN);
				const __m128i __BeginM = _mm_shuffle_epi32(__KOPM, IOP_ALL_MATCH_BEGIN);
				
				__KOPD = _mm_and_si128(__KOPD, ClearMask.xmm);
				__KM = _mm_and_si128(__KM, ClearMask.xmm);
				__KI = _mm_and_si128(__KI, ClearMask.xmm);
				
				__KM = _mm_add_epi32(__KM, __AddtoM);
				__KI = _mm_add_epi32(__KI, map->InsertionScores.xmm);
				
				__KOPM = IOP_R[iprf];
				
				const __m128i __BestOfMorI = _mm_cmpgt_epi32(__KM, __KI);
				__m128i __max   = _mm_blendv_epi8(__KI, __KM, __BestOfMorI);
				__m128i __Begin = _mm_blendv_epi8(__BeginI, __BeginM, __BestOfMorI);
				
				__m128i __KD = _mm_add_epi32(__KOPD, map->DeletionScores.xmm);
				__m128i __BestOfAll = _mm_cmpgt_epi32(__max, __KD);
				__max   = _mm_blendv_epi8(__KD, __max, __BestOfAll);
				__Begin = _mm_blendv_epi8(__BeginD, __Begin, __BestOfAll);
				
				_mm_store_si128(&IOP_W[iprf], _mm_unpackhi_epi64(__max, __Begin));
				
				__KOPD   = _mm_shuffle_epi32(__max, ALL_DELETION);
				__BeginD = _mm_shuffle_epi32(__Begin, ALL_DELETION);
			}
			
			{
				const int ReadIsN        = (Tag[iprf] == 'N') ? 1 : 0;
				const int DoWeMatch      = (Tag[iprf] == SequenceDNA) ? 1 : 0;
				const __m128i __movemask = ( DoWeMatch || ReadIsN) ? BitOnes.xmm : __Zeros;
				const __m128i __AddtoM   = _mm_blendv_epi8(map->MismatchScores.xmm, map->MatchScores.xmm, __movemask);
				
				__m128i __KI = _mm_shuffle_epi32(IOP_R[iprf], IOP_ALL_INSERTION_SCORE);
				__m128i __KM = _mm_shuffle_epi32(__KOPM, IOP_ALL_MATCH_SCORE);
				const __m128i __BeginI = _mm_shuffle_epi32(IOP_R[iprf], IOP_ALL_INSERTION_BEGIN);
				const __m128i __BeginM = _mm_shuffle_epi32(__KOPM, IOP_ALL_MATCH_BEGIN);
				
				__KOPD = _mm_and_si128(__KOPD, ClearMask.xmm);
				__KM = _mm_and_si128(__KM, ClearMask.xmm);
				__KI = _mm_and_si128(__KI, ClearMask.xmm);
				
				__KM = _mm_add_epi32(__KM, __AddtoM);
				__KI = _mm_add_epi32(__KI, map->InsertionScores.xmm);
				
				const __m128i __BestOfMorI = _mm_cmpgt_epi32(__KM, __KI);
				__m128i __max   = _mm_blendv_epi8(__KI, __KM, __BestOfMorI);
				__m128i __Begin = _mm_blendv_epi8(__BeginI, __BeginM, __BestOfMorI);
				
				__m128i __KD = _mm_add_epi32(__KOPD, map->DeletionScores.xmm);
				
				__m128i __BestOfAll = _mm_cmpgt_epi32(__max, __KD);
				__max   = _mm_blendv_epi8(__KD, __max, __BestOfAll);
				__Begin = _mm_blendv_epi8(__BeginD, __Begin, __BestOfAll);
				
				const __m128i __mask = _mm_cmpgt_epi32(__max, __BestScore);
				__BestScore         = _mm_blendv_epi8(__BestScore, __max, __mask);
				__BestLocationEnd   = _mm_blendv_epi8(__BestLocationEnd, __location, __mask);
				__BestLocationStart = _mm_blendv_epi8(__BestLocationStart, __Begin, __mask);
				
				_mm_store_si128(&IOP_W[iprf], _mm_unpackhi_epi64(__max, __Begin));
			}
		}
		else {
			for (iprf=0; iprf<MinTagLength; iprf++ ) {
				__m128i __KM = _mm_shuffle_epi32(__KOPM, IOP_ALL_MATCH_SCORE);
				__m128i __KI = _mm_shuffle_epi32(IOP_R[iprf], IOP_ALL_INSERTION_SCORE);
				const __m128i __BeginI = _mm_shuffle_epi32(IOP_R[iprf], IOP_ALL_INSERTION_BEGIN);
				const __m128i __BeginM = _mm_shuffle_epi32(__KOPM, IOP_ALL_MATCH_BEGIN);
				
				__KOPD = _mm_and_si128(__KOPD, ClearMask.xmm);
				__KM = _mm_and_si128(__KM, ClearMask.xmm);
				__KI = _mm_and_si128(__KI, ClearMask.xmm);
				
				__KM = _mm_add_epi32(__KM, map->SequenceNMatchScores.xmm);
				__KI = _mm_add_epi32(__KI, map->InsertionScores.xmm);
				
 				__KOPM = IOP_R[iprf];
				
				const __m128i __BestOfMorI = _mm_cmpgt_epi32(__KM, __KI);
				__m128i __max   = _mm_blendv_epi8(__KI, __KM, __BestOfMorI);
				__m128i __Begin = _mm_blendv_epi8(__BeginI, __BeginM, __BestOfMorI);
				
				__m128i __KD = _mm_add_epi32(__KOPD, map->DeletionScores.xmm);

				__m128i __BestOfAll = _mm_cmpgt_epi32(__max, __KD);
				__max   = _mm_blendv_epi8(__KD, __max, __BestOfAll);
				__Begin = _mm_blendv_epi8(__BeginD, __Begin, __BestOfAll);
	
				_mm_store_si128(&IOP_W[iprf], _mm_unpackhi_epi64(__max, __Begin));
				
				__KOPD   = _mm_shuffle_epi32(__max, ALL_DELETION);
				__BeginD = _mm_shuffle_epi32(__Begin, ALL_DELETION);
			}
			
			{
				__m128i __KM = _mm_shuffle_epi32(__KOPM, IOP_ALL_MATCH_SCORE);
				__m128i __KI = _mm_shuffle_epi32(IOP_R[iprf], IOP_ALL_INSERTION_SCORE);
				const __m128i __BeginI = _mm_shuffle_epi32(IOP_R[iprf], IOP_ALL_INSERTION_BEGIN);
				const __m128i __BeginM = _mm_shuffle_epi32(__KOPM, IOP_ALL_MATCH_BEGIN);
				
				__KOPD = _mm_and_si128(__KOPD, ClearMask.xmm);
				__KM = _mm_and_si128(__KM, ClearMask.xmm);
				__KI = _mm_and_si128(__KI, ClearMask.xmm);
				
				__KM = _mm_add_epi32(__KM, map->SequenceNMatchScores.xmm);
				__KI = _mm_add_epi32(__KI, map->InsertionScores.xmm);
				
				__KOPM = _mm_shuffle_epi32(IOP_R[iprf], IOP_ALL_MATCH_SCORE);
				
				const __m128i __BestOfMorI = _mm_cmpgt_epi32(__KM, __KI);
				__m128i __max   = _mm_blendv_epi8(__KI, __KM, __BestOfMorI);
				__m128i __Begin = _mm_blendv_epi8(__BeginI, __BeginM, __BestOfMorI);
				
				__m128i __KD = _mm_add_epi32(__KOPD, map->DeletionScores.xmm);

				__m128i __BestOfAll = _mm_cmpgt_epi32(__max, __KD);
				__max   = _mm_blendv_epi8(__KD, __max, __BestOfAll);
				__Begin = _mm_blendv_epi8(__BeginD, __Begin, __BestOfAll);
						
				const __m128i __mask = _mm_cmpgt_epi32(__max, __BestScore);
				__BestScore = _mm_blendv_epi8(__BestScore, __max, __mask);
				__BestLocationEnd  = _mm_blendv_epi8(__BestLocationEnd, __location, __mask);
				__BestLocationStart = _mm_blendv_epi8(__BestLocationStart, __Begin, __mask);
				
				_mm_store_si128(&IOP_W[iprf], _mm_unpackhi_epi64(__max, __Begin));
			}
		}
		
		__location = _mm_add_epi32(__location, Ones.xmm);
		
		// Swap Read and Write pointers
		const register __m128i * const ptr = IOP_W;
		IOP_W = (__m128i *) IOP_R;
		IOP_R = ptr;
		if ((iseq & 0x3FF) == 0UL) printf("Genome: %lu/%lu\r", iseq, GenomeLength);
	}
	
	__BestScore = _mm_srai_epi32(__BestScore, SCORE_SHIFT);
	const int MX = _mm_extract_epi32(__BestScore, MATCH);
	const int DX = _mm_extract_epi32(__BestScore, DELETION);
	if (MX >= DX) {
		zone->Score = MX;
		zone->End = _mm_extract_epi32(__BestLocationEnd, MATCH);
		zone->Begin = _mm_extract_epi32(__BestLocationStart, MATCH);
	}
	else {
		zone->Score = DX;
		zone->End = _mm_extract_epi32(__BestLocationEnd, DELETION);
		zone->Begin = _mm_extract_epi32(__BestLocationStart, DELETION);
	}
}
