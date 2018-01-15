
/***************************************************************************************************
                        PFTOOLS
 ***************************************************************************************************
  Jul 1, 2019 pfMapcInline.h
 ***************************************************************************************************
 (C) 2019 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 ***************************************************************************************************/
#include <stdlib.h>
#ifdef BUILD_LIBRARY
#define PF_STATIC_INLINE(type) PFEXPORT type
#define PF_EXTERN_INLINE(type) PFEXPORT type
#else
#define PF_STATIC_INLINE(type) static inline type __ALWAYS_INLINE
#define PF_EXTERN_INLINE(type) extern inline type __ALWAYS_INLINE
#endif

PF_STATIC_INLINE(const int *) createHeuristicScoreMatrix(const struct Map * const restrict map, const char * const restrict Sequence,
																												 const size_t Length)
{
	// Length size rounded to cache line
  const size_t Aligned_Length = (Length + 15) & ~15;
  int * restrict const TIMatch = _mm_malloc(Aligned_Length*6*sizeof(int),64);
	if (TIMatch == NULL) goto bail;
	
	register const int Match = map->MatchScore;
	register const int Mismatch = map->MismatchScore;
	
	size_t i = 0UL;
	while (i<Length) {
		TIMatch[i] = Mismatch;
		switch(Sequence[i]) {
			case 'A':
				TIMatch[1*Aligned_Length+i] = Match;
				TIMatch[2*Aligned_Length+i] = Mismatch;
				TIMatch[3*Aligned_Length+i] = Mismatch;
				TIMatch[4*Aligned_Length+i] = Mismatch;
				TIMatch[5*Aligned_Length+i] = Mismatch;
				break;
			case 'C':
				TIMatch[1*Aligned_Length+i] = Mismatch;
				TIMatch[2*Aligned_Length+i] = Match;
				TIMatch[3*Aligned_Length+i] = Mismatch;
				TIMatch[4*Aligned_Length+i] = Mismatch;
				TIMatch[5*Aligned_Length+i] = Mismatch;
			case 'G':
				TIMatch[1*Aligned_Length+i] = Mismatch;
				TIMatch[2*Aligned_Length+i] = Mismatch;
				TIMatch[3*Aligned_Length+i] = Match;
				TIMatch[4*Aligned_Length+i] = Mismatch;
				TIMatch[5*Aligned_Length+i] = Mismatch;
			case 'T':
				TIMatch[1*Aligned_Length+i] = Mismatch;
				TIMatch[2*Aligned_Length+i] = Mismatch;
				TIMatch[3*Aligned_Length+i] = Mismatch;
				TIMatch[4*Aligned_Length+i] = Match;
				TIMatch[5*Aligned_Length+i] = Mismatch;
			case 'N':
				TIMatch[1*Aligned_Length+i] = Mismatch;
				TIMatch[2*Aligned_Length+i] = Mismatch;
				TIMatch[3*Aligned_Length+i] = Mismatch;
				TIMatch[4*Aligned_Length+i] = Mismatch;
				TIMatch[5*Aligned_Length+i] = Match;
			default:
				TIMatch[1*Aligned_Length+i] = Mismatch;
				TIMatch[2*Aligned_Length+i] = Mismatch;
				TIMatch[3*Aligned_Length+i] = Mismatch;
				TIMatch[4*Aligned_Length+i] = Mismatch;
				TIMatch[5*Aligned_Length+i] = Mismatch;
		}
	}
	
bail:
	return TIMatch;
}

/*
PF_STATIC_INLINE(const int *) TransposeAndConvertMatchMatrix(const union Scores * const Matrices, const size_t Alphabet_Length,
									 const size_t Profile_Length)
{
  const size_t step = Matrices->Match.AlignStep;
  // Profile size rounded to cache line
  const size_t Aligned_Profile_Length = (Profile_Length+1 + 15) & ~15;
  int * restrict const TIMatch = _mm_malloc(Aligned_Profile_Length*Alphabet_Length*sizeof(int),64);

  if (TIMatch == NULL) return TIMatch;
  memset(TIMatch, 0, Aligned_Profile_Length*Alphabet_Length*sizeof(int));
  const register StoredIntegerFormat * restrict lMatch = Matrices->Match.Alphabet;
  for (size_t iprf = 0; iprf<Profile_Length; ++iprf) {
    for (size_t alpha=0; alpha <Alphabet_Length; ++alpha) {
      TIMatch[alpha*Aligned_Profile_Length+iprf] = (int) lMatch[alpha];
    }
    lMatch += step;
  }
  return TIMatch;
}

PF_EXTERN_INLINE(void) TransposeAndConvertMatchMatrixGivenMemory(int * const restrict TIMatch, const union Scores * const Matrices,
							     const size_t Alphabet_Length, const size_t Profile_Length,
							     const size_t Aligned_Profile_Length)
{
  const size_t step = Matrices->Match.AlignStep;

  memset(TIMatch, 0, Aligned_Profile_Length*Alphabet_Length*sizeof(int));
  const register StoredIntegerFormat * restrict lMatch = Matrices->Match.Alphabet;
  for (size_t iprf = 0; iprf<Profile_Length; ++iprf) {
    for (size_t alpha=0; alpha<Alphabet_Length; ++alpha) {
      TIMatch[alpha*Aligned_Profile_Length+iprf] = (int) lMatch[alpha];
    }
    lMatch += step;
  }
}

PF_STATIC_INLINE(float * ) TransposeAndConvertToFloatMatchMatrix(const union Scores * const Matrices, const size_t Alphabet_Length,
                                                            const size_t Profile_Length)
{
  const size_t step = Matrices->Match.AlignStep;
  // Profile size rounded to cache line boundary
  const size_t Aligned_Profile_Length = (Profile_Length+1 + 15) & ~15;
  float * restrict const TFMatch = _mm_malloc(Aligned_Profile_Length*Alphabet_Length*sizeof(float),64);

  if (TFMatch == NULL) return TFMatch;
  memset(TFMatch, 0, Aligned_Profile_Length*Alphabet_Length*sizeof(float));

  const register StoredIntegerFormat * restrict lMatch = Matrices->Match.Alphabet;
  for (size_t iprf = 0; iprf<Profile_Length; ++iprf) {
    for (size_t alpha=0; alpha<Alphabet_Length; ++alpha) {
      TFMatch[alpha*Aligned_Profile_Length+iprf] = (float) lMatch[alpha];
    }
    lMatch += step;
  }
  return TFMatch;
}

PF_EXTERN_INLINE(void) TransposeAndConvertToFloatMatchMatrixGivenMemory(float * const restrict TFMatch, const union Scores * const Matrices,
							            const size_t Alphabet_Length, const size_t Profile_Length,
							            const size_t Aligned_Profile_Length)
{
  const size_t step = Matrices->Match.AlignStep;

  memset(TFMatch, 0, Aligned_Profile_Length*Alphabet_Length*sizeof(int));

  const register StoredIntegerFormat * restrict lMatch = Matrices->Match.Alphabet;
  for (size_t iprf = 0; iprf<Profile_Length; ++iprf) {
    for (size_t alpha=0; alpha<Alphabet_Length; ++alpha) {
      TFMatch[alpha*Aligned_Profile_Length+iprf] = (float) lMatch[alpha];
    }
    lMatch += step;
  }
}*/
