#include "prf_config.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <smmintrin.h>
#include "pfStatistics.h" 

ContengencyTable_t * computeContengency(const unsigned char * restrict Data, 
                                        const unsigned char * restrict PacBioBase,
                                        const unsigned char * const restrict Threshold,
                                        const unsigned char Channel, 
                                        const size_t nData, const size_t nThreshold)
{
	/* Check if correctly aligned */
	if (((uintptr_t) Data & (uintptr_t) 0xF) || ((uintptr_t) PacBioBase & (uintptr_t) 0xF)) {
		fputs("Check the alignement of array arguments to be 16 byte aligned!\n", stderr);
		return NULL;
	}
	/* Allocate contengency tables */
	ContengencyTable_t * const restrict Results = (ContengencyTable_t *) _mm_malloc(nThreshold*sizeof(ContengencyTable_t),16);
	
	#define PRINT(__x) {\
		char i[16] __attribute__((aligned(16)));\
		_mm_store_si128(&i[0], __x);\
		for (size_t k=0;k<16;k++) fprintf(stderr,"\t%i", (int)i[k]);\
		fputc('\n',stderr);\
	}
	
	if (Results != NULL) {
// 		const size_t limitData           = nData & ~((size_t)15);
// 		const __m128i __RangeTranslate   = _mm_set1_epi8(0x80);
// 		const register __m128i __Channel = _mm_set1_epi8((Channel /*^ (0x8)*/) );
// 		fputs("Channel: ",stderr); PRINT(__Channel);fputc('\n',stderr);
// 		
// 		// NOTE: Aligned on cache block, check later if this is not spanning the entire cache level ?
// 		unsigned int * CTmem = alloca(2*nThreshold*sizeof(unsigned int) + 63);
// 		unsigned int (* const restrict CT)[2] = (unsigned int(*)[2])( (uintptr_t) CTmem & (uintptr_t) ~(63) );
// 		{
// 			size_t k = 0;
// 			register const __m128i __Zero = _mm_setzero_si128(); 
// 			do { 
// 				_mm_store_si128((__m128i*) &CT[k][0], __Zero);
// 				k += 2;
// 			} while (k < nThreshold);
// 			/*if (nThreshold & 0x1)*/ *((size_t*) CT[nThreshold-1]) = 0;
// 		}
		
		Results[0].TP = Results[0].FP = Results[0].TN = Results[0].FN = 0U;
		for (size_t k=0; k<nThreshold; k++) {
			register unsigned int TP, FP, TN, FN;
			TP = FP = TN = FN = 0U;
			for (size_t i=0; i<nData; i++) {
				const _Bool IsGoodChannel = PacBioBase[i] == Channel;
				if (Data[i] < Threshold[k]) {
					if (IsGoodChannel) 
						FN++;
					else
						TN++;
				}
				else {
					if (IsGoodChannel)
						TP++;
					else
						FP++;
				}
			}
			Results[k].TP = TP;
			Results[k].FP = FP;
			Results[k].TN = TN;
			Results[k].FN = FN;
		}
		
// 		{
// 			register size_t iData = 16; 
// 			while (iData < limitData) {
// 				fputs("Data   : ", stderr); for (size_t k=0;k<16;k++) fprintf(stderr,"\t%i", (int) Data[iData-16+k]);
// 				fputc('\n',stderr);
// 				const __m128i __Data          = _mm_xor_si128(*((__m128i*) &Data[iData-16]), __RangeTranslate);
// 				fputs("Tr Data: ", stderr); PRINT(__Data);
// 				fputs("Bases  : ", stderr);
// 				for (size_t k=0;k<16;k++) fprintf(stderr,"\t%c", (int) PacBioBase[iData-16+k]);
// 				fputc('\n',stderr);
// 				const __m128i __IsGoodChannel = _mm_cmpeq_epi8(*((__m128i*) &PacBioBase[iData-16]), __Channel);
// 				fputs("Is Good: ", stderr); PRINT(__IsGoodChannel);
// 				iData += 16;
// 				
// 				size_t iThreshold = 0;
// 				do {
// 					const __m128i __Threshold = _mm_xor_si128(_mm_set1_epi8(Threshold[iThreshold]), __RangeTranslate);
// 					fputs("thres  : ", stderr); PRINT(__Threshold);
// 					const __m128i __SignalLevelDOWN = _mm_cmplt_epi8(__Threshold, __Data);
// 					fputs("Level  : ", stderr); PRINT(__SignalLevelDOWN);
// 					
// 					
// 					int utmp = (unsigned int) _mm_movemask_epi8(_mm_and_si128(__IsGoodChannel, __SignalLevelDOWN));
// 					int vtmp = (unsigned int) _mm_movemask_epi8(_mm_andnot_si128(__IsGoodChannel, __SignalLevelDOWN)); 
// 					unsigned int countFP, countTN;
// 					__asm__ __volatile__ ("popcnt %1,%0;" : "=r"(countFP) : "r"(utmp));
// 					__asm__ __volatile__ ("popcnt %1,%0;" : "=r"(countTN) : "r"(vtmp));
// 					fprintf(stderr, "FP: %u TN: %u\n", countFP, countTN);
// 					CT[iThreshold][0] += countFP;
// 					CT[iThreshold][1] += countTN;
// 					fputs("", stderr);
// 				} while (++iThreshold < nThreshold);
// 			}
// 			
// 			const size_t diff = iData - limitData;
// 			if (diff) {
// 				const __m128i __Data          = _mm_xor_si128(_mm_loadu_si128((__m128i*) &Data[nData-16]), __RangeTranslate);
// 				const __m128i __IsGoodChannel = _mm_cmpeq_epi8(*((__m128i*) &PacBioBase[nData-16]), __Channel);
// 				unsigned int mask             = ~((1 << (unsigned int) diff) - 1);
// 				
// 				size_t iThreshold = 0;
// 				do {
// 					const __m128i __Threshold       = _mm_xor_si128(_mm_set1_epi8(Threshold[iThreshold]), __RangeTranslate);
// 					const __m128i __SignalLevelDOWN = _mm_cmplt_epi8(__Threshold, __Data);
// 					
// 					int utmp = (unsigned int) _mm_movemask_epi8(_mm_and_si128(__IsGoodChannel, __SignalLevelDOWN));
// 					int vtmp = (unsigned int) _mm_movemask_epi8(_mm_andnot_si128(__IsGoodChannel, __SignalLevelDOWN)); 
// 					utmp &= mask;
// 					vtmp &= mask;
// 					unsigned int countFP, countTN;
// 					__asm__ ("popcnt %1,%0;" : "=r"(countFP) : "r"(utmp));
// 					__asm__ ("popcnt %1,%0;" : "=r"(countTN) : "r"(vtmp));
// 					CT[iThreshold][0] += countFP;
// 					CT[iThreshold][1] += countTN;
// 					
// 				} while (++iThreshold < nThreshold);
// 			}
// 		}
/*		
		for(size_t k=0;k<nThreshold; k++) {
			Results[k].TP = nData - CT[k][0];
			Results[k].FP = CT[k][0];
			Results[k].TN = CT[k][1];
			Results[k].FN = nData - CT[k][1];
		}*/
		
	}
	
	return Results;
}

void computeHistogram(const unsigned char * restrict Data, unsigned int * const restrict Histogram,
                      const size_t nData)
{
	for (size_t k=0; k<256; k++) Histogram[k] = 0U;
	
	for (size_t k=0; k<nData; k++) {
		Histogram[Data[k]]++;
	}
}

unsigned char getIQRLimit(unsigned int * const restrict Histogram, const size_t HistogramSum,
                          const float * const restrict DecodeTable, const float coef)
{
	size_t LowerLimit = HistogramSum/4;
	size_t UpperLimit = (3*HistogramSum)/4;
	
	size_t k=0;
	size_t Sum = 0;
	
	for (; k<256; k++) {
		const size_t N = Sum + (size_t) Histogram[k];
		if (N > LowerLimit) break;
		Sum = N;
	}
	const float l25 = DecodeTable[k];
		
	for (; k<256; k++) {
		const size_t N = Sum + (size_t) Histogram[k];
		if (N > UpperLimit) break;
		Sum = N;
	}
	const float l75 = DecodeTable[k];
	
	const float IQR = l75 - l25;
	const float Limit = l75+coef*IQR;
	for (k=0; k<256; k++) if (DecodeTable[k] > Limit) break;
	
	fprintf(stderr,"IQR is %lf\nStop point would be %lf @ %zu\n", IQR, Limit, k-1);
	
	return (unsigned char) (k-1);
}
