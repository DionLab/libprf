/*******************************************************
                        PFTOOLS
 *******************************************************
  Nov 22, 2016 pfStatistics.h
 *******************************************************
 (C) 2011-2016 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#ifndef _STATISTICS_H
#define _STATISTICS_H
#include <mm_malloc.h>
#include "pfConfig.h"

/*
 ************************************************************************************************
 *                                    DEFINITIONS                                               *
 ************************************************************************************************
 */
typedef struct ContengencyTable_s {
	unsigned int TP;
	unsigned int FP;
	unsigned int TN;
	unsigned int FN;
	
} ContengencyTable_t;


/*
 ************************************************************************************************
 *                                  VARIABLES REQUIRED                                          *
 ************************************************************************************************
 */

/*
 ************************************************************************************************
 *                                 FUNCTION DECLARATIONS                                        *
 ************************************************************************************************
 */
ContengencyTable_t * computeContengency(const unsigned char * restrict Data, 
                                        const unsigned char * restrict PacBioBase,
                                        const unsigned char * const restrict Threshold,
                                        const unsigned char Channel, 
                                        const size_t nData, const size_t nThreshold);

unsigned char getIQRLimit(unsigned int * const restrict Histogram, const size_t HistogramSum,
                          const float * const restrict DecodeTable, const float coef);


void computeHistogram(const unsigned char * restrict Data, unsigned int * const restrict Histogram,
                      const size_t nData);

/*
 ************************************************************************************************
 *                                   INLINE FUNCTIONS                                           *
 ************************************************************************************************
 */
static inline __ALWAYS_INLINE void FreeContengencyTable(ContengencyTable_t * CT)
{
    _mm_free(CT);
		CT = NULL;
}



#endif /* _STATISTICS_H */
