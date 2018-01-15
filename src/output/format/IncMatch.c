#include "prf_config.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "pfProfile.h"
#include "pfCompute.h"
#include "pfOutput.h"

/*
 * WARNING: alignement does not start at 0 but 1 !!!
 *          NEEDS TO BE FIXED SOME DAY
 */

void PrintIncMatch(const struct Profile * const prf, const unsigned char * const AlignedSequence,
                   const struct Alignment * const alignment, const char * const Header,
                   const size_t SequenceLength, const float RAVE)
{ // replicates old pfsearch -zkx output
    char * buffer = calloc(OutputPrintWidth+1,sizeof(char));
    
		const char * cptr = Header;
    while ( *cptr != ' ' && *cptr != '\0' ) ++cptr;
		int HeaderLength = (int) ((uintptr_t) cptr - (uintptr_t) Header) - 1;

    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;

		const float norm_score = (RawToNormalizedFunction == NULL)
												? 0.0f
												: RawToNormalizedFunction(alignment->Score, NormCoefs, RAVE, SequenceLength);
		fprintf(stdout, ">%.*s/%i-%i motif=%s|%s norm_score=%.3f raw_score=%i motif_start=%i motif_end=%i\n",
				HeaderLength, &Header[1],
				alignment->Matrix.row.Begin,
				alignment->Matrix.row.End,
				prf->AC_Number,
				prf->Identification,
				norm_score,
				alignment->Score,
				alignment->IPMB,
				alignment->IPME
				);

		*buffer = '\0';
		size_t offset = 0;
	
		while ( offset <= SequenceLength ) {
				strncpy(buffer,AlignedSequence+offset,OutputPrintWidth);
				if ( '\0' != *buffer ) {
						fputs(buffer,stdout);
						fputc('\n',stdout);
				}
				offset += OutputPrintWidth;
		}
    
    free(buffer);
}
