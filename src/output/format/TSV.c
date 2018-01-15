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

// profile_ac start strop raw_score aln_string
void PrintTSV(const struct Profile * const prf, const unsigned char * const AlignedSequence,
              const struct Alignment * const alignment, const char * const Header,
              const size_t SequenceLength, const float RAVE)
{
    const char * id;
    const char * cptr = Header;
    
    while ( *cptr != ' ' && *cptr != '\0') ++cptr;
    const int HeaderLength = (int) ((uintptr_t) cptr - (uintptr_t) Header);
    id = Header;
		if (*id == '>' ) id++;
    
    cptr = prf->AC_Number;
    while ( *cptr != ';' && *cptr != '\0' ) ++cptr;
    const int AC_Length = (int) ((uintptr_t) cptr - (uintptr_t) prf->AC_Number);
    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;
		const float normtest = (RawToNormalizedFunction == NULL) 
													? 0.0f 
													: RawToNormalizedFunction(alignment->Score, NormCoefs, RAVE, SequenceLength);   
		
		if (! prf->isCircular) {
			fprintf(stdout, "%.*s\t%.*s\t%i\t%i\t%i\t%f\t%s\n",
							AC_Length, prf->AC_Number,
							HeaderLength,id,
							alignment->Matrix.row.Begin,
							alignment->Matrix.row.End,
							alignment->Score,
							normtest,
							AlignedSequence);
		}
		else {
			fprintf(stdout, "%.*s\t%.*s\t%i\t%i\t%i\t%i\t%f\t%s\n",
							AC_Length, prf->AC_Number,
							HeaderLength,id,
							alignment->Matrix.row.Begin,
							alignment->Matrix.row.End,
							alignment->Score,
							alignment->Cycles,
							normtest,
							AlignedSequence);
		}
}
