#include "prf_config.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "pfProfile.h"
#include "pfCompute.h"
#include "pfOutput.h"

void PrintSimple(const struct Profile * const prf, const unsigned char * const AlignedSequence,
                 const struct Alignment * const alignment, const char * const Header,
                 const size_t SequenceLength, const float RAVE)
{
  const char * cptr = Header;
  while ( *cptr != ' ' && *cptr != '\0') ++cptr;
  int HeaderLength = (intptr_t) cptr - (intptr_t) Header;
  RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
  const float * const restrict NormCoefs = prf->NormalizationCoefs;  
	const float normtest = (RawToNormalizedFunction == NULL) 
													? 0.0f 
													: RawToNormalizedFunction(alignment->Score, NormCoefs, RAVE, SequenceLength);   
	fprintf(stdout, "%.*s %i %f\n%s\n",
		HeaderLength,
		Header,
		alignment->Score,
		normtest,
		AlignedSequence);
  
}
