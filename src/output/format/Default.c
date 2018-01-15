#include "prf_config.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "pfProfile.h"
#include "pfOutput.h"

void PrintDefault(const struct Profile * const prf, const unsigned char * const AlignedSequence,
                  const struct Alignment * const alignment, const char * const Header,
                  const size_t SequenceLength, const float RAVE)
{  
	RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
	const float * const restrict NormCoefs = prf->NormalizationCoefs;
 
	const float normtest = (RawToNormalizedFunction == NULL) 
						? 0.0f 
						: RawToNormalizedFunction(alignment->Score, NormCoefs, RAVE, SequenceLength);
	
	if (prf->isCircular) {
		fprintf(stdout, "%7.3f %7i %7hu pos. %7i - %7i %c %s %s\n",
						normtest,
						alignment->Score,
					  alignment->Cycles,
						alignment->Matrix.row.Begin,
						alignment->Matrix.row.End,
						(alignment->Orientation < 0) ? '<' : '>',
						Header,
						AlignedSequence
		);
	}
	else {
		fprintf(stdout, "%7.3f %7i pos. %7i - %7i %c %s %s\n",
						normtest,
						alignment->Score,
						alignment->Matrix.row.Begin,
						alignment->Matrix.row.End,
						(alignment->Orientation < 0) ? '<' : '>',
						Header,
						AlignedSequence
		);
	}
}

