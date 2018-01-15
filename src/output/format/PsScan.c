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

void PrintPsScan(const struct Profile * const prf, const unsigned char * const AlignedSequence,
                 const struct Alignment * const alignment, const char * const Header,
                 const size_t SequenceLength, const float RAVE)
{
  char * id;
  char * buffer = calloc(OutputPrintWidth+1,sizeof(char));
  char * cptr = Header;

  RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
  const float * const restrict NormCoefs = prf->NormalizationCoefs;

	const float normtest = (RawToNormalizedFunction == NULL)
				? 0.0f
				: RawToNormalizedFunction(alignment->Score, NormCoefs, RAVE, SequenceLength);
	int level = -1;
	for ( unsigned int j = 0; j < MAXC; j++ ) { // find match level
			int icut = prf->CutOffData.Values[ j ].ICUT;
			int mcle = prf->CutOffData.Values[ j ].MCLE;
			if ( mcle == 0 && icut == 0 ) break;
			if ( alignment->Score > icut ) { level = mcle; break; }
			level = mcle -1; // p.s. with -c option a match could have a score lower than the lowest defined level...
	} // p.s. prf->CutOffData.Values follows CUT_OFF line order in profile src; highest level should come first... 

	{
		fprintf(stdout, ">%s L=%i %8.3f %6i pos. %8i -%8i [%5i, %5i] %s|%s %s\n",
			prf->Identification,
			level,
			normtest,
			alignment->Score,
			alignment->Matrix.row.Begin,
			alignment->Matrix.row.End,
			alignment->IPMB,
			alignment->IPME,
			prf->AC_Number,
			prf->Identification,
			prf->Description
		);
	}

	*buffer = '\0';
	size_t offset = 0;
	while (offset <= SequenceLength)
	{
		strncpy(buffer,AlignedSequence+offset,OutputPrintWidth);
		if ( '\0' != *buffer )
		{
			fputs(buffer,stdout);
			fputc('\n',stdout);
		}
		offset += OutputPrintWidth;
	}
  
  free(buffer);
}
