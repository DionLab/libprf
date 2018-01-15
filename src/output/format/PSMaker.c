#include "prf_config.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "pfProfile.h"
#include "pfCompute.h"
#include "pfOutput.h"

void PrintPSMaker(const struct Profile * const prf, const unsigned char * const AlignedSequence,
                  const struct Alignment * const alignment, const char * const Header,
                  const size_t SequenceLength, const float RAVE)
{ // replicates old pfsearch -x output
    char * buffer = calloc( OutputPrintWidth+1, sizeof( char ) );
    char * header = calloc( strlen( Header ) + 1, sizeof( char ) );
    strcpy( header, Header );

    char * cptr = header;
    while ( *cptr != ' ' && *cptr != '\0') ++cptr;
    *cptr = '\0';
    while ( *cptr != '|' && *cptr != '>' ) --cptr;
    ++cptr;

    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;
		const float norm_score = (RawToNormalizedFunction == NULL)
												? 0.0f
												: RawToNormalizedFunction(alignment->Score, NormCoefs, RAVE, SequenceLength);
				fprintf( stdout, ">%s %8.3f %6i pos. %8i -%8i %s\n",
						cptr,
						norm_score,
						alignment->Score,
						alignment->Matrix.row.Begin,
						alignment->Matrix.row.End,
						header+1
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
    free(header);
}
