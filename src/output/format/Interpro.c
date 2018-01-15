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

void PrintInterPro(const struct Profile * const prf, const unsigned char * const AlignedSequence,
                   const struct Alignment * const alignment, const char * const Header,
                   const size_t SequenceLength, const float RAVE)
{ // replicates old pfsearch -lxz output
    const char * id;
    char * des;
    char * buffer = calloc(OutputPrintWidth+1,sizeof(char));
    char * cptr = Header;
    
    // Extract ID (id), short ID (cptr) and description (des)
    while ( *cptr != ' ' && *cptr != '\0') ++cptr;
    des = ( *cptr == '\0' ) ? cptr : cptr + 1;
    *cptr = '\0';
    while ( *cptr != '|' && *cptr != '>' ) --cptr;
    ++cptr; 
    id = Header;
    ++id;
    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;

    for (unsigned int i=0; i<N; ++i) {   
        const float normtest = (RawToNormalizedFunction == NULL) 
	                       ? 0.0f 
	                       : RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
        int level = -1;
        for ( unsigned int j = 0; j < MAXC; j++ ) { // find match level
            int icut = prf->CutOffData.Values[ j ].ICUT;
            int mcle = prf->CutOffData.Values[ j ].MCLE;
            if ( mcle == 0 && icut == 0 ) break;	
            if ( alignment[i].Score >= icut ) { level = mcle; break; }
            level = mcle -1; // p.s. with -c option a match could have a score lower than the lowest defined level...
        } // p.s. prf->CutOffData.Values follows CUT_OFF line order in profile src; highest level should come first... 

        if (N > 1)
        {
            fprintf(stdout, ">%s_%i L=%i %8.3f %6i pos. %8i -%8i [%5i, %5i] %s %s\n",
                cptr,
                i+1,
                level,
                normtest,
                alignment[i].Score,
                alignment[i].Matrix.row.Begin,
                alignment[i].Matrix.row.End,
                alignment[i].IPMB,
                alignment[i].IPME,
                id,
                des 
               );
        }
        else
        {
            fprintf(stdout, ">%s L=%i %8.3f %6i pos. %8i -%8i [%5i, %5i] %s %s\n",
                cptr,
                level,
                normtest,
                alignment[i].Score,
                alignment[i].Matrix.row.Begin,
                alignment[i].Matrix.row.End,
                alignment[i].IPMB,
                alignment[i].IPME,
                id,
                des 
               );
        }

        *buffer = '\0';
        unsigned int offset = 0;
        unsigned int len    = strlen(&AlignedSequence[i][0]);

        while ( offset <= len ) {
            strncpy(buffer,&AlignedSequence[i][0]+offset,OutputPrintWidth);
            if ( '\0' != *buffer ) {
                fputs(buffer,stdout);
                fputc('\n',stdout);
            }
            offset += OutputPrintWidth;
        }
    }
    free(buffer);
}
