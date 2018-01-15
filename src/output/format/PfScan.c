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


void PrintPfscanLOpt(const struct Profile * const prf, const unsigned char * const AlignedSequence,
                     const struct Alignment * const alignment, const char * const Header,
                     const size_t SequenceLength, const float RAVE)
{ // replicates old pfscan with -l option output (e.g. L=0  32.064    9802 pos.       1 -     416 MF_00007|DNA primase DnaG [dnaG].)
  // could be used by ps_scan.pl (>=1.87)
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
				if ( alignment->Score >= icut ) { level = mcle; break; }
				level = mcle -1; // p.s. with -c option a match could have a score lower than the lowest defined level...
		} // p.s. prf->CutOffData.Values follows CUT_OFF line order in profile src; highest level should come first...

		fprintf(stdout, "L=%i %8.3f %6i pos. %8i - %7i %s|%s %s\n",
														level,
														normtest,
														alignment->Score,
														alignment->Matrix.row.Begin,
														alignment->Matrix.row.End,
														prf->AC_Number,
														prf->Identification,
														prf->Description
														);
    
}

void PrintPfscan(const struct Profile * const prf, const unsigned char * const AlignedSequence,
                 const struct Alignment * const alignment, const char * const Header,
                 const size_t SequenceLength, const float RAVE)
{ // = ~old pfscan with -lxz option output (e.g. >DNA_primase_DnaG_arc L=0   32.064   9802 pos.        1 -     416 [    1,    -1] MF_00007|DNA primase DnaG [dnaG].\nMKYLIRAR... )
    char * buffer = calloc(OutputPrintWidth+1,sizeof(char));
    char * cptr;
    
    // Extract ID (id)
    char * id = calloc(1+strlen(prf->Identification),sizeof(char));
    strcpy(id,prf->Identification);
    cptr = id;
    while ( *cptr != ';' && *cptr != '\0') ++cptr;
    *cptr = '\0';
    
    // Extract AC of matching profile
    char * ac  = calloc(1+strlen(prf->AC_Number),sizeof(char));
    strcpy(ac,prf->AC_Number);
    cptr = ac;
    while ( *cptr != ';' && *cptr != '\0' ) ++cptr;
    *cptr = '\0';
    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;
    const char * des = prf->Description;
    

		const float normtest = (RawToNormalizedFunction == NULL) 
													? 0.0f 
													: RawToNormalizedFunction(alignment->Score, NormCoefs, RAVE, SequenceLength);
		int level = -1;
		for ( unsigned int j = 0; j < MAXC; j++ ) { // find match level
				int icut = prf->CutOffData.Values[ j ].ICUT;
				int mcle = prf->CutOffData.Values[ j ].MCLE;
				if ( mcle == 0 && icut == 0 ) break;
				if ( alignment->Score >= icut ) { level = mcle; break; }
				level = mcle -1; // p.s. with -c option a match could have a score lower than the lowest defined level...
		} // p.s. prf->CutOffData.Values follows CUT_OFF line order in profile src; highest level should come first... 

		if (N > 1)
		{
				fprintf(stdout, ">%s_%i L=%i %8.3f %6i pos. %8i -%8i [%5i, %5i] %s|%s\n",
						id,    
						i+1,
						level,
						normtest,
						alignment->Score,
						alignment->Matrix.row.Begin,
						alignment->Matrix.row.End,
						alignment->IPMB,
						alignment->IPME,
						ac,
						des 
						);
		}
		else
		{
				fprintf(stdout, ">%s L=%i %8.3f %6i pos. %8i -%8i [%5i, %5i] %s|%s\n",
						id,
						level,
						normtest,
						alignment->Score,
						alignment->Matrix.row.Begin,
						alignment->Matrix.row.End,
						alignment->IPMB,
						alignment->IPME,
						ac,
						des 
						);
		}

		*buffer = '\0';
		unsigned int offset = 0;
		unsigned int len    = strlen(&AlignedSequence[i][0]);

		while (offset <= len)
		{
				strncpy(buffer,&AlignedSequence[i][0]+offset,OutputPrintWidth);
				if ( '\0' != *buffer )
				{
						fputs(buffer,stdout);
						fputc('\n',stdout);
				}
				offset += OutputPrintWidth;
		}
    
    free(buffer);
    free(id);
    free(ac);
}
