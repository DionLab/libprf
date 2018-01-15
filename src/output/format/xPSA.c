#include "prf_config.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "pfProfile.h"
#include "pfCompute.h"
#include "pfOutput.h"

void PrintxPSA(const struct Profile * const prf, const unsigned char * const AlignedSequence,
               const struct Alignment * const alignment, const char * const Header,
               const size_t SequenceLength, const float RAVE)
{
	//PFSCANV3 NEWOUT   >fig|83333.1.peg.4317/36-219 motif=MF_00223|FolE norm_score=33.703 raw_score=4046 level_tag=! seq_end=-4 motif_start=1 motif_end=-1
	char * buffer = calloc(OutputPrintWidth+1,sizeof(char));
	char * cptr;
	
	// Extract Header
	const char * hptr = Header;
	while ( *hptr != ' ' && *hptr != '\t' && *hptr != '\0') ++hptr;
	const int HeaderLength = (int) ((uintptr_t) hptr - (uintptr_t) Header);
	
	
	// Extract ID (id)
	char * id = calloc(1+strlen(prf->Identification),sizeof(char));
	strcpy(id,prf->Identification);
	cptr = id;
	while ( *cptr != ';' && *cptr != '\0') ++cptr;
	*cptr = '\0';
		    
  char * ac  = calloc(1+strlen(prf->AC_Number),sizeof(char));
	strcpy(ac,prf->AC_Number);
	cptr = ac;
	while ( *cptr != ';' && *cptr != '\0' ) ++cptr;
	*cptr = '\0';
  RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
	const float * const restrict NormCoefs = prf->NormalizationCoefs;
	const char * des = prf->Description;
	

	const float normtest = (RawToNormalizedFunction == NULL) ? 0.0f : RawToNormalizedFunction(alignment->Score, NormCoefs, RAVE, SequenceLength);
	int level = -1;
	const static char NoTag[] = "NA";
	const static char DashTag[] = "#";
	const char * tag = DashTag;
	for ( unsigned int j = 0; j < MAXC; j++ ) { // find match level
		int icut = prf->CutOffData.Values[ j ].ICUT;
		int mcle = prf->CutOffData.Values[ j ].MCLE;
		if ( mcle == 0 && icut == 0 ) break;
		if ( alignment->Score >= icut ) { level = mcle; tag = prf->CutOffData.Values[j].CCUT; break; }
		level = mcle -1; // p.s. with -c option a match could have a score lower than the lowest defined level...
	} // p.s. prf->CutOffData.Values follows CUT_OFF line order in profile src; highest level should come first... 
	if (tag[0] == '\0') tag = DashTag;

	fprintf(stdout,	"%.*s/%i-%i motif=%s|%s norm_score=%.3f raw_score=%i level_tag=%s seq_end=%i motif_start=%i motif_end=%i\n",
					HeaderLength,
					Header,
					alignment->Matrix.row.Begin,
					alignment->Matrix.row.End,
					ac,
					id,
					normtest,
					alignment->Score,
					tag,
					(int) alignment->Matrix.row.End - (int) SequenceLength,
					alignment->IPMB,
					alignment->IPME);

	*buffer = '\0';
	size_t offset = 0;

	while (offset <= SequenceLength)
	{
		strncpy(buffer,(char*) AlignedSequence+offset,OutputPrintWidth);
		if ( '\0' != *buffer ) {
			fputs(buffer,stdout);
			fputc('\n',stdout);
		}
		offset += OutputPrintWidth;
	}
	
	free(buffer);
	free(id);
	free(ac);
}
