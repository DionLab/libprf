#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include "pfConfig.h"
#include "pfCompute.h"
#include "pfOutput.h"

PFIMPEXP void PrintFASTA(const struct Profile * const prf, const unsigned char * const AlignedSequence,
			   const struct Alignment * const alignment, const char * const Header,
			   const size_t AlignmentLength, const float RAVE)
{
	static const char *Exfrmt[] = { ">%s revcomp\n", ">%s\n", ">%s\n" };
	
	printf((Header[0] != '>') ? Exfrmt[1+alignment->Orientation] : &Exfrmt[1+alignment->Orientation][1], Header);
	size_t i=0;

	while ((i+OutputPrintWidth) < AlignmentLength) {
		printf("%.*s\n", (int) OutputPrintWidth, &AlignedSequence[i]);
		i += OutputPrintWidth;
	}
	if (i<AlignmentLength) printf("%.*s\n", (int) (AlignmentLength) - (int) i, &AlignedSequence[i]);
// 	printf("%s\n", &SequenceText[i]);
}
