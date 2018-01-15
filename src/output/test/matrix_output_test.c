/*******************************************************
 *                        PFTOOLS
 *******************************************************
 *  Mar 17, 2016 matrix_output_test.c
 *******************************************************
 * (C) 2012-2016 Swiss Institute of Bioinformatics
 *     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include "pfCompute.h"
#include "pfOutput.h"
//#include "matrix_functions.h"

void TestOutput(union lScores * const restrict matrix,
								union lScores * const restrict rmatrix,
								const char * const SequenceText, const char * const Header,
								const struct Profile * const prf, const Compute_t * const restrict CoreCompute, 
								const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock)
{
	char AlignmentSequence[4096] __attribute__((aligned(16)));
	Alignment_t Alignment;
	const size_t prfLength = prf->Length;	
	const size_t Length = CoreCompute->GetBestAlignment(matrix, &Alignment, SeqLength, prfLength);
	if (Length < 4096) {
		CoreCompute->GetAlignmentSequence(matrix, (const unsigned char*)SequenceText, AlignmentSequence, 
												 &Alignment, SeqLength, prfLength);
		if (strncmp(AlignmentSequence, &Header[1], strlen(&Header[1])) == 0) {
			printf("Testing for '%s'\t'%s'\tOK\n", &Header[1], AlignmentSequence);
		}
		else {
			printf("Testing for '%s'\t'%s'\tFAILED\n", &Header[1], AlignmentSequence);

			// Generate a pdf of the error
// #ifdef PRF_OUTPUT_PDF
// 				PlotAllSequence = 1;
// 			char FileName[256];
// 			snprintf(FileName, 255, "Test_%s", &Header[1]);
// 			extern struct IO_PDF OptionsPDF;
// 			OptionsPDF.BaseFileName = FileName;
// 			PDFOutput(matrix, rmatrix, SequenceText, Header, prf, CoreCompute, SeqLength, (void*) &OptionsPDF, PrintLock);
// #endif
		}
	}
	else {
		fprintf(stderr, "Allocated buffer (4096) is not large enough to contain Alignment sequence (%lu)\n", Length);
	}

	
	

}

