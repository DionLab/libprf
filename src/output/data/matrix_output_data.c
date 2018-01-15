/*******************************************************
 *                        PFTOOLS
 *******************************************************
 *  Nov 3, 2015 pfplot_output_data.c
 *******************************************************
 * (C) 2012-2015 Swiss Institute of Bioinformatics
 *     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pfCompute.h"
#include "pfOutput.h"

void DataOutput(union lScores * const restrict matrix, union lScores * const restrict rmatrix,
                const unsigned char * const SequenceText, const char * const Header,
                const struct Profile * const prf, const Compute_t * const restrict CoreCompute,
                const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock)
{
	char OutputFileName[256] __attribute__((aligned(16)));
	unsigned char AlignmentSequence[4096] __attribute__((aligned(16)));
	char ComeFrom[4];
	Alignment_t Alignment;
	const size_t prfLength = prf->Length;
	
// 	ComeFrom[PRIORITY_MATCH]     = 'M';
// 	ComeFrom[PRIORITY_INSERTION] = 'I';
// 	ComeFrom[PRIORITY_DELETION]  = 'D';
// 	ComeFrom[PRIORITY_EXTRA]     = 'X';
		
	const size_t Length = CoreCompute->GetBestAlignment(matrix, &Alignment, SeqLength, prfLength);
	if (Length < 4096) {
		CoreCompute->GetAlignmentSequence(matrix, SequenceText, AlignmentSequence, 
						  &Alignment, SeqLength, prfLength);
		fprintf(stderr, "%i\t%i\t%lu\t%s\n", Alignment.Score, Alignment.Cycles, strlen((char*) AlignmentSequence), AlignmentSequence);
	}
	else {
		fprintf(stderr, "Allocated buffer (4096) is not large enough to contain Alignment sequence (%lu)\n", Length);
	}
	
	
	
// 	snprintf(OutputFileName, 255, "%s_%lu.dat", BaseFileName, SeqID);
// 	FILE* output = fopen(OutputFileName, "w");
// 	if ( output == NULL ) {
// 		fputs("Unable to open output file.", stderr);
// 		exit(1);
// 	}
// 	const union lScores * restrict pmatrix = matrix;
// 	for (size_t iseq=0; iseq<=SeqLength; ++iseq) {
// 		for (size_t iprf=0; iprf<=prfLength; ++iprf) {
// #ifdef TAG
// 			fprintf(output, "%5lu %5lu %9i %c %c %9i %c %c %9i %c %c %9i %c %c\n", iseq, iprf,
// 							TO_SCORE(pmatrix[iprf].Element[MATCH]),     ComeFrom[GET_STATE(pmatrix[iprf].Element[MATCH])], pmatrix[iprf].Element[MATCH] & CYCLE_MASK ? 'C' : '.',
// 							TO_SCORE(pmatrix[iprf].Element[INSERTION]), ComeFrom[GET_STATE(pmatrix[iprf].Element[INSERTION])],pmatrix[iprf].Element[INSERTION] & CYCLE_MASK ? 'C' : '.',
// 							TO_SCORE(pmatrix[iprf].Element[DELETION]),  ComeFrom[GET_STATE(pmatrix[iprf].Element[DELETION])],pmatrix[iprf].Element[DELETION] & CYCLE_MASK ? 'C' : '.',
// 							TO_SCORE(pmatrix[iprf].Element[DUMMY]),     ComeFrom[GET_STATE(pmatrix[iprf].Element[DUMMY])],pmatrix[iprf].Element[DUMMY] & CYCLE_MASK ? 'C' : '.');
// #else
// 			fprintf(output, "%5lu %5lu %7i %7i %7i %7i\n", iseq, iprf,
// 							TO_SCORE(pmatrix[iprf].Element[MATCH]),
// 							TO_SCORE(pmatrix[iprf].Element[INSERTION]),
// 							TO_SCORE(pmatrix[iprf].Element[DELETION]),
// 							TO_SCORE(pmatrix[iprf].Element[DUMMY]));
// #endif
// 		}
// 		pmatrix += prfLength+1;
// 	}
// 	
// 	fclose(output);
// 	if (rmatrix) {
// 		snprintf(OutputFileName, 255, "%s_%lu_reverse.dat", BaseFileName, SeqID);
// 		output = fopen(OutputFileName, "w");
// 		if ( output == NULL ) {
// 			fputs("Unable to open output file.", stderr);
// 			exit(1);
// 		}
// 		pmatrix = rmatrix;
// 		for (size_t iseq=0; iseq<=SeqLength; ++iseq) {
// 			for (size_t iprf=0; iprf<=prfLength; ++iprf) {
// #ifdef TAG
// 				fprintf(output, "%5lu %5lu %9i %c %c %9i %c %c %9i %c %c %9i %c %c\n", iseq, iprf,
// 								TO_SCORE(pmatrix[iprf].Element[MATCH]),     ComeFrom[GET_STATE(pmatrix[iprf].Element[MATCH])],pmatrix[iprf].Element[MATCH] & CYCLE_MASK ? 'C' : '.',
// 								TO_SCORE(pmatrix[iprf].Element[INSERTION]), ComeFrom[GET_STATE(pmatrix[iprf].Element[INSERTION])],pmatrix[iprf].Element[INSERTION] & CYCLE_MASK ? 'C' : '.',
// 								TO_SCORE(pmatrix[iprf].Element[DELETION]),  ComeFrom[GET_STATE(pmatrix[iprf].Element[DELETION])],pmatrix[iprf].Element[DELETION] & CYCLE_MASK ? 'C' : '.',
// 								TO_SCORE(pmatrix[iprf].Element[DUMMY]),     ComeFrom[GET_STATE(pmatrix[iprf].Element[DUMMY])], pmatrix[iprf].Element[DUMMY] & CYCLE_MASK ? 'C' : '.');
// #else
// 				fprintf(output, "%5lu %5lu %7i %7i %7i %7i\n", iseq, iprf,
// 								TO_SCORE(pmatrix[iprf].Element[MATCH]),
// 								TO_SCORE(pmatrix[iprf].Element[INSERTION]),
// 								TO_SCORE(pmatrix[iprf].Element[DELETION]),
// 								TO_SCORE(pmatrix[iprf].Element[DUMMY]));
// #endif
// 			}
// 			pmatrix += prfLength+1;
// 		}
// 		fclose(output);
// 	}
}

