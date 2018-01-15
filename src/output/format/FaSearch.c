#include "prf_config.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "pfProfile.h"
#include "pfOutput.h"
#include "dna.h"

/*
 * WARNING: alignement does not start at 0 but 1 !!!
 *          NEEDS TO BE FIXED SOME DAY
 */

#define OUTPUT_WIDTH 80
void Printfasearch(const struct Profile * const prf, const unsigned char * const AlignedSequence,
                   const struct Alignment * const alignment, const char * const Header,
                   const size_t SequenceLength, const float RAVE)
{
	char AminoBuffer[OUTPUT_WIDTH];
	char DNABuffer[OUTPUT_WIDTH];
	char StateBuffer[OUTPUT_WIDTH];
	char ConcensusBuffer[OUTPUT_WIDTH];
	
	/* Get Sequence identification from header*/
	const char * SeqIDPtr = Header + 1;
	const char * const MaxLengthPtr = Header + strlen(Header);
	while (*SeqIDPtr == ' ' && SeqIDPtr < MaxLengthPtr)  SeqIDPtr++;
	int c = 0;
	while ((SeqIDPtr[c] != ' ' && SeqIDPtr[c] != '\n') && &SeqIDPtr[c] < MaxLengthPtr) c++;
	
	RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
	const float * const restrict NormCoefs = prf->NormalizationCoefs;

	/* Print header */
	const int Score = alignment->Score;
	const float BestNormalizedScore = RawToNormalizedFunction(Score, NormCoefs, RAVE, SequenceLength);
	printf("ID\t%u\nIR\t%s\nIS\t%.*s\nSR\t%i\nSN\t%lf\nCB\t%i %i %lu\nCS\t%i %i %lu\n**",
					1,//ID
					prf->Identification, //IR
					//prf->Identification, //IB
					c, SeqIDPtr, //IS
					Score, //SR
					BestNormalizedScore, //SN
					alignment->IPMB, (int)prf->Length + 1 + alignment->IPME, prf->Length, //CB
					alignment->Matrix.row.Begin, alignment->Matrix.row.End, SequenceLength //CS
				);

	if (prf->Sequence == NULL) {
		fputs("No concensus with red profile !!!\n", stderr);
		exit(1);
	}
	const int ConcensusStart = alignment->IPMB - 1;
	const char * restrict Concensus = &(prf->Sequence[ConcensusStart]);
	const int counter = alignment->Matrix.row.End - alignment->Matrix.row.Begin + 1;
	int i = 0;
	int k = 0;
	int l = 0;
	while (i < counter) {
		char State;
		int index;
		if (AlignedSequence[i] >= 'a') {
			State = 'I';
			index = (int) AlignedSequence[i] - (int) 'a';
		}
		else if (AlignedSequence[i] >= 'A') {
			State = 'M';
			index = (int) AlignedSequence[i] - (int) 'A';
		}
		else {
			State = AlignedSequence[i];
		}
		switch (State) {
			case 'M':
				DNABuffer[k]     = AminoToIUPAC[index][0];
				DNABuffer[k+1]   = AminoToIUPAC[index][1];
				DNABuffer[k+2]   = AminoToIUPAC[index][2];
				StateBuffer[k]   = 'M';
				StateBuffer[k+1] = 'M';
				StateBuffer[k+2] = 'M';
				AminoBuffer[k]   = '<';
				AminoBuffer[k+1] = AlignedSequence[i];
				AminoBuffer[k+2] = '>';
				ConcensusBuffer[k]   = '<';
				ConcensusBuffer[k+1] = Concensus[l];
				ConcensusBuffer[k+2] = '>';
				k+=3;
				l+=1;
				break;
			case 'I':
				DNABuffer[k]     = (char) ((unsigned char) AminoToIUPAC[index][0] - 'A' + 'a');
				DNABuffer[k+1]   = (char) ((unsigned char) AminoToIUPAC[index][1] - 'A' + 'a');
				DNABuffer[k+2]   = (char) ((unsigned char) AminoToIUPAC[index][2] - 'A' + 'a');
				StateBuffer[k]   = 'I';
				StateBuffer[k+1] = 'I';
				StateBuffer[k+2] = 'I';
				AminoBuffer[k]   = '<';
				AminoBuffer[k+1] = AlignedSequence[i];
				AminoBuffer[k+2] = '>';
				ConcensusBuffer[k]   = ' ';
				ConcensusBuffer[k+1] = ' ';
				ConcensusBuffer[k+2] = ' ';
				k+=3;
				break;
			case '-':
				DNABuffer[k  ]   = '-';
				DNABuffer[k+1]   = '-';
				DNABuffer[k+2]   = '-';
				StateBuffer[k]   = 'D';
				StateBuffer[k+1] = 'D';
				StateBuffer[k+2] = 'D';
				AminoBuffer[k]   = '<';
				AminoBuffer[k+1] = '-';
				AminoBuffer[k+2] = '>';
				ConcensusBuffer[k]   = '<';
				ConcensusBuffer[k+1] = Concensus[l];
				ConcensusBuffer[k+2] = '>';
				k+=3;
				l+=1;
				break;
			default:
				fprintf(stderr, "Unknown state '%c'=(%i) @ %i\n%s\n", State, (int) State, counter+i, AlignedSequence);
				exit(1);
		}	      
		i++;
		if ( k > OUTPUT_WIDTH - 4 ) {
			AminoBuffer[k] = '\0';
			DNABuffer[k]   = '\0';
			StateBuffer[k] = '\0';
			ConcensusBuffer[k] = '\0';
			printf("\nEX\t%s\nAS\t%s\nSA\t%s\nSD\t%s\n", ConcensusBuffer, StateBuffer, AminoBuffer, DNABuffer);
			k = 0;
		}
	}
	if (k>0) {
		AminoBuffer[k] = '\0';
		DNABuffer[k]   = '\0';
		StateBuffer[k] = '\0';
		ConcensusBuffer[k] = '\0';
		printf("\nEX\t%s\nAS\t%s\nSA\t%s\nSD\t%s\n", ConcensusBuffer, StateBuffer, AminoBuffer, DNABuffer);
			
	}
	/* Print footer */
	fputs("**\n//\n", stdout);
			
	
}
	
