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
void PrintOneLine(const struct Profile * const prf, const unsigned char * const AlignedSequence,
                  const struct Alignment * const alignment, const char * const Header,
                  const size_t SequenceLength, const float RAVE)
{
	const float * const restrict NormCoefs = prf->NormalizationCoefs;
	
	/* Get Sequence identification from header*/
	const char * SeqIDPtr = &Header[1];
	const char * const MaxLengthPtr = &Header[strlen(Header)];
	while (*SeqIDPtr == ' ' && SeqIDPtr < MaxLengthPtr)  SeqIDPtr++;
	int c = 0;
	while ((SeqIDPtr[c] != ' ' && SeqIDPtr[c] != '\n') && &SeqIDPtr[c] < MaxLengthPtr) c++;
	
	const int Score = alignment->Score;
	const float BestNormalizedScore = prf->RawToNormalized(Score, NormCoefs, 0.0f, SequenceLength);
	
	char Orientation;
	union URegion Profile = { alignment->Matrix.column.Begin,  alignment->Matrix.column.End };
	union URegion Sequence;
	if (prf->ReverseSequence) {
		Orientation = '-';
		Sequence.Begin = (int) SequenceLength - alignment->Matrix.row.Begin;
		Sequence.End   = (int) SequenceLength - alignment->Matrix.row.End;
	}
	else{
		Orientation = '+';
		Sequence.Begin = alignment->Matrix.row.Begin;
		Sequence.End   = alignment->Matrix.row.End;
	}
	
	/*
		* [PROFILE NAME] [TARGET NAME] [DIRECTION] [RAW SCORE] [NORMALIZED SCORE]
		* [PROFILE MATCH START] [PROFILE MATCH END] [PROFILE LENGTH]
		* [TARGET MATCH START] [TARGET MATCH STOP] [TARGET LENGTH]
		*/
	printf("%s\t%.*s\t%c\t%i\t%lf\t%i\t%i\t%lu\t%i\t%i\t%lu\n",
				prf->Identification, 		//IR
				c, SeqIDPtr, 						//IS
				Orientation ,						//orientation
				Score, 									//SR
				BestNormalizedScore, 		//SN
				Profile.Begin, Profile.End, prf->Length, 	//CB
				Sequence.Begin, Sequence.End, SequenceLength 	//CS 
				);
}
