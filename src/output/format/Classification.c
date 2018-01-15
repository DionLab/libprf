#include "prf_config.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "pfProfile.h"
#define __USE_INLINE_FUNCTIONS__
#include "pfSequence.h"
#include "pfOutput.h"

/*
 * WARNING: alignement does not start at 0 but 1 !!!
 *          NEEDS TO BE FIXED SOME DAY
 */

static int RescoreAlignment(const struct Profile * const prf, const unsigned char * const AlignedSequence,
		     const struct Alignment * const alignment, const size_t SequenceLength)
{
  int Score = 0;
  int iprf=1;
  enum VectorPosition PreviousState;
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  // Prologue
  //////////////////////////////////////////////////////////////////////////////////////////////
  if (alignment->IPMB != 1 && alignment->IPME != -1) {
      fputs("Currently we do not provide classification for other than semiglobal alignment.\n", stderr);
      exit(1);
  }
  
  const size_t AlignmentLength            = strlen(AlignedSequence);
  const char * restrict Seq               = AlignedSequence;
  const char * restrict const MaxSequence = &AlignedSequence[AlignmentLength];
  unsigned int SequenceIndex              = alignment->Matrix.row.Begin;
  const unsigned int SequenceBegin        = alignment->Matrix.row.Begin;
  const unsigned int SequenceEnd          = alignment->Matrix.row.End;
  
  /* Get first alignment state */
  if (*Seq == '-') 
    PreviousState = DELETION;
  else
    PreviousState = ( *Seq < 'a' ) ? MATCH : INSERTION;
  
  const TransitionScores * restrict pTransitions = prf->Scores.Insertion.Transitions;
  const StoredIntegerFormat * pMatch             = prf->Scores.Match.Alphabet;
  const StoredIntegerFormat * pInsertion         = prf->Scores.Insertion.Alphabet;
  const size_t AlignStep                         = prf->Scores.Match.AlignStep;
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  // First Sequence Line
  ////////////////////////////////////////////////////////////////////////////////////////////// 
  if (SequenceBegin == 1) {
    register const ScoreTuple * restrict FirstSequenceProtein = prf->Scores.Insertion.FirstSequenceProtein;
    int lScore = (int) FirstSequenceProtein->To[PreviousState];
//     printf("Alignment starts at the beginning of the sequence\n");    
    /* We need to further keep looking that we are not in the case of a multiple deletion entrance */
    if (PreviousState == DELETION) {
      pInsertion += AlignStep;
      ++pTransitions;
      ++iprf;
      while (*(++Seq) == '-') {
// 	fputc(Seq[-1], stdout);
	lScore += ((int) pMatch[_D]) + ((int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(DELETION, DELETION)]);
	++pTransitions;
	pInsertion += AlignStep;
	pMatch     += AlignStep;
	++iprf;
      }
    }
    else if (PreviousState == MATCH) {
//        pInsertion += AlignStep;
//        ++pTransitions;
    }
    Score += lScore;
//     printf(" Profile pos: %4i Score: %i\n", iprf, Score);
  }
  else {
    Score += (int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(EXTRA, PreviousState)];
//     ++pTransitions;
//     pInsertion += AlignStep;
//     printf(" Profile pos: %4i Score: %i\n", iprf, Score);
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  // Loop through the internal sequence indices
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  if (Seq < MaxSequence) { 
    do {
//       fputc(*Seq, stdout);
      if (*Seq == '-') {
	 Score += (int) pMatch[_D] + (int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(PreviousState, DELETION)]; ;
	 PreviousState = DELETION;
	 ++pTransitions;
	 pInsertion += AlignStep;
	 pMatch     += AlignStep;
	 ++iprf;
      }
      else {
	const size_t index = (size_t) TranslateCharToIndex(*Seq, prf->Alphabet_Mapping);
	if ( *Seq < 'a' ) {
	  Score   += (int) pMatch[index] + (int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(PreviousState, MATCH)];
	  PreviousState = MATCH;
	  ++pTransitions;
	  pInsertion += AlignStep;
	  pMatch     += AlignStep;
	  ++iprf;
	}
	else {
	  Score    += (int) pInsertion[index] + (int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(PreviousState, INSERTION)];
// 	  printf(" Profile pos: %4i ?I: %i i: %i\n", iprf, pTransitions->Element[TRANSITION_INDEX_FROM_TO(PreviousState, INSERTION)], pInsertion[index]);
	  PreviousState  = INSERTION;
	}
	++SequenceIndex;
      }
//       printf(" Profile pos: %4i Score: %i\n", iprf, Score);
    } while ( ++Seq < MaxSequence && SequenceIndex < SequenceLength);
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  // Epilogue
  //////////////////////////////////////////////////////////////////////////////////////////////
  if (SequenceEnd == SequenceLength) {
    fputs(": End on sequence end", stdout);
  }
  else {
//     fputc(*Seq, stdout);
    if (*Seq == '-') {
      Score += (int) pMatch[_D] + (int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(PreviousState, EXTRA)]; ;
      PreviousState = DELETION;
      ++pTransitions;
      pInsertion += AlignStep;
      pMatch     += AlignStep;
    }
    else {
      const size_t index = (size_t) TranslateCharToIndex(*Seq, prf->Alphabet_Mapping);
      if ( *Seq < 'a' ) {
	Score   += (int) pMatch[index] + (int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(PreviousState, EXTRA)];
	PreviousState = MATCH;
	++pTransitions;
	pInsertion += AlignStep;
	pMatch     += AlignStep;
      }
      else {
	Score    += (int) pInsertion[index] + (int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(PreviousState, EXTRA)];
	PreviousState  = INSERTION;
      }
    }
  }
  return Score;
}

void PrintClassification(const struct Profile * const prf, const unsigned char * const AlignedSequence,
                         const struct Alignment * const alignment, const char * const Header,
                         const size_t SequenceLength, const float RAVE)
{
  RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
  const float * const restrict NormCoefs = prf->NormalizationCoefs;
  const float normtest = (RawToNormalizedFunction == NULL)
                         ? 0.0f 
                         : RawToNormalizedFunction(alignment->Score, NormCoefs, RAVE, SequenceLength);   
  fprintf(stdout, "%s %i %f\n%s\n",
	    Header,
	    alignment->Score,
	    normtest,
	    AlignedSequence);
  struct Profile * BestPrf = NULL;
  register struct Profile * tmpPrf = (struct Profile *) prf;
  int BestScore = NLOW;
  
  while (tmpPrf->next) {
		tmpPrf = tmpPrf->next;
		const int Score = RescoreAlignment(tmpPrf, AlignedSequence, alignment, SequenceLength);
		fprintf(stdout,"%s : score %i\n", tmpPrf->Description, Score);
		if (BestScore < Score) {
			BestScore = Score;
			BestPrf = tmpPrf;
		}
  }
  
  fprintf(stdout, "Best : %s : score %i\n", BestPrf->Description, BestScore);
  
}
