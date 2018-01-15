#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include "pfConfig.h"
#include "pfOutput.h"
#include <assert.h>

static inline void __ALWAYS_INLINE PerformRevComp(char * const restrict Sequence, const size_t SeqLength)
{
	unsigned char * restrict const CharPtr = (unsigned char*) Sequence;
	unsigned char * BackPtr = &CharPtr[SeqLength-1];
	for (size_t l=0; l<SeqLength/2; ++l) {
			const unsigned char c = CharPtr[l];
			CharPtr[l] = *BackPtr;
			*BackPtr-- = c;
	}
	if (SeqLength & 0x1) {
			const unsigned char c = CharPtr[SeqLength/2];
			CharPtr[SeqLength/2] = *BackPtr;
			*BackPtr = c;
	}
	for (size_t l=0; l<SeqLength; l++) {
		switch(Sequence[l]) {
			case 'A': Sequence[l] = 'T'; break;
			case 'C': Sequence[l] = 'G'; break;
			case 'G': Sequence[l] = 'C'; break;
			case 'T': Sequence[l] = 'A'; break;
			default:;
		}
	}
}

void HighlightedFASTAOutput(union lScores * const restrict matrix,
                            union lScores * const restrict rvmatrix,
                            const char * const SequenceText, const char * const restrict Header, 
                            const struct Profile * const prf, const Compute_t * const restrict CoreCompute,  
                            const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock)
{
	assert(matrix || rvmatrix);
	const char * const restrict ProfileSequence = prf->Sequence;
	const size_t prfLength = prf->Length;
	register const int (* restrict Scores)[4] = (const int (*)[4]) matrix;
	const size_t ld = prfLength+1;
	
	const _Bool OptimalOnly = ((struct IO_Wrapper *) Options)->OptimalOnly;
	const PrintFunction Print = ((struct IO_Wrapper *) Options)->Print;
	_Bool * const restrict HasAPath = (_Bool*) calloc(1+SeqLength,sizeof(_Bool));
	char * StateSequence = NULL;
	size_t StateLength = 0;
	
	/* What matrix to start first in case of optimal requested */
	union lScores * matrixPtr;
	Alignment_t lAlignment;
	char direction;
	int step;
	{
		int Score = NLOW, rvScore = NLOW;
		if (matrix) {
			Score = CoreCompute->GetBestScore(matrix, SeqLength, prfLength);
		}
		if (rvmatrix) {
			rvScore = CoreCompute->GetBestScore(rvmatrix, SeqLength, prfLength);
		}
		
		if (Score >= rvScore) {
			matrixPtr = matrix;
			direction = '>';
			
		}
		else {
			matrixPtr = rvmatrix;
			direction = '<';
		}
	}

	int run = 0;
	SEEK: ;
	do {
		if (matrixPtr == NULL) break;
		/* Get Best Alignment */
		
		const size_t Alignment_SeqLength = CoreCompute->GetNextBestAlignment(matrixPtr, HasAPath, &lAlignment,
																																				 SeqLength, prfLength);
		if (Alignment_SeqLength == 0UL) break;
		if (OutputVerbose) {
			fprintf(stderr, "Best alignment starts at (%i,%i)=[%i,%i] and ends at (%i,%i)=[%i,%i/%lu]"
			" with score %i and length %lu\n",
			lAlignment.Matrix.column.Begin, lAlignment.Matrix.row.Begin,
			lAlignment.Region.profile.Begin, lAlignment.Region.sequence.Begin,
			lAlignment.Matrix.column.End, lAlignment.Matrix.row.End,
			lAlignment.Region.profile.End, lAlignment.Region.sequence.End, SeqLength,
			lAlignment.Score, Alignment_SeqLength);
		}
		
		/* Prevent further alignment to use this region */
		_Bool AlreadyUsed = false;
		for(int k=lAlignment.Matrix.row.Begin; k<=lAlignment.Matrix.row.End; ++k) {
			AlreadyUsed |= HasAPath[k]; 
			HasAPath[k] = true;
		}
			
		if (!AlreadyUsed && (lAlignment.Score >= prf->CutOff || OptimalOnly) ) {
			/* Get Aligned Sequence */
			if (Alignment_SeqLength > StateLength) {
				StateLength = Alignment_SeqLength; 
				StateSequence = (char*) realloc(StateSequence, StateLength*sizeof(char));
				if (StateSequence == NULL) {
					fputs("Unable to allocate memory\n", stderr);
					exit(1);
				}
			}
 			const int counter = CoreCompute->GetStateSequence(matrixPtr, StateSequence, lAlignment.Matrix.array, ld,
																												Alignment_SeqLength);

			if (OutputVerbose) fprintf(stderr, "STATE (%lu ?= %i): %s\n", Alignment_SeqLength-1, counter, StateSequence);

			unsigned int cycle = 0;
			for (int j=0;j<counter; j++) if (StateSequence[j] == 'R') cycle++;
			
			/* Compute motif start and end */
			lAlignment.IPMB = lAlignment.Region.profile.Begin;
			lAlignment.IPME = lAlignment.Region.profile.End - (int) prfLength + 1;
			
			pthread_mutex_lock(PrintLock);
			printf(">%s direction=%c, raw_score=%i motif_start=%i motif_end=%i cycle=%u\n", Header,
						 direction, lAlignment.Score, lAlignment.IPMB, lAlignment.IPME, cycle);
			unsigned int column = 0U;
			int pos = 0;
			int start = lAlignment.Region.sequence.Begin;
			const char * restrict StateSequencePtr = StateSequence;
			int StateSequenceStep = 1;
			
			if (direction == '<') {
				start = (SeqLength-1) - lAlignment.Region.sequence.End;
				StateSequencePtr = &StateSequence[counter-1];
				StateSequenceStep = -1;
			}
	
			while (pos < start) {
				fputc(SequenceText[pos++], stdout);
				if (++column == OutputPrintWidth) {fputc('\n', stdout); column = 0U;}
			}
			printf("\e[31m\e[4m");
			int i=0;
			while(i<counter) {
				switch(StateSequence[i]) {
					case 'M': fputc(SequenceText[pos++], stdout); break;
					case 'I': fputc((SequenceText[pos++]) - 'A' + 'a', stdout); break;
					case 'D': fputs("\e[97-\e[31m", stdout); break;
					case 'R': fputs("\e[92m|\e[31m", stdout); break;
					case 'X': goto NEXT;
					default: printf("\nUNKNOWN SYMBOL %c\n", StateSequence[i]); column = 0U;
				}
				if(++column == OutputPrintWidth) {fputc('\n', stdout); column = 0U;}
				NEXT:;
				++i;
			}
			printf("\e[39m\e[24m");
			
			while (pos < SeqLength) {
				fputc(SequenceText[pos++], stdout);
				if (++column == OutputPrintWidth) {fputc('\n', stdout); column = 0U;}
			}
			fputc('\n', stdout);
			pthread_mutex_unlock(PrintLock);
		}
		if (OptimalOnly) goto FIN;
	} while (lAlignment.Score >= prf->CutOff);
	
	/* swap ordering */
	memset(HasAPath, 0, (1+SeqLength)*sizeof(_Bool));
	if (matrixPtr == matrix) {
		matrixPtr = rvmatrix;
		direction = '<';
	}
	else {
		matrixPtr = matrix;
		direction = '>';
	}
	if (!run) {
		run++;
		goto SEEK;
	}
	
	FIN:
	if (StateSequence) free(StateSequence);
}
