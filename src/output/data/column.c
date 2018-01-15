#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pfCompute.h"
#include "pfOutput.h"

void ColumnOutput(union lScores * const restrict matrix, union lScores * const restrict rmatrix,
               const unsigned char * const SequenceText, const char * const Header,
               const struct Profile * const prf, const Compute_t * const restrict CoreCompute,
               const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock)
{
	char Fname[512];
	const struct IO_Data * const DATOptions = Options;
	
	if (DATOptions->IsBinary) {
			fputs("Binary output no possible for tabulated format\n", stderr);
			return;
	}
	
	if (matrix == NULL) {
			fputs("TabOutput: provided matrix is null\n", stderr);
			return;
	}
	
	int count = 0;
	int limit = (int) strlen(DATOptions->BaseFileName);
	if (limit > (512-5)) limit = 512 - 5;
	for (int i=0; i<limit; i++) Fname[count++] = DATOptions->BaseFileName[i];
	limit = (int) strlen(Header);
	if (limit>0) {
		Fname[count++] = '_';
		if (limit > (512-5-count)) limit = 512 - 5 - count;
		for (int i=0; i<limit; i++) Fname[count++] = (Header[i] == '/') ? '_' : Header[i];
	}
	
	Fname[count++] = '_';
	Fname[count++] = 'c';
	const int type = count;
	Fname[count++] = '?';
	Fname[count++] = '.';
	Fname[count++] = 'c';
	Fname[count++] = 'o';
	Fname[count++] = 'l';
	Fname[count++] = '\0';
	
	static const unsigned int Index[3] = { MATCH, INSERTION, DELETION };
	static const char AlphaType[3] = { 'm', 'i', 'd' };
	const int MINIMUM = NLOW;
	const size_t LD = GetMatrixLD(prf, CoreCompute);
	const size_t prfLength = prf->Length;
	
	for (int i=0; i<3; i++) {
		Fname[type] = AlphaType[i];
		FILE * out = fopen(Fname, "w");
		const int index = Index[i];
		if (out) {
			for (size_t iprf=0; iprf<=prfLength; iprf++) {
				const union lScores * restrict ptr = &matrix[iprf];
				for (size_t iseq=0; iseq<SeqLength; iseq++) {
					const int value = GetScore(ptr->Element[index], CoreCompute); 
					fprintf(out, "%i,", (value <= MINIMUM) ? -32767 : value);
					ptr += LD;
				}
				const int value = GetScore(ptr->Element[index], CoreCompute); 
				fprintf(out, "%i\n", (value <= MINIMUM) ? -32767 : value);
			}
			fclose(out);
		}
		else {
			fprintf(stderr, "Unable to create file %s\n", Fname);
		}
	}
}

