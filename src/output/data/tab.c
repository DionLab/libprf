#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pfCompute.h"
#include "pfOutput.h"

void TabOutput(union lScores * const restrict matrix, union lScores * const restrict rmatrix,
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
	for (int i=0; i<limit; i++) Fname[count++] =  DATOptions->BaseFileName[i];
	limit = (int) strlen(Header);
	if (limit>0) {
		Fname[count++] = '_';
		if (limit > (512-5-count)) limit = 512 - 5 - count;
		for (int i=0; i<limit; i++) Fname[count++] = (Header[i] == '/') ? '_' : Header[i];
	}
	
	Fname[count++] = '.';
	Fname[count++] = 't';
	Fname[count++] = 'a';
	Fname[count++] = 'b';
	Fname[count++] = '\0';
	
	FILE * out = fopen(Fname, "w");
	if (out) {
		fprintf(out, "#NLOW is set to %i, MLOW to %i\n#seq\tprf\tmatch\tinsertion\tdeletion\n", NLOW, NLOW/4*3);
		
		const size_t prfLength = prf->Length;
		const union lScores * restrict ptr = matrix;
		const size_t LD = GetMatrixLD(prf, CoreCompute);
		
		for (size_t iseq=0; iseq<=SeqLength; iseq++) {
			for (size_t iprf=0; iprf<=prfLength; iprf++) {
					fprintf(out, "%lu\t%lu\t%i\t%i\t%i\n", iseq, iprf,
									GetScore(ptr[iprf].Element[MATCH], CoreCompute),
									GetScore(ptr[iprf].Element[INSERTION], CoreCompute),
									GetScore(ptr[iprf].Element[DELETION], CoreCompute));
			}
			ptr += LD;
		}
		fclose(out);
	}
	else {
		fprintf(stderr, "Unable to create file %s\n", DATOptions->BaseFileName);
	}
}
