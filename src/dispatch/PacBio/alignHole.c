#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <assert.h>
#include "pfProfile.h"
#include "pfSequence.h"
#include "pfOutput.h"
#include "pfCompute.h"
#include "pfDispatch.h"
#include "pb_hdf5.h"
#include "matrix_functions.h"

int alignHoleReads(const PacBio_t * const restrict PB, const struct Profile * const restrict prf,
                   const Compute_t * const restrict CoreCompute,
                   const OutputMethod Output,
                   const void * const restrict OutputMethodOptions,
                   const char * const restrict OutputFileName,
                   const PacBioDispatchOptions_t * const restrict Options)
{

	if (Options->Selection & PB_DISPATCH_CONSENSUS) {
		fputs("Single hole analysis on consensus not implemented yet!\n", stderr);
		return ERROR;
	}
	
	if (!(PB->Content & BASECALLING)) return ERROR; 
	HoleData_t HReg = HOLE_DATA_INIT;
	const unsigned int HoleNumber = Options->ZMW[0];
	int res = getHoleRegions(PB, &HReg, HoleNumber);
	if (res < 0) return res;
	
	/* Get the HQRegion from HReg starting for the end to speed things up */
	unsigned int count = HReg.nRegions;
	int HQStart, HQStop;
	int LargestSequence = 0;
	{
		int HQindex = 0;
		_Bool HasValidHQ = true;
		HoleRegion_t * const HR = HReg.Regions;
		/* WARNING: That assumes inserts are always prior to HQregion !! */
		for (int r=0; r<count; r++) {
			//fprintf(stderr, "%i : %u - %u\n", HR[r].type, HR[r].start, HR[r].stop);
			if (HR[r].type == Insert) {
				const int lstop = HR[r].stop;
				LargestSequence = (LargestSequence > lstop) ? LargestSequence : lstop;
			}
			else if (HR[r].type == HQRegion) {
				if (HR[r].stop == 0) {
					HasValidHQ = false;
				}
				HQindex = r;
			}
		}
		
		HQStart = HR[HQindex].start;
		HQStop  = HR[HQindex].stop;
		
		if (HR[HQindex].quality < PacBioHQThreshold) {
			HQStop = 0;
		}
		
		if (((Options->Selection & PB_KEEP_INVALID) && !HasValidHQ) || (Options->Selection & PB_DISCARD_FILTER)) {
			HQStart = 0;
			HQStop  = LargestSequence;
		}
	}
	if (OutputVerbose) fprintf(stderr, "Hole %u HQ region set to [%i,%i]\n", HoleNumber, HQStart, HQStop);
	
	if (!(Options->Selection & PB_DISPATCH_SUBREADS)) {
		HReg.Regions[0].type  = Insert;
		HReg.Regions[0].start = HQStart;
		HReg.Regions[0].stop  = HQStop;
		HReg.Regions[1].type  = HQRegion;
		HReg.Regions[1].start = HQStart;
		HReg.Regions[1].stop  = HQStop;
		if (OutputVerbose) fprintf(stderr, "alignHoleReads treating hole %u as a stream of length %i, region [%i:%i]\n",
		                           HoleNumber, HQStop, HQStart, HQStop);
		count = HReg.nRegions = 2;
	}
	else if (!(Options->Selection & PB_DISCARD_FILTER)){
		LargestSequence = 0;
		HoleRegion_t * restrict Regions = HReg.Regions;
		for (unsigned int r=0; r<count; r++) {
			if (Regions[r].type == Insert) {
				if (Regions[r].start < HQStart) Regions[r].start = HQStart;
				if (Regions[r].stop > HQStop) Regions[r].stop = HQStop;
				const int SeqSize = Regions[r].stop - Regions[r].start;
				if (SeqSize > 0) 
					LargestSequence = (LargestSequence < SeqSize) ? SeqSize : LargestSequence;
				else 
					Regions[r].type = Useless;
			}
		}
	}
	LargestSequence++;
	
	HoleReads_t HReads = HOLE_READS_INIT;
	res = getHoleReads_HR(PB, &HReads, &HReg);
	if (res < 0) return res;
	
	char Header[48];
	unsigned char RevComp_Alphabet_Mapping[ALPHABET_SIZE+2] __attribute__((aligned(16)));

	const size_t prfLength = prf->Length;
	const size_t WorkSize  = prf->Length + 1;
	
	if (OutputVerbose) fprintf(stderr, "alignHoleReads on hole %u with cutoff set to %i\n", HoleNumber, prf->CutOff);
		
	/*************************************************************************/
  /*                          ALLOCATE MEMORY                              */
  /*************************************************************************/
	
	union lScores * restrict matrix   = _mm_malloc(WorkSize*LargestSequence*sizeof(union lScores), 64);
	union lScores * restrict rvmatrix = _mm_malloc(WorkSize*LargestSequence*sizeof(union lScores), 64);
	int  * restrict WORK              = _mm_malloc(2*WorkSize*sizeof(union lScores)+63,64);
	unsigned char * restrict SequenceIndex = (unsigned char*) malloc(LargestSequence*sizeof(unsigned char));
	if ( rvmatrix == NULL || matrix == NULL || WORK == NULL || SequenceIndex == NULL ) {
		res = -1; 
		goto FIN;
	}
	
	PFSequence PFSeq;
	PFSeq.ProfileIndex = SequenceIndex;

	/* Builds RevComp index */
	memcpy(RevComp_Alphabet_Mapping, prf->Alphabet_Mapping, (ALPHABET_SIZE+2)*sizeof(unsigned char));
	RevComp_Alphabet_Mapping[(int) 'A' - (int) 'A'] =  prf->Alphabet_Mapping[(int) 'T' - (int) 'A' ];
	RevComp_Alphabet_Mapping[(int) 'C' - (int) 'A'] =  prf->Alphabet_Mapping[(int) 'G' - (int) 'A' ];
	RevComp_Alphabet_Mapping[(int) 'G' - (int) 'A'] =  prf->Alphabet_Mapping[(int) 'C' - (int) 'A' ];
	RevComp_Alphabet_Mapping[(int) 'T' - (int) 'A'] =  prf->Alphabet_Mapping[(int) 'A' - (int) 'A' ];

	pthread_mutex_t PrintLock = PTHREAD_MUTEX_INITIALIZER;
	for (unsigned int r=0; r<count; r++) {
		if (HReg.Regions[r].type == Insert) {
			const size_t SeqLength = HReg.Regions[r].stop - HReg.Regions[r].start;
			char * const restrict seq =  &(HReads.Stream[HReg.Regions[r].start]);
			
			/* Copy sequence to local space index */ 
			memcpy(PFSeq.ProfileIndex, seq, SeqLength);
			PFSeq.Length = SeqLength;
			
			/* Translate into indices */
			TranslateSequenceToIndex(&PFSeq, prf->Alphabet_Mapping);
			
			/* Build matrix */
			CoreCompute->BuildMatrix(prf, PFSeq.ProfileIndex, matrix, WORK, NULL, 0, SeqLength);
			
			/* Copy sequence to local space index */ 
			memcpy(PFSeq.ProfileIndex, seq, SeqLength);
			
			/* Translate into indices */
			TranslateSequenceToIndex(&PFSeq, RevComp_Alphabet_Mapping);
			ReverseTranslatedSequence(&PFSeq);
			
			/* Build reverse matrix */
			CoreCompute->BuildMatrix(prf, PFSeq.ProfileIndex, rvmatrix, WORK, NULL, 0, SeqLength);
			
			snprintf(Header,48,"Read_%u_subread_%u", HoleNumber, r); 
			Output(matrix, rvmatrix, seq, Header, prf, CoreCompute, SeqLength, OutputMethodOptions, &PrintLock);
		}
	}
	pthread_mutex_destroy(&PrintLock);
		if (WORK) _mm_free(WORK);
		if (matrix) _mm_free(matrix);
		if (rvmatrix) _mm_free(rvmatrix);
		free(SequenceIndex);
	FIN:
		
		return res;
}
