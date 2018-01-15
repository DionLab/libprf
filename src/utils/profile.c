/*******************************************************
                        PFTOOLS
 *******************************************************
  Mai 6, 2013 output.c
 *******************************************************
 (C) 2013 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#include "prf_config.h"
#include <stdlib.h>
#include <string.h>
#define __USE_INLINE_FUNCTIONS__
#include "pfProfile.h"

/*
 * Treats the level of cutoff by seeking the index corresponding to the wanted level number,
 * then storing it into the profile Level member.
 * BE CAREFUL THIS IS THE INDEX AND NOT THE VALUE OF THE LEVEL !!!
 * 
 * The function shall the index to the level if found, -1 if not found and -2 upon multiple
 * existance of the same level.
 */
int SeekProfileLevel(const struct Profile * const prf, const int Level)
{
  int index = -1;
  const SCutOffItem * const restrict cutItems = prf->CutOffData.Values;
  _Bool levelfound = false;
  const size_t N = prf->CutOffData.JCUT;
  for (size_t icut=0; icut<N; ++icut) {
    if ( cutItems[icut].MCLE == Level ) {
      index = (int) icut;
      if (levelfound) return -2;
      levelfound = true;
    }
  }
  return index;
}
/*
 * This function returns the index of the highest priority normalization 
 * matching the given mode or -1 if not found.
 */
int SeekProfileMode(const struct Profile * const prf, const int Mode)
{
  int index = -1;
  register const int SeekMode = Mode;
  const SNormalizationItem * const restrict NormItems = &(prf->NormalizationData.Values[0]);
  
  int Priority =  -1;
  
  for (int iNormalizationMode=0; iNormalizationMode<prf->NormalizationData.JNOR; ++iNormalizationMode) {
    if ( SeekMode == NormItems[iNormalizationMode].NNOR ) {
      if (index < 0) {
				Priority = NormItems[iNormalizationMode].NNPR;
				index    = iNormalizationMode;
      } 
      else if (NormItems[iNormalizationMode].NNPR < Priority) {
				Priority = NormItems[iNormalizationMode].NNPR;
				index    = iNormalizationMode;
      }
    }
  }
  return index;
}

int SetProfileMode(struct Profile * const prf, const int Mode)
{
    const int index = SeekProfileMode(prf, Mode);
    if (index >= 0) {
      prf->ModeIndex = (short int) index;
      prf->NormalizationCoefs              = prf->NormalizationData.Values[index].RNOP;
      register const int NormalizationType = prf->NormalizationData.Values[index].MNOR;
      prf->NormalizationType               = NormalizationType;
      switch (NormalizationType) {
	case LINEAR:
	  prf->NormalizedToRaw = &N2R_1;
	  prf->RawToNormalized = &R2N_1;
	  break;
	case GLE_ZSCORE:
	  prf->NormalizedToRaw = &N2R_2;
	  prf->RawToNormalized = &R2N_2;
	  break;
	case GLE_ZSCAVE:
	  prf->NormalizedToRaw = &N2R_3;
	  prf->RawToNormalized = &R2N_3;
	  if (InitAverage(prf)<0) return -5;
	  break;
	default:
	  return -4;
      }
    }
    else {
      prf->ModeIndex          = (short int) -1;
      prf->NormalizationCoefs = NULL;
      prf->NormalizationType  = 0;;
      prf->NormalizedToRaw    = NULL;
      prf->RawToNormalized    = NULL;
    }
    const int index2 = SeekProfileMode(prf, -Mode);
    prf->HeuristicModeIndex = (index2 >= 0) ? (short int) index : -1;
    
    return index;
}

void SeekProfileLevelAndMode(int * const restrict LevelIndex, int * const restrict ModeIndexWithinLevel, int * const restrict ModeIndex,
			                       int * const restrict HeuristicModeIndex, const struct Profile * const prf, const int Level, const int Mode)
{
  int localLevelIndex = -1;
  int localModeIndexWithinLevel = -1;
  int localHeuristicModeIndex = -1;
  
  /* First test if Mode exists */
  int localModeIndex = SeekProfileMode(prf, Mode);
  if (localModeIndex < 0 ) goto FIN;
  localHeuristicModeIndex = SeekProfileMode(prf, -Mode);
  
  /* Now scan the level to see if there exists one satisfying Level and Mode */
  if ( Level < MAXC ) {
    {
      const SCutOffItem * const restrict cutItem = prf->CutOffData.Values;
      register const int M = prf->CutOffData.JCUT;
      _Bool Levelfound = false;
      for (int iCutoff=0; iCutoff<M; ++iCutoff) {
	if (cutItem[iCutoff].MCLE == Level ) {
	  localLevelIndex = iCutoff;
	  if (Levelfound) {
	    localLevelIndex = -2;
	    goto FIN;
	  }
	  Levelfound = true;
	}
      }
      if (!Levelfound) {
	localLevelIndex = -1;
	goto FIN;
      }
    }
    
    const SCutOffItem * const restrict cutItem = &prf->CutOffData.Values[localLevelIndex];
    register const int N = cutItem->JCNM;
    for (int iCutoffMode=0; iCutoffMode<N; ++iCutoffMode) {
      const int CutoffMode = cutItem->MCUT[iCutoffMode];
      if (CutoffMode == Mode) {
				if (localModeIndexWithinLevel < 0) 
					localModeIndexWithinLevel = iCutoffMode;
				else {
					localModeIndexWithinLevel = -2;
					goto FIN;
				}
      }
    }
  }
  FIN:
  *LevelIndex           = localLevelIndex;
  *ModeIndexWithinLevel = localModeIndexWithinLevel;
  *ModeIndex            = localModeIndex;
  *HeuristicModeIndex   = localHeuristicModeIndex;
//   fprintf(stderr,"Level: %i Mode within level: %i Mode: %i\n", localLevelIndex, localModeIndexWithinLevel, localModeIndex);
  return;
}
/*
 * SetProfileLevelAndMod returns:
 * --------------------------------
 *  0 : all OK
 * -1 : Level not found
 * -2 : Mutliple levels satisfy conditions
 * -3 : Mode not found
 * -4 : Unknown Mode in Level (might never exist, need to check)
 * -5 : Average initialization failed
 */
int SetProfileLevelAndMode(struct Profile * const prf, const int Level, const int Mode)
{
  int lLevel=-1, lModeWithinLevel=-1, lMode=-1, lHeuristicMode=-1;
  SeekProfileLevelAndMode(&lLevel, &lModeWithinLevel, &lMode, &lHeuristicMode, prf, Level, Mode);
  
  prf->HeuristicCutOff      = 0;
  prf->ModeIndex            = -1;
  prf->ModeIndexWithinLevel = -1;
  prf->HeuristicModeIndex   = -1;
  prf->LevelIndex           = -1;
  prf->NormalizedCutOff     = 0.0f;
  prf->NormalizationCoefs   = NULL;
  prf->NormalizationType    = 0;
  prf->NormalizedToRaw      = NULL;
  prf->RawToNormalized      = NULL;
  
  if(lMode < 0) return -3;
  
  if (lLevel >= 0) {
    if (lMode >= 0) {
      prf->LevelIndex                      = (short int) lLevel;
      prf->CutOff                          = prf->CutOffData.Values[lLevel].ICUT;
      prf->HeuristicCutOff                 = prf->CutOffData.Values[lLevel].HCUT;
      prf->ModeIndex                       = (short int) lMode;
      prf->ModeIndexWithinLevel            = (short int) lModeWithinLevel;
      prf->NormalizedCutOff                = prf->CutOffData.Values[lLevel].RCUT[lMode];
      prf->NormalizationCoefs              = prf->NormalizationData.Values[lMode].RNOP;
      register const int NormalizationType = prf->NormalizationData.Values[lMode].MNOR;
      prf->NormalizationType               = NormalizationType;
      switch (NormalizationType) {
				case LINEAR:
					prf->NormalizedToRaw = &N2R_1;
					prf->RawToNormalized = &R2N_1;
					break;
				case GLE_ZSCORE:
					prf->NormalizedToRaw = &N2R_2;
					prf->RawToNormalized = &R2N_2;
					break;
				case GLE_ZSCAVE:
					prf->NormalizedToRaw = &N2R_3;
					prf->RawToNormalized = &R2N_3;
					if (InitAverage(prf)<0) return -5;
					break;
				default:
					return -4;
      }
      if (lHeuristicMode>=0) prf->HeuristicModeIndex = lHeuristicMode;
    }
    else {
      return lMode - 2; 
    }
  }
  else {
    return lLevel;
  }
  return 0;
}

/*
 * Complementary Mode 
 * modifications for DNA, IUPAC, etc ...
 */
int ComplementAlphabet(struct Profile * const prf, const enum ComplementaryMode CplmMode)
{
#define SWAP_ALPHABET(x,y) { \
  register const unsigned char ctmp = Alphabet_Mapping[x - (unsigned char) 'A'];\
  Alphabet_Mapping[x - (unsigned char) 'A'] = Alphabet_Mapping[y - (unsigned char) 'A'];\
  Alphabet_Mapping[y - (unsigned char) 'A'] = ctmp;\
}
  if (prf->Type == PF_MATRIX) { 
    unsigned char * const restrict Alphabet_Mapping = prf->Alphabet_Mapping;
//     const char * const restrict CABC = prf->CABC;
    const size_t Alphabet_Length = prf->Alphabet_Length;
  
    // Update the mapping
    switch (CplmMode) {
      case DNA:
	SWAP_ALPHABET('A', 'T');
	SWAP_ALPHABET('C', 'G');
	return 0;
	break;
      case IUPAC:
	SWAP_ALPHABET('A', 'T');
	SWAP_ALPHABET('C', 'G');
// 	SWAP_ALPHABET('G', 'C');
// 	SWAP_ALPHABET('T', 'A');
	SWAP_ALPHABET('M', 'K');
	SWAP_ALPHABET('R', 'Y');
// 	SWAP_ALPHABET('W', 'W');
// 	SWAP_ALPHABET('S', 'S');
// 	SWAP_ALPHABET('Y', 'R');
// 	SWAP_ALPHABET('K', 'M');
	SWAP_ALPHABET('V', 'B');
	SWAP_ALPHABET('H', 'D');
// 	SWAP_ALPHABET('D', 'H');
// 	SWAP_ALPHABET('B', 'V');
// 	SWAP_ALPHABET('N', 'N');
	return 0;
	break;
      default:
	fputs("Unknown complementary mode\n", stderr);
	return 1;
    }
  }
  else {
      fputs("Complementary option may not be used on pattern profile\n", stderr);
      return 1;
  }
#undef SWAP_ALPHABET
}

/* 
 * Sequence to struct Profile
 */
int DNASeqToProfile(const char * const restrict Identification, char * const restrict Sequence,
                    const size_t length, struct Profile * const restrict prf)
{
	char DefaultMatchSymbol;
  char DefaultInsertionSymbol;
	
	/* Clean */
	memset(prf, 0, sizeof(struct Profile));
	
	/* Allocate Scores matrices */
	if (AllocateScores(&(prf->Scores), 5, length) != 0) {
		goto bail;
	}
	
	/* Set default values */
	prf->isCircular         = false;
	prf->CompleteCycleOnly  = false;
	prf->isReversed         = false;
	prf->ReverseSequence    = false;
	prf->ComplementAlphabet = false;

  /*   - disjoint mode */
  prf->DisjointData.MDIS = 1;
  strcpy(prf->DisjointData.CDIS[0], "UNIQUE\0");
  strcpy(prf->DisjointData.CDIS[1], "PROTECT\0");
  prf->DisjointData.JDIP[0] = 0;
  prf->DisjointData.JDIP[1] = 2;
  

  /*   - normalization modes */
  prf->NormalizationData.JNOR = 0;
  strcpy(prf->NormalizationData.CNOR[0], "LINEAR\0");
  strcpy(prf->NormalizationData.CNOR[1], "GLE_ZSCORE\0");
  strcpy(prf->NormalizationData.CNOR[2], "GLE_ZSCAVE\0");
  prf->NormalizationData.JNOP[0] = 2;
  prf->NormalizationData.JNOP[1] = 5;
  prf->NormalizationData.JNOP[2] = 5;

  for (int i=0; i<MAXN; ++i) {
      prf->NormalizationData.Values[i].NNOR = i;
      prf->NormalizationData.Values[i].NNPR = i;
  }

  /*   - cut-off */
  prf->CutOffData.JCUT = 0;
  
  /* Pattern */
  prf->Pattern = NULL;
	
	InitializeDefault(&(prf->Scores), &DefaultMatchSymbol, &DefaultInsertionSymbol, prf->Scores.Insertion.AlignStep*(length+1UL));
	
	/* Characteristic of sequence */
	strncpy(prf->Identification, Identification, 64); 
	prf->Length = length;
	prf->Sequence = Sequence;
	prf->Alphabet_Length = 5;
	prf->CABC[0] = 'X';
	prf->CABC[1] = 'A';
	prf->CABC[2] = 'C';
	prf->CABC[3] = 'G';
	prf->CABC[4] = 'T';
	prf->CABC[4] = 'N';
	
	prf->Alphabet_Mapping['A' - 'A'] = 1;
	prf->Alphabet_Mapping['C' - 'A'] = 2;
	prf->Alphabet_Mapping['G' - 'A'] = 3;
	prf->Alphabet_Mapping['T' - 'A'] = 4;
	prf->Alphabet_Mapping['N' - 'A'] = 5;
	
	/* Scoring */
#define SCORE_I		   0
#define SCORE_IM		 0 
#define SCORE_II		 0 
#define SCORE_ID		 NLOW
#define SCORE_IX		 0

#define SCORE_D	  	-1
#define SCORE_DM		 0
#define SCORE_DI		 NLOW
#define SCORE_DD		 0
#define SCORE_DX		 0

#define SCORE_m		  -1
#define SCORE_M		   3
#define SCORE_MM		 1
#define SCORE_MI	 -15
#define SCORE_MD		-3
#define SCORE_MX		 0 

#define SCORE_BM		 0
#define SCORE_BD		-2

	const char * restrict SeqPtr = Sequence;
	const static StoredIntegerFormat DefaultMatch[5 + ALPHABET_EXTRA_LETTERS] = {
		0, SCORE_m, SCORE_m, SCORE_m, SCORE_m, 0, SCORE_D, NLOW
	};
	StoredIntegerFormat * restrict Matches = prf->Scores.Match.Alphabet;
	for (size_t iprf=0; iprf<= length; iprf++) {
		for (size_t i=0UL; i<(5 + ALPHABET_EXTRA_LETTERS); i++) Matches[i] = DefaultMatch[i];
		switch(*Sequence) {
			case 'A':
				Matches[1] = SCORE_M;
				break;
			case 'C':
				Matches[2] = SCORE_M;
				break;
			case 'G':
				Matches[3] = SCORE_M;
				break;
			case 'T':
				Matches[4] = SCORE_M;
				break;
			case 'N':
				Matches[1] = 0;
				Matches[2] = 0;
				Matches[3] = 0;
				Matches[4] = 0;
				break;
			default:
				;
		}
	}
		
	return 0;
	bail:
		FreeScores(&(prf->Scores));
		return 1;
}
