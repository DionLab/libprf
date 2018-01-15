#include <prf_config.h>
#include <stdlib.h>
#include <xmmintrin.h>
#include "pfProfile.h"
#define __USE_INLINE_FUNCTIONS__
#include "pfSequence.h"
#include "dna.h"

/* Required Symbols in the alphabet*/
const char Required_Symbols[Required_Symbols_Size] = "ACDEFGHIKLMNPQRSTVWY";

/* Stop Symbol */
char Stop_Symbols;

/* Table of Amino Acids to DNA */
const char AminoToDNA[26][4] = {
	/* A */ "GCT\0", /* B */ "xBx\0", /* C */ "TGT\0", /* D */ "GAT\0",
	/* E */ "GAA\0", /* F */ "TTT\0", /* G */ "GGT\0", /* H */ "CAT\0",
	/* I */ "ATT\0", /* J */ "xJx\0", /* K */ "AAA\0", /* L */ "CTT\0",
	/* M */ "ATG\0", /* N */ "AAT\0", /* O */ "xOx\0", /* P */ "CCT\0",
	/* Q */ "CAA\0", /* R */ "CGT\0", /* S */ "TCT\0", /* T */ "ACT\0",
	/* U */ "xUx\0", /* V */ "GTT\0", /* W */ "TGG\0", /* X */ "xXx\0",
	/* Y */ "TAT\0", /* Z */ "xZx\0"
};

/* Table of Amino Acids to IUPAC */
const char AminoToIUPAC[26][4] = {
	/* A */ "GCN\0", /* B */ "xBx\0", /* C */ "TGY\0", /* D */ "GAY\0",
	/* E */ "GAR\0", /* F */ "TTY\0", /* G */ "GGN\0", /* H */ "CAY\0",
	/* I */ "ATH\0", /* J */ "xJx\0", /* K */ "AAR\0", /* L */ "YTN\0",
	/* M */ "ATG\0", /* N */ "AAY\0", /* O */ "xOx\0", /* P */ "CCN\0",
	/* Q */ "CAR\0", /* R */ "MGN\0", /* S */ "WCN\0", /* T */ "ACN\0",
	/* U */ "xUx\0", /* V */ "GTN\0", /* W */ "TGG\0", /* X */ "xXx\0",
	/* Y */ "TAY\0", /* Z */ "xZx\0"
};

/* Translation from DNA to protein index 
 * 
 * DNA code
 * ASCII	decimal		hexadecimal	binary         
 * 	A: 	65		41		0100 0001
 * 	C:	67		43		0100 0011
 * 	G:	71		47		0100 0111
 * 	T:	84		54		0101 0100
 * 
 * 	N:	78		4E		0100 1110
 *
 */
enum DNA { A=0, C=1, G=3, T=2 };

/* 
 * Conversion table to be filled within profile reading 64 is 3 letters identification
 * and 16 is 2 first letters in the case of unique final solution
 */
unsigned char DNATripletToAminoAcidIndex[64+16] __attribute__((aligned(16)));
char DNATripletToAminoAcid[64+16] __attribute__((aligned(16)));

static inline __ALWAYS_INLINE size_t DNAToIndex(const enum DNA One, const enum DNA Two, const enum DNA Three)
{
  return (size_t) ( (Three << 4) | (Two << 2) | One );
}

static inline __ALWAYS_INLINE size_t DNA2ToIndex(const enum DNA One, const enum DNA Two)
{
  return (size_t) ( 64 | (Two << 2) | One );
}

static inline __ALWAYS_INLINE size_t DNAToIndexReversed(const enum DNA One, const enum DNA Two, const enum DNA Three)
{
  return (size_t) ( (One << 4) | (Two << 2) | Three );
}


unsigned char TranslateDNAChar3ToAminoIndex(const unsigned char * const triplet)
{
  /* Check if N is found in the 2 first letter -> unknown for sure */
  if (triplet[0] & 0x8 || triplet[1] & 0x8) {
    return (char) 0;
  } else {
    register unsigned char t0 = triplet[0] >> 1;
    register unsigned char t1 = triplet[1] >> 1;
    register unsigned char t2 = (triplet[2] & 0x8) ? 0x4 : triplet[2] >> 1;
    t0 &= 0x3;
    t1 &= 0x3;
    t2 &= 0x7;
    
    t2 <<= 4;
    t1 <<= 2;
    const size_t index = (size_t) ( (t2 | t1 | t0) & 0x7F );
    return DNATripletToAminoAcidIndex[index];
  }
}

unsigned char TranslateDNAChar3ReversedToAminoIndex(const unsigned char * const triplet)
{
  /* Check if N is found in the 2 first letter -> unknown for sure */
  if (triplet[0] & 0x8 || triplet[1] & 0x8) {
    return (char) 0;
  } else {
    register unsigned char t0 = triplet[2] >> 1;
    register unsigned char t1 = triplet[1] >> 1;
    register unsigned char t2 = (triplet[0] & 0x8) ? 0x4 : triplet[2] >> 1;
    t0 &= 0x3;
    t1 &= 0x3;
    t2 &= 0x7;
    
    t2 <<= 4;
    t1 <<= 2;
    const size_t index = (size_t) ( (t2 | t1 | t0) & 0x7F );
    return DNATripletToAminoAcidIndex[index];
  }
}

/*
 * Fill in the DNA triplet to Amino Acids table 
 */
void FillInDNATripletToAminoAcidTable(const unsigned char * restrict const Alphabet)
{
#define SET_VALUE(Id, Amino, Value) { \
  const size_t index = Id; \
  DNATripletToAminoAcidIndex[index] = Value; \
  DNATripletToAminoAcid[index] = Amino;\
}
  
  register unsigned char value;
  
  /* Set everything to zero -> unknown */
  __m128i __zero = _mm_setzero_si128();
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[ 0], __zero);
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[16], __zero);
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[32], __zero);
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[48], __zero);
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[64], __zero);
  
  memset(DNATripletToAminoAcid, '?', 80*sizeof(char));
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[ 0], __zero);
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[16], __zero);
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[32], __zero);
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[48], __zero);
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[64], __zero);
  
  
  /*  Isoleucine  	I		ATT, ATC, ATA */
  value = TranslateCharToIndex('I', Alphabet);
  SET_VALUE(DNAToIndex(A,T,T), 'I', value);
  SET_VALUE(DNAToIndex(A,T,C), 'I', value);
  SET_VALUE(DNAToIndex(A,T,A), 'I', value);
  /*  Leucine  		L		CTT, CTC, CTA, CTG, TTA, TTG */
  value = TranslateCharToIndex('L', Alphabet);
  SET_VALUE(DNAToIndex(C,T,T), 'L', value);
  SET_VALUE(DNAToIndex(C,T,C), 'L', value);
  SET_VALUE(DNAToIndex(C,T,A), 'L', value);
  SET_VALUE(DNAToIndex(C,T,G), 'L', value);
  SET_VALUE(DNAToIndex(T,T,A), 'L', value);
  SET_VALUE(DNAToIndex(T,T,G), 'L', value);
  /*  Valine		V		GTT, GTC, GTA, GTG */
  value = TranslateCharToIndex('V', Alphabet);
  SET_VALUE(DNAToIndex(G,T,T), 'V', value);
  SET_VALUE(DNAToIndex(G,T,C), 'V', value);
  SET_VALUE(DNAToIndex(G,T,A), 'V', value);
  SET_VALUE(DNAToIndex(G,T,G), 'V', value);
  /*  Phenylalanine  	F		TTT, TTC */
  value = TranslateCharToIndex('F', Alphabet);
  SET_VALUE(DNAToIndex(T,T,T), 'F', value);
  SET_VALUE(DNAToIndex(T,T,C), 'F', value);
  /*  Methionine	M		ATG */
  value = TranslateCharToIndex('M', Alphabet);
  SET_VALUE(DNAToIndex(A,T,G), 'M', value);
  /*  Cysteine 		C		TGT, TGC */
  value = TranslateCharToIndex('C', Alphabet);
  SET_VALUE(DNAToIndex(T,G,T), 'C', value);
  SET_VALUE(DNAToIndex(T,G,C), 'C', value);
  /*  Alanine      	A		GCT, GCC, GCA, GCG */
  value = TranslateCharToIndex('A', Alphabet);
  SET_VALUE(DNAToIndex(G,C,T), 'A', value);
  SET_VALUE(DNAToIndex(G,C,C), 'A', value);
  SET_VALUE(DNAToIndex(G,C,A), 'A', value);
  SET_VALUE(DNAToIndex(G,C,G), 'A', value);
  /*  Glycine  		G		GGT, GGC, GGA, GGG */
  value = TranslateCharToIndex('G', Alphabet);
  SET_VALUE(DNAToIndex(G,G,T), 'G', value);
  SET_VALUE(DNAToIndex(G,G,C), 'G', value);
  SET_VALUE(DNAToIndex(G,G,A), 'G', value);
  SET_VALUE(DNAToIndex(G,G,G), 'G', value);
  /*  Proline      	P		CCT, CCC, CCA, CCG */
  value = TranslateCharToIndex('P', Alphabet);
  SET_VALUE(DNAToIndex(C,C,T), 'P', value);
  SET_VALUE(DNAToIndex(C,C,C), 'P', value);
  SET_VALUE(DNAToIndex(C,C,A), 'P', value);
  SET_VALUE(DNAToIndex(C,C,G), 'P', value);
  /*  Threonine  	T		ACT, ACC, ACA, ACG */
  value = TranslateCharToIndex('T', Alphabet);
  SET_VALUE(DNAToIndex(A,C,T), 'T', value);
  SET_VALUE(DNAToIndex(A,C,C), 'T', value);
  SET_VALUE(DNAToIndex(A,C,A), 'T', value);
  SET_VALUE(DNAToIndex(A,C,G), 'T', value);
  /*  Serine       	S		TCT, TCC, TCA, TCG, AGT, AGC */
  value = TranslateCharToIndex('S', Alphabet);
  SET_VALUE(DNAToIndex(T,C,T), 'S', value);
  SET_VALUE(DNAToIndex(T,C,C), 'S', value);
  SET_VALUE(DNAToIndex(T,C,A), 'S', value);
  SET_VALUE(DNAToIndex(T,C,G), 'S', value);
  SET_VALUE(DNAToIndex(A,G,T), 'S', value);
  SET_VALUE(DNAToIndex(A,G,C), 'S', value);
  /*  Tyrosine  	Y		TAT, TAC */
  value = TranslateCharToIndex('Y', Alphabet);
  SET_VALUE(DNAToIndex(T,A,T), 'Y', value);
  SET_VALUE(DNAToIndex(T,A,C), 'Y', value);
  /*  Tryptophan  	W		TGG */
  value = TranslateCharToIndex('W', Alphabet);
  SET_VALUE(DNAToIndex(T,G,G), 'W', value);
  /*  Glutamine  	Q		CAA, CAG */
  value = TranslateCharToIndex('Q', Alphabet);
  SET_VALUE(DNAToIndex(C,A,A), 'Q', value);
  SET_VALUE(DNAToIndex(C,A,G), 'Q', value);
  /*  Asparagine  	N		AAT, AAC */
  value = TranslateCharToIndex('N', Alphabet);
  SET_VALUE(DNAToIndex(A,A,T), 'N', value);
  SET_VALUE(DNAToIndex(A,A,C), 'N', value);
  /*  Histidine 	H		CAT, CAC */
  value = TranslateCharToIndex('H', Alphabet);
  SET_VALUE(DNAToIndex(C,A,T), 'H', value);
  SET_VALUE(DNAToIndex(C,A,C), 'H', value);
  /*  Glutamic acid  	E		GAA, GAG */
  value = TranslateCharToIndex('E', Alphabet);
  SET_VALUE(DNAToIndex(G,A,A), 'E', value);
  SET_VALUE(DNAToIndex(G,A,G), 'E', value);
  /*  Aspartic acid 	D		GAT, GAC */
  value = TranslateCharToIndex('D', Alphabet);
  SET_VALUE(DNAToIndex(G,A,T), 'D', value);
  SET_VALUE(DNAToIndex(G,A,C), 'D', value);
  /*  Lysine       	K		AAA, AAG */
  value = TranslateCharToIndex('K', Alphabet);
  SET_VALUE(DNAToIndex(A,A,A), 'K', value);
  SET_VALUE(DNAToIndex(A,A,G), 'K', value);
  /*  Arginine  	R		CGT, CGC, CGA, CGG, AGA, AGG */
  value = TranslateCharToIndex('R', Alphabet);
  SET_VALUE(DNAToIndex(C,G,T), 'R', value);
  SET_VALUE(DNAToIndex(C,G,C), 'R', value);
  SET_VALUE(DNAToIndex(C,G,A), 'R', value);
  SET_VALUE(DNAToIndex(C,G,G), 'R', value);
  SET_VALUE(DNAToIndex(A,G,A), 'R', value);
  SET_VALUE(DNAToIndex(A,G,G), 'R', value);
  /*  Stop codons	Stop		TAA, TAG, TGA */
  SET_VALUE(DNAToIndex(T,A,A), '!', _STOP);
  SET_VALUE(DNAToIndex(T,A,G), '!', _STOP);
  SET_VALUE(DNAToIndex(T,G,A), '!', _STOP);
  
  /*
   * N at last position when two first yield a unique solution
   */
  /*  Leucine  		L		CTT, CTC, CTA, CTG, TTA, TTG */
  value = TranslateCharToIndex('L', Alphabet);
  SET_VALUE(DNA2ToIndex(C,T), 'L', value); 
  /*  Valine		V		GTT, GTC, GTA, GTG */
  value = TranslateCharToIndex('V', Alphabet);
  SET_VALUE(DNA2ToIndex(G,T), 'V', value);
  /*  Alanine      	A		GCT, GCC, GCA, GCG */
  value = TranslateCharToIndex('A', Alphabet);
  SET_VALUE(DNA2ToIndex(G,C), 'A', value);
  /*  Glycine  		G		GGT, GGC, GGA, GGG */
  value = TranslateCharToIndex('G', Alphabet);
  SET_VALUE(DNA2ToIndex(G,G), 'G', value);
  /*  Proline      	P		CCT, CCC, CCA, CCG */
  value = TranslateCharToIndex('P', Alphabet);
  SET_VALUE(DNA2ToIndex(C,C), 'P', value);
  /*  Threonine  	T		ACT, ACC, ACA, ACG */
  value = TranslateCharToIndex('T', Alphabet);
  SET_VALUE(DNA2ToIndex(A,C), 'T', value);
  /*  Serine       	S		TCT, TCC, TCA, TCG, AGT, AGC */
  value = TranslateCharToIndex('S', Alphabet);
  SET_VALUE(DNA2ToIndex(T,C), 'S', value);  
  /*  Arginine  	R		CGT, CGC, CGA, CGG, AGA, AGG */
  value = TranslateCharToIndex('R', Alphabet);
  SET_VALUE(DNA2ToIndex(C,G), 'R', value); 
  /* The rst is taken care of in the zero setting */
  #undef SET_VALUE
}

void FillInDNATripletToAminoAcidTableComplementary(const unsigned char * restrict const Alphabet)
{
#define SET_VALUE(Id, Amino, Value) { \
  const size_t index = Id; \
  DNATripletToAminoAcidIndex[index] = Value; \
  DNATripletToAminoAcid[index] = Amino;\
}
  
  register unsigned char value;
  
  /* Set everything to zero -> unknown */
  __m128i __zero = _mm_setzero_si128();
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[ 0], __zero);
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[16], __zero);
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[32], __zero);
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[48], __zero);
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[64], __zero);
  
  memset(DNATripletToAminoAcid, '?', 80*sizeof(char));
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[ 0], __zero);
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[16], __zero);
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[32], __zero);
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[48], __zero);
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[64], __zero);
  
  
  /*  Isoleucine  	I		ATT, ATC, ATA */
  value = TranslateCharToIndex('I', Alphabet);
  SET_VALUE(DNAToIndex(T,A,A), 'I', value);
  SET_VALUE(DNAToIndex(T,A,G), 'I', value);
  SET_VALUE(DNAToIndex(T,A,T), 'I', value);
  /*  Leucine  		L		CTT, CTC, CTA, CTG, TTA, TTG */
  value = TranslateCharToIndex('L', Alphabet);
  SET_VALUE(DNAToIndex(G,A,A), 'L', value);
  SET_VALUE(DNAToIndex(G,A,G), 'L', value);
  SET_VALUE(DNAToIndex(G,A,T), 'L', value);
  SET_VALUE(DNAToIndex(G,A,C), 'L', value);
  SET_VALUE(DNAToIndex(A,A,T), 'L', value);
  SET_VALUE(DNAToIndex(A,A,C), 'L', value);
  /*  Valine		V		GTT, GTC, GTA, GTG */
  value = TranslateCharToIndex('V', Alphabet);
  SET_VALUE(DNAToIndex(C,A,A), 'V', value);
  SET_VALUE(DNAToIndex(C,A,G), 'V', value);
  SET_VALUE(DNAToIndex(C,A,T), 'V', value);
  SET_VALUE(DNAToIndex(C,A,C), 'V', value);
  /*  Phenylalanine  	F		TTT, TTC */
  value = TranslateCharToIndex('F', Alphabet);
  SET_VALUE(DNAToIndex(A,A,A), 'F', value);
  SET_VALUE(DNAToIndex(A,A,G), 'F', value);
  /*  Methionine	M		ATG */
  value = TranslateCharToIndex('M', Alphabet);
  SET_VALUE(DNAToIndex(T,A,C), 'M', value);
  /*  Cysteine 		C		TGT, TGC */
  value = TranslateCharToIndex('C', Alphabet);
  SET_VALUE(DNAToIndex(A,C,A), 'C', value);
  SET_VALUE(DNAToIndex(A,C,G), 'C', value);
  /*  Alanine      	A		GCT, GCC, GCA, GCG */
  value = TranslateCharToIndex('A', Alphabet);
  SET_VALUE(DNAToIndex(C,G,A), 'A', value);
  SET_VALUE(DNAToIndex(C,G,G), 'A', value);
  SET_VALUE(DNAToIndex(C,G,T), 'A', value);
  SET_VALUE(DNAToIndex(C,G,C), 'A', value);
  /*  Glycine  		G		GGT, GGC, GGA, GGG */
  value = TranslateCharToIndex('G', Alphabet);
  SET_VALUE(DNAToIndex(C,C,A), 'G', value);
  SET_VALUE(DNAToIndex(C,C,G), 'G', value);
  SET_VALUE(DNAToIndex(C,C,T), 'G', value);
  SET_VALUE(DNAToIndex(C,C,C), 'G', value);
  /*  Proline      	P		CCT, CCC, CCA, CCG */
  value = TranslateCharToIndex('P', Alphabet);
  SET_VALUE(DNAToIndex(G,G,A), 'P', value);
  SET_VALUE(DNAToIndex(G,G,G), 'P', value);
  SET_VALUE(DNAToIndex(G,G,T), 'P', value);
  SET_VALUE(DNAToIndex(G,G,C), 'P', value);
  /*  Threonine  	T		ACT, ACC, ACA, ACG */
  value = TranslateCharToIndex('T', Alphabet);
  SET_VALUE(DNAToIndex(T,G,A), 'T', value);
  SET_VALUE(DNAToIndex(T,G,G), 'T', value);
  SET_VALUE(DNAToIndex(T,G,T), 'T', value);
  SET_VALUE(DNAToIndex(T,G,C), 'T', value);
  /*  Serine       	S		TCT, TCC, TCA, TCG, AGT, AGC */
  value = TranslateCharToIndex('S', Alphabet);
  SET_VALUE(DNAToIndex(A,G,A), 'S', value);
  SET_VALUE(DNAToIndex(A,G,G), 'S', value);
  SET_VALUE(DNAToIndex(A,G,T), 'S', value);
  SET_VALUE(DNAToIndex(A,G,C), 'S', value);
  SET_VALUE(DNAToIndex(T,C,A), 'S', value);
  SET_VALUE(DNAToIndex(T,C,G), 'S', value);
  /*  Tyrosine  	Y		TAT, TAC */
  value = TranslateCharToIndex('Y', Alphabet);
  SET_VALUE(DNAToIndex(A,T,A), 'Y', value);
  SET_VALUE(DNAToIndex(A,T,G), 'Y', value);
  /*  Tryptophan  	W		TGG */
  value = TranslateCharToIndex('W', Alphabet);
  SET_VALUE(DNAToIndex(A,C,C), 'W', value);
  /*  Glutamine  	Q		CAA, CAG */
  value = TranslateCharToIndex('Q', Alphabet);
  SET_VALUE(DNAToIndex(G,T,T), 'Q', value);
  SET_VALUE(DNAToIndex(G,T,C), 'Q', value);
  /*  Asparagine  	N		AAT, AAC */
  value = TranslateCharToIndex('N', Alphabet);
  SET_VALUE(DNAToIndex(T,T,A), 'N', value);
  SET_VALUE(DNAToIndex(T,T,G), 'N', value);
  /*  Histidine 	H		CAT, CAC */
  value = TranslateCharToIndex('H', Alphabet);
  SET_VALUE(DNAToIndex(G,T,A), 'H', value);
  SET_VALUE(DNAToIndex(G,T,G), 'H', value);
  /*  Glutamic acid  	E		GAA, GAG */
  value = TranslateCharToIndex('E', Alphabet);
  SET_VALUE(DNAToIndex(C,T,T), 'E', value);
  SET_VALUE(DNAToIndex(C,T,C), 'E', value);
  /*  Aspartic acid 	D		GAT, GAC */
  value = TranslateCharToIndex('D', Alphabet);
  SET_VALUE(DNAToIndex(C,T,A), 'D', value);
  SET_VALUE(DNAToIndex(C,T,G), 'D', value);
  /*  Lysine       	K		AAA, AAG */
  value = TranslateCharToIndex('K', Alphabet);
  SET_VALUE(DNAToIndex(T,T,T), 'K', value);
  SET_VALUE(DNAToIndex(T,T,C), 'K', value);
  /*  Arginine  	R		CGT, CGC, CGA, CGG, AGA, AGG */
  value = TranslateCharToIndex('R', Alphabet);
  SET_VALUE(DNAToIndex(G,C,A), 'R', value);
  SET_VALUE(DNAToIndex(G,C,G), 'R', value);
  SET_VALUE(DNAToIndex(G,C,T), 'R', value);
  SET_VALUE(DNAToIndex(G,C,C), 'R', value);
  SET_VALUE(DNAToIndex(T,C,T), 'R', value);
  SET_VALUE(DNAToIndex(T,C,C), 'R', value);
  /*  Stop codons	Stop		TAA, TAG, TGA */
  SET_VALUE(DNAToIndex(A,T,T), '!', _STOP);
  SET_VALUE(DNAToIndex(A,T,C), '!', _STOP);
  SET_VALUE(DNAToIndex(A,C,T), '!', _STOP);
  
  /*
   * N at last position when two first yield a unique solution
   */
  /*  Leucine  		L		CTT, CTC, CTA, CTG, TTA, TTG */
  value = TranslateCharToIndex('L', Alphabet);
  SET_VALUE(DNA2ToIndex(G,A), 'L', value); 
  /*  Valine		V		GTT, GTC, GTA, GTG */
  value = TranslateCharToIndex('V', Alphabet);
  SET_VALUE(DNA2ToIndex(C,A), 'V', value);
  /*  Alanine      	A		GCT, GCC, GCA, GCG */
  value = TranslateCharToIndex('A', Alphabet);
  SET_VALUE(DNA2ToIndex(C,G), 'A', value);
  /*  Glycine  		G		GGT, GGC, GGA, GGG */
  value = TranslateCharToIndex('G', Alphabet);
  SET_VALUE(DNA2ToIndex(C,C), 'G', value);
  /*  Proline      	P		CCT, CCC, CCA, CCG */
  value = TranslateCharToIndex('P', Alphabet);
  SET_VALUE(DNA2ToIndex(G,G), 'P', value);
  /*  Threonine  	T		ACT, ACC, ACA, ACG */
  value = TranslateCharToIndex('T', Alphabet);
  SET_VALUE(DNA2ToIndex(T,G), 'T', value);
  /*  Serine       	S		TCT, TCC, TCA, TCG, AGT, AGC */
  value = TranslateCharToIndex('S', Alphabet);
  SET_VALUE(DNA2ToIndex(A,G), 'S', value);  
  /*  Arginine  	R		CGT, CGC, CGA, CGG, AGA, AGG */
  value = TranslateCharToIndex('R', Alphabet);
  SET_VALUE(DNA2ToIndex(G,C), 'R', value); 
  /* The rst is taken care of in the zero setting */
  #undef SET_VALUE
}

void FillInDNATripletToAminoAcidTableReversed(const unsigned char * restrict const Alphabet)
{
  #define SET_VALUE(Id, Amino, Value) { \
  const size_t index = Id; \
  DNATripletToAminoAcidIndex[index] = Value; \
  DNATripletToAminoAcid[index] = Amino;\
}
  
  register unsigned char value;
  
  /* Set everything to zero -> unknown */
  __m128i __zero = _mm_setzero_si128();
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[ 0], __zero);
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[16], __zero);
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[32], __zero);
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[48], __zero);
  _mm_store_si128((__m128i*) &DNATripletToAminoAcidIndex[64], __zero);
  
  memset(DNATripletToAminoAcid, '?', 80*sizeof(char));
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[ 0], __zero);
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[16], __zero);
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[32], __zero);
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[48], __zero);
//   _mm_store_si128((__m128i*) &DNATripletToAminoAcid[64], __zero);
  
  
  /*  Isoleucine  	I		ATT, ATC, ATA */
  value = TranslateCharToIndex('I', Alphabet);
  SET_VALUE(DNAToIndexReversed(A,T,T), 'I', value);
  SET_VALUE(DNAToIndexReversed(A,T,C), 'I', value);
  SET_VALUE(DNAToIndexReversed(A,T,A), 'I', value);
  /*  Leucine  		L		CTT, CTC, CTA, CTG, TTA, TTG */
  value = TranslateCharToIndex('L', Alphabet);
  SET_VALUE(DNAToIndexReversed(C,T,T), 'L', value);
  SET_VALUE(DNAToIndexReversed(C,T,C), 'L', value);
  SET_VALUE(DNAToIndexReversed(C,T,A), 'L', value);
  SET_VALUE(DNAToIndexReversed(C,T,G), 'L', value);
  SET_VALUE(DNAToIndexReversed(T,T,A), 'L', value);
  SET_VALUE(DNAToIndexReversed(T,T,G), 'L', value);
  /*  Valine		V		GTT, GTC, GTA, GTG */
  value = TranslateCharToIndex('V', Alphabet);
  SET_VALUE(DNAToIndexReversed(G,T,T), 'V', value);
  SET_VALUE(DNAToIndexReversed(G,T,C), 'V', value);
  SET_VALUE(DNAToIndexReversed(G,T,A), 'V', value);
  SET_VALUE(DNAToIndexReversed(G,T,G), 'V', value);
  /*  Phenylalanine  	F		TTT, TTC */
  value = TranslateCharToIndex('F', Alphabet);
  SET_VALUE(DNAToIndexReversed(T,T,T), 'F', value);
  SET_VALUE(DNAToIndexReversed(T,T,C), 'F', value);
  /*  Methionine	M		ATG */
  value = TranslateCharToIndex('M', Alphabet);
  SET_VALUE(DNAToIndexReversed(A,T,G), 'M', value);
  /*  Cysteine 		C		TGT, TGC */
  value = TranslateCharToIndex('C', Alphabet);
  SET_VALUE(DNAToIndexReversed(T,G,T), 'C', value);
  SET_VALUE(DNAToIndexReversed(T,G,C), 'C', value);
  /*  Alanine      	A		GCT, GCC, GCA, GCG */
  value = TranslateCharToIndex('A', Alphabet);
  SET_VALUE(DNAToIndexReversed(G,C,T), 'A', value);
  SET_VALUE(DNAToIndexReversed(G,C,C), 'A', value);
  SET_VALUE(DNAToIndexReversed(G,C,A), 'A', value);
  SET_VALUE(DNAToIndexReversed(G,C,G), 'A', value);
  /*  Glycine  		G		GGT, GGC, GGA, GGG */
  value = TranslateCharToIndex('G', Alphabet);
  SET_VALUE(DNAToIndexReversed(G,G,T), 'G', value);
  SET_VALUE(DNAToIndexReversed(G,G,C), 'G', value);
  SET_VALUE(DNAToIndexReversed(G,G,A), 'G', value);
  SET_VALUE(DNAToIndexReversed(G,G,G), 'G', value);
  /*  Proline      	P		CCT, CCC, CCA, CCG */
  value = TranslateCharToIndex('P', Alphabet);
  SET_VALUE(DNAToIndexReversed(C,C,T), 'P', value);
  SET_VALUE(DNAToIndexReversed(C,C,C), 'P', value);
  SET_VALUE(DNAToIndexReversed(C,C,A), 'P', value);
  SET_VALUE(DNAToIndexReversed(C,C,G), 'P', value);
  /*  Threonine  	T		ACT, ACC, ACA, ACG */
  value = TranslateCharToIndex('T', Alphabet);
  SET_VALUE(DNAToIndexReversed(A,C,T), 'T', value);
  SET_VALUE(DNAToIndexReversed(A,C,C), 'T', value);
  SET_VALUE(DNAToIndexReversed(A,C,A), 'T', value);
  SET_VALUE(DNAToIndexReversed(A,C,G), 'T', value);
  /*  Serine       	S		TCT, TCC, TCA, TCG, AGT, AGC */
  value = TranslateCharToIndex('S', Alphabet);
  SET_VALUE(DNAToIndexReversed(T,C,T), 'S', value);
  SET_VALUE(DNAToIndexReversed(T,C,C), 'S', value);
  SET_VALUE(DNAToIndexReversed(T,C,A), 'S', value);
  SET_VALUE(DNAToIndexReversed(T,C,G), 'S', value);
  SET_VALUE(DNAToIndexReversed(A,G,T), 'S', value);
  SET_VALUE(DNAToIndexReversed(A,G,C), 'S', value);
  /*  Tyrosine  	Y		TAT, TAC */
  value = TranslateCharToIndex('Y', Alphabet);
  SET_VALUE(DNAToIndexReversed(T,A,T), 'Y', value);
  SET_VALUE(DNAToIndexReversed(T,A,C), 'Y', value);
  /*  Tryptophan  	W		TGG */
  value = TranslateCharToIndex('W', Alphabet);
  SET_VALUE(DNAToIndexReversed(T,G,G), 'W', value);
  /*  Glutamine  	Q		CAA, CAG */
  value = TranslateCharToIndex('Q', Alphabet);
  SET_VALUE(DNAToIndexReversed(C,A,A), 'Q', value);
  SET_VALUE(DNAToIndexReversed(C,A,G), 'Q', value);
  /*  Asparagine  	N		AAT, AAC */
  value = TranslateCharToIndex('N', Alphabet);
  SET_VALUE(DNAToIndexReversed(A,A,T), 'N', value);
  SET_VALUE(DNAToIndexReversed(A,A,C), 'N', value);
  /*  Histidine 	H		CAT, CAC */
  value = TranslateCharToIndex('H', Alphabet);
  SET_VALUE(DNAToIndexReversed(C,A,T), 'H', value);
  SET_VALUE(DNAToIndexReversed(C,A,C), 'H', value);
  /*  Glutamic acid  	E		GAA, GAG */
  value = TranslateCharToIndex('E', Alphabet);
  SET_VALUE(DNAToIndexReversed(G,A,A), 'E', value);
  SET_VALUE(DNAToIndexReversed(G,A,G), 'E', value);
  /*  Aspartic acid 	D		GAT, GAC */
  value = TranslateCharToIndex('D', Alphabet);
  SET_VALUE(DNAToIndexReversed(G,A,T), 'D', value);
  SET_VALUE(DNAToIndexReversed(G,A,C), 'D', value);
  /*  Lysine       	K		AAA, AAG */
  value = TranslateCharToIndex('K', Alphabet);
  SET_VALUE(DNAToIndexReversed(A,A,A), 'K', value);
  SET_VALUE(DNAToIndexReversed(A,A,G), 'K', value);
  /*  Arginine  	R		CGT, CGC, CGA, CGG, AGA, AGG */
  value = TranslateCharToIndex('R', Alphabet);
  SET_VALUE(DNAToIndexReversed(C,G,T), 'R', value);
  SET_VALUE(DNAToIndexReversed(C,G,C), 'R', value);
  SET_VALUE(DNAToIndexReversed(C,G,A), 'R', value);
  SET_VALUE(DNAToIndexReversed(C,G,G), 'R', value);
  SET_VALUE(DNAToIndexReversed(A,G,A), 'R', value);
  SET_VALUE(DNAToIndexReversed(A,G,G), 'R', value);
  /*  Stop codons	Stop		TAA, TAG, TGA */
  SET_VALUE(DNAToIndexReversed(T,A,A), '!', _STOP);
  SET_VALUE(DNAToIndexReversed(T,A,G), '!', _STOP);
  SET_VALUE(DNAToIndexReversed(T,G,A), '!', _STOP);
  
  /*
   * N at last position when two first yield a unique solution
   */
  /*  Leucine  		L		CTT, CTC, CTA, CTG, TTA, TTG */
  value = TranslateCharToIndex('L', Alphabet);
  SET_VALUE(DNA2ToIndex(C,T), 'L', value); 
  /*  Valine		V		GTT, GTC, GTA, GTG */
  value = TranslateCharToIndex('V', Alphabet);
  SET_VALUE(DNA2ToIndex(G,T), 'V', value);
  /*  Alanine      	A		GCT, GCC, GCA, GCG */
  value = TranslateCharToIndex('A', Alphabet);
  SET_VALUE(DNA2ToIndex(G,C), 'A', value);
  /*  Glycine  		G		GGT, GGC, GGA, GGG */
  value = TranslateCharToIndex('G', Alphabet);
  SET_VALUE(DNA2ToIndex(G,G), 'G', value);
  /*  Proline      	P		CCT, CCC, CCA, CCG */
  value = TranslateCharToIndex('P', Alphabet);
  SET_VALUE(DNA2ToIndex(C,C), 'P', value);
  /*  Threonine  	T		ACT, ACC, ACA, ACG */
  value = TranslateCharToIndex('T', Alphabet);
  SET_VALUE(DNA2ToIndex(A,C), 'T', value);
  /*  Serine       	S		TCT, TCC, TCA, TCG, AGT, AGC */
  value = TranslateCharToIndex('S', Alphabet);
  SET_VALUE(DNA2ToIndex(T,C), 'S', value);  
  /*  Arginine  	R		CGT, CGC, CGA, CGG, AGA, AGG */
  value = TranslateCharToIndex('R', Alphabet);
  SET_VALUE(DNA2ToIndex(C,G), 'R', value); 
  /* The rst is taken care of in the zero setting */
}

