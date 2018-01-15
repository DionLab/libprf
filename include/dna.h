
/*
 * 20 Amino acids, their single-letter data-base codes (SLC), and their corresponding DNA codons
 * taken from http://www.cbs.dtu.dk/courses/27619/codon.html
 * ==============================================================================  
 *  Amino Acid		SLC		DNA codons
 * ------------------------------------------------------------------------------
 *  Isoleucine  	I			ATT, ATC, ATA
 *  Leucine  			L			CTT, CTC, CTA, CTG, TTA, TTG
 *  Valine				V			GTT, GTC, GTA, GTG
 *  Phenylalanine	F			TTT, TTC
 *  Methionine		M			ATG
 *  Cysteine 			C			TGT, TGC
 *  Alanine      	A			GCT, GCC, GCA, GCG
 *  Glycine  			G			GGT, GGC, GGA, GGG
 *  Proline      	P			CCT, CCC, CCA, CCG
 *  Threonine  		T			ACT, ACC, ACA, ACG
 *  Serine       	S			TCT, TCC, TCA, TCG, AGT, AGC
 *  Tyrosine  		Y			TAT, TAC
 *  Tryptophan  	W			TGG
 *  Glutamine  		Q			CAA, CAG
 *  Asparagine  	N			AAT, AAC
 *  Histidine 		H			CAT, CAC
 *  Glutamic acid	E			GAA, GAG
 *  Aspartic acid	D			GAT, GAC
 *  Lysine       	K			AAA, AAG
 *  Arginine  		R			CGT, CGC, CGA, CGG, AGA, AGG
 *  Stop codons		Stop	TAA, TAG, TGA
 * -------------------------------------------------------------------------------
 *
 * In this table, the twenty amino acids found in proteins are listed, along with
 * the single-letter code used to represent these amino acids in protein data bases.
 * The DNA codons representing each amino acid are also listed. All 64 possible
 * 3-letter combinations of the DNA coding units T, C, A and G are used either to
 * encode one of these amino acids or as one of the three stop codons that signals
 * the end of a sequence. While DNA can be decoded unambiguously, it is not possible
 * to predict a DNA sequence from its protein sequence. Because most amino acids
 * have multiple codons, a number of possible DNA sequences might represent the
 * same protein sequence.
 */

/* Required Symbols in the alphabet*/
#define Required_Symbols_Size 20
extern const char Required_Symbols[Required_Symbols_Size];

/* Stop Symbol must be define in profile */
// extern char Stop_Symbols;

/* Table of Amino Acids to DNA */
extern const char AminoToDNA[26][4];

/* Table of Amino Acids to IUPAC */
extern const char AminoToIUPAC[26][4];

/* 
 * Conversion table to be filled within profile reading 64 is 3 letters identification
 * and 16 is 2 first letters in the case of unique final solution
 */
extern unsigned char DNATripletToAminoAcidIndex[64+16] __attribute__((aligned(16)));
extern char DNATripletToAminoAcid[64+16] __attribute__((aligned(16)));

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
 */

static inline void TranslateDNACharToAminoIndex(unsigned char * const restrict Indices,
																								char * const restrict AminoAcids, const unsigned char c)
{
	const register size_t t2 = ((c >> 1) & 0x3) << 4;
	const unsigned char * const restrict AminoIndex = &DNATripletToAminoAcidIndex[t2]; 
	const char * const restrict Amino = &DNATripletToAminoAcid[t2];
	for (size_t i=0; i<16; ++i) {
		Indices[i]    = AminoIndex[i];
		AminoAcids[i] = Amino[i];
	}
}

static inline void TranslateDNAChar2ToAminoIndex(unsigned char * const restrict Indices,
																								 char * const restrict AminoAcids,
																								 const unsigned char * const duplet)
{
	const register size_t t1 = ((duplet[0] >> 1) & 0x3) << 2;
	const register size_t t2 = ((duplet[1] >> 1) & 0x3) << 4;
	const unsigned char * const restrict AminoIndex = &DNATripletToAminoAcidIndex[t2|t1]; 
	const char * const restrict Amino = &DNATripletToAminoAcid[t2|t1];
	Indices[0] = AminoIndex[0];
	Indices[1] = AminoIndex[1];
	Indices[2] = AminoIndex[2];
	Indices[3] = AminoIndex[3];
	AminoAcids[0] = Amino[0];
	AminoAcids[1] = Amino[1];
	AminoAcids[2] = Amino[2];
	AminoAcids[3] = Amino[3];
}

/* Funtion used to get index given the last DNA base */ 
static inline unsigned char TranslateDNACharToIndex ( const unsigned char * const c)
{
	return (*c >> 1) & 0x3;
}
/* Funtion used to get index given the last 2 DNA bases */
static inline unsigned char TranslateDNAChar2ToIndex ( const unsigned char * const c)
{
	unsigned char u = (c[0] >> 1) & 0x3;
	unsigned char t = (c[1] >> 1) & 0x3;
	return (t << 2) | u;
}

static inline char TranslateDNACharToAmino(const unsigned char * const triplet)
{
	/* Check if N is found in the 2 first letter -> unknown for sure */
	if (triplet[0] & 0x8 || triplet[1] & 0x8) {
		return '?';
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
		return DNATripletToAminoAcid[index];
	}
}

#ifndef PFSEQ
#define PFSEQ
typedef struct PFSequence {
	unsigned char * ProfileIndex;
	size_t Length;
} PFSequence;
#endif

static inline PFSequence * TranslateDNASequenceToInterleaveIndex(PFSequence * const Sequence, const unsigned char * restrict const Alphabet)
{
	unsigned char * restrict CharPtr = Sequence->ProfileIndex;
	size_t i, counter = 0, MaxLength = Sequence->Length;
	
	/* Remove unknown letter and move to uppercase */
	for (i=0; i<MaxLength; ++i) {
		register unsigned char value  = ( ( CharPtr[i] >= (unsigned char) 'a' ) ? CharPtr[i] - ((unsigned char) 'a' - (unsigned char) 'A') : CharPtr[i] );
		switch (value)
		{
			case (unsigned char) 'A':
			case (unsigned char) 'C':
			case (unsigned char) 'G':
			case (unsigned char) 'T':
			case (unsigned char) 'N':
				CharPtr[counter++] = value;
				break;
			default:
				;
		}
	}
	
	for (i=2; i<counter; ++i) *CharPtr++ = TranslateDNACharToIndex(CharPtr);
	/* fill up missing phase data with unknown */
	for (i=0; i<2-(counter % 3); ++i) CharPtr[counter+i] = 0;
	
	Sequence->Length = counter/3;
	return Sequence;
}

static inline PFSequence * TranslateDNASequenceToInterleave(PFSequence * const Sequence, const unsigned char * restrict const Alphabet)
{
	unsigned char * restrict CharPtr = Sequence->ProfileIndex;
	size_t i, counter = 0, MaxLength = Sequence->Length;
	
	/* Remove unknown letter and move to uppercase */
	for (i=0; i<MaxLength; ++i) {
		register unsigned char value  = ( ( CharPtr[i] >= (unsigned char) 'a' ) ? CharPtr[i] - ((unsigned char) 'a' - (unsigned char) 'A') : CharPtr[i] );
		switch (value)
		{
			case (unsigned char) 'A':
			case (unsigned char) 'C':
			case (unsigned char) 'G':
			case (unsigned char) 'T':
			case (unsigned char) 'N':
				CharPtr[counter++] = value;
				break;
			default:
				;
		}
	}
	
	//for (i=2; i<counter; ++i) *CharPtr++ = TranslateDNACharToIndex(CharPtr);
	/* fill up missing phase data with unknown */
	for (i=0; i<2-(counter % 3); ++i) CharPtr[counter+i] = 0;
	
	Sequence->Length = counter/3;
	return Sequence;
}

unsigned char TranslateDNAChar3ToAminoIndex(const unsigned char * const triplet);
unsigned char TranslateDNAChar3ReversedToAminoIndex(const unsigned char * const triplet);

void FillInDNATripletToAminoAcidTable(const unsigned char * restrict const Alphabet);
void FillInDNATripletToAminoAcidTableComplementary(const unsigned char * restrict const Alphabet);
void FillInDNATripletToAminoAcidTableReversed(const unsigned char * restrict const Alphabet);
