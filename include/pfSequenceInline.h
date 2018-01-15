/***************************************************************************************************
                        PFTOOLS
 ***************************************************************************************************
  Oct 3, 2011 pfSequenceInline.h
 ***************************************************************************************************
 (C) 2011 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 ***************************************************************************************************/
#include "pfInput.h"

static inline PFSequence * ReadSequenceIndex(Sequence_t * const Seq, const size_t index, FILE * const stream, const s_Data * const DataPtr)
{
  /* Position into file */
  fseek(stream, (long) DataPtr[index].Offset, SEEK_SET);

  /* Read into memory */
  register const size_t Size = DataPtr[index].HeaderLength + DataPtr[index].SequenceLength + 1; // + 1 accounts for \n in between
  if (fread(Seq->Data.Header, sizeof(char), Size, stream) != Size) return NULL;

  /* Bound text */
  Seq->Data.Header[DataPtr[index].HeaderLength] = '\0';
  Seq->Data.Header[Size] = '\0';

  /* Set Sequence data start */
  Seq->ProfileData.ProfileIndex = (unsigned char*) &(Seq->Data.Header[DataPtr[index].HeaderLength+1]);

  /* Set Sequence length */
  Seq->ProfileData.Length = DataPtr[index].SequenceLength;

  return &(Seq->ProfileData);
}

#ifdef PRF_USE_MMAP
static inline PFSequence * MMAP_ReadSequenceIndex(Sequence_t * const Seq, const size_t index, const char * const restrict Array,
						 const s_Data * const DataPtr, const off_t InitialArrayOffset
#ifdef MMAP_DEBUG
	                                        ,const size_t ThreadId, const size_t NodeId, const size_t length
#endif						
)
{
  /* Position into Array */
  register const size_t ArrayOffset = (size_t) ( DataPtr[index].Offset - InitialArrayOffset);
//   fseek(stream, (long) DataPtr[index].Offset, SEEK_SET);

  /* Read into memory */
  register const size_t Size = DataPtr[index].HeaderLength + DataPtr[index].SequenceLength + 1; // + 1 accounts for \n in between
#ifdef MMAP_DEBUG
  if ( Size + ArrayOffset > length) {
    fprintf(stderr,"Thread %lu from Node %lu Seq: %lu will read beyond mmap %lu > %lu\n", ThreadId, NodeId, index, Size + ArrayOffset, length);
  }
#endif
  memcpy(Seq->Data.Header, &Array[ArrayOffset], sizeof(char)*Size);
  
//   if (fread(Seq->Data.Header, sizeof(char), Size, stream) != Size) return NULL;

  /* Bound text */
  Seq->Data.Header[DataPtr[index].HeaderLength] = '\0';
  Seq->Data.Header[Size] = '\0';

  /* Set Sequence data start */
  Seq->ProfileData.ProfileIndex = (unsigned char*) &(Seq->Data.Header[DataPtr[index].HeaderLength+1]);

  /* Set Sequence length */
  Seq->ProfileData.Length = DataPtr[index].SequenceLength;

  return &(Seq->ProfileData);
}
#endif

static inline void ReadSequenceNameIndex(char * const Name, const size_t index, FILE * const stream, const s_Data * const DataPtr)
{
  /* Position into file */
  fseek(stream, (long int) DataPtr[index].Offset, SEEK_SET);

  /* Read into memory */
  if (fscanf(stream, ">%s",Name) != 1) {
    fprintf(stderr, "Read error for sequence name %lu\n", index);
  }
}

extern inline char * __ALWAYS_INLINE CleanSequence(PFSequence * const Sequence)
{
  register size_t counter = 0;
  unsigned char * restrict const CharPtr = Sequence->ProfileIndex;

  for (size_t i=0; i<Sequence->Length; ++i) {
    register const unsigned char c = ((unsigned char) CharPtr[i] >= (unsigned char) 'a' ) 
                                    ? (unsigned char) CharPtr[i] - ((unsigned char) 'a' - (unsigned char) 'A')
				    : (unsigned char) CharPtr[i];
    if ( c >= (unsigned char) 'A' && c <= (unsigned char) 'Z' ) {
      CharPtr[counter++] = c;
    } 
  }
  
  Sequence->Length = counter;
  return (char*) Sequence->ProfileIndex;
}


extern inline unsigned char __ALWAYS_INLINE TranslateCharToIndex(const unsigned char letter, const unsigned char * restrict const Alphabet)
{
  const unsigned char lletter = (unsigned char) letter;
  register size_t index = (size_t) ( ( lletter >= (unsigned char) 'a' ) ? lletter - ((unsigned char) 'a' - (unsigned char) 'A') : lletter );
  if ( index >= (size_t) 'A' && index <= (size_t) 'Z' ) {
    return Alphabet[index - (size_t) 'A'];
  } else {
    return 0;
  }
}

/* WARNING: NEED TO OPTIMIZE THIS PIECE OF JUNK */
extern inline PFSequence * __ALWAYS_INLINE TranslateSequenceToIndex(PFSequence * const Sequence, const unsigned char * restrict const Alphabet )
{
  register size_t counter = 0;
  unsigned char * restrict const CharPtr = Sequence->ProfileIndex;

  for (size_t i=0; i<Sequence->Length; ++i) {
    register size_t index = (size_t) ( ( CharPtr[i] >= (unsigned char) 'a' ) ? CharPtr[i] - ((unsigned char) 'a' - (unsigned char) 'A') : CharPtr[i] );
    if ( index >= (size_t) 'A' && index <= (size_t) 'Z' ) {
      CharPtr[counter++] = Alphabet[index - (size_t) 'A'];
    }
  }

  Sequence->Length = counter;
  return Sequence;
}

extern inline void __ALWAYS_INLINE ReverseSequence(PFSequence * const Sequence)
{
  unsigned char * restrict const CharPtr = Sequence->ProfileIndex;
  const size_t SeqLength = Sequence->Length;
  unsigned char * BackPtr = &CharPtr[SeqLength-1];

  for (size_t i=0; i<SeqLength/2; ++i) {
      const unsigned char c = CharPtr[i];
      CharPtr[i] = *BackPtr;
      *BackPtr-- = c;
  }
  if (SeqLength & 0x1) {
      const unsigned char c = CharPtr[SeqLength/2];
      CharPtr[Sequence->Length/2] = *BackPtr;
      *BackPtr = c;
  }
}

extern inline void __ALWAYS_INLINE ReverseComplementSequence(PFSequence * const Sequence)
{
  unsigned char * restrict const CharPtr = Sequence->ProfileIndex;
  const size_t SeqLength = Sequence->Length;
  unsigned char * BackPtr = &CharPtr[SeqLength-1];

  for (size_t i=0; i<SeqLength/2; ++i) {
      unsigned char c;
			switch(CharPtr[i]) {
				case 'A': c = 'T'; break;
				case 'T': c = 'A'; break;
				case 'G': c = 'C'; break;
				case 'C': c = 'G'; break;
				default: c = '?';
			}
			unsigned char d;
			switch(*BackPtr) {
				case 'A': d = 'T'; break;
				case 'T': d = 'A'; break;
				case 'G': d = 'C'; break;
				case 'C': d = 'G'; break;
				default: d = '?';
			}
      CharPtr[i] = d;
      *BackPtr-- = c;
  }
  if (SeqLength & 0x1) {
			unsigned char c;
			switch(CharPtr[SeqLength/2]) {
				case 'A': c = 'T'; break;
				case 'T': c = 'A'; break;
				case 'G': c = 'C'; break;
				case 'C': c = 'G'; break;
				default: c = '?';
			}
			unsigned char d;
			switch(*BackPtr) {
				case 'A': d = 'T'; break;
				case 'T': d = 'A'; break;
				case 'G': d = 'C'; break;
				case 'C': d = 'G'; break;
				default: d = '?';
			}
      CharPtr[SeqLength/2] = d;
      *BackPtr = c;
  }
}

extern inline void __ALWAYS_INLINE PerformRevComp(unsigned char * const restrict Sequence, const size_t SeqLength)
{
	unsigned char * restrict const CharPtr = Sequence;
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

extern inline size_t __ALWAYS_INLINE CleanSequenceTo(unsigned char * const restrict To, const unsigned char * const restrict From, const size_t SeqLength)
{
	/* Remove carriage return */
	size_t counter = 0UL;
	for (size_t l=0; l<SeqLength; ++l) {
		register const char c = From[l];
		const register size_t index = (size_t) ( ( c >= (unsigned char) 'a' ) 
																? c - ((unsigned char) 'a' - (unsigned char) 'A')
																: c );
		if ( index >= (size_t) 'A' && index <= (size_t) 'Z' ) To[counter++] = c;
	}
	To[counter] = '\0';
	return counter;
}

extern inline void __ALWAYS_INLINE RevCompIndex(PFSequence * const Sequence, const unsigned char * restrict const Mapping)
{
	unsigned char * restrict const CharPtr = (unsigned char*) Sequence->ProfileIndex;
	const size_t SeqLength = Sequence->Length;
	unsigned char * BackPtr = &CharPtr[SeqLength-1];
	for (size_t l=0; l<SeqLength/2; ++l) {
			const unsigned char c = Mapping[CharPtr[l]];
			CharPtr[l] = Mapping[*BackPtr];
			*BackPtr-- = c;
	}
	if (SeqLength & 0x1) {
			const unsigned char c = Mapping[CharPtr[SeqLength/2]];
			CharPtr[SeqLength/2] = Mapping[*BackPtr];
			*BackPtr = c;
	}
}
