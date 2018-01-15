/*******************************************************
                        PFTOOLS
 *******************************************************
  May 29, 2013 pfSequence.h
 *******************************************************
 (C) 2013 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include "pfConfig.h"

#ifndef PFSEQ
#define PFSEQ
/*
 ************************************************************************************************
 *                                    DEFINITIONS                                               *
 ************************************************************************************************
 */
typedef struct PFSequence {
  unsigned char * ProfileIndex;
  size_t Length;
} PFSequence;
#endif

typedef struct Sequence {
  union {
    char * Header;
    void * Memory;
  } Data;
  size_t Size;
  PFSequence ProfileData;
} Sequence_t;

/*
 ************************************************************************************************
 *                             SEQUENCE FUNCTION DECLARATIONS                                   *
 ************************************************************************************************
 */
#ifndef __USE_INLINE_FUNCTIONS__
char * CleanSequence(PFSequence * const Sequence_t);
PFSequence * ReadSequenceIndex(Sequence_t * const Seq, const size_t index, FILE * const stream, const s_Data * const DataPtr);
PFSequence * MMAP_ReadSequenceIndex(Sequence_t * const Seq, const size_t index, const char * const restrict Array,
				    const s_Data * const DataPtr, const off_t InitialArrayOffset
#ifdef MMAP_DEBUG
				    ,const size_t ThreadId, const size_t NodeId, const size_t length
#endif						
);
void ReadSequenceNameIndex(char * const Name, const size_t index, FILE * const stream, const s_Data * const DataPtr);
unsigned char TranslateCharToIndex(const char letter, const unsigned char * restrict const Alphabet);
PFSequence * TranslateSequenceToIndex(PFSequence * const Sequence, const unsigned char * restrict const Alphabet );
void ReverseTranslatedSequence(PFSequence * const Sequence);
void ReverseComplementSequence(PFSequence * const Sequence);
#else
#include "pfSequenceInline.h"
#endif

/*
 ************************************************************************************************
 *                                   INLINE FUNCTIONS                                           *
 ************************************************************************************************
 */
#endif /* _SEQUENCE_H */
