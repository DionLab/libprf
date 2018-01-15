/*******************************************************
                        PFTOOLS
 *******************************************************
  Apr 5, 2016 pfInput.h
 *******************************************************
 (C) 2011-2016 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#ifndef _HANDLED_INPUT_H
#define _HANDLED_INPUT_H
#include <stdbool.h>
#include "pfConfig.h"

/*
 ************************************************************************************************
 *                                    DEFINITIONS                                               *
 ************************************************************************************************
 */
#ifdef PRF_INPUT_FASTA
enum FastaOptions {
	DoFastaIndexImport = 1,
	DoFastaIndexExport = 2
};

typedef struct s_Data {
	off_t Offset;
	unsigned int HeaderLength;
	unsigned int SequenceLength;
} s_Data;

typedef struct FASTAStructure {
	char FileName[256];
	const char * restrict SequenceFile;
	char * restrict indexFileName;
	s_Data  *DataPtr;
	off_t FileSize;
	size_t SequenceCount;
	size_t MaxSequenceSize; // this includes carriage returns and header
	enum FastaOptions Options;
} FASTAStructure;
#endif

/*
 ************************************************************************************************
 *                                        VARIABLES                                             *
 ************************************************************************************************
 */

/*
 ************************************************************************************************
 *                                FASTA FUNCTION DECLARATIONS                                   *
 ************************************************************************************************
 */
#ifdef PRF_INPUT_FASTA
int isFASTA(const char * const restrict File);
int OpenFASTAStructure(const char * const FileName, FASTAStructure * const FaS);
void CloseFASTAStructure(FASTAStructure * const FaS);
int AnalyzeFASTAStructure(const int FileDescriptor, FASTAStructure * const FaS);
int ExportFASTAStructure(FILE* const stream, const FASTAStructure * const Info);
int ImportFASTAStructure(FILE* const stream, FASTAStructure * const Info);
#endif

/*
 ************************************************************************************************
 *                                 PACBIO FUNCTION DECLARATIONS                                 *
 ************************************************************************************************
 */
#ifdef PRF_INPUT_HDF5
#include "pb_hdf5.h"
#endif
#ifdef PRF_INPUT_PBBAM
#include "pb_bam.h"
#endif

/*
 ************************************************************************************************
 *                                   INLINE FUNCTIONS                                           *
 ************************************************************************************************
 */

#endif /* _HANDLED_INPUT_H */
