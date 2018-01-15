/*******************************************************
                        PFTOOLS
 *******************************************************
  Apr 5, 2016 pfOutput.h
 *******************************************************
 (C) 2011-2016 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#ifndef _OUTPUT_FORMAT_H
#define _OUTPUT_FORMAT_H
#include <stdbool.h>
#include <pthread.h>
#include "pfConfig.h"
#include "pfProfile.h"
#include "pfCompute.h"

/*
 ************************************************************************************************
 *                                    DEFINITIONS                                               *
 ************************************************************************************************
 */
enum SearchType { PRF=1, REGEX=2, MAP=3 };
enum OutputType { TEXT, HISTOGRAM, DENSITY, PNG, DATA, PDF, TEST};
enum WrapperRange {
	Wrapper_Alignement = 0,
	Wrapper_Source = 1,
	Wrapper_Before = 2,
	Wrapper_After = 3,
	Wrapper_InBetween = 4
};
enum Constrain {
	WITH_REVERSE=1,
	SCORE_RANGE=2,
	CYCLE_RANGE=4,
	INVERSE_SELECTION=8,
	BORDER_CLIP=16
};


typedef void (*OutputMethod)(union lScores * const restrict matrix,
                             union lScores * const restrict rmatrix,
                             const unsigned char * const SequenceText, const char * const SequenceIdentification,
                             const struct Profile * const prf,
                             const Compute_t * const restrict CoreCompute,
                             const size_t SeqLength, const void * const Options,
                             pthread_mutex_t * const restrict PrintLock);

typedef void (*PrintFunction)(const struct Profile * const, const unsigned char * const,
                              const Alignment_t * const, const char * const, const size_t,
                              const float);

struct IO_Wrapper {
	PrintFunction Print;
  _Bool OptimalOnly;
	_Bool BestOfStdAndRevComp;
	enum WrapperRange Range;
};

struct IO_Data {
  const char * BaseFileName;
  _Bool IsBinary;
  _Bool Separate;
};

struct IO_Histogram {
	_Bool CycleRatherThanScore;
	const char *BaseFileName;
};
  
#ifdef PRF_OUTPUT_GRAPHICS 
struct IO_PNG { 
  const char * BaseFileName;
  int PixelWidth;
  int PixelHeight;
  int ScaleMin;
  int ScaleMax;
  int x[2];
  int y[2];
  int compress[2];
  _Bool GivenRegion;
  _Bool GivenScale;
  _Bool ScaleOnImage;
  _Bool WithReverse;
};
#endif

#ifdef PRF_OUTPUT_PDF
struct IO_PDF {
  const char * BaseFileName;
  _Bool WithPaths;
  _Bool WithProfileScoreGraph;
  _Bool WithSequenceScoreGraph;
  _Bool WithReverse;
  _Bool GivenRegion;
  _Bool WholeSequence;
  unsigned int x[2];
  unsigned int y[2];
};
#endif

typedef struct OutputType_s {
	OutputMethod OutputFct;
	int ScoreRange[2];
	unsigned short int CycleRange[2];
	unsigned short int BorderClip[2];
	enum OutputType Type;
	enum SearchType SearchWith;
	enum Constrain Constrains;
	union dummy {
		struct IO_Wrapper Text;
		struct IO_Data Data;
		struct IO_Histogram Histogram;
#ifdef PRF_OUTPUT_GRAPHICS
		struct IO_PNG PNG;
#endif
#ifdef PRF_OUTPUT_PDF
		struct IO_PDF PDF;
#endif
	} Specific;
} OutputType_t;

#define OUTPUT_TYPE_INIT {\
	.OutputFct = NULL,\
	.ScoreRange = { 0,0},\
	.CycleRange = {0U, 0U},\
	.BorderClip = {0U, 0U},\
	.Type = TEXT,\
	.SearchWith = 0,\
	.Constrains = 0,\
	.Specific = { .Text={\
		.Print = PrintDefault,\
		.OptimalOnly = false,\
		.BestOfStdAndRevComp = false,\
		.Range = Wrapper_Alignement\
		}\
	}\
}

#define INIT_TEXT(var) {\
	(var).Print = PrintDefault;\
	(var).OptimalOnly = false;\
	(var).BestOfStdAndRevComp = false;\
	(var).Range = Wrapper_Alignement;\
}

// struct IO_Data OptionsData = {
// 	.BaseFileName = OutputFileName,
// 	.IsBinary = false,
// 	.Separate = false
// };
// #ifdef PRF_OUTPUT_GRAPHICS
// struct IO_PNG  OptionsPNG  = {
// 	.BaseFileName = OutputFileName,
// 	.PixelWidth = 1,
// 	.PixelHeight = 1,
// 	.ScaleMin = 0,
// 	.ScaleMax = 0,
// 	.x = { 0, 0},
// 	.y = { 0, 0},
// 	.GivenRegion = false,
// 	.GivenScale = false,
// 	.ScaleOnImage = false,
// 	.WithReverse = false,
// 	.compress = {1,1}
// };
// #endif
// #ifdef PRF_OUTPUT_PDF
// struct IO_PDF OptionsPDF  = { 
// 	.BaseFileName = OutputFileName,
// 	.WithPaths = true, 
// 	.WithProfileScoreGraph = false,
// 	.WithSequenceScoreGraph = false,
// 	.WithReverse = false,
// 	.GivenRegion = false,
// 	.x = {0,0},
// 	.y = {0,0}
// };
// #endif


/*
 ************************************************************************************************
 *                                  VARIABLES REQUIRED                                          *
 ************************************************************************************************
 */
extern _Bool RegexCycleCount;
unsigned int OutputPrintWidth __attribute__((weak))= 80U;
_Bool OutputVerbose __attribute__((weak)) = false;

/*
 ************************************************************************************************
 *                                PRINT FUNCTION DECLARATIONS                                   *
 ************************************************************************************************
 */
#ifdef PRF_OUTPUT_FORMAT_SIMPLE
PFIMPEXP void PrintSimple(const struct Profile * const prf, const unsigned char * const AlignedSequence,
			  const Alignment_t * const alignment, const char * const Header,
			  const size_t SequenceLength, const float RAVE);
#endif
PFIMPEXP void PrintDefault(const struct Profile * const prf, const unsigned char * const AlignedSequence,
			   const Alignment_t * const alignment, const char * const Header,
			   const size_t SequenceLength, const float RAVE);
#if 0
PFIMPEXP void PrintPfscanLOpt(const struct Profile * const prf, const unsigned char * const AlignedSequence,
															const Alignment_t * const alignment, const char * const Header,
															const size_t SequenceLength, const float RAVE);
#endif
#ifdef PRF_OUTPUT_FORMAT_INTERPRO
PFIMPEXP void PrintInterpro(const struct Profile * const prf, const unsigned char * const AlignedSequence,
			    const Alignment_t * const alignment, const char * const Header,
			    const size_t SequenceLength, const float RAVE);
#endif
#ifdef PRF_OUTPUT_FORMAT_PFSCAN
PFIMPEXP void PrintPfscan(const struct Profile * const prf, const unsigned char * const AlignedSequence,
			  const Alignment_t * const alignment, const char * const Header,
			  const size_t SequenceLength, const float RAVE);
#endif
#ifdef PRF_OUTPUT_FORMAT_TSV
PFIMPEXP void PrintTSV(const struct Profile * const prf, const unsigned char * const AlignedSequence,
		       const Alignment_t * const alignment, const char * const Header,
		       const size_t SequenceLength, const float RAVE);
#endif
#ifdef PRF_OUTPUT_FORMAT_CLASSIFICATION
PFIMPEXP void PrintClassification(const struct Profile * const prf, const unsigned char * const AlignedSequence,
                         const struct Alignment * const alignment, const char * const Header,
                         const size_t SequenceLength, const float RAVE);
#endif
#ifdef PRF_OUTPUT_FORMAT_INCMATCH
PFIMPEXP void PrintIncMatch(const struct Profile * const prf, const unsigned char * const AlignedSequence,
			    const Alignment_t * const alignment, const char * const Header,
			    const size_t SequenceLength, const float RAVE);
#endif
#ifdef PRF_OUTPUT_FORMAT_PSMAKER
PFIMPEXP void PrintPSMaker(const struct Profile * const prf, const unsigned char * const AlignedSequence,
			   const Alignment_t * const alignment, const char * const Header,
			   const size_t SequenceLength, const float RAVE);
#endif
#ifdef PRF_OUTPUT_FORMAT_XPSA
PFIMPEXP void PrintxPSA( const struct Profile * const prf, const unsigned char * const AlignedSequence,
												 const Alignment_t * const alignment, const char * const Header,
												 const size_t SequenceLength, const float RAVE);
#endif
#ifdef PRF_OUTPUT_FORMAT_PSSCAN
PFIMPEXP void PrintPsScan(const struct Profile * const prf, const unsigned char * const AlignedSequence,
			  const Alignment_t * const alignment, const char * const Header,
			  const size_t SequenceLength, const float RAVE);
#endif
#ifdef PRF_OUTPUT_FORMAT_ONELINE
PFIMPEXP void PrintOneLine( const struct Profile * const prf, const unsigned char * const AlignedSequence,
			    const Alignment_t * const alignment, const char * const Header,
			    const size_t SequenceLength, const float RAVE);
#endif
#ifdef PRF_OUTPUT_FORMAT_FASEARCH
PFIMPEXP void Printfasearch(const struct Profile * const prf, const unsigned char * const AlignedSequence,
			   const Alignment_t * const alignment, const char * const Header,
			   const size_t SequenceLength, const float RAVE);
#endif
#ifdef PRF_OUTPUT_FORMAT_FASTA
PFIMPEXP void PrintFASTA(const struct Profile * const prf, const unsigned char * const AlignedSequence,
			   const Alignment_t * const alignment, const char * const Header,
			   const size_t SequenceLength, const float RAVE);
#endif
/*
 ************************************************************************************************
 *                                 OUTPUT FUNCTION DECLARATIONS                                 *
 ************************************************************************************************
 */

void WrapperOutput(union lScores * const restrict matrix,
                   union lScores * const restrict rmatrix,
                   const unsigned char * const SequenceText, const char * const SequenceIdentification,
                   const struct Profile * const prf, const Compute_t * const restrict CoreCompute,
                   const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock);

void TestOutput(union lScores * const restrict matrix,
                union lScores * const restrict rmatrix,
                const unsigned char * const SequenceText, const char * const SequenceIdentification,
                const struct Profile * const prf, const Compute_t * const restrict CoreCompute,
                const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock);

void DataOutput(union lScores * const restrict matrix,
                union lScores * const restrict rmatrix,
                const unsigned char * const SequenceText, const char * const SequenceIdentification,
                const struct Profile * const prf, const Compute_t * const restrict CoreCompute,
                const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock);

void TabOutput(union lScores * const restrict matrix, union lScores * const restrict rmatrix,
               const unsigned char * const SequenceText, const char * const Header,
               const struct Profile * const prf, const Compute_t * const restrict CoreCompute,
               const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock);

void RowOutput(union lScores * const restrict matrix, union lScores * const restrict rmatrix,
               const unsigned char * const SequenceText, const char * const Header,
               const struct Profile * const prf, const Compute_t * const restrict CoreCompute,
               const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock);

void ColumnOutput(union lScores * const restrict matrix, union lScores * const restrict rmatrix,
               const unsigned char * const SequenceText, const char * const Header,
               const struct Profile * const prf, const Compute_t * const restrict CoreCompute,
               const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock);

void HighlightedFASTAOutput(union lScores * const restrict matrix,
                            union lScores * const restrict rvmatrix,
                            const unsigned char * const SequenceText, const char * const restrict Header, 
                            const struct Profile * const prf, const Compute_t * const restrict CoreCompute,  
                            const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock);

#ifdef PRF_OUTPUT_PDF
void PDFOutput(union lScores * const restrict matrix,
               union lScores * const restrict rmatrix,
               const unsigned char * const SequenceText, const char * const SequenceIdentification,
               const struct Profile * const prf, const Compute_t * const restrict CoreCompute,
               const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock);
#endif

#ifdef PRF_OUTPUT_GRAPHICS
void PNGOutput(union lScores * const restrict matrix,
               union lScores * const restrict rmatrix,
               const unsigned char * const SequenceText, const char * const SequenceIdentification,
               const struct Profile * const prf, const Compute_t * const restrict CoreCompute,
               const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock);
#endif

int GetDispatchThreadIndex(const OutputType_t * const restrict output,
                           const struct Profile * const restrict prf);

/*
 ************************************************************************************************
 *                                   INLINE FUNCTIONS                                           *
 ************************************************************************************************
 */

#endif /* _OUTPUT_FORMAT_H */
