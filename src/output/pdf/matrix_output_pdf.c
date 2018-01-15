/*******************************************************
                        PFTOOLS
 *******************************************************
  Nov 3, 2015 pfplot_output_pdf.c
 *******************************************************
 (C) 2012-2015 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#include "prf_config.h"
#include "pfConfig.h"
#ifdef PRF_OUTPUT_PDF
#include <stdlib.h>
#include <stdio.h>
#include <smmintrin.h>
#include "pfCompute.h"
#include "pfOutput.h"
#include <hpdf.h>
#include <setjmp.h>

#define PlotAllSequence ((struct IO_PDF*)Options)->WholeSequence
enum StatePriority {
  PRIORITY_MATCH     = 1,
  PRIORITY_INSERTION = 2,
  PRIORITY_DELETION  = 0,
  PRIORITY_EXTRA     = 3
};


#define GET_STATE(x) (((x) & StateMask) >> StateShift)
#define TO_SCORE(x) GetScore((x), CoreCompute)
#define CYCLE_MASK 0x0


#ifdef HPDF_DLL
void  __stdcall
#else
void
#endif
error_handler  (HPDF_STATUS   error_no,
                HPDF_STATUS   detail_no,
                void         *user_data)
{
    printf ("ERROR: error_no=%04X, detail_no=%u\n", (HPDF_UINT) error_no, (HPDF_UINT) detail_no);
    longjmp(*((jmp_buf*) user_data), 1);
}
#define BORDER_SIZE 	10.0f
#define RIGHT_MARGIN 	50.0f
#define LEFT_MARGIN 	20.0f
#define TOP_MARGIN 	50.0f
#define BOTTOM_MARGIN 	50.0f
#define CELL_SIZE 	50.0f

#define SCORE_FONT_SIZE 	5
#define PLOT_FONT_SIZE 		5
#define SEQUENCE_FONT_SIZE 	8
#define MAX(x,y) ((x)>(y)?(x):(y))

static const HPDF_REAL cosinus270 = 0.0f; //cosf(-90.0f/180.0f*3.14159265f);
static const HPDF_REAL sinus270   = -1.0f; //sinf(-90.0f/180.0f*3.14159265f);
static const HPDF_REAL cosinus315 = 0.707106781f; //cosf(-45.0f/180.0f*3.14159265f);
static const HPDF_REAL sinus315   = 0.707106781f; //sinf(-45.0f/180.0f*3.14159265f);

void PDFOutput(union lScores * const restrict matrix,
               union lScores * const restrict rmatrix,
               const unsigned char * const SequenceText, const char * const Header,
               const struct Profile * const prf, const Compute_t * const restrict CoreCompute,
               const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock)
{
    HPDF_Doc  pdf;
    HPDF_Font font;
    HPDF_Page page;
		jmp_buf env;
    char fname[256] __attribute__((aligned(16)));
    char TextScore[64];
    char ComeFrom[4];
    const char * const restrict ProfileSequence = prf->Sequence;
    const size_t prfLength = prf->Length;
    register const int (* restrict Scores)[4];
    const size_t ld = GetMatrixLD(prf, CoreCompute);
		const size_t SeqID = 0UL;
    ComeFrom[PRIORITY_MATCH]     = 'M';
    ComeFrom[PRIORITY_INSERTION] = 'I';
    ComeFrom[PRIORITY_DELETION]  = 'D';
    ComeFrom[PRIORITY_EXTRA]     = 'X';
    
    if ( ((struct IO_PDF *)Options)->WithReverse && (rmatrix == NULL)) {
			fputs("Reverse cannot be done without reverse matrix!\n", stderr);
			((struct IO_PDF *)Options)->WithReverse = false;
		}
		
		const struct IO_PDF * const PDFOptions = Options;
    snprintf(fname, 255, "%s_seq%lu.pdf", PDFOptions->BaseFileName, SeqID);
    if (OutputVerbose) fprintf(stderr, "Working on %s having %lu characters\n", fname, SeqLength);

    pdf = HPDF_New (error_handler, (void*) &env);
    if (!pdf) {
        printf ("error: cannot create PdfDoc object\n");
        return;
    }

    if (setjmp(env)) {
        HPDF_Free (pdf);
        return;
    }

    /* create default-font */
    font = HPDF_GetFont (pdf, "Helvetica", NULL);

    /* add a new page object. */
    page = HPDF_AddPage (pdf);

    /*********************************************************************************************************************************
     *                                                         GET BEST ALIGNMENT                                                                  *********************************************************************************************************************************/
    Alignment_t lAlignment;
    const size_t Alignment_SeqLength = CoreCompute->GetBestAlignment(matrix, &lAlignment, SeqLength, prfLength);
		if (OutputVerbose) {
			fprintf(stderr, "Best alignment starts at (%i,%i) and ends at (%i,%i) with value %i, length %lu\n",
		    lAlignment.Region.sequence.Begin, lAlignment.Region.profile.Begin, 
		    lAlignment.Region.sequence.End, lAlignment.Region.profile.End, 
		    lAlignment.Score, Alignment_SeqLength);
		}
//     int iprf = lAlignment.Matrix.column.End;
//     int index_c = lAlignment.Matrix.row.End;
//     int Alignment_stop = lAlignment.Region.sequence.End;
    
    HPDF_REAL (*Alignment)[3] = (HPDF_REAL (*)[3]) malloc(3*(1+Alignment_SeqLength)*sizeof(HPDF_REAL));
    if (Alignment == NULL) {
			fputs("unable to allocate memory for alignment\n", stderr);
			goto END;
		}
    
    /*********************************************************************************************************************************
     *                                                      GET BEST REVERSED ALIGNMENT                                                            *********************************************************************************************************************************/
    Alignment_t lrAlignment;
    size_t rAlignment_SeqLength = 0U;
		HPDF_REAL (*rAlignment)[3] = NULL;
    if (PDFOptions->WithReverse) {
			rAlignment_SeqLength = CoreCompute->GetBestAlignment(rmatrix, &lrAlignment, SeqLength, prfLength);

      if (OutputVerbose) {
				fprintf(stderr, "Best reverse alignment starts at (%i,%i) and ends at (%i,%i) with value %i\n"
					"Alignment length %lu\n",
					lrAlignment.Region.profile.Begin, lrAlignment.Region.sequence.Begin, 
					lrAlignment.Region.profile.End, lrAlignment.Region.sequence.End, 
					lrAlignment.Score, rAlignment_SeqLength);
			}

			rAlignment = (HPDF_REAL (*)[3]) malloc(3*(1+rAlignment_SeqLength)*sizeof(HPDF_REAL));
      if (rAlignment == NULL) {
				fputs("unable to allocate memory for reverse alignment\n", stderr);
				free(rAlignment);
				goto END;
      }
    }

    /*********************************************************************************************************************************
     *                                                           GET BEST PATH                                                                *********************************************************************************************************************************/
    char * Alignment_Sequence=NULL, *Alignment_SequencePtr;
		if (PDFOptions->WithPaths) {}
		
    /*********************************************************************************************************************************
     *                                                         GET BEST REVERSED PATH                                                              *********************************************************************************************************************************/
    char * rAlignment_Sequence = NULL;
    if (PDFOptions->WithReverse && PDFOptions->WithPaths) { }
    
    /*********************************************************************************************************************************
     *                                                              PDF PAGE SETUP                                                                
     *********************************************************************************************************************************/
    
    /* Compute the real sequence region */
		size_t seqRegion_min = (size_t) lAlignment.Region.sequence.Begin;
    size_t seqRegion_max = (size_t) lAlignment.Region.sequence.End + 1LU;
    if (PDFOptions->WithReverse) {
      if (seqRegion_min > rAlignment[rAlignment_SeqLength-1][1]) seqRegion_min = rAlignment[rAlignment_SeqLength-1][1];
      if (seqRegion_max < rAlignment[0][1]) seqRegion_max = rAlignment[0][1];
    }
    seqRegion_min = (PlotAllSequence) ? 0 : seqRegion_min;
    seqRegion_max = (PlotAllSequence) ? SeqLength : seqRegion_max;
    
    size_t prfRegion_min = (size_t) 0;
    size_t prfRegion_max = (size_t) prfLength;
    
    if (PDFOptions->GivenRegion) {
			seqRegion_min = PDFOptions->y[0];
			seqRegion_max = PDFOptions->y[1];
			prfRegion_min = PDFOptions->x[0];
			prfRegion_max = PDFOptions->x[1];
    }
    
    if (seqRegion_max > SeqLength) seqRegion_max = SeqLength;
    if (prfRegion_max > ld) prfRegion_max = ld;
    
    if (OutputVerbose) {
			fprintf(stderr, "Plotting sequence region [%lu,%lu]\n", seqRegion_min, seqRegion_max);
	    fprintf(stderr, "Plotting profile region [%lu,%lu]\n", prfRegion_min, prfRegion_max);
		}
		const size_t region_size = (size_t) (seqRegion_max - seqRegion_min + 1UL);
    const size_t prfRegion_size = (size_t) (prfRegion_max - prfRegion_min + 1UL);
    
    /* Set to landscape mode */
    HPDF_REAL PageWidth = 2.0f*BORDER_SIZE + LEFT_MARGIN + RIGHT_MARGIN + CELL_SIZE*(prfRegion_size+1);
    HPDF_REAL PageHeight = 2.0f*BORDER_SIZE + TOP_MARGIN + BOTTOM_MARGIN + CELL_SIZE*region_size;
    if (OutputVerbose) fprintf(stderr, "PDF page size set to [%lf,%lf]\n", PageWidth, PageHeight);
    HPDF_Page_SetHeight(page, PageHeight);
    HPDF_Page_SetWidth(page, PageWidth);

    /* print the lines of the page. */
    HPDF_Page_SetLineWidth (page, 0.0);

    /* print the title of the page (with positioning center). */
    HPDF_Page_SetFontAndSize (page, font, 10);
    snprintf(TextScore, 64, "%s", PDFOptions->BaseFileName);
    HPDF_REAL tw = HPDF_Page_TextWidth (page, TextScore);
    HPDF_Page_BeginText (page);
    HPDF_Page_MoveTextPos (page, (HPDF_Page_GetWidth(page) - tw) / 2, HPDF_Page_GetHeight (page) - BORDER_SIZE - 5);
    HPDF_Page_ShowText (page, TextScore);
    HPDF_Page_EndText (page);

		if (PDFOptions->WithPaths) {
			HPDF_Page_BeginText (page);
			HPDF_Page_MoveTextPos (page, BORDER_SIZE, HPDF_Page_GetHeight (page) - BORDER_SIZE - 15);
			HPDF_Page_ShowText (page, "Standard alignment:");
			HPDF_Page_EndText (page);
			tw = HPDF_Page_TextWidth (page,"Standard alignment:");
			
			HPDF_Page_SetRGBFill (page, 1.0, 0.0, 0.0);
			HPDF_Page_BeginText (page);
			HPDF_Page_MoveTextPos (page, BORDER_SIZE + tw + 2, HPDF_Page_GetHeight (page) - BORDER_SIZE - 15);
			HPDF_Page_ShowText (page, Alignment_SequencePtr);
			HPDF_Page_EndText (page);
			HPDF_Page_SetRGBFill (page, 0.0, 1.0, 0.0);
		}
#if 0
    HPDF_Page_BeginText (page);
    HPDF_Page_MoveTextPos (page, BORDER_SIZE, HPDF_Page_GetHeight (page) - BORDER_SIZE - 27);
    HPDF_Page_ShowText (page, "Reversed alignment:");
    HPDF_Page_EndText (page);
    tw = HPDF_Page_TextWidth (page,"Reversed alignment:");

    HPDF_Page_BeginText (page);
    HPDF_Page_MoveTextPos (page, BORDER_SIZE + tw + 2, HPDF_Page_GetHeight (page) - BORDER_SIZE - 27);
    HPDF_Page_ShowText (page, rAlignment_Sequence);
    HPDF_Page_EndText (page);
    HPDF_Page_SetRGBFill (page, 0.0, 0.0, 0.0);
#endif
    HPDF_Page_SetFontAndSize (page, font, SCORE_FONT_SIZE);


    /* Compute the cell size */
    HPDF_REAL cell_width  = (HPDF_Page_GetWidth(page) - (HPDF_REAL) (2*BORDER_SIZE + LEFT_MARGIN + RIGHT_MARGIN) ) / ((HPDF_REAL) (prfRegion_size));
    HPDF_REAL cell_height = (HPDF_Page_GetHeight(page) - (HPDF_REAL) (2*BORDER_SIZE + TOP_MARGIN + BOTTOM_MARGIN)) / ((HPDF_REAL) region_size);
    if (OutputVerbose) fprintf(stderr, "Cell size: %lf x %lf\n", cell_width, cell_height);
    
    /*********************************************************************************************************************************
		 *                                                          BORDER SHADING                                                                     *********************************************************************************************************************************/
    {
      HPDF_Page_SetRGBFill (page, 0.9, 0.9, 0.9);
      
      // Left column
      if (prfRegion_min == 0) {
				HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN);
				HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN, BORDER_SIZE+BOTTOM_MARGIN);
				HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN + cell_width, BORDER_SIZE+BOTTOM_MARGIN);
				HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN + cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN);
				HPDF_Page_Fill(page);
      }
      // Right column
      if (PDFOptions->WithReverse && prfRegion_max >= prfLength) {
				HPDF_Page_MoveTo (page, HPDF_Page_GetWidth(page) - BORDER_SIZE - RIGHT_MARGIN, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN);
				HPDF_Page_LineTo (page, HPDF_Page_GetWidth(page) - BORDER_SIZE - RIGHT_MARGIN, BORDER_SIZE+BOTTOM_MARGIN);
				HPDF_Page_LineTo (page, HPDF_Page_GetWidth(page) - BORDER_SIZE - RIGHT_MARGIN - cell_width, BORDER_SIZE+BOTTOM_MARGIN);
				HPDF_Page_LineTo (page, HPDF_Page_GetWidth(page) - BORDER_SIZE - RIGHT_MARGIN - cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN);
				HPDF_Page_Fill(page);
      }
      
      // Potential first row
      if (seqRegion_min == 0) {
				HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN);
				HPDF_Page_LineTo (page, HPDF_Page_GetWidth(page) - BORDER_SIZE - RIGHT_MARGIN, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN);
				HPDF_Page_LineTo (page, HPDF_Page_GetWidth(page) - BORDER_SIZE - RIGHT_MARGIN, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN-cell_height);
				HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN-cell_height);
				HPDF_Page_Fill(page);
      }
      
      HPDF_Page_SetRGBFill (page, 0.0, 0.0, 0.0);
    }   
    /*********************************************************************************************************************************
     *                                                          SCORE SHADING                                                                      *********************************************************************************************************************************/
    {
      // Compute the maximum values for color scaling
      __m128i __max  = _mm_setzero_si128();      
      for (size_t i=0; i<ld*SeqLength; ++i) __max  = _mm_max_epi32(matrix[i].xmm, __max);
      if (PDFOptions->WithReverse) {
				__m128i __rmax = _mm_setzero_si128();
				for (size_t i=0; i<ld*SeqLength; ++i) {
						__rmax = _mm_max_epi32(rmatrix[i].xmm, __rmax);
				}
				__max  = _mm_max_epi32(__rmax, __max);
      }
#ifdef TAG
      // Clean extra bits
      __m128i __sign = _mm_cmplt_epi32(__max, _mm_setzero_si128());
      __max = _mm_srai_epi32(__max, 2);
      __max = _mm_or_si128(__max,_mm_and_si128(__sign, __ClearMask));
#endif
      const int ColorMax = _mm_extract_epi32(__max,MATCH);
      const int ColorMin = -20;
      if (OutputVerbose) fprintf(stderr, "Palette [%i,%i]\n", ColorMin, ColorMax);
      
      union lScores * restrict pmatrix = &matrix[ld];            
      for (size_t id_r=seqRegion_min>0?seqRegion_min:1; id_r<=seqRegion_max; ++id_r) {
				for (size_t id_c=(prfRegion_min<1)?1:prfRegion_min; id_c<=prfRegion_max; ++id_c) {
					const int value = TO_SCORE(pmatrix[id_c].Element[MATCH]);
					if (value>ColorMin) 
					{
						size_t column = id_c - prfRegion_min;
						size_t row    = id_r - seqRegion_min; 
						//HPDF_Page_SetRGBFill (page, 1.0, 1.0, 1.0 - (HPDF_REAL)(value-ColorMin)/((HPDF_REAL)(ColorMax-ColorMin)));
						
						HPDF_Page_SetRGBFill (page, 0.0, 1.0, 0.0);
						HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN + column*cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN - row*cell_height);
						HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN + (column+1)*cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN - row*cell_height);
						HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN + column*cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN - (row+1)*cell_height);
						HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN + column*cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN - row*cell_height);
						HPDF_Page_Fill(page);
						 if (OutputVerbose) fprintf(stderr, "Position %lu %lu : %lu %lu\t%i\n", id_r, id_c, row, column,  value);
						goto STOP;
					}
				}
				pmatrix += ld;
			}
      STOP:
#if 0
      size_t Nrow;
      size_t StartRow;
      union lScores * restrict pmatrix;
      if (seqRegion_min == 0) {
	Nrow = region_size-1;
	StartRow = 1;
	pmatrix = &matrix[ld];
      }
      else {
	Nrow = region_size;
	StartRow = 0;
	pmatrix = &matrix[seqRegion_min*ld];
      }
      
    for (size_t row=StartRow; row<=Nrow; ++row) {
	/* Avoid the first one as it is the border */
	++pmatrix;
	for (size_t column=1; column<=prfLength; ++column) {
	  const int value = TO_SCORE(pmatrix->Element[MATCH]);
#if !defined(_BEST_IS_NEGATIVE_)	  
	  if (value>ColorMin) 
#else
	  if (value<ColorMax)   
#endif
	  {
	    HPDF_Page_SetRGBFill (page, 1.0, 1.0, 1.0 - (HPDF_REAL)(value-ColorMin)/((HPDF_REAL)(ColorMax-ColorMin)));
	    HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN + column*cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN - row*cell_height);
	    HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN + (column+1)*cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN - row*cell_height);
	    HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN + column*cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN - (row+1)*cell_height);
	    HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN + column*cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN - row*cell_height);
	    HPDF_Page_Fill(page);
	  }
	  ++pmatrix;
	}
      }
      
//       if (PDFOptions->WithReverse) {
// 	const union lScores * restrict prmatrix = &rmatrix[/*SeqLength*/seqRegion_max*ld + prfLength];
// 	for (size_t row=StartRow; row<=Nrow; ++row) {
// 	  /* Avoid the first one as it is the border */ 
// 	    --prmatrix;
// 	    const int rvalue = TO_SCORE(prmatrix->Element[MATCH]);
// #if !defined(_BEST_IS_NEGATIVE_)
// 	    if (rvalue>ColorMin) 
// #else
// 	    if (rvalue<ColorMax) 
// #endif
// 	    {
// 	      HPDF_Page_SetRGBFill (page, 1.0, 1.0, 1.0 - (HPDF_REAL)(rvalue-ColorMin)/((HPDF_REAL)(ColorMax-ColorMin)));
// 	      HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN + column*cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN - (row+1)*cell_height);
// 	      HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN + (column+1)*cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN - (row+1)*cell_height);
// 	      HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN + (column+1)*cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN - row*cell_height);
// 	      HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN + column*cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN - (row+1)*cell_height);
// 	      HPDF_Page_Fill(page);
// 	    }
// 	    --prmatrix;
// 	  }
// 	}
//       }
#endif
      HPDF_Page_SetRGBFill (page, 0.0, 0.0, 0.0);
    }
    /*********************************************************************************************************************************
     *                                                               SQUARES                                                                       *********************************************************************************************************************************/
    {
      HPDF_Page_SetRGBStroke (page, 0.0, 0.0, 0.0);
      HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN,BORDER_SIZE+BOTTOM_MARGIN);
      HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN,HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN);
      HPDF_Page_Stroke (page);
      size_t column = 0;
      for (size_t id=prfRegion_min; id<=prfRegion_max; ++id, ++column) {
				HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN + column*cell_width, BORDER_SIZE + BOTTOM_MARGIN);
				HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN + column*cell_width, HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN);
				HPDF_Page_Stroke (page);

				TextScore[0] = ProfileSequence[id-1]; TextScore[1] = '\0';
				HPDF_Page_SetFontAndSize (page, font, SEQUENCE_FONT_SIZE);
				HPDF_REAL tw = HPDF_Page_TextWidth (page, TextScore);
				HPDF_Page_BeginText (page);
				HPDF_Page_MoveTextPos (page, BORDER_SIZE + LEFT_MARGIN + column*cell_width + 0.5f*(cell_width-tw),
						HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN + 3*SCORE_FONT_SIZE/2);
				HPDF_Page_ShowText (page, TextScore);
				HPDF_Page_EndText (page);

				HPDF_Page_SetFontAndSize (page, font, SCORE_FONT_SIZE);
				snprintf(TextScore,16,"%lu/%lu",id,prfLength-id+1);
				tw = HPDF_Page_TextWidth (page, TextScore);
				HPDF_Page_BeginText (page);
				HPDF_Page_MoveTextPos (page, BORDER_SIZE + LEFT_MARGIN + column*cell_width + 0.5f*(cell_width-tw),
						HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN + SCORE_FONT_SIZE/2);
				
				HPDF_Page_ShowText (page, TextScore);
				HPDF_Page_EndText (page);
      }
      HPDF_Page_MoveTo (page, HPDF_Page_GetWidth(page) - BORDER_SIZE - RIGHT_MARGIN - cell_width,BORDER_SIZE+BOTTOM_MARGIN);
      HPDF_Page_LineTo (page, HPDF_Page_GetWidth(page) - BORDER_SIZE - RIGHT_MARGIN - cell_width,HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN);
      HPDF_Page_Stroke (page);
      
      HPDF_Page_MoveTo (page, HPDF_Page_GetWidth(page) - BORDER_SIZE - RIGHT_MARGIN,BORDER_SIZE+BOTTOM_MARGIN);
      HPDF_Page_LineTo (page, HPDF_Page_GetWidth(page) - BORDER_SIZE - RIGHT_MARGIN,HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN);
      HPDF_Page_Stroke (page);

      // Sequence letters
      HPDF_REAL angle      = -90.0f/180.0f*3.14159265f;
      HPDF_REAL cosinus270 = cosf(angle);
      HPDF_REAL sinus270   = sinf(angle);
      size_t Endrow = region_size;
      if (seqRegion_min == 0) {
				Endrow--;
				HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN,                         BORDER_SIZE + BOTTOM_MARGIN + Endrow*cell_height);
				HPDF_Page_LineTo (page, HPDF_Page_GetWidth(page)-BORDER_SIZE-RIGHT_MARGIN, BORDER_SIZE + BOTTOM_MARGIN + Endrow*cell_height);
				HPDF_Page_Stroke (page);
      }
      const unsigned char * restrict Text = &SequenceText[seqRegion_max-1]; // removing the first row offset from score matrix
      for (size_t row=0; row<Endrow; ++row) {
				HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN,                         BORDER_SIZE + BOTTOM_MARGIN + row*cell_height);
				HPDF_Page_LineTo (page, HPDF_Page_GetWidth(page)-BORDER_SIZE-RIGHT_MARGIN, BORDER_SIZE + BOTTOM_MARGIN + row*cell_height);
				HPDF_Page_Stroke (page);

				HPDF_Page_SetFontAndSize (page, font, SEQUENCE_FONT_SIZE);
				HPDF_Page_BeginText (page);
				HPDF_Page_MoveTextPos (page, BORDER_SIZE, BORDER_SIZE + BOTTOM_MARGIN + (HPDF_REAL)(row+0.5)*cell_height/*+SEQUENCE_FONT_SIZE/2*/);
				TextScore[0] = Text[-row]; TextScore[1] = '\0';
				HPDF_Page_ShowText (page, TextScore);
				HPDF_Page_EndText (page);

				snprintf(TextScore,16,"%lu/%lu",seqRegion_max-row,1+SeqLength-(seqRegion_max-row));
				const HPDF_REAL tw = HPDF_Page_TextWidth (page, TextScore);

				HPDF_Page_SetFontAndSize (page, font, SCORE_FONT_SIZE);
				HPDF_Page_BeginText (page);
				HPDF_Page_SetTextMatrix (page, cosinus270, sinus270, -sinus270, cosinus270,
							/*x + 2*/ BORDER_SIZE+SEQUENCE_FONT_SIZE,
							/*y-cell_height+tw+1 + SCORE_FONT_SIZE*/ BORDER_SIZE + BOTTOM_MARGIN + (HPDF_REAL)(row+0.5)*cell_height+0.5f*tw-SCORE_FONT_SIZE/2);

			// 	HPDF_Page_MoveTextPos (page, , BORDER_SIZE + BOTTOM_MARGIN + (HPDF_REAL)(row+0.5)*cell_height-SCORE_FONT_SIZE/2);
				HPDF_Page_ShowText (page, TextScore);
				HPDF_Page_EndText (page);
      }
      HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN,HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN);
      HPDF_Page_LineTo (page, HPDF_Page_GetWidth(page) - BORDER_SIZE - RIGHT_MARGIN,HPDF_Page_GetHeight(page)-BORDER_SIZE-TOP_MARGIN);
      HPDF_Page_Stroke (page);
    }
    
    /*********************************************************************************************************************************
     *                                                                PATHS                                                                        *********************************************************************************************************************************/
    if (PDFOptions->WithPaths) {
      /* overlay the path */
      HPDF_Page_SetLineWidth (page, (HPDF_REAL) 3.0);
      HPDF_Page_SetRGBStroke (page, 1.0, 0.0, 0.0);
      HPDF_Page_SetRGBFill (page, 1.0, 0.0, 0.0);

      HPDF_Page_Circle(page, BORDER_SIZE + LEFT_MARGIN + (Alignment[1][0] + 0.5)*(cell_width),
			     BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - Alignment[1][1] + 0.5 )*cell_height, 8.0);
      HPDF_Page_Fill(page);
      if (Alignment_SeqLength>1) {
				HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN + (Alignment[0][0] + 0.5)*(cell_width),
							BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - Alignment[0][1] + 0.5 )*cell_height);
				for (size_t i=1; i<Alignment_SeqLength; ++i) {
					HPDF_Page_LineTo(page, BORDER_SIZE + LEFT_MARGIN + (Alignment[i][0] + 0.5)*(cell_width),
							BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - Alignment[i][1] + 0.5 )*cell_height);
				}
				HPDF_Page_Stroke(page);
      }
      HPDF_Page_MoveTo(page, BORDER_SIZE + LEFT_MARGIN + (Alignment[Alignment_SeqLength-1][0] + 0.5)*(cell_width) -5,
			     BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - Alignment[Alignment_SeqLength-1][1] + 0.5 )*cell_height - 5);
      HPDF_Page_LineTo(page, BORDER_SIZE + LEFT_MARGIN + (Alignment[Alignment_SeqLength-1][0] + 0.5)*(cell_width) +5,
			     BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - Alignment[Alignment_SeqLength-1][1] + 0.5 )*cell_height - 5);
      HPDF_Page_LineTo(page, BORDER_SIZE + LEFT_MARGIN + (Alignment[Alignment_SeqLength-1][0] + 0.5)*(cell_width) +5,
			     BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - Alignment[Alignment_SeqLength-1][1] + 0.5 )*cell_height + 5);
      HPDF_Page_LineTo(page, BORDER_SIZE + LEFT_MARGIN + (Alignment[Alignment_SeqLength-1][0] + 0.5)*(cell_width) -5,
			     BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - Alignment[Alignment_SeqLength-1][1] + 0.5 )*cell_height + 5);

      HPDF_Page_Fill(page);
      HPDF_Page_SetRGBStroke (page, 0.0, 0.0, 0.0);
      HPDF_Page_SetRGBFill (page, 0.0, 0.0, 0.0);


      /* Overlay the reverse path */
      if (PDFOptions->WithReverse) {
				HPDF_Page_SetLineWidth (page, (HPDF_REAL) 2.0);
				HPDF_Page_SetRGBStroke (page, 0.0, 1.0, 0.0);
				HPDF_Page_SetRGBFill (page, 0.0, 1.0, 0.0);
				HPDF_Page_Circle(page, BORDER_SIZE + LEFT_MARGIN + (rAlignment[1][0] + 0.5)*(cell_width),
									BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - rAlignment[1][1] + 0.5 )*cell_height, 4.0);
				HPDF_Page_Fill(page);

				if (rAlignment_SeqLength > 1) {
					HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN + (rAlignment[0][0] + 0.5)*(cell_width),
								BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - rAlignment[0][1] + 0.5 )*cell_height);

					for (size_t i=1; i<rAlignment_SeqLength; ++i) {
						HPDF_Page_LineTo(page, BORDER_SIZE + LEFT_MARGIN + (rAlignment[i][0] + 0.5)*(cell_width),
								BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - rAlignment[i][1] + 0.5 )*cell_height);
					}
					HPDF_Page_Stroke(page);
				}
				HPDF_Page_MoveTo(page, BORDER_SIZE + LEFT_MARGIN + (rAlignment[rAlignment_SeqLength-1][0] + 0.5)*(cell_width) -3,
									BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - rAlignment[rAlignment_SeqLength-1][1] + 0.5 )*cell_height - 3);
				HPDF_Page_LineTo(page, BORDER_SIZE + LEFT_MARGIN + (rAlignment[rAlignment_SeqLength-1][0] + 0.5)*(cell_width) +3,
									BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - rAlignment[rAlignment_SeqLength-1][1] + 0.5 )*cell_height - 3);
				HPDF_Page_LineTo(page, BORDER_SIZE + LEFT_MARGIN + (rAlignment[rAlignment_SeqLength-1][0] + 0.5)*(cell_width) +3,
									BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - rAlignment[rAlignment_SeqLength-1][1] + 0.5 )*cell_height + 3);
				HPDF_Page_LineTo(page, BORDER_SIZE + LEFT_MARGIN + (rAlignment[rAlignment_SeqLength-1][0] + 0.5)*(cell_width) -3,
									BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - rAlignment[rAlignment_SeqLength-1][1] + 0.5 )*cell_height + 3);
				HPDF_Page_Fill(page);
				HPDF_Page_SetRGBStroke (page, 0.0, 0.0, 0.0);
				HPDF_Page_SetRGBFill (page, 0.0, 0.0, 0.0);
      }
    }

    /*********************************************************************************************************************************
     *                                                           OVERLAY SCORES                                                                    *********************************************************************************************************************************/
    HPDF_REAL diag         = sqrtf(cell_height*cell_height + cell_width*cell_width);
    HPDF_REAL cosinus_diag = cell_width / diag;
    HPDF_REAL sinus_diag   = cell_height / diag;

    Scores = (const int (*)[4]) &matrix[seqRegion_min*ld+prfRegion_min];
 
    HPDF_REAL y = HPDF_Page_GetHeight(page) - BORDER_SIZE - TOP_MARGIN;

    const int MLOW = NLOW/4*3;
		const unsigned int StateShift = CoreCompute->StateShift;
		const unsigned int StateMask = CoreCompute->StateMask;
    for (size_t row=0; row<region_size; ++row) {
      HPDF_REAL x = BORDER_SIZE+LEFT_MARGIN;
      for (size_t column=0; column<prfRegion_size; ++column) {
				int value = TO_SCORE(Scores[column][DELETION]);
				int origin = ComeFrom[GET_STATE(Scores[column][DELETION])];
				if (value < NLOW) {
					snprintf(TextScore, 16, "NLOW %c", origin);
				}
				else if (value < MLOW) {
					snprintf(TextScore, 16, "MLOW %c", origin);
				}
				else {
					snprintf(TextScore, 16, "%i %c", value, origin);
				}
				if (Scores[column][DELETION] & CYCLE_MASK)
					HPDF_Page_SetRGBFill (page, 1.0, 0.0, 0.0);
				HPDF_REAL tw = HPDF_Page_TextWidth (page, TextScore);
				HPDF_Page_BeginText (page);
				HPDF_Page_MoveTextPos (page, x + (cell_width-tw) - 1 - SCORE_FONT_SIZE, y - 1 - SCORE_FONT_SIZE );
				HPDF_Page_ShowText (page, TextScore);
				HPDF_Page_EndText (page);
				if (Scores[column][DELETION] & CYCLE_MASK)
					HPDF_Page_SetRGBFill (page, 0.0, 0.0, 0.0);
				
				value = TO_SCORE(Scores[column][MATCH]);
				origin = ComeFrom[GET_STATE(Scores[column][MATCH])];
				if (value < NLOW) {
					snprintf(TextScore, 16, "NLOW %c", origin);
				}
				else if (value < MLOW) {
					snprintf(TextScore, 16, "MLOW %c", origin);
				}
				else {
					snprintf(TextScore, 16, "%i %c", value, origin);
				}

				if (Scores[column][MATCH] & CYCLE_MASK)
					HPDF_Page_SetRGBFill (page, 1.0, 0.0, 0.0);
				HPDF_Page_BeginText (page);
			// 	HPDF_Page_MoveTextPos (page, x + 1, y - 1 - SCORE_FONT_SIZE );
				HPDF_Page_SetTextMatrix (page, cosinus_diag, -sinus_diag, sinus_diag, cosinus_diag,
							x + SCORE_FONT_SIZE/2*cosinus_diag, y - SCORE_FONT_SIZE*sinus_diag);
				HPDF_Page_ShowText (page, TextScore);
				HPDF_Page_EndText (page);
				if (Scores[column][MATCH] & CYCLE_MASK)
					HPDF_Page_SetRGBFill (page, 0.0, 0.0, 0.0);
				
				value = TO_SCORE(Scores[column][INSERTION]);
				origin = ComeFrom[GET_STATE(Scores[column][INSERTION])];
				if (value < NLOW) {
					snprintf(TextScore, 16, "NLOW %c", origin);
				}
				else if (value < MLOW) {
					snprintf(TextScore, 16, "MLOW %c", origin);
				}
				else {
					snprintf(TextScore, 16, "%i %c", value, origin);
				}

				if (Scores[column][INSERTION] & CYCLE_MASK)
					HPDF_Page_SetRGBFill (page, 1.0, 0.0, 0.0);
				tw = HPDF_Page_TextWidth (page, TextScore);
				HPDF_Page_BeginText (page);
			// 	HPDF_Page_MoveTextPos (page, x + 1, y - cell_height + 1 );
				HPDF_Page_SetTextMatrix (page, cosinus270, sinus270, -sinus270, cosinus270, x + 2, y-cell_height+tw+1 + SCORE_FONT_SIZE);
				HPDF_Page_ShowText (page, TextScore);
				HPDF_Page_EndText (page);
				if (Scores[column][INSERTION] & CYCLE_MASK)
					HPDF_Page_SetRGBFill (page, 0.0, 0.0, 0.0);
				
				value = TO_SCORE(Scores[column][EXTRA]);
				origin = ComeFrom[GET_STATE(Scores[column][EXTRA])];
				if (value < NLOW) {
					snprintf(TextScore, 16, "NLOW %c", origin);
				}
				else if (value < MLOW) {
					snprintf(TextScore, 16, "MLOW %c", origin);
				}
				else {
					snprintf(TextScore, 16, "%i %c", value, origin);
				}
				if (Scores[column][EXTRA] & CYCLE_MASK)
					HPDF_Page_SetRGBFill (page, 1.0, 0.0, 0.0);
				tw = HPDF_Page_TextWidth (page, TextScore);
				HPDF_Page_BeginText (page);
			// 	HPDF_Page_MoveTextPos (page, x + (cell_width-tw) - 1, y - cell_height + 1 );
				HPDF_Page_SetTextMatrix (page, cosinus_diag, sinus_diag, -sinus_diag, cosinus_diag,
							x + 0.5*cell_width - 0.5*(tw+SCORE_FONT_SIZE)*cosinus_diag,
							y - 0.5*cell_height - 0.5*(tw-SCORE_FONT_SIZE)*sinus_diag);
				HPDF_Page_ShowText (page, TextScore);
				HPDF_Page_EndText (page);
				if (Scores[column][EXTRA] & CYCLE_MASK)
					HPDF_Page_SetRGBFill (page, 0.0, 0.0, 0.0);
				x += cell_width;
      }
      Scores += ld;
      y -= cell_height;
    }

    /*********************************************************************************************************************************
     *                                                      OVERLAY REVERSED SCORES                                                                *********************************************************************************************************************************/
    if (PDFOptions->WithReverse) {
      Scores = (const int (*)[4]) &rmatrix[(SeqLength-seqRegion_max+1)*ld];
      y = BORDER_SIZE + BOTTOM_MARGIN;
      for (int row=region_size-1; row>=0; row--) {
				HPDF_REAL x = HPDF_Page_GetWidth(page) - BORDER_SIZE - RIGHT_MARGIN - cell_width;
				for (int column=prfLength; column>=0; column--) {
					int value = TO_SCORE(Scores[0][DELETION]);
					if (value < -50000) {
						snprintf(TextScore, 16, "LOW %c",
							ComeFrom[Scores[0][DELETION] & 0x3]
							);
					}
					else {
						snprintf(TextScore, 16, "%i %c",
							value,
							ComeFrom[Scores[0][DELETION] & 0x3]
							);
					}
	
					HPDF_Page_BeginText (page);
					HPDF_Page_MoveTextPos (page, x + 1 + SCORE_FONT_SIZE, y + 1);
					HPDF_Page_ShowText (page, TextScore);
					HPDF_Page_EndText (page);

					value = TO_SCORE(Scores[0][MATCH]);
					if (value < -50000) {
						snprintf(TextScore, 16, "LOW %c",
							ComeFrom[Scores[0][MATCH] & 0x3]
							);
					}
					else {
						snprintf(TextScore, 16, "%i %c",
							value,
							ComeFrom[Scores[0][MATCH] & 0x3]
							);
					}
					HPDF_REAL tw = HPDF_Page_TextWidth (page, TextScore);
					HPDF_Page_BeginText (page);
				// 	HPDF_Page_MoveTextPos (page, x + 1, y - 1 - SCORE_FONT_SIZE );
					HPDF_Page_SetTextMatrix (page, cosinus_diag, -sinus_diag, sinus_diag, cosinus_diag,
								x + cell_width - tw - SCORE_FONT_SIZE/2*cosinus_diag, y + (SCORE_FONT_SIZE/2+tw)*sinus_diag);
					HPDF_Page_ShowText (page, TextScore);
					HPDF_Page_EndText (page);

					value = TO_SCORE(Scores[0][INSERTION]);
					if (value < -50000) {
						snprintf(TextScore, 16, "LOW %c",
							ComeFrom[Scores[0][INSERTION] & 0x3]
							);
					}
					else {
						snprintf(TextScore, 16, "%i %c",
							value,
							ComeFrom[Scores[0][INSERTION] & 0x3]
							);
					}
					HPDF_Page_BeginText (page);
				// 	HPDF_Page_MoveTextPos (page, x + 1, y - cell_height + 1 );
					HPDF_Page_SetTextMatrix (page, cosinus270, sinus270, -sinus270, cosinus270, x + cell_width - SCORE_FONT_SIZE, y + cell_height - SCORE_FONT_SIZE);
					HPDF_Page_ShowText (page, TextScore);
					HPDF_Page_EndText (page);

				// 	if (TO_SCORE(Scores[0][EXTRA]) >= 450 ) { fprintf(stderr," Look somewhere\n"); HPDF_Page_SetRGBFill(page, 1.0, 0.0, 0.0);}
					value = TO_SCORE(Scores[0][EXTRA]);
					if (value < -50000) {
						snprintf(TextScore, 16, "LOW %c",
							ComeFrom[Scores[0][EXTRA] & 0x3]
							);
					}
					else {
						snprintf(TextScore, 16, "%i %c",
							value,
							ComeFrom[Scores[0][EXTRA] & 0x3]
							);
					}
					tw = HPDF_Page_TextWidth (page, TextScore);
					HPDF_Page_BeginText (page);
				// 	HPDF_Page_MoveTextPos (page, x + (cell_width-tw) - 1, y - cell_height + 1 );
					HPDF_Page_SetTextMatrix (page, cosinus_diag, sinus_diag, -sinus_diag, cosinus_diag,
								x + 0.5*(cell_width  - (tw-2*SCORE_FONT_SIZE)*cosinus_diag),
								y + 0.5*(cell_height - (tw+2*SCORE_FONT_SIZE)*sinus_diag));
					HPDF_Page_ShowText (page, TextScore);
					HPDF_Page_EndText (page);
				// 	if (TO_SCORE(Scores[0][EXTRA]) >= 450 ) HPDF_Page_SetRGBFill(page, 0.0, 0.0, 0.0);

					++Scores;
					x -= cell_width;
				}
				y += cell_height;
      }
    }

    /*********************************************************************************************************************************
     *                                                      PLOT SCORES ALONG PATH                                                                 *********************************************************************************************************************************/
    /* get maxima and minima */
    HPDF_REAL minimum, maximum;
    minimum = Alignment[0][2];
    maximum = Alignment[0][2];
    if (OutputVerbose) fprintf(stderr, "Alignment %i : %f %f %f\n", 0, Alignment[0][0], Alignment[0][1], Alignment[0][2]);
    for (size_t i=1; i<Alignment_SeqLength; ++i) {
//       fprintf(stderr, "Alignment %lu : %f %f %f\n", i,Alignment[i][0], Alignment[i][1], Alignment[i][2]);
      minimum = (Alignment[i][2] < minimum) ? Alignment[i][2] : minimum;
      maximum = (Alignment[i][2] > maximum) ? Alignment[i][2] : maximum;
    }
    if (PDFOptions->WithReverse) {
      for (size_t i=0; i<rAlignment_SeqLength; ++i) {
				//       fprintf(stderr, "Alignment %lu : %f %f %f\n", i,rAlignment[i][0], rAlignment[i][1], rAlignment[i][2]);
				minimum = (rAlignment[i][2] < minimum) ? rAlignment[i][2] : minimum;
				maximum = (rAlignment[i][2] > maximum) ? rAlignment[i][2] : maximum;
      }
    }
    HPDF_REAL amplitude = maximum - minimum;
    if (OutputVerbose) fprintf(stderr, "Min: %f Max: %f Amplitude: %f\n", minimum, maximum, amplitude);

    /* plot profile scores */
    HPDF_Page_SetLineWidth (page, (HPDF_REAL) 0.0);
    HPDF_REAL x_offset = HPDF_Page_GetWidth(page) - 0.9*RIGHT_MARGIN - BORDER_SIZE;
    HPDF_REAL y_offset = HPDF_Page_GetHeight(page) - TOP_MARGIN - BORDER_SIZE;

    /* Y axis */
    if (PDFOptions->WithSequenceScoreGraph) {
      HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN, BORDER_SIZE);
      HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN, BORDER_SIZE + BOTTOM_MARGIN*0.9);
      HPDF_Page_Stroke(page);

      HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN, BORDER_SIZE + BOTTOM_MARGIN*0.9);
      HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN - 2, BORDER_SIZE + BOTTOM_MARGIN*0.9 - 6);
      HPDF_Page_LineTo (page, BORDER_SIZE + LEFT_MARGIN + 2, BORDER_SIZE + BOTTOM_MARGIN*0.9 - 6);
      HPDF_Page_Fill(page);
    }

    /* X axis */
    if (PDFOptions->WithProfileScoreGraph) {
      if ( maximum*minimum < 0.0) {
				HPDF_REAL ytmp = BORDER_SIZE - minimum/amplitude*BOTTOM_MARGIN*0.9;
				HPDF_REAL xtmp = HPDF_Page_GetWidth(page) - RIGHT_MARGIN - BORDER_SIZE;
				HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN, ytmp);
				HPDF_Page_LineTo (page, xtmp, ytmp);
				HPDF_Page_Stroke(page);

				HPDF_Page_MoveTo (page, xtmp, ytmp);
				HPDF_Page_LineTo (page, xtmp - 6, ytmp + 2);
				HPDF_Page_LineTo (page, xtmp - 6 , ytmp - 2);
				HPDF_Page_Fill(page);

				HPDF_Page_BeginText (page);
				HPDF_Page_MoveTextPos (page, BORDER_SIZE + LEFT_MARGIN - 5 , ytmp);
				HPDF_Page_ShowText (page, "0");
				HPDF_Page_EndText (page);

				HPDF_Page_BeginText (page);
				snprintf(TextScore, 16, "%i", (int) minimum);
				HPDF_Page_MoveTextPos (page, BORDER_SIZE + LEFT_MARGIN - 5 - HPDF_Page_TextWidth (page, TextScore) , BORDER_SIZE);
				HPDF_Page_ShowText (page, TextScore);
				HPDF_Page_EndText (page);

				HPDF_Page_BeginText (page);
				snprintf(TextScore, 16, "%i", (int) maximum);
				HPDF_Page_MoveTextPos (page, BORDER_SIZE + LEFT_MARGIN - 5 - HPDF_Page_TextWidth (page, TextScore) , BORDER_SIZE + 0.9*BOTTOM_MARGIN - PLOT_FONT_SIZE);
				HPDF_Page_ShowText (page, TextScore);
				HPDF_Page_EndText (page);
      }

      if (Alignment_SeqLength>1) {
				HPDF_Page_SetLineWidth (page, (HPDF_REAL) 1.0);

				HPDF_Page_SetRGBStroke (page, 1.0, 0.0, 0.0);
				HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN + (Alignment[0][0] + 0.5)*(cell_width),
							BORDER_SIZE + (Alignment[0][2] - minimum)/amplitude*BOTTOM_MARGIN*0.9);

				for (size_t i=1; i<Alignment_SeqLength; ++i) {
					HPDF_Page_LineTo(page, BORDER_SIZE + LEFT_MARGIN + (Alignment[i][0] + 0.5)*(cell_width),
							BORDER_SIZE + (Alignment[i][2] - minimum)/amplitude*BOTTOM_MARGIN*0.9);

				}
				HPDF_Page_Stroke(page);
				HPDF_Page_SetRGBStroke (page, 0.0, 0.0, 0.0);


				HPDF_Page_SetRGBFill (page, 0.0, 0.0, 0.0);
				for (size_t i=0; i<Alignment_SeqLength; ++i) {
					HPDF_Page_Circle(page, BORDER_SIZE + LEFT_MARGIN + (Alignment[i][0] + 0.5)*(cell_width),
							BORDER_SIZE + (Alignment[i][2] - minimum)/amplitude*BOTTOM_MARGIN*0.9, 1.0);
				}
				HPDF_Page_Fill(page);

				/* plot Sequence scores */
				HPDF_Page_SetLineWidth (page, (HPDF_REAL) 0.0);
      }

      if (PDFOptions->WithReverse) {
				if (rAlignment_SeqLength>1) {
					HPDF_Page_SetLineWidth (page, (HPDF_REAL) 1.0);

					HPDF_Page_SetRGBStroke (page, 0.0, 1.0, 0.0);
					HPDF_Page_MoveTo (page, BORDER_SIZE + LEFT_MARGIN + (rAlignment[0][0] + 0.5)*(cell_width),
								BORDER_SIZE + (rAlignment[0][2] - minimum)/amplitude*BOTTOM_MARGIN*0.9);

					for (size_t i=1; i<rAlignment_SeqLength; ++i) {
						HPDF_Page_LineTo(page, BORDER_SIZE + LEFT_MARGIN + (rAlignment[i][0] + 0.5)*(cell_width),
								BORDER_SIZE + (rAlignment[i][2] - minimum)/amplitude*BOTTOM_MARGIN*0.9);

					}
					HPDF_Page_Stroke(page);
					HPDF_Page_SetRGBStroke (page, 0.0, 0.0, 0.0);


					HPDF_Page_SetRGBFill (page, 0.0, 0.0, 0.0);
					for (size_t i=0; i<rAlignment_SeqLength; ++i) {
						HPDF_Page_Circle(page, BORDER_SIZE + LEFT_MARGIN + (rAlignment[i][0] + 0.5)*(cell_width),
								BORDER_SIZE + (rAlignment[i][2] - minimum)/amplitude*BOTTOM_MARGIN*0.9, 1.0);
					}
					HPDF_Page_Fill(page);

					/* plot Sequence scores */
					HPDF_Page_SetLineWidth (page, (HPDF_REAL) 0.0);
				}
      }

      /* X axis */
      HPDF_Page_MoveTo (page, x_offset, y_offset);
      HPDF_Page_LineTo (page, x_offset+0.9*RIGHT_MARGIN, y_offset);
      HPDF_Page_Stroke(page);

      HPDF_Page_MoveTo (page, x_offset+0.9*RIGHT_MARGIN,     y_offset);
      HPDF_Page_LineTo (page, x_offset+0.9*RIGHT_MARGIN - 6, y_offset + 2);
      HPDF_Page_LineTo (page, x_offset+0.9*RIGHT_MARGIN - 6, y_offset - 2);
      HPDF_Page_Fill(page);
    }

    if (PDFOptions->WithSequenceScoreGraph) {
      /* Y axis */
      if ( maximum*minimum < 0.0) {
				HPDF_REAL xtmp = x_offset - minimum/amplitude*RIGHT_MARGIN*0.9;
				HPDF_REAL ytmp = y_offset;
				HPDF_Page_MoveTo (page, xtmp, ytmp);
				HPDF_Page_LineTo (page, xtmp, BORDER_SIZE + BOTTOM_MARGIN);
				HPDF_Page_Stroke(page);

				HPDF_Page_MoveTo (page, xtmp,     BORDER_SIZE + BOTTOM_MARGIN);
				HPDF_Page_LineTo (page, xtmp - 2, BORDER_SIZE + BOTTOM_MARGIN + 6);
				HPDF_Page_LineTo (page, xtmp + 2, BORDER_SIZE + BOTTOM_MARGIN + 6);
				HPDF_Page_Fill(page);

				HPDF_Page_BeginText (page);
				HPDF_Page_MoveTextPos (page, xtmp , y_offset + 2);
				HPDF_Page_ShowText (page, "0");
				HPDF_Page_EndText (page);

				HPDF_Page_BeginText (page);
				snprintf(TextScore, 16, "%i", (int) minimum);
				HPDF_Page_MoveTextPos (page, HPDF_Page_GetWidth(page) - BORDER_SIZE - RIGHT_MARGIN,  y_offset + 3);
				HPDF_Page_ShowText (page, TextScore);
				HPDF_Page_EndText (page);

				HPDF_Page_BeginText (page);
				snprintf(TextScore, 16, "%i", (int) maximum);
				HPDF_Page_MoveTextPos (page, HPDF_Page_GetWidth(page) - BORDER_SIZE - HPDF_Page_TextWidth (page, TextScore) ,  y_offset + 3);
				HPDF_Page_ShowText (page, TextScore);
				HPDF_Page_EndText (page);
      }
      else {

      }

      if (Alignment_SeqLength>1) {
				HPDF_Page_SetLineWidth (page, (HPDF_REAL) 1.0);
				HPDF_Page_SetRGBStroke (page, 1.0, 0.0, 0.0);
				HPDF_Page_MoveTo (page, x_offset + (Alignment[0][2] - minimum)/amplitude*RIGHT_MARGIN*0.9,
							BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - Alignment[0][1] + 0.5)*(cell_height));

				for (size_t i=1; i<Alignment_SeqLength; ++i) {
					HPDF_Page_LineTo(page, x_offset + (Alignment[i][2] - minimum)/amplitude*RIGHT_MARGIN*0.9,
							BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - Alignment[i][1] + 0.5)*(cell_height));

				}
				HPDF_Page_Stroke(page);

				HPDF_Page_SetRGBStroke (page, 0.0, 0.0, 0.0);

				HPDF_Page_SetRGBFill (page, 0.0, 0.0, 0.0);
				for (size_t i=0; i<Alignment_SeqLength; ++i) {
					HPDF_Page_Circle(page, x_offset + (Alignment[i][2] - minimum)/amplitude*RIGHT_MARGIN*0.9,
							BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL) seqRegion_max - Alignment[i][1] + 0.5)*(cell_height), 1.0);
				}
				HPDF_Page_Fill(page);
      }
      if (PDFOptions->WithReverse) {
				if (rAlignment_SeqLength>1) {
					HPDF_Page_SetLineWidth (page, (HPDF_REAL) 1.0);
					HPDF_Page_SetRGBStroke (page, 0.0, 1.0, 0.0);
					HPDF_Page_MoveTo (page, x_offset + (rAlignment[0][2] - minimum)/amplitude*RIGHT_MARGIN*0.9,
								BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - rAlignment[0][1] + 0.5)*(cell_height));

					for (size_t i=1; i<rAlignment_SeqLength; ++i) {
						HPDF_Page_LineTo(page, x_offset + (rAlignment[i][2] - minimum)/amplitude*RIGHT_MARGIN*0.9,
								BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL)seqRegion_max - rAlignment[i][1] + 0.5)*(cell_height));

					}
					HPDF_Page_Stroke(page);

					HPDF_Page_SetRGBStroke (page, 0.0, 0.0, 0.0);

					HPDF_Page_SetRGBFill (page, 0.0, 0.0, 0.0);
					for (size_t i=0; i<rAlignment_SeqLength; ++i) {
						HPDF_Page_Circle(page, x_offset + (rAlignment[i][2] - minimum)/amplitude*RIGHT_MARGIN*0.9,
								BORDER_SIZE + BOTTOM_MARGIN + ((HPDF_REAL) seqRegion_max - rAlignment[i][1] + 0.5)*(cell_height), 1.0);
					}
					HPDF_Page_Fill(page);
				}
      }
    }


    /* save the document to a file */
    HPDF_SaveToFile (pdf, fname);

END:
    /* clean up */
    HPDF_Free (pdf);
    if (Alignment) free(Alignment);
		if (rAlignment) free(rAlignment);
		if (Alignment_Sequence) free(Alignment_Sequence);
}
#endif /* PRF_OUTPUT_PDF */
