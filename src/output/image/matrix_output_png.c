/*******************************************************
                        PFTOOLS
 *******************************************************
  May 7, 2012 pfplot_output_png.c
 *******************************************************
 (C) 2012-2015 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <smmintrin.h>
#include "pfCompute.h"
#include "pfOutput.h"

#include <gd.h>
#include <gdfonts.h>
// #define _LANDSCAPE_
#ifdef _LANDSCAPE_
#include "landscape_palette" 
#define PALETTE_SIZE 448
static void CreatePalette( gdImagePtr * const restrict Images, int (*colors)[PALETTE_SIZE] )
{
    for(size_t i = 0; i < PALETTE_SIZE; ++i) {
      const unsigned char r = (Landscape[i] >> 16 ) & 0xFF;
      const unsigned char g = (Landscape[i] >>  8 ) & 0xFF;
      const unsigned char b = (Landscape[i]       ) & 0xFF;
      
      colors[MATCH][i]     = gdImageColorAllocate(Images[MATCH],     r, g, b);
      colors[INSERTION][i] = gdImageColorAllocate(Images[INSERTION], r, g, b);
      colors[DELETION][i]  = gdImageColorAllocate(Images[DELETION],  r, g, b);
      colors[EXTRA][i]     = gdImageColorAllocate(Images[EXTRA],     r, g, b);
    }
}

// static const __m128 __CutoffOffset     = { 128.0f, 128.0f, 128.0f, 128.0f};
// static const __m128 __BelowCutoffRange = { 127.5f, 127.5f, 127.5f, 127.5f};
static const __m128 __CutoffOffset     = {  64.0f,  64.0f,  64.0f,  64.0f};
static const __m128 __BelowCutoffRange = {  63.5f,  63.5f,  63.5f,  63.5f};
static const __m128 __CutoffRange      = { 319.5f, 319.5f, 319.5f, 319.5f};
#else
#define PALETTE_SIZE 1280
static void CreatePalette( gdImagePtr * const restrict Images, int (*colors)[1280] )
{
    for(size_t i = 0; i < 256; ++i) {
      colors[MATCH][i]     = gdImageColorAllocate(Images[MATCH],     (255-i), (255-i), (255-i));
      colors[INSERTION][i] = gdImageColorAllocate(Images[INSERTION], (255-i), (255-i), (255-i));
      colors[DELETION][i]  = gdImageColorAllocate(Images[DELETION],  (255-i), (255-i), (255-i));
      colors[EXTRA][i]     = gdImageColorAllocate(Images[EXTRA],     (255-i), (255-i), (255-i));
      // Line 1: red = 0 ; green = 0 -> 255 ; blue = 255
      colors[MATCH][i+256]     = gdImageColorAllocate(Images[MATCH],     0,  i, 255);
      colors[INSERTION][i+256] = gdImageColorAllocate(Images[INSERTION], 0,  i, 255);
      colors[DELETION][i+256]  = gdImageColorAllocate(Images[DELETION],  0,  i, 255);
      colors[EXTRA][i+256]     = gdImageColorAllocate(Images[EXTRA],     0,  i, 255);
      // Line 2: red = 0 ; green = 255 ; blue = 255->0
      colors[MATCH][i+512]     = gdImageColorAllocate(Images[MATCH],     0, 255, (255-i));
      colors[INSERTION][i+512] = gdImageColorAllocate(Images[INSERTION], 0, 255, (255-i));
      colors[DELETION][i+512]  = gdImageColorAllocate(Images[DELETION],  0, 255, (255-i));
      colors[EXTRA][i+512]     = gdImageColorAllocate(Images[EXTRA],     0, 255, (255-i));
      // Line 3: red = 0->255 ; green = 255 ; blue = 0
      colors[MATCH][i+768]     = gdImageColorAllocate(Images[MATCH],     i, 255, 0);
      colors[INSERTION][i+768] = gdImageColorAllocate(Images[INSERTION], i, 255, 0);
      colors[DELETION][i+768]  = gdImageColorAllocate(Images[DELETION],  i, 255, 0);
      colors[EXTRA][i+768]     = gdImageColorAllocate(Images[EXTRA],     i, 255, 0);
      // Line 4: red = 255 ; green = 255 -> 0 ; blue = 0
      colors[MATCH][i+1024]     = gdImageColorAllocate(Images[MATCH],     255, (255-i), 0);
      colors[INSERTION][i+1024] = gdImageColorAllocate(Images[INSERTION], 255, (255-i), 0);
      colors[DELETION][i+1024]  = gdImageColorAllocate(Images[DELETION],  255, (255-i), 0);
      colors[EXTRA][i+1024]     = gdImageColorAllocate(Images[EXTRA],     255, (255-i), 0);
    }
}

static const __m128 __CutoffOffset     = {  256.0f,  256.0f,  256.0f,  256.0f};
static const __m128 __BelowCutoffRange = {  255.5f,  255.5f,  255.5f,  255.5f};
static const __m128 __CutoffRange      = { 1023.5f, 1023.5f, 1023.5f, 1023.5f};
#endif

#define _mm_neg_ps(x) _mm_xor_ps(x, _mm_set1_epi32(0x80000000))
/* This get the extrema values while getting rid of the extra information abaout coming from */
static void GetExtremumScores(__m128i * restrict Maximum, __m128i * restrict Minimum,
                              union lScores * const restrict matrix, 
															const Compute_t * const restrict CoreCompute,
															const size_t SeqLength,
                              const size_t matrix_lda)
{
  union lScores * restrict pmatrix = matrix;
  __m128i __Max = _mm_set1_epi32(NLOW);
  __m128i __Min = _mm_set1_epi32((-NLOW));
  register const __m128i __Zero = _mm_setzero_si128();
  register const __m128i __Left = _mm_set1_epi32(CoreCompute->LeftNegativeScoreMask);

  for (size_t iseq=0; iseq<=SeqLength; ++iseq) {
    for (size_t iprf=0; iprf<matrix_lda; ++iprf) {
      __m128i __mask = _mm_cmplt_epi32(pmatrix[iprf].xmm, __Zero);
      __mask = _mm_and_si128(__mask, __Left);
      const __m128i __value = _mm_or_si128( __mask, _mm_srai_epi32(pmatrix[iprf].xmm,CoreCompute->ScoreShift));
      __Max = _mm_max_epi32(__Max, __value);
      __Min = _mm_min_epi32(__Min, __value);
      _mm_store_si128(&pmatrix[iprf].xmm, __value);
    }
    pmatrix += matrix_lda;
  }
  *Maximum = __Max;
  *Minimum = __Min;
}

static void AddScores(union lScores * const matrix, const union lScores * const rmatrix,
                      const size_t SeqLength, const size_t PrfLength, const size_t matrix_lda,
                      __m128i * restrict Maximum, __m128i * restrict Minimum )
{
  register __m128i * restrict pmatrix  = &matrix[matrix_lda+1].xmm;
  const register __m128i * restrict prmatrix = &rmatrix[(1+SeqLength)*matrix_lda-1].xmm;
  
  if ( (Maximum != NULL) && (Minimum != NULL)) {
    __m128i __Max = _mm_set1_epi32(NLOW);
    __m128i __Min = _mm_set1_epi32(16384);
    for (size_t row=0; row<SeqLength; row++) { 
      for (size_t col=0; col<PrfLength; ++col) {
	__m128i __x      = *pmatrix; //_mm_srai_epi32(matrix[i].xmm,2);
	const __m128i __y = *prmatrix; //_mm_srai_epi32(pmatrix->xmm,2);
	__x = _mm_add_epi32(__x, __y);
	__Max = _mm_max_epi32(__Max, __x);
	__Min = _mm_min_epi32(__Min, __x);
	_mm_store_si128(pmatrix, __x);
	--prmatrix;
	++pmatrix;
      }
      pmatrix += 1;
      prmatrix -= 1;
    }
    *Maximum = __Max;
    *Minimum = __Min;
  }
  else {
    for (size_t row=0; row<SeqLength; row++) { 
      for (size_t col=0; col<PrfLength; ++col) {
	__m128i __x      = *pmatrix; //_mm_srai_epi32(matrix[i].xmm,2);
	const __m128i __y = *prmatrix; //_mm_srai_epi32(pmatrix->xmm,2);
	__x = _mm_add_epi32(__x, __y);
	_mm_store_si128(pmatrix, __x);
	--prmatrix;
	++pmatrix;
      }
      pmatrix += 1;
      prmatrix -= 1;
    }
  }
  if ( (uintptr_t) prmatrix != (uintptr_t) &rmatrix[matrix_lda-1]) {
      fprintf(stderr,"We do not end at beginning of rmatrix 0x%lx != 0x%lx\n%li 16 bytes vector differences\n",
	      (uintptr_t) prmatrix, (uintptr_t) &rmatrix[matrix_lda-1],
	      ((intptr_t) prmatrix - (intptr_t) &rmatrix[matrix_lda-1])/16);
  }
}


void PNGOutput(union lScores * const restrict matrix,
               union lScores * const restrict rmatrix,
               const unsigned char * const SequenceText, const char * const restrict Header,
               const struct Profile * const prf, const Compute_t * const restrict CoreCompute,
               const size_t SeqLength, const void * const Options, pthread_mutex_t * const restrict PrintLock)
{
  char OutputFileName[256] __attribute__((aligned(16)));
  char StdOutputFileName[256] __attribute__((aligned(16)));
  char ReverseOutputFileName[256] __attribute__((aligned(16)));
  char CombinedOutputFileName[256] __attribute__((aligned(16)));
  /* Declare color indexes */
  int colors[4][PALETTE_SIZE] __attribute__((aligned(16)));
  /* Declare the images */
  gdImagePtr im[4];

  union lScores Maximum[3] __attribute__((aligned(16)));
  union lScores Minimum[3] __attribute__((aligned(16)));
  union { __m128 xmmf; float Element[4]; } Scale;
  __m128i __ColorMinValue = _mm_set1_epi32(-200);
  const size_t prfLength = prf->Length;
  const size_t matrix_lda = prf->Length + 1;
	const size_t SeqID = 0U;

  union lScores * const Matrices[3] = { matrix, rmatrix, matrix};
  
  int ImageWidth, ImageHeight;
  int RealImageWidth, RealImageHeight;
  int ImageOffsetX, ImageOffsetY;
  
  if (!((struct IO_PNG*) Options)->GivenRegion) {
    ImageWidth   = (int) prfLength;
    ImageHeight  = (int) SeqLength;
    ImageOffsetX = 0;
    ImageOffsetY = 0;
    strncpy(StdOutputFileName, ((struct IO_PNG*) Options)->BaseFileName, 255);
    snprintf(ReverseOutputFileName, 255, "%s_reverse", ((struct IO_PNG*) Options)->BaseFileName);
    snprintf(CombinedOutputFileName, 255, "%s_combined", ((struct IO_PNG*) Options)->BaseFileName);
  }
  else {
    ImageWidth   = ((struct IO_PNG*) Options)->x[1] - ((struct IO_PNG*) Options)->x[0] + 1;
    ImageHeight  = ((struct IO_PNG*) Options)->y[1] - ((struct IO_PNG*) Options)->y[0] + 1;
    ImageOffsetX = ((struct IO_PNG*) Options)->x[0];
    ImageOffsetY = ((struct IO_PNG*) Options)->y[0];
    snprintf(StdOutputFileName, 255, "%s_[%i,%i][%i,%i]", ((struct IO_PNG*) Options)->BaseFileName,
	     ((struct IO_PNG*) Options)->x[0], ((struct IO_PNG*) Options)->x[1],
	     ((struct IO_PNG*) Options)->y[0], ((struct IO_PNG*) Options)->y[1]
	    );
    snprintf(ReverseOutputFileName, 255, "%s_[%i,%i][%i,%i]_reverse", ((struct IO_PNG*) Options)->BaseFileName,
	     ((struct IO_PNG*) Options)->x[0], ((struct IO_PNG*) Options)->x[1],
	     ((struct IO_PNG*) Options)->y[0], ((struct IO_PNG*) Options)->y[1]
	    );
    snprintf(CombinedOutputFileName, 255, "%s_[%i,%i][%i,%i]_combined", ((struct IO_PNG*) Options)->BaseFileName,
	     ((struct IO_PNG*) Options)->x[0], ((struct IO_PNG*) Options)->x[1],
	     ((struct IO_PNG*) Options)->y[0], ((struct IO_PNG*) Options)->y[1]
	     );
  }
  
  if (((struct IO_PNG*) Options)->compress[0] != 1 && ((struct IO_PNG*) Options)->compress[1] != 1) {
      RealImageWidth  = ImageWidth;
      RealImageHeight = ImageHeight;
      ImageWidth  = ImageWidth  / ((struct IO_PNG*) Options)->compress[0];
      ImageHeight = ImageHeight / ((struct IO_PNG*) Options)->compress[1];
  }
  
  const int PNG_CELL_X_SIZE = ((struct IO_PNG*) Options)->PixelWidth;
  const int PNG_CELL_Y_SIZE = ((struct IO_PNG*) Options)->PixelHeight;
  const int ScaledImageWidth  = ImageWidth*PNG_CELL_X_SIZE;
  const int ScaledImageHeight = ImageHeight*PNG_CELL_Y_SIZE; 
  const char *FileName[3] = { StdOutputFileName, ReverseOutputFileName, CombinedOutputFileName };
  const _Bool ScaleOnImage = ((struct IO_PNG*) Options)->ScaleOnImage;
  const int ScaleWidth = 100;
  const size_t upto = ( rmatrix == NULL) ? 1 : 3; 
  
  for (size_t mat=0; mat<upto; ++mat) {
    FILE *pngoutM, *pngoutI, *pngoutD,  *pngoutX;
    if (ScaleOnImage) {
      snprintf(OutputFileName, 255, "%s_match_%lu.png", FileName[mat], SeqID);
      pngoutM = fopen(OutputFileName, "wb");
      snprintf(OutputFileName, 255, "%s_insertion_%lu.png", FileName[mat], SeqID);
      pngoutI = fopen(OutputFileName, "wb");
      snprintf(OutputFileName, 255, "%s_deletion_%lu.png", FileName[mat], SeqID);
      pngoutD = fopen(OutputFileName, "wb");
      snprintf(OutputFileName, 255, "%s_extra_%lu.png", FileName[mat], SeqID);
      pngoutX = fopen(OutputFileName, "wb");
    } 
    else {
      snprintf(OutputFileName, 255, "%s_match_%lu_scale.png", FileName[mat], SeqID);
      pngoutM = fopen(OutputFileName, "wb");
      snprintf(OutputFileName, 255, "%s_insertion_%lu_scale.png", FileName[mat], SeqID);
      pngoutI = fopen(OutputFileName, "wb");
      snprintf(OutputFileName, 255, "%s_deletion_%lu_scale.png", FileName[mat], SeqID);
      pngoutD = fopen(OutputFileName, "wb");
      snprintf(OutputFileName, 255, "%s_extra_%lu_scale.png", FileName[mat], SeqID);
      pngoutX = fopen(OutputFileName, "wb");
    }
    if ( pngoutM == NULL || pngoutI == NULL || pngoutD == NULL || pngoutX == NULL ) {
      if (pngoutM != NULL) fclose(pngoutM);
      if (pngoutI != NULL) fclose(pngoutI);
      if (pngoutD != NULL) fclose(pngoutD);
      if (pngoutX != NULL) fclose(pngoutX);
      return;
    }

    /* Allocate the images */
    if (ScaleOnImage) {
      im[MATCH]     = gdImageCreateTrueColor(ScaleWidth + ScaledImageWidth, ScaledImageHeight);
      im[INSERTION] = gdImageCreateTrueColor(ScaleWidth + ScaledImageWidth, ScaledImageHeight);
      im[DELETION]  = gdImageCreateTrueColor(ScaleWidth + ScaledImageWidth, ScaledImageHeight);
      im[EXTRA]     = gdImageCreateTrueColor(ScaleWidth + ScaledImageWidth, ScaledImageHeight);
    }
    else {
      im[MATCH]     = gdImageCreateTrueColor(ScaleWidth, PALETTE_SIZE/2);
      im[INSERTION] = gdImageCreateTrueColor(ScaleWidth, PALETTE_SIZE/2);
      im[DELETION]  = gdImageCreateTrueColor(ScaleWidth, PALETTE_SIZE/2);
      im[EXTRA]     = gdImageCreateTrueColor(ScaleWidth, PALETTE_SIZE/2);
    }
    /*
     * Allocate the colors.
     */
    CreatePalette(im, colors);

    /* Add reverse to standard */
    if (mat < 2) {
      GetExtremumScores( &Maximum[mat].xmm, &Minimum[mat].xmm, Matrices[mat], CoreCompute, SeqLength, matrix_lda);
    }
    else {
      AddScores(matrix, rmatrix, SeqLength, prfLength, matrix_lda, &Maximum[mat].xmm, &Minimum[mat].xmm);
    }
    
    if (OutputVerbose) {
			fprintf(stderr, "Maxima: %i\t%i\t%i\t%i\n", Maximum[mat].Element[0], Maximum[mat].Element[1],
			        Maximum[mat].Element[2], Maximum[mat].Element[3]);
			fprintf(stderr, "Minima: %i\t%i\t%i\t%i\n", Minimum[mat].Element[0], Minimum[mat].Element[1],
			        Minimum[mat].Element[2], Minimum[mat].Element[3]);
    }
    
    /* Draw the scale */
    int lScaleWidth;

    if (ScaleOnImage) {
      int ScaleSize = PALETTE_SIZE;
      int colorstep = 1;
      {
				while (ScaleSize > ScaledImageHeight) { ScaleSize >>= 1; colorstep <<=1;}
				fprintf(stderr, "Scale size reduced to %i, color step of %i\n", ScaleSize, colorstep);
      }
      gdImageFilledRectangle(im[MATCH],     0, 0, ScaleWidth, ScaledImageHeight, colors[MATCH][0]);
      gdImageFilledRectangle(im[INSERTION], 0, 0, ScaleWidth, ScaledImageHeight, colors[INSERTION][0]);
      gdImageFilledRectangle(im[DELETION],  0, 0, ScaleWidth, ScaledImageHeight, colors[DELETION][0]);
      gdImageFilledRectangle(im[EXTRA],     0, 0, ScaleWidth, ScaledImageHeight, colors[EXTRA][0]);
      const int Bottom = (ScaledImageHeight + ScaleSize ) / 2;
      int lineY  = Bottom;
      for (int i=0; i<PALETTE_SIZE; i+=colorstep) {
				gdImageLine(im[MATCH],     20, lineY, 80, lineY, colors[MATCH][i]);
				gdImageLine(im[INSERTION], 20, lineY, 80, lineY, colors[INSERTION][i]);
				gdImageLine(im[DELETION],  20, lineY, 80, lineY, colors[DELETION][i]);
				gdImageLine(im[EXTRA],     20, lineY, 80, lineY, colors[EXTRA][i]);
				--lineY;
      }
      const gdFontPtr FontPtr = gdFontGetSmall();
      sprintf(OutputFileName, "%i", Minimum[mat].Element[MATCH]);
      gdImageString(im[MATCH], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[MATCH][255]);
      sprintf(OutputFileName, "%i", prf->CutOff);
      gdImageString(im[MATCH], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/colorstep - FontPtr->h/2, (unsigned char*)OutputFileName, colors[MATCH][0]);
      sprintf(OutputFileName, "%i", Maximum[mat].Element[MATCH]);
      gdImageString(im[MATCH], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), (int) SeqLength - Bottom + 2, (unsigned char*)OutputFileName, colors[MATCH][0]);

      sprintf(OutputFileName, "%i", Minimum[mat].Element[INSERTION]);
      gdImageString(im[INSERTION], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[INSERTION][255]);
      sprintf(OutputFileName, "%i", prf->CutOff);
      gdImageString(im[INSERTION], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/colorstep - FontPtr->h/2, (unsigned char*)OutputFileName, colors[INSERTION][0]);
      sprintf(OutputFileName, "%i", Maximum[mat].Element[INSERTION]);
      gdImageString(im[INSERTION], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), (int) SeqLength - Bottom + 2, (unsigned char*)OutputFileName, colors[INSERTION][0]);

      sprintf(OutputFileName, "%i", Minimum[mat].Element[DELETION]);
      gdImageString(im[DELETION], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[DELETION][255]);
      sprintf(OutputFileName, "%i", prf->CutOff);
      gdImageString(im[DELETION], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/colorstep - FontPtr->h/2, (unsigned char*)OutputFileName, colors[DELETION][0]);
      sprintf(OutputFileName, "%i", Maximum[mat].Element[DELETION]);
      gdImageString(im[DELETION], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), (int) SeqLength - Bottom + 2, (unsigned char*)OutputFileName, colors[DELETION][0]);

      sprintf(OutputFileName, "%i", Minimum[mat].Element[EXTRA]);
      gdImageString(im[EXTRA], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[EXTRA][255]);
      sprintf(OutputFileName, "%i", prf->CutOff);
      gdImageString(im[EXTRA], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/colorstep - FontPtr->h/2, (unsigned char*)OutputFileName, colors[EXTRA][0]);
      sprintf(OutputFileName, "%i", Maximum[mat].Element[EXTRA]);
      gdImageString(im[EXTRA], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), (int) SeqLength - Bottom + 2, (unsigned char*)OutputFileName, colors[EXTRA][0]);

      lScaleWidth = ScaleWidth;
    } 
    else {
      gdImageFilledRectangle(im[MATCH],     0, 0, ScaleWidth, PALETTE_SIZE/2, colors[MATCH][0]);
      gdImageFilledRectangle(im[INSERTION], 0, 0, ScaleWidth, PALETTE_SIZE/2, colors[INSERTION][0]);
      gdImageFilledRectangle(im[DELETION],  0, 0, ScaleWidth, PALETTE_SIZE/2, colors[DELETION][0]);
      gdImageFilledRectangle(im[EXTRA],     0, 0, ScaleWidth, PALETTE_SIZE/2, colors[EXTRA][0]);
      const int Bottom = PALETTE_SIZE/2 - 1;
      int lineY  = Bottom;
      for (int i=0; i<PALETTE_SIZE; i+=2) {
				gdImageLine(im[MATCH], 0, lineY, ScaleWidth, lineY, colors[MATCH][i]);
				gdImageLine(im[INSERTION], 0, lineY, ScaleWidth, lineY, colors[INSERTION][i]);
				gdImageLine(im[DELETION], 0, lineY, ScaleWidth, lineY, colors[DELETION][i]);
				gdImageLine(im[EXTRA], 0, lineY, ScaleWidth, lineY, colors[EXTRA][i]);
				--lineY;
      }
      const gdFontPtr FontPtr = gdFontGetSmall();
      sprintf(OutputFileName, "%i", Minimum[mat].Element[MATCH]);
      gdImageString(im[MATCH], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[MATCH][255]);
      sprintf(OutputFileName, "%i", prf->CutOff);
      gdImageString(im[MATCH], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/2 - FontPtr->h/2, (unsigned char*)OutputFileName, colors[MATCH][0]);
      sprintf(OutputFileName, "%i", Maximum[mat].Element[MATCH]);
      gdImageString(im[MATCH], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), 2, (unsigned char*)OutputFileName, colors[MATCH][0]);

      sprintf(OutputFileName, "%i", Minimum[mat].Element[INSERTION]);
      gdImageString(im[INSERTION], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[INSERTION][255]);
      sprintf(OutputFileName, "%i", prf->CutOff);
      gdImageString(im[INSERTION], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/2 - FontPtr->h/2, (unsigned char*)OutputFileName, colors[INSERTION][0]);
      sprintf(OutputFileName, "%i", Maximum[mat].Element[INSERTION]);
      gdImageString(im[INSERTION], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), 2, (unsigned char*)OutputFileName, colors[INSERTION][0]);

      sprintf(OutputFileName, "%i", Minimum[mat].Element[DELETION]);
      gdImageString(im[DELETION], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[DELETION][255]);
      sprintf(OutputFileName, "%i", prf->CutOff);
      gdImageString(im[DELETION], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/2 - FontPtr->h/2, (unsigned char*)OutputFileName, colors[DELETION][0]);
      sprintf(OutputFileName, "%i", Maximum[mat].Element[DELETION]);
      gdImageString(im[DELETION], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), 2, (unsigned char*)OutputFileName, colors[DELETION][0]);

      sprintf(OutputFileName, "%i", Minimum[mat].Element[EXTRA]);
      gdImageString(im[EXTRA], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[EXTRA][255]);
      sprintf(OutputFileName, "%i", prf->CutOff);
      gdImageString(im[EXTRA], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/2 - FontPtr->h/2, (unsigned char*)OutputFileName, colors[EXTRA][0]);
      sprintf(OutputFileName, "%i", Maximum[mat].Element[EXTRA]);
      gdImageString(im[EXTRA], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), 2, (unsigned char*)OutputFileName, colors[EXTRA][0]);

      /* Output the scale image to the disk file in PNG format. */
      gdImagePng(im[MATCH],     pngoutM);
      gdImagePng(im[INSERTION], pngoutI);
      gdImagePng(im[DELETION],  pngoutD);
      gdImagePng(im[EXTRA],     pngoutX);


      /* Close the files. */
      fclose(pngoutM);
      fclose(pngoutI);
      fclose(pngoutD);
      fclose(pngoutX);

      /* Destroy the image in memory. */
      gdImageDestroy(im[MATCH]);
      gdImageDestroy(im[INSERTION]);
      gdImageDestroy(im[DELETION]);
      gdImageDestroy(im[EXTRA]);

      /* Compose file names */
      snprintf(OutputFileName, 255, "%s_match_%lu.png", FileName[mat], SeqID);
      pngoutM = fopen(OutputFileName, "wb");
      snprintf(OutputFileName, 255, "%s_insertion_%lu.png", FileName[mat], SeqID);
      pngoutI = fopen(OutputFileName, "wb");
      snprintf(OutputFileName, 255, "%s_deletion_%lu.png", FileName[mat], SeqID);
      pngoutD = fopen(OutputFileName, "wb");
      snprintf(OutputFileName, 255, "%s_extra_%lu.png", FileName[mat], SeqID);
      pngoutX = fopen(OutputFileName, "wb");

      if ( pngoutM == NULL || pngoutI == NULL || pngoutD == NULL || pngoutX == NULL ) {
			if (pngoutM != NULL) fclose(pngoutM);
			if (pngoutI != NULL) fclose(pngoutI);
			if (pngoutD != NULL) fclose(pngoutD);
			if (pngoutX != NULL) fclose(pngoutX);
	return;
      }

      /* Allocate the images */
      im[MATCH]     = gdImageCreateTrueColor(ScaledImageWidth, ScaledImageHeight);
      im[INSERTION] = gdImageCreateTrueColor(ScaledImageWidth, ScaledImageHeight);
      im[DELETION]  = gdImageCreateTrueColor(ScaledImageWidth, ScaledImageHeight);
      im[EXTRA]     = gdImageCreateTrueColor(ScaledImageWidth, ScaledImageHeight);

      /* Allocate the colors */
      CreatePalette(im, colors);

      /* Correct shift */
      lScaleWidth = 0;
    }
    
    /* Compute scaling and offset */
    Maximum[mat].xmmf = _mm_cvtepi32_ps(Maximum[mat].xmm);
    Minimum[mat].xmmf = _mm_cvtepi32_ps(Minimum[mat].xmm);
    
    float cutoff = (float) prf->CutOff;
//     if (mat == 2) cutoff *= 2.0f;
    __m128 __CutoffValue    = _mm_set1_ps(cutoff);
    const __m128 __TooSmall = _mm_cmplt_ps(Maximum[mat].xmmf, __CutoffValue);
    __CutoffValue = _mm_blendv_ps(__CutoffValue, Maximum[mat].xmmf, __TooSmall);

    const __m128 __AboveCutoffScale = _mm_div_ps(__CutoffRange, _mm_sub_ps(Maximum[mat].xmmf, __CutoffValue));
    const __m128 __BelowCutoffScale = _mm_div_ps(__BelowCutoffRange, _mm_sub_ps(__CutoffValue, Minimum[mat].xmmf));

    /* Adjust Min/Max if required */
    if (((struct IO_PNG*) Options)->GivenScale) {
			Minimum[mat].xmm = _mm_set1_epi32(((struct IO_PNG*) Options)->ScaleMin);
			Maximum[mat].xmm = _mm_set1_epi32(((struct IO_PNG*) Options)->ScaleMax);
    }
    
    /* Draw the color boxes */
    const union lScores * restrict pmatrix = &((Matrices[mat])[(1+ImageOffsetY)*matrix_lda+3+ImageOffsetX]);
    int y = 0;
    
    if (((struct IO_PNG*) Options)->compress[0] == 1 && ((struct IO_PNG*) Options)->compress[1] == 1) {
      for (int iseq=0; iseq<ImageHeight; ++iseq) {
				const int y1 = y+PNG_CELL_Y_SIZE-1;
				int x = lScaleWidth;
				for (int iprf=0; iprf<ImageWidth; ++iprf) {
					/* Get the value in float format */
					__m128 __value = _mm_cvtepi32_ps(pmatrix[iprf].xmm);
					/* Set default for below cutoff value */
					__m128 __Scale  = __BelowCutoffScale;
					__m128 __Offset = _mm_setzero_ps();
					__m128 __Subs   = Minimum[mat].xmmf;
					/* Is it above cutoff ? update parameters*/
					__m128 __Mask = _mm_cmpgt_ps(__value, __CutoffValue);
					__Scale       = _mm_blendv_ps(__Scale, __AboveCutoffScale, __Mask);
					__Offset      = _mm_blendv_ps(__Offset, __CutoffOffset, __Mask);
					__Subs        = _mm_blendv_ps(__Subs, __CutoffValue, __Mask);

					__value = _mm_add_ps(__Offset,_mm_mul_ps(_mm_sub_ps(__value, __Subs), __Scale));
					union lScores pixelColor = { .xmm = _mm_cvtps_epi32(__value)};
					pixelColor.xmm = _mm_max_epi32(pixelColor.xmm, _mm_setzero_si128());
					pixelColor.xmm = _mm_min_epi32(pixelColor.xmm, _mm_set1_epi32(PALETTE_SIZE-1));
					
					gdImageFilledRectangle(im[MATCH],     x, y, x+PNG_CELL_X_SIZE-1, y1, colors[MATCH][pixelColor.Element[MATCH]]);
					gdImageFilledRectangle(im[INSERTION], x, y, x+PNG_CELL_X_SIZE-1, y1, colors[INSERTION][pixelColor.Element[INSERTION]]);
					gdImageFilledRectangle(im[DELETION],  x, y, x+PNG_CELL_X_SIZE-1, y1, colors[DELETION][pixelColor.Element[DELETION]]);
					gdImageFilledRectangle(im[EXTRA],     x, y, x+PNG_CELL_X_SIZE-1, y1, colors[EXTRA][pixelColor.Element[EXTRA]]);
					x += PNG_CELL_X_SIZE;
				}
				y += PNG_CELL_Y_SIZE;
				pmatrix += matrix_lda;
      }
    }
    else {
      /* Compression is active */
      const int compression_x = ((struct IO_PNG*) Options)->compress[0];
      const int compression_y = ((struct IO_PNG*) Options)->compress[1];
      const int yStep = compression_y*matrix_lda;
      for (int iseq=0; iseq<RealImageHeight; iseq += compression_y) {
				const int y1 = y+PNG_CELL_Y_SIZE-1;
				int x = lScaleWidth;
				for (int iprf=0; iprf<RealImageWidth; iprf += compression_x) {
					/* Get the value in float format */
					int xMax = compression_x;
					if (compression_x + iprf >= RealImageWidth) xMax = RealImageWidth - iprf + 1;
					int yMax = compression_y;
					if (compression_y + iseq >= RealImageHeight) yMax = RealImageHeight - iseq + 1;
					const __m128i * restrict ptr = &(pmatrix[iprf].xmm);
					__m128i __valueI = _mm_set1_epi32(NLOW);
					for (int icompress=0; icompress<yMax; icompress++) {
						for (int jcompress=0; jcompress<xMax; jcompress++) {
							__valueI = _mm_max_epi32(__valueI, ptr[jcompress]);
						}
						ptr += matrix_lda; 
					}
					
					__m128 __value = _mm_cvtepi32_ps(__valueI);
					
					/* Set default for below cutoff value */
					__m128 __Scale  = __BelowCutoffScale;
					__m128 __Offset = _mm_setzero_ps();
					__m128 __Subs   = Minimum[mat].xmmf;
					/* Is it above cutoff ? update parameters*/
					__m128 __Mask = _mm_cmpgt_ps(__value, __CutoffValue);
					__Scale       = _mm_blendv_ps(__Scale, __AboveCutoffScale, __Mask);
					__Offset      = _mm_blendv_ps(__Offset, __CutoffOffset, __Mask);
					__Subs        = _mm_blendv_ps(__Subs, __CutoffValue, __Mask);

					__value = _mm_add_ps(__Offset,_mm_mul_ps(_mm_sub_ps(__value, __Subs), __Scale));
					union lScores pixelColor = { .xmm = _mm_cvtps_epi32(__value)};
					pixelColor.xmm = _mm_max_epi32(pixelColor.xmm, _mm_setzero_si128());
					pixelColor.xmm = _mm_min_epi32(pixelColor.xmm, _mm_set1_epi32(PALETTE_SIZE-1));
					
					gdImageFilledRectangle(im[MATCH],     x, y, x+PNG_CELL_X_SIZE-1, y1, colors[MATCH][pixelColor.Element[MATCH]]);
					gdImageFilledRectangle(im[INSERTION], x, y, x+PNG_CELL_X_SIZE-1, y1, colors[INSERTION][pixelColor.Element[INSERTION]]);
					gdImageFilledRectangle(im[DELETION],  x, y, x+PNG_CELL_X_SIZE-1, y1, colors[DELETION][pixelColor.Element[DELETION]]);
					gdImageFilledRectangle(im[EXTRA],     x, y, x+PNG_CELL_X_SIZE-1, y1, colors[EXTRA][pixelColor.Element[EXTRA]]);
					x += PNG_CELL_X_SIZE;
				}
				y += PNG_CELL_Y_SIZE;
				pmatrix += yStep;
      }
    }

    /* Output the image to the disk file in PNG format. */
    gdImagePng(im[MATCH],     pngoutM);
    gdImagePng(im[INSERTION], pngoutI);
    gdImagePng(im[DELETION],  pngoutD);
    gdImagePng(im[EXTRA],     pngoutX);


    /* Close the files. */
    fclose(pngoutM);
    fclose(pngoutI);
    fclose(pngoutD);
    fclose(pngoutX);

    /* Destroy the image in memory. */
    gdImageDestroy(im[MATCH]);
    gdImageDestroy(im[INSERTION]);
    gdImageDestroy(im[DELETION]);
    gdImageDestroy(im[EXTRA]);
  }
}


// static void GetExtremumScores(__m128i * restrict Maximum, __m128i * restrict Minimum,
// 			      union lScores * const matrix, const size_t SeqLength, const size_t prfLength)
// {
//   union lScores * restrict pmatrix = matrix;
//   __m128i __Max = _mm_set1_epi32(NLOW);
//   __m128i __Min = _mm_set1_epi32(16384);
//   register const __m128i __Zero = _mm_setzero_si128();
//   register const __m128i __Left = _mm_set1_epi32(0xC0000000);
// 
//   for (size_t iseq=0; iseq<=SeqLength; ++iseq) {
//     for (size_t iprf=0; iprf<=prfLength; ++iprf) {
//       __m128i __mask = _mm_cmplt_epi32(pmatrix[iprf].xmm, __Zero);
//       __mask = _mm_and_si128(__mask, __Left);
//       const __m128i __value = _mm_or_si128( __mask, _mm_srai_epi32(pmatrix[iprf].xmm,2));
//       __Max = _mm_max_epi32(__Max, __value);
//       __Min = _mm_min_epi32(__Min, __value);
//       _mm_store_si128(&pmatrix[iprf].xmm, __value);
//     }
//     pmatrix += prfLength+1;
//   }
//   *Maximum = __Max;
//   *Minimum = __Min;
// }
// 
// static void AddScores(union lScores * const matrix, const union lScores * const rmatrix,
// 		      const size_t SeqLength, const size_t prfLength,
// 		      __m128i * restrict Maximum, __m128i * restrict Minimum )
// {
//   const size_t fullsize = (1+SeqLength)*(1+prfLength);
//   const union lScores * restrict pmatrix = &rmatrix[fullsize];
//   __m128i __Max = _mm_set1_epi32(NLOW);
//   __m128i __Min = _mm_set1_epi32(16384);
// 
//   for (size_t i=0; i<fullsize; ++i) {
//     __m128i __x      = matrix[i].xmm; //_mm_srai_epi32(matrix[i].xmm,2);
//     const __m128i __y = pmatrix->xmm; //_mm_srai_epi32(pmatrix->xmm,2);
//     __x = _mm_add_epi32(__x, __y);
//     __Max = _mm_max_epi32(__Max, __x);
//     __Min = _mm_min_epi32(__Min, __x);
//     _mm_store_si128(&matrix[i].xmm, __x);
//     --pmatrix;
//   }
//   if ( (size_t) pmatrix != (size_t) rmatrix) {
//       fputs("We do not end at beginning of rmatrix\n", stderr);
//       exit(1);
//   }
//   *Maximum = __Max;
//   *Minimum = __Min;
// }
// 
// 
// 
// void PNGOutput(union lScores * const matrix, union lScores * const rmatrix,
// 	       union lScores * const restrict LengthMatrix, union lScores * const restrict rLengthMatrix,
// 	       const char * const SequenceText, const struct Profile * const prf, const char * const BaseFileName,
// 	       const size_t SeqLength, const void * const Options, const size_t SeqID)
// {
//   char OutputFileName[256] __attribute__((aligned(16)));
//   char StdOutputFileName[256] __attribute__((aligned(16)));
//   char ReverseOutputFileName[256] __attribute__((aligned(16)));
//   char CombinedOutputFileName[256] __attribute__((aligned(16)));
//   /* Declare color indexes */
//   int colors[4][PALETTE_SIZE] __attribute__((aligned(16)));
//   /* Declare the images */
//   gdImagePtr im[4];
// 
//   union lScores Maximum[3] __attribute__((aligned(16)));
//   union lScores Minimum[3] __attribute__((aligned(16)));
//   union { __m128 xmmf; float Element[4]; } Scale;
//   __m128i __ColorMinValue = _mm_set1_epi32(-200);
//   const size_t prfLength = prf->Length;
//   const size_t matrix_lda = prfLength + 1;
// 
//   union lScores * const Matrices[3] = { matrix, rmatrix, matrix};
//   
//   int ImageWidth, ImageHeight;
//   int ImageOffsetX, ImageOffsetY;
//   
//   if (!((struct IO_PNG*) Options)->GivenRegion) {
//     ImageWidth   = (int) prfLength;
//     ImageHeight  = (int) SeqLength;
//     ImageOffsetX = 0;
//     ImageOffsetY = 0;
//     strncpy(StdOutputFileName, BaseFileName, 255);
//     snprintf(ReverseOutputFileName, 255, "%s_reverse", BaseFileName);
//     snprintf(CombinedOutputFileName, 255, "%s_combined", BaseFileName);
//   }
//   else {
//     ImageWidth   = ((struct IO_PNG*) Options)->x[1] - ((struct IO_PNG*) Options)->x[0] + 1;
//     ImageHeight  = ((struct IO_PNG*) Options)->y[1] - ((struct IO_PNG*) Options)->y[0] + 1;
//     ImageOffsetX = ((struct IO_PNG*) Options)->x[0];
//     ImageOffsetY = ((struct IO_PNG*) Options)->y[0];
//     snprintf(StdOutputFileName, 255, "%s_[%i,%i][%i,%i]", BaseFileName,
// 	     ((struct IO_PNG*) Options)->x[0], ((struct IO_PNG*) Options)->x[1],
// 	     ((struct IO_PNG*) Options)->y[0], ((struct IO_PNG*) Options)->y[1]
// 	    );
//     snprintf(ReverseOutputFileName, 255, "%s_[%i,%i][%i,%i]_reverse", BaseFileName,
// 	     ((struct IO_PNG*) Options)->x[0], ((struct IO_PNG*) Options)->x[1],
// 	     ((struct IO_PNG*) Options)->y[0], ((struct IO_PNG*) Options)->y[1]
// 	    );
//     snprintf(CombinedOutputFileName, 255, "%s_[%i,%i][%i,%i]_combined", BaseFileName,
// 	     ((struct IO_PNG*) Options)->x[0], ((struct IO_PNG*) Options)->x[1],
// 	     ((struct IO_PNG*) Options)->y[0], ((struct IO_PNG*) Options)->y[1]
// 	     );
//   }
//   
//   const int PNG_CELL_X_SIZE = ((struct IO_PNG*) Options)->PixelWidth;
//   const int PNG_CELL_Y_SIZE = ((struct IO_PNG*) Options)->PixelHeight;
//   const int ScaledImageWidth  = ImageWidth*PNG_CELL_X_SIZE;
//   const int ScaledImageHeight = ImageHeight*PNG_CELL_Y_SIZE; 
//   const char *FileName[3] = { StdOutputFileName, ReverseOutputFileName, CombinedOutputFileName };
//   const _Bool ScaleOnImage = ((struct IO_PNG*) Options)->ScaleOnImage;
//   const int ScaleWidth = 100;
//   const size_t upto = ( rmatrix == NULL) ? 1 : 3; 
//   
//   for (size_t mat=0; mat<upto; ++mat) {
//     FILE *pngoutM, *pngoutI, *pngoutD,  *pngoutX;
//     if (ScaleOnImage) {
//       snprintf(OutputFileName, 255, "%s_match_%lu.png", FileName[mat], SeqID);
//       pngoutM = fopen(OutputFileName, "wb");
//       snprintf(OutputFileName, 255, "%s_insertion_%lu.png", FileName[mat], SeqID);
//       pngoutI = fopen(OutputFileName, "wb");
//       snprintf(OutputFileName, 255, "%s_deletion_%lu.png", FileName[mat], SeqID);
//       pngoutD = fopen(OutputFileName, "wb");
//       snprintf(OutputFileName, 255, "%s_extra_%lu.png", FileName[mat], SeqID);
//       pngoutX = fopen(OutputFileName, "wb");
//     } 
//     else {
//       snprintf(OutputFileName, 255, "%s_match_%lu_scale.png", FileName[mat], SeqID);
//       pngoutM = fopen(OutputFileName, "wb");
//       snprintf(OutputFileName, 255, "%s_insertion_%lu_scale.png", FileName[mat], SeqID);
//       pngoutI = fopen(OutputFileName, "wb");
//       snprintf(OutputFileName, 255, "%s_deletion_%lu_scale.png", FileName[mat], SeqID);
//       pngoutD = fopen(OutputFileName, "wb");
//       snprintf(OutputFileName, 255, "%s_extra_%lu_scale.png", FileName[mat], SeqID);
//       pngoutX = fopen(OutputFileName, "wb");
//     }
//     if ( pngoutM == NULL || pngoutI == NULL || pngoutD == NULL || pngoutX == NULL ) {
//       if (pngoutM != NULL) fclose(pngoutM);
//       if (pngoutI != NULL) fclose(pngoutI);
//       if (pngoutD != NULL) fclose(pngoutD);
//       if (pngoutX != NULL) fclose(pngoutX);
//       return;
//     }
// 
//     /* Allocate the images */
//     if (ScaleOnImage) {
//       im[MATCH]     = gdImageCreateTrueColor(ScaleWidth + ScaledImageWidth, ScaledImageHeight);
//       im[INSERTION] = gdImageCreateTrueColor(ScaleWidth + ScaledImageWidth, ScaledImageHeight);
//       im[DELETION]  = gdImageCreateTrueColor(ScaleWidth + ScaledImageWidth, ScaledImageHeight);
//       im[EXTRA]     = gdImageCreateTrueColor(ScaleWidth + ScaledImageWidth, ScaledImageHeight);
//     }
//     else {
//       im[MATCH]     = gdImageCreateTrueColor(ScaleWidth, PALETTE_SIZE/2);
//       im[INSERTION] = gdImageCreateTrueColor(ScaleWidth, PALETTE_SIZE/2);
//       im[DELETION]  = gdImageCreateTrueColor(ScaleWidth, PALETTE_SIZE/2);
//       im[EXTRA]     = gdImageCreateTrueColor(ScaleWidth, PALETTE_SIZE/2);
//     }
//     /*
//      * Allocate the colors.
//      */
//     CreatePalette(im, colors);
// 
//     /* Add reverse to standard */
//     if (mat < 2) {
//       GetExtremumScores( &Maximum[mat].xmm, &Minimum[mat].xmm, Matrices[mat], SeqLength, prfLength);
//     }
//     else {
//       AddScores(matrix, rmatrix, SeqLength, prfLength, &Maximum[mat].xmm, &Minimum[mat].xmm);
//     }
//     
//     /* Draw the scale */
//     int lScaleWidth;
// 
//     if (ScaleOnImage) {
//       int ScaleSize = PALETTE_SIZE;
//       int colorstep = 1;
//       {
// 	while (ScaleSize > ScaledImageHeight) { ScaleSize >>= 1; colorstep <<=1;}
// 	fprintf(stderr, "Scale size reduced to %i, color step of %i\n", ScaleSize, colorstep);
//       }
//       gdImageFilledRectangle(im[MATCH],     0, 0, ScaleWidth, ScaledImageHeight, colors[MATCH][0]);
//       gdImageFilledRectangle(im[INSERTION], 0, 0, ScaleWidth, ScaledImageHeight, colors[INSERTION][0]);
//       gdImageFilledRectangle(im[DELETION],  0, 0, ScaleWidth, ScaledImageHeight, colors[DELETION][0]);
//       gdImageFilledRectangle(im[EXTRA],     0, 0, ScaleWidth, ScaledImageHeight, colors[EXTRA][0]);
//       const int Bottom = (ScaledImageHeight + ScaleSize ) / 2;
//       int lineY  = Bottom;
//       for (int i=0; i<PALETTE_SIZE; i+=colorstep) {
// 	  gdImageLine(im[MATCH],     20, lineY, 80, lineY, colors[MATCH][i]);
// 	  gdImageLine(im[INSERTION], 20, lineY, 80, lineY, colors[INSERTION][i]);
// 	  gdImageLine(im[DELETION],  20, lineY, 80, lineY, colors[DELETION][i]);
// 	  gdImageLine(im[EXTRA],     20, lineY, 80, lineY, colors[EXTRA][i]);
// 	  --lineY;
//       }
//       const gdFontPtr FontPtr = gdFontGetSmall();
//       sprintf(OutputFileName, "%i\0", Minimum[mat].Element[MATCH]);
//       gdImageString(im[MATCH], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[MATCH][255]);
//       sprintf(OutputFileName, "%i\0", prf->CutOff);
//       gdImageString(im[MATCH], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/colorstep - FontPtr->h/2, (unsigned char*)OutputFileName, colors[MATCH][0]);
//       sprintf(OutputFileName, "%i\0", Maximum[mat].Element[MATCH]);
//       gdImageString(im[MATCH], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), (int) SeqLength - Bottom + 2, (unsigned char*)OutputFileName, colors[MATCH][0]);
// 
//       sprintf(OutputFileName, "%i\0", Minimum[mat].Element[INSERTION]);
//       gdImageString(im[INSERTION], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[INSERTION][255]);
//       sprintf(OutputFileName, "%i\0", prf->CutOff);
//       gdImageString(im[INSERTION], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/colorstep - FontPtr->h/2, (unsigned char*)OutputFileName, colors[INSERTION][0]);
//       sprintf(OutputFileName, "%i\0", Maximum[mat].Element[INSERTION]);
//       gdImageString(im[INSERTION], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), (int) SeqLength - Bottom + 2, (unsigned char*)OutputFileName, colors[INSERTION][0]);
// 
//       sprintf(OutputFileName, "%i\0", Minimum[mat].Element[DELETION]);
//       gdImageString(im[DELETION], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[DELETION][255]);
//       sprintf(OutputFileName, "%i\0", prf->CutOff);
//       gdImageString(im[DELETION], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/colorstep - FontPtr->h/2, (unsigned char*)OutputFileName, colors[DELETION][0]);
//       sprintf(OutputFileName, "%i\0", Maximum[mat].Element[DELETION]);
//       gdImageString(im[DELETION], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), (int) SeqLength - Bottom + 2, (unsigned char*)OutputFileName, colors[DELETION][0]);
// 
//       sprintf(OutputFileName, "%i\0", Minimum[mat].Element[EXTRA]);
//       gdImageString(im[EXTRA], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[EXTRA][255]);
//       sprintf(OutputFileName, "%i\0", prf->CutOff);
//       gdImageString(im[EXTRA], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/colorstep - FontPtr->h/2, (unsigned char*)OutputFileName, colors[EXTRA][0]);
//       sprintf(OutputFileName, "%i\0", Maximum[mat].Element[EXTRA]);
//       gdImageString(im[EXTRA], FontPtr, 20 + 60/2 - (strlen(OutputFileName) * FontPtr->w / 2), (int) SeqLength - Bottom + 2, (unsigned char*)OutputFileName, colors[EXTRA][0]);
// 
//       lScaleWidth = ScaleWidth;
//     } 
//     else {
//       gdImageFilledRectangle(im[MATCH],     0, 0, ScaleWidth, PALETTE_SIZE/2, colors[MATCH][0]);
//       gdImageFilledRectangle(im[INSERTION], 0, 0, ScaleWidth, PALETTE_SIZE/2, colors[INSERTION][0]);
//       gdImageFilledRectangle(im[DELETION],  0, 0, ScaleWidth, PALETTE_SIZE/2, colors[DELETION][0]);
//       gdImageFilledRectangle(im[EXTRA],     0, 0, ScaleWidth, PALETTE_SIZE/2, colors[EXTRA][0]);
//       const int Bottom = PALETTE_SIZE/2 - 1;
//       int lineY  = Bottom;
//       for (int i=0; i<PALETTE_SIZE; i+=2) {
// 	  gdImageLine(im[MATCH], 0, lineY, ScaleWidth, lineY, colors[MATCH][i]);
// 	  gdImageLine(im[INSERTION], 0, lineY, ScaleWidth, lineY, colors[INSERTION][i]);
// 	  gdImageLine(im[DELETION], 0, lineY, ScaleWidth, lineY, colors[DELETION][i]);
// 	  gdImageLine(im[EXTRA], 0, lineY, ScaleWidth, lineY, colors[EXTRA][i]);
// 	  --lineY;
//       }
//       const gdFontPtr FontPtr = gdFontGetSmall();
//       sprintf(OutputFileName, "%i\0", Minimum[mat].Element[MATCH]);
//       gdImageString(im[MATCH], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[MATCH][255]);
//       sprintf(OutputFileName, "%i\0", prf->CutOff);
//       gdImageString(im[MATCH], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/2 - FontPtr->h/2, (unsigned char*)OutputFileName, colors[MATCH][0]);
//       sprintf(OutputFileName, "%i\0", Maximum[mat].Element[MATCH]);
//       gdImageString(im[MATCH], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), 2, (unsigned char*)OutputFileName, colors[MATCH][0]);
// 
//       sprintf(OutputFileName, "%i\0", Minimum[mat].Element[INSERTION]);
//       gdImageString(im[INSERTION], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[INSERTION][255]);
//       sprintf(OutputFileName, "%i\0", prf->CutOff);
//       gdImageString(im[INSERTION], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/2 - FontPtr->h/2, (unsigned char*)OutputFileName, colors[INSERTION][0]);
//       sprintf(OutputFileName, "%i\0", Maximum[mat].Element[INSERTION]);
//       gdImageString(im[INSERTION], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), 2, (unsigned char*)OutputFileName, colors[INSERTION][0]);
// 
//       sprintf(OutputFileName, "%i\0", Minimum[mat].Element[DELETION]);
//       gdImageString(im[DELETION], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[DELETION][255]);
//       sprintf(OutputFileName, "%i\0", prf->CutOff);
//       gdImageString(im[DELETION], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/2 - FontPtr->h/2, (unsigned char*)OutputFileName, colors[DELETION][0]);
//       sprintf(OutputFileName, "%i\0", Maximum[mat].Element[DELETION]);
//       gdImageString(im[DELETION], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), 2, (unsigned char*)OutputFileName, colors[DELETION][0]);
// 
//       sprintf(OutputFileName, "%i\0", Minimum[mat].Element[EXTRA]);
//       gdImageString(im[EXTRA], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 2 - FontPtr->h, (unsigned char*)OutputFileName, colors[EXTRA][255]);
//       sprintf(OutputFileName, "%i\0", prf->CutOff);
//       gdImageString(im[EXTRA], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), Bottom - 256/2 - FontPtr->h/2, (unsigned char*)OutputFileName, colors[EXTRA][0]);
//       sprintf(OutputFileName, "%i\0", Maximum[mat].Element[EXTRA]);
//       gdImageString(im[EXTRA], FontPtr, ScaleWidth/2 - (strlen(OutputFileName) * FontPtr->w / 2), 2, (unsigned char*)OutputFileName, colors[EXTRA][0]);
// 
//       /* Output the scale image to the disk file in PNG format. */
//       gdImagePng(im[MATCH],     pngoutM);
//       gdImagePng(im[INSERTION], pngoutI);
//       gdImagePng(im[DELETION],  pngoutD);
//       gdImagePng(im[EXTRA],     pngoutX);
// 
// 
//       /* Close the files. */
//       fclose(pngoutM);
//       fclose(pngoutI);
//       fclose(pngoutD);
//       fclose(pngoutX);
// 
//       /* Destroy the image in memory. */
//       gdImageDestroy(im[MATCH]);
//       gdImageDestroy(im[INSERTION]);
//       gdImageDestroy(im[DELETION]);
//       gdImageDestroy(im[EXTRA]);
// 
//       /* Compose file names */
//       snprintf(OutputFileName, 255, "%s_match_%lu.png", FileName[mat], SeqID);
//       pngoutM = fopen(OutputFileName, "wb");
//       snprintf(OutputFileName, 255, "%s_insertion_%lu.png", FileName[mat], SeqID);
//       pngoutI = fopen(OutputFileName, "wb");
//       snprintf(OutputFileName, 255, "%s_deletion_%lu.png", FileName[mat], SeqID);
//       pngoutD = fopen(OutputFileName, "wb");
//       snprintf(OutputFileName, 255, "%s_extra_%lu.png", FileName[mat], SeqID);
//       pngoutX = fopen(OutputFileName, "wb");
// 
//       if ( pngoutM == NULL || pngoutI == NULL || pngoutD == NULL || pngoutX == NULL ) {
// 	if (pngoutM != NULL) fclose(pngoutM);
// 	if (pngoutI != NULL) fclose(pngoutI);
// 	if (pngoutD != NULL) fclose(pngoutD);
// 	if (pngoutX != NULL) fclose(pngoutX);
// 	return;
//       }
// 
//       /* Allocate the images */
//       im[MATCH]     = gdImageCreateTrueColor((int) prfLength*PNG_CELL_X_SIZE, (int) SeqLength*PNG_CELL_Y_SIZE);
//       im[INSERTION] = gdImageCreateTrueColor((int) prfLength*PNG_CELL_X_SIZE, (int) SeqLength*PNG_CELL_Y_SIZE);
//       im[DELETION]  = gdImageCreateTrueColor((int) prfLength*PNG_CELL_X_SIZE, (int) SeqLength*PNG_CELL_Y_SIZE);
//       im[EXTRA]     = gdImageCreateTrueColor((int) prfLength*PNG_CELL_X_SIZE, (int) SeqLength*PNG_CELL_Y_SIZE);
// 
//       /* Allocate the colors */
//       CreatePalette(im, colors);
// 
//       /* Correct shift */
//       lScaleWidth = 0;
//     }
//     
//     /* Compute scaling and offset */
//     Maximum[mat].xmmf                 = _mm_cvtepi32_ps(Maximum[mat].xmm);
//     Minimum[mat].xmmf                 = _mm_cvtepi32_ps(Minimum[mat].xmm);
//     
//     float cutoff = (float) prf->CutOff;
// //     if (mat == 2) cutoff *= 2.0f;
//     const __m128 __CutoffValue      = _mm_set1_ps(cutoff);
// 
//     const __m128 __AboveCutoffScale = _mm_div_ps(__CutoffRange, _mm_sub_ps(Maximum[mat].xmmf, __CutoffValue));
//     const __m128 __BelowCutoffScale = _mm_div_ps(__BelowCutoffRange, _mm_sub_ps(__CutoffValue, Minimum[mat].xmmf));
// 
//     /* Adjust Min/Max if required */
//     if (((struct IO_PNG*) Options)->GivenScale) {
// 	Minimum[mat].xmm = _mm_set1_epi32(((struct IO_PNG*) Options)->ScaleMin);
// 	Maximum[mat].xmm = _mm_set1_epi32(((struct IO_PNG*) Options)->ScaleMax);
//     }
//     
//     /* Draw the color boxes */
//     const union lScores * restrict pmatrix = &((Matrices[mat])[(1+ImageOffsetY)*matrix_lda+1+ImageOffsetX]);
//     int y = 0;
//     for (int iseq=0; iseq<ImageHeight; ++iseq) {
//       const int y1 = y+PNG_CELL_Y_SIZE-1;
//       int x = lScaleWidth;
//       for (int iprf=0; iprf<ImageWidth; ++iprf) {
// 	/* Get the value in float format */
// 	__m128 __value = _mm_cvtepi32_ps(pmatrix[iprf].xmm);
// 	/* Set default for below cutoff value */
// 	__m128 __Scale  = __BelowCutoffScale;
// 	__m128 __Offset = _mm_setzero_ps();
// 	__m128 __Subs   = Minimum[mat].xmmf;
// 	/* Is it above cutoff ? update parameters*/
// 	__m128 __Mask = _mm_cmpgt_ps(__value, __CutoffValue);
// 	__Scale       = _mm_blendv_ps(__Scale, __AboveCutoffScale, __Mask);
// 	__Offset      = _mm_blendv_ps(__Offset, __CutoffOffset, __Mask);
// 	__Subs        = _mm_blendv_ps(__Subs, __CutoffValue, __Mask);
// 
// 	__value = _mm_add_ps(__Offset,_mm_mul_ps(_mm_sub_ps(__value, __Subs), __Scale));
// 	union lScores pixelColor = { xmm: _mm_cvtps_epi32(__value)};
// 	pixelColor.xmm = _mm_max_epi32(pixelColor.xmm, _mm_setzero_si128());
// 	pixelColor.xmm = _mm_min_epi32(pixelColor.xmm, _mm_set1_epi32(PALETTE_SIZE-1));
// 	
// 	gdImageFilledRectangle(im[MATCH],     x, y, x+PNG_CELL_X_SIZE-1, y1, colors[MATCH][pixelColor.Element[MATCH]]);
// 	gdImageFilledRectangle(im[INSERTION], x, y, x+PNG_CELL_X_SIZE-1, y1, colors[INSERTION][pixelColor.Element[INSERTION]]);
// 	gdImageFilledRectangle(im[DELETION],  x, y, x+PNG_CELL_X_SIZE-1, y1, colors[DELETION][pixelColor.Element[DELETION]]);
// 	gdImageFilledRectangle(im[EXTRA],     x, y, x+PNG_CELL_X_SIZE-1, y1, colors[EXTRA][pixelColor.Element[EXTRA]]);
// 	x += PNG_CELL_X_SIZE;
//       }
//       y += PNG_CELL_Y_SIZE;
//       pmatrix += matrix_lda;
//     }
// 
//     /* Output the image to the disk file in PNG format. */
//     gdImagePng(im[MATCH],     pngoutM);
//     gdImagePng(im[INSERTION], pngoutI);
//     gdImagePng(im[DELETION],  pngoutD);
//     gdImagePng(im[EXTRA],     pngoutX);
// 
// 
//     /* Close the files. */
//     fclose(pngoutM);
//     fclose(pngoutI);
//     fclose(pngoutD);
//     fclose(pngoutX);
// 
//     /* Destroy the image in memory. */
//     gdImageDestroy(im[MATCH]);
//     gdImageDestroy(im[INSERTION]);
//     gdImageDestroy(im[DELETION]);
//     gdImageDestroy(im[EXTRA]);
//   }
// }
