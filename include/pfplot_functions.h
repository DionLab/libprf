/*******************************************************
                        PFTOOLS
 *******************************************************
  Oct 19, 2012 pfplot_functions.h
 *******************************************************
 (C) 2011 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/

#include <xmmintrin.h>
#include "pfProfile.h"


typedef void (*FuncPtr)( __m128i, __m128i * restrict const, const size_t, const size_t, const size_t);

extern int PlotAllSequence;
