/*******************************************************
                        PFTOOLS
 *******************************************************
  Oct 19, 2012 pfplot_functions.h
 *******************************************************
 (C) 2011 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/

#ifndef _MATRIX_FUNCTION_H
#define _MATRIX_FUNCTION_H
#include <xmmintrin.h>

#define TAG

/* Priority */
enum StatePriority {
  PRIORITY_MATCH     = 1,
  PRIORITY_INSERTION = 2,
  PRIORITY_DELETION  = 0,
  PRIORITY_EXTRA     = 3
};


/*
  * INTEGER FORMAT  -------------------------------------
  *
  * |xxxxxxxxxxxxxxxxxxxxxxxxxxxxx|xx|x|
  *                              |  | |
	*                              |  | + Cycle
  *                              |  +-- State
  *                              +----- Score
  */
# define SCORE_SHIFT               3
# define SCORE_MASK                0xFFFFFFF8
# define STATE_SHIFT               1
# define STATE_MASK                0x00000006
# define CYCLE_SHIFT               0
# define CYCLE_MASK                0x00000001
# define LEFT_NEGATIVE_SCORE_MASK  0xE0000000
# define CLEAN_STATE               0xFFFFFFF8
// # define CLEAN_PHASE               0xFFFFFFFC
// # define CLEAN_STATE_AND_PHASE     0xFFFFFFF0


# ifdef TAG
#  define TO_SCORE(x) (int) ((x < 0) ? (((unsigned int) x)>>SCORE_SHIFT) | LEFT_NEGATIVE_SCORE_MASK : ((unsigned int) x)>>SCORE_SHIFT)
#  define GET_STATE(x) ((x & STATE_MASK) >> STATE_SHIFT)
# else
#  define TO_SCORE(x) x
#  define GET_STATE(x) x
# endif

typedef void (*FuncPtr)( __m128i, __m128i * restrict const, const size_t, const size_t, const size_t);

extern int PlotAllSequence;

/*
 * Output functions
 */

/*
 * Structures
 */


/*
 * Global variables
 */

/*
 * Available functions
 */


#endif /*_MATRIX_FUNCTION_H*/
