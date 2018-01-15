
/*******************************************************
                        PFTOOLS
 *******************************************************
  Oct 12, 2012 xali1_print_sse41.c
 *******************************************************
 (C) 2012 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <smmintrin.h>
#include <assert.h>
#include "pfCompute.h"
#include "pfOutput.h"
#include "sse41_inline_fcts.h"

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
  * |xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|xx|
  *                               |  |
  *                               |  +-- State
  *                               +----- Score
  */
# define SCORE_SHIFT               2
# define SCORE_MASK                0xFFFFFFFC
# define STATE_SHIFT               0
# define STATE_MASK                0x00000003
# define LEFT_NEGATIVE_SCORE_MASK  0xC0000000
# define CLEAN_STATE               0xFFFFFFF3
// # define CLEAN_PHASE               0xFFFFFFFC
// # define CLEAN_STATE_AND_PHASE     0xFFFFFFF0
#define MATCH_UPPER_BITS     0xA000
#define INSERTION_UPPER_BITS 0xB000
#define DELETION_UPPER_BITS  0xC000

#define CIGAR_SIZE 124

struct data_s {
    uint32_t score;
    uint16_t position;
    uint16_t start;
    uint16_t cigar[CIGAR_SIZE];
} __attribute__((aligned(16)));

typedef struct Cell_s {
    struct data_s Match;
    struct data_s Insertion;
} Cell_t;

#define TO_SCORE(x) (int)(((x)<0) ? ((((unsigned int)(x)) >> SCORE_SHIFT) | LEFT_NEGATIVE_SCORE_MASK ) : (((unsigned int)x) >> SCORE_SHIFT))
#define MAX(a,b) ((a>b) ? a : b)

void horizontal_sse41(const struct Profile * const restrict prf,
                             const unsigned char * const restrict Sequence,
                             struct data_s * const restrict matrix, Cell_t * const restrict WORK,
                             const size_t BSEQ, const size_t LSEQ)
/*
 * WARNING: for SSE version, 
 *          WORK should be aligned on cache size (64b) and 4 times the ((profile size + 1) + (profile_size+1))*sizeof(Cell_t) + 63 to align to cache line
 *          matrix should be of size at least (sequence size+1)*sizeof(Cell_t) and aligned on 16b.
 */
{
  Cell_t * restrict IOP_W, * restrict IOP_R;
  const __m128i __MatchMask     = _mm_set1_epi32(PRIORITY_MATCH);
  const __m128i __InsertionMask = _mm_set1_epi32(PRIORITY_INSERTION);
  const __m128i __DeletionMask  = _mm_set1_epi32(PRIORITY_DELETION);
  const __m128i __ExtraMask     = _mm_set1_epi32(PRIORITY_EXTRA);
  
  register const TransitionScores * const restrict Transitions = prf->Scores.Insertion.Transitions;
  const StoredIntegerFormat * const restrict Match             = prf->Scores.Match.Alphabet;
  const StoredIntegerFormat * const restrict Insertion         = prf->Scores.Insertion.Alphabet;
  const size_t AlignStep                                       = prf->Scores.Match.AlignStep;
  const size_t prfLength = prf->Length;
  
  struct data_s * restrict KOPD = (struct data_s *) (((uintptr_t) (WORK + 2*(prfLength+1)) + 63UL) & ~(63UL));
  IOP_W = WORK;
  /*
   * Initialize Insertion and Match Entrance Line using FirstSequenceProtein
   */  
  {
        register const StoredIntegerFormat * restrict lMatch = (const StoredIntegerFormat *) &Match[_D];
        register const ScoreTuple * restrict FirstSequenceProtein = prf->Scores.Insertion.FirstSequenceProtein;
        /*
        * PROFILE COLUMN 0 entrance
        */
        IOP_W[0].Match.score        = (int) FirstSequenceProtein[0].To[MATCH];
        IOP_W[0].Match.start        = 0;
        IOP_W[0].Match.position     = 0;
        IOP_W[0].Match.cigar[0]     = MATCH_UPPER_BITS;
        IOP_W[0].Insertion.score    = (int) FirstSequenceProtein[0].To[INSERTION];
        IOP_W[0].Insertion.start    = 0;
        IOP_W[0].Insertion.position = 0;
        IOP_W[0].Insertion.cigar[0] = INSERTION_UPPER_BITS;
        KOPD->score                 = (int) FirstSequenceProtein[0].To[DELETION];
        KOPD->start                 = 0;
        KOPD->position              = 0;
        KOPD->cigar[0]              = DELETION_UPPER_BITS;
        
        FirstSequenceProtein++;
        register const TransitionScores (* restrict pTransitions) = &Transitions[1];
            
        /*
        * LOOP THROUGH THE REST OF THE PROFILE
        */
        for (size_t iprf=1; iprf<=prfLength; ++iprf ) {
            register const int KD = KOPD->score + (int) *lMatch;
            lMatch += AlignStep;

            // Transform KD into a vector
            __m128i __KD = _mm_set1_epi32(KD);
            // Load Transitions and  Convert signed WORD into signed DWORD
            __m128i __TransitionsD = LoadStoredIntegerVector(&(pTransitions->From[DELETION]));

            // Add KD to Transitions
            __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

            // Move to next profile transitions
            pTransitions++;

            // Load FirstSequenceProtein
            __m128i __FirstSequenceProtein = LoadStoredIntegerVector(&(FirstSequenceProtein[0]));

            // Move to next profile First Sequence
            FirstSequenceProtein++;
            
            // Paste index in lowest 2 bits
            __TransitionsD = _mm_slli_epi32(__TransitionsD, SCORE_SHIFT);
            __FirstSequenceProtein = _mm_slli_epi32(__FirstSequenceProtein, SCORE_SHIFT);
            __TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);
            __FirstSequenceProtein = _mm_or_si128(__FirstSequenceProtein, __ExtraMask);
            
            // Get maximum ( this is SSE 4.1 )
            __m128i __max = _mm_max_epi32(__TransitionsD, __FirstSequenceProtein);

            // Match 
            const int match = _mm_extract_epi32(__max, MATCH);
            IOP_W[iprf].Match.score = (match >> SCORE_SHIFT);
            switch(match & STATE_MASK) {
                case (PRIORITY_EXTRA): 
                    IOP_W[iprf].Match.start = (int) iprf;
                    IOP_W[iprf].Match.position = 0;
                    IOP_W[iprf].Match.cigar[0] = MATCH_UPPER_BITS;
                    break;
                case (PRIORITY_DELETION): 
                    IOP_W[iprf].Match.start = KOPD->start;
                    IOP_W[iprf].Match.position = 1;
                    IOP_W[iprf].Match.cigar[1] = MATCH_UPPER_BITS;
                    break;
            }

            // Insertion
            const int insertion = _mm_extract_epi32(__max, INSERTION);
            IOP_W[iprf].Insertion.score = (insertion >> SCORE_SHIFT);
            switch(insertion & STATE_MASK) {
                case (PRIORITY_EXTRA): 
                    IOP_W[iprf].Insertion.start = (int) iprf;
                    IOP_W[iprf].Insertion.position = 0;
                    IOP_W[iprf].Insertion.cigar[0] = INSERTION_UPPER_BITS;
                    break;
                case (PRIORITY_DELETION): 
                    IOP_W[iprf].Insertion.start = KOPD->start;
                    IOP_W[iprf].Insertion.position = 1;
                    IOP_W[iprf].Insertion.cigar[1] = INSERTION_UPPER_BITS;
                    break;
            }
        
            // Deletion
            const int deletion = _mm_extract_epi32(__max, DELETION);
            KOPD->score = (deletion >> SCORE_SHIFT);
            switch(deletion & STATE_MASK) {
                case (PRIORITY_EXTRA): 
                    KOPD->start = (int) iprf;
                    KOPD->position = 0;
                    KOPD->cigar[0] = DELETION_UPPER_BITS;
                    break;
                case (PRIORITY_DELETION): 
                    KOPD->cigar[0]++;
                    break;
            }
        }
  }

  // Swap and assign Read and write pointers
  IOP_R = IOP_W;
  IOP_W = &WORK[prfLength + 1];

  /*
   * LOOP THROUGH THE SEQUENCE STRING
   */
  for ( size_t iseq=BSEQ; iseq < LSEQ-1; ++iseq) {
    register const size_t j1 = (size_t) Sequence[iseq];
    register const StoredIntegerFormat * restrict lInsertion = Insertion;
    /*
     * PROFILE COLUMN 0 entrance
     */
    {
        register const int KI = IOP_R[0].Insertion.score + (int) lInsertion[j1];

        // Transform KI into a vector
        __m128i __KI = _mm_set1_epi32(KI);
        // Load Transitions and Convert signed WORD into signed DWORD
        __m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[0].From[INSERTION]));
        // Add KI to Transition
        __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

        // Load Transitions and Convert signed WORD into signed DWORD
        __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[0].From[EXTRA]));

        // Insert lScore into __TransitionsX not necessary as profile loading store NLOW there automatically
        //__TransitionsX = _mm_insert_epi32(__TransitionsX, NLOW, DUMMY);
        
        // Paste index in lowest 2 bits
        __TransitionsI = _mm_slli_epi32(__TransitionsI, SCORE_SHIFT);
        __TransitionsX = _mm_slli_epi32(__TransitionsX, SCORE_SHIFT);
        __TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);
        __TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);

        // Get maximum ( this is SSE 4.1 )
        __m128i __max = _mm_max_epi32(__TransitionsI, __TransitionsX);
        
        // Match
        const int match = _mm_extract_epi32(__max, MATCH);
        IOP_W[0].Match.score = (match >> SCORE_SHIFT);
        switch(match & STATE_MASK) {
            case (PRIORITY_EXTRA): 
                IOP_W[0].Match.start    = 0;
                IOP_W[0].Match.position = 0;
                IOP_W[0].Match.cigar[0] = MATCH_UPPER_BITS;
                break;
            case (PRIORITY_INSERTION): 
                IOP_W[0].Match.start    = IOP_R[0].Insertion.start;
                IOP_W[0].Match.position = 1;
                IOP_W[0].Match.cigar[0] = IOP_R[0].Insertion.cigar[0];
                IOP_W[0].Match.cigar[1] = MATCH_UPPER_BITS;
                break;
        }
      
        // Insertion
        const int insertion = _mm_extract_epi32(__max, INSERTION);
        IOP_W[0].Match.score = (insertion >> SCORE_SHIFT);
        switch(insertion & STATE_MASK) {
            case (PRIORITY_EXTRA): 
                IOP_W[0].Insertion.start    = 0;
                IOP_W[0].Insertion.position = 0;
                IOP_W[0].Insertion.cigar[0] = INSERTION_UPPER_BITS;
                break;
            case (PRIORITY_INSERTION):
                IOP_W[0].Insertion.start    = IOP_R[0].Insertion.start;
                IOP_W[0].Insertion.position = IOP_R[0].Insertion.position;
                IOP_W[0].Insertion.cigar[0] = IOP_R[0].Insertion.cigar[0] + 1;
                break;
        }

        // Deletion
        KOPD++;
        const int deletion = _mm_extract_epi32(__max, DELETION);
        KOPD->score = (deletion >> SCORE_SHIFT);
        switch(deletion & STATE_MASK) {
            case (PRIORITY_EXTRA): 
                KOPD->start = 0;
                KOPD->position = 0;
                KOPD->cigar[0] = DELETION_UPPER_BITS;
                break;
            case (PRIORITY_INSERTION):
                KOPD->start = 0;
                KOPD->position = 1;
                KOPD->cigar[0] = IOP_R[0].Insertion.cigar[0];
                KOPD->cigar[1] = DELETION_UPPER_BITS;
                break;
        }
    }

    lInsertion += AlignStep;
    register const StoredIntegerFormat * restrict lMatch = Match;

    /*
     * LOOP THROUGH THE REST OF THE PROFILE
     */
    for (size_t iprf=1; iprf<=prfLength; ++iprf ) {
        const int KM = IOP_R[iprf-1].Match.score   + (int) lMatch[j1];
        const int KI = IOP_R[iprf].Insertion.score + (int) lInsertion[j1];
        const int KD = KOPD->score                 + (int) lMatch[_D];

        lMatch     += AlignStep;
        lInsertion += AlignStep;

        // Transform KM into a vector
        __m128i __KM = _mm_set1_epi32(KM);
        // Load Transitions and Convert signed WORD into signed DWORD
        __m128i __TransitionsM = LoadStoredIntegerVector(&(Transitions[iprf].From[MATCH]));
        // Add KM to Transition
        __TransitionsM = _mm_add_epi32(__TransitionsM, __KM);


        // Transform KI into a vector
        __m128i __KI = _mm_set1_epi32(KI);
        // Load Transitions and Convert signed WORD into signed DWORD
        __m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[iprf].From[INSERTION]));
        // Add KI to Transition
        __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

        // Paste index in lowest 2 bits
        __TransitionsM = _mm_slli_epi32(__TransitionsM, SCORE_SHIFT);
        __TransitionsI = _mm_slli_epi32(__TransitionsI, SCORE_SHIFT);
        __TransitionsM = _mm_or_si128(__TransitionsM, __MatchMask);
        __TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);

        // Get maximum ( this is SSE 4.1 )
        __m128i __max1 = _mm_max_epi32(__TransitionsM, __TransitionsI);      

        // Load Transitions and Convert signed WORD into signed DWORD
        __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));
        // Insert lscore into TransitionX not necessary as profile loading should have already done it
        // __TransitionsX = _mm_insert_epi32(__TransitionsX, lScore, DUMMY);

        // Transform KD into a vector
        __m128i __KD = _mm_set1_epi32(KD);
        // Load Transitions and Convert signed WORD into signed DWORD
        __m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
        // Add KD to Transition
        __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

        // Paste index in lowest 2 bits
        __TransitionsX = _mm_slli_epi32(__TransitionsX, SCORE_SHIFT);
        __TransitionsD = _mm_slli_epi32(__TransitionsD, SCORE_SHIFT);
        __TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);
        __TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);

        // Get maximum ( this is SSE 4.1 )
        __m128i __max2 = _mm_max_epi32(__TransitionsD, __TransitionsX);
        __max1 = _mm_max_epi32(__max1, __max2);
        
        // Match
        const int match = _mm_extract_epi32(__max1, MATCH);
        IOP_W[iprf].Match.score = (match >> SCORE_SHIFT);
        switch(match & STATE_MASK) {
            case (PRIORITY_MATCH):
                IOP_W[iprf].Match.start    = IOP_R[iprf-1].Match.start;
                IOP_W[iprf].Match.position = IOP_R[iprf-1].Match.position;
                { uint16_t i=0; do { IOP_W[iprf].Match.cigar[i] = IOP_R[iprf-1].Match.cigar[i]; } while (++i < IOP_R[iprf-1].Match.position); }
                IOP_W[iprf].Match.cigar[IOP_W[iprf].Match.position]++;
                break;
            case (PRIORITY_INSERTION): 
                IOP_W[iprf].Match.start    = IOP_R[iprf].Insertion.start;
                IOP_W[iprf].Match.position = IOP_R[iprf].Insertion.position + 1;
                { uint16_t i=0; do { IOP_W[iprf].Match.cigar[i] = IOP_R[iprf].Insertion.cigar[i]; } while (++i < IOP_R[iprf].Insertion.position); }
                IOP_W[iprf].Match.cigar[IOP_R[iprf].Insertion.position + 1] = MATCH_UPPER_BITS;
                break;
            case (PRIORITY_DELETION):
                IOP_W[iprf].Match.start    = KOPD->start;
                IOP_W[iprf].Match.position = KOPD->position + 1;
                { uint16_t i=0; do { IOP_W[iprf].Match.cigar[i] = KOPD->cigar[i]; } while (++i < KOPD->position); }
                IOP_W[iprf].Match.cigar[KOPD->position + 1] = MATCH_UPPER_BITS;
            case (PRIORITY_EXTRA): 
                IOP_W[iprf].Match.start    = (int) iprf;
                IOP_W[iprf].Match.position = 0;
                IOP_W[iprf].Match.cigar[0] = MATCH_UPPER_BITS;
                break;
        }
      
        // Insertion
        const int insertion = _mm_extract_epi32(__max1, INSERTION);
        IOP_W[iprf].Insertion.score = (insertion >> SCORE_SHIFT);
        switch(insertion & STATE_MASK) {
            case (PRIORITY_MATCH):
                IOP_W[iprf].Insertion.start    = IOP_R[iprf-1].Match.start;
                IOP_W[iprf].Insertion.position = IOP_R[iprf-1].Match.position + 1;
                { uint16_t i=0; do { IOP_W[iprf].Insertion.cigar[i] = IOP_R[iprf-1].Match.cigar[i]; } while (++i < IOP_R[iprf-1].Match.position); }
                IOP_W[iprf].Insertion.cigar[IOP_R[iprf-1].Match.position+1] = INSERTION_UPPER_BITS;
                break;
            case (PRIORITY_INSERTION):
                IOP_W[iprf].Insertion.start    = IOP_R[iprf].Insertion.start;
                IOP_W[iprf].Insertion.position = IOP_R[iprf].Insertion.position;
                { uint16_t i=0; do { IOP_W[iprf].Insertion.cigar[i] = IOP_R[iprf].Insertion.cigar[i]; } while (++i < IOP_R[iprf].Insertion.position); }
                IOP_W[iprf].Insertion.cigar[IOP_R[iprf].Insertion.position]++;
                break;
            case (PRIORITY_DELETION):
                IOP_W[iprf].Insertion.start    = KOPD->start;
                IOP_W[iprf].Insertion.position = KOPD->position + 1;
                { uint16_t i=0; do { IOP_W[iprf].Insertion.cigar[i] = KOPD->cigar[i]; } while (++i < KOPD->position); }
                IOP_W[iprf].Insertion.cigar[KOPD->position + 1] = MATCH_UPPER_BITS;
                break;
            case (PRIORITY_EXTRA): 
                IOP_W[iprf].Insertion.start    = (int) iprf;
                IOP_W[iprf].Insertion.position = 0;
                IOP_W[iprf].Insertion.cigar[0] = INSERTION_UPPER_BITS;
                break;
        }

        // Deletion
        const int deletion = _mm_extract_epi32(__max1, DELETION);
        KOPD->score = (deletion >> SCORE_SHIFT);
        switch(deletion & STATE_MASK) {
            case (PRIORITY_MATCH):
                KOPD->start = IOP_R[iprf-1].Match.start;
                KOPD->position = IOP_R[iprf-1].Match.position + 1;
                { uint16_t i=0; do { KOPD->cigar[i] = IOP_R[iprf-1].Match.cigar[i]; } while (++i < IOP_R[iprf-1].Match.position); }
                KOPD->cigar[IOP_R[iprf-1].Match.position + 1] = DELETION_UPPER_BITS;
                break;
            case (PRIORITY_INSERTION):
                KOPD->start = IOP_R[iprf].Insertion.start;
                KOPD->position = IOP_R[iprf].Insertion.position + 1;
                { uint16_t i=0; do { KOPD->cigar[i] = IOP_R[iprf].Insertion.cigar[i]; } while (++i < IOP_R[iprf].Insertion.position); }
                KOPD->cigar[IOP_R[iprf].Insertion.position + 1] = DELETION_UPPER_BITS;
                break;
            case (PRIORITY_DELETION):
                KOPD->cigar[KOPD->position]++;
                break;
            case (PRIORITY_EXTRA): 
                KOPD->start = 0;
                KOPD->position = 0;
                KOPD->cigar[0] = DELETION_UPPER_BITS;
                break;
        }
    }

    // Swap Read and Write pointers
    const register Cell_t * const ptr = IOP_W;
    IOP_W = (Cell_t*) IOP_R;
    IOP_R = (Cell_t*) ptr;
  }

  /*
   * Last position on the Sequence using LastSequenceProtein
   */
  {
    register const StoredIntegerFormat * restrict lInsertion = Insertion;
    const int j1 = (int) Sequence[LSEQ-1];
    /*
     * PROFILE COLUMN 0 entrance
     */
    int KI = IOP_R[0].Insertion.score + (int) lInsertion[j1];
    KOPD++;
    if ((KI + (int) Transitions[0].Element[_ID]) >= (int) Transitions[0].Element[_XD]) {
        KOPD->score    = (KI + (int) Transitions[0].Element[_ID]);
        KOPD->start    = IOP_R[0].Insertion.start;
        KOPD->position = IOP_R[0].Insertion.position + 1;
        { uint16_t i=0; do { KOPD->cigar[i] = IOP_R[0].Insertion.cigar[i]; } while (++i < IOP_R[0].Insertion.position); }
        KOPD->cigar[IOP_R[0].Insertion.position + 1] = DELETION_UPPER_BITS;
    }
    else {
        KOPD->score    = (int) Transitions[0].Element[_XD];
        KOPD->start    = 0;
        KOPD->position = 0;
        KOPD->cigar[0] = DELETION_UPPER_BITS;
    }
    register const ScoreTuple * const restrict LastSequenceProtein = prf->Scores.Insertion.LastSequenceProtein;
    register const StoredIntegerFormat * restrict lMatch = Match;
    lInsertion += AlignStep;

    /*
     * LOOP THROUGH THE REST OF THE PROFILE
     */
    for (size_t iprf=1; iprf<=prfLength; ++iprf) {
        const int KM = IOP_R[iprf-1].Match.score   + lMatch[j1];
        KI           = IOP_R[iprf].Insertion.score + lInsertion[j1];
        const int KD = KOPD->score                 + lMatch[_D];

        lMatch     += AlignStep;
        lInsertion += AlignStep;

        // Transform KM into a vector
        __m128i __KM = _mm_set1_epi32(KM);
        // Load Transitions and Convert signed WORD into signed DWORD
        __m128i __TransitionsM = LoadStoredIntegerVector(&(Transitions[iprf].From[MATCH]));
        // Insert LastProteinSequence
        __TransitionsM = _mm_insert_epi32(__TransitionsM, (int) LastSequenceProtein[iprf].From[MATCH], DUMMY);
        // Add KM to Transition
        __TransitionsM = _mm_add_epi32(__TransitionsM, __KM);


        // Transform KI into a vector
        __m128i __KI = _mm_set1_epi32(KI);
        // Load Transitions and Convert signed WORD into signed DWORD
        __m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[iprf].From[INSERTION]));
        // Insert LastProteinSequence
        __TransitionsI = _mm_insert_epi32(__TransitionsI, (int) LastSequenceProtein[iprf].From[INSERTION], DUMMY);
        // Add KI to Transition
        __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

        // Paste index in lowest 2 bits
        __TransitionsM = _mm_slli_epi32(__TransitionsM, SCORE_SHIFT);
        __TransitionsI = _mm_slli_epi32(__TransitionsI, SCORE_SHIFT);
        __TransitionsM = _mm_or_si128(__TransitionsM, __MatchMask);
        __TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);

        // Get maximum ( this is SSE 4.1 )
        __m128i __max1 = _mm_max_epi32(__TransitionsM, __TransitionsI);

        // Load Transitions and Convert signed WORD into signed DWORD
        __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));
        // Insert lscore into TransitionX not necessary as profile loading should have already done it
        // __TransitionsX = _mm_insert_epi32(__TransitionsX, lScore, DUMMY);

        // Transform KD into a vector
        __m128i __KD = _mm_set1_epi32(KD);
        // Load Transitions and Convert signed WORD into signed DWORD
        __m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
        // Insert LastProteinSequence
        __TransitionsD = _mm_insert_epi32(__TransitionsD, (int) LastSequenceProtein[iprf].From[DELETION], DUMMY);
        // Add KD to Transition
        __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

        // Paste index in lowest 2 bits
        __TransitionsX = _mm_slli_epi32(__TransitionsX, SCORE_SHIFT);
        __TransitionsD = _mm_slli_epi32(__TransitionsD, SCORE_SHIFT);
        __TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);
        __TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);

        // Get maximum ( this is SSE 4.1 )
        __m128i __max2 = _mm_max_epi32(__TransitionsD, __TransitionsX);
        __max1 = _mm_max_epi32(__max1, __max2);

        // Best to end
        const int best = _mm_extract_epi32(__max1, EXTRA);
        matrix[iprf-1].score = (best >> SCORE_SHIFT);
        switch(best & STATE_MASK) {
            case (PRIORITY_MATCH):
                matrix[iprf-1].start    = IOP_R[iprf-1].Match.start;
                matrix[iprf-1].position = IOP_R[iprf-1].Match.position;
                { uint16_t i=0; do { matrix[iprf-1].cigar[i] = IOP_R[iprf-1].Match.cigar[i]; } while (++i < IOP_R[iprf-1].Match.position); }
                matrix[iprf-1].cigar[IOP_W[iprf].Match.position] = IOP_R[iprf-1].Match.cigar[IOP_R[iprf-1].Match.position] + 1;
                break;
            case (PRIORITY_INSERTION): 
                matrix[iprf-1].start    = IOP_R[iprf].Insertion.start;
                matrix[iprf-1].position = IOP_R[iprf].Insertion.position;
                { uint16_t i=0; do { matrix[iprf].cigar[i] = IOP_R[iprf].Insertion.cigar[i]; } while (++i < IOP_R[iprf].Insertion.position); }
                matrix[iprf-1].cigar[IOP_W[iprf].Match.position] = IOP_R[iprf].Insertion.cigar[IOP_R[iprf].Insertion.position] + 1;
                break;
            case (PRIORITY_DELETION):
                matrix[iprf-1].start    = KOPD->start;
                matrix[iprf-1].position = KOPD->position;
                { uint16_t i=0; do { IOP_W[iprf].Match.cigar[i] = KOPD->cigar[i]; } while (++i < KOPD->position); }
                matrix[iprf-1].cigar[KOPD->position] = KOPD->cigar[KOPD->position] + 1;
        }

        // Deletion
        const int deletion = _mm_extract_epi32(__max1, DELETION);
        KOPD->score = (deletion >> SCORE_SHIFT);
        switch(deletion & STATE_MASK) {
            case (PRIORITY_MATCH):
                KOPD->start = IOP_R[iprf-1].Match.start;
                KOPD->position = IOP_R[iprf-1].Match.position + 1;
                { uint16_t i=0; do { KOPD->cigar[i] = IOP_R[iprf-1].Match.cigar[i]; } while (++i < IOP_R[iprf-1].Match.position); }
                KOPD->cigar[IOP_R[iprf-1].Match.position + 1] = DELETION_UPPER_BITS;
                break;
            case (PRIORITY_INSERTION):
                KOPD->start = IOP_R[iprf].Insertion.start;
                KOPD->position = IOP_R[iprf].Insertion.position + 1;
                { uint16_t i=0; do { KOPD->cigar[i] = IOP_R[iprf].Insertion.cigar[i]; } while (++i < IOP_R[iprf].Insertion.position); }
                KOPD->cigar[IOP_R[iprf].Insertion.position + 1] = DELETION_UPPER_BITS;
                break;
            case (PRIORITY_DELETION):
                KOPD->cigar[KOPD->position]++;
                break;
            case (PRIORITY_EXTRA): 
                KOPD->start = 0;
                KOPD->position = 0;
                KOPD->cigar[0] = DELETION_UPPER_BITS;
                break;
        }
    }
  }
}


#undef MAX
