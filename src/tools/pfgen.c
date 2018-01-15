#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <errno.h>
#define _GNU_SOURCE
#include <string.h>
#define __USE_INLINE_FUNCTIONS__
#include "pfProfile.h"
#include "pfIOInline.h"

#include <getopt.h>
#define HEADER "%----------------------------------------------------------------------------%\n"\
	       "|                              PFDUMP v" PF_VERSION "                                   |\n"\
	       "%----------------------------------------------------------------------------%\n"\
	       "Built on " __DATE__ " at " __TIME__ ".\n"
 
				 
/* default value for header is 4 digits*/
// #define INT_FORMAT "%4i"
// #define NLOW_FORMAT "NLOW"
// #define MLOW_FORMAT "MLOW"
// #define SPACE "  "
//#define LINE "-----"

#define INT_FORMAT "%6i"
#define NLOW_FORMAT "  NLOW"
#define MLOW_FORMAT "  MLOW"
#define SPACE "    "
#define LINE "-------"

#ifdef PRF_USE_AFFINITY
#include <pthread.h>
// cpu_set_t * Thread_masks[2] = {0,0};						/* Define variables to hold thread affinity mask */
// unsigned int Thread_count[2] = {0,0};
pthread_attr_t * restrict threads_attr = NULL;
#endif

static const char opt_to_test[]     = "CcrAapmibtTh";
static unsigned int out_match       = 0;
static unsigned int out_insertion   = 0;
static unsigned int out_boundaries  = 0;
static unsigned int out_transitions = 0;
static unsigned int isReversed      = 0;
static unsigned int out_alphabet    = 0;
static unsigned int out_clean       = 0;
static unsigned int out_sequence    = 0; 

static unsigned int UseColor = 0;
static unsigned int out_profile = 0;

static const struct option long_options[] =
{
	/*
	 * These options set a flag. 
	 */
	
	/* 
	 * These options don't set a flag. We distinguish them by their indices. 
	 */
	{"help",          	no_argument,   	0,	'h'},	
	{"all",     				no_argument,		0,	'A'},
	{"alphabet",				no_argument,		0,	'a'},
	{"color",						no_argument,		0,	'c'},
	{"reverse",					no_argument,		0,	'r'},
	{"profile",					no_argument,		0,	'p'},
	{"match-table",			no_argument,		0,	'm'},
	{"insertion-table",	no_argument,		0,	'i'},
	{"boundary-table",	no_argument,		0,	'b'},
	{"sequence",				no_argument, 		0,	's'}, 
	{"transition-table",no_argument,		0,	't'},
	{"tables",					no_argument,		0,	'T'},
	{"clean",						no_argument,		0,	'C'},
	{0, 0, 0, 0}
};

static void __attribute__((noreturn)) Usage(FILE * stream)
{
  fputs(
	" pfdump [options] profile database\n\n"
	" Options:\n"
	"   --color             [-c] : output profile reading with colors\n"
	"   --reverse           [-r] : reverse profile\n"
	"   --all               [-A] : output everything\n"
	"   --tables            [-T] : output score tables\n" 
	"   --alphabet          [-a] : output alphabet\n"
	"   --profile           [-p] : output profile reading\n"
	"   --sequence          [-s] : output profile query sequence\n"
	"   --match-table       [-m] : output match table\n" 
	"   --insertion-table   [-i] : output insertion table\n"
	"   --boundary-table    [-b] : output boundary table\n"
	"   --transition-table  [-t] : output transition table\n"
	"   --clean             [-C] : output profile without comments\n"
	"   --help              [-h] : output command help\n",
	stream);
  exit(0);
}

int main(int argc, char *argv[])
{
  struct Profile prf, * rprf, * pprf;
  char * ProfileFile;
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // OPTIONS
  ////////////////////////////////////////////////////////////////////////////////////////////////
 
  while (1) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    const int c = getopt_long (argc, argv, opt_to_test, long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) break;
    switch (c) {
      case 'c': UseColor = 1; break;
      case 'r': isReversed = 1; break;
      case 'p': out_profile = 1; break;
      case 'm': out_match = 1; break;
      case 'i': out_insertion = 1; break;
      case 'b': out_boundaries = 1; break;
      case 't': out_transitions = 1; break;
      case 'a': out_alphabet = 1; break;
			case 's': out_sequence = 1; break;
      case 'A':
				out_profile     = 1;
				out_match       = 1;
				out_insertion   = 1;
				out_boundaries  = 1;
				out_transitions = 1;
				out_alphabet    = 1;
				out_sequence    = 1;
				break;
						case 'T':
				out_match       = 1;
				out_insertion   = 1;
				out_boundaries  = 1;
				out_transitions = 1;
				break;
						case 'C':
					out_clean = 1;
				break;
      case 'h':
      default:
	Usage(stdout);
    }
  }

  if (optind == argc) {
    fputs("Error in given options\n", stderr);
    Usage(stderr);
  } else {
    ProfileFile = argv[optind];
  }
  
  if (ReadProfile(ProfileFile, &prf, false, false) < 0) {
     FreeProfile(&prf, false);
    return 1;
  }
  
  if (prf.Type == PF_PATTERN) {
    FreeProfile(&prf, false);
    return 0;
  }
  
  if (out_clean) {
    const char Footer[] = "CC   Tests output with pfdump\n";
    WriteProfile(ProfileFile, &prf, Footer, stdout);
    FreeProfile(&prf, false);
    return 0;
  }
  
  if (isReversed) {
    rprf = ReverseProfile(&prf);
    if (rprf == NULL) {fputs("Error while reversing profile\n",stderr); exit(1);}
    pprf = rprf;
  } else {
    pprf = &prf;
  }
  PrepareExtraTable(pprf);
  
  if (out_alphabet) {
    puts("Alphabet Mapping");
    for (size_t i=0; i<ALPHABET_SIZE; ++i) {
      printf("Map %c=%2u\t", (char) ((unsigned char) 'A' + (unsigned char) i), (unsigned int) pprf->Alphabet_Mapping[i]);
      if ((i+1) % 8 == 0 ) puts("");
    }
    puts("\n");
  }
  const size_t prfLength = pprf->Length + 1;
  
	if (out_sequence) {
		printf("SEQ: \"%s\"", pprf->Sequence);
	}
	
  struct SMatch * const Match         = &pprf->Scores.Match;
  struct SInsertion * const Insertion = &pprf->Scores.Insertion;
  
  if (out_match) {
    printf("Match matrix with alignment %lu\n\n",Match->AlignStep );
    printf("    | ");
    for (size_t alpha=0; alpha<ALPHABET_MEMORY_SIZE; ++alpha) {
      printf(SPACE "%2lu ", alpha);
    }
    fputs("\n", stdout);
    printf("    | ");
    for (size_t alpha=0; alpha<ALPHABET_MEMORY_SIZE; ++alpha) {
      fputs(LINE, stdout);
    }
    fputs("\n", stdout);
    
    for (size_t iprf=0; iprf<prfLength; ++iprf) {
      const StoredIntegerFormat * MatchLine = &Match->Alphabet[iprf*Match->AlignStep];
      printf("%3lu | ", iprf+1);
      for (size_t alpha=0; alpha<ALPHABET_MEMORY_SIZE; ++alpha) {
				if (MatchLine[alpha] == NLOW) {
					printf(NLOW_FORMAT " ");
				} else {
					printf(INT_FORMAT " ", MatchLine[alpha]);
				}
      }
      fputs("\n", stdout);
    }
  }
//   puts("Transpose Match matrix");
//   const int * TIMatch = TransposeAndConvertMatchMatrix(&(pprf->Scores.Match), pprf->Alphabet_Length, pprf->Length);
//   const int * MatchLine = TIMatch;
//   const size_t Aligned_Profile_Length = (prfLength+1 + 15) & ~15;
//   for (size_t alpha=0; alpha<pprf->Alphabet_Length; ++alpha) {
//       printf("%3lu | ", alpha);
//       for (size_t iprf=0; iprf<prfLength; ++iprf) {
// 	if (MatchLine[iprf] == NLOW) {
// 	  printf("NLOW ");
// 	} else {
// 	  printf("%4i ", MatchLine[iprf]);
// 	}
//       }
//       fputs("\n", stdout);
//       MatchLine += Aligned_Profile_Length;
//   }
  if (out_insertion) {
    struct SInsertion Insertion = pprf->Scores.Insertion;
    printf("\nInsertion alphabet matrix with alignment %lu\n\n",Insertion.AlignStep );
    printf("    | ");
    for (size_t alpha=0; alpha<ALPHABET_MEMORY_SIZE; ++alpha) {
      printf(SPACE "%2lu ", alpha);
    }
    fputs("\n", stdout);
    printf("    | ");
    for (size_t alpha=0; alpha<ALPHABET_MEMORY_SIZE; ++alpha) {
      fputs(LINE, stdout);
    }
    fputs("\n", stdout);

    for (size_t iprf=0; iprf<prfLength; ++iprf) {
      StoredIntegerFormat * InsertionLine = &Insertion.Alphabet[iprf*Insertion.AlignStep];
      printf("%3lu | ", iprf+1);
      for (size_t alpha=0; alpha<ALPHABET_MEMORY_SIZE; ++alpha) {
				if (InsertionLine[alpha] == NLOW) {
					printf(NLOW_FORMAT " ");
				} else {
					printf(INT_FORMAT " ", InsertionLine[alpha]);
				}
      }
      fputs("\n", stdout);
    }
  }
  
  if (out_boundaries) {
    printf("\nInsertion boundaries matrix with alignment %i\n\n", INSERTION_BOUNDARIES_SIZE );
    char Header[] = "     " SPACE "_B0" SPACE "_B1" SPACE "_E0" SPACE "_E1" SPACE "_BM" SPACE "_BI" SPACE "_BD" SPACE "_BE" SPACE "_ME" SPACE "_IE" SPACE "_DE\n";
    
    fputs(Header, stdout);
    fputs("    |", stdout);
    for (size_t i=0; i<strlen(Header)-5; ++i) fputc('-', stdout);
    fputs("\n", stdout);
    for (size_t iprf=0; iprf<prfLength; ++iprf) {
      StoredIntegerFormat * InsertionLine = &Insertion->Boundaries[iprf*INSERTION_BOUNDARIES_SIZE];
      printf("%3lu | ", iprf);
      for (size_t alpha=0; alpha<INSERTION_BOUNDARIES_SIZE; ++alpha) {
				if (InsertionLine[alpha] == NLOW) {
					printf(NLOW_FORMAT " ");
				} else {
					printf(INT_FORMAT " ", InsertionLine[alpha]);
				}
      }
      fputs("\n", stdout);
    }
  }

  if (out_transitions) {
    const TransitionScores * restrict InsertionLine = pprf->Scores.Insertion.Transitions;
    const ScoreTuple * const restrict FirstSeq      = pprf->Scores.Insertion.FirstSequenceProtein;
    const ScoreTuple * const restrict LastSeq       = pprf->Scores.Insertion.LastSequenceProtein;
    
#define PRINT_VALUE(x) { \
      if (x == NLOW) {\
	printf(NLOW_FORMAT " ");\
      } else if ( x == NLOW/4*3) {\
	printf(MLOW_FORMAT " ");\
      } else {\
	printf(INT_FORMAT " ", x);\
      }\
    }
    
    printf("\nInsertion transition matrix with alignment %i\n\n", INSERTION_TRANSITIONS_SIZE );
    char Header2[] = "    |" SPACE "MATCH " SPACE SPACE "   " SPACE "    | "\
										SPACE "INSERTION" SPACE SPACE SPACE "   | "\
										SPACE "DELETION  " SPACE  SPACE SPACE "  | "\
										SPACE "EXTRA   " SPACE SPACE " |"\
										SPACE "Lst" SPACE "Seq" SPACE "OUT" " |" \
										SPACE "Fst" SPACE "Seq" SPACE "IN " " |\n";
    fputs(Header2,stdout);
    fputs("    |", stdout);
    for (size_t i=0; i<strlen(Header2)-7; ++i) fputc('-', stdout);
    fputs("|\n", stdout);
    char Header3[] = "    |" SPACE "_MM" SPACE "_MI" SPACE "_MD" SPACE "_MX" " |"\
                             SPACE "_IM" SPACE "_II" SPACE "_ID" SPACE "_IX" " |"\
                             SPACE "_DM" SPACE "_DI" SPACE "_DD" SPACE "_DX" " |"\
                             SPACE "_XM" SPACE "_XI" SPACE "_XD" " |"\
                             SPACE "_MY" SPACE "_IY" SPACE "_DY" " |"\
                             SPACE "_YM" SPACE "_YI" SPACE "_YD" " |\n";
    fputs(Header3,stdout);
    fputs("    |", stdout);
    for (size_t i=0; i<strlen(Header2)-7; ++i) fputc('-', stdout);
    fputs("|\n", stdout);

    for (size_t iprf=0; iprf<prfLength; ++iprf) {
      printf("%3lu | ", iprf);
      PRINT_VALUE(InsertionLine[iprf].From[MATCH].To[MATCH]);
      PRINT_VALUE(InsertionLine[iprf].From[MATCH].To[INSERTION]);
      PRINT_VALUE(InsertionLine[iprf].From[MATCH].To[DELETION]);
      PRINT_VALUE(InsertionLine[iprf].From[MATCH].To[EXTRA]);
      fputs("| ", stdout);
      PRINT_VALUE(InsertionLine[iprf].From[INSERTION].To[MATCH]);
      PRINT_VALUE(InsertionLine[iprf].From[INSERTION].To[INSERTION]);
      PRINT_VALUE(InsertionLine[iprf].From[INSERTION].To[DELETION]);
      PRINT_VALUE(InsertionLine[iprf].From[INSERTION].To[EXTRA]);
      fputs("| ", stdout);
      PRINT_VALUE(InsertionLine[iprf].From[DELETION].To[MATCH]);
      PRINT_VALUE(InsertionLine[iprf].From[DELETION].To[INSERTION]);
      PRINT_VALUE(InsertionLine[iprf].From[DELETION].To[DELETION]);
      PRINT_VALUE(InsertionLine[iprf].From[DELETION].To[EXTRA]);
      
      fputs("| ", stdout);
      PRINT_VALUE(InsertionLine[iprf].From[EXTRA].To[MATCH]);
      PRINT_VALUE(InsertionLine[iprf].From[EXTRA].To[INSERTION]);
      PRINT_VALUE(InsertionLine[iprf].From[EXTRA].To[DELETION]);
      fputs("| ", stdout);
      
      PRINT_VALUE(LastSeq[iprf].From[MATCH]);
      PRINT_VALUE(LastSeq[iprf].From[INSERTION]);
      PRINT_VALUE(LastSeq[iprf].From[DELETION]);
      
      fputs("| ",stdout);
      
      PRINT_VALUE(FirstSeq[iprf].To[MATCH]);
      PRINT_VALUE(FirstSeq[iprf].To[INSERTION]);
      PRINT_VALUE(FirstSeq[iprf].To[DELETION]);
      
      fputs("|\n", stdout);
    }
  }
  
  FreeProfile(pprf, isReversed);
  return 0;
}


