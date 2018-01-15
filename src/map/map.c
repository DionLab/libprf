#include "prf_config.h"
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "pfMap.h"
#include "json.h"

static const char Alphabet[] = "XACGTN";
static const unsigned char Alphabet_Mapping[32] = {
	1, 0, 2, 0, 0, 0, 3, 0,
	0, 0, 0, 0, 0, 5, 0, 0,
	0, 0, 0, 4, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0
};

static const struct Map DefaultMap = {
	.RowInit = { 
		.elem = {
			( _BD << SCORE_SHIFT) | (PRIORITY_EXTRA << STATE_SHIFT),
			(NLOW << SCORE_SHIFT) | (PRIORITY_EXTRA << STATE_SHIFT),
			( _BM << SCORE_SHIFT) | (PRIORITY_EXTRA << STATE_SHIFT),
			(NLOW << SCORE_SHIFT) | (PRIORITY_EXTRA << STATE_SHIFT)
		}
	},
	.MatchScores = {
		.elem = {
			((_M+_MD) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (1 << MATCH_SHIFT),
			((_M+_MX) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (1 << MATCH_SHIFT),
			((_M+_MM) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (1 << MATCH_SHIFT),
			((_M+_MI) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (1 << MATCH_SHIFT)
		}
	},
	.MismatchScores = {
		.elem = {
			((_m+_MD) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (0 << MATCH_SHIFT),
			((_m+_MX) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (0 << MATCH_SHIFT),
			((_m+_MM) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (0 << MATCH_SHIFT),
			((_m+_MI) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (0 << MATCH_SHIFT)
		}
	},
	.SequenceNMatchScores = {
		.elem = {
			((_MD) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (1 << MATCH_SHIFT),
			((_MX) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (1 << MATCH_SHIFT),
			((_MM) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (1 << MATCH_SHIFT),
			((_MI) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (1 << MATCH_SHIFT)
		}
	},
	.InsertionScores = {
		.elem = {
			((NLOW  ) << SCORE_SHIFT) + (PRIORITY_INSERTION << STATE_SHIFT) + (0 << MATCH_SHIFT),
			((_I+_IX) << SCORE_SHIFT) + (PRIORITY_INSERTION << STATE_SHIFT) + (0 << MATCH_SHIFT),
			((_I+_IM) << SCORE_SHIFT) + (PRIORITY_INSERTION << STATE_SHIFT) + (0 << MATCH_SHIFT),
			((_I+_II) << SCORE_SHIFT) + (PRIORITY_INSERTION << STATE_SHIFT) + (0 << MATCH_SHIFT)
		}
	},
	.DeletionScores = {
		.elem = {
			((_D+_DD) << SCORE_SHIFT) + (PRIORITY_DELETION << STATE_SHIFT) + (0 << MATCH_SHIFT),
			((_D+_DX) << SCORE_SHIFT) + (PRIORITY_DELETION << STATE_SHIFT) + (0 << MATCH_SHIFT),
			((_D+_DM) << SCORE_SHIFT) + (PRIORITY_DELETION << STATE_SHIFT) + (0 << MATCH_SHIFT),
			((NLOW  ) << SCORE_SHIFT) + (PRIORITY_DELETION << STATE_SHIFT) + (0 << MATCH_SHIFT)
		}
	},
	.Alphabet = Alphabet,
	.Alphabet_Mapping = Alphabet_Mapping,
	.MatchScore = _M,
	.MismatchScore = _m,
	.NScore = 0
};

const struct Map * getDefaultMap()
{
	return &DefaultMap;
}

int dumpDefaultMap(const char * FileName)
{
	FILE* const restrict out = fopen(FileName, "w");
	if (out) {
		fprintf(out,"{\n" 
		             "/* --------- Aligner scores ---------- */\n"
								 "\t\"Begin\": {\n"
								 "\t/* Begin scores to */\n"
								 "\t\t\"Match\": %i,\n"
								 "\t\t\"Insertion\": %i,\n"
								 "\t},\n"
		             "\t\"Scores\": {\n"
		             "\t/* State scores */\n"
		             "\t\t\"Match\": %i,\n"
		             "\t\t\"MisMatch\": %i,\n"
		             "\t\t\"Insertion\": %i,\n"
		             "\t\t\"Deletion\": %i\n"
		             "\t},\n"
		             "/* State transition scores, note that minimal value (NLOW is %i) */\n"
		             "\t\"Transitions\": {\n"
		             "\t\t\"MM\": %i, \"MI\": %i, \"MD\": %i,\n"
		             "\t\t\"IM\": %i, \"II\": %i, \"ID\": %i,\n"
		             "\t\t\"DM\": %i, \"DI\": %i, \"DD\": %i\n"
		             "\t}\n"
						     "}\n",
					       _BM, _BD,
		             _M, _m, _I, _D, NLOW,
					       _MM, _MI, _MD,
		             _IM, _II, _ID,
		             _DM, _DI, _DD);
		fclose(out);
		return 0;
	}
	else {
		fprintf(stderr, "Unable to generate JSON file %s\n", FileName);
		return -1;
	}
}

int loadMap(const char * const restrict FileName, struct Map * const restrict map)
{
	int err = -1;
	union {
		struct elem {
			int m;
			int M;
			int In;
			int D;
		} elem;
		__m128i xmm;
	} Scores;
	union {
		struct To {
			int D;
			int X;
			int M;
			int In;
		} To;
		__m128i xmm;
	} FromM, FromI, FromD;
	int BM, BD;
	
	memset(map, 0 , sizeof(struct Map));
// 	map->Alphabet = Alphabet;
	
	json_object * const jBase = json_object_from_file(FileName);
	if (! jBase) {
		fprintf(stderr,"Error loading JSON file %s\n", FileName);
		goto bail;
	}
	
	json_object * jElement;
	
#define LOAD(string, where) {\
	json_object_object_get_ex(jElement, string, &jSubElement);\
	if (jSubElement) {\
		if (json_object_is_type(jSubElement, json_type_int)) {\
			where = json_object_get_int(jSubElement);\
		}\
		else {\
			fputs("Unable to read correct format for " string "\n", stderr);\
			goto bail;\
		}\
	}\
	else {\
		fprintf(stderr, "Unable to find " string " within JSON file %s\n", FileName);\
		goto bail;\
	}\
}
	
	json_object_object_get_ex(jBase, "Begin", &jElement);
	if (jElement == NULL) {
		fprintf(stderr, "Unable to find Begin within JSON file %s\n", FileName);
		goto bail;
	}
	else {
		json_object * jSubElement;
		LOAD("Match", BM);
		LOAD("Deletion", BD);
	}
	
	json_object_object_get_ex(jBase, "Scores", &jElement);
	if (jElement == NULL) {
		fprintf(stderr, "Unable to find Scores within JSON file %s\n", FileName);
		goto bail;
	}
	else {
		json_object * jSubElement;
		LOAD("Match",  Scores.elem.M);
		LOAD("MisMatch", Scores.elem.m);
		LOAD("Insertion", Scores.elem.In);
		LOAD("Deletion", Scores.elem.D);
	}
	
	json_object_object_get_ex(jBase, "Transitions", &jElement);
	if (jElement == NULL) {
		fprintf(stderr, "Unable to find Transitions within JSON file %s\n", FileName);
		goto bail;
	}
	else {
		json_object * jSubElement;
		LOAD("MM", FromM.To.M);
		LOAD("MI", FromM.To.In);
		LOAD("MD", FromM.To.D);
		
		LOAD("IM", FromI.To.M);
		LOAD("II", FromI.To.In);
		LOAD("ID", FromI.To.D);
		
		LOAD("DM", FromD.To.M);
		LOAD("DI", FromD.To.In);
		LOAD("DD", FromD.To.D);
	}
#undef LOAD

	map->RowInit.elem[MATCH] = BM;
	map->RowInit.elem[DELETION] = BD;
	map->RowInit.elem[INSERTION] = NLOW;
	map->RowInit.elem[EXTRA] = NLOW;
	
	__m128i __m = _mm_shuffle_epi32(Scores.xmm, 0b00000000);
	__m128i __M = _mm_shuffle_epi32(Scores.xmm, 0b01010101);
	__m128i __I = _mm_shuffle_epi32(Scores.xmm, 0b10101010);
	__m128i __D = _mm_shuffle_epi32(Scores.xmm, 0b11111111);
	
	FromM.To.X = NLOW;
	FromI.To.X = NLOW;
	FromD.To.X = NLOW;
	
	__m128i __Transitionm = _mm_add_epi32(FromM.xmm, __m);
	__m128i __SequenceN   = FromM.xmm;
	FromM.xmm = _mm_add_epi32(FromM.xmm, __M);
	FromI.xmm = _mm_add_epi32(FromI.xmm, __I);
	FromD.xmm = _mm_add_epi32(FromD.xmm, __D);
	
	__Transitionm = _mm_slli_epi32(__Transitionm, SCORE_SHIFT);
	__SequenceN   = _mm_slli_epi32(__SequenceN, SCORE_SHIFT);
	FromM.xmm     = _mm_slli_epi32(FromM.xmm, SCORE_SHIFT);
	FromI.xmm     = _mm_slli_epi32(FromI.xmm, SCORE_SHIFT);
	FromD.xmm     = _mm_slli_epi32(FromD.xmm, SCORE_SHIFT);
	
	map->MismatchScores.xmm  = _mm_or_si128(__Transitionm, _mm_set1_epi32((PRIORITY_MATCH<<STATE_SHIFT) || (0 << MATCH_SHIFT)));
	map->MatchScores.xmm     = _mm_or_si128(FromM.xmm, _mm_set1_epi32((PRIORITY_MATCH << STATE_SHIFT) || (1 << MATCH_SHIFT)));
	map->SequenceNMatchScores.xmm = _mm_or_si128(__SequenceN, _mm_set1_epi32((PRIORITY_MATCH << STATE_SHIFT) || (1 << MATCH_SHIFT)));
	map->InsertionScores.xmm = _mm_or_si128(FromI.xmm, _mm_set1_epi32((PRIORITY_INSERTION << STATE_SHIFT) || (0 << MATCH_SHIFT)));
	map->DeletionScores.xmm  = _mm_or_si128(FromD.xmm, _mm_set1_epi32((PRIORITY_DELETION << STATE_SHIFT) || (0 << MATCH_SHIFT)));
	
err = 0;
	
	bail:;
	return err;
}

