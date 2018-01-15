#include "prf_config.h"
#include <stdlib.h>
#include "pfProfile.h"
#include "pfOutput.h"
/****************************************************************************************************
 * From the following enums, we define the ordering as follow
 * 
enum SearchType { PRF=1, REGEX=2 };
enum OutputType { TEXT, HISTOGRAM, DENSITY, PNG, DATA, PDF };
enum Constrain { WITH_REVERSE=1, SCORE_RANGE=2, CYCLE_RANGE=4, INVERSE_SELECTION=8, BORDER_CLIP=16};
enum WrapperRange { Wrapper_Alignement = 0, Wrapper_Source = 1, Wrapper_Before = 2, Wrapper_After = 3,
	                  Wrapper_Outer = 4};

	MATRIX PRF:
		TEXT:																								: 160 cases
			- WITH_REVERSE
			- SCORE_RANGE
			- CYCLE_RANGE if circular profile only
			- INVERSE_SELECTION
			- BORDER_CLIP
			- WRAPPER (5 cases)
			
		TEXT2_OPTIMAL_WITH_REVERSE:													:  16 cases		
			- SCORE_RANGE
			- CYCLE_RANGE if circular profile only
			- INVERSE_SELECTION
			- BORDER_CLIP
			
		HISTOGRAM:																					:  16 cases
			- WITH_REVERSE
			- SCORE_RANGE 
			- CYCLE_RANGE if circular profile only
			- BORDER_CLIP
			- HISTO_CYCLE
			
		DENSITY: if circular profile												:  2 cases
			- WITH_REVERSE
			- SCORE_RANGE required anyway
			- CYCLE_RANGE required anyway
	
		PDF, PNG, DATA:																			: 2 cases
			- WITH_REVERSE
	
	
	REGEX or PATTERN PRF:
		TEXT:																								: 32 cases
			- WITH_REVERSE
			- SCORE_RANGE
			- CYCLE_RANGE
			- INVERSE_SELECTION
			- BORDER_CLIP
		
		HISTOGRAM:																					:  6 cases
			- WITH_REVERSE
			- SCORE_RANGE 
			- CYCLE_RANGE if circular profile only
			
		DENSITY:																						:  2 cases
			- WITH_REVERSE
			- SCORE_RANGE required anyway
			- CYCLE_RANGE required anyway
*****************************************************************************************************/
#define N_PRF_TEXT								      160
#define N_PRF_TEXT_OPTIMAL_WITH_REVERSE  16
#define N_PRF_HISTOGRAM								   16
#define N_PRF_DENSITY								      2
#define N_PRF_OTHER								        2
/*****************************************************************************************************/
int GetDispatchThreadIndex(const OutputType_t * const restrict output,
                           const struct Profile * const restrict prf) 
{
	int res = (output->SearchWith & PRF) ? 0 : (N_PRF_DENSITY + N_PRF_TEXT_OPTIMAL_WITH_REVERSE + N_PRF_HISTOGRAM + N_PRF_TEXT + N_PRF_OTHER);
	
	switch(output->Type) {
		case TEXT:
			if (output->Specific.Text.BestOfStdAndRevComp)
				res += 160 + (output->Constrains >> 1);
			else
				res += 32*output->Specific.Text.Range + output->Constrains;;
			break;
		case HISTOGRAM:
			{
				if (output->Specific.Histogram.CycleRatherThanScore && !(output->Constrains & CYCLE_RANGE)) {
					fprintf(stderr, "Histogram on cycles requires --cycle-range\n");
					goto bail;
				}
				unsigned int uitmp = 0U;
				if (output->Constrains & WITH_REVERSE) uitmp |= 1; 
				if (output->Constrains & BORDER_CLIP) uitmp |= 2;
				
				if (output->Specific.Histogram.CycleRatherThanScore) {
					if (output->Constrains & SCORE_RANGE) uitmp |= 4;
					res += N_PRF_TEXT + N_PRF_TEXT_OPTIMAL_WITH_REVERSE + 8 + uitmp;
				}
				else {
					if (output->Constrains & CYCLE_RANGE) uitmp |= 4;
					res += N_PRF_TEXT + N_PRF_TEXT_OPTIMAL_WITH_REVERSE + uitmp; 
				}
			}
			break;
		case DENSITY:
			if (output->SearchWith != PRF) goto bail;
			res += N_PRF_TEXT + N_PRF_TEXT_OPTIMAL_WITH_REVERSE + N_PRF_HISTOGRAM + (output->Constrains & 0b1);
			break;
		case PNG:
		case DATA:
		case PDF: 
			if ((output->SearchWith != PRF) || (prf && !(prf->Type == PF_MATRIX))) goto bail;
			res += N_PRF_TEXT + N_PRF_TEXT_OPTIMAL_WITH_REVERSE + N_PRF_HISTOGRAM + N_PRF_DENSITY + (output->Constrains & 0b1);
			break;
		case TEST:
		default:
			fputs("Oops I missed that thread configuration!!!\n", stderr);
			exit(1);
	}
	
	return res;
	
	bail: ;
		fputs("GetDispatchThreadIndex: Impossible choice to accomplish, verify your options\n", stderr);
		return -1;
}
