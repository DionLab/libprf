#include "prf_config.h"
#include <stdlib.h>
#include "pfCompute.h"

#if !defined(__USE_WINAPI__)
pthread_mutex_t PrintLock;
#else
CRITICAL_SECTION PrintLock;
#endif

int compareAlignments(const void *a, const void *b) 
{
		if ( ((Alignment_t*)a)->Region.sequence.Begin < ((Alignment_t*)b)->Region.sequence.Begin )
			return -1;
		else {
			if ( ((Alignment_t*)a)->Region.sequence.Begin == ((Alignment_t*)b)->Region.sequence.Begin ) 
				return 0;
			else 
				return 1;
		}
}
