
/***************************************************************************************************
                                        PROFILING
 ***************************************************************************************************
  Apr 1, 2016 profiling.c
 ***************************************************************************************************
 (C) 2016 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 ***************************************************************************************************/
#define _GNU_SOURCE
#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <pthread.h>
#include <stdint.h>
#include <stdbool.h>
#include <getopt.h>
# include <unistd.h>
#ifdef PRF_USE_AFFINITY
# include <sched.h>
#endif
#include <nmmintrin.h>

#define HEADER "%----------------------------------------------------------------------------%\n"\
	       "|                               Mapping v" PRF_VERSION "                         |\n"\
	       "%----------------------------------------------------------------------------%\n"\
	       "| Built on " __DATE__ " at " __TIME__ ".                                          |\n"
#include "pfVersion.h"
#define __USE_INLINE_FUNCTIONS__
#include "pfMap.h"
#include "pfInput.h"
#include "pfSequence.h"
#include "threadpool.h"

typedef struct Reference_s {
// 	pthread_mutex_t lock;
	Sequence_t * Sequences;
	size_t Count;
} Reference_t;
 
typedef struct MapJob_s {
	struct MapJob_s * prev;
	size_t id;
} MapJob_t;
				 
#ifdef PRF_USE_AFFINITY
#define __USE_AFFINITY__
#endif
#include "system.h"
static SystemInfo System;
static const struct Map * map; 												/* Map */
static FASTAStructure FASTA_DB1, FASTA_DB2;			/* Sequence Database File */
static Reference_t References;

#ifdef PRF_USE_AFFINITY
cpu_set_t * Thread_masks[2] = {0,0};						/* Define variables to hold thread affinity mask */
unsigned int Thread_count[2] = {0,0};
pthread_attr_t * restrict threads_attr = NULL;
#endif

#define TEST_IF_AVAILABLE(TYPE) if (OutputType.Type != TYPE) {\
	fprintf(stderr,"Invalid option for this type");\
	exit(1);\
}

static const char opt_to_test[] = "hsVI:i:t:"
#ifdef PRF_USE_AFFINITY
	"012:3"
#endif
;

_Bool OutputVerbose __attribute__((weak)) = false;


static const struct option long_options[] =
{
  /*
	 * These options set a flag.
	 */

  /*
	 * These options don't set a flag. We distinguish them by their indices.
	 */
	{"help",               		no_argument,       	0,	'h'},
	{"sse2",									no_argument,				0,	's'},
	{"verbose",								no_argument,				0,	'V'},
	/* Sequence */
	/* Database indexing options */
	{"index-database-1",			required_argument,	0,	'I'},
	{"index-database-2",			required_argument,	0,	'i'},
	/* Others */
	/* SMP options*/
	{"nthreads",							required_argument,	0,	't'},
#ifdef PRF_USE_AFFINITY
	{"no-affinity",						no_argument,				0,	'0'},
	{"split", 								no_argument,				0,	'1'},
	{"thread-affinity",				required_argument,	0,	'2'},
	{"no-shared-core",				no_argument,				0,	'3'},
#endif
	{0, 0, 0, 0}
};

static void __attribute__((noreturn)) Usage(FILE * stream)
{
  fputs(
	" Mapping [options] Database1 Database2\n"
	" Options:\n"
	"  Database\n"
	"   FASTA\n"
	"     --index-database-1       [-I] : use indices stored in given file, create if missing\n"
	"     --index-database-2       [-i] : use indices stored in given file, create if missing\n\n"
	"  Optimizations\n"
// 	"   --sse2                     [-s] : enforces SSE 2 only instruction set\n"
	"   --nthreads                 [-t] : max number of threads to use\n"
#ifdef PRF_USE_AFFINITY
	"   --no-affinity                   : disable CPU affinity file\n"
	"   --thread-affinity               : file containing thread mask,\n"
	"                                     one row for one thread\n"
	"   --no-shared-core                : Prevent core resource sharing\n"
	"   --split                         : if both SSE 2 & 4.1 are available,\n"
	"                                     split half-half using linked resources\n\n"
#else
	"\n"
#endif
	"  Other\n"
	"    --verbose                 [-V] : verbose on stderr\n"
	"    --help                    [-h] : output command help\n\n"
	"This is version " PRF_VERSION ".\n",
	stream);
  exit(1);
}

static void* alignment_thread(threadarg_t * const restrict _Data)
{
	int res = SUCCESS;
	Sequence_t SeqData;
	/*************************************************************************/
	/*                          GET COMMON DATA                              */
	/*************************************************************************/

	
	/*************************************************************************/
	/*                          ALLOCATE MEMORY                              */
	/*************************************************************************/
	const size_t WorkSize = FASTA_DB1.MaxSequenceSize + 1;
	SeqData.Size = FASTA_DB1.MaxSequenceSize;
	SeqData.Data.Memory = (unsigned char*) _mm_malloc(FASTA_DB1.MaxSequenceSize*sizeof(unsigned char),64);
	int * restrict WORK = _mm_malloc(2*WorkSize*sizeof(__m128i), 64);
	if ( WORK == NULL || SeqData.Data.Memory == NULL ) {
		res = -1; 
		goto FIN;
	}
	
	/*************************************************************************/
	/*                          CONFIGURE WORKER                             */
	/*************************************************************************/
		
	/* Assure all threads have been created before starting serving */
	threadpool_t* const restrict  thpool = _Data->thpool;
	
	/* Mark thread as alive (initialized) */
	pthread_mutex_lock(&thpool->thcount_lock);
	thpool->num_threads_alive += 1;
	pthread_mutex_unlock(&thpool->thcount_lock);

	// Initialize some other data
	MapJob_t * task = NULL;
	
	while(thpool->threads_keepalive) {
		bsem_wait(thpool->jobqueue_p.has_items);
		
		if (thpool->threads_keepalive){
			pthread_mutex_lock(&thpool->thcount_lock);
			thpool->num_threads_working++;
			pthread_mutex_unlock(&thpool->thcount_lock);
			
			/* Read job from queue and execute it */
			pthread_mutex_lock(&thpool->jobqueue_p.rwmutex);
			task = (MapJob_t*) jobqueue_pull(&thpool->jobqueue_p);
			pthread_mutex_unlock(&thpool->jobqueue_p.rwmutex);
			
			/***************************************************************/
			/*                     START THE TASK                          */
			/***************************************************************/
			if (task) {
				/* Get the target sequence */
				SETUP_DATABASE_ACCESS(FASTA_DB1.SequenceFile);
				PFSequence * const Seq = GET_DATABASE_SEQUENCE(&(SeqData), FASTA_DB1.DataPtr, task->id);
				UNSET_DATABASE_ACCESS();
				
				/* Transform into real Fasta sequence and split header and data */
				if (SeqData.ProfileData.Length > 0UL) {
					CleanSequence(&(SeqData.ProfileData));
					
					/////////////////////////////////////////////////////////////////////////////////////////////////////////
					// WORK
					for (size_t iRef=0UL; iRef<References.Count; iRef++) {
						Zone_t zone;
						GetMapping(map, WORK,
                       (const char*) References.Sequences[iRef].ProfileData.ProfileIndex,
						           (const char*) SeqData.ProfileData.ProfileIndex,
                       &zone,
						           SeqData.ProfileData.Length,
						           References.Sequences[iRef].ProfileData.Length);
						fprintf(stdout,"Sequence %lu : %s %u-%u %i\n", task->id, References.Sequences[iRef].Data.Header, zone.Begin, zone.End, zone.Score);
					}
				}
				
				/* Return the job memory slot */
				pthread_mutex_lock(&thpool->donequeue_p.rwmutex);
				jobqueue_push(&thpool->donequeue_p, (job_t*) task);
				pthread_mutex_unlock(&thpool->donequeue_p.rwmutex);
			}
			pthread_mutex_lock(&thpool->thcount_lock);
			thpool->num_threads_working--;
			pthread_mutex_unlock(&thpool->thcount_lock);
		}
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//SUMMARIZE;
		
	FIN:
	pthread_mutex_lock(&thpool->thcount_lock);
	thpool->num_threads_alive--;
	pthread_mutex_unlock(&thpool->thcount_lock);
			
	if (WORK) _mm_free(WORK);
	if (SeqData.Data.Memory) _mm_free(SeqData.Data.Memory);
	
	if( res != SUCCESS) {
		fprintf(stderr, "Thread %u error code %i\n", _Data->threadID, res);
		if (task) {
			pthread_mutex_lock(&thpool->thcount_lock);
			thpool->num_threads_working--;
			pthread_mutex_unlock(&thpool->thcount_lock);
		}
		thpool->threads_keepalive = 0;
	}
	
	_Data->returnState = (void*) (uintptr_t) res;
	
	return 0;
}

int main (int argc, char *argv[])
{
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // LOCAL STRUCTURES
  ////////////////////////////////////////////////////////////////////////////////////////////////
  struct timeval _t0, _t1;							/* Timing structures */

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // LOCAL DATA
  ////////////////////////////////////////////////////////////////////////////////////////////////

  size_t nCPUs=0;										    /* number of threads */
  int res;
  char * DB_1 = NULL, * DB_2 = NULL;		/* FASTA sequence file */

  enum Version ComputeVersion = SSE41;	/* Trigger SSE version to use for filter and alignment */
#ifdef PRF_USE_AFFINITY
  char buffer[128] __attribute__((aligned(16)));	/* buffer to read affinity file mask */
  _Bool noAffinity = false;												/* disable use of cpu affinity */
  _Bool split = false;
  _Bool noSharedCore = false;					/* Prevent hyperthreading or AMD compute unit to share resources */
  _Bool GivenAffinityFile = false;		/* File holding a mask for each thread */
  char * AffinityMaskFileName;				/* Name of affinity mask file provided by option m */
#endif

	////////////////////////////////////////////////////////////////////////////////////////////////
  // INITIALIZE VALUES
  ////////////////////////////////////////////////////////////////////////////////////////////////
 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // SYSTEM ARCHITECTURE ANALYSIS
  ////////////////////////////////////////////////////////////////////////////////////////////////
  getSystemInfo(&System);

  /* Check for minimum requirement */
  if (!(System.Extensions & MM_SSE2)) {
      fputs("pfrepeat requires at least a CPU capable of SSE 2.\n", stderr);
      exit(1);
  }

  /* Allow fast SSE 4.1 extensions ? */
  if (System.Extensions & MM_SSE41) {
      ComputeVersion = SSE41;
  }
  else {
      fputs("Sorry only compute methods based upon SSE 4.1 are available yet\n", stderr);
      exit(1);
  }
  
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
			case 'V':
				OutputVerbose = true;
				break;
      case 'h':
				Usage(stderr);
				break;
    }
  }
  
	if (OutputVerbose) {
   fputs(HEADER
#ifdef __USE_MMAP__
	       "| Using Linux kernel MMAP function.                                          |\n"
#endif
				 "%----------------------------------------------------------------------------%\n"
      ,stderr);
    printSystemInfo(&System);
#ifdef USE_32BIT_FORMAT
    fputs("Using 32 bit format integer for scores\n", stderr);
#endif
    if (ComputeVersion == SSE2 && (System.Extensions & MM_SSE41)) {
			fputs("Enforcing SSE 2...\n", stderr);
    }
  }
  
  optind = 1;
  while (1) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    const int c = getopt_long (argc, argv, opt_to_test, long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) break;
    switch (c) {
			/****************************************** Sequence *****************************************/
#ifdef PRF_USE_AFFINITY
      case '0':
        noAffinity = true;
        break;
      case '3':
				noSharedCore = true;
				break;
      case '1':
				if (ComputeVersion == SSE41) {
					split = true;
				} else {
					fputs("Split not possible without SSE 4.1\n", stderr);
					exit(1);
				}
				break;
      case '2':
				GivenAffinityFile = true;
				AffinityMaskFileName = optarg;
				break;
#endif
      case 'V':
				break;
      case 't':
				nCPUs = (size_t) atoi(optarg);
				break;
      case 's':
				ComputeVersion = SSE2;
				break;
      case 'I':
				FASTA_DB1.Options |= DoFastaIndexExport;
				FASTA_DB1.indexFileName = optarg;
				break;
      case 'i':
				FASTA_DB2.Options |= DoFastaIndexImport;
				FASTA_DB2.indexFileName = optarg;
				break;
      case 'h':
				Usage(stderr);
				break;
			case '?':
				fprintf(stderr, "Unknown option %s\n", argv[option_index]);
				exit(1);
			default:
				fprintf(stderr,"Option %c is unknown\n", c);
    }
  }

  if (optind >= argc) {
    fputs("Expected arguments after options\n", stderr);
    Usage(stderr);
  }
  else {
    DB_1 = argv[optind];
		DB_2 = argv[optind+1];
  }
    
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // COMMAND LINE COHERENCY
  ////////////////////////////////////////////////////////////////////////////////////////////////
  {
		_Bool AllOK = true;
		
		if (!DB_1 || !DB_2) {
			fputs("Expected arguments after options\n", stderr);
			AllOK = false;
		}
  
		if (!AllOK) Usage(stderr);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
  // THREADING ANALYSIS
  ////////////////////////////////////////////////////////////////////////////////////////////////
  /*
   * Retrieve number of cores
   */
  nCPUs = (nCPUs == 0) ? (size_t) System.nOverallCores : nCPUs;

#ifdef PRF_USE_AFFINITY
  if (noAffinity) {
    // -----------------------------------------------------------------------------
    //                        ***  NO AFFINITY ***
    // -----------------------------------------------------------------------------
    if (OutputVerbose) fputs("Thread affinity disabled\n", stderr);
    Thread_count[0] = System.nOverallCores;
    Thread_masks[0] = (cpu_set_t*) malloc(System.nOverallCores*sizeof(cpu_set_t));
    for (size_t thread=0; thread<System.nOverallCores; ++thread) {
      CPU_ZERO(&Thread_masks[0][thread]);
      for (int i=0; i<(int) System.nOverallCores; ++i) CPU_SET(i, &Thread_masks[0][thread]);
    }
  }
  else if (GivenAffinityFile) {
    // -----------------------------------------------------------------------------
    //                     ***  INPUT FILE HOLDING NASKS ***
    // -----------------------------------------------------------------------------
    if (OutputVerbose)
      fprintf(stderr,"Parsing file %s for affinity mask and number of threads\n", AffinityMaskFileName);
    FILE* in = fopen(AffinityMaskFileName, "r");
    if (in == NULL) {
			fprintf(stderr, "Cannot open thread affinity file %s.\n", optarg);
			exit(1);
    }
    size_t lines = 0;
    while (!feof(in)) {
			int num = fread(buffer, sizeof(char), 64, in);
			for (unsigned int i=0; i<num; i++)
		    if (buffer[i] == '\n') lines++;
    }
    rewind(in);
    if (lines != 0) {
			if (lines > System.nOverallCores) lines = System.nOverallCores;
			Thread_masks[0] = (cpu_set_t*) malloc(lines*sizeof(cpu_set_t));
			for (size_t i=0; i<lines; i++) {
			    if (fscanf(in, "%s\n", buffer) == 1) {
			      const size_t tmp_size = strlen(buffer) - 1;
			      CPU_ZERO(&Thread_masks[0][i]);
			      for (int j=tmp_size; j>=0; j--) {
				  if (buffer[j] != '0') CPU_SET(j, &Thread_masks[0][i]);
			      }
			    }
			}
			Thread_count[0] = lines;
			if (OutputVerbose) fprintf(stderr,"Found %2lu threads affinity masks.",nCPUs);
    }
    else {
			if (OutputVerbose) printf("Cannot understand cpu mask, keep on normally\n");
    }
    fclose(in);
  }
  else if ( split ) {
    // -----------------------------------------------------------------------------
    //                 ***  HALF SSE 2 HALF SSE 4.1 HYPERTHREADING***
    // -----------------------------------------------------------------------------
    Thread_count[0] = getMasks(&System, -1, -1, 1, &Thread_masks[0]);
    if (Thread_count[0] == 0) {
      fputs("No potential affinity mask found !!!\n", stderr);
      exit(0);
    }
    Thread_count[1] = getMasks(&System, -1, -1, 2, &Thread_masks[1]);
    if (Thread_count[1] == 0) {
      fputs("No potential affinity mask found with hyperthreading !!!\n", stderr);
      exit(0);
    }
    if (OutputVerbose)
      fprintf(stderr, "%u threads will use SSE 4.1 and %u SSE 2\n", Thread_count[0], Thread_count[1]);
  }
  else if (noSharedCore) {
    if (OutputVerbose)
      fputs("No sharing of core resources will be used: Intel Hyperthreading or AMD Compute Unit\n", stderr);
    Thread_count[0] = getMasks(&System, -1, -1, 1, &Thread_masks[0]);
    if (Thread_count[0] == 0) {
      fputs("No potential affinity mask found !!!\n", stderr);
      exit(0);
    }
  }
  else {
    // -----------------------------------------------------------------------------
    //                        *** OPERATING SYSTEM CHOICE ***
    // -----------------------------------------------------------------------------
    Thread_count[0] = getMasks(&System, -1, -1, -1, &Thread_masks[0]);
    if (Thread_count[0] == 0) {
      fputs("No potential affinity mask found !!!\n", stderr);
      exit(0);
    }
  }

  {
    register size_t total = (size_t) (Thread_count[0] + Thread_count[1]);
    if (nCPUs > total) nCPUs = total;
  }
 
  threads_attr = (pthread_attr_t*) alloca(nCPUs*sizeof(pthread_attr_t));
  {
    register const cpu_set_t * current = &Thread_masks[0][0];
    for (size_t i=0; i<nCPUs; ++i) {
      pthread_attr_init(&threads_attr[i]);
      if (i == (size_t) Thread_count[0]) current = &Thread_masks[1][0];
      pthread_attr_setaffinity_np(&threads_attr[i], sizeof(cpu_set_t), current);
      ++current;
    }
  }
#endif
  
  if (OutputVerbose) {
			fprintf(stderr, "Job dispatched over %lu cores.\n", nCPUs);
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////
  // LOAD REFERENCE
  ////////////////////////////////////////////////////////////////////////////////////////////////
	if (isFASTA(DB_2)) {
		if (OpenFASTAStructure(DB_2, &FASTA_DB2) == 0) {
			SETUP_DATABASE_ACCESS(FASTA_DB2.SequenceFile);
			
			References.Count = FASTA_DB2.SequenceCount;
			References.Sequences = (Sequence_t*) malloc(FASTA_DB2.SequenceCount*sizeof(Sequence_t));
			if (References.Sequences == NULL) {
				fprintf(stderr, "Unable to allocate memory for Reference sequence\n");
				goto bail;
			}
			
			for (size_t iRef=0UL; iRef<FASTA_DB2.SequenceCount; iRef++) {
                const size_t Size = FASTA_DB2.DataPtr[iRef].HeaderLength + FASTA_DB2.DataPtr[iRef].SequenceLength + 2UL;
				/* Allocate memory if first pass */
				References.Sequences[iRef].Data.Memory = _mm_malloc(Size, 64);
				if (References.Sequences[iRef].Data.Memory == NULL) {
					fprintf(stderr, "Unable to allocate memory for Reference sequence\n");
					goto bail;
				}
				References.Sequences[iRef].Size = Size;
				
				/* Get Reference Sequence */
				PFSequence * const RefSeq = GET_DATABASE_SEQUENCE(&(References.Sequences[iRef]), FASTA_DB2.DataPtr, iRef);
				
				CleanSequence(RefSeq);
			}
			UNSET_DATABASE_ACCESS();
			CloseFASTAStructure(&FASTA_DB2);
		}
		else {
			fprintf(stderr, "Error opening FASTA file %s\n", DB_2);
			res = 1;
			goto bail;
		}
	}
	else {
		fprintf(stderr, "Reference database %s is not a FASTA file!\n", DB_2);
		goto bail;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////
  // LOAD SEQUENCES FILE FOR MEMORY IDEAS
  ////////////////////////////////////////////////////////////////////////////////////////////////
	if ( isFASTA(DB_1)) {
		if (OpenFASTAStructure(DB_1, &FASTA_DB1) != 0) {
			fprintf(stderr, "Error opening FASTA file %s\n", DB_1);
			res = 1;
			goto bail;
		}
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////
  // STARTING THREADPOOL
  ////////////////////////////////////////////////////////////////////////////////////////////////
  map = getDefaultMap();
  threadpool_t * const pool = createThreadPool(alignment_thread,
#ifdef PRF_USE_AFFINITY
	                                             Thread_masks[0],
#endif
	                                             sizeof(MapJob_t), NULL, nCPUs, 0);
	
	if (pool == NULL) {
		fprintf(stderr, "Unable to create thread pool\n");
		goto bail;
	}
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // ANALYSIS BASED ON DATABASE
  ////////////////////////////////////////////////////////////////////////////////////////////////
	/* ---------------------------------  FASTA SEQUENCES --------------------------------------- */
	if ( isFASTA(DB_1)) {
				
		/* No more threads than sequences in DB file */
		{
			const size_t MaxJobs = FASTA_DB1.SequenceCount * FASTA_DB2.SequenceCount; 
			if (nCPUs > MaxJobs) nCPUs = MaxJobs;
		}

		for (size_t iRef=0UL; iRef<FASTA_DB1.SequenceCount; iRef++) {
			/* Get space element */
		reclaim:
			pthread_mutex_lock(&pool->donequeue_p.rwmutex);
			MapJob_t * const task = (MapJob_t*) jobqueue_pull(&pool->donequeue_p);
			pthread_mutex_unlock(&pool->donequeue_p.rwmutex);
			if (task == NULL) {
                sleep(1);
                goto reclaim;
            }
			task->id = iRef;
			
			/* Adding new task to jobqueue */
			pthread_mutex_lock(&pool->jobqueue_p.rwmutex);
			jobqueue_push(&pool->jobqueue_p, (job_t*) task);
			pthread_mutex_unlock(&pool->jobqueue_p.rwmutex);
            
            if (OutputVerbose) fprintf(stdout,"Master sent job %lu/%lu to queue\n", iRef, FASTA_DB1.SequenceCount);
		}
		
		goto Input_Found;
	}

	/* ------------------------------------  UNDEFINED ------------------------------------------ */
	
	fputs("Database files %s could not be determined\n.", stderr);
	res = 1;
	
	Input_Found: ;
	
	////////////////////////////////////////////////////////////////////////////////////////////////
  // WAIT FOR QUEUE TO BE DONE
  ////////////////////////////////////////////////////////////////////////////////////////////////
    if (OutputVerbose) 
	thpool_wait(pool);

	////////////////////////////////////////////////////////////////////////////////////////////////
  // CLEANLY CLOSE
  ////////////////////////////////////////////////////////////////////////////////////////////////
	
	CloseFASTAStructure(&FASTA_DB1);
	
  /* Free Memory */
	// References 
	
bail:;
	destroyThreadPool(pool);
	
#ifdef PRF_USE_AFFINITY
  if (Thread_masks[0]) free(Thread_masks[0]);
  if (Thread_masks[1]) free(Thread_masks[1]);
#endif
	
	freeSystemInfo(&System);

  return res;
}
