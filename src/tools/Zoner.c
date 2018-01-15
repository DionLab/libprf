#define _GNU_SOURCE
#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/mman.h> 
#include <pthread.h>
#include <stdint.h>
#include <getopt.h>
#ifdef PRF_USE_AFFINITY
# include <unistd.h>
# include <sched.h>
#endif
#include <nmmintrin.h>

#define HEADER "%----------------------------------------------------------------------------%\n"\
	       "|                              Zoner v" PRF_VERSION "                          |\n"\
	       "%----------------------------------------------------------------------------%\n"\
	       "| Built on " __DATE__ " at " __TIME__ ".                                          |\n"
#include "pfVersion.h"
#include "pfMap.h"
#include "pfSequence.h"
#include "pfInput.h"
#include "pfOutput.h"
#include "pfCompute.h"
#include "pfDispatchExt.h"

#include "threadpool.h"

#ifdef PRF_USE_AFFINITY
#define __USE_AFFINITY__
#endif
#include "system.h"

typedef struct COMMON {
	const struct Map * map;
	const char * Genome;
	size_t GenomeLength;
	pthread_mutex_t PrintLock;
} common_t;

union universal_job_s {
	fasta_job_t FASTA;
#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
	pb_job_t PacBio;
#endif
} universal_job_t;
 
static SystemInfo System;
#ifdef PRF_USE_AFFINITY
cpu_set_t * Thread_masks[2] = {0,0};						/* Define variables to hold thread affinity mask */
unsigned int Thread_count[2] = {0,0};
pthread_attr_t * restrict threads_attr = NULL;
#endif

#define TEST_IF_AVAILABLE(TYPE) if (OutputType.Type != TYPE) {\
	fprintf(stderr,"Invalid option for this type");\
	exit(1);\
}

static const char opt_to_test[] = "d:g:G:ac:H:I:i:t:hs0:o:Vw!:X:FD:em,Q"
#ifdef PRF_USE_AFFINITY
	"012:3"
#endif
#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
	"jJK:yp:L:qYTM"
#endif
;

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
	/* Profile */
	{"cutoff",								required_argument,	0,	'c'},
	{"optimal", 							no_argument,				0,	'a'},
	{"bestof",								no_argument,				0,	','},
	{"json-prf",							required_argument,	0,	'G'},
	{"dump-json",							required_argument,	0,	'd'},
	{"genome",								required_argument,	0,	'g'},
	/* Sequence */
	{"with-revcomp",					no_argument,				0,	'e'},
	/* Database indexing options */
	{"create-index-database",	required_argument,	0,	'I'},
	{"use-index-database",		required_argument,	0,	'i'},
	{"stream-fasta",					no_argument,				0,	'.'},
#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
	/* PacBio */
	{"subreads",							no_argument,				0,	'j'},
	{"consensus",							no_argument,				0,	'J'},
	{"polymerase",						no_argument,				0,	'x'},
	{"best-subread",					no_argument,				0,	'M'},
	{"invalid-zmw",						no_argument,				0,	'O'},
	{"zmw-number",						required_argument,	0,	'K'},
	{"not-only-sequencing",		no_argument,				0,	'y'},
	{"min-filter-score",			required_argument,	0,	'p'},
	{"max-filter-score",			required_argument,	0,	'L'},
	{"keep-invalid-zmw",			no_argument,				0,	'q'},
	{"discard-filter",				no_argument,				0,	'Y'},
	{"test-filter",						no_argument,				0,	'T'},
#endif
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

int OutputVerbose = 0;
static const struct timespec polling_interval = { .tv_sec = 0, .tv_nsec = 400};

static void __attribute__((noreturn)) Usage(FILE * stream)
{
  fputs(
	" Zoner --genome file [options] [Database, - ]\n"
	" Options:\n"
	"  Alignment:\n"
	"    --json-prf file                : json containing scores\n"
	"    --dump-json file               : dump default values into json file\n\n"
	"  Sequence\n"
	"    --with-revcomp                 : test also reverse complement\n"
	"    --bestof                       : when used with --with-revcomp and --optimal,\n"
	"                                     output only the best of both orientation\n\n"
	"  Database\n"
	"   FASTA\n"
	"    Fasta format allows to be piped with the '-' symbol as database file name\n\n"
#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM) 
#if !defined(PRF_INPUT_HDF5)
	"   PacBio BAM\n"
#elif !defined(PRF_INPUT_PBBAM)
	"   PacBio HDF5\n"
#else
	"   PacBio HDF5 or BAM\n"
#endif
	"     --subreads                    : work on subreads (default)\n"
	"     --consensus                   : work on consensus\n"
	"     --polymerase                  : work on entire polymerase data\n"
	"     --invalid-zmw                 : work on invalid zmw only\n"
	"     --best-subread                : keep only best subread of zmw\n"
	"     --zmw-number i                : limit to zmw i\n"
	"     --not-only-sequencing         : compute on any zmw type\n"
	"     --min-filter-score i          : minimal High Quality region score i to keep data\n"
	"     --max-filter-score i          : maximum High Quality region score i to keep data\n"
	"     --keep-invalid-zmw            : keep zmw removed by PacBio\n"
	"     --discard-filter              : does not filter or trim according to HQ region\n"
	"     --test-filter                 : output regions based upon the above filtering\n\n"
#endif
	" Optimizations\n"
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

static void * FastaThread(threadarg_t* const restrict _Data) 
{
		int res = -1;
		/*************************************************************************/
		/*                          GET COMMON DATA                              */
		/*************************************************************************/
		common_t * const common = _Data->common;
		
		/* For profile search we need to have */
		const struct Map * const restrict map = common->map;
		const char * const restrict Genome = common->Genome;
		const size_t GenomeLength = common->GenomeLength;
		
		/* For output we need */
// 		pthread_mutex_t * const restrict PrintLock = &(common->PrintLock);
		
		/*************************************************************************/
    /*                          ALLOCATE MEMORY                              */
    /*************************************************************************/
		size_t MemorySize = 2048UL*1024UL*1024UL;
		void * WORK = mmap(NULL, MemorySize, PROT_READ | PROT_WRITE,
                           MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_HUGETLB,
                           -1, 0);
		if (WORK == MAP_FAILED) {
			fprintf(stderr, "%s: Unable to allocate %lu Mbytes Huge TLB, going for standard!\n", __FUNCTION__, MemorySize >> 20);
			WORK = mmap(NULL, MemorySize, PROT_READ | PROT_WRITE,
						MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE,
                        -1, 0);
			if (WORK == MAP_FAILED) {
				fprintf(stderr, "%s: Unable to allocate %lu Mbytes standard memory, quitting!\n", __FUNCTION__, MemorySize >> 20);
				goto FIN;
			}
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
		fasta_job_t * task = NULL;
		
		while(thpool->threads_keepalive) {
			bsem_wait(thpool->jobqueue_p.has_items);
			
			if (thpool->threads_keepalive) {
				pthread_mutex_lock(&thpool->thcount_lock);
				thpool->num_threads_working++;
				pthread_mutex_unlock(&thpool->thcount_lock);
				
				/* Read job from queue and execute it */
				pthread_mutex_lock(&thpool->jobqueue_p.rwmutex);
				task = (fasta_job_t*) jobqueue_pull(&thpool->jobqueue_p);
				pthread_mutex_unlock(&thpool->jobqueue_p.rwmutex);
				
				/***************************************************************/
				/*                     START THE TASK                          */
				/***************************************************************/
				if (task) {
					/* Transform into real Fasta sequence and split header and data */
					Sequence_t * const SeqData = &(task->Fasta_Sequence);
					if (SeqData->ProfileData.Length > 0UL) {
						size_t SeqLength = 0UL;
						{
							unsigned char * restrict ptr = SeqData->ProfileData.ProfileIndex;
							const unsigned char * const limit = ptr + SeqData->ProfileData.Length; 
							while ( ptr < limit) {
								register unsigned char c = *ptr++;
								c = (c >= 'a') ? c - ('a' - 'A') : c;
								if (c >= 'A' && c <= 'Z') 
									SeqData->ProfileData.ProfileIndex[SeqLength++] = c;
							}
							SeqData->ProfileData.Length = SeqLength;
						}
						
						/* Check memory allocated size */
						while (2*sizeof(__m128i)*SeqLength >= MemorySize) {
							munmap(WORK, MemorySize);
							MemorySize *= 2UL;
							WORK = mmap(NULL, MemorySize, PROT_READ | PROT_WRITE,
										MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_HUGETLB,
										-1, 0);
							if (WORK == MAP_FAILED) {
								fprintf(stderr, "%s: Unable to allocate %lu Mbytes Huge TLB, going for standard!\n", __FUNCTION__, MemorySize >> 20);
								WORK = mmap(NULL, MemorySize, PROT_READ | PROT_WRITE,
											MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE,
											-1, 0);
								if (WORK == MAP_FAILED) goto FIN;
							}
						}

						/////////////////////////////////////////////////////////////////////////////////////////////////////////
						// WORK
						Zone_t Alignment;
						struct timeval t0, t1;
						gettimeofday(&t0, NULL);
						GetMapping(map, WORK, Genome, SeqData->ProfileData.ProfileIndex, &Alignment, SeqLength, GenomeLength);
						gettimeofday(&t1, NULL);
						const double elapse_time = 1.0E9*(double) (t1.tv_sec - t0.tv_sec) + 1000.0 * (double) (t1.tv_usec - t0.tv_usec);
						const double speed = (double) (SeqLength*GenomeLength)/elapse_time;
						fprintf(stdout,"%s: Score = %i, range = [%u,%u], speed = %5.2lf gcups\n", SeqData->Data.Header, Alignment.Score, Alignment.Begin, Alignment.End, speed); 
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

		res = 0;
		munmap(WORK, MemorySize);

    FIN:
    pthread_mutex_lock(&thpool->thcount_lock);
    thpool->num_threads_alive--;
    pthread_mutex_unlock(&thpool->thcount_lock);
		
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

#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
static void * PacBioThread(threadarg_t* const restrict _Data)
{
		int res = -1;
		/*************************************************************************/
		/*                          GET COMMON DATA                              */
		/*************************************************************************/
		common_t * const common = _Data->common;
		
		/* For profile search we need to have */
		const struct Map * const restrict map = common->map;
		const char * const restrict Genome = common->Genome;
		const size_t GenomeLength = common->GenomeLength;
		
		/* For output we need */
// 		pthread_mutex_t * const restrict PrintLock = &(common->PrintLock);
		
		/*************************************************************************/
    /*                          ALLOCATE MEMORY                              */
    /*************************************************************************/
		size_t MemorySize = 2048UL*1024UL*1024UL;
		void * WORK = mmap(NULL, MemorySize, PROT_READ | PROT_WRITE,
											 MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_HUGETLB,
	                     -1, 0);
		if (WORK == MAP_FAILED) {
			fprintf(stderr, "%s: Unable to allocate %lu Mbytes Huge TLB, going for standard!\n", __FUNCTION__, MemorySize >> 20);
			WORK = mmap(NULL, MemorySize, PROT_READ | PROT_WRITE,
									MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE,
	                -1, 0);
			if (WORK == MAP_FAILED) goto FIN;
		}
		
				
		/*************************************************************************/
		/*                          CONFIGURE WORKER                             */
		/*************************************************************************/
		/* Assure all threads have been created before starting serving */
    threadpool_t* const restrict thpool = _Data->thpool;
       
    /* Mark thread as alive (initialized) */
    pthread_mutex_lock(&thpool->thcount_lock);
    thpool->num_threads_alive += 1;
    pthread_mutex_unlock(&thpool->thcount_lock);
    
    // Initialize some other data
    pb_job_t * task = NULL;
			
		while(thpool->threads_keepalive) {
			bsem_wait(thpool->jobqueue_p.has_items);
			
			if (thpool->threads_keepalive){
			    pthread_mutex_lock(&thpool->thcount_lock);
			    thpool->num_threads_working++;
			    pthread_mutex_unlock(&thpool->thcount_lock);
			    
			    /* Read job from queue and execute it */
			    pthread_mutex_lock(&thpool->jobqueue_p.rwmutex);
			    task = (pb_job_t*) jobqueue_pull(&thpool->jobqueue_p);
			    pthread_mutex_unlock(&thpool->jobqueue_p.rwmutex);
			    
			    /***************************************************************/
			    /*                     START THE TASK                          */
			    /***************************************************************/
			    if (task) {
						const HoleData_t * const HReg = &(task->HReg);
						const unsigned int HoleNumber = HReg->HoleNumber;
				    const unsigned int nRegions   = HReg->nRegions;
// 						printf("Thread %lu working on hole %u, %u regions found\n", ((ThreadPool_Arg_t*)_Data)->ID, HoleNumber, nRegions);
						restrictReadstoHQ(HReg);
						
						for (int i=0; i<nRegions; i++) {
							if (HReg->Regions[i].type == Insert) {
								const int lSeqLength = (int) HReg->Regions[i].stop - (int) HReg->Regions[i].start;
								if (lSeqLength <= 1) continue;
								
								/* Check memory allocated size */
								while (2*sizeof(__m128i)*lSeqLength >= MemorySize) {
									munmap(WORK, MemorySize);
									MemorySize *= 2UL;
									WORK = mmap(NULL, MemorySize, PROT_READ | PROT_WRITE,
															MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_HUGETLB,
															-1, 0);
									if (WORK == MAP_FAILED) {
										fprintf(stderr, "%s: Unable to allocate %lu Mbytes Huge TLB, going for standard!\n", __FUNCTION__, MemorySize >> 20);
										WORK = mmap(NULL, MemorySize, PROT_READ | PROT_WRITE,
																MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE,
																-1, 0);
										if (WORK == MAP_FAILED) goto FIN;
									}
								}
								
								unsigned char * const restrict SequenceText = &(task->Reads.Stream[HReg->Regions[i].start]);
								/////////////////////////////////////////////////////////////////////////////////////////////////////////
								// WORK
								Zone_t Alignment;
								GetMapping(map, WORK, Genome, SequenceText, &Alignment, GenomeLength, lSeqLength);
								fprintf(stdout,"ZMW %u-%i/%u: Score = %i, range = [%u,%u]\n", HoleNumber, 1+i, nRegions, Alignment.Score, Alignment.Begin, Alignment.End); 
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
     
    munmap(WORK, MemorySize);
    res = 0; 
		
    FIN:
    pthread_mutex_lock(&thpool->thcount_lock);
    thpool->num_threads_alive--;
    pthread_mutex_unlock(&thpool->thcount_lock);
		    
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
#endif

int main (int argc, char *argv[])
{
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // LOCAL STRUCTURES
  ////////////////////////////////////////////////////////////////////////////////////////////////
	struct Map map;                           /* Mapping */
  const char * restrict GenomeFile = NULL;  /* Genome */
  const char * restrict jsonPrf = NULL;			/* JSON small prf */
  FASTAStructure FASTA = { 0 };							/* Sequence Database File */
#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
  PacBio_t * restrict PBS = NULL;						/* PacBio HDF5 datafile */
#endif
  struct timeval _t0, _t1;									/* Timing structures */

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // LOCAL DATA
  ////////////////////////////////////////////////////////////////////////////////////////////////

  size_t nCPUs=0;										/* number of threads */
  int res;
  int Cutoff = -1;									/* Default Coverage cutoff from command line, if not zero then enforces that value */
  char * DB = NULL;									/* FASTA sequence file */
#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
  int isPacBio = 0;									/* Are we crunching HDF5 file from PacBio ?*/
  PacBioDispatchOptions_t PacBioDispatchOptions = PB_DEFAULT_OPTIONS;
	unsigned int PacBioZmWNumber;
#endif
  size_t * shares = 0;
  size_t MaxProfileLength = 0;			/* Maxiumum profile length when several are given */
  struct ThreadArrayData *threads_array_arg = NULL;
	_Bool FastaStream = false;

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
  if (!(System.Extensions & MM_SSE41)) {
      fputs("Sorry only compute methods based upon SSE 4.1 are available yet\n", stderr);
      exit(1);
  }
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // OPTIONS
  ////////////////////////////////////////////////////////////////////////////////////////////////
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
  }
  
  optind = 1;
  while (1) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    const int c = getopt_long (argc, argv, opt_to_test, long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) break;
    switch (c) {
			/****************************************** Profile  *****************************************/
			case 'c':
				Cutoff = atoi(optarg);
				break;
			case 'g':
				GenomeFile = optarg;
				break;
			case 'G':
				jsonPrf = optarg;
				break;
			case 'd':
				if (dumpDefaultMap(optarg) != 0) {
					fprintf(stderr, "Unable to create %s\n", optarg);
					exit(1);
				}
				else
					exit(0);
				break;
			/****************************************** Sequence *****************************************/
#ifdef PRF_USE_AFFINITY
      case '0':
        noAffinity = true;
        break;
      case '3':
				noSharedCore = true;
				break;
      case '2':
				GivenAffinityFile = true;
				AffinityMaskFileName = optarg;
				break;
#endif
      case 'V':
				OutputVerbose = true;
				break;
      case 't':
				nCPUs = (size_t) atoi(optarg);
				break;
      case 'I':
				FASTA.Options |= DoFastaIndexExport;
				FASTA.indexFileName = optarg;
				break;
      case 'i':
				FASTA.Options |= DoFastaIndexImport;
				FASTA.indexFileName = optarg;
				break;
			case '.':
				FastaStream = true;
				break;
				/************************************** TEXT specific options ******************************/
#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
			case 'j':
				PacBioDispatchOptions.Selection |= PB_DISPATCH_SUBREADS;
				break;
			case 'J':
				PacBioDispatchOptions.Selection |= PB_DISPATCH_CONSENSUS;
				break;
			case 'x':
				PacBioDispatchOptions.Selection |= PB_DISPATCH_ZMW;
				break;
			case 'O':
				PacBioDispatchOptions.Selection |= PB_DISPATCH_INVALID;
				break;
			case 'K':
				{
					const int itmp = atoi(optarg);
					if (itmp < 0 ) {
						fprintf(stderr, "PacBio hole number must be positive or null (%i)\n",itmp);
						exit(1);
					}
					else {
						PacBioZmWNumber = (unsigned int) itmp;
						PacBioDispatchOptions.ZMW = &PacBioZmWNumber;
						PacBioDispatchOptions.nZMW = 1;
					}
				}
				break;
			case 'p':
				PacBioDispatchOptions.minReadAccuracy = atoi(optarg);
				PacBioDispatchOptions.Selection |= PB_HAS_FILTER_SCORE;
				break;
			case 'L':
				PacBioDispatchOptions.maxReadAccuracy = atoi(optarg);
				PacBioDispatchOptions.Selection |= PB_HAS_FILTER_SCORE;
				break;
			case 'q':
				PacBioDispatchOptions.Selection |= PB_KEEP_INVALID;
				break;
			case 'Y':
				PacBioDispatchOptions.Selection |= PB_DISCARD_FILTER;
				break;
			case 'y':
				PacBioDispatchOptions.Selection |= PB_NOT_ONLY_SEQUENCING;
				break;
			case 'T':
				PacBioDispatchOptions.Selection |= PB_TEST_OUTPUT;
				break;
			case 'M':
				PacBioDispatchOptions.Selection |= PB_BEST_SUBREAD_OF_ZMW;
				break;
#endif
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
    DB = argv[optind];
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // COMMAND LINE COHERENCY
  ////////////////////////////////////////////////////////////////////////////////////////////////
  {
		_Bool AllOK = true;
		
		if (GenomeFile == NULL) {
			fputs("Expecting both --genome and --json-prf\n", stderr);
			AllOK = false;
		}

		if (!DB) {
			fputs("Expected arguments after options\n", stderr);
			AllOK = false;
		}
  
		if (!AllOK) Usage(stderr);
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////
  // COMMON VARIABLE SETTINGS
  ////////////////////////////////////////////////////////////////////////////////////////////////
  common_t common = {
		.PrintLock = PTHREAD_MUTEX_INITIALIZER,
	};  

	
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // LOAD GENOME
  ////////////////////////////////////////////////////////////////////////////////////////////////
	int genFD = open(GenomeFile, O_RDONLY);
	if (genFD < 0) {
		fprintf(stderr, "Unable to oipen Genome file %s\n", GenomeFile);
		goto bail;
	}
	
	struct stat st;
	if (fstat(genFD, &st) < 0) {
		perror("fstat");
		goto bail;
	}
  
  common.GenomeLength = st.st_size;
	common.Genome = mmap(NULL, st.st_size,  PROT_READ, MAP_PRIVATE, genFD, 0);
	if (common.Genome == MAP_FAILED) {
		fprintf(stderr, "Unable to map genome file %s into RAM\n", GenomeFile);
		goto bail;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////
  // LOAD MAP
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (jsonPrf == NULL) {
		common.map = getDefaultMap();
	}
	else {
		if (loadMap(jsonPrf, &map) != 0) {
			fprintf(stderr, "Error loading json score file %s\n", jsonPrf);
			goto bail;
		}
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
#ifdef USE_PACBIO 
		if (!(PacBioDispatchOptions.nZMW == 1U))
#endif
			fprintf(stderr, "Job dispatched over %lu cores.\n", nCPUs);
	}
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // ANALYSIS BASED ON DATABASE
  ////////////////////////////////////////////////////////////////////////////////////////////////
	/* -------------------------------- FASTA SEQUENCES STREAM ---------------------------------- */
#ifdef PRF_INPUT_FASTA
	if (( DB[0] == '-' && DB[1] == '\0')) {
		 threadpool_t * const restrict thpool = createThreadPool(FastaThread,
#ifdef PRF_USE_AFFINITY
                                                                 NULL,
#endif
																sizeof(universal_job_t),
																(void*) &common,
																nCPUs,0);
		if ( !thpool ) {
			res = -3;
			goto bail;
		}
		
		res = dispatchStreamFASTAExt(stdin, thpool);
		destroyThreadPool(thpool);
		goto Input_Found;
	}
	/* ---------------------------------  FASTA SEQUENCES --------------------------------------- */
	else if (isFASTA(DB) ) {
		FILE *fd = fopen(DB, "r");
		if (fd == NULL) {
			fprintf(stderr, "Error opening file %s\n", DB);
			goto bail;
		}
		threadpool_t * const restrict thpool = createThreadPool(FastaThread,
#ifdef PRF_USE_AFFINITY
																NULL,
#endif
																sizeof(universal_job_t),
																(void*) &common,
																nCPUs,0);
		if ( !thpool ) {
			res = -3;
			goto bail;
		}
		res = dispatchStreamFASTAExt(fd, thpool);
		fclose(fd);
		destroyThreadPool(thpool);
		goto Input_Found;
	}
#endif

	/* -------------------------------  PAC BIO HDF5 SEQUENCES ---------------------------------- */
#ifdef PRF_INPUT_HDF5
	if ( isPacBioH5(DB) ) {
		PBS = PacBioOpen(DB);
		if (PBS == NULL) {
			fprintf(stderr, "Error opening Pac Bio file %s\n", DB);
			goto bail;
		}
		if (!(PBS->Content & BASECALLING)) {
			fputs("Error Pac Bio file is not containing Basecalling data\n", stderr);
			goto bail;
		}
		
		gettimeofday(&_t0,0);
		const int tres = (PacBioDispatchOptions.Selection & PB_DISPATCH_CONSENSUS) ? IndexConsensus(PBS) : IndexRegions(PBS);
		gettimeofday(&_t1,0);
		if (tres != SUCCESS) {
			fprintf(stderr,"PacBio error indexing regions, code %i\n", tres);
			goto bail;
		}
		
		if (OutputVerbose) {
			fprintf(stderr,
							"Base file name : %s\n"
							"Directory      : %s\n"
							"Parts          : %u\n",
							PBS->BaseFileName, PBS->Directory, PBS->nParts);
			for (unsigned int i=0; i<PBS->nParts; i++) {
				fprintf(stderr, "               : %s\n", PBS->PartsFileName[i]);
			}
			
			fprintf(stderr, "ZMW number    : %u\n", PBS->nHoles);
			const double T = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr,"Indexing holes took %lf seconds.\n", T);
		}
		
		if (PopulateHoleStatus(PBS) != SUCCESS) {
			fprintf(stderr, "ZMW status not populated correctly...\n");
			res = 1;
			goto bail;
		}
		
		if (PopulateHoleCoordinates(PBS) != SUCCESS) {
			fprintf(stderr, "ZMW coordinates not populated correctly...\n");
			res = 1;
			goto bail;
		}
		
		if (PacBioDispatchOptions.ZMW == NULL) {
			if (!(PacBioDispatchOptions.Selection & PB_NOT_ONLY_SEQUENCING)) {
				PacBioDispatchOptions.nZMW = getSequencingHoles(PBS, &(PacBioDispatchOptions.ZMW));
				if (PacBioDispatchOptions.nZMW > 0) {
					if (OutputVerbose) fprintf(stderr, "Found %i sequencing holes\n", PacBioDispatchOptions.nZMW);
				}
				else {
					fprintf(stderr, "getSequencingHoles returned %i\n", PacBioDispatchOptions.nZMW);
					goto bail;
				}
			}
			else {
				PacBioDispatchOptions.nZMW = PBS->nHoles;
				PacBioDispatchOptions.ZMW = (unsigned int*) malloc(PacBioDispatchOptions.nZMW*sizeof(unsigned int));
				if (!PacBioDispatchOptions.ZMW) {
					fprintf(stderr, "Unable to allocate memory for list of sequencing holes\n");
					goto bail;
				}
				for(int i=0; i<PacBioDispatchOptions.nZMW; i++) PacBioDispatchOptions.ZMW[i] = i;
			}
		}
		
		if (nCPUs > PacBioDispatchOptions.nZMW) nCPUs = PacBioDispatchOptions.nZMW;
		threadpool_t * const restrict thpool = createThreadPool(PacBioThread,
#ifdef PRF_USE_AFFINITY
																														NULL,
#endif
																														sizeof(universal_job_t),
																														(void*) &common,
																														nCPUs,0);
		if ( !thpool ) {
			res = -3;
			goto bail;
		}
		res = dispatchPacBioExt(PBS, thpool, &PacBioDispatchOptions);
		destroyThreadPool(thpool);
		if ((PacBioDispatchOptions.nZMW != 1) && PacBioDispatchOptions.ZMW) free(PacBioDispatchOptions.ZMW);
		if (PacBioClose(PBS) != SUCCESS) {
			fprintf(stderr, "Error closing Pac Bio...\n");
		}
		
		goto Input_Found;
	}
#endif
	/* -------------------------------  PAC BIO BAM SEQUENCES ----------------------------------- */
#ifdef PRF_INPUT_PBBAM
	if (isPacBioBAM(DB)) {
		PacBioBAM_t * const PBBAM = OpenPacBioBAM(DB);
		if (PBBAM) {
			threadpool_t * const restrict thpool = createThreadPool(PacBioThread,
#ifdef PRF_USE_AFFINITY
																	NULL,
#endif
																	sizeof(universal_job_t),
																	(void*) &common,
																	nCPUs,0);
			if ( !thpool ) {
				res = -3;
				goto bail;
			}
			res = dispatchPacBioBAMExt(PBBAM, thpool, &PacBioDispatchOptions);
			destroyThreadPool(thpool);
			ClosePacBioBAM(PBBAM);
		}
		else {
			fprintf(stderr, "Error in OpenPacBioBAM\n");
			res = 1;
		}
		goto Input_Found;
	}
#endif
	/* ------------------------------------  UNDEFINED ------------------------------------------ */
	
	fprintf(stderr, "Database file %s could not be determined\n.", DB);
	res = 1;
	Input_Found:;

	////////////////////////////////////////////////////////////////////////////////////////////////
  // CLEANLY CLOSE
  ////////////////////////////////////////////////////////////////////////////////////////////////
  
  /* Free Memory */
bail:;

	if (common.Genome != MAP_FAILED) munmap((void*) common.Genome, common.GenomeLength);

#ifdef PRF_USE_AFFINITY
  if (Thread_masks[0]) free(Thread_masks[0]);
  if (Thread_masks[1]) free(Thread_masks[1]);
#endif
	
	freeSystemInfo(&System);
	pthread_mutex_destroy(&(common.PrintLock));

  return res;
}

