#define _GNU_SOURCE
#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include <assert.h>
#include "pfCompute.h"
#include "pfOutput.h"
#include "pfDispatch.h"

int PacBioHQThreshold = 0;
static const struct timespec polling_interval = { .tv_sec =0, .tv_nsec = 400L};

static int thpool_add_subread_work(threadpool_t * const restrict thpool_p, const PacBio_t * const restrict PB, 
                                   const PacBioDispatchOptions_t * const restrict Options)
{
	int res = ERROR;
	if (!PB) return -1;
	if (!PB->BaseCalling.Regions) return -2;

	pb_job_t * restrict newJob;
	const size_t N = Options->nZMW;
	
	
	/* Get a job slot from the pool donequeue */
	again:;
	pthread_mutex_lock(&thpool_p->donequeue_p.rwmutex);
	newJob = (pb_job_t*) jobqueue_pull(&thpool_p->donequeue_p);
	pthread_mutex_unlock(&thpool_p->donequeue_p.rwmutex);
	
	if (!newJob) {
		if (thpool_p->num_threads_alive > 0) {
			/* Wait a bit and try again */
			nanosleep(&polling_interval, NULL);
			goto again;
		}
		else {
			return -1;
		}
	}
	
	for (size_t iJob = 0; iJob < N; iJob++) {
		const unsigned int HoleNumber = Options->ZMW[iJob];
		RegionIndex_t * const restrict Region = &(PB->BaseCalling.Regions[HoleNumber]);
		
		/* Get the number of regions in hole number */
		const int RegionSize = Region->stop - Region->start + 1;

		/* Verify available memory */ 
		if (newJob->HReg.nAllocatedRegions < RegionSize) {
			newJob->HReg.Regions = (HoleRegion_t * restrict) realloc(newJob->HReg.Regions, RegionSize*sizeof(HoleRegion_t));
			if (newJob->HReg.Regions == NULL) { res= -4; goto bail; }
			newJob->HReg.nAllocatedRegions = RegionSize;
		}
		memset(newJob->HReg.Regions, 0, RegionSize*sizeof(HoleRegion_t));
		
		hsize_t dims[2];
		const size_t FileNum = PB->HoleFileIndex[HoleNumber];

		hid_t dataset = H5Dopen2(PB->Parts[FileNum], "PulseData/Regions", H5P_DEFAULT); if (dataset < 0) goto bail;
		hid_t dataspace = H5Dget_space(dataset); if (dataspace < 0) goto bail;
		H5Sget_simple_extent_dims(dataspace, dims, NULL);
		
		hsize_t outdims[2] = { RegionSize, 4 };
		hid_t mspace_id = H5Screate_simple(2, outdims, NULL);
		
		hsize_t hyperslab_start[2] = { Region->start, 1};
		hsize_t hyperslab_count[2] = { RegionSize, 4};
		
		H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
		if (H5Dread(dataset, H5T_NATIVE_INT, mspace_id, dataspace, H5P_DEFAULT, newJob->HReg.Regions) < 0) goto bail;
		
		
		H5Sclose(mspace_id);
		H5Sclose(dataspace);
		H5Dclose(dataset);
		
		/* Get the sequence range */
		int max = 0;
		{
			_Bool HasValidHQ = true;
			int HQindex = 0;
			HoleRegion_t * const HR = newJob->HReg.Regions;
			/* WARNING: That assumes inserts are always prior to HQregion !! */
			for (int r=0; r<RegionSize; r++) {
				if (HR[r].type == Insert) {
					const int lstop = HR[r].stop;
					max = (max > lstop) ? max : lstop;
				}
				else if (HR[r].type == HQRegion) {
					if (HR[r].stop == 0) {
						HasValidHQ = false;
					}
					HQindex = r;
				}
			}
			
			/* 1. Do we keep regular */
			if (HasValidHQ && (Options->Selection & PB_DISPATCH_INVALID)) continue;
			
			/* 2. Do we keep invalid if flag s such ?*/
			else if (!HasValidHQ && !(Options->Selection & PB_DISPATCH_INVALID) && !(Options->Selection & PB_KEEP_INVALID)) continue; 

			/* 3. Filter */
			if (HasValidHQ && !(Options->Selection & PB_DISCARD_FILTER) && (Options->Selection & PB_HAS_FILTER_SCORE)) {
					if (HR[HQindex].quality < Options->minReadAccuracy || HR[HQindex].quality > Options->maxReadAccuracy ) continue;
			}

			/* 4. Correct regions in case of invalid */
			if (!HasValidHQ) {
					HR[HQindex].stop = max;
			}
            
			
			if (((Options->Selection & PB_KEEP_INVALID) && !HasValidHQ) || (Options->Selection & PB_DISCARD_FILTER)) {
				HR[HQindex].start = 0;
				HR[HQindex].stop = max;
			}
		}
		
		const size_t msize = max;// - min + 2;
		if (newJob->Reads.Size < msize) {
			newJob->Reads.Stream = (unsigned char*) realloc(newJob->Reads.Stream, msize*sizeof(unsigned char));
			if (newJob->Reads.Stream == NULL) { res = -5; goto bail; }
			newJob->Reads.Size = msize;
		}
		memset(newJob->Reads.Stream, 0, msize*sizeof(char));
		/* Read sequence string */
		
		dataset = H5Dopen2(PB->Parts[FileNum], "PulseData/BaseCalls/Basecall", H5P_DEFAULT); if (dataset < 0)  goto bail;
		dataspace = H5Dget_space(dataset); if (dataspace < 0)  goto bail;
		H5Sget_simple_extent_dims(dataspace, dims, NULL);
		
		outdims[0] = msize;
		outdims[1] =  0;
		mspace_id = H5Screate_simple(1, outdims, NULL);
		
		hyperslab_start[0] = Region->BaseCallsOffset;// + min;
		hyperslab_start[1] = 1;
		hyperslab_count[0] = msize;
		hyperslab_count[1] = 1;
		
		H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
		if (H5Dread(dataset, H5T_NATIVE_UCHAR, mspace_id, dataspace, H5P_DEFAULT, newJob->Reads.Stream) < 0)  goto bail;
		
		H5Sclose(mspace_id);
		H5Sclose(dataspace);
		H5Dclose(dataset);
		
		newJob->HReg.nRegions       = (unsigned int) RegionSize;
		newJob->HReg.HoleNumber     = HoleNumber;
		newJob->HReg.Coordinates[0] = PB->Coordinates[HoleNumber][0];
		newJob->HReg.Coordinates[1] = PB->Coordinates[HoleNumber][1];
				
		pthread_mutex_lock(&thpool_p->jobqueue_p.rwmutex);
		jobqueue_push(&thpool_p->jobqueue_p, (job_t*) newJob);
		pthread_mutex_unlock(&thpool_p->jobqueue_p.rwmutex);
		
		/* Get a job slot from the pool donequeue */
		again2:
		pthread_mutex_lock(&thpool_p->donequeue_p.rwmutex);
		newJob = (pb_job_t*) jobqueue_pull(&thpool_p->donequeue_p);
		pthread_mutex_unlock(&thpool_p->donequeue_p.rwmutex);
		
		if (!newJob) {
			if (thpool_p->num_threads_alive > 0) {
				/* Wait a bit and try again */
				nanosleep(&polling_interval, NULL);
				goto again2;
			}
			else {
				res = -1;
				goto bail;
			}
		}
	}
	
	res = SUCCESS;
	
	bail:;
	if (newJob) {
		pthread_mutex_lock(&thpool_p->donequeue_p.rwmutex);
		jobqueue_push(&thpool_p->donequeue_p, (job_t*) newJob);
		pthread_mutex_unlock(&thpool_p->donequeue_p.rwmutex);
	}
	
	return res;
}

static int thpool_add_read_work(threadpool_t * const restrict thpool_p, const PacBio_t * const restrict PB, 
                                const PacBioDispatchOptions_t * const restrict Options)
{
	int res = ERROR;
	if (!PB) return -1;
	if (!PB->BaseCalling.Regions) return -2;

	pb_job_t * restrict newJob;
	const size_t N = Options->nZMW;
	
	/* Get a job slot from the pool donequeue */
	again:
	pthread_mutex_lock(&thpool_p->donequeue_p.rwmutex);
	newJob = (pb_job_t*) jobqueue_pull(&thpool_p->donequeue_p);
	pthread_mutex_unlock(&thpool_p->donequeue_p.rwmutex);
	
	if (!newJob) {
		if (thpool_p->num_threads_alive > 0) {
			/* Wait a bit and try again */
			nanosleep(&polling_interval, NULL);
			goto again;
		}
		else {
			res = -1;
			goto bail;
		}
	}
	
	for (size_t iJob = 0; iJob < N; iJob++) {
		const unsigned int HoleNumber = Options->ZMW[iJob];
		RegionIndex_t * const restrict Region = &(PB->BaseCalling.Regions[HoleNumber]);
		
		/* Get the number of regions in hole number */
		const int RegionSize = Region->stop - Region->start + 1;
		
		if (newJob->HReg.nAllocatedRegions < RegionSize) {
			newJob->HReg.Regions = (HoleRegion_t * restrict) realloc(newJob->HReg.Regions, RegionSize*sizeof(HoleRegion_t));
			if (newJob->HReg.Regions == NULL) { res = -4; goto bail; }
			newJob->HReg.nAllocatedRegions = RegionSize;
		}
		memset(newJob->HReg.Regions, 0, RegionSize*sizeof(HoleRegion_t));
		
		hsize_t dims[2];
		const size_t FileNum = PB->HoleFileIndex[HoleNumber];

		hid_t dataset = H5Dopen2(PB->Parts[FileNum], "PulseData/Regions", H5P_DEFAULT); if (dataset < 0) goto bail;
		hid_t dataspace = H5Dget_space(dataset); if (dataspace < 0) goto bail;
		H5Sget_simple_extent_dims(dataspace, dims, NULL);
		
		hsize_t outdims[2] = { RegionSize, 4 };
		hid_t mspace_id = H5Screate_simple(2, outdims, NULL);
		
		hsize_t hyperslab_start[2] = { Region->start, 1};
		hsize_t hyperslab_count[2] = { RegionSize, 4};
		
		H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
		if (H5Dread(dataset, H5T_NATIVE_INT, mspace_id, dataspace, H5P_DEFAULT, newJob->HReg.Regions) < 0) goto bail;
			
		H5Sclose(mspace_id);
		H5Sclose(dataspace);
		H5Dclose(dataset);
		
		/* Get the sequence range */
		int max = 0;
		{
			_Bool HasValidHQ = true;
			int HQindex = 0;
			HoleRegion_t * const HR = newJob->HReg.Regions;
			/* WARNING: That assumes inserts are always prior to HQregion !! */
			for (int r=0; r<RegionSize; r++) {
				if (HR[r].type == Insert) {
					const int lstop = HR[r].stop;
					max = (max > lstop) ? max : lstop;
				}
				else if (HR[r].type == HQRegion) {
					if (HR[r].stop == 0) {
						HasValidHQ = false;
					}
					HQindex = r;
				}
			}
			
			/* 1. Do we keep regular */
			if (HasValidHQ && (Options->Selection & PB_DISPATCH_INVALID)) continue;
			
			/* 2. Do we keep invalid if flag s such ?*/
			else if (!HasValidHQ && !(Options->Selection & PB_DISPATCH_INVALID) && !(Options->Selection & PB_KEEP_INVALID)) continue; 
			
			/* 3. Filter */
			if (HasValidHQ && !(Options->Selection & PB_DISCARD_FILTER) && (Options->Selection & PB_HAS_FILTER_SCORE)) {
				if (HR[HQindex].quality < Options->minReadAccuracy || HR[HQindex].quality > Options->maxReadAccuracy ) continue;
			}
					
			/* 4. Correct regions in case of invalid */
			if (!HasValidHQ) {
				HR[HQindex].stop = max;
			}
			
			/* 5. Neglect HQ region */
			if (Options->Selection & PB_DISCARD_HQREGION) {
				HR[HQindex].start = 0;
				if (max > HR[HQindex].stop) HR[HQindex].stop = max;
			}
			
			/* Hack the regions to get only one big stretch */
			HR[0].type = Insert;
			HR[0].start = 0;
			HR[0].stop = max;
			HR[1].type = HQRegion;
			HR[1].start = HR[HQindex].start;
			HR[1].stop =  HR[HQindex].stop;
			newJob->HReg.nRegions = 2;
		}
		
		const size_t msize = max;// - min + 2;
		if (newJob->Reads.Size < msize) {
			newJob->Reads.Stream = (unsigned char*) realloc(newJob->Reads.Stream, msize*sizeof(unsigned char));
			if (newJob->Reads.Stream == NULL) { res = -5; goto bail; }
			newJob->Reads.Size = msize;
		}
		memset(newJob->Reads.Stream, 0, msize*sizeof(char));
		/* Read sequence string */
		
		dataset = H5Dopen2(PB->Parts[FileNum], "PulseData/BaseCalls/Basecall", H5P_DEFAULT); if (dataset < 0) goto bail;;
		dataspace = H5Dget_space(dataset); if (dataspace < 0) goto bail;;
		H5Sget_simple_extent_dims(dataspace, dims, NULL);
		
		outdims[0] = msize;
		outdims[1] =  0;
		mspace_id = H5Screate_simple(1, outdims, NULL);
		
		hyperslab_start[0] = Region->BaseCallsOffset;// + min;
		hyperslab_start[1] = 1;
		hyperslab_count[0] = msize;
		hyperslab_count[1] = 1;
		
		H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
		if (H5Dread(dataset, H5T_NATIVE_UCHAR, mspace_id, dataspace, H5P_DEFAULT, newJob->Reads.Stream) < 0) goto bail;
		
		H5Sclose(mspace_id);
		H5Sclose(dataspace);
		H5Dclose(dataset);
		
		newJob->HReg.HoleNumber     = HoleNumber;
		newJob->HReg.Coordinates[0] = PB->Coordinates[HoleNumber][0];
		newJob->HReg.Coordinates[1] = PB->Coordinates[HoleNumber][1];
		
		pthread_mutex_lock(&thpool_p->jobqueue_p.rwmutex);
		jobqueue_push(&thpool_p->jobqueue_p, (job_t*) newJob);
		pthread_mutex_unlock(&thpool_p->jobqueue_p.rwmutex);
		
		/* Get a job slot from the pool donequeue */
		again2:
		pthread_mutex_lock(&thpool_p->donequeue_p.rwmutex);
		newJob = (pb_job_t*) jobqueue_pull(&thpool_p->donequeue_p);
		pthread_mutex_unlock(&thpool_p->donequeue_p.rwmutex);
		
		if (!newJob) {
			if (thpool_p->num_threads_alive > 0) {
				/* Wait a bit and try again */
				nanosleep(&polling_interval, NULL);
				goto again2;
			}
			else {
				res = -1;
				goto bail;
			}
		}
	}
	
	res = SUCCESS;
	
	bail:;
	
	if (newJob) {
		pthread_mutex_lock(&thpool_p->donequeue_p.rwmutex);
		jobqueue_push(&thpool_p->donequeue_p, (job_t*) newJob);
		pthread_mutex_unlock(&thpool_p->donequeue_p.rwmutex);
	}
	
	return res;
}

static int thpool_add_consensus_work(threadpool_t * const restrict thpool_p, const PacBio_t * const restrict PB, 
                                     const PacBioDispatchOptions_t * const restrict Options)
{
	if (!PB) return -1;
	if (!PB->Content & CONSENSUS_BASECALLING) return -2;

	/* Get consensus read length */
	hsize_t lengthdims;
	hid_t lengthdataset = H5Dopen2(PB->Parts[0], "PulseData/ConsensusBaseCalls/ZMW/NumEvent", H5P_DEFAULT); 
	if (lengthdataset < 0) return ERROR;
	hid_t lengthdataspace = H5Dget_space(lengthdataset); if (lengthdataspace < 0) return ERROR;
	H5Sget_simple_extent_dims(lengthdataspace, &lengthdims, NULL);
	
	if (lengthdims != PB->nHoles) return ERROR;
	
	int * const restrict ConsensusBaseCallingLength = (int * const restrict) malloc(PB->nHoles*sizeof(int));
	if (ConsensusBaseCallingLength == NULL) return ERROR;
	
	if (H5Dread(lengthdataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ConsensusBaseCallingLength) < 0) return ERROR;
	
	H5Dclose(lengthdataset);
	H5Sclose(lengthdataspace);
	
	pb_job_t * restrict newJob;
	const size_t N = Options->nZMW;
	/* Process holes */
	for (size_t iJob = 0; iJob < N; iJob++) {
		const unsigned int HoleNumber = Options->ZMW[iJob];
		RegionIndex_t * const restrict Region = &(PB->ConsensusBaseCalling.Regions[HoleNumber]);
		
		/* Get the number of regions in hole number */
		const int RegionSize = 2; /* Read + HQ region */
		
		/* Get a job slot from the pool donequeue */
		again:
		pthread_mutex_lock(&thpool_p->donequeue_p.rwmutex);
		newJob = (pb_job_t*) jobqueue_pull(&thpool_p->donequeue_p);
		pthread_mutex_unlock(&thpool_p->donequeue_p.rwmutex);
		
		if (!newJob) {
			if (thpool_p->num_threads_alive > 0) {
				/* Wait a bit and try again */
				nanosleep(&polling_interval, NULL);
				goto again;
			}
			else {
				return -1;
			}
		}
		
		/* Verify available memory */ 
		if (newJob->HReg.nAllocatedRegions < RegionSize) {
			newJob->HReg.Regions = (HoleRegion_t * restrict) realloc(newJob->HReg.Regions, RegionSize*sizeof(HoleRegion_t));
			if (newJob->HReg.Regions == NULL) return -4;
			newJob->HReg.nAllocatedRegions = RegionSize;
		}
		memset(newJob->HReg.Regions, 0, RegionSize*sizeof(HoleRegion_t));
		hsize_t dims[2];
		const size_t FileNum = PB->HoleFileIndex[HoleNumber];

		/* Build fake regions containing read and HQ*/
		{
			HoleRegion_t * const HR = newJob->HReg.Regions;
			HR[0].type = Insert;
			HR[0].start = 0;
			HR[0].stop = (int) ConsensusBaseCallingLength[HoleNumber];
			HR[0].quality = 0;
			HR[1].type = HQRegion;
			HR[1].start = 0;
			HR[1].stop = (int) ConsensusBaseCallingLength[HoleNumber];
			HR[1].quality = 1000; /* dummy perfect for the time being */
		}
		
		newJob->HReg.nRegions       = 2U;
		newJob->HReg.HoleNumber     = HoleNumber;
		newJob->HReg.Coordinates[0] = PB->Coordinates[HoleNumber][0];
		newJob->HReg.Coordinates[1] = PB->Coordinates[HoleNumber][1];
// 			newJob->PB                 = PB;
		
		const size_t msize = (size_t) ConsensusBaseCallingLength[HoleNumber];
		if (newJob->Reads.Size < msize) {
			newJob->Reads.Stream = (unsigned char*) realloc(newJob->Reads.Stream, msize*sizeof(unsigned char));
			if (newJob->Reads.Stream == NULL) return -5;
			newJob->Reads.Size = msize;
		}
		memset(newJob->Reads.Stream, 0, msize*sizeof(char));
		/* Read sequence string */
		
		const hid_t dataset = H5Dopen2(PB->Parts[FileNum], "PulseData/ConsensusBaseCalls/Basecall", H5P_DEFAULT); if (dataset < 0) return ERROR;
		const hid_t dataspace = H5Dget_space(dataset); if (dataspace < 0) return ERROR;
		H5Sget_simple_extent_dims(dataspace, dims, NULL);
		
		hsize_t outdims[2];
		outdims[0] = msize;
		outdims[1] =  0;
		const hid_t mspace_id = H5Screate_simple(1, outdims, NULL);
		
		hsize_t hyperslab_start[2] = { Region->BaseCallsOffset, 1};
		hsize_t hyperslab_count[2] = { msize, 1};
		
		H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
		if (H5Dread(dataset, H5T_NATIVE_UCHAR, mspace_id, dataspace, H5P_DEFAULT, newJob->Reads.Stream) < 0) return ERROR;
		
		H5Sclose(mspace_id);
		H5Sclose(dataspace);
		H5Dclose(dataset);
		
		pthread_mutex_lock(&thpool_p->jobqueue_p.rwmutex);
		jobqueue_push(&thpool_p->jobqueue_p, (job_t*) newJob);
		pthread_mutex_unlock(&thpool_p->jobqueue_p.rwmutex);
	}
	
	free(ConsensusBaseCallingLength);
	
	return SUCCESS;
}

int dispatchPacBio(const struct Profile * const restrict prf,
                   const PacBio_t * const restrict PB,
#ifdef PRF_CORE_PCRE
                   struct RegEx * const restrict Regex,
#endif
                   const OutputType_t * const restrict OutputType,
                   const PacBioDispatchOptions_t * Options,
                   const size_t nCPUs)
{
		enum PBContent WhichBaseCalling = (Options->Selection & PB_DISPATCH_CONSENSUS) ? CONSENSUS_BASECALLING : BASECALLING;
		if (!(PB->Content & WhichBaseCalling)) return ERROR;
		
		int res= SUCCESS;
		if (OutputVerbose) fputs("Dispatching PacBio data...\n", stderr); 
		
		/*************************************************************************/
		/*                         CONSISTENCY CHECK                             */
		/*************************************************************************/
		if ((Options->Selection & PB_BEST_SUBREAD_OF_ZMW) && !(Options->Selection & PB_DISPATCH_SUBREADS)) {
			fputs("Aggregation is only possible on subreads!\n", stderr);
			return -10;
		}
		/*************************************************************************/
    /*                         CHOOSE THREAD FUNCTION                        */
    /*************************************************************************/
#if (defined(PRF_CORE_STD) || defined(PRF_CORE_REPEAT))
#if (defined(PRF_CORE_STD) && defined(PRF_CORE_REPEAT))
		const Compute_t * CoreCompute = (prf->isCircular) ? &Repeat_sse41 : &Standard_sse41;
#elif defined(PRF_CORE_STD)
		const Compute_t * CoreCompute = &Standard_sse41;
#else
		const Compute_t * CoreCompute = &Repeat_sse41;
#endif
#else
#error "Dispatching procedures requires at lest STD or REPEAT methods"
#endif
	
		void * (*ThreadFct)(threadarg_t* const restrict) = NULL;

		if (Options->Selection & PB_TEST_OUTPUT) {
			ThreadFct = tp_pb_filtertest;
		}
		else {
			const int index = GetDispatchThreadIndex(OutputType, prf);
			if (index >= 0) {
				ThreadFct = (Options->Selection & PB_BEST_SUBREAD_OF_ZMW) ? tp_pbzp[index] : tp_pbsp[index];
			}
			else {
				fputs("Error in the choice of thread function\n", stderr);
				return -1;
			}
		}
		   
    /*************************************************************************/
    /*                       EXTRA VARIABLES FOR THREADS                     */
    /*************************************************************************/
		size_t HistogramsMemSize = 0UL;
		if (OutputType->Type == HISTOGRAM) {
			if (OutputType->Specific.Histogram.CycleRatherThanScore) {
				HistogramsMemSize = OutputType->CycleRange[1] - OutputType->CycleRange[0] + 1; 
			}
			else {
				HistogramsMemSize = OutputType->ScoreRange[1] - OutputType->ScoreRange[0] + 1; 
			}
		}
		else if(OutputType->Type == DENSITY) {
			HistogramsMemSize = (OutputType->CycleRange[1] - OutputType->CycleRange[0] + 1)\
												*(OutputType->ScoreRange[1] - OutputType->ScoreRange[0] + 1);
		}
		
		size_t * const Histograms = (size_t *) calloc(HistogramsMemSize*nCPUs, sizeof(size_t));
		if (Histograms == NULL) {
			fprintf(stderr, "Unable to allocate memory for histograms, requested size was %lu bytes.\n",
							HistogramsMemSize*nCPUs*sizeof(size_t));
			return -1;
		}
		size_t * const Missed = (size_t *) calloc(nCPUs, sizeof(size_t));
		if (Missed == NULL) {
			fprintf(stderr, "Unable to allocate memory for missed values, requested size was %lu bytes.\n",
							nCPUs*sizeof(size_t));
			if (Histograms) free(Histograms);
			return -2;
		}

		/*************************************************************************/
    /*                    SET THE THREADS COMMON VARIABLES                   */
    /*************************************************************************/
		pb_common_t common;
		common.profile = prf;
		common.Histograms = Histograms;
		common.Counters = Missed;
		common.Compute = CoreCompute;
		common.OutputType = OutputType;
#ifdef PRF_CORE_PCRE
		common.regex = Regex;
#endif
		pthread_mutex_init(&(common.PrintLock), NULL);

		/*************************************************************************/
    /*                          PREPARE THREAD POOL                          */
    /*************************************************************************/  
#ifdef PRF_USE_AFFINITY
		threadpool_t * const restrict thpool = createThreadPool(ThreadFct, NULL, sizeof(pb_job_t), (void*) &common, nCPUs,0);
#else
		threadpool_t * const restrict thpool = createThreadPool(ThreadFct, sizeof(pb_job_t), (void*) &common, nCPUs,0);
#endif
		if ( !thpool ) {
			if (Histograms) free(Histograms);
			if (Missed) free(Missed);
			return -3;
    }
         
    /*************************************************************************/
		/*                GETTING AND ANALYSING INPUT BAM FILE                   */
		/*************************************************************************/

	
		/*************************************************************************/
    /*                        ADD TASKS TO JOB QUEUE                         */
    /*************************************************************************/
		struct timeval _t0, _t1;
		gettimeofday(&_t0,0);
		
		/* Standard reads basecalling */
		if (Options->Selection & PB_DISPATCH_ZMW) {
			res = thpool_add_read_work(thpool, PB, Options);
		}
		else if (Options->Selection & PB_DISPATCH_CONSENSUS) {
			res = thpool_add_consensus_work(thpool, PB, Options);
		}
		else {
			res = thpool_add_subread_work(thpool, PB, Options);
		}
		
    if (res != 0 ) {
			fprintf(stderr, "Master failed to send jobs, error code %i\n", res);
			goto KILL_THREADS;
		}
		
		/*************************************************************************/
    /*                       WAIT FOR TASK TO FINISH                         */
    /*************************************************************************/
    thpool_wait(thpool);
		gettimeofday(&_t1,0);
		
		if (OutputVerbose) {
			const double t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr,"This took %lf seconds to crunch on %lu cores.\n", t, nCPUs);

		}
		    
    /*************************************************************************/
    /*                  TELL THREADS TO FLUSH AND TERMINATE                  */
    /*************************************************************************/


    /*************************************************************************/
    /*                          GATHER_EXTRA_DATA                            */
    /*************************************************************************/ 
		if (HistogramsMemSize != 0) {
			size_t MissedValues = Missed[0];
			const size_t * restrict hptr = &Histograms[HistogramsMemSize];
			for (size_t i=1;i<nCPUs;i++) {
				MissedValues += Missed[i];
				for (size_t j=0; j<HistogramsMemSize; j++) Histograms[j] += hptr[j];
				hptr += HistogramsMemSize;
			}
			if (MissedValues) {
				fprintf(stderr,
								"Some sequences (%lu) bear alignment that are longer than the histogram bin number\n",
								MissedValues);
			}
			char FName[256];
			hptr = Histograms;
			if (OutputType->Type == HISTOGRAM) {
				snprintf(FName, 256, "%s.histogram", OutputType->Specific.Histogram.BaseFileName);
				FILE * const out = fopen(FName, "w");
				const int RangeStart = (OutputType->Specific.Histogram.CycleRatherThanScore) ? OutputType->CycleRange[0] :\
																OutputType->ScoreRange[0];
				if (out != NULL) {
					for (size_t j=0; j<HistogramsMemSize; j++) {
						fprintf(out, "%i\t%lu\n", RangeStart+(int)j, Histograms[j]);
					}
					fclose(out);
				}
				else {
					fprintf(stderr, "Unable to create output histogram %s\n", FName);
					res = -5;
				}
			}
			else if (OutputType->Type == DENSITY) {
				snprintf(FName, 256, "%s.density", OutputType->Specific.Histogram.BaseFileName);
				FILE * const out = fopen(FName, "w");
				if (out != NULL) {
					hptr = Histograms;
					for (int i=OutputType->ScoreRange[0]; i<OutputType->ScoreRange[1]; i++) {
						for (unsigned int j=OutputType->CycleRange[0]; j<OutputType->CycleRange[1]; j++) {
							fprintf(out, "%i\t%u\t%lu\n", i, j,* hptr++);
						}
						fprintf(out, "\n");
					}
					fclose(out);
				}
				else {
					fprintf(stderr, "Unable to create output histogram %s\n", FName);
					res = -5;
				}
			}
		}
		

    /*************************************************************************/
    /*                    CHECK FOR THREADS ERROR                            */
    /*************************************************************************/

       
    /*************************************************************************/
    /*                      DESTROY THE THREAD POOL                          */
    /*************************************************************************/
    KILL_THREADS:
		destroyThreadPool(thpool);
    pthread_mutex_destroy(&(common.PrintLock));
    
    FIN0:
		if (Histograms) free(Histograms);
    FIN2:
		return res;
}

int dispatchPacBioExt(const PacBio_t * const restrict PB,
											threadpool_t * const restrict thpool,
                      const PacBioDispatchOptions_t * Options
                     )
{
		enum PBContent WhichBaseCalling = (Options->Selection & PB_DISPATCH_CONSENSUS) ? CONSENSUS_BASECALLING : BASECALLING;
		if (!(PB->Content & WhichBaseCalling)) return ERROR;
		
		int res= SUCCESS;
		if (OutputVerbose) fputs("Dispatching PacBio data...\n", stderr); 
		
		/*************************************************************************/
		/*                         CONSISTENCY CHECK                             */
		/*************************************************************************/
		if ((Options->Selection & PB_BEST_SUBREAD_OF_ZMW) && !(Options->Selection & PB_DISPATCH_SUBREADS)) {
			fputs("Aggregation is only possible on subreads!\n", stderr);
			return -10;
		}
		/*************************************************************************/
    /*                         CHOOSE THREAD FUNCTION                        */
    /*************************************************************************/
		   
    /*************************************************************************/
    /*                       EXTRA VARIABLES FOR THREADS                     */
    /*************************************************************************/

		/*************************************************************************/
    /*                    SET THE THREADS COMMON VARIABLES                   */
    /*************************************************************************/

		/*************************************************************************/
    /*                          PREPARE THREAD POOL                          */
    /*************************************************************************/  
         
    /*************************************************************************/
		/*                GETTING AND ANALYSING INPUT BAM FILE                   */
		/*************************************************************************/
	
		/*************************************************************************/
    /*                        ADD TASKS TO JOB QUEUE                         */
    /*************************************************************************/
		struct timeval _t0, _t1;
		gettimeofday(&_t0,0);
		
		/* Standard reads basecalling */
		if (Options->Selection & PB_DISPATCH_ZMW) {
			res = thpool_add_read_work(thpool, PB, Options);
		}
		else if (Options->Selection & PB_DISPATCH_CONSENSUS) {
			res = thpool_add_consensus_work(thpool, PB, Options);
		}
		else {
			res = thpool_add_subread_work(thpool, PB, Options);
		}
		
    if (res != 0 ) {
			fprintf(stderr, "Master failed to send jobs, error code %i\n", res);
			goto KILL_THREADS;
		}
		
		/*************************************************************************/
    /*                       WAIT FOR TASK TO FINISH                         */
    /*************************************************************************/
    thpool_wait(thpool);
		gettimeofday(&_t1,0);
		
		if (OutputVerbose) {
			const double t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr,"This took %lf seconds to crunch on %i cores.\n", t, thpool->num_threads);

		}
		    
    /*************************************************************************/
    /*                  TELL THREADS TO FLUSH AND TERMINATE                  */
    /*************************************************************************/


    /*************************************************************************/
    /*                          GATHER_EXTRA_DATA                            */
    /*************************************************************************/ 

    /*************************************************************************/
    /*                    CHECK FOR THREADS ERROR                            */
    /*************************************************************************/

       
    /*************************************************************************/
    /*                      DESTROY THE THREAD POOL                          */
    /*************************************************************************/
    KILL_THREADS:
    
    FIN0:
    FIN2:
		return res;
}
