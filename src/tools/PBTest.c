#include "prf_config.h" 
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <sys/time.h>
#include "pb_hdf5.h"
#include "../dispatch/PacBio/pb_threads.h"

///////////////////////////////////////////////////////////////////
// THREAD STRUCTURE DECLARATIONS
void* thread(threadarg_t* const restrict arg)
{
	int res = SUCCESS;
	//     printf("Thread started fct THREAD_FCT_NAME\n");
    
    /*************************************************************************/
    /*                          CONFIGURE WORKER                             */
    /*************************************************************************/
		
		/* Get thread pool */
		threadpool_t * const restrict thpool = arg->thpool;
		
    /* Register signal handler */
//     struct sigaction act;
//     act.sa_handler = thread_hold;
//     if (sigaction(SIGUSR1, &act, NULL) == -1) {
// 	    fprintf(stderr, "thread_do(): cannot handle SIGUSR1");
//     }
    
    /* Mark thread as alive (initialized) */
    pthread_mutex_lock(&thpool->thcount_lock);
    thpool->num_threads_alive += 1;
    pthread_mutex_unlock(&thpool->thcount_lock);
    
    
    // Get the pointers
//     const ThreadPool_Common_t * const restrict common = arg->common;
        
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
						printf("Thread working on hole %u, %u regions found\n", HoleNumber, nRegions);
						restrictReadstoHQ(HReg);
						for (int i=0; i<nRegions; i++) {
							if (HReg->Regions[i].type == Insert) {
								printf("region %i: %i-%i\n", i, HReg->Regions[i].start, HReg->Regions[i].stop);
								for (int j=HReg->Regions[i].start; j<HReg->Regions[i].stop; j+=80) {
									printf("\t%.*s\n", (j+80)<HReg->Regions[i].stop ? 80 : HReg->Regions[i].stop - j,
									                   &(task->Reads.Stream[j]));
								}
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
     
    FIN:
    pthread_mutex_lock(&thpool->thcount_lock);
    thpool->num_threads_alive--;
    pthread_mutex_unlock(&thpool->thcount_lock);
    
    if( res != SUCCESS) {
      fprintf(stderr, "Thread error code %i\n", res);
      if (task) {
				pthread_mutex_lock(&thpool->thcount_lock);
				thpool->num_threads_working--;
				pthread_mutex_unlock(&thpool->thcount_lock);
				jobqueue_clear(&thpool->jobqueue_p);
      }
      thpool->threads_keepalive = 0;
    }
    arg->returnState = (void*) (intptr_t) res;
    return NULL;
}


_Bool PacBioKeepInvalid = false;
int PacBioHQThreshold = 0; 

static int thpool_add_work(threadpool_t * const restrict thpool_p, const PacBio_t * const restrict PB, 
                           const unsigned int * const restrict holes, const size_t N)
{
	if (!PB) return -1;
	if (!PB->BaseCalling.Regions) return -2;

	pb_job_t * restrict newJob; 
	for (size_t iJob = 0; iJob < N; iJob++) {
		const unsigned int HoleNumber = holes[iJob];
		RegionIndex_t * const restrict Region = &(PB->BaseCalling.Regions[HoleNumber]);
		
		/* Get the number of regions in hole number */
		const int RegionSize = Region->stop - Region->start + 1;

		/* Get a job slot from the pool donequeue */
		pthread_mutex_lock(&thpool_p->donequeue_p.rwmutex);
		newJob = (pb_job_t*) jobqueue_pull(&thpool_p->donequeue_p);
		pthread_mutex_unlock(&thpool_p->donequeue_p.rwmutex);
		
		if (newJob == NULL) {
			return -2;
		}
		
		if (newJob->HReg.nAllocatedRegions < RegionSize) {
			newJob->HReg.Regions = (HoleRegion_t * restrict) realloc(newJob->HReg.Regions, RegionSize*sizeof(HoleRegion_t));
			if (newJob->HReg.Regions == NULL) return -4;
			newJob->HReg.nAllocatedRegions = RegionSize;
		}
		memset(newJob->HReg.Regions, 0, RegionSize*sizeof(HoleRegion_t));
		hsize_t dims[2];
		const size_t FileNum = PB->HoleFileIndex[HoleNumber];

		hid_t dataset = H5Dopen2(PB->Parts[FileNum], "PulseData/Regions", H5P_DEFAULT); if (dataset < 0) return ERROR;
		hid_t dataspace = H5Dget_space(dataset); if (dataspace < 0) return ERROR;
		H5Sget_simple_extent_dims(dataspace, dims, NULL);
		
		hsize_t outdims[2] = { RegionSize, 4 };
		hid_t mspace_id = H5Screate_simple(2, outdims, NULL);
		
		hsize_t hyperslab_start[2] = { Region->start, 1};
		hsize_t hyperslab_count[2] = { RegionSize, 4};
		
		H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
		if (H5Dread(dataset, H5T_NATIVE_INT, mspace_id, dataspace, H5P_DEFAULT, newJob->HReg.Regions) < 0) return ERROR;
		
		newJob->HReg.nRegions       = (unsigned int) RegionSize;
		newJob->HReg.HoleNumber     = HoleNumber;
		newJob->HReg.Coordinates[0] = PB->Coordinates[HoleNumber][0];
		newJob->HReg.Coordinates[1] = PB->Coordinates[HoleNumber][1];
		
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
			
			if (HR[HQindex].quality < PacBioHQThreshold) {
				HR[HQindex].stop = 0;
			}
			
			if (PacBioKeepInvalid && !HasValidHQ) {
				HR[HQindex].stop = max;
			}
		}
		
		const size_t msize = max;// - min + 2;
		if (newJob->Reads.Size < msize) {
			newJob->Reads.Stream = (unsigned char*) realloc(newJob->Reads.Stream, msize*sizeof(char));
			if (newJob->Reads.Stream == NULL) return -5;
			newJob->Reads.Size = msize;
		}
		memset(newJob->Reads.Stream, 0, msize*sizeof(char));
		/* Read sequence string */
		
		dataset = H5Dopen2(PB->Parts[FileNum], "PulseData/BaseCalls/Basecall", H5P_DEFAULT); if (dataset < 0) return ERROR;
		dataspace = H5Dget_space(dataset); if (dataspace < 0) return ERROR;
		H5Sget_simple_extent_dims(dataspace, dims, NULL);
		
		outdims[0] = msize;
		outdims[1] =  0;
		mspace_id = H5Screate_simple(1, outdims, NULL);
		
		hyperslab_start[0] = Region->BaseCallsOffset;// + min;
		hyperslab_start[1] = 1;
		hyperslab_count[0] = msize;
		hyperslab_count[1] = 1;
		
		H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
		if (H5Dread(dataset, H5T_NATIVE_UCHAR, mspace_id, dataspace, H5P_DEFAULT, newJob->Reads.Stream) < 0) return ERROR;
		
		H5Sclose(mspace_id);
		H5Sclose(dataspace);
		H5Dclose(dataset);
		
		pthread_mutex_lock(&thpool_p->jobqueue_p.rwmutex);
		jobqueue_push(&thpool_p->jobqueue_p, (job_t*) newJob);
		pthread_mutex_unlock(&thpool_p->jobqueue_p.rwmutex);
	}
	
	return SUCCESS;
}

int master(const int start, const int stop, const size_t nthreads, const PacBio_t * const restrict PB, 
           unsigned int * const restrict holes, const size_t N )
{
		int res= SUCCESS;
    /*************************************************************************/
    /*                    SET THE TREADS COMMON VARIABLES                    */
    /*************************************************************************/  
		{
    }
    
    /*************************************************************************/
    /*                          PREPARE THREAD POOL                          */
    /*************************************************************************/  
    
    // Make new thread pool
#ifdef PRF_USE_AFFINITY
    threadpool_t * const thpool = createThreadPool( thread, NULL, sizeof(pb_job_t), NULL, 2, 0);
#else
    threadpool_t * const thpool = createThreadPool( thread, sizeof(pb_job_t), NULL, 2, 0);
#endif
    if (!thpool) {
			res = -2;
			goto FIN2;
    }
    
    /*************************************************************************/
    /*                        ADD TASKS TO JOB QUEUE                         */
    /*************************************************************************/

		res = thpool_add_work(thpool, PB, holes, N);
    if (res != SUCCESS) {
			res = -5; goto FIN3;
		}
            
    /*************************************************************************/
    /*                       WAIT FOR TASK TO FINISH                         */
    /*************************************************************************/
    thpool_wait(thpool);
        
    /*************************************************************************/
    /*                      DESTROY THE THREAD POOL                          */
    /*************************************************************************/
		FIN3:;
		destroyThreadPool(thpool);
		
    FIN2:
		return res;
}

int main(int argc, char * argv[])
{
	if (argc != 3) {
		fprintf(stderr, "Usage: %s [pac bio HDF5 file] [hole number]\n", argv[0]);
		exit(1);
	}
	
	fprintf(stderr, "Testing file %s\n", argv[1]);
	PacBio_t * const restrict PBS = PacBioOpen(argv[1]);
	if (PBS == NULL) {
		fprintf(stderr, "Error opening Pac Bio file %s\n", argv[1]);
		exit(1);
	}
	
	fprintf(stderr,"Base file name : %s\n"
				 "Directory      : %s\n"
				 "Parts          : %u\n",
				PBS->BaseFileName, PBS->Directory, PBS->nParts);
	for (unsigned int i=0; i<PBS->nParts; i++) {
		fprintf(stderr,"               : %s\n", PBS->PartsFileName[i]);
	}
	
	fputs("The given file [set] contains : ", stderr);
	if (PBS->Content & CONSENSUS_BASECALLING) fputs("CONSENSUS ", stderr);
	if (PBS->Content & BASECALLING) fputs("BASECALLING ", stderr);
	if (PBS->Content & PULSE) fputs("PULSES ", stderr);
	if (PBS->Content & TRACE) fputs("TRACES ", stderr);
	fputs("\n", stderr);
	
	fprintf(stderr,"Holes          : %u\n", PBS->nHoles);
	if (PBS->Content & TRACE) {
		
		if (PopulateHoleStatus(PBS) != SUCCESS) {
				fprintf(stderr,"Hole status not populated...\n");
				goto END;
		}
		
		if (PopulateHoleCoordinates(PBS) != SUCCESS) {
			fprintf(stderr,"Hole coordinates not populated...\n");
			goto END;
		}
		
		if (PopulateDecodeTable(PBS) != SUCCESS) {
			fputs("Error getting Decoding tables and Bias\n", stderr);
			goto END;
		}
		for (unsigned int i=0; i<PBS->nParts; i++) {
			printf("%s decode table:\n", PBS->PartsFileName[i]);
			for (int l=0; l<16; l++) {
				for (int k=0; k<7; k++) {
					printf("%8.2f ", PBS->Traces.DecodeTable[i][8*l+k]);
				}
				printf("%8.2f\n", PBS->Traces.DecodeTable[i][8*l+7]);
			}
			printf("\n%s decode bias: %8.2f\n\n", PBS->PartsFileName[i], PBS->Traces.Bias[i]);
		}
		int start, stop;
		if (sscanf(argv[2], "%i:%i", &start, &stop) != 2) {
			const int holeN = atoi(argv[2]);
			
			fprintf(stderr,"Hole %i is at location [%i,%i]\n", holeN, PBS->Coordinates[holeN][0], PBS->Coordinates[holeN][1]);
			fputs("This hole is of type ", stderr);
			switch(PBS->Status[holeN]) {
				case SEQUENCING: fputs("SEQUENCING\n", stderr); break;
				case ANTIHOLE: fputs("ANTIHOLE\n", stderr); break;
				case FIDUCIAL: fputs("FIDUCIAL\n", stderr); break;
				case SUSPECT: fputs("SUSPECT\n", stderr); break;
				case ANTIMIRROR: fputs("ANTIMIRROR\n", stderr); break;
				case FDZMW: fputs("FDZMW\n", stderr); break;
				case FBZMW: fputs("FBZMW\n", stderr); break;
				case ANTIBEAMLET: fputs("ANTIBEAMLET\n", stderr); break;
				case OUTSIDEFOV: fputs("OUTSIDEFOV\n", stderr); break;
			}

			__m128* const restrict __tr = getTrace(PBS, (unsigned int) holeN);
			
			FILE* out = fopen("Trace.dat","w");
			if (out != NULL) {
				float(* const restrict ftmp)[4] = (float (* const restrict)[4]) __tr; 
				for (size_t i=0; i<PBS->Traces.Length; i++) {
					fprintf(out, "%7lu\t%10.5f\t%10.5f\t%10.5f\t%10.5f\n",i, ftmp[i][0], ftmp[i][1], ftmp[i][2], ftmp[i][3]);
				}
				fclose(out);
				fprintf(stderr, "Hole %i traces exported to Traces.dat\n", holeN);
			}
			
			_mm_free(__tr);
			
			unsigned int Stats[256][256];
			
			Traces_t * const Indices = getTraceIndex(PBS, holeN); 
			unsigned char * restrict IndicesPtr = Indices->memory;
			
			const size_t nTimes = PBS->Traces.Length;
			for (int i=0; i<4; i++) {
				memset(Stats, 0, 256*256*sizeof(unsigned int));
				register unsigned char previous = IndicesPtr[0]; 
				for (size_t t=1; t<nTimes; t++) {
					const unsigned char Id = IndicesPtr[t];
					Stats[previous][Id]++;
					previous = Id;
				}
				char FName[16];
				snprintf(FName, 16, "Stat_%i.dat", i);
				out = fopen(FName, "w");
				if (out != NULL) {
					for (int l=0; l<256; l++) {
						for (int m=0;m<256;m++) {
							fprintf(out, "%u %u %u\n", l, m, Stats[l][m]); 
						}
						fputc('\n', out);
					}
					fclose(out);
				}
				IndicesPtr += Indices->stride;
			}
			
			freeTraces(Indices);
		}
		
	} 
	else if (PBS->Content & BASECALLING) {
		IndexRegions(PBS);
		
		int start, stop;
		if (sscanf(argv[2], "%i:%i", &start, &stop) != 2) {
			HoleData_t HReg = HOLE_DATA_INIT;
			const int holeN = atoi(argv[2]);

			int count = getHoleRegions(PBS, &HReg, holeN);
			if (count>0) {
				fprintf(stderr,"Regions     type    start     stop  quality\n");
				for (int i=0; i<count; i++) {
					fprintf(stderr,"%7i%9i%9i%9i%9i\n", i, HReg.Regions[i].type, HReg.Regions[i].start, HReg.Regions[i].stop, HReg.Regions[i].quality);
				}
				
				restrictReadstoHQ(&HReg);
				fprintf(stderr,"Restricting to HQ region\nRegions     type    start     stop  quality\n");
				for (int i=0; i<count; i++) {
					fprintf(stderr,"%7i%9i%9i%9i%9i\n", i, HReg.Regions[i].type, HReg.Regions[i].start, HReg.Regions[i].stop, HReg.Regions[i].quality);
				}
			}
			
			HoleReads_t HReads = HOLE_READS_INIT;
			HoleQVs_t HQVs = HOLE_QVS_INIT;
			getHoleReads_HN(PBS, &HReads, holeN);
			getHoleQVs_HN(PBS, &HQVs, holeN);
			const int NameSize = strlen(PBS->BaseFileName) - 7;
			
			if (PopulateHoleStatus(PBS) != SUCCESS) {
				fprintf(stderr,"Hole status not populated...\n");
			}
			
			/* Get RQ */
			int l = count;
			int HQ = 0;
			while (l>=0) {
				if (HReg.Regions[l].type == HQRegion) {HQ = HReg.Regions[l].quality; break;}
				--l;
			}
			
			if (HQ >= 0 && PBS->Status[holeN] == (unsigned char) SEQUENCING) {
				for (size_t i=0; i<HReads.Size; i++) {
					HReads.Stream[i] = (HReads.Stream[i] >= 'a') ? HReads.Stream[i] - 'a' + 'A' : HReads.Stream[i];
				}
				printf("SMRT adapters:\n");
				for (int i=0; i<count; i++) {
					if (HReg.Regions[i].type == Adapter) {
	// 					if (HReg.Regions[i].stop - HReg.Regions[i].start < 500) continue;
	// 					for (int j=HReg.Regions[i].start; j<HReg.Regions[i].stop; j++) {
	// 							HReads.Stream[j] = (unsigned char) HReads.Stream[j] - 'A' + 'a';
	// 					}
						printf(">%.*s/%i/%i_%i RQ=0.%3.3i\n",
									 NameSize, PBS->BaseFileName, holeN,
									 HReg.Regions[i].start, HReg.Regions[i].stop,
									 HReg.Regions[i].quality );
						for (int j=HReg.Regions[i].start; j<HReg.Regions[i].stop; j+=80) {
							printf("%.*s\n", (j+80)<HReg.Regions[i].stop ? 80 : HReg.Regions[i].stop - j,
							                   &(HReads.Stream[j]));
						}
					}
				}
				
				printf("Subreads:\n");
				for (int i=0; i<count; i++) {
					if (HReg.Regions[i].type == Insert) {
	// 					if (HReg.Regions[i].stop - HReg.Regions[i].start < 500) continue;
	// 					for (int j=HReg.Regions[i].start; j<HReg.Regions[i].stop; j++) {
	// 							HReads.Stream[j] = (unsigned char) HReads.Stream[j] - 'A' + 'a';
	// 					}
						printf(">%.*s/%i/%i_%i RQ=0.%3.3i\n",
									 NameSize, PBS->BaseFileName, holeN,
									 HReg.Regions[i].start, HReg.Regions[i].stop,
									 HQ );
						for (int j=HReg.Regions[i].start; j<HReg.Regions[i].stop; j+=80) {
							printf("%.*s\n", (j+80)<HReg.Regions[i].stop ? 80 : HReg.Regions[i].stop - j,
							                   &(HReads.Stream[j]));
						}
					}
				}
			}
			
			freeHoleReads(&HReads);
			freeHoleQVs(&HQVs);
			
			
			if (PopulateHoleCoordinates(PBS) != SUCCESS) {
				fprintf(stderr,"Hole coordinates not populated...\n");
			}
			else {
				fprintf(stderr,"Hole %i is at location [%i,%i]\n", holeN, PBS->Coordinates[holeN][0], PBS->Coordinates[holeN][1]);
				fputs("This hole is of type ", stderr);
				switch(PBS->Status[holeN]) {
					case SEQUENCING: fputs("SEQUENCING\n", stderr); break;
					case ANTIHOLE: fputs("ANTIHOLE\n", stderr); break;
					case FIDUCIAL: fputs("FIDUCIAL\n", stderr); break;
					case SUSPECT: fputs("SUSPECT\n", stderr); break;
					case ANTIMIRROR: fputs("ANTIMIRROR\n", stderr); break;
					case FDZMW: fputs("FDZMW\n", stderr); break;
					case FBZMW: fputs("FBZMW\n", stderr); break;
					case ANTIBEAMLET: fputs("ANTIBEAMLET\n", stderr); break;
					case OUTSIDEFOV: fputs("OUTSIDEFOV\n", stderr); break;
				}
			}
			
			unsigned int * SeqHole;
			count = getSequencingHoles(PBS, &SeqHole);
			if (count > 0) {
				fprintf(stderr,"Found %i sequencing holes\n", count);
			}
			else {
				fprintf(stderr,"getSequencingHoles returned %i\n", count);
		// 		for (int i=0; i<10; i++) printf("%i hole status : %u\n", i, SeqHole[i]);
			}
			
			free(SeqHole);
		}
		else {
		if (PopulateHoleCoordinates(PBS) != SUCCESS) {
			fprintf(stderr,"Error in Hole coordinates not populated...\n");
			exit(1);
		}
		
		size_t count = stop-start+1;
		unsigned int * const holes = (unsigned int*) malloc(count*sizeof(unsigned int));
		if (holes == NULL) {
			fprintf(stderr, "Allocation error on holes\n");
			exit(1);
		}
		unsigned int * holesPtr = holes;
		for (int i=start; i<=stop;i++) *holesPtr++ = i;
		
		if (master(start, stop, 6, PBS, holes, count ) != 0 ) {
			fprintf(stderr, " Error in master\n");
		}
	}
	}
	
	END:
	if (PacBioClose(PBS) != SUCCESS) {
		fprintf(stderr, "Error closing Pac Bio...\n");
		exit(1);
	}
	exit(0);
}

