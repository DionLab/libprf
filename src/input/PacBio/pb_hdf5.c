#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <libgen.h>
#include <pthread.h>
#include "pb_hdf5.h"
#include <xmmintrin.h>

#define __USE_INLINE_FUNCTIONS__
#include "pfSequence.h"

const char *const HoleStatusStrings[] = { 
	"SEQUENCING",
	"ANTIHOLE",
	"FIDUCIAL",
	"SUSPECT",
	"ANTIMIRROR",
	"FDZMW",
	"FBZMW",
	"ANTIBEAMLET",
	"OUTSIDEFOV"
};

static const char * WhereHoleStatus[] = {"PulseData/BaseCalls/ZMW/HoleStatus", "TraceData/HoleStatus" };
static const char * WhereHoleCoord[]  = {"PulseData/BaseCalls/ZMW/HoleXY", "TraceData/HoleXY" };

int isPacBioH5(const char * const restrict BaseFileName)
{
	int res = 0;
	const int fid = open(BaseFileName, O_RDONLY);
	if (fid < 0) return fid;
	char buffer[4];
	if (read(fid, buffer, 4) != 4) {
		res = -2;
		goto FIN;
	}
	
	res  = ( (unsigned char) buffer[0] == 0x89 );
	res &= ( (unsigned char) buffer[1] == 0x48 );
	res &= ( (unsigned char) buffer[2] == 0x44 );
	res &= ( (unsigned char) buffer[3] == 0x46 );
	
	FIN:
	close(fid);
	return res;
}

PacBio_t* PacBioOpen(const char * const restrict BaseFileName)
{
	/* Allocate memory to hold the structure */
	PacBio_t * restrict PB = (PacBio_t*) malloc(sizeof(PacBio_t));
	if (PB == NULL) return NULL;
	memset(PB, 0, sizeof(PacBio_t));
	
	/* Check a few size requirements */
	{
		const size_t L = strlen(BaseFileName);
		char * const ctmp = alloca(L+1);
		char * ptr = strncpy(ctmp, BaseFileName, L); ctmp[L] = '\0';
		char * const directory = dirname(ptr);
		{
			if (strlen(directory) >= DIRECTORY_SIZE) {
					fprintf(stderr, "Directory name of %s exceeds allowed size of %u\n",BaseFileName, DIRECTORY_SIZE);
					free(PB);
					return NULL;
			}
		}
		strncpy(PB->Directory, directory, DIRECTORY_SIZE);
		
		ptr = strncpy(ctmp, BaseFileName, L); ctmp[L] = '\0';
		char * const name = basename(ptr);
		{
			if (strlen(name) >= DIRECTORY_SIZE) {
					fprintf(stderr, "File name of %s exceeds allowed size of %u\n",BaseFileName, FILENAME_SIZE);
					free(PB);
					return NULL;
			}
		}
		strncpy(PB->BaseFileName, name, FILENAME_SIZE);
		PB->BaseFileName[FILENAME_SIZE-1] = '\0';
	}
	
	/* 
	 * Open PacBio base file and get data out of it 
	 */
	hid_t file_id = H5Fopen(BaseFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id <0 ) return NULL;
	
	PB->Base = file_id;
	
	H5E_auto2_t old_func;
	void *old_client_data;
	H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client_data);
	H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
	_Bool isMultipart;
	{
		herr_t status = H5Gget_objinfo (file_id, "MultiPart", 0, NULL);
		isMultipart = (status == 0) ? true : false;
	}
	H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);

	/* Get multiparts filename */
	if (isMultipart) {
		hsize_t field_len;
		hid_t dataset = H5Dopen2(file_id, "MultiPart/Parts", H5P_DEFAULT); if (dataset < 0) return NULL;
		hid_t field_space = H5Dget_space(dataset); if (field_space < 0) return NULL;
		H5Sget_simple_extent_dims(field_space, &field_len, NULL);
		
		PB->nParts = (unsigned int) field_len;
		H5Sclose(field_space);
		
		char ** PartsFileName = (char**) malloc(field_len*sizeof(char*));
		if (PartsFileName == NULL) {
			PacBioClose(PB);
			return NULL;
		}
		for (hsize_t p=0; p<field_len; p++) PartsFileName[p] = NULL;
		
		hid_t dataspace = H5Screate_simple(1, &field_len, NULL);
		hid_t datatype = H5Tvlen_create(H5T_NATIVE_CHAR);
		hvl_t * const restrict rdata = alloca(field_len*sizeof(hvl_t));
		
		if (H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata) < 0) {
			H5Dclose(dataset);
			H5Sclose(dataspace);
			H5Tclose(datatype);
			PacBioClose(PB);
			return NULL;
		}
		
		for (hsize_t p=0; p<field_len; p++) {
			const size_t length = rdata[p].len;
			const char * const restrict str = (char * const restrict) rdata[p].p;
			char * const restrict ctmp = (char*) malloc((length+1)*sizeof(char));
			if (ctmp == NULL) {
				PacBioClose(PB);
				return NULL;
			}

			for (size_t i=0; i<length; i++) ctmp[i] = str[i];
			ctmp[length] = '\0';
			PartsFileName[p] = ctmp;
		}
		H5Dvlen_reclaim(datatype, dataspace, H5P_DEFAULT, rdata);
		
		PB->PartsFileName = PartsFileName;
		
		H5Dclose(dataset);
		H5Sclose(dataspace);
		H5Tclose(datatype);
	}
	else {
		PB->nParts = 1;
		char ** PartsFileName = (char**) malloc(sizeof(char*));
		if (PartsFileName == NULL) {
			PacBioClose(PB);
			return NULL;
		}
		PartsFileName[0] = malloc(FILENAME_SIZE*sizeof(char));
		if (PartsFileName[0] == NULL) {
			PacBioClose(PB);
			return NULL;
		}
		strncpy(PartsFileName[0], PB->BaseFileName, 256);
		PB->PartsFileName = PartsFileName;
	}
	
	/* Get Holes location within files */
	unsigned int * tmp;
	if (isMultipart) {
		hsize_t field_len[2];
		hid_t dataset = H5Dopen2(file_id, "MultiPart/HoleLookup", H5P_DEFAULT); if (dataset < 0) return NULL;
		hid_t field_space = H5Dget_space(dataset); if (field_space < 0) return NULL;
		H5Sget_simple_extent_dims(field_space, field_len, NULL);
		
		PB->nHoles = (unsigned int) field_len[0];
		H5Sclose(field_space);
		
		unsigned char * const restrict ID = (unsigned char * const restrict) malloc(field_len[0]*sizeof(unsigned char));
		if (ID == NULL) {
			PacBioClose(PB);
			return NULL;
		}
		PB->HoleFileIndex = ID;
		
		tmp = (unsigned int * const restrict) malloc(field_len[1]*field_len[0]*sizeof(unsigned int));
		if (tmp == NULL) {
			PacBioClose(PB);
			return NULL;
		}
		hid_t dataspace = H5Screate_simple(1, field_len, NULL);
		if (H5Dread(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp) < 0) {
			H5Dclose(dataset);
			H5Sclose(dataspace);
			PacBioClose(PB);
			return NULL;
		}
		
		const unsigned int (* restrict pHole)[2] = (const unsigned int (*)[2]) tmp;
		for (hsize_t p=0; p<field_len[0]; p++) {
			ID[pHole[p][0]] = (unsigned char) pHole[p][1]-1;
		}
		
		H5Dclose(dataset);
		H5Sclose(dataspace);
	}
	else {
		hsize_t field_len;
		hid_t dataset = H5Dopen2(file_id, "PulseData/BaseCalls/ZMW/HoleNumber", H5P_DEFAULT); if (dataset < 0) return NULL;
		hid_t field_space = H5Dget_space(dataset); if (field_space < 0) return NULL;
		H5Sget_simple_extent_dims(field_space, &field_len, NULL);
		
		PB->nHoles = (unsigned int) field_len;
		H5Dclose(dataset);
		H5Sclose(field_space);
		
		unsigned char * const restrict ID = (unsigned char * const restrict) malloc(field_len*sizeof(unsigned char));
		if (ID == NULL) {
			PacBioClose(PB);
			return NULL;
		}
		
		memset(ID,0,field_len*sizeof(unsigned char));
		
		PB->HoleFileIndex = ID;
	}
		
	/* Allocate memory and open PacBio part files */
	{
		char FileName [FILENAME_SIZE+DIRECTORY_SIZE];
		const unsigned int nParts = PB->nParts;
		hid_t * const restrict Parts = (hid_t*) malloc(nParts*sizeof(hid_t*));
		if (Parts == NULL) {
			PacBioClose(PB);
			return NULL;
		}
		PB->Parts = Parts;
		memset(Parts, 0, sizeof(hid_t));
	
		if (isMultipart) {
			for (unsigned int p=0; p<nParts; p++) {
				sprintf(FileName, "%s/%s", PB->Directory, PB->PartsFileName[p]);
				PB->Parts[p] = H5Fopen(FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
				if (PB->Parts[p]<0 ) {
					PacBioClose(PB);
					return NULL;
				}
			}
		}
		else {
			PB->Parts[0] = file_id;
		}
	}
	
	/* Determine PacBio file type and return
	 *   - Trc have TraceData
	 *   - Pls have AnalysisConfig
	 *   - Bas have PulseData/BaseCalls
	 */
	H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client_data);
	H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
	{
		file_id = PB->Parts[0];
	  herr_t status = H5Gget_objinfo (file_id, "/TraceData", 0, NULL);
		if (status == 0) {
			PB->Content |= TRACE;
			if (isMultipart) {
				unsigned int * const HPI = (unsigned int * const restrict) malloc(PB->nHoles*sizeof(unsigned int));
				if (HPI == NULL) {
					PacBioClose(PB);
					return NULL;
				}
				PB->HoleFilePositionIndex = HPI;
				
				hsize_t field_len;
				for (unsigned int p=0; p<PB->nParts; p++) {
					hid_t dataset = H5Dopen2(PB->Parts[p], "TraceData/HoleNumber", H5P_DEFAULT); if (dataset < 0) return NULL;
					hid_t field_space = H5Dget_space(dataset); if (field_space < 0) return NULL;
					H5Sget_simple_extent_dims(field_space, &field_len, NULL);
					if (H5Dread(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp) < 0) {
						H5Dclose(dataset);
						H5Sclose(field_space);
						PacBioClose(PB);
						return NULL;
					}
					
					H5Sclose(field_space);
					H5Dclose(dataset);
					
					for (hsize_t p=0; p<field_len; p++) {
						HPI[tmp[p]] = p;
					}
				}
				free(tmp);
			}
			
			/* Get Trace length */
			hsize_t field_len[3];
			hid_t dataset = H5Dopen2(file_id, "TraceData/Traces", H5P_DEFAULT);
			hid_t field_space = H5Dget_space(dataset);
			H5Sget_simple_extent_dims(field_space, field_len, NULL);
			PB->Traces.Length = (size_t) field_len[2];
			H5Sclose(field_space);
			H5Dclose(dataset);
			
		}
		status = H5Gget_objinfo (file_id, "/AnalysisConfig", 0, NULL);
		if (status == 0) PB->Content |= PULSE;
		
		status = H5Gget_objinfo (file_id, "/PulseData/BaseCalls", 0, NULL);
		if (status == 0) PB->Content |= BASECALLING;
		
		status = H5Gget_objinfo (file_id, "/PulseData/ConsensusBaseCalls", 0, NULL);
		if (status == 0) PB->Content |= CONSENSUS_BASECALLING;
	}
	
	H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
	return PB;
}

int PacBioClose(PacBio_t * PB)
{
	if (PB->nParts > 1) H5Fclose(PB->Base);
	
	if (PB->PartsFileName) {
		for (unsigned int i=0; i<PB->nParts; i++) {
			if (PB->PartsFileName[i]) free(PB->PartsFileName[i]);
		}
	}
	
	if (PB->Parts) {
		for (unsigned int i=0; i<PB->nParts; i++) {
			if (PB->Parts[i]) H5Fclose(PB->Parts[i]);
		}
		free(PB->Parts);
	}
	
	if (PB->HoleFileIndex) free(PB->HoleFileIndex);
	
	if (PB->Status) free(PB->Status);
	if (PB->Coordinates) free(PB->Coordinates);
	
	if (PB->Content & BASECALLING) {
		if (PB->BaseCalling.Regions != NULL) {
			free(PB->BaseCalling.Regions);
		}
	}
	if (PB->Content & CONSENSUS_BASECALLING) {
		if (PB->ConsensusBaseCalling.Regions != NULL) {
			free(PB->ConsensusBaseCalling.Regions);
		}
	}
	
	if (PB->Content & TRACE) {
		if (PB->Traces.DecodeTable) free(PB->Traces.DecodeTable);
		if (PB->Traces.Bias) free(PB->Traces.Bias);
	}
	
	if (PB->HoleFilePositionIndex) free(PB->HoleFilePositionIndex);
	
	free(PB);
	PB = NULL;
	
	return SUCCESS;
}

/* Content contains BASECALLING */

int IndexRegions(PacBio_t * const restrict PB)
{
	if (!(PB->Content & BASECALLING)) return ERROR;
	
	hsize_t AllocatedSpace = 0;
	hsize_t lengthAllocatedSpace = 0;
	int * restrict Space  = NULL;
	int * restrict length = NULL;
	
	RegionIndex_t * const restrict Regions = (RegionIndex_t * const restrict) malloc(PB->nHoles*sizeof(RegionIndex_t));
	if (Regions == NULL) return ERROR;
	
	int HoleID = 0;
	for (unsigned int p=0; p<PB->nParts; p++) {
		/*
		 * Read Regions from file 
		 */
		
		/* Get dimensions */
		hsize_t dims[2];
		hid_t dataset = H5Dopen2(PB->Parts[p], "PulseData/Regions", H5P_DEFAULT); if (dataset < 0) return ERROR;
		hid_t dataspace = H5Dget_space(dataset); if (dataspace < 0) return ERROR;
		H5Sget_simple_extent_dims(dataspace, dims, NULL);
		hid_t plist = H5Dget_create_plist(dataset);
		hsize_t chunkdims[2];
		const int chunkdim = H5Pget_chunk(plist, 2, chunkdims); 
		if (dims[0] > AllocatedSpace) {
				Space = (int * restrict) realloc(Space, chunkdims[0]*chunkdims[1]*sizeof(int));
				if (Space == NULL) {
					return ERROR;
				}
				AllocatedSpace = chunkdims[0];
		}
		
		/* Get length of each regions */
		hsize_t lengthdims[2];
		hid_t lengthdataset = H5Dopen2(PB->Parts[p], "PulseData/BaseCalls/ZMW/NumEvent", H5P_DEFAULT); 
		if (lengthdataset < 0) return ERROR;
		hid_t lengthdataspace = H5Dget_space(lengthdataset); if (lengthdataspace < 0) return ERROR;
		H5Sget_simple_extent_dims(lengthdataspace, lengthdims, NULL);
		if (lengthdims[0] > lengthAllocatedSpace) {
				length = (int * restrict) realloc(length, lengthdims[0]*sizeof(int));
				if (length == NULL) {
					return ERROR;
				}
				lengthAllocatedSpace = lengthdims[0];
		}
		if (H5Dread(lengthdataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, length) < 0) return ERROR;
		
		hsize_t hyperslab_start[2] = {0,0};
		hsize_t hyperslab_count[2] = { chunkdims[0], chunkdims[1]};

// 		printf("Regions is: %lu x %lu\n", chunkdims[0], chunkdims[1]);
		
		hid_t mspace_id = H5Screate_simple(2, chunkdims, NULL);
		
		/* Read and parse data */
		hsize_t maxCount = dims[0] - chunkdims[0];
		unsigned int start = 0;

		hsize_t StartBase = 0;
		int LastBase = -1;
		int index = 0;

		while (hyperslab_start[0] < maxCount) {
//  			printf("Reading hyperslab %lu + %lu : start @ %lu\n", hyperslab_start[0], hyperslab_count[0], StartBase);
			H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
			if (H5Dread(dataset, H5T_NATIVE_INT, mspace_id, dataspace, H5P_DEFAULT, Space) < 0) return ERROR;
			
			for (hsize_t i=0; i<chunkdims[0]; i++) {
				if (Space[i*5] == HoleID) {
					LastBase = ( LastBase > Space[i*5+3]) ? LastBase : Space[i*5+3];
				}
				else {
					Regions[HoleID].start = start;
					Regions[HoleID].stop  = (unsigned int) (i+hyperslab_start[0]-1);
					
					Regions[HoleID].BaseCallsOffset = StartBase;
					// WARNING: we should maybe Pac for negative number here !!!
					Regions[HoleID].BaseCallsLength = length[index];
					HoleID = Space[i*5];
					start = (unsigned int) (i+hyperslab_start[0]);
					StartBase += (hsize_t) length[index++];
					LastBase = Space[i*5+3];
				}
			}
			
			hyperslab_start[0] 			 += chunkdims[0];
		}
		H5Sclose(mspace_id);
		hyperslab_count[0] = dims[0] - hyperslab_start[0];
		mspace_id = H5Screate_simple(2, hyperslab_count, NULL);
// 		printf("Reading hyperslab %lu + %lu\n", hyperslab_start[0], hyperslab_count[0]);
		H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
		if (H5Dread(dataset, H5T_NATIVE_INT, mspace_id, dataspace, H5P_DEFAULT, Space) < 0) return ERROR;
		for (hsize_t i=0; i<hyperslab_count[0]; i++) {
			if (Space[i*5] == HoleID) {
				LastBase = ( LastBase > Space[i*5+3]) ? LastBase : Space[i*5+3];
			}
			else {
				Regions[HoleID].start = start;
				Regions[HoleID].stop  = (unsigned int) (i+hyperslab_start[0]-1);
				
				Regions[HoleID].BaseCallsOffset = StartBase;
				HoleID = Space[i*5];
				start = (unsigned int) (i+hyperslab_start[0]);
				StartBase += (hsize_t) LastBase;
				LastBase = Space[i*5+3];
			}
		}
		Regions[HoleID].start = start;
		Regions[HoleID].stop  = (unsigned int) (hyperslab_count[0]+hyperslab_start[0]-1);
		Regions[HoleID].BaseCallsOffset = StartBase;
		
		++HoleID;
		
		H5Dclose(dataset);
		H5Sclose(dataspace);
		H5Sclose(mspace_id);
		H5Pclose(plist);
		H5Dclose(lengthdataset);
		H5Sclose(lengthdataspace);
	}
		
	free(Space);
	free(length);
	PB->BaseCalling.Regions = Regions;
	
	return SUCCESS;
}

int PopulateHoleStatus(PacBio_t * const restrict PB)
{
	unsigned char * const Status = (unsigned char *) malloc(PB->nHoles*sizeof(unsigned char));
	if (Status == NULL) return ERROR;
	
	size_t offset = 0;
	const unsigned int nParts = PB->nParts;
	const char * WherePtr = (PB->Content & TRACE) ? WhereHoleStatus[1] : WhereHoleStatus[0];
	
	for (unsigned int p=0; p<nParts; p++) {
		hsize_t dims[2];
		hid_t dataset = H5Dopen2(PB->Parts[p], WherePtr, H5P_DEFAULT); if (dataset < 0) return ERROR;
		hid_t dataspace = H5Dget_space(dataset); if (dataspace < 0) return ERROR;
		H5Sget_simple_extent_dims(dataspace, dims, NULL);
		
		if (H5Dread(dataset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Status[offset]) < 0) return ERROR;
		offset += dims[0];
		
		H5Dclose(dataset);
		H5Sclose(dataspace);
	}
	PB->Status = Status;
	return SUCCESS;
}

int PopulateHoleCoordinates(PacBio_t * const restrict PB)
{
	short int * const Coordinates = (short int *) malloc(PB->nHoles*2*sizeof(short int));
	if (Coordinates == NULL) return ERROR;
	
	size_t offset = 0;
	const unsigned int nParts = PB->nParts;
	const char * WherePtr = (PB->Content & TRACE) ? WhereHoleCoord[1] : WhereHoleCoord[0];
	for (unsigned int p=0; p<nParts; p++) {
		hsize_t dims[2];
		hid_t dataset = H5Dopen2(PB->Parts[p], WherePtr, H5P_DEFAULT); if (dataset < 0) return ERROR;
		hid_t dataspace = H5Dget_space(dataset); if (dataspace < 0) return ERROR;
		H5Sget_simple_extent_dims(dataspace, dims, NULL);
		
		if (H5Dread(dataset, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Coordinates[2*offset]) < 0) return ERROR;
		offset += dims[0];
		
		H5Dclose(dataset);
		H5Sclose(dataspace);
	}
	PB->Coordinates = (short int (*)[2]) Coordinates;
	return SUCCESS;
}

int getHoleRegions(const PacBio_t * const restrict PB, HoleData_t * const restrict HReg,
                   const unsigned int HoleNumber)
{
	if (!(PB->Content & BASECALLING)) return ERROR;
	
	if (PB->BaseCalling.Regions == NULL) {
		fputs("Please index regions first!\n",stderr);
		return -1;
	}
	
	const RegionIndex_t * const restrict Region = &PB->BaseCalling.Regions[HoleNumber];
	const size_t RegionSize = Region->stop - Region->start + 1;
	HoleRegion_t * restrict Regions;
	if (HReg->nAllocatedRegions < RegionSize) {
		Regions = (HoleRegion_t * const restrict) realloc(HReg->Regions,RegionSize*sizeof(HoleRegion_t));
		if (Regions == NULL) return -2;
	}
	
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
	if (H5Dread(dataset, H5T_NATIVE_INT, mspace_id, dataspace, H5P_DEFAULT, Regions) < 0) return ERROR;
	
	HReg->nAllocatedRegions = (unsigned int) RegionSize;
	HReg->nRegions          = (unsigned int) RegionSize;
	HReg->HoleNumber        = HoleNumber;
	HReg->Regions           = Regions;
	
	H5Sclose(mspace_id);
	H5Sclose(dataspace);
	H5Dclose(dataset);
	
	return (int) RegionSize;
}

int getHoleReads_HN(const PacBio_t * const restrict PB, HoleReads_t * HRead, const unsigned int HoleNumber)
{
	if (!(PB->Content & BASECALLING)) return ERROR;

	if (PB->BaseCalling.Regions == NULL) {
		fputs("Please index regions first!\n",stderr);
		return -1;
	}
	
	const RegionIndex_t * const restrict Region = &PB->BaseCalling.Regions[HoleNumber];
	const size_t RegionSize = Region->stop - Region->start + 1;
	
	HoleRegion_t * const restrict HReg = (HoleRegion_t * const restrict) malloc(RegionSize*sizeof(HoleRegion_t));
	if (HReg == NULL) return -2;
	
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
	if (H5Dread(dataset, H5T_NATIVE_INT, mspace_id, dataspace, H5P_DEFAULT, HReg) < 0) return ERROR;
	
	H5Sclose(mspace_id);
	H5Sclose(dataspace);
	H5Dclose(dataset);
	
	/* Get the sequence range */
	int min = 0, max = 0;
	unsigned int InsertCount = 0;
	for (int r=0; r<RegionSize; r++) {
		if (HReg[r].type == Insert) {
			const int lstart = HReg[r].start;
			min = (min < lstart) ? min : lstart;
			const int lstop = HReg[r].stop;
			max = (max > lstop) ? max : lstop;
			++InsertCount;
		}
	}
	
	/* Allocate memory for the Reads */
	const size_t msize = max - min + 1;
	HRead->Stream = (unsigned char*) malloc(msize*sizeof(unsigned char));
	if (HRead->Stream == NULL) {
		free(HReg);
		return -5;
	}

	/* Read sequence string */
	
	dataset = H5Dopen2(PB->Parts[FileNum], "PulseData/BaseCalls/Basecall", H5P_DEFAULT); if (dataset < 0) return ERROR;
	dataspace = H5Dget_space(dataset); if (dataspace < 0) return ERROR;
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	
	outdims[0] = msize;
	outdims[1] =  0;
	mspace_id = H5Screate_simple(1, outdims, NULL);
	
	hyperslab_start[0] = Region->BaseCallsOffset + min;
	hyperslab_start[1] = 1;
	hyperslab_count[0] = msize;
	hyperslab_count[1] = 1;
	
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
	if (H5Dread(dataset, H5T_NATIVE_UCHAR, mspace_id, dataspace, H5P_DEFAULT, HRead->Stream) < 0) return ERROR;
	
	H5Sclose(mspace_id);
	H5Sclose(dataspace);
	H5Dclose(dataset);
	
	free(HReg);
	
	return (int) InsertCount;
}

int getHoleQVs_HN(const PacBio_t * const restrict PB, HoleQVs_t * HRead, const unsigned int HoleNumber)
{
	if (!(PB->Content & BASECALLING)) return ERROR;
	if (PB->BaseCalling.Regions == NULL) {
		fputs("Please index regions first!\n",stderr);
		return -1;
	}
	
	const RegionIndex_t * const restrict Region = &PB->BaseCalling.Regions[HoleNumber];
	const size_t RegionSize = Region->stop - Region->start + 1;
	
	HoleRegion_t * const restrict HReg = (HoleRegion_t * const restrict) malloc(RegionSize*sizeof(HoleRegion_t));
	if (HReg == NULL) return -2;
	
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
	if (H5Dread(dataset, H5T_NATIVE_INT, mspace_id, dataspace, H5P_DEFAULT, HReg) < 0) return ERROR;
	
	H5Sclose(mspace_id);
	H5Sclose(dataspace);
	H5Dclose(dataset);
		
	/* Get the sequence range */
	int min = 0, max = 0;
	unsigned int InsertCount = 0;
	for (int r=0; r<RegionSize; r++) {
		if (HReg[r].type == Insert) {
			const int lstart = HReg[r].start;
			min = (min < lstart) ? min : lstart;
			const int lstop = HReg[r].stop;
			max = (max > lstop) ? max : lstop;
			++InsertCount;
		}
	}
	
	/* Allocate memory for the Reads */
	const size_t msize = max - min + 1;
	HRead->Stream = (unsigned char*) malloc(msize*sizeof(unsigned char));
	if (HRead->Stream == NULL) {
		free(HReg);
		return -5;
	}
		
	/* Read sequence string */
	
	dataset = H5Dopen2(PB->Parts[FileNum], "PulseData/BaseCalls/QualityValue", H5P_DEFAULT); if (dataset < 0) return ERROR;
	dataspace = H5Dget_space(dataset); if (dataspace < 0) return ERROR;
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	
	outdims[0] = msize;
	outdims[1] =  0;
	mspace_id = H5Screate_simple(1, outdims, NULL);
	
	hyperslab_start[0] = Region->BaseCallsOffset + min;
	hyperslab_start[1] = 1;
	hyperslab_count[0] = msize;
	hyperslab_count[1] = 1;
	
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
	if (H5Dread(dataset, H5T_NATIVE_UCHAR, mspace_id, dataspace, H5P_DEFAULT, HRead->Stream) < 0) return ERROR;
	
	H5Sclose(mspace_id);
	H5Sclose(dataspace);
	H5Dclose(dataset);
	
	free(HReg);
	
	return (int) InsertCount;
}

int getHoleReads_HR(const PacBio_t * const restrict PB, HoleReads_t * HRead, const HoleData_t * const restrict HReg)
{
	if (!(PB->Content & BASECALLING)) return ERROR;
	if (PB->BaseCalling.Regions == NULL) {
		fputs("Please index regions first!\n",stderr);
		return -1;
	}

	const RegionIndex_t * const restrict Region = &PB->BaseCalling.Regions[HReg->HoleNumber];
	const size_t RegionSize = Region->stop - Region->start + 1;
	
	/* Get the sequence range */
	int min = 0, max = 0;
	unsigned int InsertCount = 0;
	for (int r=0; r<RegionSize; r++) {
		if (HReg->Regions[r].type == Insert) {
			const int lstart = HReg->Regions[r].start;
			min = (min < lstart) ? min : lstart;
			const int lstop = HReg->Regions[r].stop;
			max = (max > lstop) ? max : lstop;
			++InsertCount;
		}
	}
	
	/* Allocate memory for the Reads */
	const size_t msize = max - min + 1;
	HRead->Stream = (unsigned char*) malloc(msize*sizeof(unsigned char));
	if (HRead->Stream == NULL) {
		return -5;
	}
	
	/* Read sequence string */
	hsize_t dims[2];
	const size_t FileNum = PB->HoleFileIndex[HReg->HoleNumber];
	hid_t dataset = H5Dopen2(PB->Parts[FileNum], "PulseData/BaseCalls/Basecall", H5P_DEFAULT); if (dataset < 0) return ERROR;
	hid_t dataspace = H5Dget_space(dataset); if (dataspace < 0) return ERROR;
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	
	hsize_t outdims[2] = {msize, 0 };
	hid_t mspace_id = H5Screate_simple(1, outdims, NULL);
	
	hsize_t hyperslab_start[2] = { Region->BaseCallsOffset + min, 1};
	hsize_t hyperslab_count[2] = { msize, 1};
	
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
	if (H5Dread(dataset, H5T_NATIVE_UCHAR, mspace_id, dataspace, H5P_DEFAULT, HRead->Stream) < 0) return ERROR;
	
	H5Sclose(mspace_id);
	H5Sclose(dataspace);
	H5Dclose(dataset);
	
	return (int) InsertCount;
}

int getHoleQVs_HR(const PacBio_t * const restrict PB, HoleQVs_t * HRead, const HoleData_t * const restrict HReg)
{
	if (!(PB->Content & BASECALLING)) return ERROR;
		if (PB->BaseCalling.Regions == NULL) {
		fputs("Please index regions first!\n",stderr);
		return -1;
	}
	
	const RegionIndex_t * const restrict Region = &PB->BaseCalling.Regions[HReg->HoleNumber];
	const size_t RegionSize = Region->stop - Region->start + 1;
		
	/* Get the sequence range */
	int min = 0, max = 0;
	unsigned int InsertCount = 0;
	for (int r=0; r<RegionSize; r++) {
		if (HReg->Regions[r].type == Insert) {
			const int lstart = HReg->Regions[r].start;
			min = (min < lstart) ? min : lstart;
			const int lstop = HReg->Regions[r].stop;
			max = (max > lstop) ? max : lstop;
			++InsertCount;
		}
	}
	
	/* Allocate memory for the Reads */
	const size_t msize = max - min + 1;
	HRead->Stream = (unsigned char*) malloc(msize*sizeof(unsigned char));
	if (HRead->Stream == NULL) {
		return -5;
	}
	
	/* Read sequence string */
	hsize_t dims[2];
	const size_t FileNum = PB->HoleFileIndex[HReg->HoleNumber];
	hid_t dataset = H5Dopen2(PB->Parts[FileNum], "PulseData/BaseCalls/QualityValue", H5P_DEFAULT); if (dataset < 0) return ERROR;
	hid_t dataspace = H5Dget_space(dataset); if (dataspace < 0) return ERROR;
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	
	hsize_t outdims = msize;
	hid_t mspace_id = H5Screate_simple(1, &outdims, NULL);
	
	hsize_t hyperslab_start[2] = { Region->BaseCallsOffset + min, 1};
	hsize_t hyperslab_count[2] = { msize, 1};
	
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
	if (H5Dread(dataset, H5T_NATIVE_UCHAR, mspace_id, dataspace, H5P_DEFAULT, HRead->Stream) < 0) return ERROR;
	
	H5Sclose(mspace_id);
	H5Sclose(dataspace);
	H5Dclose(dataset);
	
	return (int) InsertCount;
}

int getSequencingHoles(const PacBio_t * const restrict PB, unsigned int ** SequencingHoles)
{
	if (PB->Status == NULL) return ERROR;
	
	int count = 0;
	const unsigned int nHoles = PB->nHoles;
	const unsigned char * const restrict Status = PB->Status; 
	
	for (unsigned int hole=0; hole<nHoles; hole++) {
		if (Status[hole] == (unsigned char) SEQUENCING) ++count;
	}
	
	unsigned int * restrict res = (unsigned int * const) malloc(count*sizeof(unsigned int));
	if (res == NULL) return ERROR;
	*SequencingHoles = res;
	
	for (unsigned int hole=0; hole<nHoles; hole++) {
		if (Status[hole] == (unsigned char) SEQUENCING) *res++ = hole;
	}
	
	return count;
}

int getBasecallsPulses(const PacBio_t * const restrict PB, BaseCallsPulses_t * const restrict BP, const unsigned int HoleNumber)
{
	if (!(PB->Content & BASECALLING)) goto END;
	if (PB->BaseCalling.Regions == NULL) {
		fputs("Please index regions first!\n",stderr);
		goto END;
	}
	
	const hid_t file_id  = PB->Parts[PB->HoleFileIndex[HoleNumber]]; 
	
	hsize_t hyperslab_start = PB->BaseCalling.Regions[HoleNumber].BaseCallsOffset;;
	hsize_t hyperslab_count = (hsize_t) PB->BaseCalling.Regions[HoleNumber].BaseCallsLength;
	hid_t mspace_id = H5Screate_simple(1, &hyperslab_count, NULL);
	
	/* Read PreBaseFrame */
	unsigned short int * const restrict PreBaseFrames = malloc(hyperslab_count*sizeof(unsigned short int));
	if (PreBaseFrames == NULL) goto END1;
	
	hid_t dataset = H5Dopen2(file_id, "PulseData/BaseCalls/PreBaseFrames", H5P_DEFAULT); if (dataset < 0 ) goto END2;
	hid_t dataspace = H5Dget_space(dataset); if (dataspace < 0) goto CLEAN_2;
	
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &hyperslab_start, NULL, &hyperslab_count, NULL);
	if (H5Dread(dataset, H5T_NATIVE_USHORT, mspace_id, dataspace, H5P_DEFAULT, PreBaseFrames) < 0) goto CLEAN_1;
	
	H5Sclose(dataspace);
	H5Dclose(dataset);
	
	/* Read WidthInFrame */
	unsigned short int * const restrict WidthInFrames = malloc(hyperslab_count*sizeof(unsigned short int));
	if (WidthInFrames == NULL) goto END2;
	
	dataset = H5Dopen2(file_id, "PulseData/BaseCalls/WidthInFrames", H5P_DEFAULT); if (dataset < 0 ) goto END1;
	dataspace = H5Dget_space(dataset); if (dataspace < 0) goto CLEAN_2;
	
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &hyperslab_start, NULL, &hyperslab_count, NULL);
	if (H5Dread(dataset, H5T_NATIVE_USHORT, mspace_id, dataspace, H5P_DEFAULT, WidthInFrames) < 0) goto CLEAN_1;
	
	H5Sclose(dataspace);
	H5Dclose(dataset);
	
	/* Read BaseCalling */
	char * const restrict BaseCalls = malloc(hyperslab_count*sizeof(char));
	if (BaseCalls == NULL) goto END2;
	
	dataset = H5Dopen2(file_id, "PulseData/BaseCalls/Basecall", H5P_DEFAULT); if (dataset < 0 ) goto END1;
	dataspace = H5Dget_space(dataset); if (dataspace < 0) goto CLEAN_2;
	
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &hyperslab_start, NULL, &hyperslab_count, NULL);
	if (H5Dread(dataset, H5T_NATIVE_UCHAR, mspace_id, dataspace, H5P_DEFAULT, BaseCalls) < 0) goto CLEAN_1;
	
	H5Sclose(dataspace);
	H5Dclose(dataset);
	H5Sclose(mspace_id);
	
	/* Compute frames start and end */
	unsigned int * const restrict FramesStart = malloc(hyperslab_count*sizeof(unsigned int));
	unsigned int * const restrict FramesEnd   = malloc(hyperslab_count*sizeof(unsigned int));
	if (FramesStart == NULL || FramesEnd == NULL) goto END2;
	
	register unsigned int Sum = 0U;
	for (hsize_t Events=0; Events<hyperslab_count; Events++) {
		const register unsigned int utmp = PreBaseFrames[Events] + Sum;
		FramesStart[Events] = utmp;
		Sum                 = utmp + WidthInFrames[Events];
		FramesEnd[Events]   = Sum;
	}
	
	BP->FramesStart = FramesStart;
	BP->FramesEnd   = FramesEnd;
	BP->Bases       = BaseCalls;
	BP->nBases      = (size_t) hyperslab_count;
	
	free(PreBaseFrames);
	free(WidthInFrames);
	
	return SUCCESS;
	
	CLEAN_1: ;
		H5Sclose(dataspace);
	CLEAN_2: ;
		H5Dclose(dataset);
	
	END1:
		H5Sclose(mspace_id);
		
	END2:
		if (PreBaseFrames) free(PreBaseFrames);
		if (WidthInFrames) free(WidthInFrames);
		if (BaseCalls) free(BaseCalls);
		if (FramesStart) free(FramesStart);
		if (FramesEnd) free(FramesEnd);

	END:
		return ERROR;
}


/* Content contains CONSENSUS */
int IndexConsensus(PacBio_t * const restrict PB)
{
	if (!(PB->Content & CONSENSUS_BASECALLING)) return ERROR;
	const unsigned int nHoles = PB->nHoles;
	
	RegionIndex_t * const restrict Regions = (RegionIndex_t * const restrict) malloc(PB->nHoles*sizeof(RegionIndex_t));
	if (Regions == NULL) return ERROR;
	
	if (PB->nParts != 1 ) return ERROR;
	
	hid_t file_id = PB->Parts[0];
	int HoleID = 0;

	/* Get length of each regions */
	hsize_t lengthdims;
	hid_t lengthdataset = H5Dopen2(file_id, "PulseData/ConsensusBaseCalls/ZMW/NumEvent", H5P_DEFAULT); 
	if (lengthdataset < 0) return ERROR;
	hid_t lengthdataspace = H5Dget_space(lengthdataset); if (lengthdataspace < 0) return ERROR;
	H5Sget_simple_extent_dims(lengthdataspace, &lengthdims, NULL);
	
	if (lengthdims != nHoles) return ERROR;
	
	int * const restrict length = (int * const restrict) malloc(nHoles*sizeof(int));
	if (length == NULL) return ERROR;
	
	if (H5Dread(lengthdataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, length) < 0) return ERROR;
	
	H5Dclose(lengthdataset);
	H5Sclose(lengthdataspace);
	
	hsize_t Offset = 0;
	for (unsigned int iHole=0; iHole<nHoles; iHole++) {
		Regions[iHole].BaseCallsOffset = Offset;
		Regions[iHole].start = 0U;
		Regions[iHole].stop = 0U;
		Offset += (hsize_t) length[iHole];
	}
	
	free(length);
	PB->ConsensusBaseCalling.Regions = Regions;
	return SUCCESS;
}

/* Content contains TRACES */
int PopulateDecodeTable(PacBio_t * const restrict PB)
{
	if (!(PB->Content & TRACE)) goto END;
	
	const unsigned int nParts = PB->nParts;
	float (*DecodeTable)[256] = (float (*)[256]) malloc(nParts*256*sizeof(float));
	if (DecodeTable == NULL) goto END;
	
	float * Bias = (float*) malloc(nParts*sizeof(float));
	if (Bias == NULL) goto FREE_DECODE_TABLE;
	
	for (unsigned int p=0; p<nParts; p++) {
// 		hsize_t dims[2];
		hid_t dataset = H5Dopen2(PB->Parts[p], "TraceData/Codec/Decode", H5P_DEFAULT);
		if (dataset < 0) goto FREE_DECODE_TABLE;
// 		hid_t dataspace = H5Dget_space(dataset); if (dataspace < 0) goto FREE_DECODE_TABLE;
// 		H5Sget_simple_extent_dims(dataspace, dims, NULL);
		
		if (H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, DecodeTable[p]) < 0)
			goto FREE_DECODE_TABLE;
		
		H5Dclose(dataset);
// 		H5Sclose(dataspace);
		
		dataset = H5Gopen2(PB->Parts[p], "TraceData/Codec", H5P_DEFAULT);
		if (dataset < 0) goto FREE_DECODE_TABLE;
		hid_t attr = H5Aopen(dataset, "Bias", H5P_DEFAULT);
	  if (attr < 0) goto FREE_DECODE_TABLE;
		if (H5Aread(attr, H5T_NATIVE_FLOAT, &Bias[p]) < 0) goto FREE_DECODE_TABLE;
		H5Aclose(attr);
		H5Gclose(dataset);
	}
	
	PB->Traces.DecodeTable = DecodeTable;
	PB->Traces.Bias = Bias;
	return SUCCESS;
	
	FREE_DECODE_TABLE:
		free(DecodeTable);
	FREE_BIAS:
		free(Bias);
	END:
		return ERROR;
}


Traces_t * getTraceIndex(const PacBio_t * const restrict PB, const unsigned int HoleNumber)
{
	if (!(PB->Content & TRACE) || PB->Traces.DecodeTable == NULL || PB->Traces.Bias == NULL) goto END;
	
	const size_t FileNum = PB->HoleFileIndex[HoleNumber];
	hid_t dataset = H5Dopen2(PB->Parts[FileNum], "TraceData/Traces", H5P_DEFAULT);
	
	Traces_t * const Results = (Traces_t*) malloc(sizeof(Traces_t));
	if (Results == NULL) goto CLEAN;
	
	hsize_t dims[3];
	hid_t dataspace = H5Dget_space(dataset); if (dataspace < 0) goto CLEAN;
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	
	const unsigned int stride = (dims[2] + 15) & ~(15);
	Results->stride = stride;
	unsigned char * const restrict ctmp = _mm_malloc(4*stride*sizeof(unsigned char), 16);
	if (ctmp == NULL) goto CLEAN;
	
	hsize_t start = (hsize_t) PB->HoleFilePositionIndex[HoleNumber];
	

	hsize_t outdims[3] = {1, 1, dims[2]};
	hid_t mspace_id = H5Screate_simple(3, outdims, NULL);
	
	hsize_t hyperslab_start[3]  = {start,0,0};
	hsize_t hyperslab_count[3]  = {1, 1, dims[2]};
	Results->T = ctmp;
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
	if (H5Dread(dataset, H5T_NATIVE_UCHAR, mspace_id, dataspace, H5P_DEFAULT, Results->T) < 0) goto FREE1;
	
	Results->G = ctmp + stride;
	hyperslab_start[1] = 1;
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
	if (H5Dread(dataset, H5T_NATIVE_UCHAR, mspace_id, dataspace, H5P_DEFAULT, Results->G) < 0) goto FREE1;
	
	Results->A = ctmp + 2*stride;
	hyperslab_start[1] = 2;
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
	if (H5Dread(dataset, H5T_NATIVE_UCHAR, mspace_id, dataspace, H5P_DEFAULT, Results->A) < 0) goto FREE1;
	
	Results->C = ctmp + 3*stride;
	hyperslab_start[1] = 3;
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperslab_start, NULL, hyperslab_count, NULL);
	if (H5Dread(dataset, H5T_NATIVE_UCHAR, mspace_id, dataspace, H5P_DEFAULT, Results->C) < 0) goto FREE1;
	
	H5Dclose(dataset);
	H5Sclose(dataspace);
	H5Sclose(mspace_id);
	
	
	return Results;
	
	FREE1:
		H5Sclose(mspace_id);
		free(ctmp);
		H5Sclose(dataspace);
	CLEAN:
		H5Dclose(dataset);
		free(Results);
	END:
		return NULL;
}

__m128* getTrace(const PacBio_t * const restrict PB, const unsigned int HoleNumber)
{
	Traces_t * const Indices = getTraceIndex(PB, HoleNumber);
	if (Indices == NULL) goto END;
	
	const size_t stride = PB->Traces.Length;
	__m128 * const restrict Trace = (__m128*) _mm_malloc(stride*sizeof(__m128), 16);
	if (Trace == NULL) goto CLEAN;
	
	const size_t FileNum = PB->HoleFileIndex[HoleNumber];
	const __m128 __Bias = _mm_load1_ps(&(PB->Traces.Bias[FileNum]));
	const float * const restrict DecodeTable = PB->Traces.DecodeTable[FileNum];
 
	for (size_t i=0; i<stride; i++) {
		float ftmp[4] __attribute__((aligned(16)));
		
		ftmp[0] = DecodeTable[Indices->T[i]];
		ftmp[1] = DecodeTable[Indices->G[i]];
		ftmp[2] = DecodeTable[Indices->A[i]];
		ftmp[3] = DecodeTable[Indices->C[i]];
			
		_mm_store_ps((float*) &Trace[i], _mm_sub_ps(_mm_load_ps(&ftmp[0]), __Bias));
	}
	
	freeTraces(Indices);
	return Trace;
	
	CLEAN:
		freeTraces(Indices);
	END:
		return NULL;
}

__m128* DecodeTracePositionIndex(const PacBio_t * const restrict PB, Traces_t * const restrict Indices,
                                 const unsigned int HoleNumber)
{
	if (!(PB->Content & TRACE) || PB->Traces.DecodeTable == NULL || PB->Traces.Bias == NULL) goto END;
	
	const size_t Length = PB->Traces.Length;
	__m128 * const restrict Trace = (__m128*) _mm_malloc(Length*sizeof(__m128), 16);
	if (Trace == NULL) goto END;
	
	const size_t FileNum = PB->HoleFileIndex[HoleNumber];
	const __m128 __Bias = _mm_load1_ps(&(PB->Traces.Bias[FileNum]));
	const float * const restrict DecodeTable = PB->Traces.DecodeTable[FileNum];

	for (size_t i=0; i<Length; i++) {
		float ftmp[4] __attribute__((aligned(16)));
		
		ftmp[0] = DecodeTable[Indices->T[i]];
		ftmp[1] = DecodeTable[Indices->G[i]];
		ftmp[2] = DecodeTable[Indices->A[i]];
		ftmp[3] = DecodeTable[Indices->A[i]];
			
		_mm_store_ps((float*) &Trace[i], _mm_sub_ps(_mm_load_ps(&ftmp[0]), __Bias));
	}
	
	return Trace;
	
	END:
		return NULL;
}

int getHDFSummary(PacBio_t * const restrict PB, const unsigned int HQThreshold,
                  ZMWSummary_t * const restrict summary, FILE* restrict out)
{
	IndexRegions(PB);
	if (PB->Coordinates == NULL && out != NULL) {
		if (PopulateHoleCoordinates(PB) != SUCCESS) {
			fprintf(stderr,"Hole coordinates not populated...\n");
			return -1;
		}
	}
	
	if (PopulateHoleStatus(PB) != SUCCESS) {
		fprintf(stderr,"Hole status not populated...\n");
		return -1;
	}
	
	memset (summary, 0, sizeof(ZMWSummary_t));
	
	const unsigned int nHoles = PB->nHoles;
	const unsigned char * const restrict Status = PB->Status; 
	if (out == NULL) {
		for(unsigned int iHole=0; iHole< nHoles; iHole++) {
			HoleData_t HReg = HOLE_DATA_INIT;
			int HQValue = -1;
			unsigned int isValid = 0;
			if (Status[iHole] == (unsigned char) SEQUENCING) {
				summary->nSequencing++;
				const int count = getHoleRegions(PB, &HReg, iHole);
				if (count>0) {
					HoleRegion_t * restrict Regions = HReg.Regions;
					isValid = 1;
					
					/* Look for invalid HQRegion == invalid hole*/
					int k=count;
					while(--k >= 0) {
						if (Regions[k].type == HQRegion) {
							HQValue = Regions[k].quality;
							break;
						}
					}
					if (HQValue < 0) {
						fprintf(stderr, "Hole %u has invalid or no HQ region!!!\n",iHole);
					}
					else if (HQValue >= HQThreshold) {
						summary->nHQaboveThreshold++;
						restrictReadstoHQ(&HReg);
						for (int i=0; i<count; i++) {
							switch(Regions[i].type) {
								case (Insert) : summary->nWithinHQSubreads++; break;
								case (Useless): summary->nOutofHQSubreads++; break;
							}
						}
					}
					else {
						const int HasHQ = (Regions[k].stop != 0);
						if (HasHQ) {
							summary->nBelowHQThreshold++;
							restrictReadstoHQ(&HReg);
							for (int i=0; i<count; i++) {
								switch(Regions[i].type) {
									case (Insert) : summary->nBelowHQThresholdWithin++; break;
									case (Useless): summary->nBelowHQThresholdOutof++; break;
								}
							}
						}
						else {
							summary->nInvalidHQRegion++;
							isValid = 0;
							for (int i=0; i<count; i++) {
								if (Regions[i].type == Insert) {
									summary->nInvalidHQRegionSubreads++;
								}
							}
						}
					}
				}
				else {
					fprintf(stderr,"Hole %u appears to have no regions!!!\n", iHole);
				}
			}
		}
	}
	else {
		const char *QualityText[] = { "LOWQUAL", "HIGHQUAL", "INVALID" };
		short int (*Coordinates)[2] = PB->Coordinates;
		
		for(unsigned int iHole=0; iHole< nHoles; iHole++) {
			HoleData_t HReg = HOLE_DATA_INIT;
			int HQValue = -1;
			unsigned int isValid = 0;
			const char * restrict QualPtr;
			if (Status[iHole] == (unsigned char) SEQUENCING) {
				unsigned int nSubreads = 0U;
				summary->nSequencing++;
				const int count = getHoleRegions(PB, &HReg, iHole);
				if (count>0) {
					HoleRegion_t * restrict Regions = HReg.Regions;
					isValid = 1;
					
					/* Look for invalid HQRegion == invalid hole*/
					int k=count;
					while(--k >= 0) {
						if (Regions[k].type == HQRegion) {
							HQValue = Regions[k].quality;
							break;
						}
					}
					if (HQValue < 0) {
						fprintf(stderr, "Hole %u has invalid or no HQ region!!!\n",iHole);
					}
					else if (HQValue >= HQThreshold) {
						summary->nHQaboveThreshold++;
						restrictReadstoHQ(&HReg);
						QualPtr = QualityText[1];
						for (int i=0; i<count; i++) {
							switch(Regions[i].type) {
								case (Insert) : summary->nWithinHQSubreads++; nSubreads++; break;
								case (Useless): summary->nOutofHQSubreads++; nSubreads++;  break;
							}
						}
					}
					else {
						const int HasHQ = (Regions[k].stop != 0);
						if (HasHQ) {
							summary->nBelowHQThreshold++;
							restrictReadstoHQ(&HReg);
							QualPtr = QualityText[0];
							for (int i=0; i<count; i++) {
								switch(Regions[i].type) {
									case (Insert) : summary->nBelowHQThresholdWithin++; nSubreads++; break;
									case (Useless): summary->nBelowHQThresholdOutof++; nSubreads++; break;
								}
							}
						}
						else {
							summary->nInvalidHQRegion++;
							isValid = 0;
							QualPtr = QualityText[2];
							for (int i=0; i<count; i++) {
								if (Regions[i].type == Insert) {
									summary->nInvalidHQRegionSubreads++;
								}
							}
						}
					}
					unsigned int TotalLength = 0U;
					for (int i=0; i<count; i++) {
						if (Regions[i].stop > TotalLength) TotalLength = Regions[i].stop;
					}
					fprintf(out, "%u\t%hi\t%hi\t%s\t%u\t%u\t%u\n",
			          iHole, Coordinates[iHole][0], Coordinates[iHole][1], QualPtr,
						    HQValue, TotalLength, nSubreads);
				}
				else {
					fprintf(stderr,"Hole %u appears to have no regions!!!\n", iHole);
				}
			}
		}
	}
	
	
	return 0;
}
