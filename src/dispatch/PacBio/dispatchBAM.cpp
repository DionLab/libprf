#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sys/time.h>

#include "pbbam/BamFile.h"
#include "pbbam/BamReader.h"
#include "pbbam/virtual/ZmwReadStitcher.h"
#include "pbbam/PbiFilterTypes.h"
#include "pbbam/virtual/VirtualRegionType.h"
#include "pbbam/virtual/VirtualRegionTypeMap.h"

#define restrict 
#include "pb_bam.h"

extern "C" {
#include "pfProfile.h"
#include "pfOutput.h"
#include "pb_common.h"
#include "pfDispatch.h"
}

using namespace PacBio::BAM;

struct MyZmwRecord {
	std::string Sequence;
	int Accuracy;
	unsigned int HoleNumber;
	ZmwType Type;
	std::map<VirtualRegionType, std::vector<VirtualRegion> > RegionsMap;
	
	bool HasVirtualRegionType(const VirtualRegionType regionType) const {
		return RegionsMap.find(regionType) != RegionsMap.end(); 
	}
};


/// \brief Appends content of src vector to dst vector using move semantics.
///
/// \param[in]     src  Input vector that will be empty after execution
/// \param[in,out] dst  Output vector that will be appended to
///
template <typename T>
inline void MoveAppend(std::vector<T>& src, std::vector<T>& dst) noexcept
{
    if (dst.empty())
    {
        dst = std::move(src);
    }
    else
    {
        dst.reserve(dst.size() + src.size());
        std::move(src.begin(), src.end(), std::back_inserter(dst));
        src.clear();
    }
}

/// \brief Appends content of src vector to dst vector using move semantics.
///
/// \param[in]     src  Input vector via perfect forwarding
/// \param[in,out] dst  Output vector that will be appended to
///
template <typename T>
inline void MoveAppend(std::vector<T>&& src, std::vector<T>& dst) noexcept
{
    if (dst.empty())
    {
        dst = std::move(src);
    }
    else
    {
        dst.reserve(dst.size() + src.size());
        std::move(src.begin(), src.end(), std::back_inserter(dst));
        src.clear();
    }
}


static
std::string ToString(const RecordType type)
{
    static const auto lookup = std::map<RecordType, std::string>
    {
        { RecordType::ZMW,        "ZMW" },
        { RecordType::HQREGION,   "HQREGION" },
        { RecordType::SUBREAD,    "SUBREAD" },
        { RecordType::CCS,        "CCS" },
        { RecordType::SCRAP,      "SCRAP" },
        { RecordType::UNKNOWN,    "UNKNOWN" }
    };

    try {
        return lookup.at(type);
    } catch (std::exception&) {
        throw std::runtime_error("error: unknown RecordType encountered");
    }
}

static
std::string RegionToString(const VirtualRegionType type)
{
    static const auto lookup = std::map<VirtualRegionType, std::string>
    {
			{ VirtualRegionType::ADAPTER,		"ADAPTER" },
			{ VirtualRegionType::BARCODE,		"BARCODE" },
			{ VirtualRegionType::FILTERED,	"FILTERED"}, 
			{ VirtualRegionType::SUBREAD,		"SUBREAD" },
			{ VirtualRegionType::HQREGION,	"HQREGION"},
			{ VirtualRegionType::LQREGION,	"LQREGION"}, 
    };

    try {
        return lookup.at(type);
    } catch (std::exception&) {
        throw std::runtime_error("error: unknown RecordType encountered");
    }
}

static
void StitchSources(const std::vector<BamRecord>& sources, MyZmwRecord& record)
{
    const auto& firstRecord = sources[0];
    const auto& lastRecord = sources[sources.size() - 1];
		bool HasReadAccuracy = false;
		bool HasScrapZmwType = false;
    std::string& sequence = record.Sequence;
		std::map<VirtualRegionType, std::vector<VirtualRegion>>& virtualRegionsMap = record.RegionsMap;
		
    // initialize capacity
    const auto stitchedSize = lastRecord.QueryEnd() - firstRecord.QueryStart();
		sequence.clear();
    sequence.reserve(stitchedSize);
		virtualRegionsMap.clear();
		
		record.HoleNumber = firstRecord.HoleNumber();
		
//  		fprintf(stderr, "ZMW %u has %u bases:\n", record.HoleNumber, stitchedSize);
    // Stitch using tmp vars
    for(auto& b : sources)
    {
        sequence.append(b.Sequence());

//         MoveAppend(b.Qualities(), qualities);

       if (b.HasScrapRegionType())
        {
            const VirtualRegionType regionType = b.ScrapRegionType();
            virtualRegionsMap[regionType].emplace_back(regionType, b.QueryStart(), b.QueryEnd());
//  						fprintf(stderr,"\t%s\t%u\t%u",RegionToString( b.ScrapRegionType()).c_str(), b.QueryStart(), b.QueryEnd());
        }
//          else 
//  					fprintf(stderr,"\t%s\t%u\t%u",ToString( b.Type()).c_str(), b.QueryStart(), b.QueryEnd());

        if (b.HasLocalContextFlags())
        {
            std::pair<int, int> barcodes{-1, -1};
            if (b.HasBarcodes())
                barcodes = b.Barcodes();

            constexpr auto regionType = VirtualRegionType::SUBREAD;
            virtualRegionsMap[regionType].emplace_back(regionType, b.QueryStart(), b.QueryEnd(), b.LocalContextFlags(),
                barcodes.first, barcodes.second);
        }

        if (b.HasReadAccuracy() && !HasReadAccuracy) {
					record.Accuracy = (int) floorf(1000.0f*static_cast<float>(b.ReadAccuracy()));
					HasReadAccuracy = true;
				}

//  				if (b.HasReadAccuracy())
//  					fprintf(stderr, "\t%lf\n", static_cast<float>(b.ReadAccuracy()));
//  				else 
//  					fputc('\n', stderr);
					
        if (b.HasScrapZmwType())
        {
            if (!HasScrapZmwType)
                record.Type = b.ScrapZmwType();
            else if (record.Type != b.ScrapZmwType())
                throw std::runtime_error("ScrapZmwTypes do not match");
					
        }
        
    }

    assert(sequence.size() == stitchedSize);
    
    // Determine HQREGION bases on LQREGIONS
    if (virtualRegionsMap.find(VirtualRegionType::LQREGION) != virtualRegionsMap.end())
    {
        if (virtualRegionsMap[VirtualRegionType::LQREGION].size() == 1)
        {
            const auto lq = virtualRegionsMap[VirtualRegionType::LQREGION][0];
            if (lq.beginPos == 0)
                virtualRegionsMap[VirtualRegionType::HQREGION].emplace_back(
                    VirtualRegionType::HQREGION, lq.endPos, sequence.size());
            else if (lq.endPos == stitchedSize)
                virtualRegionsMap[VirtualRegionType::HQREGION].emplace_back(
                    VirtualRegionType::HQREGION, 0, lq.beginPos);
            else {
							throw std::runtime_error("Unknown HQREGION for ZMW " + std::to_string( record.HoleNumber));
						}
        }
        else
        {
            int beginPos = 0;
            for (const auto& lqregion : virtualRegionsMap[VirtualRegionType::LQREGION])
            {
                if (lqregion.beginPos - beginPos > 0)
                    virtualRegionsMap[VirtualRegionType::HQREGION].emplace_back(
                        VirtualRegionType::HQREGION, beginPos, lqregion.beginPos);
                beginPos = lqregion.endPos;
            }
        }
    }
    else
    {
        virtualRegionsMap[VirtualRegionType::HQREGION].emplace_back(
            VirtualRegionType::HQREGION, 0, sequence.size());
    }
}

static int
thpool_add_subreads(threadpool_t * const restrict thpool_p, ZmwReadStitcher& Reader, const enum PacBioSelection selection)
{
	MyZmwRecord record;
	pb_job_t * restrict newJob;
	
	record.RegionsMap[VirtualRegionType::ADAPTER] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::SUBREAD] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::HQREGION] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::LQREGION] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::BARCODE] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::FILTERED] = std::vector<VirtualRegion>();
	
	while  (Reader.HasNext()) {
		std::vector<BamRecord> records = Reader.NextRaw();
		// Sort sources by queryStart
    std::sort(records.begin(), records.end(),
              [](const BamRecord& l1, const BamRecord& l2)
              { return l1.QueryStart() < l2.QueryStart(); });
		try {
			StitchSources(records, record);
		} catch(std::runtime_error err) {
			std::cerr << err.what() << std::endl;
			continue;
		}
		
		std::map<VirtualRegionType, std::vector<VirtualRegion>>& VirtualRegionsMap = record.RegionsMap;
		unsigned int neededRegionSize = VirtualRegionsMap[VirtualRegionType::SUBREAD].size();
		if (selection & PB_DISCARD_FILTER)
			neededRegionSize += VirtualRegionsMap[VirtualRegionType::FILTERED].size();
	
		/* Check if it is worth submitting this job */
		if ( !(selection & PB_KEEP_INVALID) && neededRegionSize <= 0) continue;
		
		neededRegionSize += VirtualRegionsMap[VirtualRegionType::ADAPTER].size()
		                  + VirtualRegionsMap[VirtualRegionType::HQREGION].size();
		
		if (neededRegionSize <= 0U) {			
			fprintf(stderr, "ZMW %u has no regions\n", record.HoleNumber);
			std::cout << "This is hole " << record.HoleNumber << std::endl;
			
			std::cout << "Length " << record.Sequence.length() << " accuracy " << record.Accuracy << std::endl;
			
			std::map<VirtualRegionType, std::vector<VirtualRegion>>::const_iterator iter;
			std::cout << "Virtual region categories: " << VirtualRegionsMap.size() << std::endl; 
			for (iter = VirtualRegionsMap.begin(); iter != VirtualRegionsMap.end(); ++iter) {
				const VirtualRegionType key = iter->first;
				std::vector<VirtualRegion>::const_iterator iregions;
				std::cout << RegionToString (key) << std::endl;
				for (iregions = iter->second.begin(); iregions != iter->second.end(); ++iregions) {
					std::cout << "\t" << iregions->beginPos << " " << iregions->endPos << " " << iregions->score << std::endl;
				}
			}
			return -1;
		}
											
		/* Get the Job space slot */
		/* Get a job slot from the pool donequeue */
		again:;
		pthread_mutex_lock(&thpool_p->donequeue_p.rwmutex);
		newJob = (pb_job_t*) jobqueue_pull(&thpool_p->donequeue_p);
		pthread_mutex_unlock(&thpool_p->donequeue_p.rwmutex);
		
		if (!newJob) {
			if (thpool_p->num_threads_alive > 0) {
				/* Wait a bit and try again */
				struct timespec polling_interval = { .tv_sec = 0, .tv_nsec = 400};
				nanosleep(&polling_interval, NULL);
				goto again;
			}
			else {
				return -1;
			}
		}
		
		/* Handle regions */
		{
			newJob->HReg.nRegions = neededRegionSize;
			
			if (newJob->HReg.nAllocatedRegions < neededRegionSize) {
				newJob->HReg.Regions = (HoleRegion_t * restrict) realloc(newJob->HReg.Regions, neededRegionSize*sizeof(HoleRegion_t));
				if (newJob->HReg.Regions == NULL) return -4;
				newJob->HReg.nAllocatedRegions = neededRegionSize;
			}
			
			HoleRegion_t * restrict HR = newJob->HReg.Regions;
			if (record.HasVirtualRegionType(VirtualRegionType::ADAPTER)) {
				std::vector<VirtualRegion>::const_iterator iregions = VirtualRegionsMap[VirtualRegionType::ADAPTER].cbegin();
				while (iregions != VirtualRegionsMap[VirtualRegionType::ADAPTER].cend()) {
					HR->type = RegionType::Adapter;
					HR->start = iregions->beginPos;
					HR->stop = iregions->endPos;
					HR->quality = 0;
					HR++;
					iregions++;
				}
			}
			if (record.HasVirtualRegionType(VirtualRegionType::SUBREAD)) {
				std::vector<VirtualRegion>::const_iterator iregions = VirtualRegionsMap[VirtualRegionType::SUBREAD].cbegin();
				while (iregions != VirtualRegionsMap[VirtualRegionType::SUBREAD].cend()) {
					HR->type = RegionType::Insert;
					HR->start = iregions->beginPos;
					HR->stop = iregions->endPos;
					HR->quality = 0;
					HR++;
					iregions++;
				}
			}
			if ((selection & PB_DISCARD_FILTER) && record.HasVirtualRegionType(VirtualRegionType::FILTERED)) {
				std::vector<VirtualRegion>::const_iterator iregions = VirtualRegionsMap[VirtualRegionType::FILTERED].cbegin();
				while (iregions != VirtualRegionsMap[VirtualRegionType::FILTERED].cend()) {
					HR->type = RegionType::Insert;
					HR->start = iregions->beginPos;
					HR->stop = iregions->endPos;
					HR->quality = 0;
					HR++;
					iregions++;
				}
			}
			if (record.HasVirtualRegionType(VirtualRegionType::HQREGION)) {
				int Quality = record.Accuracy;
				if (Quality < 0) Quality = 0;
				else if (Quality > 1000) Quality = 1000;
				HR->type = RegionType::HQRegion;
				HR->quality = Quality;
				VirtualRegion& region = VirtualRegionsMap[VirtualRegionType::HQREGION].front();
				HR->start = region.beginPos;
				HR->stop = region.endPos;
			}
		}

		/* Handle Sequence */
		{
			size_t Length = record.Sequence.length()+1;
			if (Length >newJob->Reads.Size) {
				newJob->Reads.Stream = (unsigned char*) realloc(newJob->Reads.Stream, Length*sizeof(unsigned char));
				if (newJob->Reads.Stream == NULL) return -5;
				newJob->Reads.Size = Length;
			}
			Length -= 1;
			memcpy(newJob->Reads.Stream, record.Sequence.c_str(), Length*sizeof(char));
			newJob->Reads.Stream[Length] = '\0';
		}
		
		/* Handle Hole Number */
		newJob->HReg.HoleNumber = record.HoleNumber;
		
		/* Handle coordinates */
		newJob->HReg.Coordinates[0] = record.HoleNumber >> 16;
		newJob->HReg.Coordinates[1] = record.HoleNumber & 0xFFFF;
		
		pthread_mutex_lock(&thpool_p->jobqueue_p.rwmutex);
		jobqueue_push(&thpool_p->jobqueue_p, (job_t*) newJob);
		pthread_mutex_unlock(&thpool_p->jobqueue_p.rwmutex);
	}
	return 0;
}

// static int
// thpool_add_subreads(threadpool_t * const restrict thpool_p, BamReader& Reader, const enum PacBioSelection selection)
// {
// 	MyZmwRecord record;
// 	pb_job_t * restrict newJob;
// 	
// 	record.RegionsMap[VirtualRegionType::ADAPTER] = std::vector<VirtualRegion>();
// 	record.RegionsMap[VirtualRegionType::SUBREAD] = std::vector<VirtualRegion>();
// 	record.RegionsMap[VirtualRegionType::HQREGION] = std::vector<VirtualRegion>();
// 	record.RegionsMap[VirtualRegionType::LQREGION] = std::vector<VirtualRegion>();
// 	record.RegionsMap[VirtualRegionType::BARCODE] = std::vector<VirtualRegion>();
// 	record.RegionsMap[VirtualRegionType::FILTERED] = std::vector<VirtualRegion>();
// 	
// 	while  (Reader.HasNext()) {
// 		std::vector<BamRecord> records = Reader.NextRaw();
// 		// Sort sources by queryStart
//     std::sort(records.begin(), records.end(),
//               [](const BamRecord& l1, const BamRecord& l2)
//               { return l1.QueryStart() < l2.QueryStart(); });
// 		try {
// 			StitchSources(records, record);
// 		} catch(std::runtime_error err) {
// 			std::cerr << err.what() << std::endl;
// 			continue;
// 		}
// 		
// 		std::map<VirtualRegionType, std::vector<VirtualRegion>>& VirtualRegionsMap = record.RegionsMap;
// 		unsigned int neededRegionSize = VirtualRegionsMap[VirtualRegionType::SUBREAD].size();
// 		if (selection & PB_DISCARD_FILTER)
// 			neededRegionSize += VirtualRegionsMap[VirtualRegionType::FILTERED].size();
// 	
// 		/* Check if it is worth submitting this job */
// 		if ( !(selection & PB_KEEP_INVALID) && neededRegionSize <= 0) continue;
// 		
// 		neededRegionSize += VirtualRegionsMap[VirtualRegionType::ADAPTER].size()
// 		                  + VirtualRegionsMap[VirtualRegionType::HQREGION].size();
// 		
// 		if (neededRegionSize <= 0U) {			
// 			fprintf(stderr, "ZMW %u has no regions\n", record.HoleNumber);
// 			std::cout << "This is hole " << record.HoleNumber << std::endl;
// 			
// 			std::cout << "Length " << record.Sequence.length() << " accuracy " << record.Accuracy << std::endl;
// 			
// 			std::map<VirtualRegionType, std::vector<VirtualRegion>>::const_iterator iter;
// 			std::cout << "Virtual region categories: " << VirtualRegionsMap.size() << std::endl; 
// 			for (iter = VirtualRegionsMap.begin(); iter != VirtualRegionsMap.end(); ++iter) {
// 				const VirtualRegionType key = iter->first;
// 				std::vector<VirtualRegion>::const_iterator iregions;
// 				std::cout << RegionToString (key) << std::endl;
// 				for (iregions = iter->second.begin(); iregions != iter->second.end(); ++iregions) {
// 					std::cout << "\t" << iregions->beginPos << " " << iregions->endPos << " " << iregions->score << std::endl;
// 				}
// 			}
// 			return -1;
// 		}
// 											
// 		/* Get the Job space slot */
// 		/* Get a job slot from the pool donequeue */
// 		again:;
// 		pthread_mutex_lock(&thpool_p->donequeue_p.rwmutex);
// 		newJob = (pb_job_t*) jobqueue_pull(&thpool_p->donequeue_p);
// 		pthread_mutex_unlock(&thpool_p->donequeue_p.rwmutex);
// 		
// 		if (!newJob) {
// 			if (thpool_p->num_threads_alive > 0) {
// 				/* Wait a bit and try again */
// 				struct timespec polling_interval = { .tv_sec = 0, .tv_nsec = 400};
// 				nanosleep(&polling_interval, NULL);
// 				goto again;
// 			}
// 			else {
// 				return -1;
// 			}
// 		}
// 		
// 		/* Handle regions */
// 		{
// 			newJob->HReg.nRegions = neededRegionSize;
// 			
// 			if (newJob->HReg.nAllocatedRegions < neededRegionSize) {
// 				newJob->HReg.Regions = (HoleRegion_t * restrict) realloc(newJob->HReg.Regions, neededRegionSize*sizeof(HoleRegion_t));
// 				if (newJob->HReg.Regions == NULL) return -4;
// 				newJob->HReg.nAllocatedRegions = neededRegionSize;
// 			}
// 			
// 			HoleRegion_t * restrict HR = newJob->HReg.Regions;
// 			if (record.HasVirtualRegionType(VirtualRegionType::ADAPTER)) {
// 				std::vector<VirtualRegion>::const_iterator iregions = VirtualRegionsMap[VirtualRegionType::ADAPTER].cbegin();
// 				while (iregions != VirtualRegionsMap[VirtualRegionType::ADAPTER].cend()) {
// 					HR->type = RegionType::Adapter;
// 					HR->start = iregions->beginPos;
// 					HR->stop = iregions->endPos;
// 					HR->quality = 0;
// 					HR++;
// 					iregions++;
// 				}
// 			}
// 			if (record.HasVirtualRegionType(VirtualRegionType::SUBREAD)) {
// 				std::vector<VirtualRegion>::const_iterator iregions = VirtualRegionsMap[VirtualRegionType::SUBREAD].cbegin();
// 				while (iregions != VirtualRegionsMap[VirtualRegionType::SUBREAD].cend()) {
// 					HR->type = RegionType::Insert;
// 					HR->start = iregions->beginPos;
// 					HR->stop = iregions->endPos;
// 					HR->quality = 0;
// 					HR++;
// 					iregions++;
// 				}
// 			}
// 			if ((selection & PB_DISCARD_FILTER) && record.HasVirtualRegionType(VirtualRegionType::FILTERED)) {
// 				std::vector<VirtualRegion>::const_iterator iregions = VirtualRegionsMap[VirtualRegionType::FILTERED].cbegin();
// 				while (iregions != VirtualRegionsMap[VirtualRegionType::FILTERED].cend()) {
// 					HR->type = RegionType::Insert;
// 					HR->start = iregions->beginPos;
// 					HR->stop = iregions->endPos;
// 					HR->quality = 0;
// 					HR++;
// 					iregions++;
// 				}
// 			}
// 			if (record.HasVirtualRegionType(VirtualRegionType::HQREGION)) {
// 				int Quality = record.Accuracy;
// 				if (Quality < 0) Quality = 0;
// 				else if (Quality > 1000) Quality = 1000;
// 				HR->type = RegionType::HQRegion;
// 				HR->quality = Quality;
// 				VirtualRegion& region = VirtualRegionsMap[VirtualRegionType::HQREGION].front();
// 				HR->start = region.beginPos;
// 				HR->stop = region.endPos;
// 			}
// 		}
// 
// 		/* Handle Sequence */
// 		{
// 			size_t Length = record.Sequence.length()+1;
// 			if (Length >newJob->Reads.Size) {
// 				newJob->Reads.Stream = (unsigned char*) realloc(newJob->Reads.Stream, Length*sizeof(unsigned char));
// 				if (newJob->Reads.Stream == NULL) return -5;
// 				newJob->Reads.Size = Length;
// 			}
// 			Length -= 1;
// 			memcpy(newJob->Reads.Stream, record.Sequence.c_str(), Length*sizeof(char));
// 			newJob->Reads.Stream[Length] = '\0';
// 		}
// 		
// 		/* Handle Hole Number */
// 		newJob->HReg.HoleNumber = record.HoleNumber;
// 		
// 		/* Handle coordinates */
// 		newJob->HReg.Coordinates[0] = record.HoleNumber >> 16;
// 		newJob->HReg.Coordinates[1] = record.HoleNumber & 0xFFFF;
// 		
// 		pthread_mutex_lock(&thpool_p->jobqueue_p.rwmutex);
// 		jobqueue_push(&thpool_p->jobqueue_p, (job_t*) newJob);
// 		pthread_mutex_unlock(&thpool_p->jobqueue_p.rwmutex);
// 	}
// 	return 0;
// }

static int
thpool_add_zmw(threadpool_t * const restrict thpool_p, ZmwReadStitcher& Reader, const enum PacBioSelection selection)
{
	MyZmwRecord record;
	pb_job_t * restrict newJob;
	
	record.RegionsMap[VirtualRegionType::ADAPTER] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::SUBREAD] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::HQREGION] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::LQREGION] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::BARCODE] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::FILTERED] = std::vector<VirtualRegion>();
	
	
	while  (Reader.HasNext()) {
		std::vector<BamRecord> records = Reader.NextRaw();
		// Sort sources by queryStart
    std::sort(records.begin(), records.end(),
              [](const BamRecord& l1, const BamRecord& l2)
              { return l1.QueryStart() < l2.QueryStart(); });
		try {
			StitchSources(records, record);
		} catch(std::runtime_error err) {
			std::cerr << err.what() << std::endl;
			continue;
		}
		
		std::map<VirtualRegionType, std::vector<VirtualRegion>>& VirtualRegionsMap = record.RegionsMap;
		unsigned int RegionFound = VirtualRegionsMap[VirtualRegionType::SUBREAD].size();
		
		if (selection & PB_DISCARD_FILTER)
			RegionFound += VirtualRegionsMap[VirtualRegionType::FILTERED].size();
					
		/* Check if it is worth submitting this job */
		if ( !(selection & PB_KEEP_INVALID) && RegionFound < 1) continue;
		
		/* Get the Job space slot */
		/* Get a job slot from the pool donequeue */
		again:;
		pthread_mutex_lock(&thpool_p->donequeue_p.rwmutex);
		newJob = (pb_job_t*) jobqueue_pull(&thpool_p->donequeue_p);
		pthread_mutex_unlock(&thpool_p->donequeue_p.rwmutex);
		
		if (!newJob) {
			if (thpool_p->num_threads_alive > 0) {
				/* Wait a bit and try again */
				struct timespec polling_interval = { .tv_sec = 0, .tv_nsec = 400};
				nanosleep(&polling_interval, NULL);
				goto again;
			}
			else {
				return -1;
			}
		}
		
		/* Handle regions */
		{
			newJob->HReg.nRegions = 2;
			
			if (newJob->HReg.nAllocatedRegions < 2) {
				newJob->HReg.Regions = (HoleRegion_t * restrict) realloc(newJob->HReg.Regions, 2*sizeof(HoleRegion_t));
				if (newJob->HReg.Regions == NULL) return -4;
				newJob->HReg.nAllocatedRegions = 2;
			}
			
			HoleRegion_t * restrict HR = newJob->HReg.Regions;
			
			int Quality = record.Accuracy;
			if (Quality < 0) Quality = 0;
			else if (Quality > 1000) Quality = 1000;
			VirtualRegion& region = record.RegionsMap[VirtualRegionType::HQREGION].front();
			if ((selection & PB_KEEP_INVALID) && region.beginPos == region.endPos) {
				HR[1].type = RegionType::HQRegion;
				HR[1].start = 0U;
				HR[1].stop = record.Sequence.length();
				HR[1].quality = Quality;
				HR[0].type = RegionType::Insert;
				HR[0].quality = 0;
				HR[0].start = 0U;
				HR[0].stop = record.Sequence.length();
			}
			else {
				HR[1].type = RegionType::HQRegion;
				HR[1].start = region.beginPos;
				HR[1].stop = region.endPos;
				HR[1].quality = Quality;
				
				HR[0].type = RegionType::Insert;
				HR[0].quality = 0;
				if (!(selection & PB_DISCARD_FILTER)) {
					HR[0].start = region.beginPos;
					HR[0].stop = region.endPos;
				}
				else {
					HR[0].start = 0U;
					HR[0].stop = record.Sequence.length();
					HR[1].start = 0U;
					HR[1].stop = record.Sequence.length();
				}
			}
		}

		/* Handle Sequence */
		{
			size_t Length = record.Sequence.length()+1;
			if (Length >newJob->Reads.Size) {
				newJob->Reads.Stream = (unsigned char*) realloc(newJob->Reads.Stream, Length*sizeof(unsigned char));
				if (newJob->Reads.Stream == NULL) return -5;
				newJob->Reads.Size = Length;
			}
			Length -= 1;
			memcpy(newJob->Reads.Stream, record.Sequence.c_str(), Length*sizeof(char));
			newJob->Reads.Stream[Length] = '\0';
		}
		
		/* Handle Hole Number */
		newJob->HReg.HoleNumber = record.HoleNumber;
		
		/* Handle coordinates */
		newJob->HReg.Coordinates[0] = record.HoleNumber >> 16;
		newJob->HReg.Coordinates[1] = record.HoleNumber & 0xFFFF;
		
		pthread_mutex_lock(&thpool_p->jobqueue_p.rwmutex);
		jobqueue_push(&thpool_p->jobqueue_p, (job_t*) newJob);
		pthread_mutex_unlock(&thpool_p->jobqueue_p.rwmutex);
	}
	return 0;
}

static int
thpool_add_invalid_zmw(threadpool_t * const restrict thpool_p, ZmwReadStitcher& Reader, const enum PacBioSelection selection)
{
	MyZmwRecord record;
	pb_job_t * restrict newJob;
	
	record.RegionsMap[VirtualRegionType::ADAPTER] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::SUBREAD] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::HQREGION] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::LQREGION] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::BARCODE] = std::vector<VirtualRegion>();
	record.RegionsMap[VirtualRegionType::FILTERED] = std::vector<VirtualRegion>();
	
	
	while  (Reader.HasNext()) {
		std::vector<BamRecord> records = Reader.NextRaw();
		// Sort sources by queryStart
    std::sort(records.begin(), records.end(),
              [](const BamRecord& l1, const BamRecord& l2)
              { return l1.QueryStart() < l2.QueryStart(); });
		try {
			StitchSources(records, record);
		} catch(std::runtime_error err) {
			std::cerr << err.what() << std::endl;
			continue;
		}
		
		std::map<VirtualRegionType, std::vector<VirtualRegion>>& VirtualRegionsMap = record.RegionsMap;
		unsigned int RegionFound = VirtualRegionsMap[VirtualRegionType::SUBREAD].size()
		                         + VirtualRegionsMap[VirtualRegionType::FILTERED].size();
		
		if (RegionFound) continue;
		
		/* Get the Job space slot */
		/* Get a job slot from the pool donequeue */
		again:;
		pthread_mutex_lock(&thpool_p->donequeue_p.rwmutex);
		newJob = (pb_job_t*) jobqueue_pull(&thpool_p->donequeue_p);
		pthread_mutex_unlock(&thpool_p->donequeue_p.rwmutex);
		
		if (!newJob) {
			if (thpool_p->num_threads_alive > 0) {
				/* Wait a bit and try again */
				struct timespec polling_interval = { .tv_sec = 0, .tv_nsec = 400};
				nanosleep(&polling_interval, NULL);
				goto again;
			}
			else {
				return -1;
			}
		}
		
		/* Handle regions */
		{
			newJob->HReg.nRegions = 2;
			
			if (newJob->HReg.nAllocatedRegions < 2) {
				newJob->HReg.Regions = (HoleRegion_t * restrict) realloc(newJob->HReg.Regions, 2*sizeof(HoleRegion_t));
				if (newJob->HReg.Regions == NULL) return -4;
				newJob->HReg.nAllocatedRegions = 2;
			}
			
			HoleRegion_t * restrict HR = newJob->HReg.Regions;
			
			HR[1].type = RegionType::HQRegion;
			HR[1].start = 0U;
			HR[1].stop = record.Sequence.length();
			HR[1].quality = 0;
			HR[0].type = RegionType::Insert;
			HR[0].quality = 0;
			HR[0].start = 0U;
			HR[0].stop = record.Sequence.length();
		}

		/* Handle Sequence */
		{
			size_t Length = record.Sequence.length()+1;
			if (Length >newJob->Reads.Size) {
				newJob->Reads.Stream = (unsigned char*) realloc(newJob->Reads.Stream, Length*sizeof(unsigned char));
				if (newJob->Reads.Stream == NULL) return -5;
				newJob->Reads.Size = Length;
			}
			Length -= 1;
			memcpy(newJob->Reads.Stream, record.Sequence.c_str(), Length*sizeof(char));
			newJob->Reads.Stream[Length] = '\0';
		}
		
		/* Handle Hole Number */
		newJob->HReg.HoleNumber = record.HoleNumber;
		
		/* Handle coordinates */
		newJob->HReg.Coordinates[0] = record.HoleNumber >> 16;
		newJob->HReg.Coordinates[1] = record.HoleNumber & 0xFFFF;
		
		pthread_mutex_lock(&thpool_p->jobqueue_p.rwmutex);
		jobqueue_push(&thpool_p->jobqueue_p, (job_t*) newJob);
		pthread_mutex_unlock(&thpool_p->jobqueue_p.rwmutex);
	}
	return 0;
}


extern "C"
int dispatchPacBioBAM(const struct Profile * const restrict prf,
                      const PacBioBAM_t * const restrict PBBAM,
#ifdef PRF_CORE_PCRE
                      struct RegEx * const restrict Regex,
#endif
                      const OutputType_t * const restrict OutputType,
                      const PacBioDispatchOptions_t * Options,
                      const size_t nCPUs)
{
	int res= SUCCESS;
	if (OutputVerbose) fputs("Dispatching PacBioBAM data...\n", stderr); 
	
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
	void* (*ThreadFct)(threadarg_t* const restrict) = NULL;
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
	common.PrintLock = PTHREAD_MUTEX_INITIALIZER;
	common.Histograms = Histograms;
	common.Counters = Missed;
	common.Compute = CoreCompute;
	common.OutputType = OutputType;
#ifdef PRF_CORE_PCRE
	common.regex = Regex;
#endif

	
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
	PbiFilter Filter{};
	
	/* Given a list of ZMWs */
	if (Options->ZMW) {
		std::vector<int32_t> zmws;
		zmws.reserve(Options->nZMW);
		for (int i=0; i<Options->nZMW; i++) zmws.emplace_back(Options->ZMW[i]);
		Filter.Add( PbiZmwFilter{ zmws } );
	}
	
	/* Given filter score */
	if (Options->Selection & PB_HAS_FILTER_SCORE) {
		Filter.Add(PbiReadAccuracyFilter{ (float) Options->minReadAccuracy / 1000.0f, Compare::GREATER_THAN_EQUAL });
		Filter.Add(PbiReadAccuracyFilter{ (float) Options->maxReadAccuracy / 1000.0f, Compare::LESS_THAN_EQUAL });
	}

	if (!PBBAM->scraps.empty()) {
		ZmwReadStitcher Reader(PBBAM->subreads, PBBAM->scraps, Filter);
		if (Options->Selection & (PB_DISPATCH_ZMW | PB_DISPATCH_HQREGION))
			res = thpool_add_zmw(thpool, Reader, Options->Selection);
		else if (Options->Selection & PB_DISPATCH_INVALID) 
			res = thpool_add_invalid_zmw(thpool, Reader, Options->Selection);
		else
			res = thpool_add_subreads(thpool, Reader, Options->Selection);
	}
	else {
		BamReader Reader(PBBAM->subreads);
		if (Options->Selection & (PB_DISPATCH_ZMW | PB_DISPATCH_HQREGION) || (Options->Selection & PB_DISPATCH_INVALID)) {
			fputs("Without scrap file, we cannot perform invalid or HQ region\n", stderr);
			destroyThreadPool(thpool);
			pthread_mutex_destroy(&(common.PrintLock));
			if (Histograms) free(Histograms);
			if (Missed) free(Missed);
			return -4;
		}	
		else {
// 			res = thpool_add_subreads(thpool, Reader, Options->Selection);
		}
	}
	
	if (res != 0 ) {
		fprintf(stderr, "Master failed to send jobs, error code %i\n", res);
		destroyThreadPool(thpool);
		pthread_mutex_destroy(&(common.PrintLock));
		if (Histograms) free(Histograms);
		if (Missed) free(Missed);
		return -5;
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
				res = -6;
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
				res = -7;
			}
		}
	}
	

	/*************************************************************************/
	/*                    CHECK FOR THREADS ERROR                            */
	/*************************************************************************/

			
	/*************************************************************************/
	/*                      DESTROY THE THREAD POOL                          */
	/*************************************************************************/

	destroyThreadPool(thpool);
	pthread_mutex_destroy(&(common.PrintLock));
	
	if (Histograms) free(Histograms);
	if (Missed) free(Missed);
	return res;
}

extern "C"
int dispatchPacBioBAMExt(const PacBioBAM_t * const restrict PBBAM,
												 threadpool_t * const restrict thpool,
                         const PacBioDispatchOptions_t * Options)
{
	int res= SUCCESS;
	if (OutputVerbose) fputs("Dispatching PacBioBAM data...\n", stderr); 
	
	/*************************************************************************/
	/*                         CONSISTENCY CHECK                             */
	/*************************************************************************/
	if ((Options->Selection & PB_BEST_SUBREAD_OF_ZMW) && !(Options->Selection & PB_DISPATCH_SUBREADS)) {
		fputs("Aggregation is only possible on subreads!\n", stderr);
		return -10;
	}
	
	/*************************************************************************/
	/*                        ADD TASKS TO JOB QUEUE                         */
	/*************************************************************************/
	struct timeval _t0, _t1;
	gettimeofday(&_t0,0);
	PbiFilter Filter{};
	
	/* Given a list of ZMWs */
	if (Options->ZMW) {
		std::vector<int32_t> zmws;
		zmws.reserve(Options->nZMW);
		for (int i=0; i<Options->nZMW; i++) zmws.emplace_back(Options->ZMW[i]);
		Filter.Add( PbiZmwFilter{ zmws } );
	}
	
	/* Given filter score */
	if (Options->Selection & PB_HAS_FILTER_SCORE) {
		Filter.Add(PbiReadAccuracyFilter{ (float) Options->minReadAccuracy / 1000.0f, Compare::GREATER_THAN_EQUAL });
		Filter.Add(PbiReadAccuracyFilter{ (float) Options->maxReadAccuracy / 1000.0f, Compare::LESS_THAN_EQUAL });
	}

	
	ZmwReadStitcher Reader(PBBAM->subreads, PBBAM->scraps, Filter);
	if (Options->Selection & (PB_DISPATCH_ZMW | PB_DISPATCH_HQREGION))
		res = thpool_add_zmw(thpool, Reader, Options->Selection);
	else if (Options->Selection & PB_DISPATCH_INVALID) 
		res = thpool_add_invalid_zmw(thpool, Reader, Options->Selection);
	else
		res = thpool_add_subreads(thpool, Reader, Options->Selection);
	
	if (res != 0 ) {
		fprintf(stderr, "Master failed to send jobs, error code %i\n", res);
		return -4;
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
	/*                    CHECK FOR THREADS ERROR                            */
	/*************************************************************************/

			
	/*************************************************************************/
	/*                      DESTROY THE THREAD POOL                          */
	/*************************************************************************/

	return res;
}
