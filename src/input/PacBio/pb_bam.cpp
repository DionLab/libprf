#include <cstdlib>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include "pbbam/BamFile.h"
#include "pbbam/BamReader.h"
#include "pbbam/virtual/ZmwReadStitcher.h"
#include "pbbam/PbiFilterTypes.h"
#include "pbbam/virtual/VirtualRegionType.h"
#include "pbbam/virtual/VirtualRegionTypeMap.h"
#define restrict 

using namespace PacBio::BAM;

struct ZmwRecord {
	int Accuracy;
	ZmwType Type;
	std::map<VirtualRegionType, std::vector<VirtualRegion> > RegionsMap;
	
	bool HasVirtualRegionType(const VirtualRegionType regionType) const {
		return RegionsMap.find(regionType) != RegionsMap.end(); 
	}
};

struct ZmwRecordExt {
	int Accuracy;
	unsigned int x;
	unsigned int y;
	unsigned int HoleNumber;
	unsigned int Length;
	ZmwType Type;
	std::map<VirtualRegionType, std::vector<VirtualRegion> > RegionsMap;
};

extern "C" {
#include "pb_bam.h"
#include "pfOutput.h"
#include "pb_common.h"
}

using namespace PacBio::BAM;

extern "C"
int isPacBioBAM(const char * const restrict FileName)
{
	int res = 0;
	const int fid = open(FileName, O_RDONLY);
	if (fid < 0) return fid;
	char buffer[4];
	if (read(fid, buffer, 4) != 4) {
		res = -2;
		goto FIN;
	}
	
	res  = ( (unsigned char) buffer[0] == 0x1F );
	res &= ( (unsigned char) buffer[1] == 0x8B );
	res &= ( (unsigned char) buffer[2] == 0x08 );
	res &= ( (unsigned char) buffer[3] == 0x04 );
	
	FIN:
	close(fid);
	return res;
}

extern "C"
PacBioBAM_t* OpenPacBioBAM(const char * const restrict File) {
	std::string tmp = std::string(File);
	BamFile bamTmp = BamFile(tmp);
	
	if (!bamTmp.IsPacBioBAM()) {
		fprintf(stderr, "File %s is not a PacBio BAM file\n", File);
		return nullptr;
	}
	if (!bamTmp.PacBioIndexExists()) {
		fprintf(stderr, "PacBio index not found for %s\n", File);
		return nullptr;
	}
	
	PacBioBAM_t * result = new PacBioBAM_t;
	std::string::size_type pos = tmp.find("subreads");
	
	if( pos != std::string::npos) {
		result->subreads = tmp;
		result->scraps   = tmp.replace(pos, 8, "scraps");
		
		try {
			BamFile bamTmp2  = BamFile(result->scraps);
			if (!bamTmp2.IsPacBioBAM()) {
				fprintf(stderr, "File %s is not a PacBio BAM file, going as if none exists.", result->scraps.c_str());
			}
			if (!bamTmp2.PacBioIndexExists()) {
				fprintf(stderr, "PacBio index not found for %s\n, please create one", result->scraps.c_str());
				return nullptr;
			}
		}
		catch ( const std::exception& e) {
			std::cerr << "Error opening " << tmp << e.what();
		}
		if (OutputVerbose) {
				fprintf(stderr, "Pac Bio BAM files:\n\tSubreads: %s\n\tScraps:   %s\n", result->subreads.c_str(),
				        result->scraps.c_str());  
			}
		return result;
	}
	else {
		pos = tmp.find("scraps");
		if( pos == std::string::npos) {
			fprintf(stderr, "Unable to identify bam file %s as subreads or scraps.\n", File);
			delete result;
			return nullptr;
		}
		else {
			result->scraps   = tmp;
			result->subreads = tmp.replace(pos, 6, "subreads");
			BamFile bamTmp2  = BamFile(result->subreads);
		
			if (!bamTmp2.IsPacBioBAM()) {
				fprintf(stderr, "File %s is not a PacBio BAM file\n", result->subreads.c_str());
				delete result;
				return nullptr;
			}
			if (!bamTmp2.PacBioIndexExists()) {
				fprintf(stderr, "PacBio index not found for %s\n", result->subreads.c_str());
				delete result;
				return nullptr;
			}
			
			if (OutputVerbose) {
				fprintf(stderr, "Pac Bio BAM files:\n\tSubreads: %s\n\tScraps:   %s\n", result->subreads.c_str(),
				        result->scraps.c_str());  
			}
			return result;
		}
	}
}

extern "C"
void ClosePacBioBAM(const PacBioBAM_t * const restrict PBBAM) 
{
	PBBAM->subreads.empty();
	PBBAM->scraps.empty();
}

static
void StitchSources(const std::vector<BamRecord>& sources, ZmwRecord& record)
{
    const auto& firstRecord = sources[0];
    const auto& lastRecord = sources[sources.size() - 1];
		bool HasReadAccuracy = false;
		bool HasScrapZmwType = false;
		std::map<VirtualRegionType, std::vector<VirtualRegion>>& virtualRegionsMap = record.RegionsMap;
		
    // initialize capacity
    const auto stitchedSize = lastRecord.QueryEnd() - firstRecord.QueryStart();
		virtualRegionsMap.clear();
		
//  		fprintf(stderr, "ZMW %u has %u bases:\n", record.HoleNumber, stitchedSize);
    // Stitch using tmp vars
    for(auto& b : sources)
    {

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
}

static
void StitchSourcesExt(const std::vector<BamRecord>& sources, ZmwRecordExt& record)
{
    const auto& firstRecord = sources[0];
    const auto& lastRecord = sources[sources.size() - 1];
		bool HasReadAccuracy = false;
		bool HasScrapZmwType = false;
		std::map<VirtualRegionType, std::vector<VirtualRegion>>& virtualRegionsMap = record.RegionsMap;
		
    // initialize capacity
    const auto stitchedSize = lastRecord.QueryEnd() - firstRecord.QueryStart();
		virtualRegionsMap.clear();
		
//  		fprintf(stderr, "ZMW %u has %u bases:\n", record.HoleNumber, stitchedSize);
    // Stitch using tmp vars
		unsigned int Length = 0U;
    for(auto& b : sources)
    {
			 Length += b.Sequence().size();
       if (b.HasScrapRegionType())
        {
            const VirtualRegionType regionType = b.ScrapRegionType();
            virtualRegionsMap[regionType].emplace_back(regionType, b.QueryStart(), b.QueryEnd());
				}
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

        if (b.HasScrapZmwType())
        {
            if (!HasScrapZmwType)
                record.Type = b.ScrapZmwType();
            else if (record.Type != b.ScrapZmwType())
                throw std::runtime_error("ScrapZmwTypes do not match");
					
        }
    }
    
    record.Length = Length;
		{
			const unsigned int HN = firstRecord.HoleNumber();
			record.HoleNumber = HN;
			record.x = HN >> 16;
			record.y = HN & 0xFFFF;
		}
}

extern "C" {
int getBAMSummary(const PacBioBAM_t * const restrict PBBAM, const unsigned int HQThreshold,
                  ZMWSummary_t * const restrict summary, FILE* restrict out )
{
	memset (summary, 0, sizeof(ZMWSummary_t));
	ZmwReadStitcher Reader(PBBAM->subreads, PBBAM->scraps);
	
	if (out == NULL) {
		ZmwRecord record;
	
		record.RegionsMap[VirtualRegionType::ADAPTER] = std::vector<VirtualRegion>();
		record.RegionsMap[VirtualRegionType::SUBREAD] = std::vector<VirtualRegion>();
		record.RegionsMap[VirtualRegionType::LQREGION] = std::vector<VirtualRegion>();
		record.RegionsMap[VirtualRegionType::BARCODE] = std::vector<VirtualRegion>();
		record.RegionsMap[VirtualRegionType::FILTERED] = std::vector<VirtualRegion>();
		while  (Reader.HasNext()) {
			std::vector<BamRecord> records = Reader.NextRaw();
			try {
				StitchSources(records, record);
			} catch(std::runtime_error err) {
				std::cerr << err.what() << std::endl;
				continue;
			}
			
			std::map<VirtualRegionType, std::vector<VirtualRegion>>& VirtualRegionsMap = record.RegionsMap;
			
			/*	
				unsigned int nSequencing;
				unsigned int nHQaboveThreshold;
				unsigned int nBelowHQThreshold;
				unsigned int nInvalidHQRegion;
				unsigned int nWithinHQSubreads;
				unsigned int nOutofHQSubreads;
				unsigned int nBelowHQThresholdWithin;
				unsigned int nBelowHQThresholdOutof;
				unsigned int nInvalidHQRegionSubreads;
			*/
			
			const unsigned int nSubreads = VirtualRegionsMap[VirtualRegionType::SUBREAD].size();
			const unsigned int nFiltered = VirtualRegionsMap[VirtualRegionType::FILTERED].size();
			
			summary->nSequencing++;
			/* Check if it is worth submitting this job */
			if ( (nSubreads+nFiltered) < 1) summary->nInvalidHQRegion++;
			
			const int Quality = record.Accuracy;
			if (Quality < HQThreshold) {
				summary->nBelowHQThreshold++;
				summary->nBelowHQThresholdWithin += nFiltered;
			}
			else {
				summary->nHQaboveThreshold++;
				summary->nWithinHQSubreads += nSubreads;
			}
		}
		record.RegionsMap.empty();
	}
	else {
		ZmwRecordExt record;
		const char *QualityText[] = { "LOWQUAL", "HIGHQUAL", "INVALID" };
		
	
		record.RegionsMap[VirtualRegionType::ADAPTER] = std::vector<VirtualRegion>();
		record.RegionsMap[VirtualRegionType::SUBREAD] = std::vector<VirtualRegion>();
		record.RegionsMap[VirtualRegionType::LQREGION] = std::vector<VirtualRegion>();
		record.RegionsMap[VirtualRegionType::BARCODE] = std::vector<VirtualRegion>();
		record.RegionsMap[VirtualRegionType::FILTERED] = std::vector<VirtualRegion>();
		while  (Reader.HasNext()) {
			std::vector<BamRecord> records = Reader.NextRaw();
			try {
				StitchSourcesExt(records, record);
			} catch(std::runtime_error err) {
				std::cerr << err.what() << std::endl;
				continue;
			}
			const char * restrict QualPtr;
			std::map<VirtualRegionType, std::vector<VirtualRegion>>& VirtualRegionsMap = record.RegionsMap;
			
			/*	
				unsigned int nSequencing;
				unsigned int nHQaboveThreshold;
				unsigned int nBelowHQThreshold;
				unsigned int nInvalidHQRegion;
				unsigned int nWithinHQSubreads;
				unsigned int nOutofHQSubreads;
				unsigned int nBelowHQThresholdWithin;
				unsigned int nBelowHQThresholdOutof;
				unsigned int nInvalidHQRegionSubreads;
			*/
			
			const unsigned int nSubreads = VirtualRegionsMap[VirtualRegionType::SUBREAD].size();
			const unsigned int nFiltered = VirtualRegionsMap[VirtualRegionType::FILTERED].size();
			
			summary->nSequencing++;
			/* Check if it is worth submitting this job */
			const unsigned int TotalSubreads = nSubreads+nFiltered;
			if ( TotalSubreads < 1) { 
				summary->nInvalidHQRegion++;
				QualPtr = QualityText[2];
			}
			else {
				if (record.Accuracy < HQThreshold) {
					summary->nBelowHQThreshold++;
					summary->nBelowHQThresholdWithin += nFiltered;
					QualPtr = QualityText[0];
				}
				else {
					summary->nHQaboveThreshold++;
					summary->nWithinHQSubreads += nSubreads;
					QualPtr = QualityText[1];
				}
			}
			
			fprintf(out, "%u\t%u\t%u\t%s\t%u\t%u\t%u\n",
			       record.HoleNumber, record.x, record.y, QualPtr, record.Accuracy, record.Length, TotalSubreads);
		}
		record.RegionsMap.empty();
	}
		
	return 0;
	
}
}
