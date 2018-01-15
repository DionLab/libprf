#include <iostream>
#include <pbbam/BamFile.h>
#include <pbbam/BamReader.h>
#include "pbbam/virtual/ZmwReadStitcher.h"
#include "pbbam/PbiFilterTypes.h"

using namespace PacBio::BAM;

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

int main (int argc, char * argv[]) {
	
// 	BamFile BF = BamFile(std::string(argv[1]));
// 	
// 	if (BF.IsPacBioBAM()) {
// 		std::cout << "This is a PB BAM file" << std::endl;
// 		BamHeader hdr(BF.Header());
// 		std::cout << "PacBio versions is " << hdr.PacBioBamVersion() << std::endl;
// 		
// 		BamReader Reader(BF);
// 		BamRecord record;
// 		
// 		Reader.GetNext(record);
// 		
// 		std::cout << "This is hole " << record.HoleNumber() << " of type " <<  ToString(record.Type()) << std::endl;
// 		
// 		
// 	}
	
	const std::string subreads = std::string(argv[1]);
	const std::string scraps   = std::string(argv[2]);
// 	const int32_t zmwid = atoi(argv[3]);
// 		PbiFilter filter{ PbiZmwFilter{ zmwid } };
	
	PbiFilter filter{ };
	ZmwReadStitcher Reader(subreads, scraps, filter);
	
	while  (Reader.HasNext()) {
		std::vector<BamRecord> records = Reader.NextRaw();
		// Sort sources by queryStart
    std::sort(records.begin(), records.end(),
              [](const BamRecord& l1, const BamRecord& l2)
              { return l1.QueryStart() < l2.QueryStart(); });
		
		std::cout << "This is hole " << records[0].HoleNumber() << " of type " <<  ToString(records[0].Type()) << std::endl;
		
		for(auto& b : records)
    {
        if (b.HasScrapRegionType())
        {
            const VirtualRegionType regionType = b.ScrapRegionType();
            fprintf(stderr,"SCRAP\t%s\t%u\t%u",RegionToString(b.ScrapRegionType()).c_str(), b.QueryStart(), b.QueryEnd());
        }
        else 
  					fprintf(stderr,"NO SCRAP\t%s\t%u\t%u",ToString( b.Type()).c_str(), b.QueryStart(), b.QueryEnd());

        if (b.HasLocalContextFlags())
        {
//             std::pair<int, int> barcodes{-1, -1};
//             if (b.HasBarcodes()) barcodes = b.Barcodes();
						fprintf(stderr,"\tCONTEXT\tSUBREAD\t%u\t%u", b.QueryStart(), b.QueryEnd());
        }

 				if (b.HasReadAccuracy())
 					fprintf(stderr, "\t%lf\n", static_cast<float>(b.ReadAccuracy()));
 				else 
 					fputc('\n', stderr);
					
//         if (b.HasScrapZmwType())
//         {
//             if (!HasScrapZmwType)
//                 record.Type = b.ScrapZmwType();
//             else if (record.Type != b.ScrapZmwType())
//                 throw std::runtime_error("ScrapZmwTypes do not match");
// 					
//         }
        
    }
	}
	
	return 0;
}

