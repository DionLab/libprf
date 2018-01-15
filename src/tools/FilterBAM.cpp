#include <iostream>
#include <string>
#include <cstdio>
#include "pbbam/BamFile.h"
#include "pbbam/BamReader.h"
#include "pbbam/BamWriter.h"
#include <getopt.h>
#include "pfVersion.h"
using namespace std;
using namespace PacBio::BAM;

static const char opt_to_test[] = "i:o:f:ch";
static const struct option long_options[] =
{
  /*
	 * These options set a flag.
	 */
	
  /*
	 * These options don't set a flag. We distinguish them by their indices.
	 */
	{"help",    no_argument,       	0,	'h'},
	{"in",			required_argument,	0,	'i'},
	{"out",			required_argument,	0,	'o'},
	{"filter",	required_argument,	0,	'f'},
	{"check",		no_argument,				0,	'c'},
	{0, 0, 0, 0}
};

static void __attribute__((noreturn)) Usage(FILE* stream)
{
	fputs(
	" FilterBAM --in BAM file --out BAM file [--filter file]\n"
	"    --in  file                [-i] : input BAM file\n"
	"    --out file                [-o] : output BAM file\n"
	"    --filter file             [-f] : filter file containing a full\n"
	"                                     smartcell ZMW name per line\n"
	"    --check                   [-c] : print out ZMW in file\n\n"
	"  Other\n"
	"    --help                    [-h] : output command help\n\n"
	"This is version " PRF_VERSION ".\n",
	stream);
  exit(1);
}

typedef struct ReadId {
	unsigned int HoleNumber;
	int Start;
} ReadId_t;

static int compare_Index(const void * a, const void * b)
{
	const ReadId_t * const da = (const ReadId_t*) a;
	const ReadId_t * const db = (const ReadId_t*) b;
	
	if (da->HoleNumber < db->HoleNumber) {
		return -1;
	}
	else if (da->HoleNumber > db->HoleNumber) {
		return 1;
	}
	else {
		if (da->Start < db->Start) {
			return -1;
		}
		else if (da->Start > db->Start) {
			return 1;
		}
		else {
			return 0;
		}
	}
}


int main(int argc, char * argv[]) 
{
		int Check = 0;
		const char * InputBamName = NULL;
		const char * OutputBAMName = NULL;
		const char * FilterName = NULL;
	
		if (argc == 1) Usage(stderr);
		while (1) {
			/* getopt_long stores the option index here. */
			int option_index = 0;
			const int c = getopt_long (argc, argv, opt_to_test, long_options, &option_index);

			/* Detect the end of the options. */
			if (c == -1) break;
			switch (c) {
				case 'i':
					InputBamName = optarg;
					break;
				case 'o':
					OutputBAMName = optarg;
					break;
				case 'f':
					FilterName = optarg;
					break;
				case 'c':
					Check = 1;
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
		
		
		BamFile in{ InputBamName };
		BamReader inputBAM(in);
		
		if (Check == 0) {
			size_t IndexSize = 2048UL;
			ReadId_t * Index = (ReadId_t*) malloc(IndexSize*sizeof(ReadId_t));
			if (Index == NULL) {
				fputs("Unable to allocate memory for indexing filter list\n", stderr);
				return 1;
			}
		
			const BamHeader& hdr = in.Header();
			
			BamRecord record;
			BamWriter outputBAM(string(OutputBAMName), hdr);
			
			inputBAM.GetNext(record);
			size_t CellNameLength = record.FullName().find_first_of('/');
			if ( CellNameLength == string::npos) {
				cerr << "Unable to determine SMRT cell name from " << record.FullName() << endl;
				return 1;
			}
			const string CellName = record.FullName().substr(0, CellNameLength);
			
			FILE *win = fopen(FilterName,"r");
			if (win == NULL) {
				cerr << "Provide a readable white list file!" << endl;
				return 1;
			}
			
			char * ptr = NULL;
			size_t n = 0;
			unsigned int count = 0U;
			unsigned int TotalCount = 0U;
			while (!feof(win)) {
				if (getline(&ptr, &n, win) > 0) {
					string line = string(ptr);
					if (CellName.compare(line.substr(0,CellNameLength)) == 0) {
						const size_t first = line.find_first_of('/');
						if (first != string::npos) {
							const size_t second = line.find_first_of('/', 1UL+first);
							if (second != string::npos) {
								const size_t underscore = line.find_first_of('_', 1UL+second);
								if (underscore != string::npos) {
									Index[count].HoleNumber = stoi(line.substr(first+1UL, second-first-1UL));
									Index[count].Start = stoi(line.substr(second+1UL, underscore-second-1UL));
								}
								else {
									Index[count].HoleNumber = stoi(line.substr(first+1UL, second-first-1UL));
									Index[count].Start = -1;
								}
							}
							else {
								Index[count].HoleNumber = stoi(line.substr(first+1UL, second-first-1UL));
								Index[count].Start = -1;
							}
							if (++count == IndexSize) {
								IndexSize <<= 1;
								Index = (ReadId_t*) realloc(Index, sizeof(ReadId_t)*IndexSize);
								if (Index == NULL) {
									fputs("Unable to reallocate memory for indexing filter list\n", stderr);
									return 1;
								}
							}
						}
					}
	// 				else {
	// 						cout << line.substr(0,CellNameLength) << " is not " << CellName << endl;
	// 				}
				}
				TotalCount++;
			}
			free(ptr);
			fclose(win);
			
			qsort(Index, count, sizeof(ReadId_t), compare_Index);
			for (int i=0; i<count; i++) printf("%i\t%u\t%i\n", i, Index[i].HoleNumber, Index[i].Start);
			cerr << count << "/" << TotalCount << " lines in " << string(FilterName) << " satisfy SMRT cell name " << CellName << "." << endl;
			
			const ReadId_t * CurrentFilter = Index;
			const ReadId_t * const FilterLimit = &Index[count];
			do {

			Test: ;
				if (record.HoleNumber() == CurrentFilter->HoleNumber) {
					if (CurrentFilter->Start < 0) {
						outputBAM.Write(record);
					}
					else if (record.QueryStart() == CurrentFilter->Start) {
						outputBAM.Write(record);
						if ( (uintptr_t) ++CurrentFilter >= (uintptr_t) FilterLimit) break;
					}
				}
				else if (record.HoleNumber() > CurrentFilter->HoleNumber) {
					if ( (uintptr_t) ++CurrentFilter >= (uintptr_t) FilterLimit) break;
				}
			} while ( inputBAM.GetNext(record));
			
		}
		else {
			BamRecord record;
			while ( inputBAM.GetNext(record)) {
				cerr << "ZMW " << record.FullName() << endl;
			}
		}
		
		return 0;
}




