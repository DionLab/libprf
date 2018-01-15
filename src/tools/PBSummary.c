#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include "pb_hdf5.h"
#include "pb_bam.h"
#define HQ_THRESHOLD 750
#define str(x) #x
#define TO_CHAR(x) str(x)

static void __attribute__((noreturn)) Usage(FILE * stream)
{
  fputs(
	" PBSummary [options] PacBio bas hdf5 or BAM file"
	" Options:\n"
	"   -l        : output summary as a list rather than as a table\n"
	"   -o <file> : export more info of all ZMW found\n\n"
	"  Other\n"
	"   -h : output command help\n\n",
	stream);
  exit(0);
}

_Bool OutputVerbose = true;

int main(int argc, char * argv[])
{
	ZMWSummary_t Summary;
	const char * restrict inputFileName = NULL;
	const char * restrict outputFileName = NULL;
	FILE * out = NULL;
	int c;
	int asList = 0;
	opterr = 0;

  while ((c = getopt (argc, argv, "hlo:")) != -1) {
    switch (c)
      {
      case 'l':
        asList = 1;
        break;
      case 'h':
	       Usage(stdout);
				 break;
			case 'o':
					outputFileName = optarg;
					break;
      default:
				fprintf(stderr, "Unknown option %c\n", c); 
        Usage(stderr);
      }
	}
	
  if (optind >= argc) {
    fputs("Expected arguments after options\n", stderr);
    Usage(stderr);
  }
  else {
    inputFileName = argv[optind];
  }
  
	fprintf(stderr, "Testing file %s\n", inputFileName);
	if (outputFileName) {
		out = fopen(outputFileName, "w");
	}
	if (isPacBioBAM(inputFileName)) {
		PacBioBAM_t * const PBBAM = OpenPacBioBAM(inputFileName);
		if (PBBAM == NULL) {
			fprintf(stderr, "Error opening Pac Bio file %s\n", inputFileName);
			exit(1);
		}
		
		getBAMSummary(PBBAM, HQ_THRESHOLD, &Summary, out);
		ClosePacBioBAM(PBBAM);
	}
	else if (isPacBioH5(inputFileName)) {
		PacBio_t * const restrict PBS = PacBioOpen(inputFileName);
		if (PBS == NULL) {
			fprintf(stderr, "Error opening Pac Bio file %s\n", inputFileName);
			exit(1);
		}
		
		fprintf(stderr,"Base file name : %s\n"
						"Directory      : %s\n"
						"Parts          : %u\n",
						PBS->BaseFileName, PBS->Directory, PBS->nParts);
		for (unsigned int i=0; i<PBS->nParts; i++) {
			fprintf(stderr,"               : %s\n", PBS->PartsFileName[i]);
		}
		
		fprintf(stderr,"Holes          : %u\n", PBS->nHoles);
		
		if (!(PBS->Content & BASECALLING)) {
			fputs("Error Pac Bio file does not contain Basecalling data\n", stderr);
			exit(1);
		}
		
		
		getHDFSummary(PBS, HQ_THRESHOLD, &Summary, out);
		
		if (PacBioClose(PBS) != SUCCESS) {
			fprintf(stderr, "Error closing Pac Bio...\n");
			exit(1);
		}
	}
	else {
		fprintf(stderr, "The file %s is not recognized as either HDF5 nor BAM\n", inputFileName);
		exit(1);
	}
	
	if (out) fclose(out);
	
	if (!asList) {
		printf("Sequencing reads\t"
					 "Reads above HQ " TO_CHAR(HQ_THRESHOLD) "\t"
					 "Reads below HQ " TO_CHAR(HQ_THRESHOLD) "\t"
					 "Invalid reads\t"
					 "Subreads within HQ region " TO_CHAR(HQ_THRESHOLD) "\t"
					 "Subreads out of HQ region " TO_CHAR(HQ_THRESHOLD) "\t"
					 "Subreads within HQ region in reads below HQ " TO_CHAR(HQ_THRESHOLD) "\t"
					 "Subreads out of HQ region in reads below HQ " TO_CHAR(HQ_THRESHOLD) "\t"
					 "Invalid subreads\n%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
					 Summary.nSequencing, Summary.nHQaboveThreshold, Summary.nBelowHQThreshold, Summary.nInvalidHQRegion,
					 Summary.nWithinHQSubreads, Summary.nOutofHQSubreads,
					 Summary.nBelowHQThresholdWithin, Summary.nBelowHQThresholdOutof, 
					 Summary.nInvalidHQRegionSubreads);
	}
	else {
		printf("Sequencing reads : %u\n"
					 "Reads above HQ " TO_CHAR(HQ_THRESHOLD) " : %u\n"
					 "Reads below HQ " TO_CHAR(HQ_THRESHOLD) " : %u\n"
					 "Invalid reads : %u\n"
					 "Subreads within HQ region " TO_CHAR(HQ_THRESHOLD) " : %u\n"
					 "Subreads out of HQ region " TO_CHAR(HQ_THRESHOLD) " : %u\n"
					 "Subreads within HQ region in reads below HQ " TO_CHAR(HQ_THRESHOLD) " : %u\n"
					 "Subreads out of HQ region in reads below HQ " TO_CHAR(HQ_THRESHOLD) " : %u\n"
					 "Invalid subreads : %u\n",
					 Summary.nSequencing, Summary.nHQaboveThreshold, Summary.nBelowHQThreshold, Summary.nInvalidHQRegion,
					 Summary.nWithinHQSubreads, Summary.nOutofHQSubreads,
					 Summary.nBelowHQThresholdWithin, Summary.nBelowHQThresholdOutof, 
					 Summary.nInvalidHQRegionSubreads);
	}
	exit(0);
}
