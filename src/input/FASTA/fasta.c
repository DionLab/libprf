/* TODO:
 *  Use the information from the filesystem to adjust buffer size reading.
 */
#include "prf_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include "pfInput.h"

/* Beware of alignment before changing these !!! */
#define BUFFER_SIZE         1024*1024
#define CHAIN_SEGMENT_SIZE      16384
#define MAX_NUMBER_OF_CHAIN      8192

/* NOTE: For mmap reading we extend the sequence count by one and store
 *       the file size in the extra data, that is the sequence count + 1
 * 	 structure should have
 * 		off_t Offset = file size in byte;
 *		size_t HeaderLength = 0;
 *		size_t SequenceLength = 0;
 */

extern int OutputVerbose;

int AnalyzeFASTAStructure(const int FileDescriptor, FASTAStructure * const FaS)
{
  char Buffer[BUFFER_SIZE] __attribute__((aligned(4096)));
  off_t *Offset_chain_ptr[MAX_NUMBER_OF_CHAIN] ;
  size_t *Header_chain_ptr[MAX_NUMBER_OF_CHAIN];
  size_t MaxSequenceSize = 0UL;
  s_Data * restrict DataPtr = NULL;
  _Bool HeaderPending = false;
  
  /* Allocate initial CHAIN */
  register void * const data = malloc(CHAIN_SEGMENT_SIZE*(sizeof(off_t) + sizeof(size_t)));
  if (data == NULL) {
    fputs("Unable to allocate sufficient memory.\n", stderr);
    return 1;
  }
  Offset_chain_ptr[0] = (off_t*) data;
  Header_chain_ptr[0] = (size_t*) &((off_t*) data)[CHAIN_SEGMENT_SIZE];
 
  size_t ChainCount    = 1UL;
  size_t SequenceCount = 0UL;
  register off_t * restrict Offsets  = Offset_chain_ptr[0];
  register size_t * restrict Headers = Header_chain_ptr[0];
  
  /* Initial line does not begin with a '\n' */
  off_t BufferOffset  = lseek (FileDescriptor, 0, SEEK_CUR);
  size_t length       = read(FileDescriptor, Buffer, BUFFER_SIZE*sizeof(char));
  _Bool LastWasReturn = true;

  if (length > 0) {
    do {
      if (LastWasReturn) {                                              /* Tricky end of line located on */
                                                                        /* multiple of buffer size       */
        if (HeaderPending) {                                            /* Terminate pending header */
            *Headers = (size_t) ( BufferOffset - 1 - *Offsets );
            if ( ++SequenceCount == CHAIN_SEGMENT_SIZE) {               /* More allocation required or not? */
              if (++ChainCount == MAX_NUMBER_OF_CHAIN) {
                fputs("Chaining limit reached, consider increasing hard coded value.\n", stderr);
                goto FreeMemory;
              }
              {
                register void * data = malloc(CHAIN_SEGMENT_SIZE*(sizeof(off_t) + sizeof(size_t)));
                if (data == NULL) {
                  fputs("Unable to allocate sufficient memory.\n", stderr);
                  goto FreeMemory;
                }
                Offset_chain_ptr[ChainCount-1] = (off_t*) data;
                Header_chain_ptr[ChainCount-1] = (size_t*) &((off_t*) data)[CHAIN_SEGMENT_SIZE];
              }
              Offsets  = Offset_chain_ptr[ChainCount-1];
              Headers  = Header_chain_ptr[ChainCount-1];

              SequenceCount = 0;
            } else {
              ++Headers;
              ++Offsets;
            }
            HeaderPending = false;
        } else if (Buffer[0] == '>' ) {                                 /* New sequence starting */                                    
          *Offsets = BufferOffset;
          HeaderPending  = true;
        }
        LastWasReturn = false;
      }
      
      for (size_t i=1; i<length; ++i) {
        if (Buffer[i-1] == '\n') {                                      /* New Line */
          if (Buffer[i] == '>' ) {                                      /* New sequence starting */
            *Offsets = BufferOffset + (off_t) i;
            HeaderPending  = true;
          } else if (HeaderPending) {                                   /* Terminate pending header */
            *Headers = (size_t) ( BufferOffset + (off_t) i - 1 - *Offsets );
            if ( ++SequenceCount == CHAIN_SEGMENT_SIZE) {               /* More allocation required or not? */
              if (++ChainCount == MAX_NUMBER_OF_CHAIN) {
                fputs("Chaining limit reached, consider increasing hard coded value.\n", stderr);
                goto FreeMemory;
              }
              {
                register void * data = malloc(CHAIN_SEGMENT_SIZE*(sizeof(off_t) + sizeof(size_t)));
                if (data == NULL) {
                  fputs("Unable to allocate sufficient memory.\n", stderr);
                  goto FreeMemory;
                }
                Offset_chain_ptr[ChainCount-1] = (off_t*) data;
                Header_chain_ptr[ChainCount-1] = (size_t*) &((off_t*) data)[CHAIN_SEGMENT_SIZE];
              }
              Offsets  = Offset_chain_ptr[ChainCount-1];
              Headers  = Header_chain_ptr[ChainCount-1];

              SequenceCount = 0;
            } else {
              ++Headers;
              ++Offsets;
            }
            HeaderPending = false;
          }
        }
      }
      if (Buffer[length-1] == '\n') {                                   /* Tricky end of line located on */
        LastWasReturn = true;                                           /* multiple of buffer size       */
      }
      BufferOffset  = lseek (FileDescriptor, 0, SEEK_CUR);
      length = read(FileDescriptor, Buffer, BUFFER_SIZE*sizeof(char));
    } while (length > 0);

    /* Package the overall data */

    SequenceCount += (ChainCount-1)*CHAIN_SEGMENT_SIZE;
    
    DataPtr = (s_Data*) malloc((1+SequenceCount)*sizeof(s_Data));
    if ( DataPtr == NULL) {
      fputs("Unable to allocate sufficient memory\n", stderr);
      goto FreeMemory;
    }


    /* Fill up the data */
    Offsets    = Offset_chain_ptr[0];
    Headers    = Header_chain_ptr[0];
    ChainCount = 0;
    for (size_t i=0; i<SequenceCount-1; ++i) {
      const size_t index = i % CHAIN_SEGMENT_SIZE;
      if ( index != CHAIN_SEGMENT_SIZE-1 ) {
        DataPtr[i].Offset         = Offsets[index];
        DataPtr[i].HeaderLength   = (unsigned int) Headers[index];
				assert(Offsets[index+1] > Offsets[index]);
        register const size_t SequenceLength = Offsets[index+1] - Offsets[index];
        MaxSequenceSize = (SequenceLength > MaxSequenceSize) ? SequenceLength : MaxSequenceSize;
        DataPtr[i].SequenceLength = (unsigned int) (SequenceLength - Headers[index] - 1);
        
      }
      else {
        DataPtr[i].Offset       = Offsets[index];
        DataPtr[i].HeaderLength = (unsigned int) Headers[index];

        free(Offset_chain_ptr[ChainCount++]);
        Offsets = Offset_chain_ptr[ChainCount];
        Headers = Header_chain_ptr[ChainCount];
				assert(Offsets[0] > DataPtr[i].Offset);
        register const size_t SequenceLength = Offsets[0] - DataPtr[i].Offset;
        MaxSequenceSize = (SequenceLength > MaxSequenceSize) ? SequenceLength : MaxSequenceSize;
        DataPtr[i].SequenceLength = (unsigned int) (SequenceLength - DataPtr[i].HeaderLength - 1);
      }
    }

    const size_t index = (SequenceCount-1) % CHAIN_SEGMENT_SIZE;
    
    DataPtr[SequenceCount-1].Offset         = Offsets[index];
    DataPtr[SequenceCount-1].HeaderLength   = (unsigned int) Headers[index];
    register const size_t SequenceLength    = FaS->FileSize - Offsets[index];

    MaxSequenceSize = (SequenceLength > MaxSequenceSize) ? SequenceLength : MaxSequenceSize;
		assert(SequenceLength > (Headers[index] - 1));
    DataPtr[SequenceCount-1].SequenceLength = (unsigned int) (SequenceLength - Headers[index] - 1);
    free(Offsets);
  }
 
  /* Add extra file size at last DataPtr */
  DataPtr[SequenceCount].Offset = FaS->FileSize;
  DataPtr[SequenceCount].HeaderLength = 0;
  DataPtr[SequenceCount].SequenceLength = 0;
  
  FaS->DataPtr         = DataPtr;
  FaS->SequenceCount   = SequenceCount;
  FaS->MaxSequenceSize = MaxSequenceSize+1;
  
  return 0;
  
  FreeMemory:
    for (size_t i=0; i<ChainCount; ++i) {
      free(Offset_chain_ptr[i]);
    }
  return 1;
}

int ExportFASTAStructure(FILE* stream, const FASTAStructure * const FaS)
{
 
  unsigned short FileNameLength = (unsigned short) strlen(FaS->FileName);
  if (fwrite(&(FileNameLength), sizeof(unsigned short), 1, stream) != 1) return 1;
  if (fwrite(FaS->FileName, sizeof(char), (size_t) FileNameLength, stream) != (size_t) FileNameLength) return 1;
  if (fwrite(&(FaS->FileSize), sizeof(off_t), 1, stream) != 1) return 1;
  unsigned long tmp = (unsigned long) FaS->SequenceCount;
  if (fwrite(&tmp, sizeof(unsigned long), 1, stream) != 1) return 1;
  tmp = (unsigned long) FaS->MaxSequenceSize;
  if (fwrite(&tmp, sizeof(unsigned long), 1, stream) != 1) return 1;
  if (fwrite(FaS->DataPtr, sizeof(s_Data), 1+FaS->SequenceCount, stream) != 1+FaS->SequenceCount) return 1;
  
 return 0; 
}

int ImportFASTAStructure(FILE* stream, FASTAStructure * const FaS)
{
  unsigned short FileNameLength;
  unsigned long SequenceCount, MaxSequenceSize; 
  char FileName[256];
  
  
  if (fread(&FileNameLength, sizeof(unsigned short), 1, stream) != 1) return 1;
  if (fread(FileName, sizeof(char), (size_t) FileNameLength, stream) != FileNameLength) return 1;
  if (fread(&(FaS->FileSize), sizeof(off_t), 1, stream) != 1) return 1;
  if (fread(&SequenceCount, sizeof(unsigned long), 1, stream) != 1) return 1;
  FaS->SequenceCount = (size_t) SequenceCount;
  if (fread(&MaxSequenceSize, sizeof(unsigned long), 1, stream) != 1) return 1;
  FaS->MaxSequenceSize = (size_t) MaxSequenceSize;
  FaS->DataPtr = (s_Data*) malloc((1+FaS->SequenceCount)*sizeof(s_Data));
  if ( FaS->DataPtr == NULL) {
    fputs("Unable to allocate sufficient memory to hold FASTA structure\n", stderr);
    return 1;
  }
  if (fread(FaS->DataPtr, sizeof(s_Data), (1+FaS->SequenceCount), stream) != 1+FaS->SequenceCount) return 1;
  
 return 0; 
}

int isFASTA(const char * const FileName) 
{
	const int FileDescriptor = open(FileName, O_RDONLY);
	if ( FileDescriptor == -1) {
    fprintf(stderr, "Error occured while accessing file %s\n", FileName);
    perror("The error reads ");
    return 0;
  }
  
  char buffer;
	int res = 0;
	while (1) {
		if (read(FileDescriptor, &buffer, 1) != 1) {
			fprintf(stderr, "Error reading a character from %s\n", FileName);
			break;
		}
		
		if (buffer == '>') 
			res = 1;
		else if (buffer == ' ' || buffer == '\n' || buffer == '\t')
			continue;
		else 
			break;
	}
	
	close(FileDescriptor);
	return res;
}

int OpenFASTAStructure(const char * const FileName, FASTAStructure * const FaS)
{
	struct timeval _t0, _t1;
	struct stat FileStat;
	int res;
	
	{
		const size_t len = strlen(FileName);
		if (len >= 256) {
			fprintf(stderr, "Filename %s is longer (%lu) than hardcoded storage size (256)\n", FileName, len);
			return 1;
		}
		memcpy(FaS->FileName, FileName, (len+1)*sizeof(char));
	}
	
	const int FileDescriptor = open(FileName, O_RDONLY);
	if ( FileDescriptor == -1) {
    fprintf(stderr, "Error occured while accessing file %s\n", FileName);
    perror("The error reads ");
    return 1;
  }

  if (fstat(FileDescriptor, &FileStat) == -1 ) {
    fprintf(stderr, "Error occured while accessing file %s\n", FileName);
    perror("The error reads ");
    return 1;
  }
  
  if (FileStat.st_size == 0 )
    return 1;
  else
    FaS->FileSize = FileStat.st_size;
	
	gettimeofday(&_t0,0);
	if (!(FaS->Options & DoFastaIndexImport)) {
		res = AnalyzeFASTAStructure(FileDescriptor, FaS);
	}
	else {
		FILE* inIndex = fopen(FaS->indexFileName, "rb");
		if (inIndex != NULL) {
			res = ImportFASTAStructure(inIndex, FaS);
			fclose(inIndex);
		}
		else {
			fprintf(stderr,"Unable to open index file %s, will analyze database instead.\n",FaS->indexFileName);
			res = AnalyzeFASTAStructure(FileDescriptor, FaS);
		}
	}
	gettimeofday(&_t1,0);
	if (res != 0) goto bail;
	
	if (OutputVerbose) {
		const double T = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
		fprintf(stderr, "FASTA file %s indexing took %lf seconds.\n", FaS->FileName, T);
		fprintf(stderr,
					"FASTA file %s analyzed\n"
					"  - Found %lu sequences within %lu bytes\n"
					"  - Biggest sequence entry is %lu bytes\n",
					FaS->FileName, FaS->SequenceCount-1, (size_t) FaS->FileSize, FaS->MaxSequenceSize);
	}

	if ((FaS->Options & DoFastaIndexExport)) {
		FILE *io = fopen(FaS->indexFileName, "wb");
		if ( io != NULL ) {
			const int itmp = ExportFASTAStructure(io, FaS);
			if (OutputVerbose) {
				fprintf(stderr, itmp>0 ? "Export of indices failed, Pac space for %s\n": "Export of indices to file %s\n", FaS->indexFileName);
			}
			fclose(io);
		}
		else {
			if (OutputVerbose)
				fprintf(stderr, "Export of indices failed, Pac write permission for %s\n", FaS->indexFileName);
		}
	}
	
#ifdef PRF_USE_MMAP
	FaS->SequenceFile = (const char * const restrict) mmap(NULL, (size_t) FaS->FileSize, PROT_READ, MAP_PRIVATE, FileDescriptor, 0);
	if (FaS->SequenceFile == NULL) {
		fputs("Unable to map sequence file to memory\n", stderr);
		exit(1);
	}
#endif
	close(FileDescriptor);
	
	return 0;
	
	bail:
	return res;
}

void CloseFASTAStructure(FASTAStructure * const FaS)
{
#ifdef PRF_USE_MMAP
	munmap((void*) FaS->SequenceFile, FaS->FileSize);
	FaS->SequenceFile = NULL;
#endif
	free(FaS->DataPtr);
}

#ifdef _TEST
int main(int argc, char *argv[])
{
  char Buffer[2048] __attribute__((aligned(16)));
  FASTAStructure FASTA;
  
  if (argc < 2 || argc > 3) { fputs("Provide just a FASTA file\n", stderr); return 1; }
  printf("Starting analysis of %s...\n", argv[1]);

  const int res = AnalyzeFASTAStructure(argv[1], &FASTA);
  if (res != 0) {
    fputs("Error found.\n", stderr);
  } else {
    printf("FASTA file %s analyzed\n\tFound %lu sequences within %lu bytes\n\tBiggest sequence entry is %lu bytes\n",
           argv[1],
           FASTA.SequenceCount,
           FASTA.FileSize,
           FASTA.MaxSequenceSize);
    if (argc > 2) {
      const int index = (size_t) atoi(argv[2]);
      if ( index < 0) {
        sprintf(Buffer, "%s.copy\0", argv[1]);
        printf("Recreating global file for diff command in %s\n", Buffer);
        FILE* const out = fopen(Buffer, "w");
        FILE* const in  = fopen(argv[1], "r");
        for (size_t i=0; i<FASTA.SequenceCount; ++i) {
          fseek(in, (long) FASTA.DataPtr[i].Offset, SEEK_SET);
          fread(Buffer, FASTA.DataPtr[i].HeaderLength, sizeof(char), in);
          Buffer[FASTA.DataPtr[i].HeaderLength] = '\0';
          fprintf(out,"%s\n", Buffer);
          if (FASTA.DataPtr[i].SequenceLength >= 2048) {
            fprintf(stderr, "Sorry but allocated buffer cannot hold sequence length of %lu\n",FASTA.DataPtr[i].SequenceLength );
            return 1;
          }
          fread(Buffer, FASTA.DataPtr[i].SequenceLength, sizeof(char), in);
          Buffer[FASTA.DataPtr[i].SequenceLength] = '\0';
          fprintf(out, "%s\n", &Buffer[1]);
        }
        fclose(in);
        fclose(out);
      } else if (index < FASTA.SequenceCount) {
        printf("Sequence %i is at offset %li with header size of %lu and sequence size of %lu\n",
               index,
               FASTA.DataPtr[index].Offset,
               FASTA.DataPtr[index].HeaderLength,
               FASTA.DataPtr[index].SequenceLength);

        FILE* const in = fopen(argv[1], "r");
        fseek(in, (long) FASTA.DataPtr[index].Offset, SEEK_SET);
        fread(Buffer, FASTA.DataPtr[index].HeaderLength, sizeof(char), in);
        Buffer[FASTA.DataPtr[index].HeaderLength] = '\0';
        printf("Header : %s\n", Buffer);
        if (FASTA.DataPtr[index].SequenceLength < 2048) {
          fread(Buffer, FASTA.DataPtr[index].SequenceLength, sizeof(char), in);
          Buffer[FASTA.DataPtr[index].SequenceLength] = '\0';
          printf("Sequence : %s\n", Buffer);
        } else {
          fputs("Sorry but allocated buffer cannot hold sequence length\n", stderr);
        }
        fclose(in);
      } 
    }
    
  }
  
  return 0;
}
#endif
