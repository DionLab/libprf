#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sam.h>
#include <assert.h> 

static const char * const prefix[] = { "@HD", "@SQ", "@RG", "@PG", "@CO" };
static const char * const token[] = { "VN", "SO", "pb" };

void readRG(const char * const text, const size_t len)
{
	const char * restrict ptr = text;
	const uintptr_t Limit = (uintptr_t) text + len - 1;
	while ((uintptr_t) ptr < Limit) {
		// Get next tab
		const char * restrict cptr = ptr;
		while (*cptr != '\t' && *cptr != '\n') ++cptr;
		const int llen = (uintptr_t) cptr - (uintptr_t) ptr;
		printf("\tRG\t%.*s\n", llen, ptr);
		ptr += llen + 1;
	}
	
	//printf("RG: %.*s\n", (int) len, text);
	
}

void readHeader (const bam_hdr_t * const restrict hdr) 
{
	const char * restrict Text = hdr->text;
	const size_t len = hdr->l_text;
	const uintptr_t Limit = (uintptr_t) Text + len - 1;
	
	while ((uintptr_t) Text < Limit) {
		// Get next carriage return 
		const char * restrict ptr = Text;
		while ((uintptr_t) ptr < Limit && *ptr != '\n') ++ptr;
		if ((uintptr_t) ptr == Limit) break;
		
		const size_t lineLen = (uintptr_t) ptr - (uintptr_t) Text;
		
		// skip if line is not long enough to contain true values
    if (lineLen < 5) continue;

		
    if (Text[0] == '@') {
			// Token @HD
			if (Text[1] == 'H' && Text[2] == 'D') {
				const char * restrict Tokens = Text + 3;
				while ((uintptr_t) Tokens < (uintptr_t) ptr) { 
					if (*Tokens != '\t') {
						++Tokens;
					}
					else {
						if (Tokens[1] == 'V' && Tokens[2] == 'N') {
							fputs("Version = ", stdout);
							Tokens += 4;
							while (*Tokens != '\t' && *Tokens != '\n') fputc((int) *Tokens++, stdout);
							fputc('\n', stdout);
						}
						else if (Tokens[1] == 'S' && Tokens[2] == 'O') {
							fputs("Sort = ", stdout);
							Tokens += 4;
							while (*Tokens != '\t' && *Tokens != '\n') fputc((int) *Tokens++, stdout);
							fputc('\n', stdout);
						}
						else if (Tokens[1] == 'p' && Tokens[2] == 'b') {
							fputs("PacBioBamVersion = ", stdout);
							Tokens += 4;
							while (*Tokens != '\t' && *Tokens != '\n') fputc((int) *Tokens++, stdout);
							fputc('\n', stdout);
						}
					}
				}
				// Pac for required tags
// 				if (Version().empty())
// 					Version(std::string(hts_version()));
			}
			// Token @SQ
			else if (Text[1] == 'S' && Text[2] == 'Q') {
//             AddSequence(SequenceInfo::FromSam(line));
			}
			// Token RG
			else if (Text[1] == 'R' && Text[2] == 'G') {
				readRG(Text+4, lineLen-4);
			}
			// Token PG
			else if (Text[1] == 'P' && Text[2] == 'G') {
//             AddProgram(ProgramInfo::FromSam(line));
			}
			// Token CO
			else if (Text[1] == 'C' && Text[2] == 'O') {
//             AddComment(line.substr(4));
			}
    }
    Text = ptr + 1;
	}
}


int main(int argc, char *argv[]){
	
	samFile *fp_in = hts_open(argv[1],"r"); //open bam file
	if (fp_in->format.format != bam) {
		fprintf(stderr, "This is not bam file\n");
		return 1;
	}
	bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
	readHeader(bamHdr);
	bam1_t *aln = bam_init1(); //initialize an alignment
	
	char *chrom = argv[2];
	int locus = atoi(argv[3]);
	int comp ;

	printf("%s\t%d\n", chrom, locus);

	//header parse
	char  *tar = bamHdr->text ;
	uint32_t *tarlen = bamHdr->target_len ;

	printf("%s\n",tar);
	
	if (bam_read1(fp_in->fp.bgzf, aln) > 0) {
		bam_get_aux(aln);
		int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
		char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
		uint32_t len = aln->core.l_qseq; //length of the read.
		
		uint8_t *q = bam_get_seq(aln); //quality string
		uint32_t q2 = aln->core.qual ; //mapping quality
		
		
		char *qseq = (char *)malloc(len);

		for(int i=0; i< len ; i++){
			qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
		}
		
		printf("%d\t%d\t%s\n",pos,len,qseq);

// 		 if(strcmp(chrom, chr) == 0){
// 
// 			 if(locus > pos+len){
// 				 printf("%s\t%d\t%d\t%s\t%s\t%d\n",chr,pos,len,qseq,q,q2);
// 			 }
// 		 }
		assert(0);
	}
	
	bail:
	bam_destroy1(aln);
	sam_close(fp_in);
	
	return 0;
}
