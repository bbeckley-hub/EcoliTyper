/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#define _XOPEN_SOURCE 600
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compdna.h"
#include "filebuff.h"
#include "hashmap.h"
#include "makeindex.h"
#include "pherror.h"
#include "qseqs.h"
#include "seqparse.h"
#include "updateindex.h"

int (*biasPrintPtr)(FILE*, char*, unsigned char*, int) = &biasPrint;
int (*qualcheck)(CompDNA *, const int) = &lenCheck;

int biasPrint(FILE *name_out, char *format, unsigned char *name, int bias) {
	return fprintf(name_out, "%s B%d\n", name, bias);
}

int biasNoPrint(FILE *name_out, char *format, unsigned char *name, int bias) {
	return fprintf(name_out, "%s\n", name);
}

int lenCheck(CompDNA *qseq, const int minlen) {
	return minlen < qseq->seqlen;
}

int internalStopCheck1(CompDNA *qseq) {
	
	int pos, seqlen;
	unsigned codon;
	long unsigned *seq;
	
	if(qseq->seqlen % 3) {
		return 0;
	}
	
	/* valid stop: TAA, TAG and TGA  = 48, 50, 56 */
	seq = qseq->seq;
	seqlen = qseq->seqlen - 3;
	pos = 0;
	while(pos < seqlen) {
		/* get codon */
		if(((pos + 3) & 15) < 3) {
			codon = (((seq[(pos >> 5)] << ((pos & 31) << 1)) | (seq[(pos >> 5) + 1] >> (64 - ((pos & 31) << 1)))) >> 58);
		} else {
			codon = (seq[(pos >> 5)] << ((pos & 31) << 1)) >> 58;
		}
		
		/* check for stop codon */
		if(codon == 48 || codon == 50 || codon == 56) {
			return 0;
		}
		
		pos += 3;
	}
	
	return 1;
}

int internalStopCheck(CompDNA *qseq, const int minlen) {
	
	int pos, seqlen;
	unsigned codon, frames, frames_r;
	long unsigned *seq;
	
	if(qseq->seqlen < minlen) {
		return 0;
	}
	
	/* valid stop: TAA, TAG and TGA  = 48, 50, 56 */
	/* rc: TTA, CTA and TCA = 60, 28, 52*/
	frames = 0;
	frames_r = 0;
	seq = qseq->seq;
	seqlen = qseq->seqlen - 4;
	pos = 1;
	while(pos < seqlen) {
		/* get codon */
		if(((pos + 3) & 15) < 3) {
			codon = (((seq[(pos >> 5)] << ((pos & 31) << 1)) | (seq[(pos >> 5) + 1] >> (64 - ((pos & 31) << 1)))) >> 58);
		} else {
			codon = (seq[(pos >> 5)] << ((pos & 31) << 1)) >> 58;
		}
		
		/* check for stop codon */
		if(codon == 48 || codon == 50 || codon == 56) { /* forward */
			frames |= (1 << (pos % 3));
		} else if(codon == 60 || codon == 28 || codon == 52) {
			frames_r |= (1 << (pos % 3));
		}
		
		/* check if all frames has internal stop codons */
		if(frames == 7 && frames_r == 7) {
			return 0;
		}
		
		++pos;
	}
	
	/* fit to accepted frame */
	/*
	if(frames == 7) {
		comp_rc(qseq);
	}
	*/
	
	return 1;
}

int qualCheck(CompDNA *qseq, const int minlen) {
	
	unsigned start, stop, pos;
	
	if(qseq->seqlen < minlen || qseq->seqlen % 3) {
		return 0;
	}
	
	/* valid start: NTG and ATT */
	/* NTG & 15 = 0X1110 = 14 || ATT = 0x001111 = 15 */
	/* rc: CAN and AAT */
	/* CAN >> 2 = 4 || AAT = 3 */
	start = *(qseq->seq) >> 58;
	
	/* valid stop: TAA, TAG and TGA  = 48, 50, 56 */
	/* rc: TTA, CTA and TCA = 60, 28, 52*/
	pos = qseq->seqlen - 3;
	if((qseq->seqlen & 15) < 3) {
		stop = (((qseq->seq[(pos >> 5)] << ((pos & 31) << 1)) | (qseq->seq[(pos >> 5) + 1] >> (64 - ((pos & 31) << 1)))) >> 58);
	} else {
		stop = (qseq->seq[(pos >> 5)] << ((pos & 31) << 1)) >> 58;
	}
	
	/* check valid start and stop codon */
	if(((start & 15) == 14 || start == 15) && (stop == 48 || stop == 50 || stop == 56)) { /* forward */
		return internalStopCheck1(qseq);
	} else if(((stop >> 2) == 4 || stop == 3) && (start == 60 || start == 28 || start == 52)) { /* rc */
		comp_rc(qseq);
		return internalStopCheck1(qseq);
	}
	
	return 0;
}

void makeDB(HashMap *templates, int kmerindex, char **inputfiles, int fileCount, char *outputfilename, int appender, char *trans, int MinLen, int MinKlen, double homQ, double homT, unsigned **template_lengths, unsigned **template_ulengths, unsigned **template_slengths) {
	
	int fileCounter, file_len, bias, FASTQ;
	char *filename;
	unsigned char *seq;
	FILE *seq_out, *length_out, *name_out;
	Qseqs *header, *qseq;
	FileBuff *inputfile;
	CompDNA *compressor;
	
	/* allocate */
	compressor = smalloc(sizeof(CompDNA));
	allocComp(compressor, 1024);
	header = setQseqs(1024);
	qseq = setQseqs(1024);
	inputfile = setFileBuff(1024 * 1024);
	
	/* open files */
	file_len = strlen(outputfilename);
	strcat(outputfilename, 	".length.b");
	length_out = sfopen(outputfilename, "wb");
	outputfilename[file_len] = 0;
	if(appender) {
		strcat(outputfilename, 	".name");
		name_out = sfopen(outputfilename, "ab");
		outputfilename[file_len] = 0;
	} else {
		strcat(outputfilename, 	".name");
		name_out = sfopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
	}
	if(appender) {
		strcat(outputfilename, ".seq.b");
		seq_out = sfopen(outputfilename, "ab");
		outputfilename[file_len] = 0;
	} else {
		strcat(outputfilename, ".seq.b");
		seq_out = sfopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
	}
	
	fprintf(stderr, "# Updating DBs\n");
	/* iterate inputfiles */
	for(fileCounter = 0; fileCounter < fileCount; ++fileCounter) {
		/* open file */
		filename = inputfiles[fileCounter];
		/* determine filetype and open it */
		if((FASTQ = openAndDetermine(inputfile, filename)) & 2) {
			fprintf(stderr, "%s\t%s\n", "# Reading inputfile: ", filename);
			
			/* parse the file */
			while(FileBuffgetFsa(inputfile, header, qseq, trans)) {
				if(qseq->len >= compressor->size) {
					freeComp(compressor);
					allocComp(compressor, qseq->size);
				}
				bias = compDNAref(compressor, qseq->seq, qseq->len);
				
				if(qualcheck(compressor, MinLen) && update_DB(templates, compressor, templates->DB_size, MinKlen, homQ, homT, *template_ulengths, *template_slengths, header)) {
					/* Update annots */
					seq = header->seq + header->len;
					while(isspace(*--seq)) {
						*seq = 0;
					}
					if(bias > 0) {
						biasPrintPtr(name_out, "%s B%d\n", header->seq + 1, bias);
					} else {
						fprintf(name_out, "%s\n", header->seq + 1);
					}
					updateAnnotsPtr(compressor, templates->DB_size, kmerindex, seq_out, template_lengths, template_ulengths, template_slengths);
					
					fprintf(stderr, "# Added:\t%s\n", header->seq + 1);
					
					if(++(templates->DB_size) == USHRT_MAX) {
						/* convert values to unsigned */
						convertToU(templates);
					}
				} else {
					fprintf(stderr, "# Skipped:\t%s\n", header->seq + 1);
				}
			}
		} else {
			fprintf(stderr, "Unsupported format for file:\t%s\n", filename);
			errno |= 1;
			//exit(1);
		}
		
		/* close file buffer */
		if(FASTQ & 4) {
			gzcloseFileBuff(inputfile);
		} else {
			closeFileBuff(inputfile);
		}
	}
	
	/* Dump annots */
	cfwrite(&templates->DB_size, sizeof(int), 1, length_out);
	if(*template_ulengths != 0) {
		**template_ulengths = 0;
		**template_slengths = 0;
		cfwrite(*template_lengths, sizeof(unsigned), templates->DB_size, length_out);
		cfwrite(*template_slengths, sizeof(unsigned), templates->DB_size, length_out);
		cfwrite(*template_ulengths, sizeof(unsigned), templates->DB_size, length_out);
	} else {
		**template_lengths = kmerindex;
		cfwrite(*template_lengths, sizeof(unsigned), templates->DB_size, length_out);
	}
	fclose(seq_out);
	fclose(length_out);
	fclose(name_out);
	
	if(templates->n) {
		fprintf(stderr, "# Templates key-value pairs:\t%lu.\n", templates->n);// / 1048576);
	} else {
		fprintf(stderr, "DB is empty!!!\n");
		exit((errno | 1));
	}
	
	/* clean */
	freeComp(compressor);
	free(compressor);
	destroyQseqs(header);
	destroyQseqs(qseq);
	destroyFileBuff(inputfile);
}
