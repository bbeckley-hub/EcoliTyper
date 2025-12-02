/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *		http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#define _XOPEN_SOURCE 600
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "penalties.h"
#include "pherror.h"
#include "qseqs.h"
#include "runkma.h"
#include "seq2fasta.h"
#include "stdnuc.h"

int * getLengths(char *filename) {
	
	int file_len, DB_size, *template_lengths;
	FILE *infile;
	
	file_len = strlen(filename);
	strcat(filename, ".length.b");
	infile = sfopen(filename, "rb");
	filename[file_len] = 0;
	sfread(&DB_size, sizeof(int), 1, infile);
	template_lengths = calloc(DB_size, sizeof(int));
	if(!template_lengths) {
		ERROR();
	}
	sfread(template_lengths, sizeof(int), DB_size, infile);
	*template_lengths = DB_size;
	fclose(infile);
	
	return template_lengths;
}

void printFastas(char *filename, int *template_lengths) {
	
	const char bases[6] = "ACGTN-";
	int i, j, max, DB_size, file_len;
	long unsigned *compseq;
	char *seq;
	FILE *seqfile, *namefile;
	Qseqs *template_name;
	
	template_name = setQseqs(256);
	DB_size = *template_lengths;
	max = template_lengths[1];
	for(i = 2; i < DB_size; i++) {
		if(max < template_lengths[i]) {
			max = template_lengths[i];
		}
	}
	
	file_len = strlen(filename);
	strcat(filename, ".seq.b");
	seqfile = sfopen(filename, "rb");
	filename[file_len] = 0;
	strcat(filename, ".name");
	namefile = sfopen(filename, "rb");
	filename[file_len] = 0;
	
	/* Allocate stuff */
	seq = malloc(max+2);
	compseq = calloc((max >> 5) + 1, sizeof(long unsigned));
	if(!seq || !compseq) {
		ERROR();
	}
	for(i = 1; i < DB_size; ++i) {
		sfread(compseq, sizeof(long unsigned), (template_lengths[i] >> 5) + 1, seqfile);
		
		j = template_lengths[i];
		*(seq += j) = '\n';
		while(j--) {
			*--seq = bases[getNuc(compseq, j)];
		}
		
		fprintf(stdout, ">%s\n", nameLoad(template_name, namefile));
		cfwrite(seq, 1, template_lengths[i] + 1, stdout);
	}
}

int intCmpAscend(const void * a, const void * b) {
	return *((int*) a) - *((int*) b);
}

void printFastaList(char *filename, int *template_lengths, int *seqlist) {
	
	const char bases[6] = "ACGTN-";
	int i, j, n, max, DB_size, file_len;
	long unsigned skip, *compseq;
	char *seq;
	FILE *seqfile, *namefile;
	Qseqs *template_name;
	
	/* sort sequence list */
	qsort(seqlist + 1, (n = *seqlist), sizeof(int), intCmpAscend);
	
	/* skip invalid templates */
	while(n && *++seqlist <= 0) {
		--n;
	}
	
	/* open index */
	template_name = setQseqs(256);
	DB_size = *template_lengths;
	max = template_lengths[1];
	for(i = 2; i < DB_size; i++) {
		if(max < template_lengths[i]) {
			max = template_lengths[i];
		}
	}
	file_len = strlen(filename);
	strcat(filename, ".seq.b");
	seqfile = sfopen(filename, "rb");
	filename[file_len] = 0;
	strcat(filename, ".name");
	namefile = sfopen(filename, "rb");
	filename[file_len] = 0;
	
	/* Allocate stuff */
	seq = malloc(max+2);
	compseq = calloc((max >> 5) + 1, sizeof(long unsigned));
	if(!seq || !compseq) {
		ERROR();
	}
	skip = 0;
	for(i = 1; n && i < DB_size; ++i) {
		if(i == *seqlist) {
			/* get seq */
			sfseek(seqfile, skip, SEEK_CUR);
			sfread(compseq, sizeof(long unsigned), (template_lengths[i] >> 5) + 1, seqfile);
			
			j = template_lengths[i];
			*(seq += j) = '\n';
			while(j--) {
				*--seq = bases[getNuc(compseq, j)];
			}
			
			/* print seq */
			fprintf(stdout, ">%s\n", nameLoad(template_name, namefile));
			sfwrite(seq, 1, template_lengths[i] + 1, stdout);
			
			/* get next target */
			while(--n && i == *++seqlist);
			skip = 0;
		} else {
			nameSkip(namefile, max);
			skip += (((template_lengths[i] >> 5) + 1) * sizeof(long unsigned));
		}
	}
}

int * intSplit(char sep, const char *src) {
	
	int n, *seqlist;
	char *stringPtr, *intStr, *errStr;
	
	/* get number of elements */
	n = 2;
	stringPtr = (char *) src;
	while(*stringPtr) {
		if(*stringPtr++ == sep) {
			++n;
		}
	}
	seqlist = smalloc(n * sizeof(int));
	
	/* get elements */
	n = 0;
	stringPtr = (char *) src;
	intStr = stringPtr;
	while(*stringPtr) {
		if(*++stringPtr == sep) {
			*stringPtr = 0;
			++n;
			seqlist[n] = strtol(intStr, &errStr, 10);
			if(*errStr != 0 || *(intStr = ++stringPtr) == 0) {
				fprintf(stderr, "Invalid list parsed.\n");
				exit(1);
			}
		} else if(*stringPtr == 0) {
			++n;
			seqlist[n] = strtol(intStr, &errStr, 10);
			if(*errStr != 0) {
				fprintf(stderr, "Invalid list parsed.\n");
				exit(1);
			}
		}
	}
	*seqlist = n;
	
	return seqlist;
}

static void helpMessage(int status) {
	
	FILE *out;
	
	if(status) {
		out = stderr;
	} else {
		out = stdout;
	}
	fprintf(out, "kma seq2fasta prints the fasta sequence of a given kma index to stdout.\n");
	fprintf(out, "# Options are:\tDesc:\t\t\t\t\tDefault:\tRequirements:\n");
	fprintf(out, "#\t-t_db\tTemplate DB\t\t\t\tNone\t\tREQUIRED\n");
	fprintf(out, "#\t-seqs\tComma separated list of templates\tPrint entire index.\n");
	fprintf(out, "#\t-h\tShows this help message\n");
	exit(status);
}

int seq2fasta_main(int argc, char *argv[]) {
	
	int args, file_len, *template_lengths, *seqlist;
	char *filename;
	
	seqlist = 0;
	file_len = 0;
	filename = 0;
	args = 0;
	while(++args < argc) {
		if(strcmp(argv[args], "-t_db") == 0) {
			if(++args < argc) {
				file_len = strlen(argv[args]);
				filename = smalloc(file_len + 64);
				strcpy(filename, argv[args]);
			}
		} else if(strcmp(argv[args], "-seqs") == 0) {
			if(++args < argc) {
				seqlist = intSplit(',', argv[args]);
			}
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else {
			helpMessage(1);
		}
	}
	if(!filename) {
		fprintf(stderr, "Need a db\n");
		helpMessage(1);
	}
	
	/* get lengths */
	template_lengths = getLengths(filename);
	
	if(seqlist) {
		/* get sequences from list */
		printFastaList(filename, template_lengths, seqlist);
	} else {
		/* get all sequences */
		printFastas(filename, template_lengths);
	}
	
	return 0;
}
