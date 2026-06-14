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
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hashmapkma.h"
#include "pherror.h"
#include "stdnuc.h"
#include "update.h"

unsigned convertLength_014to015(char *filename) {
	
	unsigned size, file_len;
	int *lengths;
	FILE *file;
	
	file_len = strlen(filename);
	strcat(filename, ".length.b");
	file = sfopen(filename, "rb+");
	filename[file_len] = 0;
	
	sfread(&size, sizeof(unsigned), 1, file);
	lengths = smalloc(3 * size * sizeof(unsigned));
	
	file_len = fread(lengths, sizeof(unsigned), 3 * size, file);
	sfseek(file, sizeof(unsigned), SEEK_SET);
	if(file_len == size) {
		file_len = 0;
	} else if(file_len == 2 * size) {
		fprintf(stderr, "DB is old.\n");
		fprintf(stderr, "It will only work for \"-Sparse\" mapping!!!\n");
		sfwrite(lengths, sizeof(unsigned), size, file);
		sfwrite(lengths, sizeof(unsigned), 2 * size, file);
		file_len = 0;
	} else if(file_len == 3 * size) {
		sfwrite(lengths + 2 * size, sizeof(unsigned), size, file);
		sfwrite(lengths, sizeof(unsigned), 2 * size, file);
		file_len = 1;
	} else {
		fprintf(stderr, "DB is malformed.\n");
		exit(1);
	}
	
	fclose(file);
	return file_len;
}

int hashMapKMA_014to015(char *filename, unsigned prefix) {
	
	unsigned i, tmp, size, kmersize, seqsize, file_len, shifter, DB_size;
	long unsigned mask, *seq;
	FILE *file;
	HashMapKMA *dest;
	
	/* rm filename.b */
	file_len = strlen(filename);
	strcat(filename, ".b");
	remove(filename);
	filename[file_len] = 0;
	
	/* load DB */
	strcat(filename, ".comp.b");
	file = sfopen(filename, "rb");
	filename[file_len] = 0;
	
	/* load sizes */
	dest = smalloc(sizeof(HashMapKMA));
	dest->n = 0;
	dest->size = 0;
	sfread(&DB_size, sizeof(unsigned), 1, file);
	sfread(&dest->kmersize, sizeof(unsigned), 1, file);
	sfread(&dest->prefix_len, sizeof(unsigned), 1, file);
	sfread(&dest->prefix, sizeof(long unsigned), 1, file);
	sfread(&dest->size, sizeof(long unsigned), 1, file);
	
	kmersize = dest->kmersize;
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	
	/* load changed size */
	sfread(&tmp, sizeof(unsigned), 1, file);
	dest->n = tmp;
	sfread(&seqsize, sizeof(unsigned), 1, file); //seq size
	sfread(&tmp, sizeof(unsigned), 1, file);
	dest->v_index = tmp;
	sfread(&tmp, sizeof(unsigned), 1, file);
	dest->null_index = tmp;
	
	/* make checks */
	if(dest->size < dest->n || dest->n == 0) {
		fprintf(stderr, "DB is not of version 0.14\n");
		exit(1);
	}
	
	/* load arrays */
	dest->exist = smalloc(dest->size * sizeof(unsigned));
	if(dest->size != fread(dest->exist, sizeof(unsigned), dest->size, file)) {
		return 1;
	}
	
	if(mask != (dest->size - 1)) {
		seq = smalloc(seqsize * sizeof(long unsigned));
		if(seqsize != fread(seq, sizeof(long unsigned), seqsize, file)) {
			return 1;
		}
	} else {
		seq = 0;
	}
	
	dest->values = smalloc(dest->v_index * sizeof(int));
	if(dest->v_index != fread(dest->values, sizeof(int), dest->v_index, file)) {
		return 1;
	}
	
	if(mask != (dest->size - 1)) {
		dest->key_index = smalloc((dest->n + 1) * sizeof(unsigned));
		if((dest->n + 1) != fread(dest->key_index, sizeof(unsigned), dest->n + 1, file)) {
			return 1;	
		}
		
		dest->value_index = smalloc(dest->n * sizeof(unsigned));
		if(dest->n != fread(dest->value_index, sizeof(unsigned), dest->n, file)) {
			return 1;
		}
	}
	/* convert to new format */
	/* change prefix if sparse - */
	if(prefix && dest->prefix_len == 0) {
		dest->prefix = 1;
	}
	
	strcat(filename, ".comp.b");
	file = sfopen(filename, "wb");
	filename[file_len] = 0;
	sfwrite(&DB_size, sizeof(unsigned), 1, file);
	sfwrite(&dest->kmersize, sizeof(unsigned), 1, file);
	sfwrite(&dest->prefix_len, sizeof(unsigned), 1, file);
	sfwrite(&dest->prefix, sizeof(long unsigned), 1, file);
	sfwrite(&dest->size, sizeof(long unsigned), 1, file);
	sfwrite(&dest->n, sizeof(long unsigned), 1, file);
	sfwrite(&dest->v_index, sizeof(long unsigned), 1, file);
	sfwrite(&dest->null_index, sizeof(long unsigned), 1, file);
	
	/* exist */
	cfwrite(dest->exist, sizeof(unsigned), dest->size, file);
	
	/* values */
	if(DB_size < USHRT_MAX) {
		dest->values_s = (short unsigned *)(dest->values);
		for(i = 0; i < dest->v_index; ++i) {
			dest->values_s[i] = dest->values[i];
		}
		size = sizeof(short unsigned);
	} else {
		size = sizeof(unsigned);
	}
	cfwrite(dest->values, size, dest->v_index, file);
	free(dest->values);
	
	if(mask == (dest->size - 1)) {
		return 0;
	}
	
	if(dest->kmersize <= 16) {
		cfwrite(dest->key_index, sizeof(unsigned), dest->n + 1, file);
	} else {
		dest->key_index_l = realloc(dest->key_index, (dest->n + 1) * sizeof(long unsigned));
		if(dest->key_index_l) {
			dest->key_index = (unsigned *)(dest->key_index_l);
		} else {
			ERROR();
		}
		
		i = dest->n + 1;
		while(i--) {
			dest->key_index_l[i] = getKmer(seq, dest->key_index[i], shifter);
		}
		cfwrite(dest->key_index_l, sizeof(long unsigned), dest->n + 1, file);
		free(seq);
	}
	free(dest->key_index);
	
	/* value_index */
	cfwrite(dest->value_index, sizeof(unsigned), dest->n, file);
	
	return 0;
}

int index_014to015(char *filename) {
	
	unsigned prefix, file_len, returner;
	FILE *file;
	
	file_len = strlen(filename);
	
	/* change prefix if sparse - */
	prefix = convertLength_014to015(filename);
	
	returner = hashMapKMA_014to015(filename, prefix);
	
	/* check for deCon */
	strcat(filename, ".decon.b");
	file = fopen(filename, "rb");
	if(file) {
		fclose(file);
		/* change filename to: "filename.decon" */
		filename[file_len + 6] = 0;
		returner += hashMapKMA_014to015(filename, prefix);
	}
	filename[file_len] = 0;
	
	return returner;
}

static void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# KMA_update synchronises kma-indexes to the needed version.\n");
	fprintf(helpOut, "# Options are:\t\tDesc:\t\t\t\t\tRequirements:\n");
	fprintf(helpOut, "#\n");
	fprintf(helpOut, "#\t-t_db\t\tTemplate DB\t\t\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-v\t\t[XXYY], from version major version XX\n#\t\t\tto major version YY. Use minor version,\n#\t\t\tif major version is 0.\t\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-h\t\tShows this help message\n");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int update_main(int argc, char *argv[]) {
	
	unsigned args, version;
	char *filename, *error;
	
	/* set defaults */
	filename = 0;
	version = 0;
	
	args = 1;
	while(args < argc) {
		
		if(strcmp(argv[args], "-t_db") == 0) {
			if(++args < argc) {
				filename = smalloc(strlen(argv[args]) + 64);
				strcpy(filename, argv[args]);
			}
		} else if(strcmp(argv[args], "-v") == 0) {
			if(++args < argc) {
				version = strtoul(argv[args], &error, 10);
				if(*error != 0) {
					fprintf(stderr, " Invalid version specified.\n");
					exit(1);
				}
			}
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else {
			fprintf(stderr, " Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, " Printing help message:\n");
			helpMessage(1);
		}
		++args;
	}
	
	if(!filename || !version) {
		fprintf(stderr, "Insuffient amount of arguments handed!!!\n");
	} else if(version == 1415) {
		if(index_014to015(filename)) {
			fprintf(stderr, "Conversion error.\n");
			exit(1);
		}
	} else {
		fprintf(stderr, "Invalid version swifting specified.\n");
		fprintf(stderr, "Valid conversions:\n");
		fprintf(stderr, "\t%d\t%.2f -> %.2f\n", 1415, 0.14, 0.15);
		return 2;
	}
	
	return 0;
}
