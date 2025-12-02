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
#include "db.h"
#include "hashmapkma.h"
#include "pherror.h"
#include "stdstat.h"

void dbInfo(char *filename) {
	
	const char bases[6] = "ACGTN-";
	int i, filename_len, min, max;
	unsigned *values, *value_index;
	short unsigned *values_s;
	long unsigned ntcount, prefix, ntax, v_index, n, size, index;
	long unsigned *value_index_l;
	double mean, var, tmp;
	char Prefix[33];
	FILE *dbfile;
	HashMapKMA *templates;
	
	/* print DB name */
	fprintf(stdout, "DB:\t%s\n", filename);
	
	/* get # nucleotides */
	filename_len = strlen(filename);
	strcpy(filename + filename_len, ".seq.b");
	dbfile = sfopen(filename, "rb");
	sfseek(dbfile, 0, SEEK_END);
	ntcount = 4 * ftell(dbfile);
	fclose(dbfile);
	filename[filename_len] = 0;
	fprintf(stdout, "# nucleotides:\t%lu\n", ntcount);
	
	/* load database */
	strcpy(filename + filename_len, ".comp.b");
	dbfile = sfopen(filename, "rb" );
	templates = smalloc(sizeof(HashMapKMA));
	if(hashMapKMA_load(templates, dbfile, filename) == 1) {
		fprintf(stderr, "Wrong format of DB.\n");
		exit(1);
	}
	fclose(dbfile);
	filename[filename_len] = 0;
	
	/* get basic statistics */
	fprintf(stdout, "# templates:\t%d\n", templates->DB_size - 1);
	fprintf(stdout, "k:\t%d\n", templates->kmersize);
	fprintf(stdout, "m:\t%d\n", templates->mlen);
	fprintf(stdout, "hc:\t%d\n", templates->flag & 1);
	fflush(stdout);
	/* get prefix */
	if(templates->prefix_len) {
		prefix = templates->prefix;
		i = templates->prefix_len;
		Prefix[i] = 0;
		while(i--) {
			Prefix[i] = bases[prefix & 3];
			prefix >>= 2;
		}
		fprintf(stdout, "prefix:\t%s\n", Prefix);
	} else if(templates->prefix != 0) {
		fprintf(stdout, "prefix:\t-\n");
	}
	fprintf(stdout, "# uniq k-mers:\t%lu\n", templates->n);
	fprintf(stdout, "k-mer fraction covered:\t%f\n", templates->n / power(4, templates->kmersize));
	fprintf(stdout, "inferred tax size:\t%lu\n", templates->v_index);
	fflush(stdout);
	
	/* get number of unique inferred tax */
	ntax = 0;
	v_index = templates->v_index;
	if(templates->DB_size < USHRT_MAX) {
		values_s = templates->values_s;
		while(v_index) {
			++ntax;
			i = *values_s + 1;
			v_index -= i;
			values_s += i;
		}
	} else {
		values = templates->values;
		while(v_index) {
			++ntax;
			i = *values + 1;
			v_index -= i;
			values += i;
		}
	}
	fprintf(stdout, "# inferred taxids:\t%lu\n", ntax);
	fflush(stdout);
	
	/* get min, max, mean and variance of k-mer uniqueness */
	if((templates->size - 1) == templates->mask) {
		if(templates->v_index <= UINT_MAX) {
			value_index = templates->exist;
			value_index_l = 0;
		} else {
			value_index = 0;
			value_index_l = templates->exist_l;
		}
	} else {
		if(templates->n <= UINT_MAX) {
			value_index = templates->value_index;
			value_index_l = 0;
		} else {
			value_index = 0;
			value_index_l = templates->value_index_l;
		}
	}
	if(templates->DB_size < USHRT_MAX) {
		values = 0;
		values_s = templates->values_s;
	} else {
		values = templates->values;
		values_s = 0;
	}
	n = templates->n;
	size = n;
	min = templates->DB_size;
	max = 1;
	mean = 0;
	var = 0;
	while(size) {
		index = value_index ? *value_index++ : *value_index_l++;
		if(index != 1) {
			index = values ? values[index] : values_s[index];
			if(index < min) {
				min = index;
			}
			if(max < index) {
				max = index;
			}
			mean += index;
			
			/* avoid overflow on var */
			tmp = index * index;
			tmp /= n;
			var += tmp;
			
			--size;
		}
	}
	mean /= n;
	var -= mean * mean;
	
	fprintf(stdout, "k-mer co-occurence var:\t%f\n", var);
	fprintf(stdout, "k-mer co-occurence mean:\t%f\n", mean);
	fprintf(stdout, "k-mer co-occurence min:\t%d\n", min);
	fprintf(stdout, "k-mer co-occurence max:\t%d\n", max);
	fflush(stdout);
}

static void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# KMA db gives statistics on a KMA database\n");
	fprintf(helpOut, "# Options are:\t\tDesc:\t\t\t\t\tRequirements:\n");
	fprintf(helpOut, "#\n");
	fprintf(helpOut, "#\t-t_db\t\tTemplate DB\t\t\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-h\t\tShows this help message\n");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int db_main(int argc, char *argv[]) {
	
	unsigned args;
	char *filename;
	
	/* set defaults */
	filename = 0;
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-t_db") == 0) {
			if(++args < argc) {
				filename = smalloc(strlen(argv[args]) + 64);
				strcpy(filename, argv[args]);
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
	
	if(!filename) {
		fprintf(stderr, "Insuffient amount of arguments handed!!!\n");
		helpMessage(1);
	}
	
	dbInfo(filename);
	
	return 0;
}
