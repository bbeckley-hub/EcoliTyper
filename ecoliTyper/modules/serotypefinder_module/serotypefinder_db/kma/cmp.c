/* Philip T.L.C. Clausen Mar 2025 plan@dtu.dk */

/*
 * Copyright (c) 2025, Philip Clausen, Technical University of Denmark
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
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "cmp.h"
#include "hashmapkma.h"
#include "kmmap.h"
#include "pherror.h"
#include "tmp.h"
#include "version.h"

int signaturecmp(unsigned *v1, unsigned *v2, size_t type) {
	
	unsigned n;
	short unsigned *s1, *s2;
	
	if(!v1 || !v2) {
		return 1;
	} else if(type == sizeof(unsigned)) {
		n = *v1-- + 2;
		--v2;
		while(--n) {
			if(*++v1 != *++v2) {
				return 1;
			}
		}
	} else {
		s1 = (short unsigned *)(v1);
		s2 = (short unsigned *)(v2);
		n = *s1-- + 2;
		--s2;
		while(--n) {
			if(*++s1 != *++s2) {
				return 1;
			}
		}
	}
	
	return 0;
}

void printsignature(unsigned *v, size_t type) {
	
	unsigned n;
	short unsigned *s;
	
	if(!v) {
		fprintf(stderr, "%d", 0);
	} else if(type == sizeof(unsigned)) {
		n = *v;
		fprintf(stderr, "%d", n++);
		while(--n) {
			fprintf(stderr, ", %d", *++v);
		}
	} else {
		s = (short unsigned *)(v);
		n = *s;
		fprintf(stderr, "%d", n++);
		while(--n) {
			fprintf(stderr, ", %d", *++s);
		}
	}
	fprintf(stderr, "\n");
}

void printmismatch(long unsigned kmer, unsigned kmersize, unsigned *v1, unsigned *v2, size_t type) {
	
	int i;
	char Kmer[kmersize + 1], bases[4] = "ACGT";
	
	fprintf(stderr, "Kmer:\t%lu\n", kmer);
	i = kmersize;
	Kmer[i] = 0;
	while(i--) {
		Kmer[i] = bases[kmer & 3];
		kmer >>= 2;
	}
	fprintf(stderr, "%s\n", Kmer);
	printsignature(v1, type);
	printsignature(v2, type);
}

int hashMapKMA_cmp(HashMapKMA *t1, HashMapKMA *t2) {
	
	short unsigned *values_s;
	unsigned type, kmersize, *exist, *values, *v_index, *v1, *v2;
	long unsigned kmer, size, null_index, index, *exist_l, *v_index_l;
	
	/* cmp base structure, and set pointers */
	if(t1->n != t2->n || t1->v_index != t2->v_index || t1->mlen != t2->mlen || t1->kmersize != t2->kmersize || t1->flag != t2->flag || t1->prefix_len != t2->prefix_len || t1->prefix || t2->prefix || t1->DB_size != t2->DB_size) {
		fprintf(stderr, "n:\t%lu, %lu\n", t1->n, t2->n);
		fprintf(stderr, "v_index:\t%lu, %lu\n", t1->v_index, t2->v_index);
		fprintf(stderr, "mlen:\t%d, %d\n", t1->mlen, t2->mlen);
		fprintf(stderr, "kmersize:\t%d, %d\n", t1->kmersize, t2->kmersize);
		fprintf(stderr, "flag:\t%d, %d\n", t1->flag, t2->flag);
		fprintf(stderr, "prefix_len:\t%d, %d\n", t1->prefix_len, t2->prefix_len);
		fprintf(stderr, "prefix:\t%lu, %lu\n", t1->prefix, t2->prefix);
		fprintf(stderr, "DB_size:\t%d, %d\n", t1->DB_size, t2->DB_size);
		return 1;
	} else if(t1->DB_size < USHRT_MAX) {
		type = sizeof(short unsigned);
		values = 0;
		values_s = t1->values_s;
	} else {
		type = sizeof(unsigned);
		values = t1->values;
		values_s = 0;
	}
	kmersize = t1->mlen;
	
	if(t1->size == (t1->mask + 1)) { /* direct on t1 */
		if(t1->v_index <= UINT_MAX) {
			exist = t1->exist - 1;
			exist_l = 0;
		} else {
			exist = 0;
			exist_l = t1->exist_l - 1;
		}
		null_index = t1->null_index;
		size = t1->size + 1;
		kmer = 0;
		while(--size) {
			index = exist ? *++exist : *++exist_l;
			if(index != null_index) {
				v1 = values ? (values + index) : ((unsigned *)(values_s + index));
				v2 = hashMap_get(t2, kmer);
				if(signaturecmp(v1, v2, type)) {
					printmismatch(kmer, kmersize, v1, v2, type);
					return 1;
				}
			}
			++kmer;
		}
	} else { /* hashmap on t1 */
		if(t1->mlen <= 16) {
			exist = t1->key_index - 1;
			exist_l = 0;
		} else {
			exist = 0;
			exist_l = t1->key_index_l - 1;
		}
		if(t1->v_index < UINT_MAX) {
			v_index = t1->value_index - 1;
			v_index_l = 0;
		} else {
			v_index = 0;
			v_index_l = t1->value_index_l - 1;
		}
		size = t1->n + 1;
		while(--size) {
			kmer = exist ? *++exist : *++exist_l;
			index = v_index ? *++v_index : *++v_index_l;
			v1 = values ? (values + index) : ((unsigned *)(values_s + index));
			v2 = hashMap_get(t2, kmer);
			if(signaturecmp(v1, v2, type)) {
				printmismatch(kmer, kmersize, v1, v2, type);
				return 1;
			}
		}
	}
	
	return 0;
}

static void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# kma cmp compare two indexed kma databases.\n");
	fprintf(helpOut, "# %16s\t%-32s\t%s\n", "Options:", "Desc:", "Default:");
	fprintf(helpOut, "# %16s\t%-32s\t%s\n", "-t_db", "DB to compare to", "");
	fprintf(helpOut, "# %16s\t%-32s\t%s\n", "-s_db", "DB to compare with", "");
	fprintf(helpOut, "# %16s\t%-32s\t%s\n", "-tmp", "Set directory for temporary files", "");
	fprintf(helpOut, "# %16s\t%-32s\t%s\n", "-v", "Version", "");
	fprintf(helpOut, "# %16s\t%-32s\t%s\n", "-h", "Shows this help message", "");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int cmp_main(int argc, char *argv[]) {
	
	int args, t_len, s_len;
	char *templatefilename, *secondfilename;
	FILE *templatefile1, *templatefile2;
	HashMapKMA *t1, *t2;
	
	/* init */
	t_len = 0;
	s_len = 0;
	templatefilename = 0;
	secondfilename = 0;
	
	/* PARSE COMMAND LINE OPTIONS */
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-t_db") == 0) {
			++args;
			if(args < argc) {
				t_len = strlen(argv[args]);
				templatefilename = smalloc(t_len + 64);
				strcpy(templatefilename, argv[args]);
			}
		} else if(strcmp(argv[args], "-s_db") == 0) {
			++args;
			if(args < argc) {
				s_len = strlen(argv[args]);
				secondfilename = smalloc(s_len + 64);
				strcpy(secondfilename, argv[args]);
			}
		} else if(strcmp(argv[args], "-tmp") == 0) {
			if(++args < argc) {
				if(argv[args][0] != '-') {
					if(argv[args][strlen(argv[args]) - 1] != '/') {
						fprintf(stderr, "Invalid output directory specified.\n");
						exit(1);
					}
					tmpF(argv[args]);
				} else {
					--args;
				}
			}
		} else if(strcmp(argv[args], "-v") == 0) {
			fprintf(stdout, "KMA_index-%s\n", KMA_VERSION);
			exit(0);
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else {
			fprintf(stderr, "# Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, "# Printing help message:\n");
			helpMessage(1);
		}
		++args; 
	}
	
	/* check if all input was parsed */
	if(!t_len || !s_len) {
		fprintf(stderr, "Insufficient number of agruments parsed.\n");
		helpMessage(1);
	}
	if(strcmp(templatefilename, secondfilename) == 0) {
		fprintf(stderr, "Indexes to compare cannot be the same.\n");
		return 1;
	}
	
	/* merge *.comp.b */
	fprintf(stderr, "# Comparing *.comp.b\n");
	strcpy(templatefilename + t_len, ".comp.b");
	strcpy(secondfilename + s_len, ".comp.b");
	templatefile1 = sfopen(templatefilename, "rb");
	templatefile2 = sfopen(secondfilename, "rb");
	t1 = smalloc(sizeof(HashMapKMA));
	t2 = smalloc(sizeof(HashMapKMA));
	//hashMapKMA_load(t1, templatefile1, templatefilename);
	//hashMapKMA_load(t2, templatefile2, secondfilename);
	hashMapKMAmmap(t1, templatefile1);
	hashMapKMAmmap(t2, templatefile2);
	if(hashMapKMA_cmp(t1, t2)) {
		fprintf(stderr, "# Hashmaps does not match.\n");
	} else {
		fprintf(stderr, "# Hashmaps match.\n");
	}
	hashMapKMA_munmap(t1);
	hashMapKMA_munmap(t2);
	fclose(templatefile1);
	fclose(templatefile2);
	
	return 0;
}
