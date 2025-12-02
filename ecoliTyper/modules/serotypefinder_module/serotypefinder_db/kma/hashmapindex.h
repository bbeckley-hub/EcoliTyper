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

#ifndef HASHMAPINDEX
typedef struct hashMap_index HashMap_index;
struct hashMap_index {
	unsigned len; // seqlen
	unsigned size; // size of index
	int *index; // k-mer posititions in seq
	long unsigned *seq; // 2-bit sequence
	unsigned kmerindex;
};
#define HASHMAPINDEX 1
#endif

/* pointers determining how indexes a stored */
extern void (*destroyPtr)(HashMap_index *);
extern HashMap_index * (*alignLoadPtr)(HashMap_index *, int, int, int, int, long unsigned, long unsigned);

int hashMap_index_initialize(HashMap_index *dest, int len, int kmerindex);
void hashMap_index_set(HashMap_index *dest);
void hashMap_index_destroy(HashMap_index *dest);
int hashMap_index_get(const HashMap_index *dest, long unsigned key, unsigned shifter);
int hashMap_index_get_bound(const HashMap_index *dest, long unsigned key, int min, int max, unsigned shifter);
int hashMap_index_getDubPos(const HashMap_index *dest, long unsigned key, int value, unsigned shifter);
int hashMap_index_getNextDubPos(const HashMap_index *dest, long unsigned key, int min, int max, unsigned index, unsigned shifter);
void hashMap_index_add(HashMap_index *dest, long unsigned key, int newpos, unsigned shifter);
void hashMapIndex_add(HashMap_index *dest, long unsigned key, int newpos);
HashMap_index * hashMap_index_load(HashMap_index *src, int seq, int index, int len, int kmersize);
void hashMap_index_dump(HashMap_index *src, FILE *seq, FILE *index);
HashMap_index * hashMap_index_build(HashMap_index *src, int seq, int len, int kmersize);
HashMap_index * alignLoad_fly(HashMap_index *dest, int seq_in, int index_in, int len, int kmersize, long unsigned seq_index, long unsigned index_index);
HashMap_index * alignLoad_fly_mem(HashMap_index *dest, int seq_in, int index_in, int len, int kmersize, long unsigned seq_index, long unsigned index_index);
HashMap_index * alignLoad_fly_build(HashMap_index *dest, int seq_in, int index_in, int len, int kmersize, long unsigned seq_index, long unsigned index_index);
HashMap_index * alignLoad_fly_build_mem(HashMap_index *dest, int seq_in, int index_in, int len, int kmersize, long unsigned seq_index, long unsigned index_index);
HashMap_index * alignLoad_skip(HashMap_index *dest, int seq_in, int index_in, int len, int kmersize, long unsigned seq_index, long unsigned index_index);
HashMap_index * alignLoad_fly_shm(HashMap_index *dest, int seq_in, int index_in, int len, int kmersize, long unsigned seq_index, long unsigned index_index);
HashMap_index * alignLoad_shm_initial(char *templatefilename, int file_len, int seq_in, int index_in, int kmersize);
void alignClean(HashMap_index *template_index);
void alignClean_shm(HashMap_index *template_index);
