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

#ifndef HASHMAPCCI
typedef struct hashMapCCI HashMapCCI;
struct hashMapCCI {
	long unsigned mask;
	unsigned len; // seqlen
	long unsigned size; // size of index
	unsigned kmerindex;
	unsigned bandwidth;
	unsigned ndup;
	int *index; // k-mer posititions in seq / chain
	int *chain; // k-mer collision positions
	long unsigned *seq; // 2-bit sequence
	unsigned cci_size; // size of closed chains
	unsigned cci_next; // next available chain
	unsigned cci_avail; // avail chains
};
#define HASHMAPCCI 1
#endif

/* pointers determining how indexes a stored */
extern HashMapCCI * (*alignLoadPtr)(HashMapCCI *, int, int, int, long unsigned);

long unsigned hashMapCCI_initialize(HashMapCCI *dest, int len, int kmerindex);
void hashMapCCI_destroy(HashMapCCI *dest);
int hashMapCCI_get(const HashMapCCI *dest, long unsigned key, unsigned shifter);
int hashMapCCI_get_bound(const HashMapCCI *dest, long unsigned key, int min, int max, unsigned shifter);
int * hashMapCCI_getDubPos(const HashMapCCI *dest, long unsigned key, int value, unsigned shifter);
int * hashMapCCI_getNextDubPos(const HashMapCCI *dest, int *chain, long unsigned key, int min, int max, unsigned shifter);
int defragChain(HashMapCCI *dest, int size, int shifter);
int newChain(HashMapCCI *dest, int pos, int newpos, long unsigned kmer, int shifter);
int extendChain(HashMapCCI *dest, int chainpos, int newpos, long unsigned kmer, int shifter);
void hashMapCCI_add(HashMapCCI *dest, long unsigned key, int newpos, unsigned shifter);
void hashMapCCI_add_thread(HashMapCCI *dest, long unsigned key, int newpos, unsigned shifter);
HashMapCCI * hashMapCCI_load(HashMapCCI *src, int seq, int len, int kmersize);
HashMapCCI * hashMapCCI_load_thread(HashMapCCI *src, int seq, int len, int kmersize, int thread_num);
void hashMapCCI_dump(HashMapCCI *src, FILE *seq);
HashMapCCI * alignLoad_fly(HashMapCCI *dest, int seq_in, int len, int kmersize, long unsigned seq_index);
HashMapCCI * alignLoad_fly_mem(HashMapCCI *dest, int seq_in, int len, int kmersize, long unsigned seq_index);
HashMapCCI * alignLoad_skip(HashMapCCI *dest, int seq_in, int len, int kmersize, long unsigned seq_index);
