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
#include "hashmapkma.h"
#include "middlelayer.h"

#define hashMapKMA_compatible(t1, t2)(t1->kmersize == t2->kmersize && t1->prefix_len == t2->prefix_len && t1->prefix == t2->prefix && t1->mlen == t2->mlen)

long unsigned bucket_insertionsort(unsigned *key_index, unsigned *value_index, long unsigned *value_index_l, long unsigned mask_new, long unsigned mask_org, unsigned flag);
long unsigned bucket_insertionsort_l(long unsigned *key_index, unsigned *value_index, long unsigned *value_index_l, long unsigned mask_new, long unsigned mask_org, unsigned flag);
void hashMapKMA_sortbuckets(HashMapKMA *dest, HashMapKMA *src);
FILE * hashMapKMA_dumpbuckets(HashMapKMA *src);
long unsigned getV_index(long unsigned *exist, long unsigned size, long unsigned null_index, long unsigned v_indexes, MiddleLayer *middle, MiddleLayer *alternative);
unsigned * adjustV_index(long unsigned *exist_l, long unsigned size, long unsigned v_index);
long unsigned add_pairs(long unsigned *exist, long unsigned size, long unsigned null_index, long unsigned v_indexes, MiddleLayer *middle, MiddleLayer *alternative, HashMapKMA *t2);
void hashMapKMA_merge(HashMapKMA *dest, MiddleLayer *middle, MiddleLayer *alternative, FILE *tmp_1, FILE *tmp_2, long unsigned v_indexes);
unsigned loadValues1(unsigned *values, short unsigned *values_s, unsigned *values1, short unsigned *values1_s, long unsigned index);
unsigned loadValues2(unsigned *values, short unsigned *values_s, unsigned *values2, short unsigned *values2_s, long unsigned index, unsigned offset);
unsigned loadValues12(unsigned *values, short unsigned *values_s, unsigned *values1, short unsigned *values1_s, unsigned *values2, short unsigned *values2_s, long unsigned index1, long unsigned index2, unsigned offset);
void hashMapKMA_dumpmerge(HashMapKMA *src, HashMapKMA *t1, HashMapKMA *t2, MiddleLayer *middle, MiddleLayer *alternative, FILE *out);
HashMapKMA * merge_kmersignatures(HashMapKMA *t1, HashMapKMA *t2, MiddleLayer *middle, MiddleLayer *alternative);
int merge(char *templatefilename, char *templatefilename1, char *templatefilename2);
int merge_lengths(char *outname, char *inname1, char *inname2);
int cat(char *outname, char *inname1, char *inname2);
int merge_main(int argc, char *argv[]);
