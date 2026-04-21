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
#include "hashmapkma.h"
#include "hashmapkmers.h"

#ifndef HASHTABLE
typedef struct hashTable HashTable;
typedef struct hit Hit;

struct hashTable {
	long unsigned key;
	unsigned *value;
	struct hashTable *next;
};

struct hit {
	long unsigned n;
	long unsigned tot;
};
#define HASHTABLE 1
#endif

int intpos_bin(const unsigned *str1, const int str2);
HashTable * collect_Kmers(const HashMapKMA *templates, unsigned *Scores, long unsigned *Scores_tot, HashMap_kmers *foundKmers, Hit *hits);
HashTable ** collect_Kmers_deCon(const HashMapKMA *templates, unsigned *Scores, long unsigned *Scores_tot, HashMap_kmers *foundKmers, Hit *hits, int contamination);
HashTable * withDraw_Kmers(unsigned *Scores, long unsigned *Scores_tot, HashTable *kmerList, int template, Hit *hits);
Hit withDraw_Contamination(unsigned *Scores, long unsigned *Scores_tot, HashTable *kmerList, HashTable *deConTable, int template, Hit hits);
