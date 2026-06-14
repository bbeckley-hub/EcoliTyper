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

#ifndef HASHMAPKMERS
typedef struct hashTable_kmers HashTable_kmers;
typedef struct hashMap_kmers HashMap_kmers;

struct hashTable_kmers {
	long unsigned key;
	int value;
	struct hashTable_kmers *next;
};

struct hashMap_kmers {
	long unsigned size;
	long unsigned n;
	unsigned flag;
	unsigned filler; //unused
	struct hashTable_kmers **table;
};
#define HASHMAPKMERS 1
#endif

void hashMap_kmers_initialize(HashMap_kmers *dest, long unsigned newSize);
void reallocHashMap_kmers(HashMap_kmers *dest);
void hashMap_kmers_CountIndex(HashMap_kmers *dest, long unsigned key);
int hashMap_CountKmer(HashMap_kmers *dest, long unsigned key);
void emptyHash(HashMap_kmers *dest);
