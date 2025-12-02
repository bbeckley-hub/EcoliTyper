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
#include "hashtable.h"

#ifndef STRUCTHASHMAP
typedef struct hashMap HashMap;
struct hashMap {
	/* open hash structure */
	long unsigned size;		// size of DB
	long unsigned n;		// k-mers stored
	long unsigned mask;
	unsigned kmersize;		// k
	unsigned prefix_len;	// prefix length
	long unsigned prefix;	// prefix
	unsigned mlen;			// Minimizer length
	unsigned flag;			// flag to describe k-mer format: std, hom, min, minhom
	HashTable **table;		// org
	unsigned **values;		// ME
	int DB_size;
};
#define STRUCTHASHMAP 1
#endif

extern int (*hashMap_add)(HashMap *, long unsigned, unsigned);
extern unsigned * (*hashMapGet)(HashMap *, long unsigned);
extern void (*addUniqueValues)(HashMap *, long unsigned, unsigned *);
extern unsigned * (*updateValuePtr)(unsigned *, unsigned);
HashMap * hashMap_initialize(const long unsigned size, const unsigned kmersize, const unsigned mlen, const unsigned flag);
int megaMap_addKMA(HashMap *templates, long unsigned key, unsigned value);
unsigned * megaMap_getValue(HashMap *templates, long unsigned key);
void hashMap2megaMap(HashMap *templates, HashTable *table);
unsigned * updateValue(unsigned *values, unsigned value);
unsigned * updateShortValue(unsigned *valuesOrg, unsigned value);
int hashMap_addKMA(HashMap *templates, long unsigned key, unsigned value);
unsigned * hashMapGetValue(HashMap *templates, long unsigned key);
void hashMap_addUniqueValues(HashMap *dest, long unsigned key, unsigned *values);
void megaMap_addUniqueValues(HashMap *dest, long unsigned key, unsigned *values);
unsigned * HU2U(unsigned *values);
void convertToU(HashMap *templates);
