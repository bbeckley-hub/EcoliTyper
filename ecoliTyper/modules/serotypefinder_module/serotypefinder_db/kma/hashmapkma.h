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
#ifndef HASHMAPKMA
#ifdef _WIN32
typedef int key_t;
#endif

typedef struct hashMapKMA HashMapKMA;
struct hashMapKMA {
	long unsigned size;				// size of DB
	long unsigned n;				// k-mers stored
	long unsigned mask;
	long unsigned null_index;		// null value
	long unsigned v_index;			// size of values
	unsigned kmersize;				// k
	unsigned prefix_len;			// prefix length
	long unsigned prefix;			// prefix
	int DB_size;
	unsigned shmFlag;
	unsigned mlen;
	unsigned flag;
	unsigned *exist;				// size long
	long unsigned *exist_l;			// size long, big DBs
	unsigned *values;				// compressed values
	short unsigned *values_s;		// compressed values, few templates
	unsigned *key_index;			// Relative
	long unsigned *key_index_l;		// Relative, 16 < k
	unsigned *value_index;			// Relative
	long unsigned *value_index_l;	// Relative, big DBs
};
#define HASHMAPKMA 1
#endif

/* DB size dependent pointers */
extern void (*hashMapKMA_destroy)(HashMapKMA *);
extern long unsigned (*getExistPtr)(const unsigned *, const long unsigned);
extern long unsigned (*getKeyPtr)(const unsigned *, const long unsigned);
extern long unsigned (*getValueIndexPtr)(const unsigned *, const long unsigned);
extern unsigned * (*getValuePtr)(const HashMapKMA *, const long unsigned);
extern unsigned * (*hashMap_get)(const HashMapKMA *, const long unsigned);
extern int (*intpos_bin_contaminationPtr)(const unsigned *, const int);
extern int (*getSizePtr)(const unsigned *);
extern void (*hashMapKMA_addKey_ptr)(HashMapKMA *, long unsigned, long unsigned);
extern void (*hashMapKMA_addValue_ptr)(HashMapKMA *, long unsigned, long unsigned);
extern void (*hashMapKMA_addExist_ptr)(HashMapKMA *, long unsigned, long unsigned);

/* HASHMAP FUNCTIONS */
long unsigned getExist(const unsigned *exist, const long unsigned pos);
long unsigned getExistL(const unsigned *exist, const long unsigned pos);
long unsigned getKey(const unsigned *key_index, const long unsigned pos);
long unsigned getKeyL(const unsigned *key_index, const long unsigned pos);
long unsigned getValueIndex(const unsigned *value_index, const long unsigned pos);
long unsigned getValueIndexL(const unsigned *value_index, const long unsigned pos);
unsigned * getValue(const HashMapKMA *dest, const long unsigned pos);
unsigned * getValueS(const HashMapKMA *dest, const long unsigned pos);
int getSize(const unsigned *values);
int getSizeS(const unsigned *values);
int intpos_bin_contamination(const unsigned *str1, const int str2);
int intpos_bin_contamination_s(const unsigned *Str1, const int str2);
unsigned * hashMap_getGlobal(const HashMapKMA *templates, const long unsigned key);
void loadPrefix(HashMapKMA *dest, FILE *file);
unsigned * megaMap_getGlobal(const HashMapKMA *templates, const long unsigned key);
int hashMapKMA_load(HashMapKMA *dest, FILE *file, const char *filename);
void hashMapKMA_load_shm(HashMapKMA *dest, FILE *file, const char *filename);
int hashMapKMAload(HashMapKMA *dest, FILE *file);
void hashMapKMA_dump(HashMapKMA *dest, FILE *out);
void megaMapKMA_dump(HashMapKMA *dest, FILE *out);
void hashMapKMA_addKey(HashMapKMA *dest, long unsigned index, long unsigned key);
void hashMapKMA_addKeyL(HashMapKMA *dest, long unsigned index, long unsigned key);
void hashMapKMA_addValue(HashMapKMA *dest, long unsigned index, long unsigned v_index);
void hashMapKMA_addValueL(HashMapKMA *dest, long unsigned index, long unsigned v_index);
void hashMapKMA_addExist(HashMapKMA *dest, long unsigned index, long unsigned relative);
void hashMapKMA_addExistL(HashMapKMA *dest, long unsigned index, long unsigned relative);
void hashMapKMA_free(HashMapKMA *dest);
