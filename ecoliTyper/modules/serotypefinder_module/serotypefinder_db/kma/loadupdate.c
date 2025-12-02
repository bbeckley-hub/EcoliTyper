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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "loadupdate.h"
#include "pherror.h"
#include "hashmap.h"
#include "hashmapkma.h"
#include "updateindex.h"

void * memdup(const void * src, size_t size) {
	
	void *dest;
	
	dest = smalloc(size);
	memcpy(dest, src, size);
	
	return dest;
}

int CP(char *templatefilename, char *outputfilename) {
	
	int bytes, buffSize;
	char *buffer;
	FILE *file_in, *file_out;
	
	if(strcmp(templatefilename, outputfilename) == 0) {
		return 1;
	}
	
	file_in = sfopen(templatefilename, "rb");
	file_out = sfopen(outputfilename, "wb");
	buffer = smalloc((buffSize = 1048576));
	
	while((bytes = fread(buffer, 1, buffSize, file_in))) {
		cfwrite(buffer, 1, bytes, file_out);
	}
	
	fclose(file_in);
	fclose(file_out);
	free(buffer);
	
	return 0;
}

HashMap * hashMapKMA_openChains(HashMapKMA *src) {
	
	long unsigned i, key;
	unsigned *values;
	HashMap *dest;
	
	if(src->mask != src->size) {
		free(src->exist);
	}
	dest = hashMap_initialize(src->size + 1, src->kmersize, src->mlen, src->flag);
	dest->prefix = src->prefix;
	dest->prefix_len = src->prefix_len;
	dest->DB_size = src->DB_size;
	
	if(dest->size == dest->mask) {
		addUniqueValues = &megaMap_addUniqueValues;
	} else {
		addUniqueValues = &hashMap_addUniqueValues;
	}
	
	if(src->mask != src->size) {
		/* norm */
		i = src->n;
		while(i--) {
			key = getKeyPtr(src->key_index, i);
			values = getValuePtr(src, getValueIndexPtr(src->value_index, i));
			values = memdup(values, getSizePtr(values));
			addUniqueValues(dest, key, values);
		}
		free(src->key_index);
		free(src->value_index);
	} else {	
		/* mega */
		i = src->size + 1;
		while(i--) {
			if(getExistPtr(src->exist, i) != 1) {
				values = getValuePtr(src, i);
				values = memdup(values, getSizePtr(values));
				
				addUniqueValues(dest, i, values);
			}
		}
		free(src->exist);
	}
	free(src->values);
	free(src);
	
	return dest;
}

unsigned ** hashMapKMA_openValues(HashMapKMA *src) {
	
	long unsigned i, index, pos;
	unsigned *values, **Values;
	
	Values = smalloc(src->n * sizeof(unsigned *));
	index = src->n;
	if(src->mask != (src->size - 1)) {
		/* norm */
		i = src->n;
		while(i--) {
			/* get key and values */
			values = getValuePtr(src, getValueIndexPtr(src->value_index, i));
			Values[i] = memdup(values, getSizePtr(values));
			hashMapKMA_addValue_ptr(src, i, i);
		}
	} else {
		/* mega */
		i = src->size;
		getSizePtr = &getSizeS;
		while(i--) {
			if((pos = getExistPtr(src->exist, i)) != 1) {
				/* get key and values */
				values = getValuePtr(src, pos);
				Values[--index] = memdup(values, getSizePtr(values));
				hashMapKMA_addExist_ptr(src, i, index);
			} else {
				hashMapKMA_addExist_ptr(src, i, src->n);
			}
		}
	}
	free(src->values);
	src->values = 0;
	
	return Values;
}

unsigned load_DBs(char *templatefilename, char *outputfilename, unsigned **template_lengths, unsigned **template_ulengths, unsigned **template_slengths, HashMapKMA *finalDB) {
	
	int file_len, out_len, DB_size, kmerindex;
	FILE *infile;
	
	file_len = strlen(templatefilename);
	out_len = strlen(outputfilename);
	
	/* load hash */
	strcat(templatefilename, ".comp.b");
	infile = sfopen(templatefilename, "rb");
	if(hashMapKMAload(finalDB, infile)) {
		fprintf(stderr, "Wrong format of DB\n");
		exit(1);
	} else {
		finalDB->size--;
	}
	templatefilename[file_len] = 0;
	fclose(infile);
	
	/* load lengths */
	strcat(templatefilename, ".length.b");
	infile = sfopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	sfread(&DB_size, sizeof(unsigned), 1, infile);
	if(finalDB->prefix) {
		*template_lengths = smalloc((DB_size << 1) * sizeof(unsigned));
		*template_slengths = smalloc((DB_size << 1) * sizeof(unsigned));
		*template_ulengths = smalloc((DB_size << 1) * sizeof(unsigned));
		sfread(*template_lengths, sizeof(unsigned), DB_size, infile);
		sfread(*template_slengths, sizeof(unsigned), DB_size, infile);
		sfread(*template_ulengths, sizeof(unsigned), DB_size, infile);
		kmerindex = **template_slengths;
		**template_ulengths = DB_size << 1;
		**template_slengths = DB_size << 1;
	} else {
		*template_lengths = smalloc((DB_size << 1) * sizeof(unsigned));
		*template_slengths = 0;
		*template_ulengths = 0;
		sfread(*template_lengths, sizeof(unsigned), DB_size, infile);
		kmerindex = **template_lengths;
		**template_lengths = DB_size << 1;
	}
	fclose(infile);
	
	/* cp name, seq and index */
	strcat(templatefilename, ".name");
	strcat(outputfilename, ".name");
	CP(templatefilename, outputfilename);
	templatefilename[file_len] = 0;
	outputfilename[out_len] = 0;
	
	strcat(templatefilename, ".seq.b");
	strcat(outputfilename, ".seq.b");
	CP(templatefilename, outputfilename);
	templatefilename[file_len] = 0;
	outputfilename[out_len] = 0;
	
	return kmerindex;
}
