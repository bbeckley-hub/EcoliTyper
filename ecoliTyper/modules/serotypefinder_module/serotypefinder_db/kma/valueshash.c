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

#include <stdlib.h>
#include "pherror.h"
#include "valueshash.h"

long unsigned (*valuesKeyPtr)(unsigned *, int);
int (*cmpValuesPtr)(unsigned *, unsigned *, unsigned);

ValuesHash * initialize_hashValues(long unsigned size, int DB_size) {
	
	ValuesHash *dest;
	
	dest = smalloc(sizeof(ValuesHash));
	dest->n = 0;
	dest->size = size;
	dest->DB_size = DB_size;
	
	dest->table = calloc(size, sizeof(ValuesTable *));
	if(!dest->table) {
		ERROR();
	}
	
	return dest;
}

void valuesHash_destroy(ValuesHash *src) {
	
	long unsigned i;
	ValuesTable *node, *next;
	
	i = src->size;
	while(i--) {
		for(node = src->table[i]; node != 0; node = next) {
			next = node->next;
			free(node);
		}
	}
	free(src->table);
	free(src);
}

long unsigned valuesKey(unsigned *values, int DB_size) {
	
	unsigned i;
	long unsigned key;
	
	key = 0;
	for(i = 0; i <= *values; ++i) {
		key = key * DB_size + values[i];
	}
	
	return key;
}

long unsigned huValuesKey(unsigned *valuesOrg, int DB_size) {
	
	unsigned i;
	long unsigned key;
	short unsigned *values;
	
	values = (short unsigned *)(valuesOrg);
	key = 0;
	for(i = 0; i <= *values; ++i) {
		key = key * DB_size + values[i];
	}
	
	return key;
}

unsigned uSize(unsigned *values) {
	return *values + 1;
}

unsigned huSize(unsigned *valuesOrg) {
	short unsigned *values;
	values = (short unsigned *)(valuesOrg);
	return *values + 1;
}

int cmpValues(unsigned *s1, unsigned *s2, unsigned len) {
	
	if(len == 0) {
		return 1;
	}
	
	while(len--) {
		if(s1[len] != s2[len]) {
			return 0;
		}
	}
	
	return 1;
}

int cmpHuValues(unsigned *s1_org, unsigned *s2_org, unsigned len_org) {
	
	short unsigned *s1, *s2, len;
	
	len = len_org;
	--len;
	
	if(len == 0) {
		return 1;
	}
	s1 = (short unsigned *)(s1_org);
	s2 = (short unsigned *)(s2_org);
	
	while(len--) {
		if(s1[len] != s2[len]) {
			return 0;
		}
	}
	
	return 1;
}

long unsigned valuesHash_add(ValuesHash *src, unsigned *newValues, long unsigned v_index) {
	
	/* return 0 if values are new, 
	else return first index of seen value */
	
	unsigned *values;
	long unsigned index;
	ValuesTable *node;
	
	/* get index */
	index = valuesKeyPtr(newValues, src->DB_size) % src->size;
	
	/* search for values */
	for(node = src->table[index]; node != 0; node = node->next) {
		values = node->values;
		if(*values == *newValues && cmpValuesPtr(values + 1, newValues + 1, newValues[0])) { // Value exists
			return node->v_index;
		}
	}
	
	/* new values */
	++src->n;
	node = smalloc(sizeof(ValuesTable));
	node->v_index = v_index;
	node->values = newValues;
	node->next = src->table[index];
	src->table[index] = node;
	
	return v_index;
}
