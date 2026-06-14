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
#include <stdlib.h>
#include "hashmap.h"
#include "hashtable.h"
#include "pherror.h"
#include "stdstat.h"

int (*hashMap_add)(HashMap *, long unsigned, unsigned);
unsigned * (*hashMapGet)(HashMap *, long unsigned);
void (*addUniqueValues)(HashMap *, long unsigned, unsigned *);
unsigned * (*updateValuePtr)(unsigned *, unsigned);

HashMap * hashMap_initialize(const long unsigned size, const unsigned kmersize, const unsigned mlen, const unsigned flag) {
	
	HashMap *src;
	
	src = smalloc(sizeof(HashMap));
	src->size = size;
	src->n = 0;
	src->mask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	src->kmersize = kmersize;
	src->prefix_len = 0;
	src->prefix = 0;
	src->mlen = mlen;
	src->flag = flag;
	
	src->DB_size = 1;
	
	if((size - 1) == src->mask) {
		src->table = 0;
		src->values = calloc(src->size, sizeof(unsigned *));
		if(!src->values) {
			ERROR();
		}
	} else {
		src->values = 0;
		src->table = calloc(src->size, sizeof(HashTable *));
		if(!src->table) {
			ERROR();
		}
	}
	
	/* masking */
	--src->size;
	
	return src;
}

int megaMap_addKMA(HashMap *templates, long unsigned key, unsigned value) {
	
	unsigned *values;
	
	key &= templates->mask;
	values = templates->values[key];
	if(values == 0) {
		templates->values[key] = updateValuePtr(0, value);
		++templates->n;
		return 1;
	} else if((values = updateValuePtr(values, value))) {
		templates->values[key] = values;
		return 1;
	}
	
	return 0;
}

unsigned * megaMap_getValue(HashMap *templates, long unsigned key) {
	return templates->values[key & templates->mask];
}

void hashMap2megaMap(HashMap *templates, HashTable *table) {
	
	HashTable *node, *next;
	
	templates->table = 0;
	templates->size = templates->mask + 1;
	templates->values = calloc(templates->size, sizeof(unsigned *));
	if(!templates->values) {
		ERROR();
	}
	--templates->size;
	
	/* convert table */
	for(node = table; node != 0; node = next) {
		next = node->next;
		
		/* move values */
		templates->values[node->key & templates->mask] = node->value;
		
		/* clean */
		free(node);
	}
	
	/* clean */
	templates->table = 0;
	
	/* set pointers */
	hashMap_add = &megaMap_addKMA;
	hashMapGet = &megaMap_getValue;
}

unsigned * updateValue(unsigned *values, unsigned value) {
	
	if(!values) {
		values = smalloc(2 * sizeof(unsigned));
		values[0] = 1;
		values[1] = value;
	} else if(values[*values] == value) {
		return 0;
	} else {
		values[0]++;
		values = realloc(values, (values[0] + 1) * sizeof(unsigned));
		if(!values) {
			ERROR();
		}
		values[*values] = value;
	}
	
	return values;
}

unsigned * updateShortValue(unsigned *valuesOrg, unsigned value) {
	
	short unsigned *values;
	
	values = (short unsigned *)(valuesOrg);
	if(!values) {
		values = smalloc(2 * sizeof(short unsigned));
		values[0] = 1;
		values[1] = value;
	} else if(values[*values] == value) {
		return 0;
	} else {
		values[0]++;
		values = realloc(values, (values[0] + 1) * sizeof(short unsigned));
		if(!values) {
			ERROR();
		}
		values[*values] = value;
	}
	valuesOrg = (unsigned *)(values);
	
	return valuesOrg;
}

int hashMap_addKMA(HashMap *templates, long unsigned key, unsigned value) {
	
	unsigned *values, flag;
	long unsigned index, mask;
	HashTable *node, *next, *table;
	
	if(templates->flag) {
		murmur(index, key);
		index &= templates->size;
	} else {
		index = key & templates->size;
	}
	
	/* check if key exists */
	for(node = templates->table[index]; node != 0; node = node->next) {
		if(key == node->key) {
			if((values = updateValuePtr(node->value, value))) {
				node->value = values;
				return 1;
			} else {
				return 0;
			}
		}
	}
	
	/* new value check if there is space */
	if(templates->n == templates->size) {
		++templates->size;
		/* link table */
		table = 0;
		index = templates->size;
		while(index--) {
			for(node = templates->table[index]; node != 0; node = next) {
				next = node->next;
				node->next = table;
				table = node;
			}
		}
		free(templates->table);
		
		/* check for megamap */
		templates->size <<= 1;
		if((templates->mask + 1) <= (templates->size << 1)) {
			hashMap2megaMap(templates, table);
			return megaMap_addKMA(templates, key, value);
		}
		
		/* reallocate */
		templates->table = calloc(templates->size, sizeof(HashTable));
		if(!templates->table) {
			ERROR();
		}
		mask = --templates->size;
		
		flag = templates->flag;
		for(node = table; node != 0; node = next) {
			next = node->next;
			if(flag) {
				murmur(index, node->key);
				index &= mask;
			} else {
				index = node->key & mask;
			}
			node->next = templates->table[index];
			templates->table[index] = node;
		}
		
		
		if(flag) {
			murmur(index, key);
			index &= mask;
		} else {
			index = key & mask;
		}
	}
	
	/* add new value */
	node = smalloc(sizeof(HashTable));
	
	/* key */
	node->key = key;
	
	/* value */
	node->value = updateValuePtr(0, value);
	
	/* push it */
	node->next = templates->table[index];
	templates->table[index] = node;
	
	++templates->n;
	
	return 1;
}

unsigned * hashMapGetValue(HashMap *templates, long unsigned key) {
	
	long unsigned index;
	HashTable *node;
	
	if(templates->flag) {
		murmur(index, key);
		index &= templates->size;
	} else {
		index = key & templates->size;
	}
	
	for(node = templates->table[index]; node != 0; node = node->next) {
		if(key == node->key) {
			return node->value;
		}
	}
	
	return 0;
}

void hashMap_addUniqueValues(HashMap *dest, long unsigned key, unsigned *values) {
	
	long unsigned index;
	HashTable *node;
	
	node = smalloc(sizeof(HashTable));
	node->key = key;
	node->value = values;
	if(dest->flag) {
		murmur(index, key);
		index &= dest->size;
	} else {
		index = key & dest->size;
	}
	node->next = dest->table[index];
	dest->table[index] = node;
	dest->n++;
}

void megaMap_addUniqueValues(HashMap *dest, long unsigned key, unsigned *values) {
	
	dest->values[key & dest->mask] = values;
	dest->n++;
}

unsigned * HU2U(unsigned *values) {
	
	int i;
	short unsigned *hu_values;
	
	hu_values = (short unsigned *)(values);
	values = realloc(values, (hu_values[0] + 1) * sizeof(unsigned));
	if(!values) {
		ERROR();
	} else {
		hu_values = (short unsigned *)(values);
	}
	
	i = *hu_values + 1;
	while(i--) {
		values[i] = hu_values[i];
	}
	
	return values;
}

void convertToU(HashMap *templates) {
	
	long unsigned index;
	HashTable *node;
	
	/* convert values */
	index = templates->size + 1;
	if(templates->table) {
		while(index--) {
			for(node = templates->table[index]; node != 0; node = node->next) {
				node->value = HU2U(node->value);
			}
		}
	} else {
		while(index--) {
			if(templates->values[index]) {
				templates->values[index] = HU2U(templates->values[index]);
			}
		}	
	}
	
	/* set pointers */
	updateValuePtr = &updateValue;
}
