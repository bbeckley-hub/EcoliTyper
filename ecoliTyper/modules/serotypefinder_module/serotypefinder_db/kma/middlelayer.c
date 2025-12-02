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

#include <limits.h>
#include <stdlib.h>
#include "hashmapkma.h"
#include "middlelayer.h"
#include "pherror.h"

MiddleLayer * MiddleLayer_init(long unsigned size) {
	
	MiddleLayer *dest;
	
	dest = malloc(sizeof(MiddleLayer));
	if(dest) {
		dest->n = 0;
		dest->size = size;
		dest->layer = calloc(size << 1, sizeof(long unsigned));
		if(!dest->layer) {
			free(dest);
			dest = 0;
		}
	}
	
	return dest;
}

MiddleLayer * AlternativeLayer_init(long unsigned size) {
	
	MiddleLayer *dest;
	
	dest = malloc(sizeof(MiddleLayer));
	if(dest) {
		dest->n = 0;
		dest->size = size;
		dest->layer = calloc(size + (size << 1), sizeof(long unsigned));
		if(!dest->layer) {
			free(dest);
			dest = 0;
		}
	}
	
	return dest;
}

void MiddleLayer_free(MiddleLayer *src) {
	
	if(src) {
		if(src->layer) {
			free(src->layer);
		}
		free(src);
	}
}

MiddleLayer * MiddleLayer_resize(MiddleLayer *dest) {
	
	long unsigned i, *layer;
	
	/* reallocate layer */
	layer = realloc(dest->layer, (dest->size << 2) * sizeof(long unsigned));
	if(!layer) {
		ERROR();
	}
	
	/* set new part to zero */
	dest->layer = layer;
	i = (dest->size <<= 1);
	layer += i++ - 1;
	while(--i) {
		*++layer = 0;
	}
	
	return dest;
}

MiddleLayer * AlternativeLayer_resize(MiddleLayer *dest) {
	
	long unsigned i, *layer;
	
	/* reallocate layer */
	layer = realloc(dest->layer, ((dest->size << 1) + (dest->size << 2)) * sizeof(long unsigned));
	if(!layer) {
		ERROR();
	}
	
	/* set new part to zero */
	dest->layer = layer;
	i = dest->size + (dest->size << 1);
	dest->size <<= 1;
	layer += i++ - 1;
	while(--i) {
		*++layer = 0;
	}
	
	return dest;
}

MiddleLayer * MiddleLayer_readjust(MiddleLayer *dest) {
	
	long unsigned *layer;
	
	if(dest->n == dest->size) {
		return dest;
	}
	
	/* reallocate layer */
	layer = realloc(dest->layer, (dest->n << 1) * sizeof(long unsigned));
	if(layer) {
		dest->size = dest->n;
		dest->layer = layer;
	}
	
	return dest;
}

MiddleLayer * AlternativeLayer_readjust(MiddleLayer *dest) {
	
	long unsigned *layer;
	
	if(dest->n == dest->size) {
		return dest;
	}
	
	/* reallocate layer */
	layer = realloc(dest->layer, (dest->n + (dest->n << 1)) * sizeof(long unsigned));
	if(layer) {
		dest->size = dest->n;
		dest->layer = layer;
	}
	
	return dest;
}

long unsigned MiddleLayer_populate(MiddleLayer *dest, HashMapKMA *src, const long unsigned offset) {
	
	unsigned *values;
	short unsigned *values_s;
	long unsigned n, size, index, vindex, *layer;
	
	n = dest->n;
	size = dest->size;
	layer = dest->layer + (n << 1);
	if(src->DB_size < USHRT_MAX) {
		values = 0;
		values_s = src->values_s;
	} else {
		values = src->values;
		values_s = 0;
	}
	index = offset;
	vindex = index + src->v_index;
	while(index < vindex) {
		/* resize */
		if(n == size) {
			MiddleLayer_resize(dest);
			size = dest->size;
			layer = dest->layer + (n << 1);
		}
		
		/* update layer */
		*layer = index;
		layer += 2;
		++n;
		
		/* go to next signature */
		if(values) {
			index += (*values + 1);
			values += (*values + 1);
		} else {
			index += (*values_s + 1);
			values_s += (*values_s + 1);
		}
	}
	
	/* set new n */
	dest->n = n;
	
	return index;
}

long unsigned MiddleLayer_search(MiddleLayer *src, long unsigned v_index) {
	
	long unsigned mid, l1, l2, l_index, *layer;
	
	layer = src->layer;
	l1 = 0;
	l2 = src->n - 1;
	while(l1 <= l2) {
		mid = (l1 + l2) >> 1;
		l_index = layer[mid << 1]; /* remember layer contains two elements per index */
		if(l_index < v_index) {
			l1 = mid + 1;
		} else if(v_index < l_index) {
			l2 = mid - 1;
		} else {
			return mid;
		}
	}
	
	/* should not be possible */
	return src->n;
}

long unsigned AlternativeLayer_add(MiddleLayer *src, long unsigned v_index, long unsigned v_index1, long unsigned v_index2) {
	
	long unsigned index, *layer;
	
	/* check if combination is unique */
	if((index = v_index) != src->n) {
		do {
			layer = src->layer + (index + (index << 1)); /* remember layer contains three elements per index */
			if(layer[0] == v_index1 && layer[1] == v_index2) {
				return index;
			} else {
				index = layer[2];
			}
		} while(index);
		
		/* add index to new combination */
		layer[2] = src->n;
	}
	
	/* resize */
	if((src->n + 1) == src->size) {
		AlternativeLayer_resize(src);
	}
	
	/* add new combination */
	index = src->n++;
	layer = src->layer + (index + (index << 1));
	layer[0] = v_index1;
	layer[1] = v_index2;
	layer[2] = 0;
	
	/* if index == n - 1 -> increase final v_index */
	return index;
}
