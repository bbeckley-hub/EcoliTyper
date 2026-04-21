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
#include <stdlib.h>
#include "kmerlink.h"
#include "pherror.h"

KmerLink * initKmerLink(unsigned size) {
	
	KmerLink *dest;
	
	dest = smalloc(sizeof(KmerLink));
	dest->size = size;
	dest->n = 0;
	dest->list = smalloc(size * sizeof(KmerAnker));
	
	return dest;
}

void reallocKmerLink(KmerLink *src, unsigned newSize) {
	
	src->size = newSize;
	src->list = realloc(src->list, newSize * sizeof(KmerAnker));
	if(!src->list) {
		ERROR();
	}
}

KmerAnker * pushKmerLink(KmerLink *src, KmerAnker *node) {
	
	KmerAnker *dest;
	
	if(src->size == ++src->n) {
		reallocKmerLink(src, src->size << 1);
	}
	dest = src->list + (src->n - 1);
	*dest = *node;
	
	return dest;
}

KmerAnker * popKmerLink(KmerLink *src, int n) {
	
	if(n < src->n) {
		src->n -= n;
	} else {
		src->n = 0;
	}
	
	return src->list + src->n;
}

void destroyKmerLink(KmerLink *src) {
	
	free(src->list);
	free(src);
}
