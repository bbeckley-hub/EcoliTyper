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

#include "hashmapkma.h"
#ifndef MIDDLELAYER
typedef struct middlelayer MiddleLayer;
typedef struct Alternativelayer Alternativelayer;
struct middlelayer {
	long unsigned n;
	long unsigned size;
	long unsigned *layer;	/* vindex, alternative */
	//long unsigned *val_indexes;	/* vindex, alternative */
	//long unsigned *alt_indexes;	/* vindex1, vindex2, alternative */
};
struct alternativelayer {
	long unsigned n;
	long unsigned size;
	long unsigned *alt_indexes;	/* vindex1, vindex2, alternative */
};
#define MIDDLELAYER 1
#endif

MiddleLayer * MiddleLayer_init(long unsigned size);
MiddleLayer * AlternativeLayer_init(long unsigned size);
void MiddleLayer_free(MiddleLayer *src);
MiddleLayer * MiddleLayer_resize(MiddleLayer *dest);
MiddleLayer * AlternativeLayer_resize(MiddleLayer *dest);
MiddleLayer * MiddleLayer_readjust(MiddleLayer *dest);
MiddleLayer * AlternativeLayer_readjust(MiddleLayer *dest);
long unsigned MiddleLayer_populate(MiddleLayer *dest, HashMapKMA *src, const long unsigned offset);
long unsigned MiddleLayer_search(MiddleLayer *src, long unsigned v_index);
long unsigned AlternativeLayer_add(MiddleLayer *src, long unsigned v_index, long unsigned v_index1, long unsigned v_index2);
