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

#include "penalties.h"

#ifndef CHAIN
typedef struct alnPoints AlnPoints;
struct alnPoints {
	int size;
	int len;
	int *tStart;
	int *tEnd;
	int *qStart;
	int *qEnd;
	int *weight;
	int *score;
	int *next;
	Penalties *rewards;
};
#define CHAIN 1
#endif

/* pointer to chaining method */
extern int (*chainSeedsPtr)(AlnPoints *, int, int, int, unsigned *);
extern void (*trimSeedsPtr)(AlnPoints *points, int start);

/* FUNCTIONS */
AlnPoints * seedPoint_init(int size, Penalties *rewards);
void seedPoint_realloc(AlnPoints *dest, int size);
void seedPoint_free(AlnPoints *src);
int chainSeeds(AlnPoints *points, int q_len, int t_len, int kmersize, unsigned *mapQ);
int chainSeeds_circular(AlnPoints *points, int q_len, int t_len, int kmersize, unsigned *mapQ);
void trimSeeds(AlnPoints *points, int start);
void trimSeedsNoLead(AlnPoints *points, int start);
