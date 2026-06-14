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


#ifndef VALUESHASH
typedef struct valuesTable ValuesTable;
typedef struct valuesHash ValuesHash;

struct valuesTable {
	long unsigned v_index;
	unsigned *values;
	struct valuesTable *next;
};

struct valuesHash {
	long unsigned n;
	long unsigned size;
	struct valuesTable **table;
	int DB_size;
};
#define VALUESHASH 1
#endif

extern long unsigned (*valuesKeyPtr)(unsigned *, int);
extern int (*cmpValuesPtr)(unsigned *, unsigned *, unsigned);
ValuesHash * initialize_hashValues(long unsigned size, int DB_size);
void valuesHash_destroy(ValuesHash *src);
long unsigned valuesKey(unsigned *values, int DB_size);
long unsigned huValuesKey(unsigned *valuesOrg, int DB_size);
unsigned uSize(unsigned *values);
unsigned huSize(unsigned *valuesOrg);
int cmpValues(unsigned *s1, unsigned *s2, unsigned len);
int cmpHuValues(unsigned *s1_org, unsigned *s2_org, unsigned len_org);
long unsigned valuesHash_add(ValuesHash *src, unsigned *newValues, long unsigned v_index);
