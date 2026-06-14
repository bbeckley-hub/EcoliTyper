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
#include "compdna.h"
#include "hashmapkma.h"

extern int (*deConNode_ptr)(CompDNA *, HashMapKMA *, unsigned **);
extern int (*addCont)(HashMapKMA *, long unsigned, int, unsigned **);
int hashMap_addCont(HashMapKMA *dest, long unsigned key, int value, unsigned **Values);
int megaMap_addCont(HashMapKMA *dest, long unsigned index, int value, unsigned **Values);
int deConNode(CompDNA *qseq, HashMapKMA *finalDB, unsigned **Values);
int deConNode_sparse(CompDNA *qseq, HashMapKMA *finalDB, unsigned **Values);
unsigned deConDB(HashMapKMA *finalDB, char **inputfiles, int fileCount, char *trans, unsigned **Values);
