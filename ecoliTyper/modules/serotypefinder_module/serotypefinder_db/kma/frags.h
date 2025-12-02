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
#include <stdio.h>
#include "filebuff.h"
#include "qseqs.h"

#ifndef FRAG
typedef struct frag Frag;
struct frag {
	int buffer[7];
	unsigned char *qseq;
	unsigned char *header;
	struct frag *next;
};
#define FRAG 1
#endif

FILE * printFrags(Frag **alignFrags, int DB_size);
void updateAllFrag(unsigned char *qseq, int q_len, int bestHits, int best_read_score, int *best_start_pos, int *best_end_pos, int *bestTemplates, Qseqs *header, FileBuff *dest);
