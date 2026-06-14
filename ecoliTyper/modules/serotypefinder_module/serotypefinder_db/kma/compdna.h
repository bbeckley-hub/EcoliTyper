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
#include "qseqs.h"
#ifndef COMPDNA
typedef struct compDNA CompDNA;
struct compDNA {
	unsigned seqlen;
	unsigned size;
	unsigned complen;
	long unsigned *seq;
	int *N;
};
#define COMPDNA 1
#endif

void allocComp(CompDNA *compressor, unsigned size);
void reallocComp(CompDNA *compressor, unsigned size);
void freeComp(CompDNA *compressor);
CompDNA * setComp(unsigned size);
void destroyComp(CompDNA *src);
void resetComp(CompDNA *compressor);
void compDNA(CompDNA *compressor, unsigned char *seq, int seqlen);
int compDNAref(CompDNA *compressor, unsigned char *qseq, int seqlen);
void unCompDNA(CompDNA *compressor, unsigned char *seq);
void qseqCompDNA(CompDNA *compressor, Qseqs *qseq);
long unsigned binRev(long unsigned mer);
void rc_comp(CompDNA *compressor, CompDNA *compressor_rc);
void comp_rc(CompDNA *compressor);
void dumpComp(CompDNA *compressor, FILE* file);
int loadComp(CompDNA *compressor, FILE* file);
int getComp(CompDNA *compressor, FILE* file);
void dallocComp(CompDNA *compressor, unsigned size);
