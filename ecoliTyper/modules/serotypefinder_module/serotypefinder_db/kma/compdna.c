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

#include <stdlib.h>
#include <stdio.h>
#include "compdna.h"
#include "pherror.h"
#include "qseqs.h"
#include "stdnuc.h"

void allocComp(CompDNA *compressor, unsigned size) {
	
	compressor->seqlen = 0;
	if(size & 31) {
		compressor->size = (size >> 5) + 1;
		compressor->size <<= 5;
	} else {
		compressor->size = size;
	}
	
	compressor->seq = calloc(compressor->size >> 5, sizeof(long unsigned));
	compressor->N = malloc((compressor->size + 1) * sizeof(int));
	
	if(!compressor->seq || !compressor->N) {
		ERROR();
	}
	
	compressor->N[0] = 0;
	
}

void reallocComp(CompDNA *compressor, unsigned size) {
	
	if(size <= compressor->size) {
		return;
	} else if(size & 31) {
		compressor->size = (size >> 5) + 1;
		compressor->size <<= 5;
	} else {
		compressor->size = size;
	}
	
	compressor->seq = realloc(compressor->seq, (compressor->size >> 5) * sizeof(long unsigned));
	compressor->N = realloc(compressor->N, (compressor->size + 1) * sizeof(int));
	
	if(!compressor->seq || !compressor->N) {
		ERROR();
	}
}

void freeComp(CompDNA *compressor) {
	
	compressor->seqlen = 0;
	compressor->complen = 0;
	compressor->size = 0;
	
	free(compressor->seq);
	free(compressor->N);
	
}

CompDNA * setComp(unsigned size) {
	
	CompDNA *dest;
	
	dest = smalloc(sizeof(CompDNA));
	dest->seqlen = 0;
	dest->size = 
	dest->complen = 0;
	
	if(size & 31) {
		dest->size = (size >> 5) + 1;
		dest->size <<= 5;
	} else {
		dest->size = size;
	}
	
	dest->seq = calloc(dest->size >> 5, sizeof(long unsigned));
	if(!dest->seq) {
		ERROR();
	}
	dest->N = smalloc((dest->size + 1) * sizeof(int));
	*(dest->N) = 0;
	
	return dest;
}

void destroyComp(CompDNA *src) {
	
	free(src->seq);
	free(src->N);
	free(src);
}

void resetComp(CompDNA *compressor) {
	compressor->N[0] = 0;
	compressor->seqlen = 0;
	compressor->complen = 0;
}


void compDNA(CompDNA *compressor, unsigned char *seq, int seqlen) {
	
	int i, j, pos, end;
	
	/* get compressed length */
	compressor->seqlen = seqlen;
	if(seqlen & 31) {
		compressor->complen = (seqlen >> 5) + 1;
	} else {
		compressor->complen = seqlen >> 5;
	}
	
	/* reallocate */
	if(compressor->size < (compressor->complen << 5)) {
		compressor->size = (compressor->complen << 5);
		compressor->seq = realloc(compressor->seq, (compressor->size >> 5) * sizeof(long unsigned));
		compressor->N = realloc(compressor->N, (compressor->size + 1) * sizeof(int));
		if(!compressor->seq || !compressor->N) {
			ERROR();
		}
	}
	
	/* load seq */
	compressor->N[0] = 0;
	for(i = 0, pos = 0; i < seqlen; i += 32) {
		end = (i + 32 < seqlen) ? i + 32 : seqlen;
		pos = i >> 5;
		for(j = i; j < end; ++j) {
			if(seq[j] == 4) {
				compressor->seq[pos] <<= 2;
				compressor->N[0]++;
				compressor->N[compressor->N[0]] = j;
			} else {
				compressor->seq[pos] = (compressor->seq[pos] << 2) | seq[j];
			}
		}
	}
	if(seqlen & 31) {
		compressor->seq[pos] <<= (64 - ((seqlen & 31) << 1));
	}
}

int compDNAref(CompDNA *compressor, unsigned char *qseq, int seqlen) {
	
	int i, j, pos, end, bias;
	unsigned char *seq;
	
	/* trim leadin N's */
	seq = qseq;
	bias = 0;
	while(*seq == 4) {
		++seq;
		++bias;
	}
	seqlen -= bias;
	
	/* trim trailing N's */
	if(seqlen) {
		while(seq[--seqlen] == 4);
		++seqlen;
	}
	
	compressor->seqlen = seqlen;
	if(seqlen & 31) {
		compressor->complen = (seqlen >> 5) + 1;
	} else {
		compressor->complen = seqlen >> 5;
	}
	
	pos = 0;
	compressor->N[0] = 0;
	for(i = 0; i < seqlen; i += 32) {
		end = (i + 32 < seqlen) ? i + 32 : seqlen;
		pos = i >> 5;
		for(j = i; j < end; ++j) {
			if(seq[j] == 4) {
				compressor->seq[pos] <<= 2;
				compressor->N[0]++;
				compressor->N[compressor->N[0]] = j;
			} else {
				compressor->seq[pos] = (compressor->seq[pos] << 2) | seq[j];
			}
		}
	}
	if(seqlen & 31) {
		compressor->seq[pos] <<= (64 - ((seqlen & 31) << 1));
	}
	
	return bias;
}

void unCompDNA(CompDNA *compressor, unsigned char *seq) {
	
	int i;
	
	/* get nucs */
	i = compressor->seqlen;
	seq += i;
	*seq = 0;
	while(i--) {
		*--seq = getNuc(compressor->seq, i);
	}
	/*
	for(i = 0; i < compressor->seqlen; ++i) {
		seq[i] = getNuc(compressor->seq, i);
	}
	*/
	
	/* get N's */
	for(i = 1; i <= compressor->N[0]; ++i) {
		seq[compressor->N[i]] = 4;
	}
	
}

void qseqCompDNA(CompDNA *compressor, Qseqs *qseq) {
	
	qseq->len = compressor->seqlen;
	if(qseq->size <= qseq->len) {
		qseq->size = qseq->len << 1;
		free(qseq->seq);
		qseq->seq = smalloc(qseq->size);
	}
	unCompDNA(compressor, qseq->seq);
	
}

long unsigned binRev(long unsigned mer) {
	
	/* swap consecutive pairs */
	mer = ((mer >> 2) & 0x3333333333333333) | ((mer & 0x3333333333333333) << 2);
	/* swap nibbles */
	mer = ((mer >> 4) & 0x0F0F0F0F0F0F0F0F) | ((mer & 0x0F0F0F0F0F0F0F0F) << 4);
	/* swap bytes */
	mer = ((mer >> 8) & 0x00FF00FF00FF00FF) | ((mer & 0x00FF00FF00FF00FF) << 8);
	/* swap 2-bytes */
	mer = ((mer >> 16) & 0x0000FFFF0000FFFF) | ((mer & 0x0000FFFF0000FFFF) << 16);
	/* swap 4-bytes */
	return ((mer >> 32) | (mer << 32));
}

void rc_comp(CompDNA *compressor, CompDNA *compressor_rc) {
	
	int i, j, shift, r_shift;
	
	compressor_rc->seqlen = compressor->seqlen;
	compressor_rc->complen = compressor->complen;
	
	/* reallocate */
	if(compressor_rc->size < (compressor_rc->complen << 5)) {
		compressor_rc->size = (compressor_rc->complen << 5);
		compressor_rc->seq = realloc(compressor_rc->seq, (compressor_rc->size >> 5) * sizeof(long unsigned));
		compressor_rc->N = realloc(compressor_rc->N, (compressor_rc->size + 1) * sizeof(int));
		if(!compressor_rc->seq || !compressor_rc->N) {
			ERROR();
		}
	}
	
	/* reverse and complement*/
	for(i = 0, j = compressor->complen - 1; i < compressor->complen; ++i, --j) {
		compressor_rc->seq[j] = binRev(~compressor->seq[i]);
	}
	
	/* shift */
	if((compressor->seqlen & 31)) {
		shift = (((compressor->complen << 5) - compressor->seqlen) << 1);
		r_shift = 64 - shift;
		for(i = 0, j = 1; j < compressor->complen; i = j++) {
			compressor_rc->seq[i] = (compressor_rc->seq[i] << shift) | (compressor_rc->seq[j] >> r_shift);
		}
		compressor_rc->seq[i] <<= shift;
	}
	
	/* add N's */
	shift = compressor->seqlen - 1;
	compressor_rc->N[0] = compressor->N[0];
	for(i = 1, j = compressor->N[0]; j != 0; ++i, --j) {
		compressor_rc->N[i] = shift - compressor->N[j];
	}
}

void comp_rc(CompDNA *compressor) {
	
	int i, j, shift, r_shift;
	long unsigned carry;
	
	/* reverse and complement*/
	for(i = 0, j = compressor->complen - 1; i < j; ++i, --j) {
		carry = binRev(~compressor->seq[i]);
		compressor->seq[i] = binRev(~compressor->seq[j]);
		compressor->seq[j] = carry;
	}
	if(i == j) {
		compressor->seq[i] = binRev(~compressor->seq[i]);
	}
	
	/* shift */
	if((compressor->seqlen & 31)) {
		shift = (((compressor->complen << 5) - compressor->seqlen) << 1);
		r_shift = 64 - shift;
		for(i = 0, j = 1; j < compressor->complen; i = j++) {
			compressor->seq[i] = (compressor->seq[i] << shift) | (compressor->seq[j] >> r_shift);
		}
		compressor->seq[i] <<= shift;
	}
	
	/* add N's */
	i = 1;
	j = compressor->N[0];
	shift = compressor->seqlen - 1;
	for(i = 1, j = compressor->N[0]; i < j; ++i, --j) {
		r_shift = shift - compressor->N[i];
		compressor->N[i] = shift - compressor->N[j];
		compressor->N[j] = r_shift;
	}
	if(i == j) {
		compressor->N[i] = shift - compressor->N[j];
	}
	
}

void dumpComp(CompDNA *compressor, FILE* file) {
	
	sfwrite(&compressor->seqlen, sizeof(int), 1, file);
	sfwrite(&compressor->complen, sizeof(int), 1, file);
	sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, file);
	sfwrite(compressor->N, sizeof(int), compressor->N[0] + 1, file);
	
}

int loadComp(CompDNA *compressor, FILE* file) {
	
	if(fread(&compressor->seqlen, sizeof(int), 1, file)) {
		sfread(&compressor->complen, sizeof(int), 1, file);
		sfread(compressor->seq, sizeof(long unsigned), compressor->complen, file);
		sfread(compressor->N, sizeof(int), 1, file);
		sfread(compressor->N + 1, sizeof(int), compressor->N[0], file);
		return 1;
	}
	return 0;
}

int getComp(CompDNA *compressor, FILE* file) {
	
	compressor->seqlen = 0;
	if(!fread(&compressor->seqlen, sizeof(int), 1, file)) {
		return 0;
	}
	
	sfread(&compressor->complen, sizeof(int), 1, file);
	
	/* realloc */
	if(compressor->seqlen >= compressor->size) {
		free(compressor->N);
		free(compressor->seq);
		if(compressor->seqlen & 31) {
			compressor->size = (compressor->seqlen >> 5) + 1;
			compressor->size <<= 6;
		} else {
			compressor->size = compressor->seqlen << 1;
		}
		
		compressor->seq = calloc(compressor->size >> 5, sizeof(long unsigned));
		compressor->N = malloc((compressor->size + 1) * sizeof(int));
		if(!compressor->seq || !compressor->N) {
			ERROR();
		}
		compressor->N[0] = 0;
	}
	
	sfread(compressor->seq, sizeof(long unsigned), compressor->complen, file);
	sfread(compressor->N, sizeof(int), 1, file);
	sfread(compressor->N + 1, sizeof(int), compressor->N[0], file);
	
	return 1;
}

void dallocComp(CompDNA *compressor, unsigned size) {
	
	compressor->seqlen = 0;
	compressor->size = size;
	compressor->complen = 0;
	free(compressor->N);
	free(compressor->seq);
	if(size & 31) {
		compressor->size = (size >> 5) + 1;
		compressor->size <<= 5;
	} else {
		compressor->size = size;
	}
	compressor->seq = calloc(compressor->size >> 5, sizeof(long unsigned));
	if(!compressor->seq) {
		ERROR();
	}
	compressor->N = smalloc((compressor->size + 1) * sizeof(int));
	compressor->N[0] = 0;
}
