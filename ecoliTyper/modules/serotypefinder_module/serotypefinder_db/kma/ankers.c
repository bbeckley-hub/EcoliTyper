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
#include <stdlib.h>
#include "ankers.h"
#include "compdna.h"
#include "pherror.h"
#include "qseqs.h"

int (*printPtr)(int*, CompDNA*, int, const Qseqs*, const int, FILE *out) = &print_ankers;
int (*printPairPtr)(int*, CompDNA*, int, const Qseqs*, CompDNA*, int, const Qseqs*, const int flag, const int flag_r, FILE *out) = &printPair;
int (*deConPrintPtr)(int*, CompDNA*, int, const Qseqs*, const int flag, FILE *out) = &print_ankers;

int print_ankers(int *out_Tem, CompDNA *qseq, int rc_flag, const Qseqs *header, const int flag, FILE *out) {
	
	int infoSize[7];
	
	infoSize[0] = qseq->seqlen;
	infoSize[1] = qseq->complen;
	infoSize[2] = *(qseq->N);
	infoSize[3] = rc_flag;
	infoSize[4] = *out_Tem;
	infoSize[5] = header->len;
	infoSize[6] = flag;
	sfwrite(infoSize, sizeof(int), 7, out);
	sfwrite(qseq->seq, sizeof(long unsigned), qseq->complen, out);
	if(*(qseq->N)) {
		sfwrite(qseq->N + 1, sizeof(int), *(qseq->N), out);
	}
	sfwrite(out_Tem + 1, sizeof(int), *out_Tem, out);
	sfwrite(header->seq, 1, header->len, out);
	
	return 0;
}

int print_ankers_Sparse(int *out_Tem, CompDNA *qseq, int rc_flag, const Qseqs *header, const int flag, FILE *out) {
	
	int infoSize[7];
	
	infoSize[0] = qseq->seqlen;
	infoSize[1] = qseq->complen;
	infoSize[2] = *(qseq->N);
	infoSize[3] = rc_flag < 0 ? rc_flag : -rc_flag;
	infoSize[4] = *out_Tem;
	infoSize[5] = header->len;
	infoSize[6] = flag;
	sfwrite(infoSize, sizeof(int), 7, out);
	sfwrite(qseq->seq, sizeof(long unsigned), qseq->complen, out);
	if(*(qseq->N)) {
		sfwrite(qseq->N + 1, sizeof(int), *(qseq->N), out);
	}
	sfwrite(out_Tem + 1, sizeof(int), *out_Tem, out);
	sfwrite(header->seq, 1, header->len, out);
	
	return 0;
}

int find_contamination(int *out_Tem, const int contamination) {
	
	int i;
	
	i = *out_Tem + 1;
	out_Tem += i;
	while(--i) {
		if(*--out_Tem == contamination) {
			return i;
		}
	}
	
	return 0;
}

int find_contamination2(int *out_Tem, const int contamination) {
	
	int i;
	
	i = *out_Tem + 1;
	out_Tem += i;
	while(--i) {
		if(*--out_Tem == contamination) {
			return i;
		} else if(0 < *out_Tem) {
			return 0;
		}
	}
	
	return 0;
}

int deConPrint(int *out_Tem, CompDNA *qseq, int rc_flag, const Qseqs *header, const int flag, FILE *out) {
	
	int contPos;
	
	if((contPos = find_contamination(out_Tem, out_Tem[-3])) != 0) {
		out_Tem[contPos] = out_Tem[*out_Tem];
		--*out_Tem;
	}
	if((contPos = find_contamination2(out_Tem, -out_Tem[-3])) != 0) {
		out_Tem[contPos] = out_Tem[*out_Tem];
		--*out_Tem;
	}
	
	if(0 < *out_Tem) {
		return printPtr(out_Tem, qseq, rc_flag, header, flag, out);
	}
	
	return 1;
}

int deConPrintPair(int *out_Tem, CompDNA *qseq, int bestScore, const Qseqs *header, CompDNA *qseq_r, int bestScore_r, const Qseqs *header_r, const int flag, const int flag_r, FILE *out) {
	
	int contPos;
	
	if((contPos = find_contamination(out_Tem, out_Tem[-3])) != 0) {
		out_Tem[contPos] = out_Tem[*out_Tem];
		--*out_Tem;
	}
	if((contPos = find_contamination2(out_Tem, -out_Tem[-3])) != 0) {
		out_Tem[contPos] = out_Tem[*out_Tem];
		--*out_Tem;
	}
	
	if(0 < *out_Tem) {
		contPos = *out_Tem;
		*out_Tem = 0;
		printPtr(out_Tem, qseq, bestScore, header, flag, out);
		*out_Tem = contPos;
		return printPtr(out_Tem, qseq_r, bestScore_r, header_r, flag_r, out);
	}
	
	return 1;
}

int printPair(int *out_Tem, CompDNA *qseq, int bestScore, const Qseqs *header, CompDNA *qseq_r, int bestScore_r, const Qseqs *header_r, const int flag, const int flag_r, FILE *out) {
	
	int contPos;
	
	contPos = *out_Tem;
	*out_Tem = 0;
	printPtr(out_Tem, qseq, bestScore, header, flag, out);
	*out_Tem = contPos;
	printPtr(out_Tem, qseq_r, bestScore_r, header_r, flag_r, out);
	
	return 0;
}

int get_ankers(int *out_Tem, CompDNA *qseq, Qseqs *header, int *flag, FILE *inputfile) {
	
	static int infoSize[7] = {0, 0, 0, 0, 0, 0, 0};
	
	if(infoSize[0] < 0) {
		*out_Tem = -infoSize[0];
		return 0;
	} else if(fread(infoSize, sizeof(int), 7, inputfile) == 7) {
		if(infoSize[0] < 0) {
			*out_Tem = -infoSize[0];
			return 0;
		}
		qseq->seqlen = infoSize[0];
		qseq->complen = infoSize[1];
		*out_Tem = infoSize[4];
		header->len = infoSize[5];
		*flag = infoSize[6];
		
		/* reallocate */
		if(qseq->size <= qseq->seqlen) {
			free(qseq->N);
			free(qseq->seq);
			if(qseq->seqlen & 31) {
				qseq->size = (qseq->seqlen >> 5) + 1;
				qseq->size <<= 6;
			} else {
				qseq->size = qseq->seqlen << 1;
			}
			
			qseq->seq = calloc(qseq->size >> 5, sizeof(long unsigned));
			qseq->N = malloc((qseq->size + 1) * sizeof(int));
			if(!qseq->seq || !qseq->N) {
				ERROR();
			}
		}
		
		qseq->N[0] = infoSize[2];
		if(header->size <= header->len) {
			free(header->seq);
			header->size = header->len << 1;
			header->seq = malloc(header->size);
			if(!header->seq) {
				ERROR();
			}
		}
		
		sfread(qseq->seq, sizeof(long unsigned), qseq->complen, inputfile);
		sfread(qseq->N + 1, sizeof(int), qseq->N[0], inputfile);
		sfread(out_Tem + 1, sizeof(int), *out_Tem, inputfile);
		sfread(header->seq, 1, header->len, inputfile);
	} else {
		*out_Tem = -infoSize[0];
		return 0;
	}
	
	/* return score */
	return infoSize[3];
}
