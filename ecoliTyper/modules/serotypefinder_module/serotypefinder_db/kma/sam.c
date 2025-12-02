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
#include <stdio.h>
#include <stdlib.h>
#include "nw.h"
#include "pherror.h"
#include "qseqs.h"
#include "runkma.h"
#include "sam.h"
#include "threader.h"
#include "version.h"

char * makeCigar(Qseqs *Cigar, const Aln *aligned) {
	
	int len, cLen, rep;
	char op, pop, *s, *cigar;
	unsigned char *t, *q;
	
	if(Cigar->size < (aligned->len << 1) + 20) {
		Cigar->size = (aligned->len << 1) + 20;
		free(Cigar->seq);
		Cigar->seq = smalloc(Cigar->size);
	} else if(aligned->len == 0) {
		return 0;
	}
	
	len = aligned->len;
	t = aligned->t;
	s = aligned->s;
	q = aligned->q;
	cigar = (char *) Cigar->seq;
	
	if(aligned->start) {
		cLen = sprintf(cigar, "%dS", aligned->start);
	} else {
		cLen = 0;
	}
	
	rep = 1;
	if(*s == '|') {
		pop = '=';
	} else if(*t == 5) {
		pop = 'I';
	} else if(*q == 5) {
		pop = 'D';
	} else {
		pop = 'X';
	}
	++t;
	++s;
	++q;
	while(--len) {
		if(*s == '|') {
			op = '=';
		} else if(*t == 5) {
			op = 'I';
		} else if(*q == 5) {
			op = 'D';
		} else {
			op = 'X';
		}
		if(op == pop) {
			++rep;
		} else {
			cLen += sprintf(cigar + cLen, "%d%c", rep, pop);
			rep = 1;
			pop = op;
		}
		++t;
		++s;
		++q;
	}
	cLen += sprintf(cigar + cLen, "%d%c", rep, pop);
	
	if(aligned->end) {
		cLen += sprintf(cigar + cLen, "%dS", aligned->end);
	}
	Cigar->len = cLen;
	
	return cigar;
}

void saminit(Qseqs *template_name, FILE *name_file, int *template_lengths, int DB_size, char *cmd) {
	
	fprintf(stdout, "@HD\tVN:1.6\tGO:reference\n");
	if(cmd) {
		fprintf(stdout, "@PG\tID:KMA\tPN:kma\tVN:%s\tCL:%s\n", KMA_VERSION, cmd);
	} else {
		fprintf(stdout, "@PG\tID:KMA\tPN:kma\tVN:%s\n", KMA_VERSION);
	}
	while(--DB_size) {
		fprintf(stdout, "@SQ\tSN:%s\tLN:%d\n", nameLoad(template_name, name_file), *++template_lengths);
	}
	sfseek(name_file, 0, SEEK_SET);
}

int samwrite(const Qseqs *qseq, const Qseqs *header, const Qseqs *Qual, char *rname, const Aln *aligned, const int *stats) {
	
	static volatile int Lock = 0;
	volatile int *lock = &Lock;
	static Qseqs *Cigar = 0;
	int flag, pos, mapQ, pnext, tlen, size, et, score, tab;
	char *qname, *cigar, *rnext, *qual;
	unsigned char *seq;
	
	/* flag */
	/*
	1		read paired
	2		read mapped in proper pair
	4		read unmapped
	8		mate unmapped
	16		read reverse strand
	32		mate reverse strand
	64		first in pair
	128		second in pair
	256		not primary alignment
	512		read fails platform/vendor quality checks
	1024	read is PCR or optical duplicate
	2048	supplementary alignment
	
	1		template having multiple segments in sequencing
	2		each segment properly aligned according to the aligner
	4		segment unmapped
	8		next segment in the template unmapped
	16		SEQ being reverse complemented
	32		SEQ of the next segment in the template being reverse complemented
	64		the first segment in the template
	128		the last segment in the template
	256		secondary alignment
	512		not passing filters, such as platform/vendor quality controls
	1024	PCR or optical duplicate
	2048	supplementary alignment
	*/
	
	
	qname = (char *) header->seq;
	seq = qseq->seq;
	if(Qual) {
		qual = (char *) Qual->seq;
	} else {
		qual = "*";
	}
	if(aligned) {
		mapQ = 254 < aligned->mapQ ? 254 : aligned->mapQ;
		et = *stats;
		score = stats[1];
		pos = stats[2] + 1;
		tlen = stats[3] - pos;
		flag = stats[4];
	} else {
		mapQ = 0;
		et = 0;
		score = 0;
		pos = 0;
		tlen = 0;
		et = *stats;
		flag = stats[1];
		if(rname == 0) {
			rname = "*";
		}
		cigar = "*";
	}
	rnext = "*";
	pnext = 0;
	
	tab = 0;
	if(qname) {
		while(*qname) {
			if(*qname == '\t') {
				tab = -tab;
				*qname = 0;
			} else {
				++qname;
				++tab;
			}
		}
		qname = (char *) header->seq;
	}
	
	lock(lock);
	if(Cigar == 0) {
		Cigar = setQseqs(256);
	}
	if(aligned) {
		cigar = makeCigar(Cigar, aligned);
	}
	size = fprintf(stdout, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tET:i:%d\tAS:i:%d\n", qname, flag, rname, pos, mapQ, cigar, rnext, pnext, tlen, (char *) seq, qual, et, score);
	unlock(lock);
	if(tab < 0) {
		qname[-tab] = '\t';
	}
	
	return size;
}
