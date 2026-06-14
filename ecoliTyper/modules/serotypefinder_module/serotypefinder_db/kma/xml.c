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
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "nw.h"
#include "penalties.h"
#include "pherror.h"
#include "threader.h"
#include "version.h"
#include "xml.h"

void initXML(FILE *out, const char *templatefilename, const unsigned totFrags, int argc, char **argv) {
	fprintf(out, "<?xml version=\"1.0\"?>\n");
	fprintf(out, "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">\n");
	fprintf(out, "<BlastOutput>\n");
	fprintf(out, "\t<BlastOutput_program>kma</BlastOutput_program>\n");
	fprintf(out, "\t<BlastOutput_version>KMA %s</BlastOutput_version>\n", KMA_VERSION);
	fprintf(out, "\t<BlastOutput_reference>Philip T.L.C. Clausen, Frank M. Aarestrup & Ole Lund, \"Rapid and precise alignment of raw reads against redundant databases with KMA\", BMC Bioinformatics, 2018;19:307.</BlastOutput_reference>\n");
	fprintf(out, "\t<BlastOutput_db>%s</BlastOutput_db>\n", templatefilename);
	fprintf(out, "\t<BlastOutput_query-ID>kma-%s-%ld</BlastOutput_query-ID>\n", templatefilename, time(NULL));
	fprintf(out, "\t<BlastOutput_query-def>%s</BlastOutput_query-def>\n", "nucl");
	fprintf(out, "\t<BlastOutput_query-len>%d</BlastOutput_query-len>\n", totFrags);
	fprintf(out, "\t<BlastOutput_param>\n");
	fprintf(out, "\t\t<Parameters>\n");
	fprintf(out, "\t\t\t<Parameters_cmd>%s", *argv);
	while(--argc) {
		fprintf(out, " %s", *++argv);
	}
	fprintf(out, "</Parameters_cmd>\n");
	fprintf(out, "\t\t</Parameters>\n");
	fprintf(out, "\t</BlastOutput_param>\n");
	fprintf(out, "<BlastOutput_iterations>\n");
}

void capXML(FILE *out) {
	
	fprintf(out, "</BlastOutput_iterations>\n");
	fprintf(out, "</BlastOutput>\n");
}

FILE * openInitXML(const char *filename, const char *templatefilename, const unsigned totFrags, int argc, char **argv) {
	
	FILE *fileP;
	
	if(*filename == '-' && filename[1] == '-' && filename[2] == 0) {
		fileP = stdout;
	} else {
		fileP = sfopen(filename, "wb");
	}
	
	initXML(fileP, templatefilename, totFrags, argc, argv);
	
	return fileP;
}

void closeCapXML(FILE *out) {
	
	capXML(out);
	if(out && out != stdout) {
		fclose(out);
	}
}

void newIterXML(FILE *out, const int template, const int t_len, const char *template_name) {
	
	fprintf(out, "<Iteration>\n");
	fprintf(out, "\t<Iteration_iter-num>%d</Iteration_iter-num>\n", template);
	fprintf(out, "\t<Iteration_query-ID>Query_%d</Iteration_query-ID>\n", template);
	fprintf(out, "\t<Iteration_query-def>%s</Iteration_query-def>\n", template_name);
	fprintf(out, "\t<Iteration_query-len>%d</Iteration_query-len>\n", t_len);
	fprintf(out, "<Iteration_hits>\n");
}

double getEntropy(const unsigned char *aligned_assem_q, const int len) {
	
	unsigned i, *cPtr, counts[256];
	unsigned char *qPtr;
	double p, h;
	
	if(len == 0) {
		return 0;
	}
	
	/* init counts */
	cPtr = counts - 1;
	i = 257;
	while(--i) {
		*++cPtr = 0;
	}
	
	/* get counts */
	qPtr = (unsigned char *) aligned_assem_q + len;
	i = len + 1;
	while(--i) {
		++counts[*--qPtr];
	}
	
	/* calculate entropy */
	h = 0;
	i = 257;
	cPtr = counts - 1;
	while(--i) {
		if((p = *++cPtr)) {
			p /= len;
			h -= p * log2(p);
		}
	}
	
	return h;
}

void capIterXML(FILE *out, const int DB_size, const long unsigned seqsize, const int t_len, const int readCounts, const double p_value, const long read_score, const unsigned char *aligned_assem_q, const int len) {
	
	fprintf(out, "</Iteration_hits>\n");
	fprintf(out, "\t<Iteration_stat>\n");
	fprintf(out, "\t\t<Statistics>\n");
	fprintf(out, "\t\t\t<Statistics_db-num>%d</Statistics_db-num>\n", DB_size);
	fprintf(out, "\t\t\t<Statistics_db-len>%ld</Statistics_db-len>\n", seqsize);
	fprintf(out, "\t\t\t<Statistics_hsp-len>%d</Statistics_hsp-len>\n", readCounts);
	fprintf(out, "\t\t\t<Statistics_eff-space>%lu</Statistics_eff-space>\n", seqsize * t_len);
	fprintf(out, "\t\t\t<Statistics_kappa>%4.1e</Statistics_kappa>\n", p_value * read_score);
	fprintf(out, "\t\t\t<Statistics_lambda>%4.1e</Statistics_lambda>\n", p_value);
	fprintf(out, "\t\t\t<Statistics_entropy>%f</Statistics_entropy>\n", getEntropy(aligned_assem_q, len));
	fprintf(out, "\t\t</Statistics>\n");
	fprintf(out, "\t</Iteration_stat>\n");
	fprintf(out, "</Iteration>\n");
}

void hitXML(FILE *out, const int template, const unsigned char *template_name, const Aln *aligned, const AlnScore *alnStat, const Penalties *rewards, const int flag) {
	
	static volatile int Lock = 0;
	static int num = 0;
	volatile int *lock = &Lock;
	int i, Ms, MMs, W1s, Us, gap, pos, **d;
	unsigned char *t, *q;
	char *s, bases[6] = "ACGTN-";
	
	/* get stats */
	d = rewards->d;
	Ms = 0;
	MMs = 0;
	W1s = 0;
	Us = 0;
	gap = 0;
	pos = 0;
	t = aligned->t - 1;
	s = aligned->s - 1;
	q = aligned->q;
	i = aligned->len + 1;
	while(--i) {
		if(*++s == '_') {
			*s = ' ';
			if(*++t == 5 || *q == 5) {
				if(gap) {
					++Us;
				} else {
					++W1s;
					gap = 1;
				}
			} else {
				++MMs;
				if(0 < d[*t][*q]) {
					++pos;
				}
				gap = 1;
			}
		} else {
			++Ms;
			++t;
			if(0 < d[*t][*q]) {
				++pos;
			}
			gap = 1;
		}
		*t = bases[*t];
		*q = bases[*q];
		++q;
	}
	pos += W1s * (0 < rewards->W1) + Us * (0 < rewards->U);
	
	/* print stats */
	lock(lock);
	fprintf(out, "<Hit>\n");
	fprintf(out, "\t<Hit_num>%d</Hit_num>\n", ++num);
	fprintf(out, "\t<Hit_id>gnl|BL_ORD_ID|%d</Hit_id>\n", template + 1);
	fprintf(out, "\t<Hit_def>%s</Hit_def>\n", template_name);
	fprintf(out, "\t<Hit_accession>%d</Hit_accession>\n", template);
	fprintf(out, "\t<Hit_len>%d</Hit_len>\n", aligned->len);
	fprintf(out, "\t<Hit_hsps>\n");
	fprintf(out, "\t\t<Hsp>\n");
	fprintf(out, "\t\t\t<Hsp_num>1</Hsp_num>\n");
	fprintf(out, "\t\t\t<Hsp_bit-score>%d</Hsp_bit-score>\n", aligned->score);
	fprintf(out, "\t\t\t<Hsp_score>%d</Hsp_score>\n", aligned->mapQ);
	fprintf(out, "\t\t\t<Hsp_evalue>%f</Hsp_evalue>\n", pow(10, aligned->mapQ / (-10.0)));
	fprintf(out, "\t\t\t<Hsp_query-from>%d</Hsp_query-from>\n", ((flag & 16) ? (aligned->end) : (aligned->start)) + 1);
	fprintf(out, "\t\t\t<Hsp_query-to>%d</Hsp_query-to>\n", ((flag & 16) ? (aligned->start) : (aligned->end)) + 1);
	fprintf(out, "\t\t\t<Hsp_hit-from>%d</Hsp_hit-from>\n", alnStat->pos + 1);
	fprintf(out, "\t\t\t<Hsp_hit-to>%d</Hsp_hit-to>\n", alnStat->pos + alnStat->len - alnStat->tGaps + 1);
	fprintf(out, "\t\t\t<Hsp_query-frame>%d</Hsp_query-frame>\n", aligned->start % 3);
	fprintf(out, "\t\t\t<Hsp_hit-frame>%d</Hsp_hit-frame>\n", alnStat->pos % 3);
	fprintf(out, "\t\t\t<Hsp_identity>%d</Hsp_identity>\n", Ms);
	fprintf(out, "\t\t\t<Hsp_positive>%d</Hsp_positive>\n", pos);
	fprintf(out, "\t\t\t<Hsp_gaps>%d</Hsp_gaps>\n", W1s + Us);
	fprintf(out, "\t\t\t<Hsp_align-len>%d</Hsp_align-len>\n", aligned->len);
	fprintf(out, "\t\t\t<Hsp_qseq>%s</Hsp_qseq>\n", aligned->q);
	fprintf(out, "\t\t\t<Hsp_hseq>%s</Hsp_hseq>\n", aligned->t);
	fprintf(out, "\t\t\t<Hsp_midline>%s</Hsp_midline>\n", aligned->s);
	fprintf(out, "\t\t</Hsp>\n");
	fprintf(out, "\t</Hit_hsps>\n");
	fprintf(out, "</Hit>\n");
	unlock(lock);
}
