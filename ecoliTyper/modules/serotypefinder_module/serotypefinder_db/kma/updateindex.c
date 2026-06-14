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
#include <stdio.h>
#include <stdlib.h>
#include "compdna.h"
#include "hashmap.h"
#include "hashmapcci.h"
#include "pherror.h"
#include "qualcheck.h"
#include "stdnuc.h"
#include "updateindex.h"

int (*update_DB)(HashMap *, CompDNA *, unsigned, int, double, double, unsigned *, unsigned *, Qseqs *) = &updateDBs;
void (*updateAnnotsPtr)(CompDNA *, int, int, FILE *, unsigned **, unsigned **, unsigned **) = &updateAnnots;

int updateDBs(HashMap *templates, CompDNA *qseq, unsigned template, int MinKlen, double homQ, double homT, unsigned *template_ulengths, unsigned *template_slengths, Qseqs *header) {
	
	int i, j, end, shifter, mPos, hLen, seqend;
	unsigned kmersize, mlen, flag, cPos, iPos, n;
	long unsigned mask, mmask, kmer, cmer, hmer, *seq;
	
	if(qseq->seqlen < templates->kmersize) {
		return 0;
	}
	
	/* set parameters */
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	seq = qseq->seq;
	kmersize = templates->kmersize;
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = templates->flag;
	seqend = qseq->seqlen - kmersize + 1;
	hLen = kmersize;
	n = 0;
	
	/* iterate sequence */
	for(i = 1, j = 0; i <= qseq->N[0] && j < seqend; ++i) {
		/* init k-mer */
		getKmer_macro(kmer, seq, j, cPos, iPos, (shifter + 2));
		cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
		end = qseq->N[i];
		for(j += kmersize - 1; j < end; ++j) {
			/* update k-mer */
			kmer = updateKmer_macro(kmer, seq, j, mask);
			cmer = flag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
			
			/* update hashMap */
			hashMap_add(templates, cmer, template);
			++n;
		}
		j = end + 1;
	}
	qseq->N[0]--;
	
	return n;
}

int updateDBs_sparse(HashMap *templates, CompDNA *qseq, unsigned template, int MinKlen, double homQ, double homT, unsigned *template_ulengths, unsigned *template_slengths, Qseqs *header) {
	
	int i, j, end, rc, prefix_len, prefix_shifter, mPos, hLen, seqend;
	unsigned kmersize, shifter, mlen, flag, cPos, iPos, slen, ulen;
	long unsigned pmask, mask, mmask, prefix, pmer, kmer, cmer, hmer, *seq;
	
	if(qseq->seqlen < templates->kmersize) {
		return 0;
	}
	seq = qseq->seq;
	prefix = templates->prefix;
	prefix_len = templates->prefix_len;
	prefix_shifter = 64 - (prefix_len << 1);
	pmask = 0xFFFFFFFFFFFFFFFF >> prefix_shifter;
	kmersize = templates->kmersize;
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = templates->flag;
	seqend = qseq->seqlen - kmersize - prefix_len + 1;
	hLen = kmersize;
	
	/* test homology and length */
	if(QualCheck(templates, qseq, MinKlen, homQ, homT, template_ulengths, header)) {
		slen = 0;
		ulen = 0;
		for(rc = 0; rc < 2; ++rc) {
			/* revers complement */
			if(rc) {
				comp_rc(qseq);
			}
			
			/* iterate seq */
			qseq->N[0]++;
			qseq->N[qseq->N[0]] = qseq->seqlen;
			j = 0;
			if(prefix_len) {
				for(i = 1; i <= qseq->N[0] && j < seqend; ++i) {
					getKmer_macro(pmer, seq, j, cPos, iPos, (prefix_shifter + 2));
					end = qseq->N[i] - kmersize;
					for(j += prefix_len - 1; j < end; ++j) {
						pmer = updateKmer_macro(pmer, seq, j, pmask);
						if(pmer == prefix) {
							/* get kmer */
							getKmer_macro(kmer, seq, (j + 1), cPos, iPos, shifter);
							cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
							
							/* add kmer */
							if(hashMap_add(templates, cmer, template)) {
								++ulen;
							}
							++slen;
						}
					}
					j = qseq->N[i] + 1;
				}
				qseq->N[0]--;
			} else {
				for(i = 1; i <= qseq->N[0] && j < seqend; ++i) {
					getKmer_macro(kmer, seq, j, cPos, iPos, (shifter + 2));
					cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
					end = qseq->N[i];
					for(j += kmersize - 1;j < end; ++j) {
						/* update kmer */
						kmer = updateKmer_macro(kmer, seq, j, mask);
						cmer = flag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
						
						/* add kmer */
						if(hashMap_add(templates, cmer, template)) {
							++ulen;
						}
						++slen;
					}
					j = end + 1;
				}
				qseq->N[0]--;
			}
		}
		if(prefix_len == 0 && !prefix) {
			comp_rc(qseq);
		}
		template_slengths[template] = slen;
		template_ulengths[template] = ulen;
		return slen;
	}
	
	return 0;
}

void updateAnnots(CompDNA *qseq, int DB_size, int kmerindex, FILE *seq_out, unsigned **template_lengths, unsigned **template_ulengths, unsigned **template_slengths) {
	
	/* Dump annots */
	cfwrite(qseq->seq, sizeof(long unsigned), (qseq->seqlen >> 5) + 1, seq_out);
	
	(*template_lengths)[DB_size] = qseq->seqlen;
	if((DB_size + 1) >= **template_lengths) {
		**template_lengths <<= 1;
		*template_lengths = realloc(*template_lengths, **template_lengths * sizeof(unsigned));
		if(!*template_lengths) {
			ERROR();
		}
	}
}

void updateAnnots_sparse(CompDNA *qseq, int DB_size, int kmerindex, FILE *seq_out, unsigned **template_lengths, unsigned **template_ulengths, unsigned **template_slengths) {
	
	/* Dump annots */
	cfwrite(qseq->seq, sizeof(long unsigned), (qseq->seqlen >> 5) + 1, seq_out);
	
	(*template_lengths)[DB_size] = qseq->seqlen;
	if((DB_size + 1) >= **template_ulengths) {
		**template_ulengths <<= 1;
		*template_slengths = realloc(*template_slengths, **template_ulengths * sizeof(unsigned));
		*template_ulengths = realloc(*template_ulengths, **template_ulengths * sizeof(unsigned));
		*template_lengths = realloc(*template_lengths, **template_ulengths * sizeof(unsigned));
		if(!template_lengths || !template_slengths || !template_ulengths) {
			ERROR();
		}
	}
}
