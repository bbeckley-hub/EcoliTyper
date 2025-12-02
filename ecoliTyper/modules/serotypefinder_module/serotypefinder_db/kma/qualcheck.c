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
#include <limits.h>
#include <stdlib.h>
#include "compdna.h"
#include "hashmap.h"
#include "pherror.h"
#include "qualcheck.h"
#include "stdnuc.h"
#include "stdstat.h"

int (*QualCheck)(HashMap *templates, CompDNA *, int, double, double, unsigned *, Qseqs *) = &lengthCheck;

int lengthCheck(HashMap *templates, CompDNA *qseq, int MinKlen, double homQ, double homT, unsigned *template_ulengths, Qseqs *header) {
	
	int i, j, end, rc, thisKlen, prefix_len, prefix_shifter, cPos, iPos;
	long unsigned prefix, kmer;
	
	thisKlen = MinKlen;
	prefix = templates->prefix;
	prefix_len = templates->prefix_len;
	prefix_shifter = sizeof(long unsigned) * sizeof(long unsigned) - (prefix_len << 1);
	
	if(qseq->seqlen < templates->kmersize) {
		return 0;
	} else if(prefix_len == 0) {
		if((qseq->seqlen - templates->kmersize + 1) * 2 < MinKlen) {
			return 0;
		} else {
			comp_rc(qseq);
			return 1;
		}
	}
	
	for(rc = 0; rc < 2; ++rc) {
		/* revers complement */
		if(rc) {
			comp_rc(qseq);
		}
		
		/* iterate seq */
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1, j = 0; i <= qseq->N[0] && thisKlen != 0; ++i) {
			end = qseq->N[i] - prefix_len - templates->kmersize + 1;
			for(;j < end && thisKlen != 0; ++j) {
				getKmer_macro(kmer, qseq->seq, j, cPos, iPos, prefix_shifter);
				if(kmer == prefix) {
					--thisKlen;
				}
			}
			j = qseq->N[i] + 1;
		}
		qseq->N[0]--;
	}
	
	if(thisKlen) {
		return 0;
	} else {
		return 1;
	}
}

int queryCheck(HashMap *templates, CompDNA *qseq, int MinKlen, double homQ, double homT, unsigned *template_ulengths, Qseqs *header) {
	
	static unsigned *Scores_tot = 0, *bestTemplates = 0;
	int i, j, end, rc, template, cPos, iPos, mlen, flag, mPos, hLen;
	unsigned thisKlen, prefix_len, prefix_shifter, shifter, DB_size, *values;
	double bestQ, thisQ;
	long unsigned mmask, prefix, kmer, cmer;
	void (*updateScoreAndTemplate_ptr)(unsigned *, unsigned *, unsigned *);
	
	thisKlen = 0;
	prefix_len = templates->prefix_len;
	if(qseq->seqlen < templates->kmersize + prefix_len) {
		return 0;
	} else if(prefix_len == 0 && templates->prefix != 0) {
		prefix = 0;
	} else {
		prefix = templates->prefix;
	}
	prefix_shifter = sizeof(long unsigned) * sizeof(long unsigned) - (prefix_len << 1);
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (templates->kmersize << 1);
	DB_size = templates->DB_size;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = templates->flag;
	hLen = templates->kmersize;
	
	/* realloc */
	if(!bestTemplates) {
		bestTemplates = smalloc(1024 * sizeof(unsigned));
		Scores_tot = calloc(1024, sizeof(unsigned));
		if(!Scores_tot) {
			ERROR();
		}
		*Scores_tot = 1024;
	} else if(DB_size >= *Scores_tot) {
		free(Scores_tot);
		Scores_tot = calloc(2 * DB_size, sizeof(unsigned));
		free(bestTemplates);
		bestTemplates = malloc(2 * DB_size * sizeof(unsigned));
		if(!Scores_tot || !bestTemplates) {
			ERROR();
		}
		*Scores_tot = 2 * DB_size;
	}
	
	if(DB_size < USHRT_MAX) {
		updateScoreAndTemplate_ptr = &updateScoreAndTemplateHU;
	} else {
		updateScoreAndTemplate_ptr = &updateScoreAndTemplate;
	}
	
	/* get scores */
	bestTemplates[0] = 0;
	for(rc = 0; rc < 2; ++rc) {
		/* revers complement */
		if(rc) {
			comp_rc(qseq);
		}
		
		/* iterate seq */
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1, j = 0; i <= qseq->N[0]; ++i) {
			end = qseq->N[i] - prefix_len - templates->kmersize + 1;
			for(;j < end; ++j) {
				if(prefix_len == 0 || getKmer(qseq->seq, j, prefix_shifter) == prefix) {
					++thisKlen;
					getKmer_macro(kmer, qseq->seq, (j + prefix_len), cPos, iPos, shifter);
					cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
					
					if((values = hashMapGet(templates, cmer))) {
						updateScoreAndTemplate_ptr(Scores_tot, bestTemplates, values);
					}
				}
			}
			j = qseq->N[i] + 1;
		}
		
		qseq->N[0]--;
	}
	
	if(thisKlen < MinKlen) {
		for(i = 1; i <= *bestTemplates; ++i) {
			Scores_tot[bestTemplates[i]] = 0;
		}
		return 0;
	}
	
	/* get query cov */
	bestQ = 0;
	template = 0;
	for(i = 1; i <= *bestTemplates; ++i) {
		thisQ = 1.0 * Scores_tot[bestTemplates[i]] / thisKlen;
		if(bestQ < thisQ) {
			bestQ = thisQ;
			template = bestTemplates[i];
		}
		Scores_tot[bestTemplates[i]] = 0;
	}
	
	if(bestQ < homQ) {
		fprintf(stdout, "%s\t%d\t%f\t%d\n", header->seq + 1, DB_size, bestQ, template);
		return 1;
	} else {
		fprintf(stdout, "%s\t%d\t%f\t%d\n", header->seq + 1, template, bestQ, template);
		return 0;
	}
}

int templateCheck(HashMap *templates, CompDNA *qseq, int MinKlen, double homQ, double homT, unsigned *template_ulengths, Qseqs *header) {
	
	static unsigned *Scores = 0, *Scores_tot = 0, *bestTemplates = 0;
	static HashMap_kmers *foundKmers = 0;
	int i, j, end, rc, templateQ, templateT, cPos, iPos, mlen, mPos, hLen;
	unsigned thisKlen, prefix_len, prefix_shifter, shifter, DB_size, flag;
	unsigned *values;
	double bestQ, thisQ, bestT, thisT;
	long unsigned mmask, key, prefix, cmer;
	void (*updateScoreAndTemplate_ptr)(unsigned *, unsigned *, unsigned *);
	void (*addUscore_ptr)(unsigned *, unsigned *);
	
	thisKlen = 0;
	prefix_len = templates->prefix_len;
	if(qseq->seqlen < templates->kmersize + prefix_len) {
		return 0;
	} else if(prefix_len == 0 && templates->prefix != 0) {
		prefix = 0;
	} else {
		prefix = templates->prefix;
	}
	prefix_shifter = sizeof(long unsigned) * sizeof(long unsigned) - (prefix_len << 1);
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (templates->kmersize << 1);
	DB_size = templates->DB_size;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = templates->flag;
	hLen = templates->kmersize;
	
	if(DB_size < USHRT_MAX) {
		updateScoreAndTemplate_ptr = &updateScoreAndTemplateHU;
		addUscore_ptr = &addUscoreHU;
	} else {
		updateScoreAndTemplate_ptr = &updateScoreAndTemplate;
		addUscore_ptr = &addUscore;
	}
	
	if(foundKmers == 0) {
		foundKmers = smalloc(sizeof(HashMap_kmers));
		foundKmers->n = 0;
		foundKmers->size = templates->size;
		foundKmers->table = calloc(foundKmers->size + 1, sizeof(HashTable_kmers *));
		Scores = calloc(1024, sizeof(unsigned));
		Scores_tot = calloc(1024, sizeof(unsigned));
		bestTemplates = smalloc(1024 * sizeof(unsigned));
		if(!Scores || !Scores_tot || !foundKmers->table) {
			ERROR();
		}
		*Scores_tot = 1024;
	}
	
	/* realloc */
	if(DB_size >= *Scores_tot) {
		free(Scores);
		Scores = calloc(2 * DB_size, sizeof(unsigned));
		free(Scores_tot);
		Scores_tot = calloc(2 * DB_size, sizeof(unsigned));
		free(bestTemplates);
		bestTemplates = malloc(2 * DB_size * sizeof(unsigned));
		if(!Scores || !Scores_tot || !bestTemplates) {
			ERROR();
		}
		Scores_tot[0] = 2 * DB_size;
	}
	
	/* get scores */
	bestTemplates[0] = 0;
	for(rc = 0; rc < 2; ++rc) {
		/* revers complement */
		if(rc) {
			comp_rc(qseq);
		}
		
		/* iterate seq */
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1, j = 0; i <= qseq->N[0]; ++i) {
			end = qseq->N[i] - prefix_len - templates->kmersize + 1;
			for(;j < end; ++j) {
				if(prefix_len == 0 || getKmer(qseq->seq, j, prefix_shifter) == prefix) {
					++thisKlen;
					getKmer_macro(key, qseq->seq, (j + prefix_len), cPos, iPos, shifter);
					cmer = flag ? getCmer(key, &mPos, &hLen, shifter, mlen, mmask) : key;
					
					if((values = hashMapGet(templates, cmer))) {
						updateScoreAndTemplate_ptr(Scores_tot, bestTemplates, values);
						if(hashMap_CountKmer(foundKmers, cmer)) {
							addUscore_ptr(Scores, values);
						}
					}
				}
			}
			j = qseq->N[i] + 1;
		}
		qseq->N[0]--;
	}
	
	if(thisKlen < MinKlen) {
		return 0;
	}
	
	/* get query cov */
	bestQ = 0;
	templateQ = 0;
	bestT = 0;
	templateT = 0;
	for(i = 1; i <= *bestTemplates; ++i) {
		thisQ = 1.0 * Scores_tot[bestTemplates[i]] / thisKlen;
		if(thisQ > bestQ) {
			bestQ = thisQ;
			templateQ = bestTemplates[i];
		}
		thisT = 1.0 * Scores[bestTemplates[i]] / template_ulengths[bestTemplates[i]];
		if(thisT > bestT) {
			bestT = thisT;
			templateT = bestTemplates[i];
		}
		Scores_tot[bestTemplates[i]] = 0;
		Scores[bestTemplates[i]] = 0;
	}
	
	emptyHash(foundKmers);
	
	if(thisKlen >= MinKlen) {
		if(cmp(bestT < homT, bestQ < homQ)) {
			fprintf(stdout, "%s\t%d\t%f\t%d\t%f\t%d\n", header->seq + 1, DB_size, bestQ, templateQ, bestT, templateT);
			return 1;
		} else {
			fprintf(stdout, "%s\t%d\t%f\t%d\t%f\t%d\n", header->seq + 1, (bestT < homT) ? templateQ : templateT, bestQ, templateQ, bestT, templateT);
			return 0;
		}
	} else {
		return 0;
	}
}

void updateScoreAndTemplate(unsigned *Scores_tot, unsigned *bestTemplates, unsigned *values) {
	
	unsigned i;
	
	i = *values + 1;
	while(--i) {
		if(Scores_tot[*++values] == 0) {
			bestTemplates[0]++;
			bestTemplates[*bestTemplates] = *values;
		}
		Scores_tot[*values]++;
	}
}

void updateScoreAndTemplateHU(unsigned *Scores_tot, unsigned *bestTemplates, unsigned *values_org) {
	
	unsigned i;
	short unsigned *values;
	
	values = (short unsigned *)(values_org);
	i = *values + 1;
	while(--i) {
		if(Scores_tot[*++values] == 0) {
			bestTemplates[0]++;
			bestTemplates[*bestTemplates] = *values;
		}
		Scores_tot[*values]++;
	}
}

void addUscore(unsigned *Scores, unsigned *values) {
	
	unsigned i;
	
	i = *values + 1;
	while(--i) {
		Scores[*++values]++;
	}
}

void addUscoreHU(unsigned *Scores, unsigned *values_org) {
	
	unsigned i;
	short unsigned *values;
	
	values = (short unsigned *)(values_org);
	i = *values + 1;
	while(--i) {
		Scores[*++values]++;
	}
}
