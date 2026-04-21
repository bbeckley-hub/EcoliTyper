/* Philip T.L.C. Clausen Jun 2025 plan@dtu.dk */

/*
 * Copyright (c) 2025, Philip Clausen, Technical University of Denmark
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

#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include "assembly.h"
#include "compdna.h"
#include "ef.h"
#include "hashmapkma.h"
#include "kmmap.h"
#include "pherror.h"
#include "qseqs.h"
#include "reassign.h"
#include "runkma.h"
#include "stdnuc.h"
#include "stdstat.h"

int reassign_template(HashMapKMA *templates, FILE *name_file, int seq_in, int *template_lengths, Assem *aligned_assem, AssemInfo *matrix, Qseqs *qseq, CompDNA *tseq) {
	
	static int *bestTemplates = 0, *Score = 0;
	static long *seq_indexes = 0, *name_indexes = 0;
	static unsigned char *to2Bit = 0, *complement = 0;
	static CompDNA *qseq_comp = 0, *qseq_rc_comp = 0;
	int i, template, match;
	char *templatefilename;
	FILE *templatefile;
	
	/* init */
	if(!name_file && !seq_in && !template_lengths && !aligned_assem && !matrix && !qseq && !tseq) { /* clean up */
		/* free */
		if(bestTemplates) {
			hashMapKMA_destroy(templates);
			free(to2Bit);
			free(complement);
			freeComp(qseq_comp);
			freeComp(qseq_rc_comp);
			free(qseq_comp);
			free(qseq_rc_comp);
			free(bestTemplates);
			free(Score);
			free(seq_indexes);
			free(name_indexes);
		} else if(templates) {
			free(templates);
		}
		
		/* reset */
		to2Bit = 0;
		complement = 0;
		qseq_comp = 0;
		qseq_rc_comp = 0;
		bestTemplates = 0;
		Score = 0;
		seq_indexes = 0;
		name_indexes = 0;
		
		return 0;
	} else if(!bestTemplates) { /* init */
		/* k-mer signatures */
		templatefilename = (char *)(templates->exist);
		i = strlen(templatefilename);
		strcat(templatefilename, ".comp.b");
		templatefile = sfopen(templatefilename, "rb" );
		if(templates->shmFlag & 32) {
			hashMapKMAmmap(templates, templatefile);
		} else if(hashMapKMA_load(templates, templatefile, templatefilename) == 1) {
			fprintf(stderr, "Wrong format of DB.\n");
			exit(1);
		}
		templatefilename[i] = 0;
		fclose(templatefile);
		
		/* arrays */
		to2Bit = base2nibble_table();
		complement = complement_table();
		qseq_comp = malloc(sizeof(CompDNA));
		qseq_rc_comp = malloc(sizeof(CompDNA));
		bestTemplates = malloc(((templates->DB_size + 1) << 1) * sizeof(int));
		Score = calloc(templates->DB_size, sizeof(unsigned));
		seq_indexes = malloc((templates->DB_size + 1) * sizeof(long));
		name_indexes = malloc((templates->DB_size + 1) * sizeof(long));
		if(!qseq_comp || !qseq_rc_comp || !bestTemplates || !Score || !seq_indexes || !name_indexes) {
			ERROR();
		}
		allocComp(qseq_comp, 1024);
		allocComp(qseq_rc_comp, 1024);
		
		/* make file indexes */
		*seq_indexes = seq_in;
		seq_indexes[1] = 0;
		for(i = 2; i < templates->DB_size; ++i) {
			seq_indexes[i] = seq_indexes[i - 1] + ((template_lengths[i - 1] >> 5) + 1) * sizeof(long unsigned);
		}
		*name_indexes = ftell(name_file);
		fseek(name_file, 0, SEEK_SET);
		name_indexes[1] = 0;
		name_indexes[2] = 0;
		i = 2;
		while(i < templates->DB_size) {
			name_indexes[i]++;
			if(fgetc(name_file) == '\n') {
				name_indexes[i + 1] = name_indexes[i];
				++i;
			}
		}
		fseek(name_file, *name_indexes, SEEK_SET);
	}
	
	/* Convert consensus to qseq */
	consensus2qseq(aligned_assem->q, aligned_assem->len, qseq, to2Bit);
	compDNA(qseq_comp, qseq->seq, qseq->len);
	rc_comp(qseq_comp, qseq_rc_comp);
	
	/* Find potential reassignment candidates */
	reassign_kmers(templates, bestTemplates, Score, template_lengths, qseq_comp, qseq_rc_comp);
	
	/* Get sequence(s) of candidates and cmp to consensus in order of declining lengths */
	reassign_buildheap(bestTemplates, template_lengths);
	
	/* Go over candidates */
	match = -1;
	while(match < 0 && (template = reassign_popheap(bestTemplates, template_lengths))) {
		/* load seq of candidate */
		reassign_loadseq(tseq, template_lengths[NORM(template)], seq_in, seq_indexes[NORM(template)]);
		
		/* match candidate with consensus */
		if(templates->prefix || templates->prefix_len) {
			/* check both strands */
			if((match = reassign_matchseqs(qseq_comp, tseq)) < 0) {
				match = reassign_matchseqs(qseq_rc_comp, tseq);
			}
		} else if(0 < template) {
			/* check forwrd */
			match = reassign_matchseqs(qseq_comp, tseq);
		} else {
			/* check rc */
			match = reassign_matchseqs(qseq_rc_comp, tseq);
		}
	}
	
	if(!template) { /* no matches found */
		return 0;
	} else if(template < 0) { /* match is rc */
		/* rc org match */
		assemble_rc(matrix, assem_rc(aligned_assem, complement));
		
		/* reset direction of matches */
		template = -template;
	}
	
	/* Readjust consensus alignment and matrix to fit reassigned template */
	reassign_matrix_offset(matrix, aligned_assem, match, tseq->seq, template_lengths[template]);
	
	/* Update statistics */
	getExtendedFeatures(aligned_assem, matrix, tseq->seq, template_lengths[template], 1);
	*name_indexes = ftell(name_file);
	fseek(name_file, name_indexes[template], SEEK_SET);
	nameLoad(qseq, name_file);
	fseek(name_file, *name_indexes, SEEK_SET);
	
	return template;
}

int consensus2qseq(unsigned char *src, int len, Qseqs *dest, const unsigned char *trans) {
	
	int n, nibble;
	unsigned char *seq;
	
	/* reallocate */
	if(dest->size < len++) {
		dest->seq = realloc(dest->seq, len);
		if(!dest->seq) {
			ERROR();
		}
		dest->size = len;
	}
	
	--src;
	seq = dest->seq - 1;
	n = 0;
	while(--len) {
		nibble = trans[*++src];
		if(nibble <= 4) {
			*++seq = nibble;
			++n;
		}
	}
	dest->len = n;
	
	return n;
}

int reassign_heapify(int *bests, int *lengths, int index) {
	
	int n, root, child;
	
	/* init */
	n = *bests++;
	root = index;
	
	/* test children */
	child = (index << 1) + 1;
	if(child < n && lengths[NORM(bests[root])] < lengths[NORM(bests[child])]) {
		root = child;
	}
	++child;
	if(child < n && lengths[NORM(bests[root])] < lengths[NORM(bests[child])]) {
		root = child;
	}
	
	/* test if root has changed */
	if(root != index) {
		/* exchange root with max child */
		child = bests[index];
		bests[index] = bests[root];
		bests[root] = child;
		
		/* heapify affected sub-tree */
		return 1 + reassign_heapify(bests, lengths, root);
	}
	
	return 0;
}

int reassign_buildheap(int *bests, int *lengths) {
	
	int i, n;
	
	n = 0; /* number of re-arrangements */
	i = *bests >> 1; /* first leaf */
	while(i--) { /* go over sub-trees */
		n += reassign_heapify(bests, lengths, i);
	}
	
	return n;
}

int reassign_popheap(int *bests, int *lengths) {
	
	int template;
	
	if(!*bests) {
		return 0;
	}
	
	/* get top element */
	template = bests[1];
	/* move bottom to top */
	bests[1] = bests[(*bests)--];
	/* heapify new root */
	reassign_heapify(bests, lengths, 0);
	
	return template;
}

int reassign_kmers(const HashMapKMA *templates, int *bestTemplates, int *Score, int *template_lengths, CompDNA *qseq_fw, CompDNA *qseq_rc) {
	
	int i, j, end, seqend, rc, kmersize, mPos, mlen, hLen, reps;
	int prefix_len, prefix_shifter, *bests, *Scores, *bestTemplates_r;
	unsigned n, SU, shifter, cPos, iPos, flag, *values, *last;
	long unsigned mask, mmask, pmask, prefix, pmer, kmer, cmer, hmer, *seq;
	short unsigned *values_s;
	CompDNA *qseq;
	
	SU = (templates->DB_size < USHRT_MAX) ? 1 : 0;
	kmersize = templates->kmersize;
	prefix = templates->prefix;
	prefix_len = templates->prefix_len;
	prefix_shifter = 64 - (prefix_len << 1);
	pmask = 0xFFFFFFFFFFFFFFFF >> prefix_shifter;
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = templates->flag;
	hLen = kmersize;
	qseq = qseq_fw;
	end = qseq->seqlen;
	seqend = end - kmersize + 1;
	seq = qseq->seq;
	bests = bestTemplates;
	Scores = Score;
	*bestTemplates = 0;
	bestTemplates_r = bestTemplates + templates->DB_size;
	*bestTemplates_r = 0;
	if(prefix_len) {
		template_lengths += templates->DB_size;
		for(rc = 0; rc < 2; ++rc) {
			if(rc) { /* go over reverse complement */
				qseq = qseq_rc;
				seq = qseq->seq;
			}
			
			/* iterate seq */
			last = 0;
			reps = 0;
			j = 0;
			qseq->N[0]++;
			qseq->N[qseq->N[0]] = qseq->seqlen;
			for(i = 1; i <= qseq->N[0]; ++i) {
				/* init prefix */
				getKmer_macro(pmer, seq, j, cPos, iPos, (prefix_shifter + 2));
				end = qseq->N[i] - kmersize;
				for(j += prefix_len - 1; j < end; ++j) {
					/* update prefix */
					pmer = updateKmer_macro(pmer, seq, j, pmask);
					if(pmer == prefix) {
						/* get kmer */
						getKmer_macro(kmer, seq, (j + 1), cPos, iPos, shifter);
						cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
						/* lookup */
						if((values = hashMap_get(templates, cmer))) {
							if(values == last) {
								++reps;
							} else if(last) {
								if(SU) {
									values_s = (short unsigned *) last;
									n = *values_s + 1;
									while(--n) {
										if((Scores[*++values_s] += reps) == reps) {
											bests[++*bests] = *values_s;
										}
									}
								} else {
									n = *last + 1;
									while(--n) {
										if((Scores[*++last] += reps) == reps) {
											bests[++*bests] = *last;
										}
									}
								}
								reps = 1;
								last = values;
							} else {
								reps = 1;
								last = values;
							}
						}
					}
				}
			}
			/* get last k-mer hits */
			if(last) {
				if(SU) {
					values_s = (short unsigned *) last;
					n = *values_s + 1;
					while(--n) {
						if((Scores[*++values_s] += reps) == reps) {
							bests[++*bests] = *values_s;
						}
					}
				} else {
					n = *last + 1;
					while(--n) {
						if((Scores[*++last] += reps) == reps) {
							bests[++*bests] = *last;
						}
					}
				}
			}
			
			reps = 0;
			qseq->N[0]--;
		}
		/* evaluate scores */
		values = (unsigned *)(bests);
		n = *values + 1;
		last = values++;
		while(--n) {
			if(Scores[*values] < template_lengths[*values]) {
				--*bests;
			} else {
				*++last = *values;
			}
			Scores[*values] = 0;
			++values;
		}
	} else {
		for(rc = 0; rc < 1 + (prefix ? 0 : 1); ++rc) {
			if(rc) { /* go over reverse complement */
				bests = bestTemplates_r;
				qseq = qseq_rc;
				seq = qseq->seq;
			}
			
			/* iterate seq */
			last = 0;
			reps = 0;
			j = 0;
			qseq->N[0]++;
			qseq->N[qseq->N[0]] = qseq->seqlen;
			for(i = 1; i <= qseq->N[0] && j < seqend; ++i) {
				/* init k-mer */
				getKmer_macro(kmer, seq, j, cPos, iPos, (shifter + 2));
				cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
				end = qseq->N[i];
				for(j += kmersize - 1; j < end; ++j) {
					/* update k-mer */
					kmer = updateKmer_macro(kmer, seq, j, mask);
					cmer = flag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
					/* lookup */
					if((values = hashMap_get(templates, cmer))) {
						if(values == last) {
							++reps;
						} else if(last) {
							if(SU) {
								values_s = (short unsigned *) last;
								n = *values_s + 1;
								while(--n) {
									if((Scores[*++values_s] += reps) == reps) {
										bests[++*bests] = *values_s;
									}
								}
							} else {
								n = *last + 1;
								while(--n) {
									if((Scores[*++last] += reps) == reps) {
										bests[++*bests] = *last;
									}
								}
							}
							reps = 1;
							last = values;
						} else {
							reps = 1;
							last = values;
						}
					}
				}
			}
			/* get last k-mer hits */
			if(last) {
				if(SU) {
					values_s = (short unsigned *) last;
					n = *values_s + 1;
					while(--n) {
						if((Scores[*++values_s] += reps) == reps) {
							bests[++*bests] = *values_s;
						}
					}
				} else {
					n = *last + 1;
					while(--n) {
						if((Scores[*++last] += reps) == reps) {
							bests[++*bests] = *last;
						}
					}
				}
			}
			
			/* evaluate scores */
			values = (unsigned *)(bests);
			n = *values + 1;
			last = values++;
			while(--n) {
				if(Scores[*values] < (template_lengths[*values] - kmersize + 1)) {
					--*bests;
				} else {
					*++last = *values;
				}
				Scores[*values] = 0;
				++values;
			}
			
			reps = 0;
			qseq->N[0]--;
		}
		
		/* merge matches */
		bests = bestTemplates_r;
		n = *bests + 1;
		while(--n) {
			bestTemplates[++*bestTemplates] = -(*++bests);
		}
	}
	
	return *bestTemplates;
}

int reassign_getoffset(CompDNA *qseq, long unsigned kmer, int offset) {
	
	unsigned cPos, iPos;
	long unsigned qmer;
	
	--offset;
	qmer = kmer + 1;
	while(qmer != kmer && ++offset < qseq->seqlen) {
		getKmer_macro(qmer, qseq->seq, offset, cPos, iPos, 0);
	}
	
	return (qmer == kmer) ? offset : -1;
}

int reassign_cmpseqs(long unsigned *qseq, long unsigned *target, int len, int offset) {
	
	int complen;
	unsigned cPos, iPos;
	long unsigned kmer;
	
	if(offset < 0) {
		return -1;
	}
	
	/* cmp in chunks of 32 */
	complen = (len >> 5) + ((len & 31) ? 1 : 0);
	while(--complen) {
		getKmer_macro(kmer, qseq, offset, cPos, iPos, 0);
		if(*target != kmer) {
			return *target < kmer ? -1 : 1;
		}
		++target;
		offset += 32;
	}
	
	/* cmp remaining */
	if(len & 31) {
		complen = 64 - ((len & 31) << 1);
		getKmer_macro(kmer, qseq, offset, cPos, iPos, complen);
		kmer <<= complen;
		if(*target != kmer) {
			return *target < kmer ? -1 : 1;
		}
	}
	
	return 0;
}

int reassign_testNs(int *N, int start, int end) {
	
	int n;
	
	n = *N + 1;
	while(--n && *++N < end) {
		if(start < *N && *N < end) {
			return *N;
		}
	}
	
	return 0;
}

int reassign_matchseqs(CompDNA *consensus, CompDNA *qseq) {
	
	/* return offset of match, -1 of no match */
	int i, start, match;
	
	/* check if N's comes in the way */
	start = 0;
	for(i = 1; i <= consensus->N[0]; ++i) {
		if(qseq->seqlen <= consensus->N[i] - start) {
			/* it is possible to place qseq in consensus */
			i = consensus->N[0];
		} else if(consensus->seqlen - consensus->N[i] < qseq->seqlen) {
			/* it is not possible to place qseq in consensus */
			return -1;
		} else {
			/* check next N */
			start = consensus->N[i] + 1;
		}
	}
	
	/* cmp sequences */
	--start;
	match = 1;
	while(match && 0 <= ++start && qseq->seqlen <= (consensus->seqlen - start)) {
		/* first match in consensus */
		if(0 <= (match = reassign_getoffset(consensus, *(qseq->seq), start))) {
			start = match;
		} else {
			start = qseq->seqlen;
		}
		
		/* extend from there */
		match = reassign_cmpseqs(consensus->seq, qseq->seq, qseq->seqlen, start);
		
		/* test for N's in match */
		if(!match && (match = reassign_testNs(consensus->N, start, start + qseq->seqlen))) {
			start = match;
		}
	}
	
	return match ? -1 : start;
}

int reassign_matchseqs_rc(CompDNA *consensus, CompDNA *consensus_rc, CompDNA *candidate) {
	
	if(reassign_matchseqs(consensus, candidate)) {
		return reassign_matchseqs(consensus_rc, candidate);
	}
	
	return 0;
}

void reassign_loadseq(CompDNA *qseq, int len, int file, long offset) {
	
	long unsigned curr;
	
	reallocComp(qseq, len);
	qseq->seqlen = len;
	qseq->complen = (len >> 5) + 1;
	qseq->N[0] = 0;
	curr = lseek(file, 0, SEEK_CUR);
	if(offset != lseek(file, offset, SEEK_SET) || read(file, qseq->seq, qseq->complen * sizeof(long unsigned)) < (qseq->complen * sizeof(long unsigned))) {
		fprintf(stderr, "Corrupted *.seq.b\n");
		if(errno) {
			ERROR();
		} else {
			exit(1);
		}
	}
	/* reset file descriptor */
	lseek(file, curr, SEEK_SET);
}

int reassign_matrix_bias(AssemInfo *matrix, int bias) {
	
	int len;
	Assembly *asm_ptr, *new_ptr;
	
	/* reallocate */
	if(matrix->size <= matrix->len + bias) {
		matrix->size = matrix->len + bias + 1;
		matrix->assmb = realloc(matrix->assmb, matrix->size * sizeof(Assembly));
	} else if(bias <= 0) {
		return 0;
	}
	
	/* bias */
	asm_ptr = matrix->assmb + matrix->len;
	new_ptr = asm_ptr + bias;
	matrix->len += bias;
	len = matrix->len + 1;
	while(--len) {
		*--new_ptr = *--asm_ptr;
		new_ptr->next += bias;
	}
	
	return bias;
}

unsigned reassign_matrix_insertions(AssemInfo *matrix) {
	
	unsigned insertions;
	
	/* reallocate */
	if(matrix->size == matrix->len) {
		matrix->size += 1024;
		matrix->assmb = realloc(matrix->assmb, matrix->size * sizeof(Assem));
		if(!matrix->assmb) {
			ERROR();
		}
	}
	
	/* get next available spot */
	insertions = matrix->len++;
	matrix->assmb[insertions].next = 0;
	return insertions;
}

void reassign_matrix_offset(AssemInfo *matrix, Assem *aligned_assem, int offset, long unsigned *tseq, int t_len) {
	
	const unsigned char bases[6] = "ACGTN-";
	int pos, bias, aln_len, asm_len, new_pos, tmp;
	unsigned insertions;
	char *s, *news;
	unsigned char *t, *q, *newt, *newq;
	short unsigned *counts;
	long unsigned depthUpdate;
	Assembly *assembly, *asm_ptr, *new_ptr;
	
	/* get to offset in matrix */
	t = aligned_assem->t;
	s = aligned_assem->s;
	q = aligned_assem->q;
	pos = 0;
	assembly = matrix->assmb;
	asm_ptr = assembly;
	bias = -offset;
	while(offset) {
		if(*q != '-') {
			--offset;
		}
		
		/* advance */
		++t;
		++s;
		++q;
		pos = asm_ptr->next;
		asm_ptr = assembly + pos;
	}
	
	/* get bias from insertions, to avoid overwriting org before new is created */
	newt = t;
	newq = q;
	aln_len = t_len;
	while(*newt && aln_len) {
		if(*newt == '-' && *newq != '-') {
			++bias;
		}
		++newt;
		++newq;
		--aln_len;
	}
	/* make space for insertions */
	bias = reassign_matrix_bias(matrix, bias);
	
	/* basecall from offset */
	insertions = 0;
	new_pos = 0;
	new_ptr = 0;
	pos += bias;
	assembly = matrix->assmb;
	asm_ptr = assembly + pos;
	newt = aligned_assem->t - 1;
	news = aligned_assem->s - 1;
	newq = aligned_assem->q - 1;
	aln_len = 0;
	asm_len = 0;
	aligned_assem->depth = 0;
	aligned_assem->depthVar = 0;
	while(aln_len != t_len) {
		if(*q != '-') {
			/* update depth */
			counts = asm_ptr->counts;
			depthUpdate = counts[0] + counts[1] + counts[2] + counts[3] + counts[4] + counts[5];
			aligned_assem->depth += depthUpdate;
			aligned_assem->depthVar += depthUpdate * depthUpdate;
			
			/* update consensus */
			*++newt = bases[getNuc(tseq, aln_len)];
			*++news = '|';
			*++newq = *q;
			
			/* update next of previous node */
			if(new_ptr) {
				new_ptr->next = aln_len;
			}
			
			/* update matrix with match */
			new_pos = aln_len;
			new_ptr = assembly + aln_len++;
			*new_ptr = *asm_ptr;
			
			/* add org node to available insertions */
			if(*t == '-') {
				asm_ptr->next = insertions;
				insertions = pos;
				asm_ptr = new_ptr;
			}
		} else {
			/* move node */
			if(pos < t_len && *t != '-') {
				/* get placement of insertion */
				if(!insertions) {
					insertions = reassign_matrix_insertions(matrix);
					assembly = matrix->assmb;
					asm_ptr = assembly + pos;
					new_ptr = assembly + new_pos;
				}
				
				/* move node to insertion pos */
				tmp = assembly[insertions].next;
				assembly[insertions] = *asm_ptr;
				pos = insertions;
				asm_ptr = assembly + pos;
				insertions = tmp;
			}
			
			/* update next of previous node */
			new_ptr->next = pos;
			
			/* update matrix with (minor) insertion */
			new_pos = pos;
			new_ptr = asm_ptr;
		}
		
		/* advance */
		++asm_len;
		++t;
		++s;
		++q;
		pos = asm_ptr->next;
		asm_ptr = assembly + pos;
	}
	
	/* terminate consensus and assembly */
	*++newt = 0;
	*++news = 0;
	*++newq = 0;
	new_ptr->next = 0;
	matrix->len = asm_len;
	
	/* update aligned_assem */
	aligned_assem->len = aln_len;
	aligned_assem->aln_len = aln_len;
	aligned_assem->cover = t_len;
}
