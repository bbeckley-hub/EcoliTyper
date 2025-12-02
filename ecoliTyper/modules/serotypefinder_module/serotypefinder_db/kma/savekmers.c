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
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include "ankers.h"
#include "compdna.h"
#include "hashmapkma.h"
#include "kmeranker.h"
#include "penalties.h"
#include "pherror.h"
#include "qseqs.h"
#include "sam.h"
#include "savekmers.h"
#include "seqmenttree.h"
#include "stdnuc.h"
#include "stdstat.h"
#include "threader.h"

void (*ankerPtr)(int*, int*, int*, char*, int*, unsigned**, unsigned**, int*, CompDNA*, int, int, int, int, Qseqs*, volatile int*, FILE*) = &ankerAndClean;
int (*kmerScan)(const HashMapKMA *, const Penalties *, int*, int*, int*, int*, CompDNA*, CompDNA*, Qseqs*, int*, const int, volatile int*, FILE*) = &save_kmers_chain;
int (*save_kmers_pair)(const HashMapKMA *, const Penalties *, int*, int*, int*, int*, int*, int*, CompDNA*, CompDNA*, const Qseqs*, const Qseqs*, int*, const int, volatile int*, FILE*) = &save_kmers_unionPair;
int (*get_kmers_for_pair_ptr)(const HashMapKMA *, const Penalties *, int *, int *, int *, int *, CompDNA *, int *, int) = &get_kmers_for_pair;
int (*getMatch)(int*, int*) = &getBestMatch;
int (*getMatchSparse)(int*, int*, int, int, int, int) = &getBestMatchSparse;
int (*getSecondForce)(int*, int*, int*, int*, int*, int*) = &getSecondBestForce;
int (*getSecondPen)(int*, int*, int*, int*, int*, int*, int, int) = &getSecondBestPen;
int (*getF)(int*, int*, int*, int*, int*) = &getF_Best;
int (*getR)(int*, int*, int*, int*, int*) = &getR_Best;

int loadFsa(CompDNA *qseq, Qseqs *header, FILE *inputfile) {
	
	int buffer[4];
	
	if(fread(buffer, sizeof(int), 4, inputfile)) {
		qseq->seqlen = buffer[0];
		qseq->complen = buffer[1];
		/* if pair, header->len < 0 */
		header->len = abs(buffer[3]);
		
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
		qseq->N[0] = buffer[2];
		
		if(header->size <= header->len) {
			header->size = header->len << 1;
			free(header->seq);
			header->seq = smalloc(header->size);
		}
		sfread(qseq->seq, sizeof(long unsigned), qseq->complen, inputfile);
		sfread(qseq->N + 1, sizeof(int), qseq->N[0], inputfile);
		sfread(header->seq, 1, header->len, inputfile);
	} else {
		qseq->seqlen = 0;
		return 0;
	}
	
	return buffer[3];
}

void * save_kmers_threaded(void *arg) {
	
	static volatile int Lock[2] = {0, 0};
	static unsigned readNum = 0;
	volatile int *excludeIn = &Lock[0], *excludeOut = &Lock[1];
	KmerScan_thread *thread = arg;
	int *Score, *Score_r, *bestTemplates, *bestTemplates_r, *regionTemplates;
	int *regionScores, *extendScore, *p_readNum, *pr_readNum, *preg_readNum;
	int go, exhaustive, unmapped, sam, flag, cflag, stats[2];;
	FILE *inputfile, *out;
	HashMapKMA *templates;
	CompDNA *qseq, *qseq_r;
	Qseqs *header, *header_r, *samseq;
	Penalties *rewards;
	
	stats[0] = 0;
	templates = thread->templates;
	exhaustive = thread->exhaustive;
	bestTemplates = thread->bestTemplates;
	bestTemplates_r = thread->bestTemplates_r;
	rewards = thread->rewards;
	if(save_kmers_pair != &save_kmers_unionPair) {
		regionScores = calloc(templates->DB_size << 1, sizeof(int));
		if(!regionScores) {
			ERROR();
		}
	} else {
		regionScores = 0;
	}
	
	/* make constant part of flag */
	if(save_kmers_pair == &save_kmers_forcePair) {
		cflag = 1;
	} else {
		cflag = 0;
	}
	if(templates->prefix && templates->prefix_len == 0) {
		cflag |= 2;
	}
	
	extendScore = calloc((templates->DB_size + 1) * sizeof(int) + templates->DB_size + 1, 1);
	if(!extendScore) {
		ERROR();
	}
	
	qseq = smalloc(sizeof(CompDNA));
	qseq_r = smalloc(sizeof(CompDNA));
	header = setQseqs(256);
	header_r = setQseqs(256);
	allocComp(qseq, 1024);
	allocComp(qseq_r, 1024);
	regionTemplates = smalloc(((templates->DB_size << 1) + 4) * sizeof(int));
	Score = calloc(templates->DB_size, sizeof(int));
	Score_r = calloc(templates->DB_size, sizeof(int));
	if(!Score || !Score_r) {
		ERROR();
	}
	inputfile = thread->inputfile;
	*Score = thread->num;
	*bestTemplates++ = templates->DB_size;
	*bestTemplates_r++ = templates->DB_size;
	*regionTemplates++ = templates->DB_size;
	*bestTemplates++ = thread->num;
	*bestTemplates_r++ = thread->num;
	*regionTemplates++ = thread->num;
	*bestTemplates++ = 0;
	*bestTemplates_r++ = 0;
	*regionTemplates++ = 0;
	out = thread->out;
	if((sam = thread->sam)) {
		samseq = setQseqs(256);
	} else {
		samseq = 0;
	}
	
	p_readNum = (bestTemplates - 1);
	pr_readNum = (bestTemplates_r - 1);
	preg_readNum = (regionTemplates - 1);
	
	go = 1;
	while(go != 0) {
		/* load qseqs */
		lock(excludeIn);
		if((go = loadFsa(qseq, header, inputfile)) < 0) {
			/* PE */
			loadFsa(qseq_r, header_r, inputfile);
		}
		if(go != 0) {
			*p_readNum = ++readNum;
		} else {
			*p_readNum = readNum;
		}
		unlock(excludeIn);
		*pr_readNum = readNum;
		*preg_readNum = readNum;
		
		/* allocate memory */
		if(qseq_r->size < qseq->size && 0 < go) {
			freeComp(qseq_r);
			allocComp(qseq_r, qseq->size);
		}
		
		/* find ankers */
		if(0 < go) {
			unmapped = kmerScan(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, qseq, qseq_r, header, extendScore, exhaustive, excludeOut, out);
		} else if(go < 0) {
			unmapped = save_kmers_pair(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores, qseq, qseq_r, header, header_r, extendScore, exhaustive, excludeOut, out);
		} else {
			unmapped = 0;
		}
		
		if(sam && unmapped) {
			if(unmapped & 1) {
				/* set flag */
				flag = 4;
				if(go < 0) {
					flag |= 65;
					if((unmapped & 2) || (cflag & 1)) {
						flag |= 8;
					}
				}
				if((cflag & 2) == 0) {
					flag |= 16;
					if((flag & 8) && (unmapped & 2)) {
						flag |= 32;
					}
				}
				qseqCompDNA(qseq, samseq);
				nibble2base(samseq->seq, samseq->len);
				stats[1] = flag;
				samwrite(samseq, header, 0, 0, 0, stats);
			}
			if(unmapped & 2 || ((cflag & 1) && (unmapped & 1))) {
				/* set flag */
				flag = 4;
				if(go < 0) {
					flag |= 129;
					if(unmapped & 1) {
						flag |= 8;
					}
				}
				if((cflag & 2) == 0) {
					if(unmapped & 2) {
						flag |= 16;
					}
					if((flag & 8)) {
						flag |= 32;
					}
				}
				qseqCompDNA(qseq_r, samseq);
				nibble2base(samseq->seq, samseq->len);
				stats[1] = flag;
				samwrite(samseq, header_r, 0, 0, 0, stats);
			}
		}
	}
	
	/* clean up */
	if(save_kmers_pair != &save_kmers_unionPair) {
		free(regionScores);
	}
	if(sam) {
		destroyQseqs(samseq);
	}
	
	free(extendScore);
	freeComp(qseq);
	free(qseq);
	freeComp(qseq_r);
	free(qseq_r);
	destroyQseqs(header);
	destroyQseqs(header_r);
	free(regionTemplates - 3);
	free(Score);
	free(Score_r);
	
	return NULL;
}

int getBestMatch(int *bestTemplates, int *Score) {
	
	int i, score, bestScore, bestHits, template;
	
	bestHits = 0;
	bestScore = 0;
	for(i = 1; i <= *bestTemplates; ++i) {
		score = Score[(template = bestTemplates[i])];
		if(score > bestScore) {
			bestScore = score;
			bestHits = 1;
			bestTemplates[bestHits] = template;
		} else if(score == bestScore) {
			++bestHits;
			bestTemplates[bestHits] = template;
		}
		Score[template] = 0;
	}
	*bestTemplates = bestHits;
	
	return bestScore;
}

int getProxiMatch(int *bestTemplates, int *Score) {
	
	static long unsigned *softProxi = 0;
	static double minFrac = 1.0;
	int i, score, bestScore, proxiScore, bestHits, template, *Templates;
	
	if(Score == 0) {
		if(minFrac < 1.0) {
			softProxi = (long unsigned *) bestTemplates;
		} else {
			minFrac = *((double *) bestTemplates);
		}
		return 0;
	}
	
	/* get proximity score */
	bestScore = 0;
	Templates = bestTemplates;
	i = *Templates++ + 1;
	while(--i) {
		if(bestScore < (score = Score[*Templates++])) {
			bestScore = score;
		}
	}
	/*proxiScore = (bestScore - proxi < 1) ? 1 : (bestScore - proxi);*/
	proxiScore = minFrac * bestScore;
	
	/* get best matches within proximity */
	bestHits = 0;
	Templates = bestTemplates;
	i = *Templates++ + 1;
	while(--i) {
		score = Score[(template = *Templates++)];
		if(proxiScore <= score) {
			bestTemplates[++bestHits] = template;
			if(softProxi) {
				softProxi[template] += score;
			}
		}
		Score[template] = 0;
	}
	*bestTemplates = bestHits;
	
	return bestScore;
}

int getBestMatchSparse(int *bestTemplates, int *Score, int kmersize, int n_kmers, int M, int MM) {
	
	int i, score, bestScore, bestHits, template;
	
	bestHits = 0;
	bestScore = 0;
	for(i = 1; i <= *bestTemplates; ++i) {
		score = Score[(template = bestTemplates[i])];
		score = score * kmersize * M + (n_kmers - score) * MM;
		if(score > bestScore) {
			bestScore = score;
			bestHits = 1;
			bestTemplates[bestHits] = template;
		} else if(score == bestScore) {
			++bestHits;
			bestTemplates[bestHits] = template;
		}
		Score[template] = 0;
	}
	*bestTemplates = bestHits;
	
	return bestScore;
}

int getProxiMatchSparse(int *bestTemplates, int *Score, int kmersize, int n_kmers, int M, int MM) {
	
	static long unsigned *softProxi = 0;
	static double minFrac = 1.0;
	int i, score, bestScore, proxiScore, bestHits, template, *Templates;
	
	if(Score == 0) {
		if(minFrac < 1.0) {
			softProxi = (long unsigned *) bestTemplates;
		} else {
			minFrac = *((double *) bestTemplates);
		}
		return 0;
	}
	
	/* get proximity score */
	bestScore = 0;
	Templates = bestTemplates;
	i = *Templates++ + 1;
	while(--i) {
		score = Score[*Templates++];
		score = score * kmersize * M + (n_kmers - score) * MM;
		if(bestScore < score) {
			bestScore = score;
		}
	}
	/*proxiScore = (bestScore - proxi < 1) ? 1 : (bestScore - proxi);*/
	proxiScore = minFrac * bestScore;
	
	/* get best matches within proximity */
	bestHits = 0;
	Templates = bestTemplates;
	i = *Templates++ + 1;
	while(--i) {
		score = Score[(template = *Templates++)];
		score = score * kmersize * M + (n_kmers - score) * MM;
		if(proxiScore <= score) {
			bestTemplates[++bestHits] = template;
			if(softProxi) {
				softProxi[template] += score;
			}
		}
		Score[template] = 0;
	}
	*bestTemplates = bestHits;
	
	return bestScore;
}

static int clearScore(int *bestTemplates, int *Score) {
	
	int i;
	
	i = *bestTemplates + 1;
	while(--i) {
		Score[*++bestTemplates] = 0;
	}
	
	return 0;
}

int get_kmers_for_pair(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, int *extendScore, const int exhaustive) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	   and is the time determining step */
	int i, j, l, rc, end, HIT, gaps, score, Ms, MMs, Us, W1s, template, flag;
	int hitCounter, bestSeqCount, kmersize, mlen, SU, shifter, W1, U, M, MM;
	int n, cPos, iPos, mPos, hLen, seqend, j_u, m, mm, *bests, *Scores;
	unsigned *values, *last;
	short unsigned *values_s;
	long unsigned mask, mmask, kmer, cmer, hmer, *seq;
	char *include;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return 0;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = templates->flag;
	W1 = rewards->W1;
	U = rewards->U;
	M = rewards->M;
	MM = rewards->MM;
	*extendScore = 0;
	*bestTemplates = 0;
	*bestTemplates_r = 0;
	bests = bestTemplates;
	Scores = Score;
	bestSeqCount = 0;
	kmersize = templates->kmersize;
	include = (char *) (extendScore + (templates->DB_size + 1));
	seq = qseq->seq;
	seqend = qseq->seqlen - kmersize + 1;
	hLen = kmersize;
	
	for(rc = 0; rc < 2; ++rc) {
		if(rc) {
			bests = bestTemplates_r;
			Scores = Score_r;
			comp_rc(qseq);
		}
		/* Make quick check of the qseq */
		HIT = exhaustive;
		hitCounter = 0;
		j = 0;
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(; j < end && !HIT; j += kmersize) {
				getKmer_macro(kmer, seq, j, cPos, iPos, shifter);
				cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
				if(hashMap_get(templates, cmer)) {
					HIT = 1;
				}
			}
			j = qseq->N[i] + 1;
		}
		
		/* If deltamer qseq hits, then continue */
		if(HIT) {
			/* Scan the deltamer exhaustively, and collect scores in Score*/
			last = 0;
			gaps = 0;
			HIT = 0;
			Ms = 0;
			MMs = 0;
			Us = 0;
			W1s = 0;
			j = 0;
			for(i = 1; i <= qseq->N[0] && j < seqend; ++i) {
				/* init k-mer */
				getKmer_macro(kmer, seq, j, cPos, iPos, (shifter + 2));
				cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
				end = qseq->N[i];
				for(j_u = j + kmersize - 1; j_u < end; ++j_u) {
					/* update k-mer */
					kmer = updateKmer_macro(kmer, seq, j_u, mask);
					cmer = flag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
					/* lookup */
					if((values = hashMap_get(templates, cmer))) {
						if(values == last) {
							/*
							gaps == 0 -> Match
							gaps == kmersize -> 1 MM
							kmersize < gaps -> several mismatches or indel(s)
							gaps < kmersize -> deletion
							*/
							if(gaps == 0) {
								/* match */
								++Ms;
							} else if(mlen <= gaps && gaps <= kmersize) {
								/* snp */
								Ms += kmersize;
								++MMs;
							} else if(kmersize < gaps) {
								/* mismatch or insersion */
								Ms += kmersize;
								gaps -= (kmersize - 1); /* adjust for consecutive k-mer mismatches */
								if(gaps <= 2) {
									mm = gaps;
									m = 0;
								} else {
									mm = gaps / kmersize + (gaps % kmersize ? 1 : 0);
									mm = MAX(2, mm);
									m = MIN(gaps - mm, kmersize);
									m = MIN(m, mm);
								}
								
								/* evaluate best option */
								if((W1 + (gaps - 1) * U) <= (mm * MM + m * M)) {
									MMs += mm;
									Ms += m;
								} else {
									++W1s;
									Us += (gaps -1);
								}
							} else if(mlen != kmersize) {
								Ms += gaps;
								++MMs;
								/* unlikely deletion or random k-mer mismatch, 
								unless minimizers or homopolymer compression,
								assume better and go random zero score else. */
							} else {
								/* indel */
								Ms += gaps;
								++W1s;
								Us += (kmersize - gaps);
							}
							
							HIT = j;
							gaps = 0;
						} else {
							if(last) {
								score = Ms * M + MMs * MM + Us * U + W1s * W1;
								values_s = (short unsigned *) last;
								l = SU ? (*values_s + 1) : (*last + 1);
								while(--l) {
									template = SU ? *++values_s : *++last;
									Scores[template] += score;
									extendScore[template] = HIT;
								}
								HIT = j - 1;
								last = values;
								values_s = (short unsigned *) values;
								n = SU ? *values_s : *values;
								l = n + 1;
								while(--l) {
									template = SU ? *++values_s : *++values;
									if(include[template]) {
										gaps = HIT - extendScore[template];
										if(gaps == 0) {
											/* match */
											score = M;
										} else if(mlen <= gaps && gaps <= kmersize) {
											/* snp */
											score = kmersize * M + MM;
										} else if(kmersize < gaps) {
											/* mismatch or insersion */
											gaps -= (kmersize - 1); /* adjust for consecutive k-mer mismatches */
											if(gaps <= 2) {
												mm = gaps;
												m = 0;
											} else {
												mm = gaps / kmersize + (gaps % kmersize ? 1 : 0);
												mm = MAX(2, mm);
												m = MIN(gaps - mm, kmersize);
												m = MIN(m, mm);
											}
											
											/* evaluate best option */
											if((W1 + (gaps - 1) * U) <= (mm * MM + m * M)) {
												score = kmersize * M + (mm * MM + m * M);
											} else {
												score = kmersize * M + (W1 + (gaps - 1) * U);
											}
										} else if(mlen != kmersize) {
											score = gaps * M + MM;
											/* unlikely deletion or random k-mer mismatch, 
											unless minimizers or homopolymer compression,
											assume better and go random zero score else. */
										} else {
											/* indel */
											score = gaps * M + (kmersize - gaps) * U + W1;
										}
										Scores[template] += score;
									} else {
										Scores[template] = kmersize * M;
										if(include[template] == 0) {
											include[template] = 1;
											bests[++*bests] = template;
										}
									}
								}
							} else {
								last = values;
								values_s = (short unsigned *) values;
								n = SU ? *values_s : *values;
								score = kmersize * M;
								for(l = 1; l <= n; ++l) {
									template = SU ? *++values_s : *++values;
									Scores[template] = score;
									include[template] = 1;
									bests[l] = template;
								}
								*bests = n;
							}
							
							HIT = j;
							gaps = 0;
							Ms = 0;
							MMs = 0;
							Us = 0;
							W1s = 0;
						}
						++hitCounter;
					} else {
						++gaps;
					}
					++j;
				}
				gaps += (qseq->N[i] + 1 - j); /* gap over N's */
				j = end + 1;
			}
			
			if(last) {
				score = Ms * M + MMs * MM + Us * U + W1s * W1;
				if(SU) {
					values_s = (short unsigned *) last;
					l = (*values_s) + 1;
					while(--l) {
						Scores[values_s[l]] += score;
					}
				} else {
					l = (*last) + 1;
					while(--l) {
						Scores[last[l]] += score;
					}
				}
				for(l = *bests; l != 0; --l) {
					extendScore[(template = bests[l])] = 0;
					include[template] = 0;
					if(Scores[template] < 0) {
						Scores[template] = 0;
					}
				}
			}
			
			if(bestSeqCount < hitCounter) {
				bestSeqCount = hitCounter;
			}
		}
		qseq->N[0]--;
	}
	
	return bestSeqCount;
}

int get_kmers_for_pair_count(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, int *extendScore, const int exhaustive) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	   and is the time determining step */
	int i, j, rc, end, HIT, hitCounter, bestSeqCount, reps, SU, kmersize;
	int shifter, mlen, flag, n, cPos, iPos, mPos, hLen, seqend;
	int *bests, *Scores;
	unsigned *values, *last;
	short unsigned *values_s;
	long unsigned mask, mmask, kmer, cmer, hmer, *seq;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return 0;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	*extendScore = 0;
	*bestTemplates = 0;
	*bestTemplates_r = 0;
	bests = bestTemplates;
	Scores = Score;
	bestSeqCount = 0;
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = templates->flag;
	seq = qseq->seq;
	seqend = qseq->seqlen - kmersize + 1;
	hLen = kmersize;
	
	for(rc = 0; rc < 2; ++rc) {
		if(rc) {
			bests = bestTemplates_r;
			Scores = Score_r;
			comp_rc(qseq);
		}
		/* Make quick check of the qseq */
		HIT = exhaustive;
		hitCounter = 0;
		j = 0;
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(; j < end && !HIT; j += kmersize) {
				getKmer_macro(kmer, seq, j, cPos, iPos, shifter);
				cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
				if(hashMap_get(templates, cmer)) {
					HIT = 1;
				}
			}
			j = qseq->N[i] + 1;
		}
		
		/* If deltamer qseq hits, then continue */
		if(HIT) {
			/* Scan the deltamer exhaustively, and collect scores in Score*/
			last = 0;
			reps = 0;
			j = 0;
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
						} else {
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
								hitCounter += reps;
							}
							reps = 1;
							last = values;
						}
					}
				}
				j = qseq->N[i] + 1;
			}
			
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
				hitCounter += reps;
			}
			
			reps = 0;
			if(bestSeqCount < hitCounter) {
				bestSeqCount = hitCounter;
			}
		}
		qseq->N[0]--;
	}
	
	return bestSeqCount;
}

int get_kmers_for_pair_Sparse(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, int *extendScore, const int exhaustive) {
	
	int i, j, n, end, rc, prefix_len, hitCounter, reps, n_kmers, kmersize;
	int HIT, SU, mlen, flag, cPos, iPos, mPos, hLen, seqend, *bests, *Scores;
	unsigned shifter, prefix_shifter, *values, *last;
	short unsigned *values_s;
	long unsigned pmask, mask, mmask, prefix, pmer, kmer, cmer, hmer, *seq;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return 0;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
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
	
	if(*extendScore) {
		bests = bestTemplates_r;
		Scores = Score_r;
	} else {
		bests = bestTemplates;
		Scores = Score;
	}
	
	*extendScore = 0;
	*bestTemplates = 0;
	*bestTemplates_r = 0;
	hitCounter = 0;
	n_kmers = 0;
	end = qseq->seqlen;
	if(prefix_len) {
		for(rc = 0; rc < 2; ++rc) {
			if(rc) {
				comp_rc(qseq);
			}
			j = 0;
			qseq->N[0]++;
			qseq->N[qseq->N[0]] = qseq->seqlen;
			for(i = 1; i <= qseq->N[0] && j < seqend; ++i) {
				getKmer_macro(pmer, seq, j, cPos, iPos, (prefix_shifter + 2));
				end = qseq->N[i] - kmersize;
				for(j += prefix_len - 1; j < end; ++j) {
					pmer = updateKmer_macro(pmer, seq, j, pmask);
					if(pmer == prefix) {
						/* get kmer */
						getKmer_macro(kmer, seq, (j + 1), cPos, iPos, shifter);
						cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
						if((values = hashMap_get(templates, cmer))) {
							if(SU) {
								values_s = (short unsigned *) values;
								n = *values_s + 1;
								while(--n) {
									if(Scores[*++values_s]++ == 0) {
										bests[++*bests] = *values_s;
									}
								}
							} else {
								n = *values + 1;
								while(--n) {
									if(Scores[*++values]++ == 0) {
										bests[++*bests] = *values;
									}
								}
							}
							++hitCounter;
						}
						++n_kmers;
					}
				}
				j = qseq->N[i] + 1;
			}
			qseq->N[0]--;
		}
		if(hitCounter) {
			hitCounter *= (((qseq->seqlen - kmersize + 1) << 1) / n_kmers);
		}
	} else {
		HIT = exhaustive;
		j = 0;
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(; j < end && !HIT; j += kmersize) {
				getKmer_macro(kmer, seq, j, cPos, iPos, shifter);
				cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
				if(hashMap_get(templates, cmer)) {
					HIT = 1;
				}
			}
			j = qseq->N[i] + 1;
		}
		
		if(HIT) {
			last = 0;
			reps = 0;
			j = 0;
			for(i = 1; i <= qseq->N[0] && j < seqend; ++i) {
				getKmer_macro(kmer, seq, j, cPos, iPos, (shifter + 2));
				cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
				end = qseq->N[i];
				for(j += kmersize - 1; j < end; ++j) {
					/* update kmer */
					kmer = updateKmer_macro(kmer, seq, j, mask);
					cmer = flag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
					if((values = hashMap_get(templates, cmer))) {
						if(values == last) {
							++reps;
						} else {
							if(last) {
								if(SU) {
									values_s = (short unsigned *) values;
									n = *values_s + 1;
									while(--n) {
										if((Scores[*++values_s] += reps) == reps) {
											bests[++*bests] = *values_s;
										}
									}
								} else {
									n = *values + 1;
									while(--n) {
										if((Scores[*++last] += reps) == reps) {
											bests[++*bests] = *last;
										}
									}
								}
								hitCounter += reps;
							}
							reps = 1;
							last = values;
						}
					}
				}
				j = end + 1;
			}
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
				hitCounter += reps;
			}
		}
		qseq->N[0]--;
	}
	
	return hitCounter;
}

int get_kmers_for_pair_pseoudoSparse(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, int *extendScore, const int exhaustive) {
	
	int i, j, l, n, end, template, hitCounter, gaps, Ms, MMs, Us, W1s, flag;
	int W1, U, M, MM, HIT, SU, kmersize, score, mlen, mPos, hLen, cPos, iPos;
	int seqend, j_u, m, mm, *bests, *Scores;
	unsigned shifter, *values, *last;
	short unsigned *values_s;
	char *include;
	long unsigned mask, mmask, kmer, cmer, hmer, *seq;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return 0;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	seq = qseq->seq;
	kmersize = templates->kmersize;
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = templates->flag;
	W1 = rewards->W1;
	U = rewards->U;
	M = rewards->M;
	MM = rewards->MM;
	include = (char *) (extendScore + (templates->DB_size + 1));
	hLen = kmersize;
	
	if(*extendScore) {
		bests = bestTemplates_r;
		Scores = Score_r;
	} else {
		bests = bestTemplates;
		Scores = Score;
	}
	
	*extendScore = 0;
	*bestTemplates = 0;
	*bestTemplates_r = 0;
	hitCounter = 0;
	end = qseq->seqlen;
	seqend = qseq->seqlen - kmersize + 1;
	
	HIT = exhaustive;
	j = 0;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
		end = qseq->N[i] - kmersize + 1;
		for(; j < end && !HIT; j += kmersize) {
			getKmer_macro(kmer, seq, j, cPos, iPos, shifter);
			cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
			if(hashMap_get(templates, cmer)) {
				HIT = 1;
			}
		}
		j = qseq->N[i] + 1;
	}
	
	if(HIT) {
		last = 0;
		gaps = 0;
		HIT = 0;
		Ms = 0;
		MMs = 0;
		Us = 0;
		W1s = 0;
		j = 0;
		for(i = 1; i <= qseq->N[0] && j < seqend; ++i) {
			/* init k-mer */
			getKmer_macro(kmer, seq, j, cPos, iPos, (shifter + 2));
			cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
			end = qseq->N[i];
			for(j_u = j + kmersize - 1; j_u < end; ++j_u) {
				/* update k-mer */
				kmer = updateKmer_macro(kmer, seq, j_u, mask);
				cmer = flag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
				/* lookup */
				if((values = hashMap_get(templates, cmer))) {
					if(values == last) {
						/*
						gaps == 0 -> Match
						gaps == kmersize -> 1 MM
						kmersize < gaps -> several mismatches or indel(s)
						gaps < kmersize -> deletion
						*/
						if(gaps == 0) {
							/* match */
							++Ms;
						} else if(gaps == kmersize) {
							/* snp */
							Ms += kmersize;
							++MMs;
						} else if(kmersize < gaps) {
							/* mismatch or insersion */
							Ms += kmersize;
							gaps -= (kmersize - 1); /* adjust for consecutive k-mer mismatches */
							if(gaps <= 2) {
								mm = gaps;
								m = 0;
							} else {
								mm = gaps / kmersize + (gaps % kmersize ? 1 : 0);
								mm = MAX(2, mm);
								m = MIN(gaps - mm, kmersize);
								m = MIN(m, mm);
							}
							
							/* evaluate best option */
							if((W1 + (gaps - 1) * U) <= (mm * MM + m * M)) {
								MMs += mm;
								Ms += m;
							} else {
								++W1s;
								Us += (gaps -1);
							}
						} else if(mlen != kmersize) {
							Ms += gaps;
							++MMs;
							/* unlikely deletion or random k-mer mismatch, 
							unless minimizers or homopolymer compression,
							assume better and go random zero score else. */
						} else {
							/* indel */
							Ms += gaps;
							++W1s;
							Us += (kmersize - gaps);
						}
						
						HIT = j;
						gaps = 0;
					} else {
						if(last) {
							score = Ms * M + MMs * MM + Us * U + W1s * W1;
							values_s = (short unsigned *) last;
							l = SU ? (*values_s + 1) : (*last + 1);
							while(--l) {
								template = SU ? *++values_s : *++last;
								Scores[template] += score;
								extendScore[template] = HIT;
							}
							HIT = j - 1;
							last = values;
							values_s = (short unsigned *) values;
							n = SU ? *values_s : *values;
							l = n + 1;
							while(--l) {
								template = SU ? *++values_s : *++values;
								if(include[template]) {
									gaps = HIT - extendScore[template];
									if(gaps == 0) {
										/* match */
										score = M;
									} else if(gaps == kmersize) {
										/* snp */
										score = kmersize * M + MM;
									} else if(kmersize < gaps) {
										/* mismatch or insersion */
										gaps -= (kmersize - 1); /* adjust for consecutive k-mer mismatches */
										if(gaps <= 2) {
											mm = gaps;
											m = 0;
										} else {
											mm = gaps / kmersize + (gaps % kmersize ? 1 : 0);
											mm = MAX(2, mm);
											m = MIN(gaps - mm, kmersize);
											m = MIN(m, mm);
										}
										
										/* evaluate best option */
										if((W1 + (gaps - 1) * U) <= (mm * MM + m * M)) {
											score = kmersize * M + (mm * MM + m * M);
										} else {
											score = kmersize * M + (W1 + (gaps - 1) * U);
										}
									} else if(mlen != kmersize) {
										score = gaps * M + MM;
										/* unlikely deletion or random k-mer mismatch, 
										unless minimizers or homopolymer compression,
										assume better and go random zero score else. */
									} else {
										/* indel */
										score = gaps * M + (kmersize - gaps) * U + W1;
									}
									Scores[template] += score;
								} else {
									Scores[template] = kmersize * M;
									if(include[template] == 0) {
										include[template] = 1;
										bests[++*bests] = template;
									}
								}
							}
						} else {
							last = values;
							values_s = (short unsigned *) values;
							n = SU ? *values_s : *values;
							Ms = kmersize * M;
							for(l = 1; l <= n; ++l) {
								template = SU ? *++values_s : *++values;
								Scores[template] = Ms;
								include[template] = 1;
								bests[l] = template;
							}
							*bests = n;
						}
						
						HIT = j;
						gaps = 0;
						Ms = 0;
						MMs = 0;
						Us = 0;
						W1s = 0;
					}
					++hitCounter;
				} else {
					++gaps;
				}
				++j;
			}
			gaps += (qseq->N[i] + 1 - j); /* gap over N's */
			j = end + 1;
		}
		
		if(last) {
			score = Ms * M + MMs * MM + Us * U + W1s * W1;
			if(SU) {
				values_s = (short unsigned *) last;
				l = *values_s + 1;
				while(--l) {
					Scores[*++values_s] += score;
				}
			} else {
				l = *last + 1;
				while(--l) {
					Scores[*++last] += score;
				}
			}
			for(l = *bests; l != 0; --l) {
				extendScore[(template = bests[l])] = 0;
				include[template] = 0;
				if(Scores[template] < 0) {
					Scores[template] = 0;
				}
			}
		}
	}
	qseq->N[0]--;
	
	return hitCounter;
}

void getFirstForce(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores) {
	
	int i, bestHits;
	
	bestHits = 0;
	
	for(i = 1; i <= *bestTemplates; ++i) {
		++bestHits;
		regionTemplates[bestHits] = bestTemplates[i];
		regionScores[bestHits] = Score[bestTemplates[i]];
		Score[bestTemplates[i]] = 0;
	}
	for(i = 1; i <= *bestTemplates_r; ++i) {
		++bestHits;
		regionTemplates[bestHits] = -bestTemplates_r[i];
		regionScores[bestHits] = Score_r[bestTemplates_r[i]];
		Score_r[bestTemplates_r[i]] = 0;
	}
	*regionTemplates = bestHits;
}

int getSecondBestForce(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores) {
	
	int i, score, bestHits, bestScore;
	
	bestHits = 0;
	bestScore = 0;
	for(i = 1; i <= *regionTemplates; ++i) {
		if(0 < regionTemplates[i]) {
			if((score = Score[regionTemplates[i]])) {
				score += regionScores[i];
				if(bestScore < score) {
					bestScore = score;
					bestHits = 1;
					regionTemplates[bestHits] = regionTemplates[i];
				} else if(bestScore == score) {
					++bestHits;
					regionTemplates[bestHits] = regionTemplates[i];
				}
			}
		} else {
			if((score = Score_r[-regionTemplates[i]])) {
				score += regionScores[i];
				if(bestScore < score) {
					bestScore = score;
					bestHits = 1;
					regionTemplates[bestHits] = regionTemplates[i];
				} else if(bestScore == score) {
					++bestHits;
					regionTemplates[bestHits] = regionTemplates[i];
				}
			}
		}
	}
	*regionTemplates = bestHits;
	
	for(i = *bestTemplates; i != 0; --i) {
		Score[bestTemplates[i]] = 0;
	}
	for(i = *bestTemplates_r; i != 0; --i) {
		Score_r[bestTemplates_r[i]] = 0;
	}
	
	return bestScore;
}

int getSecondProxiForce(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores) {
	
	static long unsigned *softProxi = 0;
	static double minFrac = 0.0;
	int i, score, bestHits, bestScore, proxiScore, template, *Templates;
	
	if(Score == 0) {
		minFrac = *((double *) bestTemplates);
		softProxi = (long unsigned *) bestTemplates_r;
		return 0;
	}
	
	/* get proximity score */
	bestScore = 0;
	Templates = regionTemplates;
	i = *Templates++ + 1;
	while(--i) {
		template = *Templates++;
		if(template < 0) {
			if(bestScore < (score = Score_r[-template])) {
				bestScore = score;
			}
		} else {
			if(bestScore < (score = Score[template])) {
				bestScore = score;
			}
		}
	}
	proxiScore = minFrac * bestScore;
	
	/* get best matches within proximity */
	bestHits = 0;
	Templates = regionTemplates;
	i = *Templates++ + 1;
	while(--i) {
		template = *Templates++;
		if(template < 0) {
			if(proxiScore <= (score = Score_r[-template])) {
				regionTemplates[++bestHits] = template;
				if(softProxi) {
					softProxi[-template] += score;
				}
			}
		} else {
			if(proxiScore <= (score = Score[template])) {
				regionTemplates[++bestHits] = template;
				if(softProxi) {
					softProxi[template] += score;
				}
			}
		}
	}
	
	for(i = *bestTemplates; i != 0; --i) {
		Score[bestTemplates[i]] = 0;
	}
	for(i = *bestTemplates_r; i != 0; --i) {
		Score_r[bestTemplates_r[i]] = 0;
	}
	
	return bestScore;
}

int getFirstPen(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores) {
	
	int i, score, bestScore, bestHits;
	
	bestScore = 0;
	bestHits = 0;
	/* save forward matches */
	for(i = 1; i <= *bestTemplates; ++i) {
		if(bestScore < (score = Score[bestTemplates[i]])) {
			bestScore = score;
		}
		++bestHits;
		regionTemplates[bestHits] = bestTemplates[i];
		regionScores[bestHits] = score;
		Score[bestTemplates[i]] = 0;
	}
	
	/* save reverse matches */
	for(i = 1; i <= *bestTemplates_r; ++i) {
		if(bestScore < (score = Score_r[bestTemplates_r[i]])) {
			bestScore = score;
		}
		++bestHits;
		regionTemplates[bestHits] = -bestTemplates_r[i];
		regionScores[bestHits] = score;
		Score_r[bestTemplates_r[i]] = 0;
	}
	*regionTemplates = bestHits;
	
	return bestScore;
}

int getSecondBestPen(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores, int bestScore, int PE) {
	
	int i, score, bestScore_r, compScore, bestHits, template;
	
	/* get best scoring tempates */
	bestScore_r = 0;
	for(i = 1; i <= *bestTemplates; ++i) {
		if(bestScore_r < Score[bestTemplates[i]]) {
			bestScore_r = Score[bestTemplates[i]];
		}
	}
	bestHits = *bestTemplates;
	for(i = 1; i <= *bestTemplates_r; ++i) {
		if(bestScore_r < Score_r[bestTemplates_r[i]]) {
			bestScore_r = Score_r[bestTemplates_r[i]];
		}
		++bestHits;
		bestTemplates[bestHits] = -bestTemplates_r[i];
	}
	*bestTemplates = bestHits;
	
	/* check union */
	bestHits = 0;
	if(bestScore_r) {
		compScore =  bestScore + bestScore_r - PE;
		compScore = MAX(0, compScore);
		for(i = 1; i <= *regionTemplates; ++i) {
			if(0 < regionTemplates[i]) {
				/* we got one */
				if(0 < (score = Score_r[regionTemplates[i]])) {
					score += regionScores[i];
					if(compScore < score) {
						compScore = score;
						bestHits = 1;
						regionTemplates[bestHits] = regionTemplates[i];
					} else if(compScore == score) {
						++bestHits;
						regionTemplates[bestHits] = regionTemplates[i];
					}
				}
			} else {
				/* we got one */
				if(0 < (score = Score[-regionTemplates[i]])) {
					score += regionScores[i];
					if(compScore < score) {
						compScore = score;
						bestHits = 1;
						regionTemplates[bestHits] = regionTemplates[i];
					} else if(compScore == score) {
						++bestHits;
						regionTemplates[bestHits] = regionTemplates[i];
					}
				}
			}
		}
	}
	
	/* mark as PE */
	if(bestHits) {
		*regionTemplates = -bestHits;
		/* clear scores */
		for(i = *bestTemplates; i != 0; --i) {
			if(0 < bestTemplates[i]) {
				Score[bestTemplates[i]] = 0;
			} else {
				Score_r[-bestTemplates[i]] = 0;
			}
		}
	} else {
		/* get bestHits from each as SE */
		for(i = 1; i <= *regionTemplates; ++i) {
			if(bestScore == regionScores[i]) {
				regionTemplates[++bestHits] = regionTemplates[i];
			}
		}
		*regionTemplates = bestHits;
		
		bestHits = 0;
		for(i = 1; i <= *bestTemplates; ++i) {
			if(0 < (template = bestTemplates[i])) {
				if(bestScore_r == Score[template]) {
					++bestHits;
					bestTemplates[bestHits] = template;
				}
				Score[template] = 0;
			} else {
				if(bestScore_r <= Score_r[-template]) {
					++bestHits;
					bestTemplates[bestHits] = template;
				}
				Score_r[-template] = 0;
			}
		}
		*bestTemplates = bestHits;
	}
	
	return bestScore_r;
}

int getSecondProxiPen(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores, int bestScore, int PE) {
	
	static long unsigned *softProxi = 0;
	static double minFrac = 0.0;
	int i, score, bestScore_r, compScore, proxiScore, bestHits, template;
	
	if(Score == 0) {
		minFrac = *((double *) bestTemplates);
		softProxi = (long unsigned *) bestTemplates_r;
		return 0;
	}
	
	/* get best scoring tempates */
	bestScore_r = 0;
	for(i = 1; i <= *bestTemplates; ++i) {
		if(bestScore_r < Score[bestTemplates[i]]) {
			bestScore_r = Score[bestTemplates[i]];
		}
	}
	bestHits = *bestTemplates;
	for(i = 1; i <= *bestTemplates_r; ++i) {
		if(bestScore_r < Score_r[bestTemplates_r[i]]) {
			bestScore_r = Score_r[bestTemplates_r[i]];
		}
		++bestHits;
		bestTemplates[bestHits] = -bestTemplates_r[i];
	}
	*bestTemplates = bestHits;
	
	/* check union */
	bestHits = 0;
	if(bestScore_r) {
		compScore = 0;
		for(i = 1; i <= *regionTemplates; ++i) {
			if(0 < (template = regionTemplates[i])) {
				/* we got one */
				if(0 < (score = Score_r[template])) {
					score += regionScores[i];
					if(compScore < score) {
						compScore = score;
					}
				}
			} else {
				/* we got one */
				if(0 < (score = Score[-template])) {
					score += regionScores[i];
					if(compScore < score) {
						compScore = score;
					}
				}
			}
		}
		
		/* union is better */
		if((bestScore + bestScore_r - PE) <= compScore) {
			proxiScore = minFrac * compScore;
			for(i = 1; i <= *regionTemplates; ++i) {
				if(0 < (template = regionTemplates[i])) {
					/* we got one */
					if(0 < (score = Score_r[template])) {
						score += regionScores[i];
						if(proxiScore <= score) {
							regionTemplates[++bestHits] = template;
							if(softProxi) {
								softProxi[template] += score;
							}
						}
					}
				} else {
					/* we got one */
					if(0 < (score = Score[-template])) {
						score += regionScores[i];
						if(proxiScore <= score) {
							regionTemplates[++bestHits] = template;
							if(softProxi) {
								softProxi[-template] += score;
							}
						}
					}
				}
			}
		}
	}
	
	/* mark as PE */
	if(bestHits) {
		*regionTemplates = -bestHits;
		/* clear scores */
		for(i = *bestTemplates; i != 0; --i) {
			if(0 < bestTemplates[i]) {
				Score[bestTemplates[i]] = 0;
			} else {
				Score_r[-bestTemplates[i]] = 0;
			}
		}
	} else {
		/* get bestHits from each as SE */
		proxiScore = minFrac * bestScore;
		for(i = 1; i <= *regionTemplates; ++i) {
			if(proxiScore <= regionScores[i]) {
				regionTemplates[++bestHits] = regionTemplates[i];
			}
		}
		*regionTemplates = bestHits;
		
		bestHits = 0;
		proxiScore = minFrac * bestScore_r;
		for(i = 1; i <= *bestTemplates; ++i) {
			if(0 < (template = bestTemplates[i])) {
				if(proxiScore <= (score = Score[template])) {
					++bestHits;
					bestTemplates[bestHits] = template;
					if(softProxi) {
						softProxi[template] += score;
					}
				}
				Score[template] = 0;
			} else {
				if(proxiScore <= (score = Score_r[-template])) {
					++bestHits;
					bestTemplates[bestHits] = template;
					if(softProxi) {
						softProxi[-template] += score;
					}
				}
				Score_r[-template] = 0;
			}
		}
		*bestTemplates = bestHits;
	}
	
	return bestScore_r;
}

int getF_Best(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates) {
	
	int i, score, bestScore, bestHits;
	
	bestScore = 0;
	bestHits = 0;
	
	for(i = 1; i <= *bestTemplates; ++i) {
		if(bestScore < (score = Score[bestTemplates[i]])) {
			bestScore = score;
			bestHits = 1;
			regionTemplates[bestHits] = bestTemplates[i];
		} else if(bestScore == score) {
			++bestHits;
			regionTemplates[bestHits] = bestTemplates[i];
		}
		Score[bestTemplates[i]] = 0;
	}
	for(i = 1; i <= *bestTemplates_r; ++i) {
		if(bestScore < (score = Score_r[bestTemplates_r[i]])) {
			bestScore = score;
			bestHits = 1;
			regionTemplates[bestHits] = -bestTemplates_r[i];
		} else if(bestScore == score) {
			++bestHits;
			regionTemplates[bestHits] = -bestTemplates_r[i];
		}
		Score_r[bestTemplates_r[i]] = 0;
	}
	*regionTemplates = bestHits;
	
	return bestScore;
}

int getR_Best(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates) {
	
	int i, j, score, bestScore_r, bestHits;
	
	/* get best scoring tempates */
	bestScore_r = 0;
	bestHits = 0;
	for(i = 1; i <= *bestTemplates; ++i) {
		if(bestScore_r < (score = Score[bestTemplates[i]])) {
			for(j = bestHits; j != 0; --j) {
				Score[bestTemplates[j]] = 0;
			}
			bestScore_r = score;
			bestHits = 1;
			bestTemplates[bestHits] = bestTemplates[i];
		} else if(bestScore_r == score) {
			++bestHits;
			bestTemplates[bestHits] = bestTemplates[i];
		} else {
			Score[bestTemplates[i]] = 0;
		}
	}
	for(i = 1; i <= *bestTemplates_r; ++i) {
		if(bestScore_r < (score = Score_r[bestTemplates_r[i]])) {
			for(j = bestHits; j != 0; --j) {
				if(0 < bestTemplates[j]) {
					Score[bestTemplates[j]] = 0;
				} else {
					Score_r[-bestTemplates[j]] = 0;
				}
			}
			bestScore_r = score;
			bestHits = 1;
			bestTemplates[bestHits] = -bestTemplates_r[i];
		} else if(bestScore_r == score) {
			++bestHits;
			bestTemplates[bestHits] = -bestTemplates_r[i];
		} else {
			Score_r[bestTemplates_r[i]] = 0;
		}
	}
	*bestTemplates = bestHits;
	
	/* check union */
	bestHits = 0;
	for(i = 1; i <= *regionTemplates; ++i) {
		if(0 < regionTemplates[i]) {
			/* we got one */
			if(Score_r[regionTemplates[i]]) {
				++bestHits;
				score = regionTemplates[bestHits];
				regionTemplates[bestHits] = regionTemplates[i];
				regionTemplates[i] = score;
			}
		} else {
			/* we got one */
			if(Score[-regionTemplates[i]]) {
				++bestHits;
				score = regionTemplates[bestHits];
				regionTemplates[bestHits] = regionTemplates[i];
				regionTemplates[i] = score;
			}
		}
	}
	
	/* mark as PE */
	if(bestHits) {
		*regionTemplates = -bestHits;
	}
	
	/* clear scores */
	for(i = *bestTemplates; i != 0; --i) {
		if(0 < bestTemplates[i]) {
			Score[bestTemplates[i]] = 0;
		} else {
			Score_r[-bestTemplates[i]] = 0;
		}
	}
	
	return bestScore_r;
}

int getF_Proxi(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates) {
	
	static long unsigned *softProxi = 0;
	static double minFrac = 0.0;
	int i, score, bestScore, proxiScore, bestHits, template, *Templates;
	
	if(Score == 0) {
		minFrac = *((double *) bestTemplates);
		softProxi = (long unsigned *) bestTemplates_r;
		return getR_Proxi(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates);
	}
	
	/* get proximity score */
	bestScore = 0;
	Templates = bestTemplates;
	i = *Templates++ + 1;
	while(--i) {
		if(bestScore < (score = Score[*Templates++])) {
			bestScore = score;
		}
	}
	Templates = bestTemplates_r;
	i = *Templates++ + 1;
	while(--i) {
		if(bestScore < (score = Score_r[*Templates++])) {
			bestScore = score;
		}
	}
	proxiScore = minFrac * bestScore;
	
	/* get best matches within proximity */
	bestHits = 0;
	Templates = bestTemplates;
	i = *Templates++ + 1;
	while(--i) {
		score = Score[(template = *Templates++)];
		if(proxiScore <= score) {
			regionTemplates[++bestHits] = template;
			if(softProxi) {
				softProxi[template] += score;
			}
		}
		Score[template] = 0;
	}
	Templates = bestTemplates_r;
	i = *Templates++ + 1;
	while(--i) {
		score = Score_r[(template = *Templates++)];
		if(proxiScore <= score) {
			regionTemplates[++bestHits] = -template;
			if(softProxi) {
				softProxi[template] += score;
			}
		}
		Score_r[template] = 0;
	}
	*regionTemplates = bestHits;
	
	return bestScore;
}

int getR_Proxi(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates) {
	
	static long unsigned *softProxi = 0;
	static double minFrac = 0.0;
	int i, score, proxiScore, bestScore, bestHits, template, *Templates;
	
	if(Score == 0) {
		minFrac = *((double *) bestTemplates);
		softProxi = (long unsigned *) bestTemplates_r;
		return 0;
	}
	
	/* get proximity score */
	bestScore = 0;
	Templates = bestTemplates;
	i = *Templates++ + 1;
	while(--i) {
		if(bestScore < (score = Score[*Templates++])) {
			bestScore = score;
		}
	}
	Templates = bestTemplates_r;
	i = *Templates++ + 1;
	while(--i) {
		if(bestScore < (score = Score_r[*Templates++])) {
			bestScore = score;
		}
	}
	proxiScore = minFrac * bestScore;
	
	/* get best matches within proximity */
	bestHits = 0;
	Templates = bestTemplates;
	i = *Templates++ + 1;
	while(--i) {
		score = Score[(template = *Templates++)];
		if(proxiScore <= score) {
			bestTemplates[++bestHits] = template;
			if(softProxi) {
				softProxi[template] += score;
			}
		} else {
			Score[template] = 0;
		}
	}
	Templates = bestTemplates_r;
	i = *Templates++ + 1;
	while(--i) {
		score = Score_r[(template = *Templates++)];
		if(proxiScore <= score) {
			bestTemplates[++bestHits] = -template;
			if(softProxi) {
				softProxi[template] += score;
			}
		} else {
			Score_r[template] = 0;
		}
	}
	*bestTemplates = bestHits;
	
	/* check union */
	bestHits = 0;
	for(i = 1; i <= *regionTemplates; ++i) {
		if(0 < regionTemplates[i]) {
			/* we got one */
			if(Score_r[regionTemplates[i]]) {
				++bestHits;
				score = regionTemplates[bestHits];
				regionTemplates[bestHits] = regionTemplates[i];
				regionTemplates[i] = score;
			}
		} else {
			/* we got one */
			if(Score[-regionTemplates[i]]) {
				++bestHits;
				score = regionTemplates[bestHits];
				regionTemplates[bestHits] = regionTemplates[i];
				regionTemplates[i] = score;
			}
		}
	}
	
	/* mark as PE */
	if(bestHits) {
		*regionTemplates = -bestHits;
	}
	
	/* clear scores */
	for(i = *bestTemplates; i != 0; --i) {
		if(0 < bestTemplates[i]) {
			Score[bestTemplates[i]] = 0;
		} else {
			Score_r[-bestTemplates[i]] = 0;
		}
	}
	
	return bestScore;
}

int save_kmers_Sparse(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, CompDNA *qseq_r, Qseqs *header, int *extendScore, const int exhaustive, volatile int *excludeOut, FILE *out) {
	
	int i, j, k, l, n, end, rc, prefix_len, template, hitCounter, HIT, SU;
	int M, MM, n_kmers, bestScore, reps, kmersize, mlen, tflag, flag, seqend;
	int cPos, iPos, mPos, hLen;
	unsigned shifter, prefix_shifter, *values, *last;
	short unsigned *values_s;
	long unsigned pmask, mask, mmask, prefix, pmer, kmer, cmer, hmer, *seq;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return 1;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	seq = qseq->seq;
	seqend = qseq->seqlen - kmersize + 1;
	prefix = templates->prefix;
	prefix_len = templates->prefix_len;
	prefix_shifter = 64 - (prefix_len << 1);
	pmask = 0xFFFFFFFFFFFFFFFF >> prefix_shifter;
	kmersize = templates->kmersize;
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	tflag = templates->flag;
	hLen = kmersize;
	M = rewards->M;
	MM = rewards->MM;
	*bestTemplates = 0;
	hitCounter = 0;
	n_kmers = 0;
	bestScore = 0;
	end = qseq->seqlen;
	if(prefix_len) {
		flag = 16;
		for(rc = 0; rc < 2; ++rc) {
			if(rc) {
				comp_rc(qseq);
			}
			j = 0;
			qseq->N[0]++;
			qseq->N[qseq->N[0]] = qseq->seqlen;
			for(i = 1; i <= qseq->N[0]; ++i) {
				getKmer_macro(pmer, seq, j, cPos, iPos, (prefix_shifter + 2));
				end = qseq->N[i] - kmersize;
				for(j += prefix_len - 1; j < end; ++j) {
					pmer = updateKmer_macro(pmer, seq, j, pmask);
					if(pmer == prefix) {
						/* get kmer */
						getKmer_macro(kmer, seq, (j + 1), cPos, iPos, shifter);
						cmer = tflag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
						if((values = hashMap_get(templates, cmer))) {
							if(SU) {
								values_s = (short unsigned *) values;
								n = *values_s;
								for(k = 1; k <= n; ++k) {
									if(Score[(template = values_s[k])] == 0) {
										(*bestTemplates)++;
										bestTemplates[*bestTemplates] = template;
									}
									Score[template]++;
								}
							} else {
								n = *values;
								for(k = 1; k <= n; ++k) {
									if(Score[(template = values[k])] == 0) {
										(*bestTemplates)++;
										bestTemplates[*bestTemplates] = template;
									}
									Score[template]++;
								}
							}
							++hitCounter;
						}
						++n_kmers;
					}
				}
				j = qseq->N[i] + 1;
			}
			qseq->N[0]--;
		}
		
		/* get best match(es) */
		/*if(hitCounter * kmersize > (n_kmers - hitCounter)) {*/
		if(hitCounter) {
			bestScore = getMatchSparse(bestTemplates, Score, kmersize, n_kmers, M, MM);
			
			/*
			bestHits = 0;
			for(l = 1; l <= *bestTemplates; ++l) {
				score = Score[(template = bestTemplates[l])];
				score = score * kmersize * M + (n_kmers - score) * MM;
				if(score > bestScore) {
					bestScore = score;
					bestHits = 1;
					bestTemplates[bestHits] = template;
				} else if(score == bestScore) {
					++bestHits;
					bestTemplates[bestHits] = template;
				}
				Score[template] = 0;
			}
			*bestTemplates = bestHits;
			*/
		} else {
			/*
			for(l = *bestTemplates; l != 0; --l) {
				Score[bestTemplates[l]] = 0;
			}
			*bestTemplates = 0;
			*/
			*bestTemplates = clearScore(bestTemplates, Score);
		}
		end = n_kmers - hitCounter - bestScore;
	} else {
		flag = 0;
		HIT = exhaustive;
		j = 0;
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(; j < end && !HIT; j += kmersize) {
				getKmer_macro(kmer, seq, j, cPos, iPos, shifter);
				cmer = tflag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
				if(hashMap_get(templates, cmer)) {
					HIT = 1;
				}
			}
			j = qseq->N[i] + 1;
		}
		
		if(HIT) {
			last = 0;
			reps = 0;
			j = 0;
			for(i = 1; i <= qseq->N[0] && j < seqend; ++i) {
				getKmer_macro(kmer, seq, j, cPos, iPos, (shifter + 2));
				cmer = tflag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
				end = qseq->N[i];
				for(j += kmersize - 1; j < end; ++j) {
					/* update kmer */
					kmer = updateKmer_macro(kmer, seq, j, mask);
					cmer = tflag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
					if((values = hashMap_get(templates, cmer))) {
						if(values == last) {
							++reps;
						} else {
							if(last) {
								if(SU) {
									values_s = (short unsigned *) values;
									n = *values_s;
									for(l = 1; l <= n; ++l) {
										if(Score[(template = values_s[l])]) {
											Score[template] += reps;
										} else {
											Score[template] = reps;
											bestTemplates[0]++;
											bestTemplates[*bestTemplates] = template;
										}
									}
								} else {
									n = *values;
									for(l = 1; l <= n; ++l) {
										if(Score[(template = last[l])]) {
											Score[template] += reps;
										} else {
											Score[template] = reps;
											bestTemplates[0]++;
											bestTemplates[*bestTemplates] = template;
										}
									}
								}
								hitCounter += reps;
							}
							reps = 1;
							last = values;
						}
					}
				}
				j = end + 1;
			}
			if(last) {
				if(SU) {
					values_s = (short unsigned *) last;
					n = *values_s;
					for(l = 1; l <= n; ++l) {
						if(Score[(template = values_s[l])]) {
							Score[template] += reps;
						} else {
							Score[template] = reps;
							bestTemplates[0]++;
							bestTemplates[*bestTemplates] = template;
						}
					}
				} else {
					n = *last;
					for(l = 1; l <= n; ++l) {
						if(Score[(template = last[l])]) {
							Score[template] += reps;
						} else {
							Score[template] = reps;
							bestTemplates[0]++;
							bestTemplates[*bestTemplates] = template;
						}
					}
				}
				hitCounter += reps;
			}
		}
		qseq->N[0]--;
		
		/* get best match(es) */
		/*if(hitCounter * kmersize > (end - hitCounter + kmersize)) {*/
		if(hitCounter) {
			bestScore = getMatch(bestTemplates, Score);
		} else {
			*bestTemplates = clearScore(bestTemplates, Score);
		}
		end = qseq->seqlen + 1 - bestScore;
	}
	
	i = 1;
	if(kmersize <= bestScore || bestScore * kmersize > end) {
		lock(excludeOut);
		i = deConPrintPtr(bestTemplates, qseq, bestScore, header, flag, out);
		unlock(excludeOut);
	}
	
	return i;
}

int save_kmers_pseuodeSparse(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, CompDNA *qseq_r, Qseqs *header, int *extendScore, const int exhaustive, volatile int *excludeOut, FILE *out) {
	
	int i, j, l, n, end, template, hitCounter, gaps, Ms, MMs, Us, W1s, flag;
	int HIT, SU, score, bestScore, kmersize, mlen, mPos, hLen, cPos, iPos;
	int W1, U, M, MM, seqend, j_u, m, mm;
	unsigned shifter, *values, *last;
	short unsigned *values_s;
	char *include;
	long unsigned mask, mmask, kmer, cmer, hmer, *seq;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return 1;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	seq = qseq->seq;
	seqend = qseq->seqlen - kmersize + 1;
	kmersize = templates->kmersize;
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = templates->flag;
	W1 = rewards->W1;
	U = rewards->U;
	M = rewards->M;
	MM = rewards->MM;
	*bestTemplates = 0;
	hitCounter = 0;
	bestScore = 0;
	end = qseq->seqlen;
	HIT = exhaustive;
	include = (char *) (extendScore + (templates->DB_size + 1));
	hLen = kmersize;
	
	j = 0;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
		end = qseq->N[i] - kmersize + 1;
		for(; j < end && !HIT; j += kmersize) {
			getKmer_macro(kmer, seq, j, cPos, iPos, shifter);
			cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
			if(hashMap_get(templates, cmer)) {
				HIT = 1;
			}
		}
		j = qseq->N[i] + 1;
	}
	
	if(HIT) {
		last = 0;
		gaps = 0;
		HIT = 0;
		Ms = 0;
		MMs = 0;
		Us = 0;
		W1s = 0;
		j = 0;
		for(i = 1; i <= qseq->N[0] && j < seqend; ++i) {
			/* init k-mer */
			getKmer_macro(kmer, seq, j, cPos, iPos, (shifter + 2));
			cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
			end = qseq->N[i];
			for(j_u = j + kmersize - 1; j_u < end; ++j_u) {
				/* update k-mer */
				kmer = updateKmer_macro(kmer, seq, j_u, mask);
				cmer = flag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
				/* lookup */
				if((values = hashMap_get(templates, cmer))) {
					if(values == last) {
						/*
						gaps == 0 -> Match
						gaps == kmersize -> 1 MM
						kmersize < gaps -> several mismatches or indel(s)
						gaps < kmersize -> deletion
						*/
						if(gaps == 0) {
							/* match */
							++Ms;
						} else if(gaps == kmersize) {
							/* snp */
							Ms += kmersize;
							++MMs;
						} else if(kmersize < gaps) {
							/* mismatch or insersion */
							Ms += kmersize;
							gaps -= (kmersize - 1); /* adjust for consecutive k-mer mismatches */
							if(gaps <= 2) {
								mm = gaps;
								m = 0;
							} else {
								mm = gaps / kmersize + (gaps % kmersize ? 1 : 0);
								mm = MAX(2, mm);
								m = MIN(gaps - mm, kmersize);
								m = MIN(m, mm);
							}
							
							/* evaluate best option */
							if((W1 + (gaps - 1) * U) <= (mm * MM + m * M)) {
								MMs += mm;
								Ms += m;
							} else {
								++W1s;
								Us += (gaps -1);
							}
						} else if(mlen != kmersize) {
							Ms += gaps;
							++MMs;
							/* unlikely deletion or random k-mer mismatch, 
							unless minimizers or homopolymer compression,
							assume better and go random zero score else. */
						} else {
							/* indel */
							Ms += gaps;
							++W1s;
							Us += (kmersize - gaps);
						}
						
						HIT = j;
						gaps = 0;
					} else {
						if(last) {
							score = Ms * M + MMs * MM + Us * U + W1s * W1;
							values_s = (short unsigned *) last;
							l = SU ? (*values_s + 1) : (*last + 1);
							while(--l) {
								template = SU ? *++values_s : *++last;
								Score[template] += score;
								extendScore[template] = HIT;
							}
							HIT = j - 1;
							last = values;
							values_s = (short unsigned *) values;
							n = SU ? *values_s : *values;
							l = n + 1;
							while(--l) {
								template = SU ? *++values_s : *++values;
								if(include[template]) {
									gaps = HIT - extendScore[template];
									if(gaps == 0) {
										/* match */
										score = M;
									} else if(gaps == kmersize) {
										/* snp */
										score = kmersize * M + MM;
									} else if(kmersize < gaps) {
										/* mismatch or insersion */
										gaps -= (kmersize - 1); /* adjust for consecutive k-mer mismatches */
										if(gaps <= 2) {
											mm = gaps;
											m = 0;
										} else {
											mm = gaps / kmersize + (gaps % kmersize ? 1 : 0);
											mm = MAX(2, mm);
											m = MIN(gaps - mm, kmersize);
											m = MIN(m, mm);
										}
										
										/* evaluate best option */
										if((W1 + (gaps - 1) * U) <= (mm * MM + m * M)) {
											score = kmersize * M + (mm * MM + m * M);
										} else {
											score = kmersize * M + (W1 + (gaps - 1) * U);
										}
									} else if(mlen != kmersize) {
										score = gaps * M + MM;
										/* unlikely deletion or random k-mer mismatch, 
										unless minimizers or homopolymer compression,
										assume better and go random zero score else. */
									} else {
										/* indel */
										score = gaps * M + (kmersize - gaps) * U + W1;
									}
									Score[template] += score;
								} else {
									Score[template] = kmersize * M;
									if(include[template] == 0) {
										include[template] = 1;
										bestTemplates[++*bestTemplates] = template;
									}
								}
							}
						} else {
							last = values;
							values_s = (short unsigned *) values;
							n = SU ? *values_s : *values;
							Ms = kmersize * M;
							for(l = 1; l <= n; ++l) {
								template = SU ? *++values_s : *++values;
								Score[template] = Ms;
								include[template] = 1;
								bestTemplates[l] = template;
							}
							*bestTemplates = n;
						}
						
						HIT = j;
						gaps = 0;
						Ms = 0;
						MMs = 0;
						Us = 0;
						W1s = 0;
					}
					++hitCounter;
				} else {
					++gaps;
				}
				++j;
			}
			gaps += (qseq->N[i] + 1 - j); /* gap over N's */
			j = end + 1;
		}
		if(last) {
			score = Ms * M + MMs * MM + Us * U + W1s * W1;
			if(SU) {
				values_s = (short unsigned *) last;
				l = (*values_s) + 1;
				while(--l) {
					Score[values_s[l]] += score;
				}
			} else {
				l = (*last) + 1;
				while(--l) {
					Score[last[l]] += score;
				}
			}
			for(l = *bestTemplates; l != 0; --l) {
				extendScore[(template = bestTemplates[l])] = 0;
				include[template] = 0;
				if(Score[template] < 0) {
					Score[template] = 0;
				}
			}
		}
	}
	qseq->N[0]--;
	
	/* get best match(es) */
	/*clearScore(bestTemplates, extendScore);*/
	/*if(hitCounter * kmersize > (end - hitCounter + kmersize)) {*/
	if(hitCounter) {
		bestScore = getMatch(bestTemplates, Score);
		/*
		bestHits = 0;
		for(l = 1; l <= *bestTemplates; ++l) {
			if(Score[(template = bestTemplates[l])] > bestScore) {
				bestScore = Score[template];
				bestHits = 1;
				bestTemplates[bestHits] = template;
			} else if(Score[template] == bestScore) {
				++bestHits;
				bestTemplates[bestHits] = template;
			}
			Score[template] = 0;
			extendScore[template] = 0;
		}
		*bestTemplates = bestHits;
		*/
	} else {
		*bestTemplates = clearScore(bestTemplates, Score);
		/*
		for(l = *bestTemplates; l != 0; --l) {
			extendScore[(template = bestTemplates[l])] = 0;
			Score[template] = 0;
		}
		*bestTemplates = 0;
		*/
	}
	end = qseq->seqlen + 1 - bestScore;
	
	i = 1;
	if(kmersize <= bestScore || bestScore * kmersize > end) {
		lock(excludeOut);
		i = deConPrintPtr(bestTemplates, qseq, bestScore, header, 0, out);
		unlock(excludeOut);
	}
	
	return i;
}

int save_kmers(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, CompDNA *qseq_r, Qseqs *header, int *extendScore, const int exhaustive, volatile int *excludeOut, FILE *out) {
	
	int i, j, j_u, l, end, HIT, gaps, score, Ms, MMs, Us, W1s, W1, U, M, MM;
	int template, hitCounter, bestScore, bestScore_r, kmersize, mPos, hLen;
	int seqend, m, mm, mlen;
	unsigned *values, *last, n, SU, shifter, cPos, iPos, flag;
	short unsigned *values_s;
	long unsigned mask, mmask, kmer, cmer, hmer, *seq;
	char *include;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return 1;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	bestScore = 0;
	bestScore_r = 0;
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = templates->flag;
	W1 = rewards->W1;
	U = rewards->U;
	M = rewards->M;
	MM = rewards->MM;
	include = (char *) (extendScore + (templates->DB_size + 1));
	hLen = kmersize;
	
	/* reverse complement qseq */
	rc_comp(qseq, qseq_r);
	
	/* Search forward strand */
	/* Make quick check of the qseq */
	seq = qseq->seq;
	seqend = qseq->seqlen - kmersize + 1;
	HIT = exhaustive;
	j = 0;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
		end = qseq->N[i] - kmersize + 1;
		for(; j < end && !HIT; j += kmersize) {
			getKmer_macro(kmer, seq, j, cPos, iPos, shifter);
			cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
			if(hashMap_get(templates, cmer)) {
				HIT = 1;
			}
		}
		j = qseq->N[i] + 1;
	}
	
	/* If deltamer qseq hits, then continue */
	if(HIT) {
		/* Scan the deltamer exhaustively, and collect scores in Score*/
		hitCounter = 0;
		*bestTemplates = 0;
		last = 0;
		gaps = 0;
		HIT = 0;
		Ms = 0;
		MMs = 0;
		Us = 0;
		W1s = 0;
		j = 0;
		end = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && j < seqend; ++i) {
			/* init k-mer */
			getKmer_macro(kmer, seq, j, cPos, iPos, (shifter + 2));
			cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
			end = qseq->N[i];
			for(j_u = j + kmersize - 1; j_u < end; ++j_u) {
				/* update k-mer */
				kmer = updateKmer_macro(kmer, seq, j_u, mask);
				cmer = flag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
				/* lookup */
				if((values = hashMap_get(templates, cmer))) {
					if(values == last) {
						/*
						gaps == 0 -> Match
						gaps == kmersize -> 1 MM
						kmersize < gaps -> several mismatches or indel(s)
						gaps < kmersize -> deletion
						*/
						if(gaps == 0) {
							/* match */
							++Ms;
						} else if(mlen <= gaps && gaps <= kmersize) {
							/* snp */
							Ms += kmersize;
							++MMs;
						} else if(kmersize < gaps) {
							/* mismatch or insersion */
							Ms += kmersize;
							gaps -= (kmersize - 1); /* adjust for consecutive k-mer mismatches */
							if(gaps <= 2) {
								mm = gaps;
								m = 0;
							} else {
								mm = gaps / kmersize + (gaps % kmersize ? 1 : 0);
								mm = MAX(2, mm);
								m = MIN(gaps - mm, kmersize);
								m = MIN(m, mm);
							}
							
							/* evaluate best option */
							if((W1 + (gaps - 1) * U) <= (mm * MM + m * M)) {
								MMs += mm;
								Ms += m;
							} else {
								++W1s;
								Us += (gaps - 1);
							}
						} else if(mlen != kmersize) {
							Ms += gaps;
							++MMs;
							/* unlikely deletion or random k-mer mismatch, 
							unless minimizers or homopolymer compression,
							assume better and go random zero score else. */
						} else {
							/* indel */
							Ms += gaps;
							++W1s;
							Us += (kmersize - gaps);
						}
						
						HIT = j;
						gaps = 0;
					} else {
						if(last) {
							score = Ms * M + MMs * MM + Us * U + W1s * W1;
							values_s = (short unsigned *) last;
							l = SU ? (*values_s + 1) : (*last + 1);
							while(--l) {
								template = SU ? *++values_s : *++last;
								Score[template] += score;
								extendScore[template] = HIT;
							}
							HIT = j - 1;
							last = values;
							values_s = (short unsigned *) values;
							n = SU ? *values_s : *values;
							l = n + 1;
							while(--l) {
								template = SU ? *++values_s : *++values;
								if(include[template]) {
									gaps = HIT - extendScore[template];
									if(gaps == 0) {
										/* match */
										score = M;
									} else if(mlen <= gaps && gaps <= kmersize) {
										/* snp */
										score = gaps * M + MM;
									} else if(kmersize < gaps) {
										/* mismatch or insersion */
										gaps -= (kmersize - 1); /* adjust for consecutive k-mer mismatches */
										if(gaps <= 2) {
											mm = gaps;
											m = 0;
										} else {
											mm = gaps / kmersize + (gaps % kmersize ? 1 : 0);
											mm = MAX(2, mm);
											m = MIN(gaps - mm, kmersize);
											m = MIN(m, mm);
										}
										
										/* evaluate best option */
										if((W1 + (gaps - 1) * U) <= (mm * MM + m * M)) {
											score = kmersize * M + (mm * MM + m * M);
										} else {
											score = kmersize * M + (W1 + (gaps - 1) * U);
										}
									} else if(mlen != kmersize) {
										score = gaps * M + MM;
										/* unlikely deletion or random k-mer mismatch, 
										unless minimizers or homopolymer compression,
										assume better and go random zero score else. */
									} else {
										/* indel */
										score = gaps * M + (kmersize - gaps) * U + W1;
									}
									Score[template] += score;
								} else {
									/* consider opening penalty */
									/*
									gaps = j;
									if(gaps == 0) {
										Score[template] = score;
									} else {
										if(gaps <= 2) {
											mm = gaps;
											m = 0;
										} else {
											mm = gaps / kmersize + (gaps % kmersize ? 1 : 0);
											mm = MAX(2, mm);
											m = MIN(gaps - mm, kmersize);
											m = MIN(m, mm);
										}
										// evaluate best option
										open = ((W1 + (gaps - 1) * U) < Wl) ? Wl : (W1 + (gaps - 1) * U);
										open = open < (mm * MM + m * M) ? (mm * MM + m * M) : open;
										Score[template] = score + open;
									}
									*/
									Score[template] = kmersize * M;
									if(include[template] == 0) {
										include[template] = 1;
										bestTemplates[++*bestTemplates] = template;
									}
								}
							}
						} else {
							last = values;
							values_s = (short unsigned *) values;
							n = SU ? *values_s : *values;
							score = kmersize * M;
							
							/* consider opening penalty */
							/*
							gaps = j;
							if(gaps > 0) {
								if(gaps <= 2) {
									mm = gaps;
									m = 0;
								} else {
									mm = gaps / kmersize + (gaps % kmersize ? 1 : 0);
									mm = MAX(2, mm);
									m = MIN(gaps - mm, kmersize);
									m = MIN(m, mm);
								}
								// evaluate best option
								open = ((W1 + (gaps - 1) * U) < Wl) ? Wl : (W1 + (gaps - 1) * U);
								open = open < (mm * MM + m * M) ? (mm * MM + m * M) : open;
								score += open;
							}
							*/
							
							for(l = 1; l <= n; ++l) {
								template = SU ? *++values_s : *++values;
								Score[template] = score;
								include[template] = 1;
								bestTemplates[l] = template;
							}
							*bestTemplates = n;
						}
						
						HIT = j;
						gaps = 0;
						Ms = 0;
						MMs = 0;
						Us = 0;
						W1s = 0;
					}
					++hitCounter;
				} else {
					++gaps;
				}
				++j;
			}
			gaps += (qseq->N[i] + 1 - j); /* gap over N's */
			j = end + 1;
		}
		if(last) {
			score = Ms * M + MMs * MM + Us * U + W1s * W1;
			if(SU) {
				values_s = (short unsigned *) last;
				l = (*values_s) + 1;
				while(--l) {
					Score[(template = values_s[l])] += score;
					extendScore[template] = HIT;
				}
			} else {
				l = (*last) + 1;
				while(--l) {
					Score[(template = last[l])] += score;
					extendScore[template] = HIT;
				}
			}
			
			for(l = *bestTemplates; l != 0; --l) {
				template = bestTemplates[l];
				/* consider closing penalties */
				/*
				HIT = qseq->seqlen - kmersize;
				gaps = HIT - extendScore[template];
				if(gaps > 0) {
					if(gaps <= 2) {
						mm = gaps;
						m = 0;
					} else {
						mm = gaps / kmersize + (gaps % kmersize ? 1 : 0);
						mm = MAX(2, mm);
						m = MIN(gaps - mm, kmersize);
						m = MIN(m, mm);
					}
					// evaluate best option
					open = ((W1 + (gaps - 1) * U) < Wl) ? Wl : (W1 + (gaps - 1) * U);
					open = open < (mm * MM + m * M) ? (mm * MM + m * M) : open;
					Score[template] += open;
				}
				*/
				
				/* clear up */
				extendScore[template] = 0;
				include[template] = 0;
				if(Score[template] < 0) {
					Score[template] = 0;
				}
			}
		}
		
		/* get best match(es) */
		/*clearScore(bestTemplates, extendScore);*/
		/*if(hitCounter * kmersize > (end - hitCounter + kmersize)) {*/
		if(hitCounter) {
			bestScore = getMatch(bestTemplates, Score);
			/*
			bestHits = 0;
			for(l = 1; l <= *bestTemplates; ++l) {
				if(Score[(template = bestTemplates[l])] > bestScore) {
					bestScore = Score[template];
					bestHits = 1;
					bestTemplates[bestHits] = template;
				} else if(Score[template] == bestScore) {
					++bestHits;
					bestTemplates[bestHits] = template;
				}
				Score[template] = 0;
				extendScore[template] = 0;
			}
			*bestTemplates = bestHits;
			*/
		} else {
			/*
			for(l = *bestTemplates; l != 0; --l) {
				extendScore[(template = bestTemplates[l])] = 0;
				Score[template] = 0;
			}
			*/
			*bestTemplates = clearScore(bestTemplates, Score);
		}
	}
	qseq->N[0]--;
	
	/* search rc strand */
	/* Make quick check of the qseq */
	seq = qseq_r->seq;
	HIT = exhaustive;
	j = 0;
	qseq_r->N[0]++;
	qseq_r->N[qseq_r->N[0]] = qseq_r->seqlen;
	for(i = 1; i <= qseq_r->N[0] && !HIT; ++i) {
		end = qseq_r->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			getKmer_macro(kmer, seq, j, cPos, iPos, shifter);
			cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
			if(hashMap_get(templates, cmer)) {
				HIT = 1;
			}
		}
		j = qseq_r->N[i] + 1;
	}
	
	/* If deltamer qseq hits, then continue */
	if(HIT) {
		/* Scan the deltamer exhaustively, and collect scores in Score*/
		hitCounter = 0;
		*bestTemplates_r = 0;
		last = 0;
		gaps = 0;
		HIT = 0;
		Ms = 0;
		MMs = 0;
		Us = 0;
		W1s = 0;
		j = 0;
		end = qseq_r->seqlen;
		for(i = 1; i <= qseq_r->N[0] && j < seqend; ++i) {
			/* init k-mer */
			getKmer_macro(kmer, seq, j, cPos, iPos, (shifter + 2));
			cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
			end = qseq_r->N[i];
			for(j_u = j + kmersize - 1; j_u < end; ++j_u) {
				/* update k-mer */
				kmer = updateKmer_macro(kmer, seq, j_u, mask);
				cmer = flag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
				/* lookup */
				if((values = hashMap_get(templates, cmer))) {
					if(values == last) {
						/*
						gaps == 0 -> Match
						gaps == kmersize -> 1 MM
						kmersize < gaps -> several mismatches or indel(s)
						gaps < kmersize -> deletion
						*/
						if(gaps == 0) {
							/* match */
							++Ms;
						} else if(mlen <= gaps && gaps <= kmersize) {
							/* snp */
							Ms += kmersize;
							++MMs;
						} else if(kmersize < gaps) {
							/* mismatch or insersion */
							Ms += kmersize;
							gaps -= (kmersize - 1); /* adjust for consecutive k-mer mismatches */
							if(gaps <= 2) {
								mm = gaps;
								m = 0;
							} else {
								mm = gaps / kmersize + (gaps % kmersize ? 1 : 0);
								mm = MAX(2, mm);
								m = MIN(gaps - mm, kmersize);
								m = MIN(m, mm);
							}
							
							/* evaluate best option */
							if((W1 + (gaps - 1) * U) <= (mm * MM + m * M)) {
								MMs += mm;
								Ms += m;
							} else {
								++W1s;
								Us += (gaps - 1);
							}
						} else if(mlen != kmersize) {
							Ms += gaps;
							++MMs;
							/* unlikely deletion or random k-mer mismatch, 
							unless minimizers or homopolymer compression,
							assume better and go random zero score else. */
						} else {
							/* indel */
							Ms += gaps;
							++W1s;
							Us += (kmersize - gaps);
						}
						
						HIT = j;
						gaps = 0;
					} else {
						if(last) {
							score = Ms * M + MMs * MM + Us * U + W1s * W1;
							values_s = (short unsigned *) last;
							l = SU ? (*values_s + 1) : (*last + 1);
							while(--l) {
								template = SU ? *++values_s : *++last;
								Score_r[template] += score;
								extendScore[template] = HIT;
							}
							HIT = j - 1;
							last = values;
							values_s = (short unsigned *) values;
							n = SU ? *values_s : *values;
							l = n + 1;
							while(--l) {
								template = SU ? *++values_s : *++values;
								if(include[template]) {
									gaps = HIT - extendScore[template];
									if(gaps == 0) {
										/* match */
										score = M;
									} else if(mlen <= gaps && gaps <= kmersize) {
										/* snp */
										score = kmersize * M + MM;
									} else if(kmersize < gaps) {
										/* mismatch or insersion */
										gaps -= (kmersize - 1); /* adjust for consecutive k-mer mismatches */
										if(gaps <= 2) {
											mm = gaps;
											m = 0;
										} else {
											mm = gaps / kmersize + (gaps % kmersize ? 1 : 0);
											mm = MAX(2, mm);
											m = MIN(gaps - mm, kmersize);
											m = MIN(m, mm);
										}
										
										/* evaluate best option */
										if((W1 + (gaps - 1) * U) <= (mm * MM + m * M)) {
											score = kmersize * M + (mm * MM + m * M);
										} else {
											score = kmersize * M + (W1 + (gaps - 1) * U);
										}
									} else if(mlen != kmersize) {
										score = gaps * M + MM;
										/* unlikely deletion or random k-mer mismatch, 
										unless minimizers or homopolymer compression,
										assume better and go random zero score else. */
									} else {
										/* indel */
										score = gaps * M + (kmersize - gaps) * U + W1;
									}
									Score_r[template] += score;
								} else {
									Score_r[template] = kmersize * M;
									if(include[template] == 0) {
										include[template] = 1;
										bestTemplates_r[++*bestTemplates_r] = template;
									}
								}
							}
						} else {
							last = values;
							values_s = (short unsigned *) values;
							n = SU ? *values_s : *values;
							score = kmersize * M;
							for(l = 1; l <= n; ++l) {
								template = SU ? *++values_s : *++values;
								Score_r[template] = score;
								include[template] = 1;
								bestTemplates_r[l] = template;
							}
							*bestTemplates_r = n;
						}
						
						HIT = j;
						gaps = 0;
						Ms = 0;
						MMs = 0;
						Us = 0;
						W1s = 0;
					}
					++hitCounter;
				} else {
					++gaps;
				}
				++j;
			}
			gaps += (qseq_r->N[i] + 1 - j); /* gap over N's */
			j = end + 1;
		}
		if(last) {
			score = Ms * M + MMs * MM + Us * U + W1s * W1;
			if(SU) {
				values_s = (short unsigned *) last;
				l = (*values_s) + 1;
				while(--l) {
					Score_r[values_s[l]] += score;
				}
			} else {
				l = (*last) + 1;
				while(--l) {
					Score_r[last[l]] += score;
				}
			}
			
			for(l = *bestTemplates_r; l != 0; --l) {
				template = bestTemplates_r[l];
				
				/* clear up */
				extendScore[template] = 0;
				include[template] = 0;
				if(Score_r[template] < 0) {
					Score_r[template] = 0;
				}
			}
		}
		
		/* get best match(es) */
		/*clearScore(bestTemplates_r, extendScore);*/
		/*if(hitCounter * kmersize > (end - hitCounter + kmersize)) {*/
		if(hitCounter) {
			bestScore_r = getMatch(bestTemplates_r, Score_r);
			
			/*
			bestHits = 0;
			for(l = 1; l <= *bestTemplates_r; ++l) {
				if(Score_r[(template = bestTemplates_r[l])] > bestScore_r) {
					bestScore_r = Score_r[template];
					bestHits = 1;
					bestTemplates_r[bestHits] = template;
				} else if(Score_r[template] == bestScore_r) {
					++bestHits;
					bestTemplates_r[bestHits] = template;
				}
				Score_r[template] = 0;
				extendScore[template] = 0;
			}
			*bestTemplates_r = bestHits;
			*/
		} else {
			/*
			for(l = *bestTemplates_r; l != 0; --l) {
				extendScore[(template = bestTemplates_r[l])] = 0;
				Score_r[template] = 0;
			}
			*/
			*bestTemplates_r = clearScore(bestTemplates_r, Score_r);
		}
	}
	qseq_r->N[0]--;
	
	/* Validate best match */
	i = 1;
	if(bestScore > 0 || bestScore_r > 0) {
		end = qseq->seqlen + 1;
		//if((bestScore >= bestScore_r && bestScore * kmersize > (end - bestScore)) || (bestScore < bestScore_r && bestScore_r * kmersize > (end - bestScore_r))) {
		if(kmersize <= bestScore || kmersize <= bestScore_r) {
			if(bestScore > bestScore_r) {
				lock(excludeOut);
				i = deConPrintPtr(bestTemplates, qseq, bestScore, header, 0, out);
				unlock(excludeOut);
			} else if(bestScore < bestScore_r) {
				lock(excludeOut);
				i = deConPrintPtr(bestTemplates_r, qseq_r, bestScore_r, header, 16, out);
				unlock(excludeOut);
			} else {
				/* merge */
				for(i = 1; i <= *bestTemplates_r; ++i) {
					bestTemplates[0]++;
					bestTemplates[*bestTemplates] = -bestTemplates_r[i];
				}
				lock(excludeOut);
				i = deConPrintPtr(bestTemplates, qseq, -bestScore, header, 0, out);
				unlock(excludeOut);
			}
		}
	}
	
	return i;
}

int save_kmers_count(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, CompDNA *qseq_r, Qseqs *header, int *extendScore, const int exhaustive, volatile int *excludeOut, FILE *out) {
	
	int i, j, end, hitCounter, bestScore, bestScore_r, reps, n, mPos, hLen;
	int seqend;
	unsigned kmersize, mlen, shifter, flag, SU, HIT, iPos, cPos;
	unsigned *values, *last;
	short unsigned *values_s;
	long unsigned mask, mmask, kmer, cmer, hmer, *seq;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return 1;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	bestScore = 0;
	bestScore_r = 0;
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = templates->flag;
	hLen = kmersize;
	
	/* reverse complement qseq */
	rc_comp(qseq, qseq_r);
	
	/* Search forward strand */
	/* Make quick check of the qseq */
	seq = qseq->seq;
	seqend = qseq->seqlen - kmersize + 1;
	HIT = exhaustive;
	j = 0;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
		end = qseq->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			getKmer_macro(kmer, seq, j, cPos, iPos, shifter);
			cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
			if(hashMap_get(templates, cmer)) {
				HIT = 1;
			}
		}
		j = qseq->N[i] + 1;
	}
	
	/* If deltamer qseq hits, then continue */
	if(HIT) {
		/* Scan the deltamer exhaustively, and collect scores in Score*/
		hitCounter = 0;
		*bestTemplates = 0;
		last = 0;
		reps = 0;
		j = 0;
		end = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && j < seqend; ++i) {
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
					} else {
						if(last) {
							if(SU) {
								values_s = (short unsigned *) last;
								n = *values_s + 1;
								while(--n) {
									if((Score[*++values_s] += reps) == reps) {
										bestTemplates[++*bestTemplates] = *values_s;
									}
								}
							} else {
								n = *last + 1;
								while(--n) {
									if((Score[*++last] += reps) == reps) {
										bestTemplates[++*bestTemplates] = *last;
									}
								}
							}
							hitCounter += reps;
						}
						reps = 1;
						last = values;
					}
				}
			}
			j = end + 1;
		}
		if(last) {
			if(SU) {
				values_s = (short unsigned *) last;
				n = *values_s + 1;
				while(--n) {
					if((Score[*++values_s] += reps) == reps) {
						bestTemplates[++*bestTemplates] = *values_s;
					}
				}
			} else {
				n = *last + 1;
				while(--n) {
					if((Score[*++last] += reps) == reps) {
						bestTemplates[++*bestTemplates] = *last;
					}
				}
			}
			hitCounter += reps;
		}
		reps = 0;
		
		/* get best match(es) */
		/*if(hitCounter * kmersize > (end - hitCounter + kmersize)) {*/
		if(hitCounter) {
			
			bestScore = getMatch(bestTemplates, Score);
			
			/*
			bestHits = 0;
			for(l = 1; l <= *bestTemplates; ++l) {
				if(Score[(template = bestTemplates[l])] > bestScore) {
					bestScore = Score[template];
					bestHits = 1;
					bestTemplates[bestHits] = template;
				} else if(Score[template] == bestScore) {
					++bestHits;
					bestTemplates[bestHits] = template;
				}
				Score[template] = 0;
			}
			*bestTemplates = bestHits;
			*/
			
		} else {
			/*
			for(l = *bestTemplates; l != 0; --l) {
				Score[bestTemplates[l]] = 0;
			}
			*/
			*bestTemplates = clearScore(bestTemplates, Score);
		}
	}
	qseq->N[0]--;
	
	/* search rc strand */
	/* Make quick check of the qseq */
	seq = qseq_r->seq;
	HIT = exhaustive;
	j = 0;
	qseq_r->N[0]++;
	qseq_r->N[qseq_r->N[0]] = qseq_r->seqlen;
	for(i = 1; i <= qseq_r->N[0] && !HIT; ++i) {
		end = qseq_r->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			getKmer_macro(kmer, seq, j, cPos, iPos, shifter);
			cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
			if(hashMap_get(templates, cmer)) {
				HIT = 1;
			}
		}
		j = qseq_r->N[i] + 1;
	}
	
	/* If deltamer qseq hits, then continue */
	if(HIT) {
		/* Scan the deltamer exhaustively, and collect scores in Score*/
		hitCounter = 0;
		*bestTemplates_r = 0;
		last = 0;
		reps = 0;
		j = 0;
		end = qseq_r->seqlen;
		for(i = 1; i <= qseq_r->N[0] && j < seqend; ++i) {
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
					} else {
						if(last) {
							if(SU) {
								values_s = (short unsigned *) last;
								n = *values_s + 1;
								while(--n) {
									if((Score_r[*++values_s] += reps) == reps) {
										bestTemplates_r[++*bestTemplates_r] = *values_s;
									}
								}
							} else {
								n = *last + 1;
								while(--n) {
									if((Score_r[*++last] += reps) == reps) {
										bestTemplates_r[++*bestTemplates_r] = *last;
									}
								}
							}
							hitCounter += reps;
						}
						last = values;
						reps = 1;
					}
				}
			}
			j = end + 1;
		}
		if(last) {
			if(SU) {
				values_s = (short unsigned *) last;
				n = *values_s + 1;
				while(--n) {
					if((Score_r[*++values_s] += reps) == reps) {
						bestTemplates_r[++*bestTemplates_r] = *values_s;
					}
				}
			} else {
				n = *last + 1;
				while(--n) {
					if((Score_r[*++last] += reps) == reps) {
						bestTemplates_r[++*bestTemplates_r] = *last;
					}
				}
			}
			hitCounter += reps;
		}
		reps = 0;
		
		/* get best match(es) */
		/*if(hitCounter * kmersize > (end - hitCounter + kmersize)) {*/
		if(hitCounter) {
			bestScore_r = getMatch(bestTemplates_r, Score_r);
			
			/*
			bestHits = 0;
			for(l = 1; l <= *bestTemplates_r; ++l) {
				if(Score_r[(template = bestTemplates_r[l])] > bestScore_r) {
					bestScore_r = Score_r[template];
					bestHits = 1;
					bestTemplates_r[bestHits] = template;
				} else if(Score_r[template] == bestScore_r) {
					++bestHits;
					bestTemplates_r[bestHits] = template;
				}
				Score_r[template] = 0;
			}
			*bestTemplates_r = bestHits;
			*/
		} else {
			/*
			for(l = *bestTemplates_r; l != 0; --l) {
				Score_r[bestTemplates_r[l]] = 0;
			}
			*/
			*bestTemplates_r = clearScore(bestTemplates_r, Score_r);
		}
	}
	qseq_r->N[0]--;
	
	/* Validate best match */
	i = 1;
	if(bestScore > 0 || bestScore_r > 0) {
		end = qseq->seqlen + 1;
		//if((bestScore >= bestScore_r && bestScore * kmersize > (end - bestScore)) || (bestScore < bestScore_r && bestScore_r * kmersize > (end - bestScore_r))) {
		if(kmersize <= bestScore || kmersize <= bestScore_r) {
			if(bestScore > bestScore_r) {
				lock(excludeOut);
				i = deConPrintPtr(bestTemplates, qseq, bestScore, header, 0, out);
				unlock(excludeOut);
			} else if(bestScore < bestScore_r) {
				lock(excludeOut);
				i = deConPrintPtr(bestTemplates_r, qseq_r, bestScore_r, header, 16, out);
				unlock(excludeOut);
			} else {
				/* merge */
				for(i = 1; i <= *bestTemplates_r; ++i) {
					bestTemplates[0]++;
					bestTemplates[*bestTemplates] = -bestTemplates_r[i];
				}
				lock(excludeOut);
				i = deConPrintPtr(bestTemplates, qseq, -bestScore, header, 0, out);
				unlock(excludeOut);
			}
		}
	}
	return i;
}

int save_kmers_unionPair(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores, CompDNA *qseq, CompDNA *qseq_r, const Qseqs *header, const Qseqs *header_r, int *extendScore, const int exhaustive, volatile int *excludeOut, FILE *out) {
	
	int i, bestScore, bestScore_r, hitCounter, kmersize, flag, flag_r, rev;
	
	/* get_kmers_for_pair, returns a positive number if templates are found.
	zero otherwise */
	kmersize = templates->kmersize;
	if(templates->prefix_len == 0 && templates->prefix != 0) {
		rev = 0;
	} else {
		rev = 1;
	}
	
	/* get forward */
	/*if((hitCounter = get_kmers_for_pair_ptr(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, qseq, extendScore, exhaustive)) && (qseq->seqlen - hitCounter - kmersize) < hitCounter * kmersize) {*/
	if((hitCounter = get_kmers_for_pair_ptr(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, qseq, extendScore, exhaustive))) {
		/* got hits */
		bestScore = getF(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates);
		
		if(kmersize < bestScore && bestScore * kmersize < (qseq->seqlen - bestScore)) {
			bestScore = 0;
		}
	} else {
		bestScore = 0;
		if(hitCounter) {
			i = *bestTemplates + 1;
			while(--i) {
				Score[bestTemplates[i]] = 0;
			}
			*bestTemplates = 0;
			i = *bestTemplates_r + 1;
			while(--i) {
				Score_r[bestTemplates_r[i]] = 0;
			}
			*bestTemplates_r = 0;
		}
	}
	*extendScore = 1;
	
	/* get reverse */
	/*if((hitCounter = get_kmers_for_pair_ptr(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, qseq_r, extendScore, exhaustive)) && (qseq_r->seqlen - hitCounter - kmersize) < hitCounter * kmersize) {*/
	if((hitCounter = get_kmers_for_pair_ptr(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, qseq_r, extendScore, exhaustive))) {
		if(bestScore) {
			bestScore_r = getR(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates);
		} else {
			bestScore_r = getF(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates);
		}
		if(kmersize < bestScore_r && bestScore_r * kmersize < (qseq_r->seqlen - bestScore_r)) {
			bestScore_r = 0;
			*regionTemplates = abs(*regionTemplates);
		}
	} else {
		bestScore_r = 0;
		i = *bestTemplates + 1;
		while(--i) {
			Score[bestTemplates[i]] = 0;
		}
		*bestTemplates = 0;
		i = *bestTemplates_r + 1;
		while(--i) {
			Score_r[bestTemplates_r[i]] = 0;
		}
		*bestTemplates_r = 0;
	}
	
	flag = 65;
	flag_r = 129;
	if(0 < bestScore && 0 < bestScore_r) {
		if(*regionTemplates < 0) {
			flag |= 2;
			flag_r |= 2;
			
			*regionTemplates = -(*regionTemplates);
			if(0 < regionTemplates[1]) {
				if(rev) {
					flag |= 32;
					flag_r |= 16;
					comp_rc(qseq);
				} else {
					flag |= 16;
					flag_r |= 32;
					comp_rc(qseq_r);
				}
				lock(excludeOut);
				i = printPairPtr(regionTemplates, qseq, bestScore, header, qseq_r, bestScore_r, header_r, flag, flag_r, out);
				unlock(excludeOut);
			} else {
				if(rev) {
					flag |= 16;
					flag_r |= 32;
					comp_rc(qseq_r);
				} else {
					flag |= 32;
					flag_r |= 16;
					comp_rc(qseq);
				}
				for(i = *regionTemplates; i != 0; --i) {
					regionTemplates[i] = -regionTemplates[i];
				}
				lock(excludeOut);
				i = printPairPtr(regionTemplates, qseq_r, bestScore_r, header_r, qseq, bestScore, header, flag_r, flag, out);
				unlock(excludeOut);
			}
			if(i) {
				i = 3;
			} else {
				i = 0;
			}
		} else {
			if(0 < regionTemplates[1]) {
				if(rev) {
					comp_rc(qseq);
				}
				if(regionTemplates[*regionTemplates] < 0) {
					bestScore = -bestScore;
				}
			} else {
				if(rev) {
					flag |= 16;
					flag_r |= 32;
				}
				for(i = 1; i <= *regionTemplates; ++i) {
					regionTemplates[i] = -regionTemplates[i];
				}
			}
			if(0 < bestTemplates[1]) {
				if(rev) {
					comp_rc(qseq_r);
				}
				if(bestTemplates[*bestTemplates] < 0) {
					bestScore_r = -bestScore_r;
				}
			} else {
				if(rev) {
					flag |= 32;
					flag_r |= 16;
				}
				for(i = 1; i <= *bestTemplates; ++i) {
					bestTemplates[i] = -bestTemplates[i];
				}
			}
			lock(excludeOut);
			i = deConPrintPtr(regionTemplates, qseq, bestScore, header, flag, out);
			if(deConPrintPtr(bestTemplates, qseq_r, bestScore_r, header_r, flag_r, out)) {
				i += 2;
			}
			unlock(excludeOut);
		}
		return i;
	} else if(bestScore) {
		if(rev) {
			flag |= 8;
			flag |= 32;
		}
		if(0 < regionTemplates[1]) {
			if(rev) {
				comp_rc(qseq);
			}
			if(regionTemplates[*regionTemplates] < 0) {
				bestScore = -bestScore;
			}
		} else {
			if(rev) {
				flag |= 16;
			}
			for(i = 1; i <= *regionTemplates; ++i) {
				regionTemplates[i] = -regionTemplates[i];
			}
		}
		lock(excludeOut);
		i = deConPrintPtr(regionTemplates, qseq, bestScore, header, flag, out);
		unlock(excludeOut);
		if(i == 0) {
			return 2;
		}
	} else if(bestScore_r) {
		if(rev) {
			flag_r |= 8;
			flag_r |= 32;
		}
		if(0 < regionTemplates[1]) {
			if(rev) {
				comp_rc(qseq_r);
			}
			if(regionTemplates[*regionTemplates] < 0) {
				bestScore_r = -bestScore_r;
			}
		} else {
			if(rev) {
				flag_r |= 16;
			}
			for(i = 1; i <= *regionTemplates; ++i) {
				regionTemplates[i] = -regionTemplates[i];
			}
		}
		lock(excludeOut);
		i = deConPrintPtr(regionTemplates, qseq_r, bestScore_r, header_r, flag_r, out);
		unlock(excludeOut);
		if(i == 0) {
			return 1;
		}
	}
	return 3;
}

int save_kmers_penaltyPair(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores, CompDNA *qseq, CompDNA *qseq_r, const Qseqs *header, const Qseqs *header_r, int *extendScore, const int exhaustive, volatile int *excludeOut, FILE *out) {
	
	int i, bestScore, bestScore_r, compScore, hitCounter, hitCounter_r;
	int kmersize, flag, flag_r, rev;
	
	/* get_kmers_for_pair, returns a positive number if templates are found.
	zero otherwise */
	kmersize = templates->kmersize;
	if(templates->prefix_len == 0 && templates->prefix != 0) {
		rev = 0;
	} else {
		rev = 1;
	}
	
	/* get forward */
	if((hitCounter = get_kmers_for_pair_ptr(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, qseq, extendScore, exhaustive))) {
		/* got hits */
		bestScore = getFirstPen(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores);
	} else {
		bestScore = 0;
	}
	*extendScore = 1;
	
	/* get reverse */
	if((hitCounter_r = get_kmers_for_pair_ptr(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, qseq_r, extendScore, exhaustive))) {
		if(0 < bestScore) {
			bestScore_r = getSecondPen(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores, bestScore, rewards->PE);
		} else {
			bestScore_r = getF(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates);
		}
	} else {
		bestScore_r = 0;
	}
	
	flag = 65;
	flag_r = 129;
	if(0 < bestScore && 0 < bestScore_r) {
		if(*regionTemplates < 0) {
			flag |= 2;
			flag_r |= 2;
			compScore = MIN((hitCounter + hitCounter_r), (bestScore + bestScore_r));
			if(kmersize <= compScore || (qseq->seqlen + qseq_r->seqlen - compScore - (kmersize << 1)) < compScore * kmersize) {
				*regionTemplates = -(*regionTemplates);
				if(0 < regionTemplates[1]) {
					if(rev) {
						flag |= 32;
						flag_r |= 16;
						comp_rc(qseq);
					} else {
						flag |= 16;
						flag_r |= 32;
						comp_rc(qseq_r);
					}
					lock(excludeOut);
					i = printPairPtr(regionTemplates, qseq, bestScore, header, qseq_r, bestScore_r, header_r, flag, flag_r, out);
					unlock(excludeOut);
				} else {
					if(rev) {
						flag |= 16;
						flag_r |= 32;
						comp_rc(qseq_r);
					} else {
						flag |= 32;
						flag_r |= 16;
						comp_rc(qseq);
					}
					for(i = *regionTemplates; i != 0; --i) {
						regionTemplates[i] = -regionTemplates[i];
					}
					lock(excludeOut);
					i = printPairPtr(regionTemplates, qseq_r, bestScore_r, header_r, qseq, bestScore, header, flag_r, flag, out);
					unlock(excludeOut);
				}
				if(i) {
					i = 3;
				} else {
					i = 0;
				}
			} else {
				i = 3;
			}
		} else {
			hitCounter = MIN(hitCounter, bestScore);
			hitCounter = kmersize <= hitCounter || (qseq->seqlen - hitCounter - kmersize) < hitCounter * kmersize;
			if(hitCounter) {
				if(0 < regionTemplates[1]) {
					if(rev) {
						comp_rc(qseq);
					}
					if(regionTemplates[*regionTemplates] < 0) {
						bestScore = -bestScore;
					}
				} else {
					if(rev) {
						flag |= 16;
						flag_r |= 32;
					}
					for(i = *regionTemplates; i != 0; --i) {
						regionTemplates[i] = -regionTemplates[i];
					}
				}
			}
			hitCounter_r = MIN(hitCounter_r, bestScore_r);
			hitCounter_r = kmersize <= hitCounter_r || (qseq_r->seqlen - hitCounter_r - kmersize) < hitCounter_r * kmersize;
			if(hitCounter_r) {
				if(0 < bestTemplates[1]) {
					if(rev) {
						comp_rc(qseq_r);
					}
					if(bestTemplates[*bestTemplates] < 0) {
						bestScore_r = -bestScore_r;
					}
				} else {
					if(rev) {
						flag |= 32;
						flag_r |= 16;
					}
					for(i = *bestTemplates; i != 0; --i) {
						bestTemplates[i] = -bestTemplates[i];
					}
				}
			}
			
			if(hitCounter) {
				lock(excludeOut);
				i = deConPrintPtr(regionTemplates, qseq, bestScore, header, flag, out);
				unlock(excludeOut);
			} else {
				i = 1;
			}
			if(hitCounter_r) {
				lock(excludeOut);
				if(deConPrintPtr(bestTemplates, qseq_r, bestScore_r, header_r, flag_r, out)) {
					i += 2;
				}
				unlock(excludeOut);
			} else {
				i += 2;
			}
		}
		return i;
	} else if(0 < bestScore) {
		hitCounter = MIN(hitCounter, bestScore);
		if(kmersize <= hitCounter || (qseq->seqlen - hitCounter - kmersize) < hitCounter * kmersize) {
			if(rev) {
				flag |= 8;
				flag |= 32;
			}
			if(0 < regionTemplates[1]) {
				if(rev) {
					comp_rc(qseq);
				}
				if(regionTemplates[*regionTemplates] < 0) {
					bestScore = -bestScore;
				}
			} else {
				if(rev) {
					flag |= 16;
				}
				for(i = *regionTemplates; i != 0; --i) {
					regionTemplates[i] = -regionTemplates[i];
				}
			}
			lock(excludeOut);
			i = deConPrintPtr(regionTemplates, qseq, bestScore, header, flag, out);
			unlock(excludeOut);
		} else {
			i = 1;
		}
		if(i == 0) {
			return 2;
		}
	} else if(0 < bestScore_r) {
		hitCounter_r = MIN(hitCounter_r, bestScore_r);
		if(kmersize <= hitCounter_r || (qseq_r->seqlen - hitCounter_r - kmersize) < hitCounter_r * kmersize) {
			if(rev) {
				flag_r |= 8;
				flag_r |= 32;
			}
			if(0 < regionTemplates[1]) {
				if(rev) {
					comp_rc(qseq_r);
				}
				if(regionTemplates[*regionTemplates] < 0) {
					bestScore_r = -bestScore_r;
				}
			} else {
				if(rev) {
					flag_r |= 16;
				}
				for(i = 1; i <= *regionTemplates; ++i) {
					regionTemplates[i] = -regionTemplates[i];
				}
			}
			lock(excludeOut);
			i = deConPrintPtr(regionTemplates, qseq_r, bestScore_r, header_r, flag_r, out);
			unlock(excludeOut);
		} else {
			i = 1;
		}
		if(i == 0) {
			return 1;
		}
	}
	return 3;
}

int save_kmers_forcePair(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores, CompDNA *qseq, CompDNA *qseq_r, const Qseqs *header, const Qseqs *header_r, int *extendScore, const int exhaustive, volatile int *excludeOut, FILE *out) {
	
	int i, bestScore, hitCounter, hitCounter_r, kmersize, flag, flag_r, rev;
	
	/* get_kmers_for_pair, returns a positive number if templates are found.
	zero otherwise */
	kmersize = templates->kmersize;
	if(templates->prefix_len == 0 && templates->prefix != 0) {
		rev = 0;
	} else {
		rev = 1;
	}
	
	/* get forward */
	if((hitCounter = get_kmers_for_pair_ptr(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, qseq, extendScore, exhaustive))) {
		/* got hits */
		getFirstForce(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores);
	} else {
		return 1;
	}
	*extendScore = 1;
	
	/* get reverse */
	/*if((hitCounter_r = get_kmers_for_pair_ptr(templates, rewards, bestTemplates_r, bestTemplates, Score_r, Score, qseq_r, extendScore, exhaustive)) && 
		(qseq->seqlen + qseq_r->seqlen - hitCounter - hitCounter_r - (kmersize << 1)) < (hitCounter + hitCounter_r) * kmersize && 
	(bestScore = getSecondForce(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores))) {*/
	if((hitCounter_r = get_kmers_for_pair_ptr(templates, rewards, bestTemplates_r, bestTemplates, Score_r, Score, qseq_r, extendScore, exhaustive)) && 
	(bestScore = getSecondForce(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores))) {
		
		if(kmersize <= bestScore || (qseq->seqlen + qseq_r->seqlen - bestScore) < bestScore * kmersize) {
			flag = 67;
			flag_r = 131;
			
			if(regionTemplates[*regionTemplates] < 0) {
				bestScore = -bestScore;
			}
			if(0 < regionTemplates[1]) {
				if(rev) {
					flag |= 32;
					flag_r |= 16;
					comp_rc(qseq);
				} else {
					flag |= 16;
					flag_r |= 32;
					comp_rc(qseq_r);
				}
				lock(excludeOut);
				i = printPairPtr(regionTemplates, qseq, bestScore, header, qseq_r, bestScore, header_r, flag, flag_r, out);
				unlock(excludeOut);
			} else {
				if(rev) {
					flag |= 16;
					flag_r |= 32;
					comp_rc(qseq_r);
				} else {
					flag |= 32;
					flag_r |= 16;
					comp_rc(qseq);
				}
				for(i = *regionTemplates; i != 0; --i) {
					regionTemplates[i] = -regionTemplates[i];
				}
				lock(excludeOut);
				i = printPairPtr(regionTemplates, qseq_r, bestScore, header_r, qseq, bestScore, header, flag_r, flag, out);
				unlock(excludeOut);
			}
		} else {
			i = 1;
		}
		if(i == 0) {
			return 0;
		}
	} else if(hitCounter || hitCounter_r) {
		i = *bestTemplates + 1;
		while(--i) {
			Score[bestTemplates[i]] = 0;
		}
		*bestTemplates = 0;
		i = *bestTemplates_r + 1;
		while(--i) {
			Score_r[bestTemplates_r[i]] = 0;
		}
		*bestTemplates_r = 0;
	}
	return 3;
}

int save_kmers_HMM(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, CompDNA *qseq_r, Qseqs *headerOrg, int *extendScore, const int exhaustive, volatile int *excludeOut, FILE *out) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	and is the time determining step */
	static int **RegionTemplates, **TmpNs, *Sizes, *template_lengths, minLen = 0;
	static unsigned ***tVF_scores, ***tVR_scores;
	int i, j, k, l, N, n, i_r, j_r, j_u, seqlen, seqend, end, HIT, SU;
	int hitCounter, template, bestScore, bestHits, start, stop, returner;
	int start_cut, end_cut, deCon, DB_size, mPos, hLen, mPosR, hLenR, Ncheck;
	int *regionTemplates;
	unsigned num, shifter, kmersize, mlen, flag, iPos, cPos;
	unsigned *values, *last, *rlast, **VF_scores, **VR_scores;
	short unsigned *values_s;
	long unsigned mask, mmask, kmer, cmer, hmer, kmerR, cmerR, hmerR;
	int *tmpNs, reps, rreps;
	double Ms, Ns, Ms_prev, Ns_prev, HMM_param[8];
	char *include;
	Qseqs *header;
	
	if(qseq == 0 || qseq->seqlen < (kmersize = templates->kmersize)) {
		if(minLen == 0) {
			/* initial allocate */
			template_lengths = bestTemplates_r;
			n = *bestTemplates;
			Sizes = smalloc(n * sizeof(int));
			tVF_scores = smalloc(n * sizeof(int **));
			tVR_scores = smalloc(n * sizeof(int **));
			RegionTemplates = smalloc(n * sizeof(int*));
			TmpNs = smalloc(n * sizeof(int *));
			i = n;
			while(i--) {
				RegionTemplates[i] = smalloc(((templates->DB_size << 1) + 4) * sizeof(int));
				Sizes[i] = 256;
				TmpNs[i] = smalloc(256 * sizeof(int));
				tVF_scores[i] = calloc(256, sizeof(unsigned *));
				tVR_scores[i] = calloc(256, sizeof(unsigned *));
				if(!tVF_scores[i] || !tVR_scores[i]) {
					ERROR();
				}
			}
			
			/* set minLen */
			/*
			minLen = *template_lengths >> 1;
			i = templates->DB_size;
			while(--i) {
				if(minLen > ++*template_lengths >> 1)) {
					minLen = *template_lengths >> 1;
				}
			}
			minLen = exhaustive < minLen ? minLen : exhaustive;
			*/
			minLen = exhaustive;
		}
		return 1;
	} else if((DB_size = templates->DB_size) < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	returner = 1;
	deCon = deConPrintPtr == &deConPrint;
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = templates->flag;
	header = (Qseqs *) headerOrg;
	include = (char *) (extendScore + (templates->DB_size + 1));
	hLen = kmersize;
	kmer = 0;
	cmer = 0;
	kmerR = 0;
	cmerR = 0;
	
	/* calculate HMM parameters */
	HMM_param[0] = log(1 - pow(0.25, kmersize));
	HMM_param[1] = log(pow(0.25, kmersize));
	HMM_param[2] = log(1 - pow(0.25, kmersize - 1) * 0.75);
	HMM_param[3] = log(pow(0.25, kmersize - 1) * 0.75);
	HMM_param[4] = log(1 - 1.0 / kmersize * 0.75 * 0.25);
	HMM_param[5] = log(1.0 / kmersize * 0.75 * 0.25);
	HMM_param[6] = log(0.75);
	HMM_param[7] = log(0.25);
	
	if(Sizes[*Score] < qseq->size) {
		Sizes[*Score] = qseq->size;
		free(tVF_scores[*Score]);
		free(tVR_scores[*Score]);
		free(TmpNs[*Score]);
		TmpNs[*Score] = smalloc(qseq->size * sizeof(int));
		tVF_scores[*Score] = calloc(qseq->size, sizeof(unsigned *));
		tVR_scores[*Score] = calloc(qseq->size, sizeof(unsigned *));
		if(!tVF_scores[*Score] || !tVR_scores[*Score]) {
			ERROR();
		}
	}
	regionTemplates = RegionTemplates[*Score];
	*regionTemplates++ = templates->DB_size;
	*regionTemplates++ = *Score;
	*regionTemplates++ = 0;
	VF_scores = tVF_scores[*Score];
	VR_scores = tVR_scores[*Score];
	tmpNs = TmpNs[*Score];
	
	/* reverse complement qseq */
	rc_comp(qseq, qseq_r);
	seqlen = qseq->seqlen;
	seqend = seqlen - kmersize + 1;
	i = 0;
	i_r = seqlen - kmersize;
	N = 1;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = seqlen;
	while(N <= qseq->N[0] && 0 <= i_r) {
		/* find a seed */
		HIT = 0;
		if(exhaustive) {
			/* init k-mer */
			if(i < seqend) {
				getKmer_macro(kmer, qseq->seq, i, cPos, iPos, (shifter + 2));
				cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
				getKmer_macro(kmerR, qseq_r->seq, (i_r + 1), cPos, iPos, (shifter + 2));
				cmerR = flag ? initCmer(kmerR, &mPosR, &hmerR, &hLenR, shifter + 2, kmersize, mlen, mmask) : kmerR;
				kmerR <<= 2;
			}
			end = qseq->N[N];
			i += kmersize - 1;
			while(i < end && !HIT) {
				/* update k-mer */
				kmer = updateKmer_macro(kmer, qseq->seq, i, mask);
				cmer = flag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
				kmerR = updateKmerR_macro(kmerR, qseq_r->seq, i_r, shifter);
				cmerR = flag ? updateCmerR(cmerR, &mPosR, &hmerR, &hLenR, kmerR, kmersize, mlen, mmask) : kmerR;
				
				/* lookup */
				if(hashMap_get(templates, cmer) || hashMap_get(templates, cmerR)) {
					HIT = 1;
				} else {
					++i;
					--i_r;
				}
			}
			i -= kmersize - 1;
			end -= kmersize - 1;
		} else {
			end = qseq->N[N] - kmersize + 1;
			while(i < end && !HIT) {
				getKmer_macro(kmer, qseq->seq, i, cPos, iPos, shifter);
				cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
				getKmer_macro(kmerR, qseq_r->seq, i_r, cPos, iPos, shifter);
				cmerR = flag ? getCmerR(kmerR, &mPosR, &hLenR, shifter, mlen, mmask) : kmerR;
				if(hashMap_get(templates, cmer) || hashMap_get(templates, cmerR)) {
					HIT = 1;
				} else {
					i += kmersize;
					i_r -= kmersize;
				}
			}
		}
		
		/* evaluate seed */
		if(HIT) {
			/* set scores attr */
			bestScore = 0;
			*bestTemplates = 0;
			hitCounter = 1;
			
			/* save seed */
			VF_scores[i] = hashMap_get(templates, cmer);
			VR_scores[i] = hashMap_get(templates, cmerR);
			
			/* init HMM */
			Ms_prev = HMM_param[7] + HMM_param[2];
			Ns_prev = HMM_param[5] + HMM_param[0];
			Ms = 0;
			Ns = 0;
			
			/* extend backward */
			j = i - 1;
			j_r = i_r + 1;
			n = N - 1;
			Ncheck = (n > 0) ? -1 : qseq->N[n];
			if(0 <= j) {
				getKmer_macro(kmer, qseq->seq, (j + 1), cPos, iPos, (shifter + 2));
				cmer = flag ? initCmerR(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
				kmer <<= 2;
				getKmer_macro(kmerR, qseq_r->seq, j_r, cPos, iPos, (shifter + 2));
				cmerR = flag ? initCmer(kmerR, &mPosR, &hmerR, &hLenR, shifter + 2, kmersize, mlen, mmask) : kmerR;
			}
			
			j_u = j_r + kmersize - 1;
			while(j >= 0) {
				if(j == Ncheck) {
					k = j;
					while(k >= kmersize && k < (j - kmersize)) {
						/* update next N check */
						if(k == Ncheck) {
							j = Ncheck;
							--n;
							Ncheck = (n > 0) ? -1 : qseq->N[n];
						}
						/* Non match seq */
						if(Ns_prev + HMM_param[0] + HMM_param[4] >= Ms_prev + HMM_param[3] + HMM_param[4]) {
							Ns = Ns_prev + HMM_param[0] + HMM_param[4];
						} else {
							Ns = Ms_prev + HMM_param[3] + HMM_param[4];
						}
						
						/* Match seq */
						if(Ms_prev + HMM_param[2] + HMM_param[6] >= Ns_prev + HMM_param[1] + HMM_param[6]) {
							Ms = Ms_prev + HMM_param[2] + HMM_param[6];
						} else {
							Ms = Ns_prev + HMM_param[1] + HMM_param[6];
							break;
						}
						--k;
						Ns_prev = Ns;
						Ms_prev = Ms;
					}
					
					if(k >= kmersize && k < (j - kmersize)) {
						j = k - 1;
						break;
					} else {
						j = k - 1;
						j_r = seqlen - kmersize - k + 1;
						j_u = j_r + kmersize - 1;
					}
					if(j != Ncheck) {
						getKmer_macro(kmer, qseq->seq, (j + 1), cPos, iPos, (shifter + 2));
						cmer = flag ? initCmerR(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
						kmer <<= 2;
						getKmer_macro(kmerR, qseq_r->seq, j_r, cPos, iPos, (shifter + 2));
						cmerR = flag ? initCmer(kmerR, &mPosR, &hmerR, &hLenR, shifter + 2, kmersize, mlen, mmask) : kmerR;
					}
				} else {
					/* update k-mer */
					kmer = updateKmerR_macro(kmer, qseq->seq, j, shifter);
					cmer = flag ? updateCmerR(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
					kmerR = updateKmer_macro(kmerR, qseq_r->seq, j_u, mask);
					cmerR = flag ? updateCmer(cmerR, &mPosR, &hmerR, &hLenR, kmerR, kmersize, mlen, mmask) : kmerR;
					
					/* lookup and save */
					VF_scores[j] = hashMap_get(templates, cmer);
					VR_scores[j] = hashMap_get(templates, cmerR);
					
					/* HMM */
					if(VF_scores[j] || VR_scores[j]) {
						++hitCounter;
						
						/* Non match seq */
						if(Ns_prev + HMM_param[0] + HMM_param[5] >= Ms_prev + HMM_param[3] + HMM_param[5]) {
							Ns = Ns_prev + HMM_param[0] + HMM_param[5];
						} else {
							Ns = Ms_prev + HMM_param[3] + HMM_param[5];
						}
						
						/* Match seq */
						if(Ms_prev + HMM_param[2] + HMM_param[7] >= Ns_prev + HMM_param[1] + HMM_param[7]) {
							Ms = Ms_prev + HMM_param[2] + HMM_param[7];
						} else {
							Ms = Ns_prev + HMM_param[1] + HMM_param[7];
							--j;
							break;
						}
					} else {
						/* Non match seq */
						if(Ns_prev + HMM_param[0] + HMM_param[4] >= Ms_prev + HMM_param[3] + HMM_param[4]) {
							Ns = Ns_prev + HMM_param[0] + HMM_param[4];
						} else {
							Ns = Ms_prev + HMM_param[3] + HMM_param[4];
						}
						
						/* Match seq */
						if(Ms_prev + HMM_param[2] + HMM_param[6] >= Ns_prev + HMM_param[1] + HMM_param[6]) {
							Ms = Ms_prev + HMM_param[2] + HMM_param[6];
						} else {
							Ms = Ns_prev + HMM_param[1] + HMM_param[6];
							--j;
							break;
						}
					}
					--j;
					++j_u;
					++j_r;
				}
				Ns_prev = Ns;
				Ms_prev = Ms;
			}
			start = j + 1;
			
			/* init HMM */
			Ms_prev = HMM_param[7] + HMM_param[2];
			Ns_prev = HMM_param[5] + HMM_param[0];
			
			/* extend forward */
			j = i + 1;
			j_r = i_r - 1;
			Ncheck = qseq->N[N] - kmersize + 1;
			
			/* init k-mer */
			if(j != Ncheck) {
				getKmer_macro(kmer, qseq->seq, j, cPos, iPos, (shifter + 2));
				cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
				getKmer_macro(kmerR, qseq_r->seq, (j_r + 1), cPos, iPos, (shifter + 2));
				cmerR = flag ? initCmerR(kmerR, &mPosR, &hmerR, &hLenR, shifter + 2, kmersize, mlen, mmask) : kmerR;
				kmerR <<= 2;
			}
			j_u = j + kmersize - 1;
			while(j < seqend) {
				if(j == Ncheck) {
					k = j;
					while(k < seqend && k < (j + kmersize)) {
						/* update next N check */
						if(k == Ncheck) {
							j = Ncheck;
							++N;
							Ncheck = (N == qseq->N[0]) ? seqlen : qseq->N[N] - kmersize + 1;
						}
						
						/* Non match seq */
						if(Ns_prev + HMM_param[0] + HMM_param[4] >= Ms_prev + HMM_param[3] + HMM_param[4]) {
							Ns = Ns_prev + HMM_param[0] + HMM_param[4];
						} else {
							Ns = Ms_prev + HMM_param[3] + HMM_param[4];
						}
						
						/* Match seq */
						if(Ms_prev + HMM_param[2] + HMM_param[6] >= Ns_prev + HMM_param[1] + HMM_param[6]) {
							Ms = Ms_prev + HMM_param[2] + HMM_param[6];
						} else {
							Ms = Ns_prev + HMM_param[1] + HMM_param[6];
							break;
						}
						++k;
						Ns_prev = Ns;
						Ms_prev = Ms;
					}
					
					if(k < seqend && k < (j + kmersize)) {
						j = k;
						break;
					} else {
						j = k;
						j_r = seqlen - kmersize - k;
					}
					if(j != Ncheck) {
						getKmer_macro(kmer, qseq->seq, j, cPos, iPos, (shifter + 2));
						cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
						getKmer_macro(kmerR, qseq_r->seq, (j_r + 1), cPos, iPos, (shifter + 2));
						cmerR = flag ? initCmerR(kmerR, &mPosR, &hmerR, &hLenR, shifter + 2, kmersize, mlen, mmask) : kmerR;
						kmerR <<= 2;
					}
					j_u = j + kmersize - 1;
				} else {
					/* update k-mer */
					kmer = updateKmer_macro(kmer, qseq->seq, j_u, mask);
					cmer = flag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
					kmerR = updateKmerR_macro(kmerR, qseq_r->seq, j_r, shifter);
					cmerR = flag ? updateCmerR(cmerR, &mPosR, &hmerR, &hLenR, kmerR, kmersize, mlen, mmask) : kmerR;
					
					/* lookup and save */
					VF_scores[j] = hashMap_get(templates, cmer);
					VR_scores[j] = hashMap_get(templates, cmerR);
					
					/* HMM */
					if(VF_scores[j] || VR_scores[j]) {
						++hitCounter;
						
						/* Non match seq */
						if(Ns_prev + HMM_param[0] + HMM_param[5] >= Ms_prev + HMM_param[3] + HMM_param[5]) {
							Ns = Ns_prev + HMM_param[0] + HMM_param[5];
						} else {
							Ns = Ms_prev + HMM_param[3] + HMM_param[5];
						}
						
						/* Match seq */
						if(Ms_prev + HMM_param[2] + HMM_param[7] >= Ns_prev + HMM_param[1] + HMM_param[7]) {
							Ms = Ms_prev + HMM_param[2] + HMM_param[7];
						} else {
							Ms = Ns_prev + HMM_param[1] + HMM_param[7];
							++j;
							break;
						}
					} else {
						/* Non match seq */
						if(Ns_prev + HMM_param[0] + HMM_param[4] >= Ms_prev + HMM_param[3] + HMM_param[4]) {
							Ns = Ns_prev + HMM_param[0] + HMM_param[4];
						} else {
							Ns = Ms_prev + HMM_param[3] + HMM_param[4];
						}
						
						/* Match seq */
						if(Ms_prev + HMM_param[2] + HMM_param[6] >= Ns_prev + HMM_param[1] + HMM_param[6]) {
							Ms = Ms_prev + HMM_param[2] + HMM_param[6];
						} else {
							Ms = Ns_prev + HMM_param[1] + HMM_param[6];
							++j;
							break;
						}
					}
					++j;
					++j_u;
					--j_r;
				}
				Ns_prev = Ns;
				Ms_prev = Ms;
			}
			stop = j + kmersize - 1;
			
			/* evaluate hit */
			/*if(hitCounter > 0 && (hitCounter * kmersize > (stop - start - hitCounter + kmersize)) 
				&& ((stop - start) > minLen || start == 0 || stop == seqlen)) {*/
			if(hitCounter > 0 && ((stop - start) > minLen || start == 0 || stop == seqlen)) {
				if(deCon) {
					if(SU) {
						for(k = start; k < j; ++k) {
							if(((values_s = (short unsigned *) VF_scores[k]) && values_s[*values_s] == DB_size) 
								|| ((values_s = (short unsigned *) VR_scores[k]) && values_s[*values_s] == DB_size)) {
								--hitCounter;
							}
						}
					} else {
						for(k = start; k < j; ++k) {
							if(((values = VF_scores[k]) && values[*values] == DB_size)
								|| ((values = VR_scores[k]) && values[*values] == DB_size)) {
								--hitCounter;
							}
						}
					}
				}
				
				/* accept hit */
				if(hitCounter > 0) {
					/* gain total scores and mapping templates for this region */
					*bestTemplates = 0;
					*bestTemplates_r = 0;
					
					last = 0;
					reps = 0;
					rlast = 0;
					rreps = 0;
					for(k = start; k < j; ++k) {
						/* forward */
						if(VF_scores[k]) {
							if(VF_scores[k] == last) {
								++reps;
							} else {
								if(last) {
									if(SU) {
										values_s = (short unsigned *) last;
										num = *values_s + 1;
										while(--num) {
											if((Score[*++values_s] += reps) == reps) {
												bestTemplates[++*bestTemplates] = *values_s;
											}
										}
									} else {
										n = *last + 1;
										while(--n) {
											if((Score[*++last] += reps) == reps) {
												bestTemplates[++*bestTemplates] = *last;
											}
										}
									}
								}
								reps = 1;
								last = VF_scores[k];
							}
						}
						
						/* rc */
						if(VR_scores[k]) {
							if(VR_scores[k] == rlast) {
								++rreps;
							} else {
								if(rlast) {
									if(SU) {
										values_s = (short unsigned *) rlast;
										num = *values_s + 1;
										while(--num) {
											if((Score_r[*++values_s] += rreps) == rreps) {
												bestTemplates_r[++*bestTemplates_r] = *values_s;
											}
										}
									} else {
										num = *rlast + 1;
										while(--num) {
											if((Score_r[*++rlast] += rreps) == rreps) {
												bestTemplates_r[++*bestTemplates_r] = *rlast;
											}
										}
									}
								}
								rreps = 1;
								rlast = VR_scores[k];
							}
						}
						
					}
					if(last) {
						if(SU) {
							values_s = (short unsigned *) last;
							num = *values_s + 1;
							while(--num) {
								if((Score[*++values_s] += reps) == reps) {
									bestTemplates[++*bestTemplates] = *values_s;
								}
							}
						} else {
							num = *last + 1;
							while(--num) {
								if((Score[*++last] += reps) == reps) {
									bestTemplates[++*bestTemplates] = *last;
								}
							}
						}
					}
					if(rlast) {
						if(SU) {
							values_s = (short unsigned *) rlast;
							num = *values_s + 1;
							while(--num) {
								if((Score_r[*++values_s] += rreps) == rreps) {
									bestTemplates_r[++*bestTemplates_r] = *values_s;
								}
							}
						} else {
							num = *rlast + 1;
							while(--num) {
								if((Score_r[*++rlast] += rreps) == rreps) {
									bestTemplates_r[++*bestTemplates_r] = *rlast;
								}
							}
						}
					}
					
					/* cut out template hits */
					while(HIT != 0) {
						/* get best score */
						bestScore = 0;
						bestHits = 0;
						/* forward */
						for(k = 1; k <= *bestTemplates; ++k) {
							template = bestTemplates[k];
							if(Score[template] > bestScore) {
								bestScore = Score[template];
								regionTemplates[(bestHits = 1)] = template;
							} else if(Score[template] == bestScore) {
								if(Score[template]) {
									regionTemplates[++bestHits] = template;
								} else {
									bestTemplates[k] = bestTemplates[(*bestTemplates)--];
									--k;
								}
							}
						}
						
						/* rc */
						for(k = 1; k <= *bestTemplates_r; ++k) {
							template = bestTemplates_r[k];
							if(Score_r[template] > bestScore) {
								bestScore = Score_r[template];
								regionTemplates[(bestHits = 1)] = -template;
							} else if(Score_r[template] == bestScore) {
								if(bestScore) {
									regionTemplates[++bestHits] = -template;
								} else {
									bestTemplates_r[k] = bestTemplates_r[(*bestTemplates_r)--];
									--k;
								}
							}
						}
						*regionTemplates = bestHits;
						
						
						if(bestScore > 0) {
							bestHits = *regionTemplates;
							/* find limits of match */
							start_cut = j;
							for(k = 1; k <= bestHits; ++k) {
								template = (regionTemplates[k] > 0) ? regionTemplates[k] : -regionTemplates[k];
								for(l = start; l < start_cut; ++l) {
									if(VR_scores[l] && intpos_bin_contaminationPtr(VR_scores[l], template) != -1) {
										start_cut = l;
									}
									if(VF_scores[l] && intpos_bin_contaminationPtr(VF_scores[l], template) != -1) {
										start_cut = l;
									}
								}
							}
							end_cut = start_cut;
							for(k = 1; k <= bestHits; ++k) {
								template = (regionTemplates[k] > 0) ? regionTemplates[k] : -regionTemplates[k];
								for(l = j; l > end_cut; --l) {
									if(VR_scores[l] && intpos_bin_contaminationPtr(VR_scores[l], template) != -1) {
										end_cut = l;
									}
									if(VF_scores[l] && intpos_bin_contaminationPtr(VF_scores[l], template) != -1) {
										end_cut = l;
									}
								}
							}
							
							/* evaluate best hit */
							if(bestScore * kmersize > (end_cut - start_cut - bestScore + kmersize)) {
								/* check for hits on rc */
								HIT = (regionTemplates[*regionTemplates] > 0) ? 1 : -1;
								/* print */
								if(start != 0 && j != qseq->seqlen) {
									ankerAndClean(regionTemplates, Score, Score_r, include, template_lengths, VF_scores, VR_scores, tmpNs, qseq, HIT, bestScore, start_cut, end_cut, header, excludeOut, out);
								} else {
									ankerPtr(regionTemplates, Score, Score_r, include, template_lengths, VF_scores, VR_scores, tmpNs, qseq, HIT, bestScore, start_cut, end_cut, header, excludeOut, out);
								}
								returner = 0;
							} else {
								/* clear scores */
								for(k = 1; k <= *bestTemplates; ++k) {
									Score[bestTemplates[k]] = 0;
								}
								for(k = 1; k <= *bestTemplates_r; ++k) {
									Score_r[bestTemplates_r[k]] = 0;
								}
								HIT = 0;
							}
						} else {
							/* clear scores */
							for(k = 1; k <= *bestTemplates; ++k) {
								Score[bestTemplates[k]] = 0;
							}
							for(k = 1; k <= *bestTemplates_r; ++k) {
								Score_r[bestTemplates_r[k]] = 0;
							}
							
							HIT = 0;
						}
					}
				}
			}
			
			/* clear scores */
			for(k = start; k < j; ++k) {
				VF_scores[k] = 0;
				VR_scores[k] = 0;
			}
			
			i = stop + 1;
			i_r = seqlen - kmersize - i;
		} else {
			++N;
		}
	}
	
	return returner;
}

void ankerAndClean(int *regionTemplates, int *Score, int *Score_r, char *include, int *template_lengths, unsigned **VF_scores, unsigned **VR_scores, int *tmpNs, CompDNA *qseq, int HIT, int bestScore, int start_cut, int end_cut, Qseqs *header, volatile int *excludeOut, FILE *out) {
	
	static long unsigned *softProxi = 0;
	static double minFrac = 0.0;
	int k, l, bestHits, bestHitsCov, DB_size, end, score, proxiScore;
	int template, *Templates;
	unsigned *values, n, SU;
	short unsigned *values_s;
	double thisCov, bestCov;
	CompDNA tmpQseq;
	
	if(Score_r == 0) {
		minFrac = *((double *) regionTemplates);
		softProxi = (long unsigned *) Score;
		return;
	} else if((DB_size = regionTemplates[-3]) < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	/* mark hits */
	bestHits = *regionTemplates + 1;
	Templates = regionTemplates;
	while(--bestHits) {
		include[abs(*++Templates)] = 1;
	}
	
	/* get best cov */
	bestHits = *regionTemplates;
	Templates = regionTemplates;
	bestHitsCov = template_lengths[abs(*++Templates)];
	while(--bestHits) {
		if(template_lengths[abs(*++Templates)] < bestHitsCov) {
			bestHitsCov = template_lengths[abs(*Templates)];
		}
	}
	
	/* make sure cuts isn't random seeds */
	if(minFrac) {
		proxiScore = minFrac * bestScore;
		bestCov = 1.0 * proxiScore / bestHitsCov;
		bestHits = *regionTemplates;
		end = end_cut - 92;
		for(k = start_cut + 92; k <= end; ++k) {
			if(VF_scores[k]) {
				if(SU) {
					values_s = (short unsigned *) VF_scores[k];
					n = *values_s;
					for(l = 1; l <= n; ++l) {
						if(include[(template = values_s[l])] == 0 && template != DB_size && 
							(proxiScore <= (score = Score[template]) || bestCov * template_lengths[template] <= score)) {
							include[template] = 1;
							regionTemplates[++bestHits] = template;
							if(softProxi) {
								softProxi[template] += score;
							}
						}
						Score[template]--;
					}
				} else {
					values = VF_scores[k];
					n = *values;
					for(l = 1; l <= n; ++l) {
						if(include[(template = values[l])] == 0 && template != DB_size &&
							(proxiScore <= (score = Score[template]) || bestCov * template_lengths[template] <= score)) {
							include[template] = 1;
							regionTemplates[++bestHits] = template;
							if(softProxi) {
								softProxi[template] += score;
							}
						}
						Score[template]--;
					}
				}
				VF_scores[k] = 0;
			}
			if(VR_scores[k]) {
				if(SU) {
					values_s = (short unsigned *) VR_scores[k];
					n = *values_s;
					for(l = 1; l <= n; ++l) {
						if(include[(template = values_s[l])] == 0 && template != DB_size &&
							(proxiScore <= (score = Score_r[template]) || bestCov * template_lengths[template] <= score)) {
							include[template] = 1;
							regionTemplates[bestHits] = -template;
							if(softProxi) {
								softProxi[template] += score;
							}
						}
						Score_r[template]--;
					}
				} else {
					values = VR_scores[k];
					n = *values;
					for(l = 1; l <= n; ++l) {
						if(include[(template = values[l])] == 0 && template != DB_size &&
							(proxiScore <= (score = Score_r[template]) || bestCov * template_lengths[template] <= score)) {
							include[template] = 1;
							regionTemplates[bestHits] = -template;
							if(softProxi) {
								softProxi[template] += score;
							}
						}
						Score_r[template]--;
					}
				}
				VR_scores[k] = 0;
			}
		}
		*regionTemplates = bestHits;
	} else {
		bestCov = 1.0 * bestScore / bestHitsCov;
		bestHits = *regionTemplates;
		end = end_cut - 92;
		for(k = start_cut + 92; k <= end; ++k) {
			if(VF_scores[k]) {
				if(SU) {
					values_s = (short unsigned *) VF_scores[k];
					n = *values_s;
					for(l = 1; l <= n; ++l) {
						template = values_s[l];
						if(include[template] == 0 && template != DB_size) {
							thisCov = 1.0 * Score[template] / template_lengths[template];
							if(thisCov > bestCov) {
								include[template] = 1;
								bestCov = thisCov;
								bestHits = *regionTemplates + 1;
								regionTemplates[bestHits] = template;
							} else if(thisCov == bestCov) {
								include[template] = 1;
								regionTemplates[++bestHits] = template;
							}
						}
						Score[template]--;
					}
				} else {
					values = VF_scores[k];
					n = *values;
					for(l = 1; l <= n; ++l) {
						template = values[l];
						if(include[template] == 0 && template != DB_size) {
							thisCov = 1.0 * Score[template] / template_lengths[template];
							if(thisCov > bestCov) {
								include[template] = 1;
								bestCov = thisCov;
								bestHits = *regionTemplates + 1;
								regionTemplates[bestHits] = template;
							} else if(thisCov == bestCov) {
								include[template] = 1;
								regionTemplates[++bestHits] = template;
							}
						}
						Score[template]--;
					}
				}
				VF_scores[k] = 0;
			}
			if(VR_scores[k]) {
				if(SU) {
					values_s = (short unsigned *) VR_scores[k];
					n = *values_s;
					for(l = 1; l <= n; ++l) {
						template = values_s[l];
						if(include[template] == 0 && template != DB_size) {
							thisCov = 1.0 * Score_r[template] / template_lengths[template];
							if(thisCov > bestCov) {
								include[template] = 1;
								HIT = -1;
								bestCov = thisCov;
								bestHits = *regionTemplates + 1;
								regionTemplates[bestHits] = -template;
							} else if(thisCov == bestCov) {
								include[template] = 1;
								HIT = -1;
								regionTemplates[++bestHits] = -template;
							}
						}
						Score_r[template]--;
					}
				} else {
					values = VR_scores[k];
					n = *values;
					for(l = 1; l <= n; ++l) {
						template = values[l];
						if(include[template] == 0 && template != DB_size) {
							thisCov = 1.0 * Score_r[template] / template_lengths[template];
							if(thisCov > bestCov) {
								include[template] = 1;
								HIT = -1;
								bestCov = thisCov;
								bestHits = *regionTemplates + 1;
								regionTemplates[bestHits] = -template;
							} else if(thisCov == bestCov) {
								include[template] = 1;
								HIT = -1;
								regionTemplates[++bestHits] = -template;
							}
						}
						Score_r[template]--;
					}
				}
				VR_scores[k] = 0;
			}
		}
		*regionTemplates = bestHits;
	}
	
	/* clear nearest templates on both sides of match */
	end = (qseq->seqlen < (start_cut + 92)) ? qseq->seqlen : (start_cut + 92);
	start_cut = ((start_cut - 92) < 0) ? 0 : (start_cut - 92);
	for(k = start_cut; k < end; ++k) {
		if(VF_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VF_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score[values_s[l]]--;
				}
			} else {
				values = VF_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score[values[l]]--;
				}
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VR_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score_r[values_s[l]]--;
				}
			} else {
				values = VR_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score_r[values[l]]--;
				}
			}
			VR_scores[k] = 0;
		}
	}
	end = (end_cut - 92) < 0 ? 0 : (end_cut - 92);
	end_cut = ((end_cut + 92) > qseq->seqlen) ? qseq->seqlen : (end_cut + 92);
	for(k = end_cut; k > end; --k) {
		if(VF_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VF_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score[values_s[l]]--;
				}
			} else {
				values = VF_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score[values[l]]--;
				}
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VR_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score_r[values_s[l]]--;
				}
			} else {
				values = VR_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score_r[values[l]]--;
				}
			}
			VR_scores[k] = 0;
		}
	}
	
	/* unmark hits */
	bestHits = *regionTemplates + 1;
	Templates = regionTemplates;
	while(--bestHits) {
		include[abs(*++Templates)] = 0;
	}
	
	
	/* modify limits of match seq */
	start_cut = ((start_cut - 92) < 0) ? 0 : (start_cut - 92);
	end_cut = ((end_cut + 92) > qseq->seqlen) ? qseq->seqlen : (end_cut + 92);
	start_cut = (start_cut >> 5) << 5;
	end_cut = ((end_cut >> 5) << 5) + 32;
	end_cut = (end_cut < qseq->seqlen) ? end_cut : qseq->seqlen;
	tmpQseq.seqlen = (end_cut - start_cut);
	tmpQseq.seq = qseq->seq + (start_cut >> 5);
	tmpQseq.N = tmpNs;
	
	for(k = 1, l = 0; k < qseq->N[0]; ++k) {
		if(start_cut <= qseq->N[k]) {
			++l;
			tmpQseq.N[l] = qseq->N[k] - start_cut;
			if(tmpQseq.N[l] >= tmpQseq.seqlen) {
				--l;
				k = qseq->N[0];
			}
		}
	}
	
	/* trim trailing gaps */
	--tmpQseq.seqlen;
	while(tmpQseq.N[l] == tmpQseq.seqlen && l != 0) {
		--tmpQseq.seqlen;
		--l;
	}
	++tmpQseq.seqlen;
	tmpQseq.complen = (tmpQseq.seqlen >> 5) + 1;
	tmpQseq.N[0] = l;
	
	/* modify header */
	l = header->len;
	if(header->size <= l + 22) {
		header->size += 32;
		header->seq = realloc(header->seq, header->size);
		if(!header->seq) {
			ERROR();
		}
	}
	header->len += sprintf((char *) header->seq + l - 1, "\t%d\t%d", start_cut, end_cut);
	
	lock(excludeOut);
	deConPrintPtr(regionTemplates, &tmpQseq, HIT * bestScore, header, 0, out);
	unlock(excludeOut);
	
	header->seq[l] = 0;
	header->len = l;
}

void ankerAndClean_MEM(int *regionTemplates, int *Score, int *Score_r, char *include, int *template_lengths, unsigned **VF_scores, unsigned **VR_scores, int *tmpNs, CompDNA *qseq, int HIT, int bestScore, int start_cut, int end_cut, Qseqs *header, volatile int *excludeOut, FILE *out) {
	
	static long unsigned *softProxi = 0;
	static double minFrac = 0.0;
	int k, l, end, DB_size, SU, bestHits, proxiScore, template, *Templates;
	unsigned n, *values;
	short unsigned *values_s;
	CompDNA tmpQseq;
	
	if(Score_r == 0) {
		minFrac = *((double *) regionTemplates);
		softProxi = (long unsigned *) Score;
		return;
	} else if((DB_size = regionTemplates[-3]) < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	if(minFrac) {
		/* mark hits */
		bestHits = *regionTemplates + 1;
		Templates = regionTemplates;
		while(--bestHits) {
			include[abs(*++Templates)] = 1;
		}
		proxiScore = minFrac * bestScore;
		bestHits = *regionTemplates;
		end = end_cut - 92;
		for(k = start_cut + 92; k <= end; ++k) {
			if(VF_scores[k]) {
				if(SU) {
					values_s = (short unsigned *) VF_scores[k];
					n = *values_s;
					for(l = 1; l <= n; ++l) {
						template = values_s[l];
						if(include[template] == 0 && proxiScore <= Score[template] && template != DB_size) {
							include[template] = 1;
							regionTemplates[++bestHits] = template;
							if(softProxi) {
								softProxi[template] += Score[template];
							}
						}
						Score[template]--;
					}
				} else {
					values = VF_scores[k];
					n = *values;
					for(l = 1; l <= n; ++l) {
						template = values[l];
						if(include[template] == 0 && proxiScore <= Score[template] && template != DB_size) {
							include[template] = 1;
							regionTemplates[++bestHits] = template;
							if(softProxi) {
								softProxi[template] += Score[template];
							}
						}
						Score[template]--;
					}
				}
				VF_scores[k] = 0;
			}
			if(VR_scores[k]) {
				if(SU) {
					values_s = (short unsigned *) VR_scores[k];
					n = *values_s;
					for(l = 1; l <= n; ++l) {
						template = values_s[l];
						if(include[template] == 0 && proxiScore <= Score_r[template] && template != DB_size) {
							include[template] = 1;
							regionTemplates[++bestHits] = -template;
							if(softProxi) {
								softProxi[template] += Score_r[template];
							}
						}
						Score_r[template]--;
					}
				} else {
					values = VR_scores[k];
					n = *values;
					for(l = 1; l <= n; ++l) {
						template = values[l];
						if(include[template] == 0 && proxiScore <= Score_r[template] && template != DB_size) {
							include[template] = 1;
							regionTemplates[++bestHits] = -template;
							if(softProxi) {
								softProxi[template] += Score_r[template];
							}
						}
						Score_r[template]--;
					}
				}
				VR_scores[k] = 0;
			}
		}
		*regionTemplates = bestHits;
		
		/* clean up scores */
		end = (qseq->seqlen < (start_cut + 92)) ? qseq->seqlen : (start_cut + 92);
		for(k = ((start_cut - 92) < 0) ? 0 : (start_cut - 92); k < end; ++k) {
			if(VF_scores[k]) {
				if(SU) {
					values_s = (short unsigned *) VF_scores[k];
					l = (*values_s) + 1;
					while(--l) {
						Score[*++values_s]--;
					}
				} else {
					values = VF_scores[k];
					l = (*values) + 1;
					while(--l) {
						Score[*++values]--;
					}
				}
				VF_scores[k] = 0;
			}
			if(VR_scores[k]) {
				if(SU) {
					values_s = (short unsigned *) VR_scores[k];
					l = (*values_s) + 1;
					while(--l) {
						Score_r[*++values_s]--;
					}
				} else {
					values = VR_scores[k];
					l = (*values) + 1;
					while(--l) {
						Score_r[*++values]--;
					}
				}
				VR_scores[k] = 0;
			}
		}
		end = ((end_cut + 92) > qseq->seqlen) ? qseq->seqlen : (end_cut + 92);
		for(k = ((end_cut - 92) < 0) ? 0 : (end_cut - 92); k < end; ++k) {
			if(VF_scores[k]) {
				if(SU) {
					values_s = (short unsigned *) VF_scores[k];
					l = (*values_s) + 1;
					while(--l) {
						Score[*++values_s]--;
					}
				} else {
					values = VF_scores[k];
					l = (*values) + 1;
					while(--l) {
						Score[*++values]--;
					}
				}
				VF_scores[k] = 0;
			}
			if(VR_scores[k]) {
				if(SU) {
					values_s = (short unsigned *) VR_scores[k];
					l = (*values_s) + 1;
					while(--l) {
						Score_r[*++values_s]--;
					}
				} else {
					values = VR_scores[k];
					l = (*values) + 1;
					while(--l) {
						Score_r[*++values]--;
					}
				}
				VR_scores[k] = 0;
			}
		}
		start_cut = ((start_cut - 92) < 0) ? 0 : (start_cut - 92);
		end_cut = ((end_cut + 92) > qseq->seqlen) ? qseq->seqlen : (end_cut + 92);
		/* unmark hits */
		bestHits = *regionTemplates + 1;
		Templates = regionTemplates;
		while(--bestHits) {
			include[abs(*++Templates)] = 0;
		}
	} else {
		/* clean up scores */
		start_cut = ((start_cut - 92) < 0) ? 0 : (start_cut - 92);
		end_cut = ((end_cut + 92) > qseq->seqlen) ? qseq->seqlen : (end_cut + 92);
		for(k = start_cut; k < end_cut; ++k) {
			if(VF_scores[k]) {
				if(SU) {
					values_s = (short unsigned *) VF_scores[k];
					l = (*values_s) + 1;
					while(--l) {
						Score[*++values_s]--;
					}
				} else {
					values = VF_scores[k];
					l = (*values) + 1;
					while(--l) {
						Score[*++values]--;
					}
				}
				VF_scores[k] = 0;
			}
			if(VR_scores[k]) {
				if(SU) {
					values_s = (short unsigned *) VR_scores[k];
					l = (*values_s) + 1;
					while(--l) {
						Score_r[*++values_s]--;
					}
				} else {
					values = VR_scores[k];
					l = (*values) + 1;
					while(--l) {
						Score_r[*++values]--;
					}
				}
				VR_scores[k] = 0;
			}
		}
	}
	
	/* modify limits of match seq */
	start_cut = (start_cut >> 5) << 5;
	end_cut = ((end_cut >> 5) << 5) + 32;
	end_cut = (end_cut < qseq->seqlen) ? end_cut : qseq->seqlen;
	tmpQseq.seqlen = (end_cut - start_cut);
	tmpQseq.seq = qseq->seq + (start_cut >> 5);
	tmpQseq.N = tmpNs;
	
	for(k = 1, l = 0; k < qseq->N[0]; ++k) {
		if(start_cut <= qseq->N[k]) {
			++l;
			tmpQseq.N[l] = qseq->N[k] - start_cut;
			if(tmpQseq.N[l] >= tmpQseq.seqlen) {
				--l;
				k = qseq->N[0];
			}
		}
	}
	
	/* trim trailing gaps */
	--tmpQseq.seqlen;
	while(tmpQseq.N[l] == tmpQseq.seqlen && l != 0) {
		--tmpQseq.seqlen;
		--l;
	}
	++tmpQseq.seqlen;
	tmpQseq.complen = (tmpQseq.seqlen >> 5) + 1;
	tmpQseq.N[0] = l;
	
	/* modify header */
	l = header->len;
	if(header->size <= l + 22) {
		header->size += 32;
		header->seq = realloc(header->seq, header->size);
		if(!header->seq) {
			ERROR();
		}
	}
	header->len += sprintf((char *) header->seq + l - 1, "\t%d\t%d", start_cut, end_cut);
	
	lock(excludeOut);
	deConPrintPtr(regionTemplates, &tmpQseq, HIT * bestScore, header, 0, out);
	unlock(excludeOut);
	
	header->seq[l] = 0;
	header->len = l;
}

int save_kmers_chain(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, CompDNA *qseq_r, Qseqs *header, int *extendScore, const int exhaustive, volatile int *excludeOut, FILE *out) {
	
	/* return 0 for match, 1 otherwise */
	static int *Sizes = 0, *template_lengths = 0, minlen = 0;
	static double coverT = 0.5, mrs = 0.5;
	static KmerAnker **tVF_scores, **tVR_scores;
	static SeqmentTree **tSeqments;
	int i, j, j_u, rc, Wl, W1, U, M, MM, Ms, MMs, Us, W1s, score, gaps, HIT;
	int start, end, pos, shifter, kmersize, mlen, cover, len, template, test;
	int cStart, cStart_r, cPos, iPos, len_len, mPos, hLen, flag, seqend, SU;
	int VF_start, VR_start, *bests;
	unsigned DB_size, hitCounter, hitCounter_r, ties, ties_len;
	unsigned *values, *last;
	short unsigned *values_s;
	long unsigned mask, mmask, kmer, cmer, hmer, *seq;
	char *include;
	double score_len;
	KmerAnker *V_score, *V_scores, *VF_scores, *VR_scores, *tmp_score;
	KmerAnker *best_score, *best_score_r, *best_score_len, *best_score_len_r;
	SeqmentTree *chainSeqments;
	
	if(qseq == 0) {
		if(Sizes) {
			free(Sizes);
			Sizes = 0;
			i = *bestTemplates;
			while(i--) {
				free(tVF_scores[i]);
			}
			free(tVF_scores);
		} else {
			/* set coverT */
			coverT = *((double *)(bestTemplates_r));
			mrs = *((double *)(Score));
			template_lengths = Score_r;
			minlen = exhaustive;
			i = *bestTemplates;
			Sizes = smalloc(i * sizeof(int));
			tVF_scores = smalloc(2 * i * sizeof(KmerAnker *));
			tVR_scores = tVF_scores + i;
			tSeqments = calloc(i, sizeof(SeqmentTree *));
			if(!tSeqments) {
				ERROR();
			}
			while(i--) {
				Sizes[i] = 1024;
				tVF_scores[i] = calloc(2048, sizeof(KmerAnker));
				if(!tVF_scores[i]) {
					ERROR();
				}
				tVR_scores[i] = tVF_scores[i] + 1024;
				tSeqments[i] = initializeSeqmentTree(64);
			}
		}
		return 0;
	} else if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return 1;
	} else if((DB_size = templates->DB_size) < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = templates->flag;
	M = rewards->M;
	MM = rewards->MM;
	U = rewards->U;
	W1 = rewards->W1;
	Wl = rewards->Wl;
	hLen = kmersize;
	
	include = (char *) (extendScore + (templates->DB_size + 1));
	values = 0;
	values_s = 0;
	rc_comp(qseq, qseq_r);
	
	if(Sizes[*Score] < qseq->size) {
		Sizes[*Score] = qseq->size;
		free(tVF_scores[*Score]);
		tVF_scores[*Score] = calloc(2 * qseq->size, sizeof(KmerAnker));
		if(!tVF_scores[*Score]) {
			ERROR();
		} else {
			tVR_scores[*Score] = tVF_scores[*Score] + qseq->size;
		}
	}
	VF_scores = tVF_scores[*Score];
	VR_scores = tVR_scores[*Score];
	chainSeqments = tSeqments[*Score];
	chainSeqments->n = 0;
	
	/* get forward ankers */
	hitCounter = 0;
	V_scores = VF_scores;
	V_scores->start = 0;
	V_scores->end = 0;
	V_scores->values = 0;
	V_scores->descend = 0;
	HIT = exhaustive;
	j = 0;
	seqend = qseq->seqlen - kmersize + 1;
	seq = qseq->seq;
	++*(qseq->N);
	qseq->N[*(qseq->N)] = qseq->seqlen;
	for(i = 1; i <= *(qseq->N) && !HIT; ++i) {
		end = qseq->N[i] - kmersize + 1;
		for(; j < end && !HIT; j += kmersize) {
			getKmer_macro(kmer, seq, j, cPos, iPos, shifter);
			cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
			if(hashMap_get(templates, cmer)) {
				HIT = 1;
			}
		}
		j = qseq->N[i] + 1;
	}
	
	if(HIT) {
		V_score = V_scores;
		Ms = 0;
		MMs = 0;
		Us = 0;
		W1s = 0;
		gaps = 0;
		last = 0;
		j = 0;
		for(i = 1; i <= *(qseq->N) && j < seqend; ++i) {
			/* init k-mer */
			getKmer_macro(kmer, seq, j, cPos, iPos, (shifter + 2));
			cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
			end = qseq->N[i];
			for(j_u = j + kmersize - 1; j_u < end; ++j_u) {
				/* update k-mer */
				kmer = updateKmer_macro(kmer, seq, j_u, mask);
				cmer = flag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
				
				/* lookup */
				if((values = hashMap_get(templates, cmer))) {
					if(values == last) {
						/*
						gaps == 0 -> Match
						gaps == kmersize -> 1 MM
						kmersize < gaps -> several mismatches or indel(s)
						gaps < kmersize -> deletion
						*/
						if(gaps == 0) {
							/* match */
							++Ms;
						} else if(gaps == kmersize) {
							/* snp */
							Ms += kmersize;
							++MMs;
						} else {
							/* close anker */
							if(last) {
								/* update and link between ankers */
								V_score->weight = Ms * M + MMs * MM + Us * U + W1s * W1;
								V_score->end = j - gaps + kmersize;
								V_score->descend = V_score + 1;
								++V_score;
							}
							V_score->start = j;
							V_score->values = values;
							V_score->descend = 0;
							last = values;
							Ms = kmersize;
							MMs = 0;
							Us = 0;
							W1s = 0;
							++hitCounter;
						}
					} else {
						if(last) {
							/* update and link between ankers */
							V_score->weight = Ms * M + MMs * MM + Us * U + W1s * W1;
							V_score->end = j - gaps + kmersize;
							V_score->descend = V_score + 1;
							++V_score;
						}
						V_score->start = j;
						V_score->values = values;
						V_score->descend = 0;
						last = values;
						Ms = kmersize;
						MMs = 0;
						Us = 0;
						W1s = 0;
						++hitCounter;
					}
					gaps = 0;
				} else {
					++gaps;
				}
				++j;
			}
			gaps += (qseq->N[i] + 1 - j); /* gap over N's */
			j = end + 1;
		}
		if(last) {
			/* update anker */
			V_score->weight = Ms * M + MMs * MM + Us * U + W1s * W1;
			V_score->end = qseq->seqlen - gaps;
		}
	}
	
	/* get rc ankers, in forward notation */
	hitCounter_r = 0;
	V_scores = VR_scores;
	V_scores->start = 0;
	V_scores->end = 0;
	V_scores->values = 0;
	V_scores->descend = 0;
	HIT = exhaustive;
	j = 0;
	seq = qseq_r->seq;
	++*(qseq_r->N);
	qseq_r->N[*(qseq_r->N)] = qseq_r->seqlen;
	rc = qseq->seqlen - kmersize;
	for(i = 1; i <= *(qseq_r->N) && !HIT; ++i) {
		end = qseq_r->N[i] - kmersize + 1;
		for(; j < end && !HIT; j += kmersize) {
			getKmer_macro(kmer, seq, j, cPos, iPos, shifter);
			cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
			
			if(hashMap_get(templates, cmer)) {
				HIT = 1;
			}
		}
		j = qseq_r->N[i] + 1;
	}
	
	if(HIT) {
		V_score = V_scores;
		Ms = 0;
		MMs = 0;
		Us = 0;
		W1s = 0;
		gaps = 0;
		last = 0;
		j = 0;
		rc = qseq->seqlen - kmersize;
		for(i = 1; i <= *(qseq->N) && j < seqend; ++i) {
			/* init k-mer */
			getKmer_macro(kmer, seq, (rc + 1), cPos, iPos, (shifter + 2));
			cmer = flag ? initCmerR(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
			kmer <<= 2;
			end = qseq->N[i] - kmersize + 1;
			while(j < end) {
				/* update k-mer */
				kmer = updateKmerR_macro(kmer, seq, rc, shifter);
				cmer = flag ? updateCmerR(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
				
				/* lookup */
				if((values = hashMap_get(templates, cmer))) {
					if(values == last) {
						/*
						gaps == 0 -> Match
						gaps == kmersize -> 1 MM
						kmersize < gaps -> several mismatches or indel(s)
						gaps < kmersize -> deletion
						*/
						if(gaps == 0) {
							/* match */
							++Ms;
						} else if(gaps == kmersize) {
							/* snp */
							Ms += kmersize;
							++MMs;
						} else {
							/* close anker */
							if(last) {
								/* update and link between ankers */
								V_score->weight = Ms * M + MMs * MM + Us * U + W1s * W1;
								V_score->end = j - gaps + kmersize;
								V_score->descend = V_score + 1;
								++V_score;
							}
							V_score->start = j;
							V_score->values = values;
							V_score->descend = 0;
							last = values;
							Ms = kmersize;
							MMs = 0;
							Us = 0;
							W1s = 0;
							++hitCounter_r;
						}
					} else {
						if(last) {
							/* update and link between ankers */
							V_score->weight = Ms * M + MMs * MM + Us * U + W1s * W1;
							V_score->end = j - gaps + kmersize;
							V_score->descend = V_score + 1;
							++V_score;
						}
						V_score->start = j;
						V_score->values = values;
						V_score->descend = 0;
						last = values;
						Ms = kmersize;
						MMs = 0;
						Us = 0;
						W1s = 0;
						++hitCounter_r;
					}
					gaps = 0;
				} else {
					++gaps;
				}
				++j;
				--rc;
			}
			gaps += (qseq->N[i] + 1 - j); /* gap over N's */
			j = qseq->N[i] + 1;
			rc = qseq->seqlen - j;
		}
		if(last) {
			/* update anker */
			V_score->weight = Ms * M + MMs * MM + Us * U + W1s * W1;
			V_score->end = qseq->seqlen - gaps;
		}
	}
	--*(qseq->N);
	--*(qseq_r->N);
	
	/* no matches */
	if(!hitCounter && !hitCounter_r) {
		return 1;
	}
	
	/* make chains */
	V_score = VF_scores;
	V_score->score = 0;
	best_score = 0;
	best_score_r = V_score;
	best_score_len = 0;
	best_score_len_r = V_score;
	HIT = hitCounter + 1;
	bests = bestTemplates;
	*bestTemplates = 0;
	*bestTemplates_r = 0;
	ties = 0;
	ties_len = 0;
	rc = 2;
	while(rc--) {
		if(rc == 0) {
			V_score = VR_scores;
			V_score->score = 0;
			V_score->score_len = 0;
			V_score->len_len = 1;
			HIT = hitCounter_r + 1;
			bests = bestTemplates_r;
			best_score = best_score_r;
			best_score_r = V_score;
			best_score_len = best_score_len_r;
			best_score_len_r = V_score;
		}
		
		*bests = 0;
		while(--HIT) {
			/*
			Save pos of best score, and chain it (0 if local).
			Score is always current / w.r.t. the anker you are on.
			extendScore contains endpos of last hit anker.
			*/
			start = V_score->start;
			end = V_score->end;
			
			/* chain anker */
			V_score->score = 0;
			V_score->score_len = 0;
			V_score->len_len = 1;
			if(SU) {
				values_s = (short unsigned *) V_score->values;
				i = *values_s + 1;
				values_s += i;
			} else {
				values = V_score->values;
				i = *values + 1;
				values += i;
			}
			while(--i) {
				template = SU ? *--values_s : *--values;
				score = Score[template];
				pos = extendScore[template];
				gaps = start - pos;
				
				/* extend chain */
				if(!include[template]) {
					include[template] = 1;
					bests[++*bests] = template;
					if(start) {
						score = W1 + (start - 1) * U;
						score = V_score->weight + MAX(Wl, score);
					} else {
						score = V_score->weight;
					}
				} else {
					if(gaps == -kmersize) {
						/* Perfect match */
						score += V_score->weight - (kmersize - 1) * M;
					} else if(gaps == 0) {
						score += V_score->weight + MM;
					} else if(0 < gaps) {
						/* mismatch or insersion */
						if(gaps <= 2) {
							MMs = gaps;
							Ms = 0;
						} else {
							MMs = gaps / kmersize + (gaps % kmersize ? 1 : 0);
							MMs = MAX(2, MMs);
							Ms = MIN(gaps - MMs, kmersize);
							Ms = MIN(Ms, MMs);
						}
						
						/* evaluate best option */
						if((W1 + (gaps - 1) * U) <= (MMs * MM + Ms * M)) {
							score += V_score->weight + Ms * M + MMs * MM;
						} else {
							score += V_score->weight + (W1 + (gaps - 1) * U);
						}
					} else if(mlen != kmersize) {
						/* snp */
						score += V_score->weight + gaps * M + MM;
					} else {
						/* indel */
						score += V_score->weight + gaps * M - (gaps + 1) * U + W1;
					}
					
					/* verify extension */
					if(score < 0) {
						test = start ? (W1 + (start - 1) * U) : 0;
						test = MAX(test, Wl);
						if(score < test + V_score->weight) {
							score = test + V_score->weight;
						}
					}
				}
				
				/* update Scores */
				if(V_score->score < score) {
					V_score->score = score;
				}
				/* length scores */
				if(qseq->seqlen < (len_len = template_lengths[template])) {
					len_len = qseq->seqlen;
				}
				score_len = score;
				if(V_score->len_len != len_len) {
					score_len /= len_len;
					score_len *= V_score->len_len;
				}
				if(V_score->score_len < score_len || (V_score->score_len == score_len && V_score->score_len < score)) {
					V_score->score_len = score;
					V_score->len_len = len_len;
				}
				/* array scores */
				Score[template] = score;
				extendScore[template] = end;
			}
			
			/* Mark last best length corrected hit */
			score_len = V_score->score;
			if(V_score->len_len != best_score_len_r->len_len) {
				score_len /= V_score->len_len;
				score_len *= best_score_len_r->len_len;
			}
			if(best_score_len_r->score_len < score_len) {
				best_score_len_r = V_score;
				ties_len = 0;
			} else if(best_score_len_r->score_len == score_len) {
				if(best_score_len_r->score_len < V_score->score_len) {
					best_score_len_r = V_score;
					ties_len = 0;
				} else if(best_score_len_r->score_len == V_score->score_len) {
					/* first hit on rc is likely last hit on forward */
					best_score_len_r = V_score;
					++ties_len;
				}
			}
			/* Mark last best hit */
			if(best_score_r->score < V_score->score) {
				best_score_r = V_score;
				ties = 0;
			} else if(best_score_r->score == V_score->score) {
				if(best_score_r->score_len < V_score->score_len) {
					best_score_r = V_score;
					ties = 0;
				} else {
					/* first hit on rc is likely last hit on forward */
					best_score_r = V_score;
					++ties;
				}
			}
			
			++V_score;
		}
		
		/* clear Score arrays */
		clearScore(bests, Score);
		clearScore(bests, extendScore);
		i = *bests + 1;
		while(--i) {
			include[*++bests] = 0;
		}
	}
	
	/* no good hits */
	if(best_score->score < kmersize && best_score_r->score < kmersize) {
		return 1;
	}
	
	/* prune hits */
	VF_start = VF_scores ? VF_scores->start : 0;
	VR_start = VR_scores ? VR_scores->start : 0;
	if(!(VF_scores = pruneAnkers(VF_scores, kmersize))) {
		best_score->score = 0;
	}
	if(!(VR_scores = pruneAnkers(VR_scores, kmersize))) {
		best_score_r->score = 0;
	}
	
	/* get best chain, 
	start/end, score, template(s)
	and zero out ankers */
	*include = SU;
	*bestTemplates = 0;
	*bestTemplates_r = 0;
	if(kmerAnkerScore != ankerScore) {
		ties = ties_len;
		ties_len = 0;
		best_score = best_score_len;
		best_score_len = 0;
		best_score_r = best_score_len_r;
		best_score_len_r = 0;
	}
	
	/* proxi might be off */
	cStart = -1;
	cStart_r = -1;
	if(!best_score->score || !best_score_r->score) {
		if(best_score->score) {
			tmp_score = getChainTemplates(best_score, rewards, template_lengths, qseq->seqlen, kmersize, mlen, bestTemplates, Score, extendScore, include);
			cStart = tmp_score->start;
			start = cStart;
			len = best_score->end - start;
			rc = 1;
		} else {
			tmp_score = getChainTemplates(best_score_r, rewards, template_lengths, qseq->seqlen, kmersize, mlen, bestTemplates_r, Score, extendScore, include);
			cStart_r = tmp_score->start;
			start = cStart_r;
			len = best_score_r->end - start;
			rc = 2;
		}
	} else {
		tmp_score = getChainTemplates(best_score, rewards, template_lengths, qseq->seqlen, kmersize, mlen, bestTemplates, Score, extendScore, include);
		cStart = tmp_score->start;
		tmp_score = getChainTemplates(best_score_r, rewards, template_lengths, qseq->seqlen, kmersize, mlen, bestTemplates_r, Score, extendScore, include);
		cStart_r = tmp_score->start;;
		rc = chooseChain(best_score, best_score_r, cStart, cStart_r, &start, &len);
	}
	
	score = MAX(best_score->score, best_score_r->score);
	if(len < minlen || score < kmersize) {
		*include = 0;
		return 1;
	}
	
	/* get best chains */
	while(best_score || best_score_r) {
		/* get equal ankers inside chain,
		that lies under the threshold */
		if(ties) {
			if(rc & 1) {
				score = best_score->score;
				V_score = best_score;
				while((V_score = getTieAnker(start < VF_start ? VF_start : start, V_score, best_score))) {
					/*
					1. tie anker -> tie_len == best_len
					2. Best match is last on seq -> tie_start < best_start
					1. && 2. -> cover = |tie_start - best_end| / len
					*/
					if((V_score->end - start) < coverT * len) {
						/* not a co-match / overlap is insufficient ->
						no more equal matches */
						V_score = 0;
					} else { /* Anker is equal, update with templates */
						/* Mark current templates to avoid double hitting */
						i = *bestTemplates + 1;
						bests = bestTemplates;
						while(--i) {
							include[*++bests] = 1;
							Score[*bests] = 0;
							extendScore[*bests] = 0;
						}
						/* update best template candidates */
						template = *bests;
						*bests = 0;
						getChainTemplates(V_score, rewards, template_lengths, qseq->seqlen, kmersize, mlen, bests, Score, extendScore, include);
						/* update bests */
						*bestTemplates += *bests;
						*bests = template;
					}
				}
				/* Unmark */
				i = *bestTemplates + 1;
				bests = bestTemplates;
				while(--i) {
					include[*++bests] = 0;
					Score[*bests] = 0;
					extendScore[*bests] = 0;
				}
			}
			if(rc & 2) {
				score = best_score_r->score;
				V_score = best_score_r;
				while((V_score = getTieAnker(start < VR_start ? VR_start : start, V_score, best_score_r))) {
					/*
					1. tie anker -> tie_len == best_len
					2. Best match is last on seq -> tie_start < best_start
					1. && 2. -> cover = |tie_start - best_end| / len
					*/
					if((V_score->end - start) < coverT * len ) {
						/* not a co-match / overlap is insufficient ->
						no more equal matches */
						V_score = 0;
					} else { /* Anker is equal, update with templates */
						/* Mark current templates to avoid double hitting */
						i = *bestTemplates_r + 1;
						bests = bestTemplates_r;
						while(--i) {
							include[*++bests] = 1;
							Score[*bests] = 0;
							extendScore[*bests] = 0;
						}
						/* update best template candidates */
						template = *bests;
						*bests = 0;
						getChainTemplates(V_score, rewards, template_lengths, qseq->seqlen, kmersize, mlen, bests, Score, extendScore, include);
						/* update bests */
						*bestTemplates_r += *bests;
						*bests = template;
					}
				}
				/* Unmark */
				i = *bestTemplates_r + 1;
				bests = bestTemplates_r;
				while(--i) {
					include[*++bests] = 0;
					Score[*bests] = 0;
					extendScore[*bests] = 0;
				}
			}
		}
		
		/* check mapping length, and remove insufficient hits */
		if(rc & 1) {
			if(!mrchain(bestTemplates, template_lengths, qseq->seqlen, len)) {
				rc ^= 1;
			}
		}
		if(rc & 2) {
			if(!mrchain(bestTemplates_r, template_lengths, qseq->seqlen, len)) {
				rc ^= 2;
			}
		}
		
		if(rc) {
			/* use segment trees to mark "used" regions of query.
			Allow for quick check of intersection with chosen chains. 
			And insert anker-boundaries in header*/
			growSeqmentTree(chainSeqments, start, start + len);
			if(rc & 1) {
				//growSeqmentTree(chainSeqments, start, best_score->end);
				insertKmerBound(header, start, start + len);
			} else {
				//growSeqmentTree(chainSeqments, start, best_score_r->end);
				insertKmerBound(header, qseq->seqlen - best_score_r->end, qseq->seqlen - start);
			}
			
			/* print match */
			if(rc & 1) {
				if(rc & 2) {
					j = *bestTemplates;
					for(i = 1; i <= *bestTemplates_r; ++i) {
						bestTemplates[++j] = -bestTemplates_r[i];
					}
					*bestTemplates += *bestTemplates_r;
					best_score->score = -best_score->score;
					best_score_r->score = 0;
					*bestTemplates_r = 0;
				}
				lock(excludeOut);
				i = deConPrintPtr(bestTemplates, qseq, best_score->score, header, 0, out);
				unlock(excludeOut);
				best_score->score = 0;
				*bestTemplates = 0;
			} else {
				lock(excludeOut);
				i = deConPrintPtr(bestTemplates_r, qseq_r, best_score_r->score, header, 0, out);
				unlock(excludeOut);
				best_score_r->score = 0;
				*bestTemplates_r = 0;
			}
			/* remove k-mer bound */
			header->len -= (2 * sizeof(int) + 1);
		}
		
		/* get next match */
		ties = 0;
		rc = 0;
		if(best_score) {
			if(best_score->score && (tmp_score = getChainTemplates(best_score, rewards, template_lengths, qseq->seqlen, kmersize, mlen, bestTemplates, Score, extendScore, include))) {
				cStart = tmp_score->start;
				cover = queSeqmentTree(chainSeqments->root, cStart, best_score->end);
				
				/* verify chain */
				len = best_score->end - cStart;
				if(minlen <= len && cover <= coverT * len && mrs * len <= best_score->score) {
					/* get chain */
					rc = 1;
				} else {
					/* silence anker */
					best_score->score = 0;
				}
			} else {
				best_score->score = 0;
			}
			while(best_score && best_score->score == 0) {
				/* find best anker */
				if((best_score = getBestAnker(&VF_scores, &ties, template_lengths))) {
					if(kmersize < best_score->score && (tmp_score = getChainTemplates(best_score, rewards, template_lengths, qseq->seqlen, kmersize, mlen, bestTemplates, Score, extendScore, include))) {
						/* check chain */
						cStart = tmp_score->start;
						cover = queSeqmentTree(chainSeqments->root, cStart, best_score->end);
						
						/* verify chain */
						len = best_score->end - cStart;
						if(minlen <= len && cover <= coverT * len && mrs * len <= best_score->score) {
							/* get chain */
							rc = 1;
						} else {
							/* silence anker */
							best_score->score = 0;
						}
					} else {
						/* silence anker */
						best_score->score = 0;
					}
				}
			}
		}
		
		if(best_score_r) {
			if(best_score_r->score && (tmp_score = getChainTemplates(best_score_r, rewards, template_lengths, qseq->seqlen, kmersize, mlen, bestTemplates_r, Score, extendScore, include))) {
				cStart_r = tmp_score->start;
				cover = queSeqmentTree(chainSeqments->root, cStart_r, best_score_r->end);
				//cover = queSeqmentTree(chainSeqments->root, qseq->seqlen - best_score_r->end, qseq->seqlen - cStart_r);
				
				/* verify chain */
				len = best_score_r->end - cStart_r;
				if(minlen <= len && cover <= coverT * len && mrs * len <= best_score_r->score) {
					/* get chain */
					rc |= 2;
				} else {
					/* silence anker */
					best_score_r->score = 0;
				}
			} else {
				best_score_r->score = 0;
			}
			while(best_score_r && best_score_r->score == 0) {
				/* find best anker */
				if((best_score_r = getBestAnker(&VR_scores, &ties, template_lengths))) {
					if(kmersize < best_score_r->score && (tmp_score = getChainTemplates(best_score_r, rewards, template_lengths, qseq->seqlen, kmersize, mlen, bestTemplates_r, Score, extendScore, include))) {
						/* check chain */
						cStart_r = tmp_score->start;
						cover = queSeqmentTree(chainSeqments->root, cStart_r, best_score_r->end);
						//cover = queSeqmentTree(chainSeqments->root, qseq->seqlen - best_score_r->end, qseq->seqlen - cStart_r);
						
						/* verify chain */
						len = best_score_r->end - cStart_r;
						if(minlen <= best_score_r->end - cStart_r && cover <= coverT * len && mrs * len <= best_score_r->score) {
							/* get chain */
							rc |= 2;
						} else {
							/* silence anker */
							best_score_r->score = 0;
						}
					} else {
						/* silence anker */
						best_score_r->score = 0;
					}
				}
			}
		}
		
		if(!best_score && !best_score_r) {
			*include = 0;
			return 0;
		} else if(best_score && best_score_r) {
			rc = chooseChain(best_score, best_score_r, cStart, cStart_r, &start, &len);
		} else if(best_score) {
			rc = 1;
			start = cStart;
			len = best_score->end - start;
		} else {
			rc = 2;
			start = cStart_r;
			len = best_score_r->end - start;
		}
	}
	*include = 0;
	
	return 1;
}

int save_kmers_sparse_chain(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, CompDNA *qseq_r, Qseqs *header, int *extendScore, const int exhaustive, volatile int *excludeOut, FILE *out) {
	
	/* return 0 for match, 1 otherwise */
	static int *Sizes = 0, *template_lengths = 0, minlen = 0;
	static double coverT = 0.5, mrs = 0.5;
	static KmerAnker **tVF_scores;
	static SeqmentTree **tSeqments;
	int i, j, j_u, shifter, prefix_shifter, kmersize, DB_size, pos, start;
	int Wl, W1, U, M, MM, Ms, MMs, Us, W1s, end, len, gaps, template, score;
	int len_len, mlen, tflag, cPos, iPos, mPos, hLen, seqend, test, VF_start;
	int *bests;
	unsigned SU, HIT, cover, hitCounter, ties, ties_len, prefix_len, flag;
	unsigned *values, *last;
	short unsigned *values_s;
	long unsigned mask, mmask, prefix, kmer, cmer, hmer, *seq;
	char *include;
	double score_len;
	KmerAnker *V_score, *VF_scores, *tmp_score, *best_score, *best_score_len;
	SeqmentTree *chainSeqments;
	
	if(qseq == 0) {
		if(Sizes) {
			free(Sizes);
			Sizes = 0;
			i = *bestTemplates;
			while(i--) {
				free(tVF_scores[i]);
			}
			free(tVF_scores);
		} else {
			/* set coverT */
			coverT = *((double *)(bestTemplates_r));
			mrs = *((double *)(Score));
			template_lengths = Score_r;
			minlen = exhaustive;
			i = *bestTemplates;
			Sizes = smalloc(i * sizeof(int));
			tVF_scores = smalloc(i * sizeof(KmerAnker *));
			tSeqments = calloc(i, sizeof(SeqmentTree *));
			if(!tSeqments) {
				ERROR();
			}
			while(i--) {
				Sizes[i] = 1024;
				tVF_scores[i] = calloc(1024, sizeof(KmerAnker));
				if(!tVF_scores[i]) {
					ERROR();
				}
				tSeqments[i] = initializeSeqmentTree(64);
			}
		}
		return 0;
	} else if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return 1;
	} else if((DB_size = templates->DB_size) < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	tflag = templates->flag;
	prefix = templates->prefix;
	prefix_len = templates->prefix_len;
	prefix_shifter = 64 - (prefix_len << 1);
	M = rewards->M;
	MM = rewards->MM;
	U = rewards->U;
	W1 = rewards->W1;
	Wl = rewards->Wl;
	include = (char *) (extendScore + (templates->DB_size + 1));
	hitCounter = 0;
	values = 0;
	values_s = 0;
	last = 0;
	start = 0;
	seq = qseq->seq;
	seqend = qseq->seqlen - kmersize + 1;
	hLen = kmersize;
	
	if(Sizes[*Score] < qseq->size) {
		Sizes[*Score] = qseq->size << 1;
		free(tVF_scores[*Score]);
		tVF_scores[*Score] = calloc(qseq->size, sizeof(KmerAnker));
		if(!tVF_scores[*Score]) {
			ERROR();
		}
	}
	VF_scores = tVF_scores[*Score];
	chainSeqments = tSeqments[*Score];
	chainSeqments->n = 0;
	V_score = VF_scores;
	V_score->start = 0;
	V_score->end = 0;
	V_score->values = 0;
	V_score->descend = 0;
	if(prefix_len) {
		flag = 16;
		rc_comp(qseq, qseq_r);
		++*(qseq->N);
		qseq->N[*(qseq->N)] = qseq->seqlen;
		i = 0;
		j = qseq->seqlen - kmersize - prefix_len;
		for(gaps = 1; gaps <= *(qseq->N); ++gaps) {
			V_score->end = i;
			end = qseq->N[gaps] - kmersize - prefix_len + 1;
			while(i < end) {
				getKmer_macro(kmer, seq, i, cPos, iPos, prefix_shifter);
				if(kmer == prefix) {
					/* get kmer */
					getKmer_macro(kmer, seq, (i + prefix_len), cPos, iPos, shifter);
					cmer = tflag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
					if((values = hashMap_get(templates, cmer))) {
						if(values == last) {
							/* update end */
							V_score->end = i;
						} else if(last) {
							/* find mid point between ankers */
							tmp_score = V_score++;
							tmp_score->end = (tmp_score->end + i) >> 1;
							
							/* start new anker */
							V_score->start = tmp_score->end + 1;
							V_score->end = i;
							V_score->values = values;
							
							/* finish last anker */
							tmp_score->end += kmersize + prefix_len;
							tmp_score->weight = (tmp_score->end - tmp_score->start) * M;
							tmp_score->descend = V_score;
							
							/* count hit */
							++hitCounter;
							last = values;
						} else {
							/* make hit anker */
							V_score->start = V_score->end ? ((V_score->end + i) >> 1) : 0;
							V_score->end = i;
							V_score->values = values;
							last = values;
						}
					} else if(last) {
						/* finish last anker */
						tmp_score = V_score++;
						tmp_score->end = ((tmp_score->end + i) >> 1) + kmersize + prefix_len;
						tmp_score->weight = (tmp_score->end - tmp_score->start) * M;
						tmp_score->descend = V_score;
						
						/* start new anker */
						V_score->end = i;
						V_score->values = 0;
						
						++hitCounter;
						last = 0;
					} else {
						V_score->end = i;
					}
				} else if(getKmer(qseq_r->seq, j, prefix_shifter) == prefix) { /* f-hit -> rc-hit -> else */
					/* get kmer */
					getKmer_macro(kmer, qseq_r->seq, (j + prefix_len), cPos, iPos, shifter);
					cmer = tflag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
					if((values = hashMap_get(templates, cmer))) {
						if(values == last) {
							/* update end */
							V_score->end = i;
						} else if(last) {
							/* find mid point between ankers */
							tmp_score = V_score++;
							tmp_score->end = (tmp_score->end + i) >> 1;
							
							/* start new anker */
							V_score->start = tmp_score->end + 1;
							V_score->end = i;
							V_score->values = values;
							
							/* finish last anker */
							tmp_score->end += kmersize + prefix_len;
							tmp_score->weight = (tmp_score->end - tmp_score->start) * M;
							tmp_score->descend = V_score;
							
							/* count hit */
							++hitCounter;
							last = values;
						} else {
							/* make hit anker */
							V_score->start = V_score->end ? ((V_score->end + i) >> 1) : 0;
							V_score->end = i;
							V_score->values = values;
							last = values;
						}
					} else if(last) {
						/* finish last anker */
						tmp_score = V_score++;
						tmp_score->end = ((tmp_score->end + i) >> 1) + kmersize + prefix_len;
						tmp_score->weight = (tmp_score->end - tmp_score->start) * M;
						tmp_score->descend = V_score;
						
						/* start new anker */
						V_score->end = i;
						V_score->values = 0;
						
						++hitCounter;
						last = 0;
					} else {
						V_score->end = i;
					}
				}
				++i;
				--j;	
			}
			
			if(last) {
				/* finish last anker */
				tmp_score = V_score++;
				tmp_score->end = i;
				tmp_score->weight = (tmp_score->end - tmp_score->start) * M;
				tmp_score->descend = V_score;
				
				/* start new anker */
				V_score->values = 0;
				
				++hitCounter;
				last = 0;
			}
			i = qseq->N[gaps] + 1;
			j = qseq->seqlen - kmersize - prefix_len - i;
		}
		
		if(hitCounter) {
			--V_score;
			V_score->descend = 0;
		}
		--*(qseq->N);
		
		/* adjust the kmersize for chaining */
		kmersize += prefix_len - 1;
	} else {
		flag = 0;
		HIT = exhaustive;
		j = 0;
		++*(qseq->N);
		qseq->N[*(qseq->N)] = qseq->seqlen;
		for(i = 1; i <= *(qseq->N) && !HIT; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(; j < end && !HIT; j += kmersize) {
				getKmer_macro(kmer, seq, j, cPos, iPos, shifter);
				cmer = tflag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
				if(hashMap_get(templates, cmer)) {
					HIT = 1;
				}
			}
			j = qseq->N[i] + 1;
		}
		
		if(HIT) {
			Ms = 0;
			MMs = 0;
			Us = 0;
			W1s = 0;
			gaps = 0;
			last = 0;
			j = 0;
			for(i = 1; i <= *(qseq->N) && j < seqend; ++i) {
				/* init k-mer */
				getKmer_macro(kmer, seq, j, cPos, iPos, (shifter + 2));
				cmer = tflag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
				end = qseq->N[i];
				for(j_u = j + kmersize - 1; j_u < end; ++j_u) {
					/* update kmer */
					kmer = updateKmer_macro(kmer, seq, j_u, mask);
					cmer = tflag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
					
					/* lookup */
					if((values = hashMap_get(templates, cmer))) {
						if(values == last) {
							/*
							gaps == 0 -> Match
							gaps == kmersize -> 1 MM
							kmersize < gaps -> several mismatches or indel(s)
							gaps < kmersize -> deletion
							*/
							if(gaps == 0) {
								/* match */
								++Ms;
							} else if(gaps == kmersize) {
								/* snp */
								Ms += kmersize;
								++MMs;
							} else {
								if(last) {
									/* update and link between ankers */
									V_score->weight = Ms * M + MMs * MM + Us * U + W1s * W1;
									V_score->end = j - gaps + kmersize;
									V_score->descend = V_score + 1;
									++V_score;
								}
								V_score->start = j;
								V_score->values = values;
								V_score->descend = 0;
								last = values;
								Ms = kmersize;
								MMs = 0;
								Us = 0;
								W1s = 0;
								++hitCounter;
							}
						} else {
							if(last) {
								/* update and link between ankers */
								V_score->weight = Ms * M + MMs * MM + Us * U + W1s * W1;
								V_score->end = j - gaps + kmersize;
								V_score->descend = V_score + 1;
								++V_score;
							}
							V_score->start = j;
							V_score->values = values;
							V_score->descend = 0;
							last = values;
							Ms = kmersize;
							MMs = 0;
							Us = 0;
							W1s = 0;
							++hitCounter;
						}
						gaps = 0;
					} else {
						++gaps;
					}
					++j;
				}
				gaps += (qseq->N[i] + 1 - j); /* gap over N's */
				j = end + 1;
			}
			if(last) {
				/* update anker */
				V_score->weight = Ms * M + MMs * MM + Us * U + W1s * W1;
				V_score->end = qseq->seqlen - gaps;
			}
		}
		--*(qseq->N);
	}
	
	/* no matches */
	if(!hitCounter) {
		return 1;
	}
	
	/* make chains */
	V_score = VF_scores;
	V_score->score = 0;
	V_score->score_len = 0;
	V_score->len_len = 1;
	best_score = V_score;
	best_score_len = V_score;
	HIT = hitCounter + 1;
	bests = bestTemplates;
	*bests = 0;
	ties = 0;
	ties_len = 0;
	while(--HIT) {
		/*
		Save pos of best score, and chain it (0 if local).
		Score is always current / w.r.t. the anker you are on.
		extendScore contains endpos of last hit anker.
		*/
		start = V_score->start;
		end = V_score->end;
		
		/* chain anker */
		V_score->score = 0;
		V_score->score_len = 0;
		V_score->len_len = 1;
		if(SU) {
			values_s = (short unsigned *) V_score->values;
			i = *values_s + 1;
			values_s += i;
		} else {
			values = V_score->values;
			i = *values + 1;
			values += i;
		}
		while(--i) {
			template = SU ? *--values_s : *--values;
			score = Score[template];
			pos = extendScore[template];
			gaps = start - pos;
			
			/* extend chain */
			if(!include[template]) {
				include[template] = 1;
				bests[++*bests] = template;
				if(start) {
					score = W1 + (start - 1) * U;
					score = V_score->weight + MAX(Wl, score);
				} else {
					score = V_score->weight;
				}
			} else {
				if(gaps == -kmersize) {
					/* Perfect match */
					score += V_score->weight - (kmersize - 1) * M;
				} else if(0 < gaps) {
					/* mismatch or insersion */
					if(gaps <= 2) {
						MMs = gaps;
						Ms = 0;
					} else {
						MMs = gaps / kmersize + (gaps % kmersize ? 1 : 0);
						MMs = MAX(2, MMs);
						Ms = MIN(gaps - MMs, kmersize);
						Ms = MIN(Ms, MMs);
					}
					
					/* evaluate best option */
					if((W1 + (gaps - 1) * U) <= (MMs * MM + Ms * M)) {
						score += V_score->weight + Ms * M + MMs * MM;
					} else {
						score += V_score->weight + (W1 + (gaps - 1) * U);
					}
				} else if(mlen != kmersize) {
					/* snp */
					score += V_score->weight + gaps * M + MM;
				} else {
					/* indel */
					score += V_score->weight + gaps * M - (gaps + 1) * U + W1;
				}
				
				/* verify extension */
				if(score < 0) {
					test = start ? (W1 + (start - 1) * U) : 0;
					test = MAX(test, Wl);
					if(score <= test + V_score->weight) {
						score = test + V_score->weight;
					}
				}
			}
			
			/* update Scores */
			if(V_score->score < score) {
				V_score->score = score;
			}
			/* length scores */
			if(qseq->seqlen < (len_len = template_lengths[template])) {
				len_len = qseq->seqlen;
			}
			score_len = score;
			if(V_score->len_len != len_len) {
				score_len /= len_len;
				score_len *= V_score->len_len;
			}
			if(V_score->score_len < score_len || (V_score->score_len == score_len && V_score->score_len < score)) {
				V_score->score_len = score;
				V_score->len_len = len_len;
			}
			/* array scores */
			Score[template] = score;
			extendScore[template] = end;
		}
		
		/* Mark last best length corrected hit */
		score_len = V_score->score;
		if(V_score->len_len != best_score_len->len_len) {
			score_len /= V_score->len_len;
			score_len *= best_score_len->len_len;
		}
		if(best_score_len->score_len < score_len) {
			best_score_len = V_score;
			ties_len = 0;
		} else if(best_score_len->score_len == score_len) {
			if(best_score_len->score_len < V_score->score_len) {
				best_score_len = V_score;
				ties_len = 0;
			} else if(best_score_len->score_len == V_score->score_len) {
				/* first hit on rc is likely last hit on forward */
				best_score_len = V_score;
				++ties_len;
			}
		}
		/* Mark last best hit */
		if(best_score->score < V_score->score) {
			best_score = V_score;
			ties = 0;
		} else if(best_score->score == V_score->score) {
			if(best_score->score_len < V_score->score_len) {
				best_score = V_score;
				ties = 0;
			} else {
				/* first hit on rc is likely last hit on forward */
				best_score = V_score;
				++ties;
			}
		}
		
		++V_score;
	}
	
	/* clear Score arrays */
	clearScore(bests, Score);
	clearScore(bests, extendScore);
	i = *bests + 1;
	while(--i) {
		include[*++bests] = 0;
	}
	
	/* no good hits */
	if(best_score->score < kmersize) {
		return 1;
	}
	
	/* prune hits */
	VF_start = VF_scores ? VF_scores->start : 0;
	if(!(VF_scores = pruneAnkers(VF_scores, kmersize))) {
		best_score->score = 0;
		hitCounter = 0;
	}
	
	/* get best chain, 
	start/end, score, template(s)
	and zero out ankers */
	*include = SU;
	*bestTemplates = 0;
	tmp_score = getChainTemplates(best_score, rewards, template_lengths, qseq->seqlen, kmersize, mlen, bestTemplates, Score, extendScore, include);
	score = best_score->score;
	start = tmp_score->start;
	len = best_score->end - start;
	
	/* verify best anker */
	if(len < minlen || score < kmersize) {
		*include = 0;
		return 1;
	}
	if(kmerAnkerScore != ankerScore) {
		ties = ties_len;
		ties_len = 0;
	}
	
	/* get best chains */
	while(best_score) {
		/* get equal ankers inside chain,
		that lies under the threshold */
		if(ties) {
			score = best_score->score;
			V_score = best_score;
			while((V_score = getTieAnker(start < VF_start ? VF_start : start, V_score, best_score))) {
				/*
				1. tie anker -> tie_len == best_len
				2. Best match is last on seq -> tie_start < best_start
				1. && 2. -> cover = |tie_start - best_end| / len
				*/
				if((V_score->end - start) <= coverT * len ) {
					/* not a co-match / overlap is insufficient ->
					no more equal matches */
					V_score = 0;
				} else { /* Anker is equal, update with templates */
					/* Mark current templates to avoid double hitting */
					i = *bestTemplates + 1;
					bests = bestTemplates;
					while(--i) {
						include[*++bests] = 1;
						Score[*bests] = 0;
						extendScore[*bests] = 0;
					}
					/* update best template candidates */
					template = *bests;
					*bests = 0;
					getChainTemplates(V_score, rewards, template_lengths, qseq->seqlen, kmersize, mlen, bests, Score, extendScore, include);
					/* update */
					*bestTemplates += *bests;
					*bests = template;
				}
			}
			/* unmark */
			i = *bestTemplates + 1;
			bests = bestTemplates;
			while(--i) {
				include[*++bests] = 0;
				Score[*bests] = 0;
				extendScore[*bests] = 0;
			}
		}
		
		/* check mapping length, and remove insufficient hits */
		if(mrchain(bestTemplates, template_lengths, qseq->seqlen, len)) {
			/* use segment trees to mark "used" regions of query.
			Allow for quick check of intersection with chosen chains. */
			growSeqmentTree(chainSeqments, start, best_score->end);
			
			/* insert anker-boundaries */
			insertKmerBound((Qseqs *) header, start, best_score->end);
			
			/* print match */
			lock(excludeOut);
			i = deConPrintPtr(bestTemplates, qseq, best_score->score, header, flag, out);
			unlock(excludeOut);
			
			/* remove k-mer bound */
			header->len -= (2 * sizeof(int) + 1);
		}
		
		/* get next match */
		ties = 0;
		best_score->score = 0;
		*bestTemplates = 0;
		while(best_score->score == 0) {
			/* find best anker */
			if((best_score = getBestAnker(&VF_scores, &ties, template_lengths))) {
				if(kmersize < best_score->score && (tmp_score = getChainTemplates(best_score, rewards, template_lengths, qseq->seqlen, kmersize, mlen, bestTemplates, Score, extendScore, include))) {
					/* check chain */
					start = tmp_score->start;
					cover = queSeqmentTree(chainSeqments->root, start, best_score->end);
					len = best_score->end - start;
					
					/* verify chain */
					if(len < minlen || coverT * len < cover || best_score->score < mrs * len) {
						/* silence anker */
						best_score->score = 0;
					}
				} else {
					/* silence anker */
					best_score->score = 0;
				}
			} else {
				*include = 0;
				return 0;
			}
		}
	}
	*include = 0;
	
	return 1;
}
