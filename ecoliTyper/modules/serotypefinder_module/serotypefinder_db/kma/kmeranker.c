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
#include <stdlib.h>
#include <math.h>
#include "kmeranker.h"
#include "penalties.h"

KmerAnker * (*getChainTemplates)(KmerAnker*, const Penalties*, const int*, const int, const int, const int, int*, int*, int*, char*) = &getBestChainTemplates;
int (*kmerAnkerScore)(KmerAnker*) = &ankerScore;
const int (*testExtension)(const int, const int, const int) = &testExtensionScore;
const int (*proxiTestBest)(const double, const int, const int, const int, const int) = &proxiTestBestScore;
KmerAnker * (*getBestAnker)(KmerAnker**, unsigned*, const int*) = &getBestAnkerScore;
KmerAnker * (*getTieAnker)(int, KmerAnker*, const KmerAnker*) = &getTieAnkerScore;

/* helper functions */
int ankerScore(KmerAnker *src) {
	return src->score;
}

int ankerScoreLen(KmerAnker *src) {
	return src->score_len;
}

const int testExtensionScore(const int q_len, const int t_len, const int best_len) {
	return 1;
}

const int testExtensionScoreLen(const int q_len, const int t_len, const int best_len) {
	return (q_len < t_len ? q_len : t_len) == best_len;
}

const int proxiTestBestScore(const double proxiScore, const int score, const int q_len, const int t_len, const int best_len) {
	return (proxiScore <= score);
}

const int proxiTestBestScoreLen(const double proxiScore, const int score, const int q_len, const int t_len, const int best_len) {
	return (proxiScore / best_len * (q_len < t_len ? q_len : t_len) <= score) || proxiScore <= score;
}

const int mrchain(int *bestTemaples, const int *template_lengths, const int q_len, const int maplen) {
	
	static double mrc = 0;
	int i, n, *bests;
	
	if(mrc) {
		if(q_len < mrc * maplen) {
			bests = bestTemaples;
			i = *bests + 1;
			n = 0;
			while(--i) {
				if(mrc * maplen <= template_lengths[*++bests]) {
					++n;
					*++bestTemaples = *bests;
				}
			}
			bestTemaples -= n;
			return (*bestTemaples = n);
		}
	} else if(!maplen) {
		mrc = *((double *)(bestTemaples));
	}
	
	return 1;
}

KmerAnker * getBestChainTemplates(KmerAnker *src, const Penalties *rewards, const int *template_lengths, const int q_len, const int kmersize, const int mlen, int *bests, int *Score, int *extendScore, char *include) {
	
	/* get set of best templates and silences the chain, except for the initial anker */
	int i, j, template, score, tmpScore, bestScore, Wl, W1, U, M, MM, Ms, MMs;
	int gaps, start, end, pos, SU, nextAnker, target_len;
	unsigned *values;
	short unsigned *values_s;
	KmerAnker *node, *prev;
	
	/* mark template candidates */
	nextAnker = 0;
	if(!src) {
		return 0;
	} else if((SU = *include)) {
		values_s = (short unsigned *) src->values;
		i = *values_s;
		*bests = i;
		++i;
		values_s += i;
		bests += i;
		while(--i) {
			//if((include[(*--bests = *--values_s)] ^= 1)) {
			if(++include[(*--bests = *--values_s)] == 1) {
				nextAnker = 1;
			}
		}
	} else {
		values = src->values;
		i = *values;
		*bests = i;
		++i;
		values += i;
		bests += i;
		while(--i) {
			//if((include[(*--bests = *--values)] ^= 1)) {
			if(++include[(*--bests = *--values)] == 1) {
				nextAnker = 1;
			}
		}
	}
	--bests;
	M = rewards->M;
	MM = rewards->MM;
	U = rewards->U;
	W1 = rewards->W1;
	Wl = rewards->Wl;
	bestScore = kmerAnkerScore(src);
	values_s = 0;
	values = 0;
	prev = src;
	target_len = src->len_len;
	
	/* get chaining scores */
	for(node = src; nextAnker; --node) {
		if(SU) {
			values_s = (short unsigned *) node->values;
			i = *values_s + 1;
			values_s += i;
		} else {
			values = node->values;
			i = *values + 1;
			values += i;
		}
		start = node->start;
		end = node->end;
		
		while(--i) {
			template = SU ? *--values_s : *--values;
			if(include[template]) {
				score = Score[template];
				pos = extendScore[template];
				gaps = pos - end;
				
				/* extend chain */
				if(pos == 0) {
					score = node->weight;
				} else {
					if(gaps == -kmersize) {
						/* Perfect match */
						score += node->weight - (kmersize - 1) * M;
					} else if(gaps == 0) {
						score += node->weight + MM;
					} else if(0 < gaps) {
						/* mismatch or insersion */
						if(gaps <= 2) {
							MMs = gaps;
							Ms = 0;
						} else {
							MMs = gaps / kmersize + ((gaps % kmersize) ? 1 : 0);
							MMs = 2 < MMs ? MMs : 2;
							Ms = (gaps - MMs) < kmersize ? (gaps - MMs) : kmersize;
							Ms = Ms < MMs ? Ms : MMs;
						}
						
						/* evaluate best option */
						if((W1 + (gaps - 1) * U) <= (MMs * MM + Ms * M)) {
							score += node->weight + Ms * M + MMs * MM;
						} else {
							score += node->weight + (W1 + (gaps - 1) * U);
						}
					} else if(mlen != kmersize) {
						/* snp */
						score += node->weight + gaps * M + MM;
					} else {
						/* indel */
						score += node->weight + gaps * M - (gaps + 1) * U + W1;
					}
					
					/* mark as used */
					node->score = 0;
				}
				
				/* verify extension */
				if(bestScore <= score) {
					/* test if anker is finalising the chain */
					if(node->start) {
						tmpScore = W1 + (node->start - 1) * U;
						tmpScore = score + (Wl < tmpScore ? tmpScore : Wl);
					} else {
						tmpScore = score;
					}
					if(tmpScore == bestScore && testExtension(q_len, template_lengths[template], target_len)) {
						score = bestScore;
						nextAnker = 0;
						prev = node;
					}
				}
				
				extendScore[template] = start;
				Score[template] = score;
			}
		}
	}
	
	/* get best templates */
	j = 0;
	for(i = 1; i <= *bests; ++i) {
		template = bests[i];
		if(include[template] == 1 && proxiTestBest(bestScore, Score[template], q_len, template_lengths[template], target_len)) {
			bests[++j] = template;
		}
		
		/* clear Score arrays */
		Score[template] = 0;
		include[template] = 0;
		extendScore[template] = 0;
	}
	*bests = j;
	
	return j ? prev : 0;
}

KmerAnker * getProxiChainTemplates(KmerAnker *src, const Penalties *rewards, const int *template_lengths, const int q_len, const int kmersize, const int mlen, int *bests, int *Score, int *extendScore, char *include) {
	
	/* get set of best templates and silences the chain, except for the initial anker */
	static long unsigned *softProxi = 0;
	static double minFrac = 0.0;
	int i, j, template, score, tmpScore, bestScore, Wl, W1, U, M, MM, Ms, MMs;
	int gaps, start, end, pos, SU, nextAnker, target_len;
	unsigned *values;
	short unsigned *values_s;
	double proxiScore;
	KmerAnker *node, *prev;
	
	/* mark template candidates */
	if(!src) {
		minFrac = rewards ? *((double *) rewards) : 0;
		softProxi = template_lengths ? (long unsigned *) template_lengths : 0;
		return 0;
	}
	
	SU = *include;
	*bests = 0;
	M = rewards->M;
	MM = rewards->MM;
	U = rewards->U;
	W1 = rewards->W1;
	Wl = rewards->Wl;
	bestScore = kmerAnkerScore(src);
	proxiScore = minFrac * bestScore;
	values_s = 0;
	values = 0;
	nextAnker = 1;
	prev = src;
	target_len = src->len_len;
	
	/* get chaining scores */
	for(node = src; nextAnker; --node) {
		if(SU) {
			values_s = (short unsigned *) node->values;
			i = *values_s + 1;
			values_s += i;
		} else {
			values = node->values;
			i = *values + 1;
			values += i;
		}
		start = node->start;
		end = node->end;
		
		while(--i) {
			template = SU ? *--values_s : *--values;
			score = Score[template];
			pos = extendScore[template];
			gaps = pos - end;
			
			/* extend chain */
			if(pos == 0) {
				/* new template */
				score = node->weight;
				bests[++*bests] = template;
			} else {
				if(gaps == -kmersize) {
					/* Perfect match */
					score += node->weight - (kmersize - 1) * M;
				} else if(gaps == 0) {
					score += node->weight + MM;
				} else if(0 < gaps) {
					/* mismatch or insersion */
					if(gaps <= 2) {
						MMs = gaps;
						Ms = 0;
					} else {
						MMs = gaps / kmersize + ((gaps % kmersize) ? 1 : 0);
						MMs = 2 < MMs ? MMs : 2;
						Ms = (gaps - MMs) < kmersize ? (gaps - MMs) : kmersize;
						Ms = Ms < MMs ? Ms : MMs;
					}
					
					/* evaluate best option */
					if((W1 + (gaps - 1) * U) <= (MMs * MM + Ms * M)) {
						score += node->weight + Ms * M + MMs * MM;
					} else {
						score += node->weight + (W1 + (gaps - 1) * U);
					}
				} else if(mlen != kmersize) {
					/* snp */
					score += node->weight + gaps * M + MM;
				} else {
					/* indel */
					score += node->weight + gaps * M - (gaps + 1) * U + W1;
				}
				
				/* mark as used */
				node->score = 0;
			}
			
			/* verify extension */
			if(bestScore <= score) {
				/* test if anker is finalising the chain */
				if(node->start) {
					tmpScore = W1 + (node->start - 1) * U;
					tmpScore = score + (Wl < tmpScore ? tmpScore : Wl);
				} else {
					tmpScore = score;
				}
				if(tmpScore == bestScore && testExtension(q_len, template_lengths[template], target_len)) {
					score = bestScore;
					nextAnker = 0;
					prev = node;
				}
			}
			
			extendScore[template] = start;
			Score[template] = score;
		}
	}
	
	/* get best templates */
	j = 0;
	for(i = 1; i <= *bests; ++i) {
		template = bests[i];
		if(!include[template] && proxiTestBest(proxiScore, Score[template], q_len, template_lengths[template], target_len)) {
			bests[++j] = template;
			if(softProxi) {
				softProxi[template] += Score[template];
			}
		}
		
		/* clear Score arrays */
		Score[template] = 0;
		extendScore[template] = 0;
		include[template] = 0;
	}
	*bests = j;
	
	return prev;
}

KmerAnker * pruneAnkers(KmerAnker *V_score, int kmersize) {
	
	KmerAnker *node, *prev;
	
	if(!V_score) {
		return 0;
	}
	while(V_score->score < kmersize && (V_score = V_score->descend));
	
	if(!V_score) {
		return 0;
	}
	
	prev = V_score;
	node = prev;
	while((node = node->descend)) {
		if(kmersize <= node->score) {
			prev->descend = node;
			prev = node;
		}
	}
	prev->descend = 0;
	
	return V_score;
}

KmerAnker * getBestAnkerScore(KmerAnker **src, unsigned *ties, const int *template_lengths) {
	
	KmerAnker *node, *prev, *best;
	
	*ties = 0;
	prev = *src;
	while(prev && prev->score == 0) {
		prev = prev->descend;
	}
	*src = prev;
	if(!prev) {
		return 0;
	}
	best = prev;
	node = prev->descend;
	while(node) {
		if(node->score) {
			if(best->score < node->score) {
				best = node;
				*ties = 0;
			} else if(best->score == node->score) {
				best = node;
				++*ties;
			}
			prev->descend = node;
			prev = node;
		}
		node = node->descend;
	}
	prev->descend = 0;
			
	return best;
}

KmerAnker * getBestAnkerScoreLen(KmerAnker **src, unsigned *ties, const int *template_lengths) {
	
	double score_len;
	KmerAnker *node, *prev, *best;
	
	*ties = 0;
	prev = *src;
	while(prev && prev->score == 0) {
		prev = prev->descend;
	}
	*src = prev;
	if(!prev) {
		return 0;
	}
	best = prev;
	node = prev->descend;
	while(node) {
		if(node->score) {
			
			score_len = node->score_len;
			if(node->len_len != best->len_len) {
				score_len /= node->len_len;
				score_len *= best->len_len;
			}
			
			if(best->score_len < score_len) {
				best = node;
				*ties = 0;
			} else if(best->score_len == score_len) {
				if(best->score_len < node->score_len) {
					best = node;
					*ties = 0;
				} else if(best->score_len == node->score_len) {
					best = node;
					++*ties;
				}
			}
			
			prev->descend = node;
			prev = node;
		}
		node = node->descend;
	}
	prev->descend = 0;
			
	return best;
}

KmerAnker * getTieAnkerScore(int stop, KmerAnker *src, const KmerAnker *bestScore) {
	
	if(!src || src->start <= stop) {
		return 0;
	}
	
	/* search downwards */
	while(stop < (--src)->start) {
		if(src->score == bestScore->score) {
			return src;
		}
	}
	
	return 0;
}

KmerAnker * getTieAnkerScoreLen(int stop, KmerAnker *src, const KmerAnker *bestScore) {
	
	if(!src || src->start <= stop) {
		return 0;
	}
	
	/* search downwards */
	while(stop < (--src)->start) {
		if(src->score_len == bestScore->score_len && src->len_len == bestScore->len_len) {
			return src;
		}
	}
	
	return 0;
}

int chooseChain(const KmerAnker *best_score, const KmerAnker *best_score_r, int cStart, int cStart_r, int *Start, int *Len) {
	
	static double coverT = 0.5, proxi = 1.0;
	int rc, start, end;
	
	if(!best_score) {
		coverT = *((double *)(Start));
		proxi = *((double *)(Len));
		proxi = proxi < 0 ? -proxi : proxi;
		return 0;
	}
	
	/* check scores */
	if(proxi == 1.0) {
		rc = best_score_r->score < best_score->score ? 1 : best_score->score < best_score_r->score ? 2 : 3;
	} else if(best_score_r->score <= best_score->score) {
		/* forward best, check proxi reverse */
		rc = (proxi * best_score->score <= best_score_r->score) ? 3 : 1;
	} else {
		/* reverse best, check proxi forward */
		rc = (proxi * best_score_r->score <= best_score->score) ? 3 : 2;
	}
	
	/* get coordinates */
	if(rc == 1) {
		//tmp_score = getChainTemplates(best_score, rewards, template_lengths, qseq->seqlen, kmersize, bestTemplates, Score, extendScore, include);
		start = cStart;
		end = best_score->end;
	} else if(rc == 2) {
		//tmp_score = getChainTemplates(best_score_r, rewards, template_lengths, qseq->seqlen, kmersize, bestTemplates_r, Score, extendScore, include);
		start = cStart_r;
		end = best_score_r->end;
	} else {
		/* check overlap */
		/*
		tmp_score = getChainTemplates(best_score, rewards, template_lengths, qseq->seqlen, kmersize, bestTemplates, Score, extendScore, include);
		cStart = tmp_score->start;
		tmp_score = getChainTemplates(best_score_r, rewards, template_lengths, qseq->seqlen, kmersize, bestTemplates_r, Score, extendScore, include);
		cStart_r = tmp_score->start;;
		*/
		if(best_score->end < cStart_r) { /* no overlaps */
			start = cStart;
			end = best_score->end;
			rc = 1;
		} else if(best_score_r->end < cStart) {
			start = cStart_r;
			end = best_score_r->end;
			rc = 2;
		} else if(cStart <= cStart_r && best_score_r->end <= best_score->end) { /* contained or complete overlaps */
			start = cStart;
			end = best_score->end;
		} else if(cStart_r <= cStart && best_score->end <= best_score_r->end) {
			start = cStart_r;
			end = best_score_r->end;
		} else if(best_score_r->end < best_score->end) { /* check partial overlaps */
			start = best_score->end - cStart;
			end = best_score_r->end - cStart_r;
			end = start < end ? start : end;
			start = cStart_r;
			if(coverT * end <= best_score_r->end - cStart) {
				end = best_score->end;
			} else {
				end = best_score_r->end;
				rc = 2;
			}
		} else {
			start = best_score->end - cStart;
			end = best_score_r->end - cStart_r;
			end = start < end ? start : end;
			start = cStart;
			if(coverT * end <= best_score->end - cStart_r) {
				end = best_score_r->end;
			} else {
				end = best_score->end;
				rc = 1;
			}
		}
	}
	
	*Start = start;
	*Len = end - start;
	return rc;
}

