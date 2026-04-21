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

#include <math.h>
#include <stdlib.h>
#include "chain.h"
#include "penalties.h"
#include "pherror.h"
#include "stdstat.h"

int (*chainSeedsPtr)(AlnPoints *, int, int, int, unsigned *) = &chainSeeds;
void (*trimSeedsPtr)(AlnPoints *points, int start) = &trimSeeds;

AlnPoints * seedPoint_init(int size, Penalties *rewards) {
	
	AlnPoints *dest;
	
	dest = smalloc(sizeof(AlnPoints));
	dest->len = 0;
	dest->size = size;
	size *= sizeof(int);
	dest->tStart = smalloc(size);
	dest->tEnd = smalloc(size);
	dest->qStart = smalloc(size);
	dest->qEnd = smalloc(size);
	dest->weight = smalloc(size);
	dest->score = smalloc(size);
	dest->next = smalloc(size);
	dest->rewards = rewards;
	
	return dest;
}

void seedPoint_realloc(AlnPoints *dest, int size) {
	
	dest->size = size;
	size *= sizeof(int);
	dest->tStart = realloc(dest->tStart, size);
	dest->tEnd = realloc(dest->tEnd, size);
	dest->qStart = realloc(dest->qStart, size);
	dest->qEnd = realloc(dest->qEnd, size);
	dest->weight = realloc(dest->weight, size);
	dest->score = realloc(dest->score, size);
	dest->next = realloc(dest->next, size);
	if(!dest->tStart || !dest->tEnd || !dest->qStart || !dest->qEnd || !dest->weight || !dest->score || !dest->next) {
		ERROR();
	}
}

void seedPoint_free(AlnPoints *src) {
	
	free(src->tStart);
	free(src->tEnd);
	free(src->qStart);
	free(src->qEnd);
	free(src->weight);
	free(src->score);
	free(src->next);
	src->rewards = 0;
	free(src);
}

int chainSeeds(AlnPoints *points, int q_len, int t_len, int kmersize, unsigned *mapQ) {
	
	int i, j, nMems, weight, gap, score, bestScore, secondScore, bestPos;
	int tStart, tEnd, qStart, qEnd, tGap, qGap, nMin, W1, U, M, MM, Ms, MMs;
	Penalties *rewards;
	
	rewards = points->rewards;
	W1 = rewards->W1;
	U = rewards->U;
	M = rewards->M;
	MM = rewards->MM;
	nMems = points->len;
	i = nMems;
	bestPos = i - 1;
	bestScore = 0;
	secondScore = 0;
	points->score[i] = 0;
	points->next[i] = 0;
	
	while(i--) {
		weight = points->weight[i] * M;
		points->next[i] = 0;
		tEnd = points->tEnd[i];
		qEnd = points->qEnd[i];
		
		/* get stop score */
		gap = MIN(t_len - tEnd, q_len - qEnd);
		Ms = gap;
		
		/* gap */
		if(--gap) {
			gap *= U;
			gap += W1;
		} else if(gap == 0) {
			gap = W1;
		} else {
			gap = 0;
		}
		
		/* match */
		if(Ms == 2) {
			MMs = 2;
			Ms = 0;
		} else {
			MMs = Ms / kmersize + (Ms % kmersize ? 1 : 0);
			MMs = MAX(2, MMs);
			Ms = MIN(Ms - MMs, kmersize);
			Ms = MIN(Ms, MMs);
		}
		Ms = Ms * M + MMs * MM;
		score = weight + (Ms < gap ? gap : Ms);
		
		/* 128 is the bandwidth */
		nMin = MIN(nMems, i + 128);
		
		/* find best link */
		for(j = i + 1; j < nMin; ++j) {
			/* check compability */
			if(qEnd < points->qStart[j]) {
				if(tEnd < points->tStart[j]) { /* full compatability */
					tGap = points->tStart[j] - tEnd;
					qGap = points->qStart[j] - qEnd;
					
					/* calculate score for this chaining */
					if((gap = abs(tGap - qGap))) {
						--gap;
						gap *= U;
						gap += W1;
					}
					//gap += ((MIN(tGap, qGap)) * pM + weight + points->score[j]);
					Ms = MIN(tGap, qGap);
					if(Ms == 2) {
						MMs = 2;
						Ms = 0;
					} else {
						MMs = Ms / kmersize + (Ms % kmersize ? 1 : 0);
						MMs = MAX(2, MMs);
						Ms = MIN(Ms - MMs, kmersize);
						Ms = MIN(Ms, MMs);
					}
					
					gap += weight + points->score[j] + Ms * M + MMs * MM;
					
					/* check if score is max */
					if(score <= gap) {
						score = gap;
						points->next[i] = j;
					}
				} else if(kmersize <= points->tEnd[j] - tEnd) { /* semi compatability */
					/* calculate score for this chaining */
					if((gap = (points->qStart[j] - qEnd))) {
						--gap;
						gap *= U;
						gap += W1;
					}
					/* cut mem score */	
					gap += (weight + points->score[j] - (points->tStart[j] - tEnd) * M);
					
					/* check if score is max */
					if(score < gap) {
						score = gap;
						points->next[i] = j;
					}
				}
			} else if(kmersize <= points->qEnd[j] - qEnd) {
				/* cut in query is reflected in template */
				tStart = points->tStart[j] + qEnd - points->qStart[j];
				if(tEnd < tStart) { /* full compatability */
					tGap = tStart - tEnd;
					
					/* calculate score for this chaining */
					if((gap = tGap)) {
						--gap;
						gap *= U;
						gap += W1;
					}
					/* cut mem score */
					gap += (weight + points->score[j] - (tStart - tEnd) * M);
					
					/* check if score is max */
					if(score < gap) {
						score = gap;
						points->next[i] = j;
					}
				} /* if the modified tStart was invalid, then the entire mem is lost. */
			}
		}
		/* update seed */
		if(points->next[i]) {
			points->weight[i] += (points->weight[points->next[i]] - kmersize + 1);
		} else {
			points->weight[i] -= (kmersize - 1);
		}
		points->score[i] = score;
		
		/* penalize start */
		tStart = points->tStart[i];
		qStart = points->qStart[i];
		gap = MIN(tStart, qStart);
		Ms = gap;
		
		/* gap */
		if(0 < --gap) {
			gap *= U;
			gap += W1;
		} else if(gap == 0) {
			gap = W1;
		} else {
			gap = 0;
		}
		
		/* match */
		if(Ms == 2) {
			MMs = 2;
			Ms = 0;
		} else {
			MMs = Ms / kmersize + (Ms % kmersize ? 1 : 0);
			MMs = MAX(2, MMs);
			Ms = MIN(Ms - MMs, kmersize);
			Ms = MIN(Ms, MMs);
		}
		Ms = Ms * M + MMs * MM;
		score += Ms < gap ? gap : Ms;
		
		/* update bestScore */
		if(bestScore <= score) {
			if(points->next[i] != bestPos) {
				secondScore = bestScore;
			}
			bestScore = score;
			bestPos = i;
		} else if(secondScore <= score && points->next[i] != bestPos) {
			secondScore = bestScore;
		}
	}
	
	/* calculate mapping quality */
	*mapQ = 0 < bestScore ? (ceil(40 * (1 - 1.0 * secondScore / bestScore) * MIN(1, points->weight[bestPos] / 10.0) * log(bestScore))) : 0;
	points->score[bestPos] = bestScore;
	
	return bestPos;
}

int chainSeeds_circular(AlnPoints *points, int q_len, int t_len, int kmersize, unsigned *mapQ) {
	
	int i, j, nMems, weight, gap, score, bestScore, secondScore, bestPos;
	int tStart, tEnd, qStart, qEnd, tGap, qGap, nMin, W1, U, M, MM, Ms, MMs;
	Penalties *rewards;
	
	rewards = points->rewards;
	W1 = rewards->W1;
	U = rewards->U;
	M = rewards->M;
	MM = rewards->MM;
	nMems = points->len;
	i = nMems;
	bestPos = i - 1;
	bestScore = 0;
	secondScore = 0;
	points->score[i] = 0;
	points->next[i] = 0;
	
	while(i--) {
		weight = points->weight[i] * M;
		points->next[i] = 0;
		tEnd = points->tEnd[i];
		qEnd = points->qEnd[i];
		
		/* get stop score */
		gap = MIN(t_len - tEnd, q_len - qEnd);
		Ms = gap;
		
		/* gap */
		if(--gap) {
			gap *= U;
			gap += W1;
		} else if(gap == 0) {
			gap = W1;
		} else {
			gap = 0;
		}
		
		/* match */
		if(Ms == 2) {
			MMs = 2;
			Ms = 0;
		} else {
			MMs = Ms / kmersize + (Ms % kmersize ? 1 : 0);
			MMs = MAX(2, MMs);
			Ms = MIN(Ms - MMs, kmersize);
			Ms = MIN(Ms, MMs);
		}
		Ms = Ms * M + MMs * MM;
		score = weight + (Ms < gap ? gap : Ms);
		
		/* 128 is the bandwidth */
		nMin = MIN(nMems, i + 128);
		
		/* find best link */
		for(j = i + 1; j < nMin; ++j) {
			/* check compability */
			if(qEnd < points->qStart[j]) {
				tStart = points->tStart[j];
				if(tEnd < tStart) { /* full compatability */
					tGap = tStart - tEnd;
					qGap = points->qStart[j] - qEnd;
					
					/* calculate score for this chaining */
					if((gap = abs(tGap - qGap))) {
						--gap;
						gap *= U;
						gap += W1;
					}
					//gap += ((MIN(tGap, qGap)) * pM + weight + points->score[j]);
					Ms = MIN(tGap, qGap);
					if(Ms == 2) {
						MMs = 2;
						Ms = 0;
					} else {
						MMs = Ms / kmersize + (Ms % kmersize ? 1 : 0);
						MMs = MAX(2, MMs);
						Ms = MIN(Ms - MMs, kmersize);
						Ms = MIN(Ms, MMs);
					}
					
					gap += weight + points->score[j] + Ms * M + MMs * MM;
					
					/* check if score is max */
					if(score <= gap) {
						score = gap;
						points->next[i] = j;
					}
				} else if(kmersize <= points->tEnd[j] - tEnd) { /* semi compatability */
					/* calculate score for this chaining */
					if((gap = (points->qStart[j] - qEnd))) {
						--gap;
						gap *= U;
						gap += W1;
					}
					/* cut mem score*/	
					gap += (weight + points->score[j] - (tStart - tEnd) * M);
					
					/* check if score is max */
					if(score < gap) {
						score = gap;
						points->next[i] = j;
					}
				} else if(points->tEnd[j] < points->tStart[i]) { /* circular joinning, full compability */
					tGap = t_len - tEnd + tStart;
					qGap = points->qStart[j] - qEnd;
					
					/* calculate score for this chaining */
					if((gap = abs(tGap - qGap))) {
						--gap;
						gap *= U;
						gap += W1;
					}
					//gap += ((MIN(tGap, qGap)) * pM + weight + points->score[j]);
					Ms = MIN(tGap, qGap);
					if(Ms == 2) {
						MMs = 2;
						Ms = 0;
					} else {
						MMs = Ms / kmersize + (Ms % kmersize ? 1 : 0);
						MMs = MAX(2, MMs);
						Ms = MIN(Ms - MMs, kmersize);
						Ms = MIN(Ms, MMs);
					}
					
					gap += weight + points->score[j] + Ms * M + MMs * MM;
					
					/* check if score is max */
					if(score < gap) {
						score = gap;
						points->next[i] = j;
					}
				}
			} else if(kmersize <= points->qEnd[j] - qEnd) {
				/* cut in query is reflected in template */
				tStart = points->tStart[j] + qEnd - points->qStart[j];
				if(tEnd < tStart) { /* full compatability */
					tGap = tStart - tEnd;
					
					/* calculate score for this chaining */
					if((gap = tGap)) {
						--gap;
						gap *= U;
						gap += W1;
					}
					/* cut mem score */
					gap += (weight + points->score[j] - (tStart - tEnd) * M);
					
					/* check if score is max */
					if(score <= gap) {
						score = gap;
						points->next[i] = j;
					}
				} else {
					if(t_len < tStart) {
						tStart -= t_len;
					}
					if(tStart != tEnd && points->tEnd[j] < tStart) {
						tGap = t_len - tEnd + tStart;
					
						/* calculate score for this chaining */
						if((gap = tGap)) {
							--gap;
							gap *= U;
							gap += W1;
						}
						gap += (weight + points->score[j] - (tEnd - tStart) * M);
						
						/* check if score is max */
						if(score < gap) {
							score = gap;
							points->next[i] = j;
						}
					}
				}/* if the modified tStart was invalid, then the entire mem is lost. */
			}
		}
		/* update seed */
		if(points->next[i]) {
			points->weight[i] += (points->weight[points->next[i]] - kmersize + 1);
		} else {
			points->weight[i] -= (kmersize - 1);
		}
		points->score[i] = score;
		
		/* penalize start */
		tStart = points->tStart[i];
		qStart = points->qStart[i];
		gap = MIN(tStart, qStart);
		Ms = gap;
		
		/* gap */
		if(0 < --gap) {
			gap *= U;
			gap += W1;
		} else if(gap == 0) {
			gap = W1;
		} else {
			gap = 0;
		}
		
		/* match */
		if(Ms == 2) {
			MMs = 2;
			Ms = 0;
		} else {
			MMs = Ms / kmersize + (Ms % kmersize ? 1 : 0);
			MMs = MAX(2, MMs);
			Ms = MIN(Ms - MMs, kmersize);
			Ms = MIN(Ms, MMs);
		}
		Ms = Ms * M + MMs * MM;
		score += Ms < gap ? gap : Ms;
		
		/* update bestScore */
		if(bestScore <= score) {
			if(points->next[i] != bestPos) {
				secondScore = bestScore;
			}
			bestScore = score;
			bestPos = i;
		} else if(secondScore <= score && points->next[i] != bestPos) {
			secondScore = bestScore;
		}
	}
	
	/* calculate mapping quality */
	*mapQ = 0 < bestScore ? ceil(40 * (1 - 1.0 * secondScore / bestScore) * MIN(1, points->weight[bestPos] / 10.0) * log(bestScore)) : 0;
	points->score[bestPos] = bestScore;
	
	return bestPos;
}

void trimSeeds(AlnPoints *points, int start) {
	
	/* trim the start of each seed */
	static int ts = 0;
	int len;
	
	if(!points) {
		ts = start;
		return;
	} else if(!ts) {
		return;
	}
	
	if(points->qStart[start]) {
		/* iterate seeds on best chain */
		do {
			/* trim seed */
			len = points->qEnd[start] - points->qStart[start];
			if(len < ts) {
				/* ensure at least one nucleotide remains in seed */
				points->tStart[start] += --len;
				points->qStart[start] += len;
			} else {
				points->tStart[start] += ts;
				points->qStart[start] += ts;
			}
		} while((start = points->next[start]));
	} else {
		/* iterate seeds on best chain */
		while((start = points->next[start])) {
			/* trim seed */
			len = points->qEnd[start] - points->qStart[start];
			if(len < ts) {
				/* ensure at least one nucleotide remains in seed */
				points->tStart[start] += --len;
				points->qStart[start] += len;
			} else {
				points->tStart[start] += ts;
				points->qStart[start] += ts;
			}
		}
	}
}

void trimSeedsNoLead(AlnPoints *points, int start) {
	
	/* trim the start of each seed */
	static int ts = 0;
	int len;
	
	if(!points) {
		ts = start;
		return;
	} else if(!ts) {
		return;
	}
	
	/* iterate seeds on best chain */
	while((start = points->next[start])) {
		/* trim seed */
		len = points->qEnd[start] - points->qStart[start];
		if(len < ts) {
			/* ensure at least one nucleotide remains in seed */
			points->tStart[start] += --len;
			points->qStart[start] += len;
		} else {
			points->tStart[start] += ts;
			points->qStart[start] += ts;
		}
	}
}
