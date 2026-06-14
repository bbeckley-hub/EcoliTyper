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
#include <string.h>
#include "align.h"
#include "chain.h"
#include "compdna.h"
#include "hashmapcci.h"
#include "nw.h"
#include "stdnuc.h"
#include "stdstat.h"

AlnScore (*leadTailAlnPtr)(Aln *, Aln *, const long unsigned*, const unsigned char*, int, int, int, const int, NWmat *) = &leadTailAln;
void (*trailTailAlnPtr)(Aln *, Aln *, AlnScore *, const long unsigned *, const unsigned char *, int, int, int, int, const int, NWmat *) = &trailTailAln;

AlnScore skipLeadAln(Aln *aligned, Aln *Frag_align, const long unsigned *tseq, const unsigned char *qseq, int t_e, int t_len, int q_e, const int bandwidth, NWmat *matrices) {
	
	AlnScore Stat;
	
	/* initialize */
	Stat.len = 0;
	Stat.score = 0;
	Stat.match = 0;
	Stat.tGaps = 0;
	Stat.qGaps = 0;
	Stat.pos = t_e;
	
	if(aligned) {
		aligned->start = q_e;
	}
	
	return Stat;
}

AlnScore leadTailAln(Aln *aligned, Aln *Frag_align, const long unsigned *tseq, const unsigned char *qseq, int t_e, int t_len, int q_e, const int bandwidth, NWmat *matrices) {
	
	int bias, band, t_s, q_s;
	AlnScore Stat, NWstat;
	
	/* initialize */
	Stat.len = 0;
	Stat.score = 0;
	Stat.match = 0;
	Stat.tGaps = 0;
	Stat.qGaps = 0;
	Stat.pos = t_e;
	
	if(q_e) {
		/* get boundaries */
		t_s = 0;
		q_s = 0;
		if((q_e << 1) < t_e || (q_e + bandwidth) < t_e) { // big leading template gap, cut down
			//t_s = t_e - MIN(bandwidth, (q_e << 1));
			t_s = t_e - (q_e + (q_e < bandwidth ? q_e : bandwidth));
		} else if((t_e << 1) < q_e || (t_e + bandwidth) < q_e) { // big leading query gap, cut down
			//q_s = q_e - MIN(bandwidth, (t_e << 1));
			q_s = q_e - (t_e + (t_e < bandwidth ? t_e : bandwidth));
		}
		
		/* align */
		if(t_e - t_s > 0 && q_e - q_s > 0) {
			band = abs(t_e - t_s - q_e + q_s) + bandwidth;
			if(q_e - q_s <= band || t_e - t_s <= band) {// || abs(t_e - t_s - q_e - q_s) >= 32) {
				if(Frag_align) {
					NWstat = NW(tseq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, Frag_align, matrices, t_len);
				} else {
					NWstat = NW_score(tseq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, matrices, t_len);
				}
			} else if(Frag_align) {
				NWstat = NW_band(tseq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, Frag_align, band, matrices, t_len);
				//NWstat = NW(tseq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, Frag_align, matrices, t_len);
			} else {
				NWstat = NW_band_score(tseq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, band, matrices, t_len);
				//NWstat = NW_score(tseq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, matrices, t_len);
			}
			
			if(Frag_align) {
				/* trim leading gaps */
				bias = 0;
				if(t_s == 0) {
					while(bias < NWstat.len && (Frag_align->t[bias] == 5 || Frag_align->q[bias] == 5)) {
						if(Frag_align->t[bias] == 5) {
							--NWstat.tGaps;
							++(Frag_align->start);
						} else {
							--NWstat.qGaps;
						}
						++bias;
					}
					NWstat.len -= bias;
					/*if(bias) {
						NWstat.score -= (W1 + (bias - 1) * U);
					}*/
				}
				
				memcpy(aligned->t, Frag_align->t + bias, NWstat.len);
				memcpy(aligned->s, Frag_align->s + bias, NWstat.len);
				memcpy(aligned->q, Frag_align->q + bias, NWstat.len);
				aligned->start = q_s + Frag_align->start;
			}
			Stat.pos -= (NWstat.len - NWstat.tGaps);
			Stat.score = NWstat.score;
			Stat.len = NWstat.len;
			Stat.match = NWstat.match;
			Stat.tGaps = NWstat.tGaps;
			Stat.qGaps = NWstat.qGaps;
		} else if(aligned) {
			aligned->start = q_s;
		}
	}
	
	return Stat;
}

void skipTrailAln(Aln *aligned, Aln *Frag_align, AlnScore *Stat, const long unsigned *tseq, const unsigned char *qseq, int t_s, int t_len, int q_s, int q_len, const int bandwidth, NWmat *matrices) {
	if(aligned) {
		aligned->end = q_len - q_s;
		Frag_align->end = 0;
	}
}

void trailTailAln(Aln *aligned, Aln *Frag_align, AlnScore *Stat, const long unsigned *tseq, const unsigned char *qseq, int t_s, int t_len, int q_s, int q_len, const int bandwidth, NWmat *matrices) {
	
	int band, bias, q_e, t_e;
	AlnScore NWstat;
	
	/* Get intervals in query and template to align */
	q_e = q_len;
	t_e = t_len;
	if(((q_len - q_s) << 1) < (t_len - t_s) || (q_len - q_s + bandwidth) < (t_len - t_s)) { // big trailing template gap, cut down
		//t_e = t_s + MIN(bandwidth, ((q_len - q_s) << 1));
		t_e = q_len - q_s;
		t_e = t_s + (t_e + (t_e < bandwidth ? t_e : bandwidth));
	} else if(((t_len - t_s) << 1) < (q_len - q_s) || (t_len - t_s + bandwidth) < (q_len - q_s)) { // big leading query gap, cut down
		//q_e = q_s + MIN(bandwidth, ((t_len - t_s) << 1));
		q_e = t_len - t_s;
		q_e = q_s + (q_e + (q_e < bandwidth ? q_e : bandwidth));
	}
	
	/* align trailing gap */
	if(t_e - t_s > 0 && q_e - q_s > 0) {
		band = abs(t_e - t_s - q_e + q_s) + bandwidth;
		if(q_e - q_s <= band || t_e - t_s <= band) {//|| abs(t_e - t_s - q_e - q_s) >= 32) {
			if(Frag_align) {
				NWstat = NW(tseq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, Frag_align, matrices, t_len);
			} else {
				NWstat = NW_score(tseq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, matrices, t_len);
			}
		} else if(Frag_align) {
			NWstat = NW_band(tseq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, Frag_align, band, matrices, t_len);
			//NWstat = NW(tseq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, Frag_align, matrices, t_len);
		} else {
			NWstat = NW_band_score(tseq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, band, matrices, t_len);
		}
		
		if(Frag_align) {
			/* trim trailing gaps */
			if(t_e == t_len) {
				bias = NWstat.len - 1;
				while(bias && (Frag_align->t[bias] == 5 || Frag_align->q[bias] == 5)) {
					if(Frag_align->t[bias] == 5) {
						--NWstat.tGaps;
						++(Frag_align->end);
					} else {
						--NWstat.qGaps;
					}
					--bias;
				}
				++bias;
				
				if(bias != NWstat.len) {
					//NWstat.score -= (W1 + (NWstat.len - bias) * U);
					NWstat.len = bias;
				}
			}
			
			memcpy(aligned->t + Stat->len, Frag_align->t, NWstat.len);
			memcpy(aligned->s + Stat->len, Frag_align->s, NWstat.len);
			memcpy(aligned->q + Stat->len, Frag_align->q, NWstat.len);
		}
		
		Stat->score += NWstat.score;
		Stat->len += NWstat.len;
		Stat->match += NWstat.match;
		Stat->tGaps += NWstat.tGaps;
		Stat->qGaps += NWstat.qGaps;
	} else if(Frag_align) {
		Frag_align->end = 0;
	}
	
	if(aligned) {
		aligned->end = q_len - q_e + Frag_align->end;
	}
}

AlnScore KMA(const HashMapCCI *template_index, const unsigned char *qseq, int q_len, int q_start, int q_end, Aln *aligned, Aln *Frag_align, int min, int max, int mq, double scoreT, AlnPoints *points, NWmat *matrices) {
	
	int i, j, k, bias, prev, start, stop, t_len, value, end, mem_count, band;
	int t_l, t_s, t_e, q_s, q_e, score, shifter, kmersize, U, M, *seeds, **d;
	long unsigned key, mask;
	unsigned char nuc;
	AlnScore Stat, NWstat;
	Penalties *rewards;
	
	/* Extract indexes and template sequence */
	rewards = points->rewards;
	U = rewards->U;
	M = rewards->M;
	d = rewards->d;
	t_len = template_index->len;
	kmersize = template_index->kmerindex;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	key = 0;
	
	/* circular, skip boundaries */
	if(min < max) {
		min = 0;
		max = t_len;
	}
	aligned->start = 0;
	aligned->end = 0;
	
	/* find seeds */
	if(points->len) {
		mem_count = points->len;
	} else {
		mem_count = 0;
		i = q_start;
		while(i < q_end) {
			end = charpos(qseq, 4, i, q_len);
			if(end == -1) {
				end = q_end;
			}
			
			if(i < end - kmersize) {
				key = makeKmer(qseq, i, kmersize - 1);
				i += (kmersize - 1);
			} else {
				i = end + 1;
			}
			
			while(i < end) {
				key = ((key << 2) | qseq[i]) & mask;
				value = hashMapCCI_get_bound(template_index, key, min, max, shifter);
				
				if(value == 0) {
					++i;
				} else if(0 < value) {
					i -= (kmersize - 1);
					
					/* backseed for overlapping seeds */
					prev = value - 2;
					for(j = i - 1; 0 <= j && 0 <= prev && qseq[j] == getNuc(template_index->seq, prev); --j) {
						--prev;
					}
					
					/* get start positions */
					points->qStart[mem_count] = j + 1;
					points->tStart[mem_count] = prev + 2;
					
					/* skip k-mer bases */
					value += (kmersize - 1);
					i += kmersize;
					
					/* extend */
					while(i < end && value < t_len && qseq[i] == getNuc(template_index->seq, value)) {
						++i;
						++value;
					}
					
					/* get end positions */
					points->qEnd[mem_count] = i;
					points->tEnd[mem_count] = value + 1;
					
					/* calculate weight */
					points->weight[mem_count] = (points->qEnd[mem_count] - points->qStart[mem_count]);
					++mem_count;
					
					/* realloc seeding points */
					if(mem_count == points->size) {
						seedPoint_realloc(points, points->size << 1);
					}
					
					/* update position */
					if(i < end - kmersize) {
						key = makeKmer(qseq, i, kmersize - 1);
						i += (kmersize - 1);
					} else {
						i = end + 1;
					}
				} else {
					/* move counter back */
					i -= (kmersize - 1);
					
					/* get position in hashmap */
					seeds = hashMapCCI_getDubPos(template_index, key, value, shifter);
					
					/* get all mems */
					bias = i;
					while(seeds) {
						/* get mem info */
						k = i;
						/* backseed for overlapping seeds */
						value = abs(*seeds);
						prev = value - 2;
						for(j = k - 1; 0 <= j && 0 <= prev && qseq[j] == getNuc(template_index->seq, prev); --j) {
							--prev;
						}
						
						/* get start positions */
						points->qStart[mem_count] = j + 1;
						points->tStart[mem_count] = prev + 2;
						
						/* skip k-mer bases */
						value += (kmersize - 1);
						k += kmersize;
						
						/* extend */
						while(k < end && value < t_len && qseq[k] == getNuc(template_index->seq, value)) {
							++k;
							++value;
						}
						
						/* get end positions */
						points->qEnd[mem_count] = k;
						points->tEnd[mem_count] = value + 1;
						
						/* calculate weight */
						points->weight[mem_count] = (points->qEnd[mem_count] - points->qStart[mem_count]);
						++mem_count;
						
						/* realloc seeding points */
						if(mem_count == points->size) {
							seedPoint_realloc(points, points->size << 1);
						}
						
						if(bias < k) {
							bias = k;
						}
						
						seeds = hashMapCCI_getNextDubPos(template_index, seeds, key, min, max, shifter);
					}
					i = bias + 1;
					
					/* update position */
					if(i < end - kmersize) {
						key = makeKmer(qseq, i, kmersize - 1);
						i += (kmersize - 1);
					} else {
						i = end + 1;
					}
					
				}
			}
			i = end + 1;
		}
	}
	aligned->mapQ = 0;
	
	if(mem_count) {
		points->len = mem_count;
	} else {
		Stat.score = 0;
		Stat.len = 1;
		Stat.match = 0;
		Stat.tGaps = 0;
		Stat.qGaps = 0;
		Stat.pos = 0;
		aligned->s[0] = 0;
		aligned->len = 0;
		points->len = 0;
		return Stat;
	}
	
	/* get best seed chain, returns best starting point */
	start = chainSeedsPtr(points, q_len, t_len, kmersize, &aligned->mapQ);
	score = points->score[start];
	if(aligned->mapQ < mq || score < kmersize) {
		Stat.score = 0;
		Stat.len = 1;
		Stat.match = 0;
		Stat.tGaps = 0;
		Stat.qGaps = 0;
		Stat.pos = 0;
		aligned->s[0] = 0;
		aligned->len = 0;
		points->len = 0;
		return Stat;
	}
	
	/* trim seeds */
	trimSeedsPtr(points, start);
	
	/* align leading tail */
	Stat = leadTailAlnPtr(aligned, Frag_align, template_index->seq, qseq, points->tStart[start] - 1, t_len, points->qStart[start], template_index->bandwidth, matrices);
	
	/* piece seeds together */
	stop = 1;
	while(stop) {
		/* MEM */
		q_s = points->qStart[start];
		end = points->qEnd[start] - q_s;
		memcpy(aligned->t + Stat.len, qseq + q_s, end);
		memset(aligned->s + Stat.len, '|', end);
		memcpy(aligned->q + Stat.len, qseq + q_s, end);
		Stat.len += end;
		Stat.match += end;
		end = points->qEnd[start];
		for(i = points->qStart[start]; i < end; ++i) {
			nuc = qseq[i];
			Stat.score += d[nuc][nuc];
		}
		
		/* join MEMs */
		if(points->next[start]) {
			/* get positions between seed-extends */
			q_s = points->qEnd[start];
			t_s = points->tEnd[start] - 1;
			start = points->next[start];
			
			/* check if next MEM is a semi match */
			if(points->qStart[start] < q_s) {
				points->tStart[start] += (q_s - points->qStart[start]);
				points->qStart[start] = q_s;
			}
			t_e = points->tStart[start] - 1;
			
			if(t_e < t_s) {
				if(t_s <= points->tEnd[start]) {
					points->qStart[start] += (t_s - t_e);
					t_e = t_s;
					t_l = t_e - t_s;
				} else {
					/* circular joining */
					Frag_align->pos = t_len;
					t_l = t_len - t_s + t_e;
				}
			} else {
				t_l = t_e - t_s;
			}
			q_e = points->qStart[start];
			
			/* piece seed-extends together */
			if(abs(t_l - q_e + q_s) * U > q_len * M || t_l > q_len || q_e - q_s > (q_len >> 1)) {
				/* gap is too big to give a positive score */
				Stat.score = 0;
				Stat.len = 1;
				Stat.match = 0;
				Stat.tGaps = 0;
				Stat.qGaps = 0;
				aligned->s[0] = 0;
				aligned->len = 0;
				points->len = 0;
				return Stat;
			}
			if((t_l > 0 || q_e - q_s > 0)) {
				band = abs(t_l - q_e + q_s) + template_index->bandwidth;
				if(q_e - q_s <= band || t_l <= band) {// || abs(t_e - t_s - q_e - q_s) >= 32) {
					NWstat = NW(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, Frag_align, matrices, t_len);
				} else {
					NWstat = NW_band(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, Frag_align, band, matrices, t_len);
					//NWstat = NW(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, Frag_align, matrices, t_len);
				}
				memcpy(aligned->t + Stat.len, Frag_align->t, NWstat.len);
				memcpy(aligned->s + Stat.len, Frag_align->s, NWstat.len);
				memcpy(aligned->q + Stat.len, Frag_align->q, NWstat.len);
				Stat.score += NWstat.score;
				Stat.len += NWstat.len;
				Stat.match += NWstat.match;
				Stat.tGaps += NWstat.tGaps;
				Stat.qGaps += NWstat.qGaps;
			}
		} else {
			stop = 0;
		}
	}
	
	/* align trailing tail */
	trailTailAlnPtr(aligned, Frag_align, &Stat, template_index->seq, qseq, points->tEnd[start] - 1, t_len, points->qEnd[start], q_len, template_index->bandwidth, matrices);
	aligned->s[Stat.len] = 0;
	aligned->len = Stat.len;
	points->len = 0;
	
	return Stat;
}

AlnScore KMA_score(const HashMapCCI *template_index, const unsigned char *qseq, int q_len, int q_start, int q_end, const CompDNA *qseq_comp, int mq, double scoreT, AlnPoints *points, NWmat *matrices) {
	
	int i, j, k, l, bias, prev, start, stop, t_len, value, end, U, M, band;
	int t_l, t_s, t_e, q_s, q_e, mem_count, score, kmersize, *seeds, **d;
	unsigned mapQ, shifter, cPos, iPos;
	long unsigned key;
	unsigned char nuc;
	AlnScore Stat, NWstat;
	Penalties *rewards;
	
	/* Extract indexes and template sequence */
	rewards = points->rewards;
	U = rewards->U;
	M = rewards->M;
	d = rewards->d;
	t_len = template_index->len;
	kmersize = template_index->kmerindex;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	
	/* find seeds */
	if(points->len) {
		mem_count = points->len;
	} else {
		mem_count = 0;
		j = q_start;
		for(i = 1; i <= qseq_comp->N[0]; ++i) {
			if(i != qseq_comp->N[0]) {
				end = qseq_comp->N[i] - kmersize + 1;
			} else {
				end = q_end - kmersize + 1;
			}
			while(j < end) {
				getKmer_macro(key, qseq_comp->seq, j, cPos, iPos, shifter);
				value = hashMapCCI_get(template_index, key, shifter);
				
				if(value == 0) {
					++j;
				} else if(0 < value) {
					/* backseed for ambiguos seeds */
					prev = value - 2;
					for(k = j - 1; 0 <= k && 0 <= prev && qseq[k] == getNuc(template_index->seq, prev); --k) {
						--prev;
					}
					
					/* get start positions */
					points->qStart[mem_count] = k + 1;
					points->tStart[mem_count] = prev + 2;
					
					
					/* skip k-mer bases */
					value += (kmersize - 1);
					j += kmersize;
					
					/* extend */
					end += (kmersize - 1);
					while(j < end && value < t_len && qseq[j] == getNuc(template_index->seq, value)) {
						++j;
						++value;
					}
					end -= (kmersize - 1);
					
					/* get end positions */
					points->qEnd[mem_count] = j;
					points->tEnd[mem_count] = value + 1;
					
					/* calculate weight */
					points->weight[mem_count] = (points->qEnd[mem_count] - points->qStart[mem_count]);
					++mem_count;
					
					/* realloc seeding points */
					if(mem_count == points->size) {
						seedPoint_realloc(points, points->size << 1);
					}
				} else {
					/* get position in hashmap */
					seeds = hashMapCCI_getDubPos(template_index, key, value, shifter);
					
					/* get all mems */
					bias = j;
					while(seeds) {
						/* get mem info */
						l = j;
						/* backseed for overlapping seeds */
						value = abs(*seeds);
						prev = value - 2;
						for(k = l - 1; 0 <= k && 0 <= prev && qseq[k] == getNuc(template_index->seq, prev); --k) {
							--prev;
						}
						
						/* get start positions */
						points->qStart[mem_count] = k + 1;
						points->tStart[mem_count] = prev + 2;
						
						/* skip k-mer bases */
						value += (kmersize - 1);
						l += kmersize;
						
						/* extend */
						end += (kmersize - 1);
						while(l < end && value < t_len && qseq[l] == getNuc(template_index->seq, value)) {
							++l;
							++value;
						}
						end -= (kmersize - 1);
						
						/* get end positions */
						points->qEnd[mem_count] = l;
						points->tEnd[mem_count] = value + 1;
						
						/* calculate weight */
						points->weight[mem_count] = (points->qEnd[mem_count] - points->qStart[mem_count]);
						++mem_count;
						
						/* realloc seeding points */
						if(mem_count == points->size) {
							seedPoint_realloc(points, points->size << 1);
						}
						
						if(bias < l) {
							bias = l;
						}
						
						seeds = hashMapCCI_getNextDubPos(template_index, seeds, key, 0, t_len, shifter);
					}
					j = bias + 1;
				}
			}
			j = qseq_comp->N[i] + 1;
		}
	}
	mapQ = 0;
	
	if(mem_count) {
		points->len = mem_count;
	} else {
		Stat.score = 0;
		Stat.len = 1;
		Stat.match = 0;
		Stat.tGaps = 0;
		Stat.qGaps = 0;
		Stat.pos = 0;
		points->len = 0;
		return Stat;
	}
	
	/* get best seed chain, returns best starting point */
	start = chainSeedsPtr(points, q_len, t_len, kmersize, &mapQ);
	score = points->score[start];
	if(mapQ < mq || score < kmersize) {
		Stat.score = 0;
		Stat.len = 1;
		Stat.match = 0;
		Stat.tGaps = 0;
		Stat.qGaps = 0;
		Stat.pos = 0;
		points->len = 0;
		return Stat;
	}
	
	/* align leading tail */
	Stat = leadTailAlnPtr(0, 0, template_index->seq, qseq, points->tStart[start] - 1, t_len, points->qStart[start], template_index->bandwidth, matrices);
	
	/* piece seeds together */
	stop = 1;
	while(stop) {
		/* MEM */
		q_s = points->qStart[start];
		end = points->qEnd[start] - q_s;
		Stat.len += end;
		Stat.match += end;
		end = points->qEnd[start];
		for(i = points->qStart[start]; i < end; ++i) {
			nuc = qseq[i];
			Stat.score += d[nuc][nuc];
		}
		
		/* join MEMs */
		if(points->next[start]) {
			/* get positions between seed-extends */
			q_s = points->qEnd[start];
			t_s = points->tEnd[start] - 1;
			start = points->next[start];
			
			/* check if next MEM is a semi match, or a circular joining */
			if(points->qStart[start] < q_s) {
				points->tStart[start] += (q_s - points->qStart[start]);
				points->qStart[start] = q_s;
			}
			t_e = points->tStart[start] - 1;
			
			if(t_e < t_s) {
				if(t_s <= points->tEnd[start]) {
					points->qStart[start] += (t_s - t_e);
					t_e = t_s;
					t_l = t_e - t_s;
				} else {
					/* circular joining */
					t_l = t_len - t_s + t_e;
				}
			} else {
				t_l = t_e - t_s;
			}
			q_e = points->qStart[start];
			
			/* piece seed-extends together */
			if(abs(t_l - q_e + q_s) * U > q_len * M || t_l > q_len || q_e - q_s > (q_len >> 1)) {
				/* gap is too big to give a positive score */
				Stat.score = 0;
				Stat.len = 1;
				Stat.match = 0;
				Stat.tGaps = 0;
				Stat.qGaps = 0;
				points->len = 0;
				return Stat;
			}
			if((t_l > 0 || q_e - q_s > 0)) {
				band = abs(t_l - q_e + q_s) + template_index->bandwidth;
				if(q_e - q_s <= band || t_l <= band) {
					NWstat = NW_score(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, matrices, t_len);
				} else {
					NWstat = NW_band_score(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, band, matrices, t_len);
				}
				Stat.score += NWstat.score;
				Stat.len += NWstat.len;
				Stat.match += NWstat.match;
				Stat.tGaps += NWstat.tGaps;
				Stat.qGaps += NWstat.qGaps;
			}
		} else {
			stop = 0;
		}
	}
	
	/* align trailing tail */
	trailTailAlnPtr(0, 0, &Stat, template_index->seq, qseq, points->tEnd[start] - 1, t_len, points->qEnd[start], q_len, template_index->bandwidth, matrices);
	points->len = 0;
	
	return Stat;
}

int preseed(const HashMapCCI *template_index, unsigned char *qseq, int q_len) {
	
	static int exhaustive = 1;
	int i, shifter, kmersize, len;
	
	if(exhaustive) {
		exhaustive = q_len;
		return 0;
	}
	
	kmersize = template_index->kmerindex;
	len = template_index->len;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	for(i = 0; i < q_len; i += kmersize) {
		if(hashMapCCI_get_bound(template_index, makeKmer(qseq, i, kmersize), 0, len, shifter)) {
			return 0;
		}
	}
	
	return i;
}

void intcpy(int *dest, int *src, int size) {
	
	*dest = *src;
	while(--size) {
		*++dest = *++src;
	}
}

int anker_rc(const HashMapCCI *template_index, unsigned char *qseq, int q_len, int q_start, int q_end, AlnPoints *points) {
	
	static int one2one = 0;
	int i, j, k, rc, end, score, score_r, value, t_len, prev, bias;
	int bestScore, mem_count, totMems, shifter, kmersize, *seeds;
	long unsigned key, mask;
	
	if(!template_index) {
		one2one = q_len;
		return 0;
	}
	
	t_len = template_index->len;
	kmersize = template_index->kmerindex;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	key = 0;
	
	/* find seeds */
	bestScore = 0;
	score = 0;
	score_r = 0;
	mem_count = 0;
	totMems = 0;
	points->len = 0;
	for(rc = 0; rc < 2; ++rc) {
		if(rc) {
			strrc(qseq, q_len);
			score = score_r;
			points->len = mem_count;
			i = q_len - q_start;
			q_start = q_len - q_end;
			q_end = i;
			i = q_start;
		}
		i = q_start ? q_start : preseed(template_index, qseq, q_end - q_start);
		score_r = 0;
		mem_count = 0;
		
		while(i < q_end) {
			end = charpos(qseq, 4, i, q_len);
			if(end == -1) {
				end = q_end;
			}
			
			if(i < end - kmersize) {
				key = makeKmer(qseq, i, kmersize - 1);
				i += (kmersize - 1);
			} else {
				i = end + 1;
			}
			
			while(i < end) {
				key = ((key << 2) | qseq[i]) & mask;
				value = hashMapCCI_get(template_index, key, shifter);
				
				if(value == 0) {
					++i;
				} else if(0 < value) {
					i -= (kmersize - 1);
					
					/* backseed for ambiguos seeds */
					prev = value - 2;
					for(j = i - 1; 0 <= j && 0 <= prev && qseq[j] == getNuc(template_index->seq, prev); --j) {
						--prev;
						++score_r;
					}
					
					/* get start positions */
					points->qStart[totMems] = j + 1;
					points->tStart[totMems] = prev + 2;
					
					/* skip k-mer bases */
					value += (kmersize - 1);
					i += kmersize;
					score_r += kmersize;
					
					/* extend */
					while(i < end && value < t_len && qseq[i] == getNuc(template_index->seq, value)) {
						++i;
						++value;
						++score_r;
					}
					
					/* get end positions */
					points->qEnd[totMems] = i;
					points->tEnd[totMems] = value + 1;
					
					/* calculate weight */
					points->weight[totMems] = (points->tEnd[totMems] - points->tStart[totMems]);
					++mem_count;
					++totMems;
					
					/* realloc seeding points */
					if(totMems == points->size) {
						seedPoint_realloc(points, points->size << 1);
					}
					
					/* update position */
					if(i < end - kmersize) {
						key = makeKmer(qseq, i, kmersize - 1);
						i += (kmersize - 1);
					} else {
						i = end + 1;
					}
				} else {
					/* move counter back */
					i -= (kmersize - 1);
					score_r += kmersize;
					
					/* get position in hashmap */
					seeds = hashMapCCI_getDubPos(template_index, key, value, shifter);
					
					/* get all mems */
					bias = i;
					while(seeds) {
						/* get mem info */
						k = i;
						/* backseed for overlapping seeds */
						value = abs(*seeds);
						prev = value - 2;
						for(j = k - 1; 0 <= j && 0 <= prev && qseq[j] == getNuc(template_index->seq, prev); --j) {
							--prev;
						}
						
						/* get start positions */
						points->qStart[totMems] = j + 1;
						points->tStart[totMems] = prev + 2;
						
						/* skip k-mer bases */
						value += (kmersize - 1);
						k += kmersize;
						
						/* extend */
						while(k < end && value < t_len && qseq[k] == getNuc(template_index->seq, value)) {
							++k;
							++value;
						}
						
						/* get end positions */
						points->qEnd[totMems] = k;
						points->tEnd[totMems] = value + 1;
						
						/* calculate weight */
						points->weight[totMems] = (points->qEnd[totMems] - points->qStart[totMems]);
						++mem_count;
						++totMems;
						
						/* realloc seeding points */
						if(totMems == points->size) {
							seedPoint_realloc(points, points->size << 1);
						}
						
						if(bias < k) {
							bias = k;
						}
						
						seeds = hashMapCCI_getNextDubPos(template_index, seeds, key, 0, t_len, shifter);
					}
					
					/* add best anker score */
					score_r += (bias - i);
					i = bias + 1;
					
					/* update position */
					if(i < end - kmersize) {
						key = makeKmer(qseq, i, kmersize - 1);
						i += (kmersize - 1);
					} else {
						i = end + 1;
					}
				}
			}
			i = end + 1;
		}
		
		if(bestScore < score_r) {
			bestScore = score_r;
		}
	}
	
	if(one2one && bestScore < kmersize && bestScore * kmersize < (q_len - kmersize - bestScore)) {
		bestScore = 0;
		points->len = 0;
	} else if(bestScore == score) {
		/* forward match */
		strrc(qseq, q_len);
	} else {
		/* apply sam-flag */
		/* reverse match */
		/*
		if(flag | 64) { // Forward read
			flag ^= 32;
		}
		flag ^= 16;
		*/
		/* move mems down */
		if(points->len) {
			intcpy(points->tStart, points->tStart + points->len, mem_count);
			intcpy(points->tEnd, points->tEnd + points->len, mem_count);
			intcpy(points->qStart, points->qStart + points->len, mem_count);
			intcpy(points->qEnd, points->qEnd + points->len, mem_count);
			intcpy(points->weight, points->weight + points->len, mem_count);
		}
		points->len = mem_count;
	}
	
	return bestScore;
}

int anker_rc_comp(const HashMapCCI *template_index, unsigned char *qseq, unsigned char *qseq_r, CompDNA *qseq_comp, CompDNA *qseq_r_comp, int q_start, int q_end, AlnPoints *points) {
	
	static int one2one = 0;
	int i, j, k, rc, end, score, score_r, value, t_len, q_len, prev;
	int bestScore, mem_count, totMems, shifter, kmersize, bias, *Ns, *seeds;
	unsigned cPos, iPos;
	long unsigned key, mask, *seq;
	
	if(!template_index) {
		one2one = *((int *)(qseq_r));
		return 0;
	}
	
	q_len = qseq_comp->seqlen;
	t_len = template_index->len;
	kmersize = template_index->kmerindex;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	key = 0;
	
	/* find seeds */
	bestScore = 0;
	score = 0;
	score_r = 0;
	mem_count = 0;
	totMems = 0;
	points->len = 0;
	seq = qseq_comp->seq;
	Ns = qseq_comp->N;
	Ns[*Ns] = q_len;
	for(rc = 0; rc < 2; ++rc) {
		if(rc) {
			qseq = qseq_r;
			seq = qseq_r_comp->seq;
			score = score_r;
			points->len = mem_count;
			Ns = qseq_r_comp->N;
			Ns[*Ns] = q_len;
			i = q_len - q_start;
			q_start = q_len - q_end;
			q_end = i;
			i = q_start;
		} else if(q_start) {
			i = q_start;
		} else {
			i = preseed(template_index, qseq, q_end - q_start);
		}
		score_r = 0;
		mem_count = 0;
		
		while(i < q_end) {
			end = *++Ns - kmersize + 1;
			while(i < end) {
				getKmer_macro(key, seq, i, cPos, iPos, shifter);
				value = hashMapCCI_get(template_index, key, shifter);
				
				if(value == 0) {
					++i;
				} else if(0 < value) {
					/* backseed for ambiguos seeds */
					prev = value - 2;
					for(j = i - 1; 0 <= j && 0 <= prev && qseq[j] == getNuc(template_index->seq, prev); --j) {
						--prev;
						++score_r;
					}
					
					/* get start positions */
					points->qStart[totMems] = j + 1;
					points->tStart[totMems] = prev + 2;
					
					/* skip k-mer bases */
					value += (kmersize - 1);
					i += kmersize;
					score_r += kmersize;
					
					/* extend */
					while(i < end && value < t_len && qseq[i] == getNuc(template_index->seq, value)) {
						++i;
						++value;
						++score_r;
					}
					
					/* get end positions */
					points->qEnd[totMems] = i;
					points->tEnd[totMems] = value + 1;
					
					/* calculate weight */
					points->weight[totMems] = (points->tEnd[totMems] - points->tStart[totMems]);
					++mem_count;
					++totMems;
					
					/* realloc seeding points */
					if(totMems == points->size) {
						seedPoint_realloc(points, points->size << 1);
					}
					
					/* update position */
					++i;
				} else {
					score_r += kmersize;
					
					/* get position in hashmap */
					seeds = hashMapCCI_getDubPos(template_index, key, value, shifter);
					
					/* get all mems */
					bias = i;
					while(seeds) {
						/* get mem info */
						k = i;
						/* backseed for overlapping seeds */
						value = abs(*seeds);
						prev = value - 2;
						for(j = k - 1; 0 <= j && 0 <= prev && qseq[j] == getNuc(template_index->seq, prev); --j) {
							--prev;
						}
						
						/* get start positions */
						points->qStart[totMems] = j + 1;
						points->tStart[totMems] = prev + 2;
						
						/* skip k-mer bases */
						value += (kmersize - 1);
						k += kmersize;
						
						/* extend */
						while(k < end && value < t_len && qseq[k] == getNuc(template_index->seq, value)) {
							++k;
							++value;
						}
						
						/* get end positions */
						points->qEnd[totMems] = k;
						points->tEnd[totMems] = value + 1;
						
						/* calculate weight */
						points->weight[totMems] = (points->qEnd[totMems] - points->qStart[totMems]);
						++mem_count;
						++totMems;
						
						/* realloc seeding points */
						if(totMems == points->size) {
							seedPoint_realloc(points, points->size << 1);
						}
						
						if(bias < k) {
							bias = k;
						}
						
						seeds = hashMapCCI_getNextDubPos(template_index, seeds, key, 0, t_len, shifter);
					}
					/* add best anker score */
					score_r += (bias - i);
					i = bias + 1;
				}
			}
			i = end + kmersize;
		}
		
		if(bestScore < score_r) {
			bestScore = score_r;
		}
	}
	
	if(one2one && bestScore < kmersize && bestScore * kmersize < (q_len - kmersize - bestScore)) {
		bestScore = 0;
		points->len = 0;
	} else if(bestScore == score) {
		return bestScore;
	} else {
		/* move mems down */
		if(points->len) {
			intcpy(points->tStart, points->tStart + points->len, mem_count);
			intcpy(points->tEnd, points->tEnd + points->len, mem_count);
			intcpy(points->qStart, points->qStart + points->len, mem_count);
			intcpy(points->qEnd, points->qEnd + points->len, mem_count);
			intcpy(points->weight, points->weight + points->len, mem_count);
		}
		points->len = mem_count;
		return -bestScore;
	}
	
	return bestScore;
}
