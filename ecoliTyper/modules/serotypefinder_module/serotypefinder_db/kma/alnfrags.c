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
#include "align.h"
#include "alnfrags.h"
#include "ankers.h"
#include "chain.h"
#include "compdna.h"
#include "frags.h"
#include "filebuff.h"
#include "hashmapcci.h"
#include "pherror.h"
#include "qseqs.h"
#include "sam.h"
#include "stdnuc.h"
#include "stdstat.h"
#include "threader.h"
#include "updatescores.h"
#define exchange(src1, src2, tmp) tmp = src1; src1 = src2; src2 = tmp;
#define mrcheck(mrc, Stat, q_len, t_len) ((mrc * q_len <= Stat.len - Stat.qGaps) || (mrc * t_len <= Stat.len - Stat.tGaps))

int (*alnFragsPE)(HashMapCCI**, int*, int*, int, double, double, double, int, CompDNA*, CompDNA*, CompDNA*, CompDNA*, unsigned char*, unsigned char*, unsigned char*, unsigned char*, Qseqs*, Qseqs*, int, int*, int*, long unsigned*, long unsigned*, int*, int*, int*, int*, int*, int*, int, long*, FILE*, AlnPoints*, NWmat*, volatile int*, volatile int*) = alnFragsUnionPE;

int alnFragsSE_old(HashMapCCI **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, double mrc, int minlen, int rc_flag, CompDNA *qseq_comp, CompDNA *qseq_r_comp, unsigned char *qseq, unsigned char *qseq_r, int q_len, int kmersize, Qseqs *header, int *bestTemplates, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int *flag, int *best_read_score, int seq_in, long *seq_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, bestHits, aln_len, start, end, Wl, arc, rc;
	int q_start, q_end, t_len, *qBoundPtr;
	double score, bestScore;
	AlnScore alnStat;
	
	/* reverse complement seq */
	if(rc_flag < 0) {
		if(qseq_r_comp->size < qseq_comp->size) {
			freeComp(qseq_r_comp);
			allocComp(qseq_r_comp, qseq_comp->size);
		}
		rc_comp(qseq_comp, qseq_r_comp);
		unCompDNA(qseq_r_comp, qseq_r);
		qseq_r_comp->N[0]++;
		qseq_r_comp->N[qseq_r_comp->N[0]] = q_len;
	}
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = q_len;
	
	bestScore = 0;
	*best_read_score = 0;
	bestHits = 0;
	Wl = -NWmatrices->rewards->Wl;
	arc = points->len;
	points->len = 0;
	
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		/* check if index DB is loaded */
		lock(excludeDB);
		if(template >= 0 && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(templates_index[template], seq_in, template_lengths[template], kmersize, seq_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(templates_index[-template], seq_in, template_lengths[-template], kmersize, seq_indexes[-template]);
		}
		unlock(excludeDB);
		
		/* q-bound */
		if(2 * sizeof(int) + 1 < header->len && header->seq[header->len - 2 * sizeof(int) - 1] == 0) {
			qBoundPtr = (int*) (header->seq + (header->len - 2 * sizeof(int)));
			q_start = *qBoundPtr;
			q_end = *++qBoundPtr;
		} else {
			q_start = 0;
			q_end = q_len;
		}
		
		/* check ankers */
		if(arc) {
			rc = anker_rc_comp(templates_index[abs(template)], qseq, qseq_r, qseq_comp, qseq_r_comp, q_start, q_end, points);
			if(rc < 0) {
				if(0 < template) {
					template = -template;
				}
				alnStat = KMA_score(templates_index[-template], qseq_r, q_len, q_len - q_end, q_len - q_start, qseq_r_comp, mq, scoreT, points, NWmatrices);
			} else if(rc) {
				if(template < 0) {
					template = -template;
				}
				alnStat = KMA_score(templates_index[template], qseq, q_len, q_start, q_end, qseq_comp, mq, scoreT, points, NWmatrices);
			} else {
				alnStat.score = 0;
				alnStat.pos = 0;
				alnStat.len = 0;
				alnStat.tGaps = 0;
				alnStat.qGaps = 0;
				points->len = 0;
			}
		} else {
			if(template < 0) {
				alnStat = KMA_score(templates_index[-template], qseq_r, q_len, q_len - q_end, q_len - q_start, qseq_r_comp, mq, scoreT, points, NWmatrices);
			} else {
				alnStat = KMA_score(templates_index[template], qseq, q_len, q_start, q_end, qseq_comp, mq, scoreT, points, NWmatrices);
			}
		}
		
		/* get read score */
		aln_len = alnStat.len;
		start = alnStat.pos;
		end = start + aln_len - alnStat.tGaps;
		t_len = template_lengths[abs(template)];
		if(t_len < end) {
			end -= t_len;
		}
		
		/* penalty for non complete mapping */
		read_score = alnStat.score;
		/* full gene award */
		if(start == 0) {
			read_score += Wl;
		}
		if(end == t_len) {
			read_score += Wl;
		}
		//read_score += (((start != 0) + (end != template_lengths[abs(template)])) * Wl);
		
		/* Get normed score check read coverage */
		read_score = alnStat.score;
		if(minlen <= aln_len && mrcheck(mrc, alnStat, q_len, t_len)) {
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
			score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT) {
			if(score > bestScore) { // save as best match
				bestScore = score;
				*best_read_score = read_score;
				*bestTemplates = template;
				*best_start_pos = start;
				*best_end_pos = end;
				bestHits = 1;
			} else if(score == bestScore && read_score > *best_read_score) { // save as best match
				bestScore = score;
				*best_read_score = read_score;
				*bestTemplates = template;
				*best_start_pos = start;
				*best_end_pos = end;
				bestHits = 1;
			} else if(score == bestScore && read_score == *best_read_score) { // update best match
				bestTemplates[bestHits] = template;
				best_start_pos[bestHits] = start;
				best_end_pos[bestHits] = end;
				++bestHits;
			}
		}
	}
	
	if(*best_read_score > kmersize) {
		lock(excludeOut);
		update_Scores_old(qseq, q_len, bestHits, *best_read_score, best_start_pos, best_end_pos, bestTemplates, header, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
		unlock(excludeOut);
		return 0;
	} else {
		*flag |= 4;
	}
	
	return 1;
}

int alnFragsUnionPE_old(HashMapCCI **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, double mrc, double minFrac, int minlen, CompDNA *qseq_comp, CompDNA *qseq_r_comp, CompDNA *qseq_fr_comp, CompDNA *qseq_rr_comp, unsigned char *qseq, unsigned char *qseq_r, unsigned char *qseq_fr, unsigned char *qseq_rr, Qseqs *header, Qseqs *header_r, int kmersize, int *bestTemplates, int *bestTemplates_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int *flag, int *flag_r, int *best_read_score, int *best_read_score_r, int seq_in, long *seq_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, Wl;
	int compScore, bestHits, bestHits_r, aln_len, start, end, arc, rc, t_len;
	double score;
	AlnScore alnStat;
	
	/* unpack qseqs */
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = qseq_comp->seqlen;
	unCompDNA(qseq_r_comp, qseq_r);
	qseq_r_comp->N[0]++;
	qseq_r_comp->N[qseq_r_comp->N[0]] = qseq_r_comp->seqlen;
	if((arc = points->len)) {
		qseq_fr_comp->N[0]++;
		qseq_fr_comp->N[qseq_fr_comp->N[0]] = qseq_fr_comp->seqlen;
		qseq_rr_comp->N[0]++;
		qseq_rr_comp->N[qseq_rr_comp->N[0]] = qseq_rr_comp->seqlen;
		points->len = 0;
	}
	Wl = -NWmatrices->rewards->Wl;
	start = 0;
	end = 0;
	score = 0;
	*best_read_score = 0;
	*best_read_score_r = 0;
	compScore = 0;
	bestHits = 0;
	rc = 1;
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		/* check if index DB is loaded, and rc */
		if(template < 0) {
			if(rc) {
				qseq_comp->N[0]--;
				comp_rc(qseq_comp);
				qseq_comp->N[0]++;
				qseq_r_comp->N[0]--;
				comp_rc(qseq_r_comp);
				qseq_r_comp->N[0]++;
				
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
				
				rc = 0;
			}
		}
		
		lock(excludeDB);
		if(0 <= template && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(templates_index[template], seq_in, template_lengths[template], kmersize, seq_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(templates_index[-template], seq_in, template_lengths[-template], kmersize, seq_indexes[-template]);
		}
		unlock(excludeDB);
		
		/* align qseqs */
		template = abs(template);
		if(arc) {
			rc = anker_rc_comp(templates_index[template], qseq, qseq_fr, qseq_comp, qseq_fr_comp, 0, qseq_comp->seqlen,  points);
			if(rc < 0) {
				/* rc */
				alnStat = KMA_score(templates_index[template], qseq_fr, qseq_fr_comp->seqlen, 0, qseq_fr_comp->seqlen, qseq_fr_comp, mq, scoreT, points, NWmatrices);
			} else if(rc) {
				/* forward */
				matched_templates[t_i] = -matched_templates[t_i];
				alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, 0, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
			} else {
				alnStat.score = 0;
				alnStat.pos = 0;
				alnStat.len = 0;
				alnStat.tGaps = 0;
				alnStat.qGaps = 0;
				points->len = 0;
			}
		} else {
			alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, 0, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
		}
		
		/* get read score */
		aln_len = alnStat.len;
		read_score = alnStat.score;
		t_len = template_lengths[abs(template)];
		if(minlen <= aln_len && 0 < read_score && mrcheck(mrc, alnStat, qseq_comp->seqlen, t_len)) {
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.tGaps;
			if(start == 0) {
				read_score += Wl;
			}
			if(end == t_len) {
				read_score += Wl;
			}
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* adjust to save all, and unite later */
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT && *best_read_score <= read_score) {
			*best_read_score = read_score;
			bestTemplates[t_i] = read_score;
			best_start_pos[t_i] = start;
			best_end_pos[t_i] = end;
		} else {
			bestTemplates[t_i] = 0;
			best_start_pos[t_i] = -1;
			best_end_pos[t_i] = -1;
		}
		
		/* align qseqs */
		template = abs(template);
		if(arc) {
			if(rc < 0) {
				/* rc */
				alnStat = KMA_score(templates_index[template], qseq_rr, qseq_rr_comp->seqlen, 0, qseq_rr_comp->seqlen, qseq_rr_comp, mq, scoreT, points, NWmatrices);
			} else if(rc) {
				/* forward */
				alnStat = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, 0, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
			} else {
				alnStat.score = 0;
				alnStat.pos = 0;
				alnStat.len = 0;
				alnStat.tGaps = 0;
				alnStat.qGaps = 0;
			}
			rc = 1;
		} else {
			alnStat = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, 0, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
		}
		
		
		/* get read score */
		aln_len = alnStat.len;
		read_score = alnStat.score;
		if(minlen <= aln_len && 0 < read_score && mrcheck(mrc, alnStat, qseq_r_comp->seqlen, t_len)) {
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.tGaps;
			if(start == 0) {
				read_score += Wl;
			}
			if(end == t_len) {
				read_score += Wl;
			}
			
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* adjust to sae all, and unite later */
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT && *best_read_score_r <= read_score) {
			*best_read_score_r = read_score;
			bestTemplates_r[t_i] = read_score;
			if(bestTemplates[t_i]) {
				/* Handle negative insertsizes caused by trimming,
				user stupidity or sample error. */
				if(start < best_start_pos[t_i]) {
					best_start_pos[t_i] = start;
				} else {
					best_end_pos[t_i] = end;
				}
			} else {
				best_start_pos[t_i] = start;
				best_end_pos[t_i] = end;
			}
		} else {
			bestTemplates_r[t_i] = 0;
		}
		
		read_score += bestTemplates[t_i];
		if(best_start_pos[t_i] == 0) {
			read_score += Wl;
		}
		if(best_end_pos[t_i] == t_len) {
			read_score += Wl;
		}
		
		
		if(compScore < read_score) {
			compScore = read_score;
		}
	}
	
	/* get rc flag */
	if(arc) {
		rc = 0;
		for(t_i = 1; t_i <= *matched_templates && !rc; ++t_i) {
			rc = matched_templates[t_i] < 0;
		}
	}
	
	if(*best_read_score && *best_read_score_r) {
		/* both matched */
		if(compScore && compScore == (*best_read_score + *best_read_score_r)) {
			/* proper pair */
			bestHits = 0;
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(compScore == (bestTemplates[t_i] + bestTemplates_r[t_i])) {
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				}
			}
			/* check direction of qseqs */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
				lock(excludeOut);
				update_Scores_pe_old(qseq_r, qseq_r_comp->seqlen, qseq, qseq_comp->seqlen, bestHits, compScore, best_start_pos, best_end_pos, bestTemplates, header_r, header, *flag_r, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
				unlock(excludeOut);
			} else {
				if(!rc) {
					strrc(qseq, qseq_comp->seqlen);
					strrc(qseq_r, qseq_r_comp->seqlen);
					*flag ^= 48;
					*flag_r ^= 48;
				}
				lock(excludeOut);
				update_Scores_pe_old(qseq, qseq_comp->seqlen, qseq_r, qseq_r_comp->seqlen, bestHits, compScore, best_start_pos, best_end_pos, bestTemplates, header, header_r, *flag, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
				unlock(excludeOut);
			}
		} else {
			/* unmaided pair */
			bestHits = 0;
			bestHits_r = 0;
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(bestTemplates[t_i] == *best_read_score) {
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				} else if(bestTemplates_r[t_i] == *best_read_score_r) {
					bestTemplates_r[bestHits_r] = matched_templates[t_i];
					best_start_pos[bestHits_r] = best_start_pos[t_i];
					best_end_pos[bestHits_r] = best_end_pos[t_i];
					++bestHits_r;
				}
			}
			
			/* check direction of qseqs */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
			} else if(!rc) {
				strrc(qseq, qseq_comp->seqlen);
				*flag ^= 16;
				*flag_r ^= 32;
			}
			if(*bestTemplates_r < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates_r[t_i] = -bestTemplates_r[t_i];
				}
			} else if(!rc) {
				strrc(qseq_r, qseq_r_comp->seqlen);
				*flag ^= 32;
				*flag_r ^= 16;
			}
			if(*flag & 2) {
				*flag ^= 2;
				*flag_r ^= 2;
			}
			lock(excludeOut);
			update_Scores_old(qseq, qseq_comp->seqlen, bestHits, *best_read_score, best_start_pos, best_end_pos, bestTemplates, header, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
			update_Scores_old(qseq_r, qseq_r_comp->seqlen, bestHits_r, *best_read_score_r, best_start_pos, best_end_pos, bestTemplates_r, header_r, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
		}
		return 0;
	} else if(*best_read_score) {
		bestHits = 0;
		for(t_i = 1; t_i <= *matched_templates; ++t_i) {
			if(bestTemplates[t_i] == *best_read_score) {
				bestTemplates[bestHits] = matched_templates[t_i];
				best_start_pos[bestHits] = best_start_pos[t_i];
				best_end_pos[bestHits] = best_end_pos[t_i];
				++bestHits;
			}
		}
		/* check direction of qseqs */
		if(*bestTemplates < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates[t_i] = -bestTemplates[t_i];
			}
		} else if(!rc) {
			strrc(qseq, qseq_comp->seqlen);
			*flag ^= 16;
			*flag_r ^= 32;
		}
		*flag |= 8;
		*flag_r ^= 4;
		if(*flag & 2) {
			*flag ^= 2;
			*flag_r ^= 2;
		}
		
		lock(excludeOut);
		update_Scores_old(qseq, qseq_comp->seqlen, bestHits, *best_read_score, best_start_pos, best_end_pos, bestTemplates, header, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
		unlock(excludeOut);
		
		return 2;
	} else if(*best_read_score_r) {
		bestHits_r = 0;
		for(t_i = 1; t_i <= *matched_templates; ++t_i) {
			if(bestTemplates_r[t_i] == *best_read_score_r) {
				bestTemplates_r[bestHits_r] = matched_templates[t_i];
				best_start_pos[bestHits_r] = best_start_pos[t_i];
				best_end_pos[bestHits_r] = best_end_pos[t_i];
				++bestHits_r;
			}
		}
		/* check direction of qseqs */
		if(*bestTemplates_r < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates_r[t_i] = -bestTemplates_r[t_i];
			}
		} else if(!rc) {
			strrc(qseq_r, qseq_r_comp->seqlen);
			*flag ^= 32;
			*flag_r ^= 16;
		}
		*flag_r |= 8;
		*flag ^= 4;
		if(*flag_r & 2) {
			*flag ^= 2;
			*flag_r ^= 2;
		}
		lock(excludeOut);
		update_Scores_old(qseq_r, qseq_r_comp->seqlen, bestHits_r, *best_read_score_r, best_start_pos, best_end_pos, bestTemplates_r, header_r, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
		unlock(excludeOut);
		return 1;
	}
	
	return 3;
}

int alnFragsPenaltyPE_old(HashMapCCI **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, double mrc, int minlen, CompDNA *qseq_comp, CompDNA *qseq_r_comp, CompDNA *qseq_fr_comp, CompDNA *qseq_rr_comp, unsigned char *qseq, unsigned char *qseq_r, unsigned char *qseq_fr, unsigned char *qseq_rr, Qseqs *header, Qseqs *header_r, int kmersize, int *bestTemplates, int *bestTemplates_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int *flag, int *flag_r, int *best_read_score, int *best_read_score_r, int seq_in, long *seq_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, Wl, PE;
	int compScore, bestHits, bestHits_r, aln_len, start, end, t_len, arc, rc;
	double score;
	AlnScore alnStat;
	
	/* unpack qseqs */
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = qseq_comp->seqlen;
	unCompDNA(qseq_r_comp, qseq_r);
	qseq_r_comp->N[0]++;
	qseq_r_comp->N[qseq_r_comp->N[0]] = qseq_r_comp->seqlen;
	if((arc = points->len)) {
		qseq_fr_comp->N[0]++;
		qseq_fr_comp->N[qseq_fr_comp->N[0]] = qseq_fr_comp->seqlen;
		qseq_rr_comp->N[0]++;
		qseq_rr_comp->N[qseq_rr_comp->N[0]] = qseq_rr_comp->seqlen;
		points->len = 0;
	}
	Wl = -NWmatrices->rewards->Wl;
	PE = NWmatrices->rewards->PE;
	start = 0;
	end = 0;
	score = 0;
	*best_read_score = 0;
	*best_read_score_r = 0;
	compScore = 0;
	bestHits = 0;
	rc = 1;
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		/* check if index DB is loaded, and rc */
		if(template < 0) {
			if(rc) {
				qseq_comp->N[0]--;
				comp_rc(qseq_comp);
				qseq_comp->N[0]++;
				qseq_r_comp->N[0]--;
				comp_rc(qseq_r_comp);
				qseq_r_comp->N[0]++;
				
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
				
				rc = 0;
			}
		}
		
		lock(excludeDB);
		if(0 <= template && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(templates_index[template], seq_in, template_lengths[template], kmersize, seq_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(templates_index[-template], seq_in, template_lengths[-template], kmersize, seq_indexes[-template]);
		}
		unlock(excludeDB);
		
		/* align qseqs */
		template = abs(template);
		if(arc) {
			rc = anker_rc_comp(templates_index[template], qseq, qseq_fr, qseq_comp, qseq_fr_comp, 0, qseq_comp->seqlen, points);
			if(rc < 0) {
				/* rc */
				alnStat = KMA_score(templates_index[template], qseq_fr, qseq_fr_comp->seqlen, 0, qseq_fr_comp->seqlen, qseq_fr_comp, mq, scoreT, points, NWmatrices);
			} else if(rc) {
				/* forward */
				matched_templates[t_i] = -matched_templates[t_i];;
				alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, 0, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
			} else {
				alnStat.score = 0;
				alnStat.pos = 0;
				alnStat.len = 0;
				alnStat.tGaps = 0;
				alnStat.qGaps = 0;
				points->len = 0;
			}
		} else {
			alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, 0, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
		}
		
		/* get read score */
		aln_len = alnStat.len;
		read_score = alnStat.score;
		t_len = template_lengths[abs(template)];
		if(minlen <= aln_len && 0 < read_score && mrcheck(mrc, alnStat, qseq_comp->seqlen, t_len)) {
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.tGaps;
			if(start == 0) {
				read_score += Wl;
			}
			if(end == t_len) {
				read_score += Wl;
			}
			
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* adjust to see all, and unite later */
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT && *best_read_score <= read_score) {
			*best_read_score = read_score;
			bestTemplates[t_i] = read_score;
			best_start_pos[t_i] = start;
			best_end_pos[t_i] = end;
		} else {
			bestTemplates[t_i] = 0;
			best_start_pos[t_i] = -1;
			best_end_pos[t_i] = -1;
		}
		
		if(arc) {
			if(rc < 0) {
				/* rc */
				alnStat = KMA_score(templates_index[template], qseq_rr, qseq_rr_comp->seqlen, 0, qseq_rr_comp->seqlen, qseq_rr_comp, mq, scoreT, points, NWmatrices);
			} else if(rc) {
				/* forward */
				alnStat = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, 0, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
			} else {
				alnStat.score = 0;
				alnStat.pos = 0;
				alnStat.len = 0;
				alnStat.tGaps = 0;
				alnStat.qGaps = 0;
			}
			rc = 1;
		} else {
			alnStat = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, 0, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
		}
		
		/* get read score */
		aln_len = alnStat.len;
		read_score = alnStat.score;
		if(minlen <= aln_len && 0 < read_score) {
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.tGaps;
			if(start == 0) {
				read_score += Wl;
			}
			if(end == t_len) {
				read_score += Wl;
			}
			
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* adjust to save all, and unite later */
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT && *best_read_score_r <= read_score) {
			*best_read_score_r = read_score;
			bestTemplates_r[t_i] = read_score;
			if(bestTemplates[t_i]) {
				/* Handle negative insertsizes caused by trimming,
				user stupidity or sample error. */
				if(start < best_start_pos[t_i]) {
					best_start_pos[t_i] = start;
				} else {
					best_end_pos[t_i] = end;
				}
			} else {
				best_start_pos[t_i] = start;
				best_end_pos[t_i] = end;
			}
		} else {
			bestTemplates_r[t_i] = 0;
			if(bestTemplates[t_i] != 0) {
				best_start_pos[t_i] = -1;
				best_end_pos[t_i] = -1;
			}
		}
		
		read_score += (bestTemplates[t_i] + PE);
		if(compScore < read_score) {
			compScore = read_score;
		}
	}
	
	/* get rc flag */
	if(arc) {
		rc = 0;
		for(t_i = 1; t_i <= *matched_templates && !rc; ++t_i) {
			rc = matched_templates[t_i] < 0;
		}
	}
	
	if(*best_read_score && *best_read_score_r) {
		/* both matched */
		if(compScore && (*best_read_score + *best_read_score_r) <= compScore) {
			/* proper pair */
			compScore -= PE;
			bestHits = 0;
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(compScore == (bestTemplates[t_i] + bestTemplates_r[t_i])) {
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				}
			}
			/* check direction of qseqs */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
				lock(excludeOut);
				update_Scores_pe_old(qseq_r, qseq_r_comp->seqlen, qseq, qseq_comp->seqlen, bestHits, compScore, best_start_pos, best_end_pos, bestTemplates, header_r, header, *flag_r, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
				unlock(excludeOut);
			} else {
				if(!rc) {
					strrc(qseq, qseq_comp->seqlen);
					strrc(qseq_r, qseq_r_comp->seqlen);
					*flag ^= 48;
					*flag_r ^= 48;
				}
				lock(excludeOut);
				update_Scores_pe_old(qseq, qseq_comp->seqlen, qseq_r, qseq_r_comp->seqlen, bestHits, compScore, best_start_pos, best_end_pos, bestTemplates, header, header_r, *flag, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
				unlock(excludeOut);
			}
		} else {
			/* unmaided pair */
			bestHits = 0;
			bestHits_r = 0;
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(bestTemplates[t_i] == *best_read_score) {
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				} else if(bestTemplates_r[t_i] == *best_read_score_r) {
					bestTemplates_r[bestHits_r] = matched_templates[t_i];
					best_start_pos[bestHits_r] = best_start_pos[t_i];
					best_end_pos[bestHits_r] = best_end_pos[t_i];
					++bestHits_r;
				}
			}
			/* check direction of qseqs */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
			} else if(!rc) {
				strrc(qseq, qseq_comp->seqlen);
				*flag ^= 16;
				*flag_r ^= 32;
			}
			if(*bestTemplates_r < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates_r[t_i] = -bestTemplates_r[t_i];
				}
			} else if(!rc) {
				strrc(qseq_r, qseq_r_comp->seqlen);
				*flag ^= 32;
				*flag_r ^= 16;
			}
			if(*flag & 2) {
				*flag ^= 2;
				*flag_r ^= 2;
			}
			lock(excludeOut);
			update_Scores_old(qseq, qseq_comp->seqlen, bestHits, *best_read_score, best_start_pos, best_end_pos, bestTemplates, header, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
			update_Scores_old(qseq_r, qseq_r_comp->seqlen, bestHits_r, *best_read_score_r, best_start_pos, best_end_pos, bestTemplates_r, header_r, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
		}
		return 0;
	} else if(*best_read_score) {
		bestHits = 0;
		for(t_i = 1; t_i <= *matched_templates; ++t_i) {
			if(bestTemplates[t_i] == *best_read_score) {
				bestTemplates[bestHits] = matched_templates[t_i];
				best_start_pos[bestHits] = best_start_pos[t_i];
				best_end_pos[bestHits] = best_end_pos[t_i];
				++bestHits;
			}
		}
		/* check direction of qseqs */
		if(*bestTemplates < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates[t_i] = -bestTemplates[t_i];
			}
		} else if(!rc) {
			strrc(qseq, qseq_comp->seqlen);
			*flag ^= 16;
			*flag_r ^= 32;
		}
		*flag |= 8;
		*flag_r ^= 4;
		if(*flag & 2) {
			*flag ^= 2;
			*flag_r ^= 2;
		}
		
		lock(excludeOut);
		update_Scores_old(qseq, qseq_comp->seqlen, bestHits, *best_read_score, best_start_pos, best_end_pos, bestTemplates, header, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
		unlock(excludeOut);
		
		return 2;
	} else if(*best_read_score_r) {
		bestHits_r = 0;
		for(t_i = 1; t_i <= *matched_templates; ++t_i) {
			if(bestTemplates_r[t_i] == *best_read_score_r) {
				bestTemplates_r[bestHits_r] = matched_templates[t_i];
				best_start_pos[bestHits_r] = best_start_pos[t_i];
				best_end_pos[bestHits_r] = best_end_pos[t_i];
				++bestHits_r;
			}
		}
		/* check direction of qseqs */
		if(*bestTemplates_r < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates_r[t_i] = -bestTemplates_r[t_i];
			}
		} else if(!rc) {
			strrc(qseq_r, qseq_r_comp->seqlen);
			*flag ^= 32;
			*flag_r ^= 16;
		}
		*flag_r |= 8;
		*flag ^= 4;
		if(*flag_r & 2) {
			*flag ^= 2;
			*flag_r ^= 2;
		}
		lock(excludeOut);
		update_Scores_old(qseq_r, qseq_r_comp->seqlen, bestHits_r, *best_read_score_r, best_start_pos, best_end_pos, bestTemplates_r, header_r, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
		unlock(excludeOut);
		return 1;
	}
	
	return 3;
}

int alnFragsForcePE_old(HashMapCCI **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, double mrc, int minlen, CompDNA *qseq_comp, CompDNA *qseq_r_comp, CompDNA *qseq_fr_comp, CompDNA *qseq_rr_comp, unsigned char *qseq, unsigned char *qseq_r, unsigned char *qseq_fr, unsigned char *qseq_rr, Qseqs *header, Qseqs *header_r, int kmersize, int *bestTemplates, int *bestTemplates_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int *flag, int *flag_r, int *best_read_score, int *best_read_score_r, int seq_in, long *seq_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, bestHits, aln_len, Wl, start, end, arc, rc;
	int t_len;
	double score, bestScore;
	AlnScore alnStat, alnStat_r;
	
	/* unpack qseqs */
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = qseq_comp->seqlen;
	unCompDNA(qseq_r_comp, qseq_r);
	qseq_r_comp->N[0]++;
	qseq_r_comp->N[qseq_r_comp->N[0]] = qseq_r_comp->seqlen;
	if((arc = points->len)) {
		qseq_fr_comp->N[0]++;
		qseq_fr_comp->N[qseq_fr_comp->N[0]] = qseq_fr_comp->seqlen;
		qseq_rr_comp->N[0]++;
		qseq_rr_comp->N[qseq_rr_comp->N[0]] = qseq_rr_comp->seqlen;
		points->len = 0;
	}
	Wl = -NWmatrices->rewards->Wl;
	start = 0;
	end = 0;
	score = 0;
	bestScore = 0;
	*best_read_score = 0;
	bestHits = 0;
	rc = 1;
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		/* check if index DB is loaded, and rc */
		if(template < 0) {
			if(rc) {
				qseq_comp->N[0]--;
				comp_rc(qseq_comp);
				qseq_comp->N[0]++;
				qseq_r_comp->N[0]--;
				comp_rc(qseq_r_comp);
				qseq_r_comp->N[0]++;
				
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
				
				rc = 0;
			}
		}
		
		lock(excludeDB);
		if(0 <= template && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(templates_index[template], seq_in, template_lengths[template], kmersize, seq_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(templates_index[-template], seq_in, template_lengths[-template], kmersize, seq_indexes[-template]);
		}
		unlock(excludeDB);
		
		/* align qseqs */
		template = abs(template);
		if(arc) {
			rc = anker_rc_comp(templates_index[template], qseq, qseq_fr, qseq_comp, qseq_fr_comp, 0, qseq_comp->seqlen, points);
			if(rc < 0) {
				/* rc */
				alnStat = KMA_score(templates_index[template], qseq_fr, qseq_fr_comp->seqlen, 0, qseq_fr_comp->seqlen, qseq_fr_comp, mq, scoreT, points, NWmatrices);
			} else if(rc) {
				/* forward */
				matched_templates[t_i] = -matched_templates[t_i];;
				alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, 0, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
			} else {
				alnStat.score = 0;
				alnStat.pos = 0;
				alnStat.len = 0;
				alnStat.tGaps = 0;
				alnStat.qGaps = 0;
				points->len = 0;
			}
		} else {
			alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, 0, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
		}
		
		t_len = template_lengths[abs(template)];
		if(0 < alnStat.score && minlen <= alnStat.len && mrcheck(mrc, alnStat, qseq_comp->seqlen, t_len)) {
			if(arc) {
				if(rc < 0) {
					/* rc */
					alnStat_r = KMA_score(templates_index[template], qseq_rr, qseq_rr_comp->seqlen, 0, qseq_rr_comp->seqlen, qseq_rr_comp, mq, scoreT, points, NWmatrices);
				} else if(rc) {
					/* forward */
					alnStat_r = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, 0, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
				} else {
					alnStat_r.score = 0;
					alnStat_r.pos = 0;
					alnStat_r.len = 0;
					alnStat_r.tGaps = 0;
					alnStat_r.qGaps = 0;
				}
				rc = 1;
			} else {
				alnStat_r = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, 0, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
			}
			
			/* get read score */
			if(0 < alnStat_r.score && minlen <= alnStat_r.len && mrcheck(mrc, alnStat_r, qseq_r_comp->seqlen, t_len)) {
				aln_len = alnStat.len + alnStat_r.len;
				
				/* Handle negative insertsizes caused by trimming,
				user stupidity or sample error. */
				if(alnStat.pos < alnStat_r.pos) {
					start = alnStat.pos;
					end = alnStat_r.pos + alnStat_r.len - alnStat_r.tGaps;
				} else {
					start = alnStat_r.pos;
					end = alnStat.pos + alnStat.len - alnStat.tGaps;
				}
				
				read_score = alnStat.score + alnStat_r.score;
				if(start == 0) {
					read_score += Wl;
				}
				if(end == t_len) {
					read_score += Wl;
				}
				
				score = 1.0 * read_score / aln_len;
			} else {
				read_score = 0;
			}
		} else {
			read_score = 0;
		}
		
		/* adjust to se all */
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT) {
			if(score > bestScore) { // save as best match
				bestScore = score;
				*best_read_score = read_score;
				*bestTemplates = matched_templates[t_i];
				*best_start_pos = start;
				*best_end_pos = end;
				bestHits = 1;
			} else if(score == bestScore && read_score > *best_read_score) { // save as best match
				bestScore = score;
				*best_read_score = read_score;
				*bestTemplates = matched_templates[t_i];
				*best_start_pos = start;
				*best_end_pos = end;
				bestHits = 1;
			} else if(score == bestScore && read_score == *best_read_score) { // update best match
				bestTemplates[bestHits] = matched_templates[t_i];
				best_start_pos[bestHits] = start;
				best_end_pos[bestHits] = end;
				++bestHits;
			}
		}
	}
	
	/* get best matches */
	/* get rc flag */
	if(arc) {
		rc = 0;
		for(t_i = 1; t_i <= *matched_templates && !rc; ++t_i) {
			rc = matched_templates[t_i] < 0;
		}
	}
	
	if((*best_read_score_r = *best_read_score)) {
		/* check direction of qseqs */
		if(*bestTemplates < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates[t_i] = -bestTemplates[t_i];
			}
			lock(excludeOut);
			update_Scores_pe_old(qseq_r, qseq_r_comp->seqlen, qseq, qseq_comp->seqlen, bestHits, *best_read_score, best_start_pos, best_end_pos, bestTemplates, header_r, header, *flag_r, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
		} else {
			if(!rc) {
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
				*flag ^= 48;
				*flag_r ^= 48;
			}
			lock(excludeOut);
			update_Scores_pe_old(qseq, qseq_comp->seqlen, qseq_r, qseq_r_comp->seqlen, bestHits, *best_read_score, best_start_pos, best_end_pos, bestTemplates, header, header_r, *flag, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
		}
	}
	
	return 3;
}

int alnFragsSE(HashMapCCI **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, double mrc, double minFrac, int minlen, int rc_flag, CompDNA *qseq_comp, CompDNA *qseq_r_comp, unsigned char *qseq, unsigned char *qseq_r, int q_len, int kmersize, Qseqs *header, int *bestTemplates, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *Scores, int *Lengths, int *best_start_pos, int *best_end_pos, int *flag, int *best_read_score, int seq_in, long *seq_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, bestHits, aln_len, start, end, Wl, arc, rc;
	int q_start, q_end, t_len, *qBoundPtr;
	double score, bestScore;
	AlnScore alnStat;
	
	/* reverse complement seq */
	if(rc_flag < 0) {
		if(qseq_r_comp->size < qseq_comp->size) {
			freeComp(qseq_r_comp);
			allocComp(qseq_r_comp, qseq_comp->size);
		}
		rc_comp(qseq_comp, qseq_r_comp);
		unCompDNA(qseq_r_comp, qseq_r);
		qseq_r_comp->N[0]++;
		qseq_r_comp->N[qseq_r_comp->N[0]] = q_len;
	}
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = q_len;
	
	bestScore = 0;
	*best_read_score = 0;
	bestHits = 0;
	Wl = -NWmatrices->rewards->Wl;
	arc = points->len;
	points->len = 0;
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		/* check if index DB is loaded */
		lock(excludeDB);
		if(template >= 0 && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(templates_index[template], seq_in, template_lengths[template], kmersize, seq_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(templates_index[-template], seq_in, template_lengths[-template], kmersize, seq_indexes[-template]);
		}
		unlock(excludeDB);
		
		/* q-bound */
		if(2 * sizeof(int) + 1 < header->len && header->seq[header->len - 2 * sizeof(int) - 1] == 0) {
			qBoundPtr = (int*) (header->seq + (header->len - 2 * sizeof(int)));
			q_start = *qBoundPtr;
			q_end = *++qBoundPtr;
		} else {
			q_start = 0;
			q_end = q_len;
		}
		
		/* check ankers */
		if(arc) {
			rc = anker_rc_comp(templates_index[abs(template)], qseq, qseq_r, qseq_comp, qseq_r_comp, q_start, q_end, points);
			if(rc < 0) {
				if(0 < template) {
					template = -template;
				}
				alnStat = KMA_score(templates_index[-template], qseq_r, q_len, q_len - q_end, q_len - q_start, qseq_r_comp, mq, scoreT, points, NWmatrices);
			} else if(rc) {
				if(template < 0) {
					template = -template;
				}
				alnStat = KMA_score(templates_index[template], qseq, q_len, q_start, q_end, qseq_comp, mq, scoreT, points, NWmatrices);
			} else {
				alnStat.score = 0;
				alnStat.pos = 0;
				alnStat.len = 0;
				alnStat.tGaps = 0;
				alnStat.qGaps = 0;
				points->len = 0;
			}
		} else {
			if(template < 0) {
				alnStat = KMA_score(templates_index[-template], qseq_r, q_len, q_len - q_end, q_len - q_start, qseq_r_comp, mq, scoreT, points, NWmatrices);
			} else {
				alnStat = KMA_score(templates_index[template], qseq, q_len, q_start, q_end, qseq_comp, mq, scoreT, points, NWmatrices);
			}
		}
		
		/* get read score */
		aln_len = alnStat.len;
		start = alnStat.pos;
		end = start + aln_len - alnStat.tGaps;
		t_len = template_lengths[abs(template)];
		if(t_len < end) {
			/* circular aln */
			end -= t_len;
		}
		
		/* penalty for non complete mapping */
		read_score = alnStat.score;
		/* full gene award */
		if(q_len <= aln_len || t_len <= aln_len) {
			read_score += Wl;
		} else {
			if(start == 0) {
				read_score += Wl;
			}
			if(end == t_len) {
				read_score += Wl;
			}
		}
		//read_score += (((start != 0) + (end != template_lengths[abs(template)])) * Wl);
		
		/* Get normed score check read coverage */
		if(minlen <= aln_len && mrcheck(mrc, alnStat, q_len, t_len)) {
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
			score = 0;
		}
		
		/* save best match(es) */
		if(kmersize < read_score && scoreT <= score) {
			bestTemplates[bestHits] = template;
			best_start_pos[bestHits] = start;
			best_end_pos[bestHits] = end;
			Scores[bestHits] = read_score;
			Lengths[bestHits] = aln_len;
			++bestHits;
			
			/* save best scores */
			/*
			if(*best_read_score < read_score) {
				bestScore = score;
				*best_read_score = read_score;
				bestLen = aln_len;
			} else if(*best_read_score == read_score && bestScore < score) {
				bestScore = score;
				bestLen = aln_len;
			}
			*/
			
			if(bestScore < score) {
				bestScore = score;
			}
			if(*best_read_score < read_score) {
				*best_read_score = read_score;
			}
			
			/* org */
			/*if(bestScore < score) {
				bestScore = score;
				*best_read_score = read_score;
				bestLen = aln_len;
			} else if(bestScore == score && *best_read_score < read_score) {
				*best_read_score = read_score;
				bestLen = aln_len;
			}*/
		}
	}
	
	if(*best_read_score > kmersize) {
		lock(excludeOut);
		//*matched_templates = update_Scores(qseq, q_len, minFrac, bestHits, *best_read_score, bestLen, best_start_pos, best_end_pos, bestTemplates, Scores, Lengths, header, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
		*matched_templates = update_Scores(qseq, q_len, minFrac, bestHits, *best_read_score, bestScore, best_start_pos, best_end_pos, bestTemplates, Scores, Lengths, header, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
		unlock(excludeOut);
		return 0;
	} else {
		*flag |= 4;
	}
	
	return 1;
}

int alnFragsUnionPE(HashMapCCI **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, double mrc, double minFrac, int minlen, CompDNA *qseq_comp, CompDNA *qseq_r_comp, CompDNA *qseq_fr_comp, CompDNA *qseq_rr_comp, unsigned char *qseq, unsigned char *qseq_r, unsigned char *qseq_fr, unsigned char *qseq_rr, Qseqs *header, Qseqs *header_r, int kmersize, int *bestTemplates, int *bestTemplates_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int *flag, int *flag_r, int *best_read_score, int *best_read_score_r, int seq_in, long *seq_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, Wl;
	int compScore, bestHits, bestHits_r, aln_len, start, end, t_len, arc, rc;
	double score, score_r;
	AlnScore alnStat;
	
	/* unpack qseqs */
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = qseq_comp->seqlen;
	unCompDNA(qseq_r_comp, qseq_r);
	qseq_r_comp->N[0]++;
	qseq_r_comp->N[qseq_r_comp->N[0]] = qseq_r_comp->seqlen;
	if((arc = points->len)) {
		qseq_fr_comp->N[0]++;
		qseq_fr_comp->N[qseq_fr_comp->N[0]] = qseq_fr_comp->seqlen;
		qseq_rr_comp->N[0]++;
		qseq_rr_comp->N[qseq_rr_comp->N[0]] = qseq_rr_comp->seqlen;
		points->len = 0;
	}
	Wl = -NWmatrices->rewards->Wl;
	start = 0;
	end = 0;
	score = 0;
	*best_read_score = 0;
	*best_read_score_r = 0;
	compScore = 0;
	bestHits = 0;
	rc = 1;
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		/* check if index DB is loaded, and rc */
		if(template < 0) {
			if(rc) {
				qseq_comp->N[0]--;
				comp_rc(qseq_comp);
				qseq_comp->N[0]++;
				qseq_r_comp->N[0]--;
				comp_rc(qseq_r_comp);
				qseq_r_comp->N[0]++;
				
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
				
				rc = 0;
			}
		}
		
		lock(excludeDB);
		if(0 <= template && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(templates_index[template], seq_in, template_lengths[template], kmersize, seq_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(templates_index[-template], seq_in, template_lengths[-template], kmersize, seq_indexes[-template]);
		}
		unlock(excludeDB);
		
		/* align qseqs */
		template = abs(template);
		if(arc) {
			rc = anker_rc_comp(templates_index[template], qseq, qseq_fr, qseq_comp, qseq_fr_comp, 0, qseq_comp->seqlen,  points);
			if(rc < 0) {
				/* rc */
				alnStat = KMA_score(templates_index[template], qseq_fr, qseq_fr_comp->seqlen, 0, qseq_fr_comp->seqlen, qseq_fr_comp, mq, scoreT, points, NWmatrices);
			} else if(rc) {
				/* forward */
				matched_templates[t_i] = -matched_templates[t_i];
				alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, 0, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
			} else {
				alnStat.score = 0;
				alnStat.pos = 0;
				alnStat.len = 0;
				alnStat.tGaps = 0;
				alnStat.qGaps = 0;
				points->len = 0;
			}
		} else {
			alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, 0, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
		}
		
		/* get read score */
		aln_len = alnStat.len;
		read_score = alnStat.score;
		t_len = template_lengths[abs(template)];
		if(minlen <= aln_len && 0 < read_score && mrcheck(mrc, alnStat, qseq_comp->seqlen, t_len)) {
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.tGaps;
			if(start == 0) {
				read_score += Wl;
			}
			if(end == t_len) {
				read_score += Wl;
			}
			
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* save match(es) */
		if(read_score > kmersize && score >= scoreT) {
			bestTemplates[t_i] = read_score;
			best_start_pos[t_i] = start;
			best_end_pos[t_i] = end;
			if(*best_read_score < read_score) {
				*best_read_score = read_score;
			}
		} else {
			bestTemplates[t_i] = 0;
			best_start_pos[t_i] = -1;
			best_end_pos[t_i] = -1;
		}
		
		/* align qseqs */
		template = abs(template);
		if(arc) {
			if(rc < 0) {
				/* rc */
				alnStat = KMA_score(templates_index[template], qseq_rr, qseq_rr_comp->seqlen, 0, qseq_rr_comp->seqlen, qseq_rr_comp, mq, scoreT, points, NWmatrices);
			} else if(rc) {
				/* forward */
				alnStat = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, 0, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
			} else {
				alnStat.score = 0;
				alnStat.pos = 0;
				alnStat.len = 0;
				alnStat.tGaps = 0;
				alnStat.qGaps = 0;
			}
			rc = 1;
		} else {
			alnStat = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, 0, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
		}
		
		/* get read score */
		aln_len = alnStat.len;
		read_score = alnStat.score;
		if(minlen <= aln_len && 0 < read_score && mrcheck(mrc, alnStat, qseq_r_comp->seqlen, t_len)) {
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.tGaps;
			if(start == 0) {
				read_score += Wl;
			}
			if(end == t_len) {
				read_score += Wl;
			}
			
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* save match(es) */
		if(read_score > kmersize && score >= scoreT) {
			bestTemplates_r[t_i] = read_score;
			if(bestTemplates[t_i]) {
				/* Handle negative insertsizes caused by trimming,
				user stupidity or sample error. */
				if(start < best_start_pos[t_i]) {
					best_start_pos[t_i] = start;
				} else {
					best_end_pos[t_i] = end;
				}
			} else {
				best_start_pos[t_i] = start;
				best_end_pos[t_i] = end;
			}
			if(*best_read_score_r < read_score) {
				*best_read_score_r = read_score;
			}
		} else {
			bestTemplates_r[t_i] = 0;
		}
		
		read_score += bestTemplates[t_i];
		if(compScore < read_score) {
			compScore = read_score;
		}
	}
	
	if(*best_read_score || *best_read_score_r) {
		/* get rc flag */
		if(arc) {
			rc = 0;
			for(t_i = 1; t_i <= *matched_templates && !rc; ++t_i) {
				rc = matched_templates[t_i] < 0;
			}
		}
		
		/* check pairing */
		bestHits = 0;
		if(*best_read_score && *best_read_score_r) {
			/* check union */
			score = minFrac < 0 ? -minFrac : minFrac;
			score_r = *best_read_score_r * score;
			score *= *best_read_score;
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(score <= bestTemplates[t_i] && score_r <= bestTemplates_r[t_i]) {
					bestTemplates_r[bestHits] = bestTemplates[t_i] + bestTemplates_r[t_i];
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				}
			}
			
			if(bestHits) {
				/* proper pair */
				/* check direction of qseqs */
				if(*bestTemplates < 0) {
					for(t_i = 0; t_i < bestHits; ++t_i) {
						bestTemplates[t_i] = -bestTemplates[t_i];
					}
					lock(excludeOut);
					*matched_templates = update_Scores_pe(qseq_r, qseq_r_comp->seqlen, qseq, qseq_comp->seqlen, minFrac, bestHits, *best_read_score + *best_read_score_r, best_start_pos, best_end_pos, bestTemplates, bestTemplates_r, header_r, header, *flag_r, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
					unlock(excludeOut);
				} else {
					if(!rc) {
						strrc(qseq, qseq_comp->seqlen);
						strrc(qseq_r, qseq_r_comp->seqlen);
						*flag ^= 48;
						*flag_r ^= 48;
					}
					lock(excludeOut);
					*matched_templates = update_Scores_pe(qseq, qseq_comp->seqlen, qseq_r, qseq_r_comp->seqlen, minFrac, bestHits, *best_read_score + *best_read_score_r, best_start_pos, best_end_pos, bestTemplates, bestTemplates_r, header, header_r, *flag, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
					unlock(excludeOut);
				}
			} else {
				/* unmaided pair */
				bestHits_r = 0;
				
				/* move reverse up and forward down */
				/*
				1. match is forward, move down (no overwrite of reverse) -> increase counter
				2. match is reverse, swap up (new entry might be reverse or forward -> decrease end
				3. neither is a match -> increase counter.
				*/
				score = minFrac < 0 ? -minFrac : minFrac;
				score_r = *best_read_score_r * score;
				score *= *best_read_score;
				
				/* compute scores */
				t_i = 1;
				end = *matched_templates;
				while(t_i <= end) {
					if(score <= bestTemplates[t_i]) {
						matched_templates[bestHits] = matched_templates[t_i];
						bestTemplates[bestHits] = bestTemplates[t_i];
						best_start_pos[bestHits] = best_start_pos[t_i];
						best_end_pos[bestHits] = best_end_pos[t_i];
						++bestHits;
						++t_i;
					} else if(score_r <= bestTemplates_r[t_i]) {
						exchange(matched_templates[t_i], matched_templates[end], template);
						exchange(bestTemplates_r[t_i], bestTemplates_r[end], template);
						exchange(best_start_pos[t_i], best_start_pos[end], template);
						exchange(best_end_pos[t_i], best_end_pos[end], template);
						++bestHits_r;
						--end;
					} else {
						++t_i;
					}
				}
				bestTemplates_r += end;
				
				/* check direction of qseqs */
				if(*bestTemplates < 0) {
					for(t_i = 0; t_i < bestHits; ++t_i) {
						bestTemplates[t_i] = -bestTemplates[t_i];
					}
				} else if(!rc) {
					strrc(qseq, qseq_comp->seqlen);
					*flag ^= 16;
					*flag_r ^= 32;
				}
				if(*bestTemplates_r < 0) {
					for(t_i = 0; t_i < bestHits_r; ++t_i) {
						bestTemplates_r[t_i] = -bestTemplates_r[t_i];
					}
				} else if(!rc) {
					strrc(qseq_r, qseq_r_comp->seqlen);
					*flag ^= 32;
					*flag_r ^= 16;
				}
				if(*flag & 2) {
					*flag ^= 2;
					*flag_r ^= 2;
				}
				lock(excludeOut);
				/*
				bestLen,
				Lengths
				*/
				*matched_templates = update_Scores_se(qseq, qseq_comp->seqlen, minFrac, bestHits, *best_read_score, best_start_pos, best_end_pos, matched_templates, bestTemplates, header, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
				*matched_templates += update_Scores_se(qseq_r, qseq_r_comp->seqlen, minFrac, bestHits_r, *best_read_score_r, best_start_pos + end, best_end_pos + end, matched_templates + end, bestTemplates_r, header_r, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
				unlock(excludeOut);
			}
			return 0;
		} else if(*best_read_score) {
			/* forward only */
			/* compute scores */
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(bestTemplates[t_i]) {
					bestTemplates_r[bestHits] = bestTemplates[t_i];
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				}
			}
			
			/* check direction of qseq */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
			} else if(!rc) {
				strrc(qseq, qseq_comp->seqlen);
				*flag ^= 16;
				*flag_r ^= 32;
			}
			*flag |= 8;
			*flag_r ^= 4;
			if(*flag & 2) {
				*flag ^= 2;
				*flag_r ^= 2;
			}
			
			lock(excludeOut);
			*matched_templates = update_Scores_se(qseq, qseq_comp->seqlen, minFrac, bestHits, *best_read_score, best_start_pos, best_end_pos, bestTemplates, bestTemplates_r, header, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
			
			return 2;
		} else {
			/* reverse only */
			/* compute scores */
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(bestTemplates_r[t_i]) {
					bestTemplates_r[bestHits] = bestTemplates_r[t_i];
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				}
			}
			/* check direction of qseqs */
			if(*bestTemplates_r < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates_r[t_i] = -bestTemplates_r[t_i];
				}
			} else if(!rc) {
				strrc(qseq_r, qseq_r_comp->seqlen);
				*flag ^= 32;
				*flag_r ^= 16;
			}
			*flag_r |= 8;
			*flag ^= 4;
			if(*flag_r & 2) {
				*flag ^= 2;
				*flag_r ^= 2;
			}
			lock(excludeOut);
			/*
			bestLen,
			Lengths
			*/
			*matched_templates = update_Scores_se(qseq_r, qseq_r_comp->seqlen, minFrac, bestHits, *best_read_score_r, best_start_pos, best_end_pos, bestTemplates, bestTemplates_r, header_r, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
			
			return 1;
		}
	}
	
	return 3;
}

int alnFragsPenaltyPE(HashMapCCI **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, double mrc, double minFrac, int minlen, CompDNA *qseq_comp, CompDNA *qseq_r_comp, CompDNA *qseq_fr_comp, CompDNA *qseq_rr_comp, unsigned char *qseq, unsigned char *qseq_r, unsigned char *qseq_fr, unsigned char *qseq_rr, Qseqs *header, Qseqs *header_r, int kmersize, int *bestTemplates, int *bestTemplates_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int *flag, int *flag_r, int *best_read_score, int *best_read_score_r, int seq_in, long *seq_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, Wl, PE;
	int compScore, bestHits, bestHits_r, aln_len, start, end, t_len, arc, rc;
	double score, score_r;
	AlnScore alnStat;
	
	/* unpack qseqs */
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = qseq_comp->seqlen;
	unCompDNA(qseq_r_comp, qseq_r);
	qseq_r_comp->N[0]++;
	qseq_r_comp->N[qseq_r_comp->N[0]] = qseq_r_comp->seqlen;
	if((arc = points->len)) {
		qseq_fr_comp->N[0]++;
		qseq_fr_comp->N[qseq_fr_comp->N[0]] = qseq_fr_comp->seqlen;
		qseq_rr_comp->N[0]++;
		qseq_rr_comp->N[qseq_rr_comp->N[0]] = qseq_rr_comp->seqlen;
		points->len = 0;
	}
	Wl = -NWmatrices->rewards->Wl;
	PE = NWmatrices->rewards->PE;
	start = 0;
	end = 0;
	score = 0;
	*best_read_score = 0;
	*best_read_score_r = 0;
	compScore = 0;
	bestHits = 0;
	rc = 1;
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		/* check if index DB is loaded, and rc */
		if(template < 0) {
			if(rc) {
				qseq_comp->N[0]--;
				comp_rc(qseq_comp);
				qseq_comp->N[0]++;
				qseq_r_comp->N[0]--;
				comp_rc(qseq_r_comp);
				qseq_r_comp->N[0]++;
				
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
				
				rc = 0;
			}
		}
		
		lock(excludeDB);
		if(0 <= template && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(templates_index[template], seq_in, template_lengths[template], kmersize, seq_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(templates_index[-template], seq_in, template_lengths[-template], kmersize, seq_indexes[-template]);
		}
		unlock(excludeDB);
		
		/* align qseqs */
		template = abs(template);
		if(arc) {
			rc = anker_rc_comp(templates_index[template], qseq, qseq_fr, qseq_comp, qseq_fr_comp, 0, qseq_comp->seqlen, points);
			if(rc < 0) {
				/* rc */
				alnStat = KMA_score(templates_index[template], qseq_fr, qseq_fr_comp->seqlen, 0, qseq_fr_comp->seqlen, qseq_fr_comp, mq, scoreT, points, NWmatrices);
			} else if(rc) {
				/* forward */
				matched_templates[t_i] = -matched_templates[t_i];;
				alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, 0, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
			} else {
				alnStat.score = 0;
				alnStat.pos = 0;
				alnStat.len = 0;
				alnStat.tGaps = 0;
				alnStat.qGaps = 0;
				points->len = 0;
			}
		} else {
			alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, 0, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
		}
		
		/* get read score */
		aln_len = alnStat.len;
		read_score = alnStat.score;
		t_len = template_lengths[abs(template)];
		if(minlen <= aln_len && 0 < read_score && mrcheck(mrc, alnStat, qseq_comp->seqlen, t_len)) {
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.tGaps;
			if(start == 0) {
				read_score += Wl;
			}
			if(end == t_len) {
				read_score += Wl;
			}
			
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT) {
			bestTemplates[t_i] = read_score;
			best_start_pos[t_i] = start;
			best_end_pos[t_i] = end;
			if(*best_read_score < read_score) {
				*best_read_score = read_score;
			}
		} else {
			bestTemplates[t_i] = 0;
			best_start_pos[t_i] = -1;
			best_end_pos[t_i] = -1;
		}
		
		/* align qseqs */
		template = abs(template);
		if(arc) {
			if(rc < 0) {
				/* rc */
				alnStat = KMA_score(templates_index[template], qseq_rr, qseq_rr_comp->seqlen, 0, qseq_rr_comp->seqlen, qseq_rr_comp, mq, scoreT, points, NWmatrices);
			} else if(rc) {
				/* forward */
				alnStat = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, 0, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
			} else {
				alnStat.score = 0;
				alnStat.pos = 0;
				alnStat.len = 0;
				alnStat.tGaps = 0;
				alnStat.qGaps = 0;
			}
			rc = 1;
		} else {
			alnStat = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, 0, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
		}
		
		/* get read score */
		aln_len = alnStat.len;
		read_score = alnStat.score;
		if(minlen <= aln_len && 0 < read_score && mrcheck(mrc, alnStat, qseq_r_comp->seqlen, t_len)) {
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.tGaps;
			if(start == 0) {
				read_score += Wl;
			}
			if(end == t_len) {
				read_score += Wl;
			}
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* save match(es) */
		if(read_score > kmersize && score >= scoreT) {
			bestTemplates_r[t_i] = read_score;
			if(bestTemplates[t_i]) {
				/* Handle negative insertsizes caused by trimming,
				user stupidity or sample error. */
				if(start < best_start_pos[t_i]) {
					best_start_pos[t_i] = start;
				} else {
					best_end_pos[t_i] = end;
				}
			} else {
				best_start_pos[t_i] = start;
				best_end_pos[t_i] = end;
			}
			if(*best_read_score_r < read_score) {
				*best_read_score_r = read_score;
			}
		} else {
			bestTemplates_r[t_i] = 0;
		}
		
		read_score += bestTemplates[t_i];
		if(compScore < read_score) {
			compScore = read_score;
		}
	}
	
	if(*best_read_score || *best_read_score_r) {
		/* get rc flag */
		if(arc) {
			rc = 0;
			for(t_i = 1; t_i <= *matched_templates && !rc; ++t_i) {
				rc = matched_templates[t_i] < 0;
			}
		}
		
		/* check pairing */
		bestHits = 0;
		score = minFrac < 0 ? -minFrac : minFrac;
		score *= (*best_read_score + *best_read_score_r);
		if(compScore && score <= (compScore + PE)) {
			/* proper pair */
			*best_read_score = compScore + PE;
			*best_read_score_r = *best_read_score;
			
			/* compute scores */
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(bestTemplates[t_i] && bestTemplates_r[t_i]) {
					bestTemplates_r[bestHits] = bestTemplates[t_i] + bestTemplates_r[t_i] + PE;
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				}
			}
			
			/* check direction of qseqs */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
				lock(excludeOut);
				*matched_templates = update_Scores_pe(qseq_r, qseq_r_comp->seqlen, qseq, qseq_comp->seqlen, minFrac, bestHits, *best_read_score, best_start_pos, best_end_pos, bestTemplates, bestTemplates_r, header_r, header, *flag_r, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
				unlock(excludeOut);
			} else {
				if(!rc) {
					strrc(qseq, qseq_comp->seqlen);
					strrc(qseq_r, qseq_r_comp->seqlen);
					*flag ^= 48;
					*flag_r ^= 48;
				}
				lock(excludeOut);
				*matched_templates = update_Scores_pe(qseq, qseq_comp->seqlen, qseq_r, qseq_r_comp->seqlen, minFrac, bestHits, *best_read_score, best_start_pos, best_end_pos, bestTemplates, bestTemplates_r, header, header_r, *flag, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
				unlock(excludeOut);
			}
			return 0;
		} else if(*best_read_score && *best_read_score_r) {
			/* unmaided pair */
			bestHits_r = 0;
			
			/* move reverse up and forward down */
			/*
			1. match is forward, move down (no overwrite of reverse) -> increase counter
			2. match is reverse, swap up (new entry might be reverse or forward -> decrease end
			3. neither is a match -> increase counter.
			*/
			score = minFrac < 0 ? -minFrac : minFrac;
			score_r = *best_read_score_r * score;
			score *= *best_read_score;
			
			/* compute scores */
			t_i = 1;
			end = *matched_templates;
			while(t_i <= end) {
				if(score <= bestTemplates[t_i]) {
					matched_templates[bestHits] = matched_templates[t_i];
					bestTemplates[bestHits] = bestTemplates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
					++t_i;
				} else if(score_r <= bestTemplates_r[t_i]) {
					exchange(matched_templates[t_i], matched_templates[end], template);
					exchange(bestTemplates_r[t_i], bestTemplates_r[end], template);
					exchange(best_start_pos[t_i], best_start_pos[end], template);
					exchange(best_end_pos[t_i], best_end_pos[end], template);
					++bestHits_r;
					--end;
				} else {
					++t_i;
				}
			}
			bestTemplates_r += end;
			
			/* check direction of qseqs */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
			} else if(!rc) {
				strrc(qseq, qseq_comp->seqlen);
				*flag ^= 16;
				*flag_r ^= 32;
			}
			if(*bestTemplates_r < 0) {
				for(t_i = 0; t_i < bestHits_r; ++t_i) {
					bestTemplates_r[t_i] = -bestTemplates_r[t_i];
				}
			} else if(!rc) {
				strrc(qseq_r, qseq_r_comp->seqlen);
				*flag ^= 32;
				*flag_r ^= 16;
			}
			if(*flag & 2) {
				*flag ^= 2;
				*flag_r ^= 2;
			}
			lock(excludeOut);
			/*
			bestLen,
			Lengths
			*/
			*matched_templates = update_Scores_se(qseq, qseq_comp->seqlen, minFrac, bestHits, *best_read_score, best_start_pos, best_end_pos, matched_templates, bestTemplates, header, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
			*matched_templates += update_Scores_se(qseq_r, qseq_r_comp->seqlen, minFrac, bestHits_r, *best_read_score_r, best_start_pos + end, best_end_pos + end, matched_templates + end, bestTemplates_r, header_r, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
			
			return 0;
		} else if(*best_read_score) {
			/* forward only */
			/* compute scores */
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(bestTemplates[t_i]) {
					bestTemplates_r[bestHits] = bestTemplates[t_i];
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				}
			}
			
			/* check direction of qseq */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
			} else if(!rc) {
				strrc(qseq, qseq_comp->seqlen);
				*flag ^= 16;
				*flag_r ^= 32;
			}
			*flag |= 8;
			*flag_r ^= 4;
			if(*flag & 2) {
				*flag ^= 2;
				*flag_r ^= 2;
			}
			
			lock(excludeOut);
			*matched_templates = update_Scores_se(qseq, qseq_comp->seqlen, minFrac, bestHits, *best_read_score, best_start_pos, best_end_pos, bestTemplates, bestTemplates_r, header, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
			
			return 2;
		} else {
			/* reverse only */
			/* compute scores */
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(bestTemplates_r[t_i]) {
					bestTemplates_r[bestHits] = bestTemplates_r[t_i];
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				}
			}
			/* check direction of qseqs */
			if(*bestTemplates_r < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates_r[t_i] = -bestTemplates_r[t_i];
				}
			} else if(!rc) {
				strrc(qseq_r, qseq_r_comp->seqlen);
				*flag ^= 32;
				*flag_r ^= 16;
			}
			*flag_r |= 8;
			*flag ^= 4;
			if(*flag_r & 2) {
				*flag ^= 2;
				*flag_r ^= 2;
			}
			lock(excludeOut);
			/*
			bestLen,
			Lengths
			*/
			*matched_templates = update_Scores_se(qseq_r, qseq_r_comp->seqlen, minFrac, bestHits, *best_read_score_r, best_start_pos, best_end_pos, bestTemplates, bestTemplates_r, header_r, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
			
			return 1;
		}
	}
	
	return 3;
}

int alnFragsForcePE(HashMapCCI **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, double mrc, double minFrac, int minlen, CompDNA *qseq_comp, CompDNA *qseq_r_comp, CompDNA *qseq_fr_comp, CompDNA *qseq_rr_comp, unsigned char *qseq, unsigned char *qseq_r, unsigned char *qseq_fr, unsigned char *qseq_rr, Qseqs *header, Qseqs *header_r, int kmersize, int *bestTemplates, int *bestTemplates_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int *flag, int *flag_r, int *best_read_score, int *best_read_score_r, int seq_in, long *seq_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, bestHits, aln_len, Wl, start, end, arc, rc;
	int t_len;
	double score;
	AlnScore alnStat, alnStat_r;
	
	/* unpack qseqs */
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = qseq_comp->seqlen;
	unCompDNA(qseq_r_comp, qseq_r);
	qseq_r_comp->N[0]++;
	qseq_r_comp->N[qseq_r_comp->N[0]] = qseq_r_comp->seqlen;
	if((arc = points->len)) {
		qseq_fr_comp->N[0]++;
		qseq_fr_comp->N[qseq_fr_comp->N[0]] = qseq_fr_comp->seqlen;
		qseq_rr_comp->N[0]++;
		qseq_rr_comp->N[qseq_rr_comp->N[0]] = qseq_rr_comp->seqlen;
		points->len = 0;
	}
	Wl = -NWmatrices->rewards->Wl;
	start = 0;
	end = 0;
	score = 0;
	*best_read_score = 0;
	bestHits = 0;
	rc = 1;
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		/* check if index DB is loaded, and rc */
		if(template < 0) {
			if(rc) {
				qseq_comp->N[0]--;
				comp_rc(qseq_comp);
				qseq_comp->N[0]++;
				qseq_r_comp->N[0]--;
				comp_rc(qseq_r_comp);
				qseq_r_comp->N[0]++;
				
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
				
				rc = 0;
			}
		}
		
		lock(excludeDB);
		if(0 <= template && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(templates_index[template], seq_in, template_lengths[template], kmersize, seq_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(templates_index[-template], seq_in, template_lengths[-template], kmersize, seq_indexes[-template]);
		}
		unlock(excludeDB);
		
		/* align qseqs */
		template = abs(template);
		if(arc) {
			rc = anker_rc_comp(templates_index[template], qseq, qseq_fr, qseq_comp, qseq_fr_comp, 0, qseq_comp->seqlen, points);
			if(rc < 0) {
				/* rc */
				alnStat = KMA_score(templates_index[template], qseq_fr, qseq_fr_comp->seqlen, 0, qseq_fr_comp->seqlen, qseq_fr_comp, mq, scoreT, points, NWmatrices);
			} else if(rc) {
				/* forward */
				matched_templates[t_i] = -matched_templates[t_i];;
				alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, 0, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
			} else {
				alnStat.score = 0;
				alnStat.pos = 0;
				alnStat.len = 0;
				alnStat.tGaps = 0;
				alnStat.qGaps = 0;
				points->len = 0;
			}
		} else {
			alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, 0, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
		}
		
		t_len = template_lengths[abs(template)];
		if(0 < alnStat.score && minlen <= alnStat.len && mrcheck(mrc, alnStat, qseq_comp->seqlen, t_len)) {
			if(arc) {
				if(rc < 0) {
					/* rc */
					alnStat_r = KMA_score(templates_index[template], qseq_rr, qseq_rr_comp->seqlen, 0, qseq_rr_comp->seqlen, qseq_rr_comp, mq, scoreT, points, NWmatrices);
				} else if(rc) {
					/* forward */
					alnStat_r = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, 0, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
				} else {
					alnStat_r.score = 0;
					alnStat_r.pos = 0;
					alnStat_r.len = 0;
					alnStat_r.tGaps = 0;
					alnStat_r.qGaps = 0;
				}
				rc = 1;
			} else {
				alnStat_r = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, 0, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
			}
			
			/* get read score */
			if(0 < alnStat_r.score && minlen <= alnStat_r.len && mrcheck(mrc, alnStat_r, qseq_r_comp->seqlen, t_len)) {
				aln_len = alnStat.len + alnStat_r.len;
				
				/* Handle negative insertsizes caused by trimming,
				user stupidity or sample error. */
				if(alnStat.pos < alnStat_r.pos) {
					start = alnStat.pos;
					end = alnStat_r.pos + alnStat_r.len - alnStat_r.tGaps;
				} else {
					start = alnStat_r.pos;
					end = alnStat.pos + alnStat.len - alnStat.tGaps;
				}
				
				read_score = alnStat.score + alnStat_r.score;
				if(start == 0) {
					read_score += Wl;
				}
				if(end == t_len) {
					read_score += Wl;
				}
				
				score = 1.0 * read_score / aln_len;
			} else {
				read_score = 0;
			}
		} else {
			read_score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT) {
			bestTemplates[bestHits] = template;
			best_start_pos[bestHits] = start;
			best_end_pos[bestHits] = end;
			bestTemplates_r[bestHits] = read_score;
			++bestHits;
			
			if(*best_read_score < read_score) {
				*best_read_score = read_score;
			}
		}
	}
	
	if((*best_read_score_r = *best_read_score)) {
		/* check direction of qseqs */
		if(*bestTemplates < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates[t_i] = -bestTemplates[t_i];
			}
			lock(excludeOut);
			*matched_templates = update_Scores_pe(qseq_r, qseq_r_comp->seqlen, qseq, qseq_comp->seqlen, minFrac, bestHits, *best_read_score, best_start_pos, best_end_pos, bestTemplates, bestTemplates_r, header_r, header, *flag_r, *flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
		} else {
			/* get rc flag */
			if(arc) {
				rc = 0;
				for(t_i = 0; t_i < bestHits && !rc; ++t_i) {
					rc = bestTemplates[t_i] < 0;
				}
			}
			if(!rc) {
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
				*flag ^= 48;
				*flag_r ^= 48;
			}
			lock(excludeOut);
			*matched_templates = update_Scores_pe(qseq, qseq_comp->seqlen, qseq_r, qseq_r_comp->seqlen, minFrac, bestHits, *best_read_score, best_start_pos, best_end_pos, bestTemplates, bestTemplates_r, header, header_r, *flag, *flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
		}
		return 0;
	}
	
	return 3;
}

void * alnFrags_threaded(void * arg) {
	
	static volatile int Lock[3] = {0, 0, 0};
	volatile int *excludeIn = &Lock[0], *excludeOut = &Lock[1], *excludeDB = &Lock[2];
	Aln_thread *thread = arg;
	int rc_flag, read_score, delta, seq_in, kmersize, minlen;
	int flag, flag_r, mq, sam, unmapped, best_read_score, stats[2];
	int *matched_templates, *bestTemplates, *bestTemplates_r;
	int *template_lengths, *best_start_pos, *best_end_pos, *Lengths;
	long *seq_indexes;
	long unsigned *alignment_scores, *uniq_alignment_scores;
	unsigned char *qseq_fr, *qseq_rr;
	double scoreT, mrc, minFrac;
	FILE *inputfile, *frag_out_raw;
	FileBuff *frag_out_all;
	CompDNA *qseq_comp, *qseq_r_comp, *qseq_fr_comp, *qseq_rr_comp;
	Qseqs *qseq, *qseq_r, *header, *header_r;
	AlnPoints *points;
	NWmat *NWmatrices;
	HashMapCCI **templates_index;
	
	/* get input */
	matched_templates = thread->matched_templates;
	bestTemplates = thread->bestTemplates;
	bestTemplates_r = thread->bestTemplates_r;
	best_start_pos = thread->best_start_pos;
	best_end_pos = thread->best_end_pos;
	Lengths = thread->Lengths;
	alignment_scores = thread->alignment_scores;
	uniq_alignment_scores = thread->uniq_alignment_scores;
	seq_indexes = thread->seq_indexes;
	inputfile = thread->inputfile;
	frag_out_raw = thread->frag_out_raw;
	frag_out_all = thread->frag_out_all;
	seq_in = thread->seq_in;
	qseq_comp = thread->qseq_comp;
	qseq_r_comp = thread->qseq_r_comp;
	qseq = thread->qseq;
	qseq_r = thread->qseq_r;
	header = thread->header;
	header_r = thread->header_r;
	points = thread->points;
	NWmatrices = thread->NWmatrices;
	kmersize = thread->kmersize;
	minlen = thread->minlen;
	mq = thread->mq;
	sam = thread->sam;
	scoreT = thread->scoreT;
	mrc = thread->mrc;
	minFrac = thread->minFrac;
	template_lengths = thread->template_lengths;
	templates_index = thread->templates_index;
	
	qseq_fr_comp = setComp(32);
	qseq_rr_comp = setComp(32);
	qseq_fr = smalloc(32);
	qseq_rr = smalloc(32);
	delta = qseq->size;
	read_score = 0;
	stats[0] = 0;
	//lock(excludeIn);
	lockTime(excludeIn, 65536);
	while((rc_flag = get_ankers(matched_templates, qseq_comp, header, &flag, inputfile)) != 0) {
		points->len = rc_flag < 0;
		if(*matched_templates) { // SE
			read_score = 0;
			qseq_r->len = 0;
		} else { // PE
			read_score = get_ankers(matched_templates, qseq_r_comp, header_r, &flag_r, inputfile);
			read_score = labs(read_score);
			qseq_r->len = qseq_r_comp->seqlen;
		}
		unlock(excludeIn);
		qseq->len = qseq_comp->seqlen;
		
		if(delta <= MAX(qseq->len, qseq_r->len)) {
			delta = MAX(qseq->len, qseq_r->len);
			delta <<= 1;
			qseq->size = delta;
			qseq_r->size = delta;
			free(qseq->seq);
			free(qseq_r->seq);
			qseq->seq = smalloc(delta);
			qseq_r->seq = smalloc(delta);
		}
		
		if(points->len && read_score) {
			if(qseq_fr_comp->size < delta) {
				free(qseq_fr);
				free(qseq_rr);
				qseq_fr = smalloc(delta);
				qseq_rr = smalloc(delta);
				dallocComp(qseq_fr_comp, delta);
				dallocComp(qseq_rr_comp, delta);
			}
			rc_comp(qseq_comp, qseq_fr_comp);
			rc_comp(qseq_r_comp, qseq_rr_comp);
			unCompDNA(qseq_fr_comp, qseq_fr);
			unCompDNA(qseq_rr_comp, qseq_rr);
		}
		
		if(kmersize <= qseq->len) {
			if(read_score && kmersize <= qseq_r->len) { // PE
				unmapped = alnFragsPE(templates_index, matched_templates, template_lengths, mq, scoreT, mrc, minFrac, minlen, qseq_comp, qseq_r_comp, qseq_fr_comp, qseq_rr_comp, qseq->seq, qseq_r->seq, qseq_fr, qseq_rr, header, header_r, kmersize, bestTemplates, bestTemplates_r, alignment_scores, uniq_alignment_scores, best_start_pos, best_end_pos, &flag, &flag_r, &best_read_score, &read_score, seq_in, seq_indexes, frag_out_raw, points, NWmatrices, excludeOut, excludeDB);
			} else { // SE
				unmapped = alnFragsSE(templates_index, matched_templates, template_lengths, mq, scoreT, mrc, minFrac, minlen, rc_flag, qseq_comp, qseq_r_comp, qseq->seq, qseq_r->seq, qseq->len, kmersize, header, bestTemplates, alignment_scores, uniq_alignment_scores, bestTemplates_r, Lengths, best_start_pos, best_end_pos, &flag, &best_read_score, seq_in, seq_indexes, frag_out_raw, points, NWmatrices, excludeOut, excludeDB);
				//unmapped = alnFragsSE(templates_index, matched_templates, template_lengths, mq, scoreT, mrc, minlen, rc_flag, qseq_comp, qseq_r_comp, qseq->seq, qseq_r->seq, qseq->len, kmersize, header, bestTemplates, alignment_scores, uniq_alignment_scores, best_start_pos, best_end_pos, &flag, &best_read_score, seq_in, seq_indexes, frag_out_raw, points, NWmatrices, excludeOut, excludeDB);
			}
		} else {
			unmapped = 0;
		}
		
		if(sam && unmapped) {
			if(unmapped & 1) {
				stats[1] = flag;
				nibble2base(qseq->seq, qseq->len);
				samwrite(qseq, header, 0, 0, 0, stats); 
			}
			if(unmapped & 2) {
				stats[1] = flag_r;
				nibble2base(qseq_r->seq, qseq_r->len);
				samwrite(qseq_r, header_r, 0, 0, 0, stats); 
			}
		}
		
		/* dump seq to all */
		if(frag_out_all) {
			if((unmapped & 1) == 0) {
				updateAllFrag(qseq->seq, qseq->len, *matched_templates, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_all);
			}
			if((unmapped & 2) == 0) {
				updateAllFrag(qseq_r->seq, qseq_r->len, *matched_templates, read_score, best_start_pos, best_end_pos, bestTemplates, header_r, frag_out_all);
			}
		}
		lock(excludeIn);
	}
	unlock(excludeIn);
	
	destroyComp(qseq_fr_comp);
	destroyComp(qseq_rr_comp);
	free(qseq_fr);
	free(qseq_rr);
	
	return NULL;
}
