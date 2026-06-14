/* Philip T.L.C. Clausen Oct 2021 plan@dtu.dk */

/*
 * Copyright (c) 2021, Philip Clausen, Technical University of Denmark
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

#include <stdio.h>
#include <stdlib.h>
#include "conclave.h"
#include "frags.h"
#include "limits.h"
#include "pherror.h"
#include "qseqs.h"
#include "stdnuc.h"
#include "stdstat.h"

int (*ConClavePtr)(FILE *, FILE ***, int, int, long unsigned *, unsigned *, unsigned *, long unsigned *, long unsigned *, int *, Qseqs *, Qseqs *, int *, int *, int *, Frag **) = &runConClave;
int (*ConClave2Ptr)(FILE *, FILE ***, int, int, long unsigned *, unsigned *, unsigned *, long unsigned *, long unsigned *, int *, Qseqs *, Qseqs *, int *, int *, int *, Frag **, long unsigned, double, double) = &runConClave2;

unsigned char * ustrdup(unsigned char *src, size_t n) {
	
	unsigned char *dest;
	
	dest = smalloc(n);
	memcpy(dest, src, n);
	
	return dest;
}

int runConClave(FILE *frag_in_raw, FILE ***Template_fragments, int DB_size, int maxFrag, long unsigned *w_scores, unsigned *fragmentCounts, unsigned *readCounts, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *template_lengths, Qseqs *header, Qseqs *qseq, int *bestTemplates, int *best_start_pos, int *best_end_pos, Frag **alignFrags) {
	
	int i, fileCount, sparse, bestHits, read_score, flag, start, end, bestNum;
	int bestTemplate, best_read_score, fragCount, tmp_start, tmp_end;
	int tmp_template, tmp_tmp_template, stats[5], *qBoundPtr;
	double bestScore, tmp_score;
	FILE **template_fragments;
	Frag *alignFrag;
	
	/* init */
	fileCount = 0;
	fragCount = 0;
	template_fragments = *Template_fragments;
	
	while(fread(stats, sizeof(int), 5, frag_in_raw) && stats[0] != 0) {
		qseq->len = stats[0];
		sparse = stats[1];
		bestHits = abs(sparse);
		read_score = abs(stats[2]);
		header->len = stats[3];
		flag = stats[4];
		
		sfread(qseq->seq, 1, qseq->len, frag_in_raw);
		sfread(header->seq, 1, header->len, frag_in_raw);
		sfread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
		sfread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
		sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
		
		/* Several mapped templates, choose best */
		if(bestHits > 1) {
			bestTemplate = -1;
			bestScore = 0;
			best_read_score = 0;
			bestNum = 0;
			start = 0;
			end = 0;
			/* iterate hits */
			for(i = 0; i != bestHits; ++i) {
				tmp_tmp_template = bestTemplates[i];
				tmp_start = best_start_pos[i];
				tmp_end = best_end_pos[i];
				if(tmp_tmp_template < 0) {
					tmp_template = -tmp_tmp_template;
				} else {
					tmp_template = tmp_tmp_template;
				}
				tmp_score = 1.0 * alignment_scores[tmp_template] / template_lengths[tmp_template];
				if(alignment_scores[tmp_template] > best_read_score) {
					bestTemplate = tmp_tmp_template;
					best_read_score = alignment_scores[tmp_template];
					bestScore = tmp_score;
					bestNum = uniq_alignment_scores[tmp_template];
					start = tmp_start;
					end = tmp_end;
				} else if(alignment_scores[tmp_template] == best_read_score) {
					if(tmp_score > bestScore) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template];
						bestScore = tmp_score;
						bestNum = uniq_alignment_scores[tmp_template];
						start = tmp_start;
						end = tmp_end;
					} else if(tmp_score == bestScore) {
						if(uniq_alignment_scores[tmp_template] > bestNum) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
							start = tmp_start;
							end = tmp_end;
						} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
							start = tmp_start;
							end = tmp_end;
						}
					}
				}
			}
		} else {
			bestTemplate = *bestTemplates;
			start = *best_start_pos;
			end = *best_end_pos;
		}
		
		/* reverse complement seq */
		if(bestTemplate < 0) {
			bestTemplate = -bestTemplate;
			strrc(qseq->seq, qseq->len);
			flag |= 16;
			
			/* check and invert query bounds */
			if(2 * sizeof(int) + 1 < header->len && header->seq[header->len - 2 * sizeof(int) - 1] == 0) {
				/* get q-bounds */
				qBoundPtr = (int*) (header->seq + (header->len - 2 * sizeof(int)));
				tmp_start = *qBoundPtr;
				tmp_end = *++qBoundPtr;
				/* invert */
				*qBoundPtr = qseq->len - tmp_start;
				*--qBoundPtr = qseq->len - tmp_end;
			}
		}
		w_scores[bestTemplate] += read_score;
		if(fragmentCounts) {
			fragmentCounts[bestTemplate]++;
			readCounts[bestTemplate]++;
		}
		
		/* dump frag info */
		alignFrag = smalloc(sizeof(Frag));
		alignFrag->buffer[0] = qseq->len;
		alignFrag->buffer[1] = bestHits;
		alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
		alignFrag->buffer[3] = start;
		alignFrag->buffer[4] = end;
		alignFrag->buffer[5] = header->len;
		alignFrag->buffer[6] = flag;
		alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
		alignFrag->header = ustrdup(header->seq, header->len);
		alignFrag->next = alignFrags[bestTemplate];
		alignFrags[bestTemplate] = alignFrag;
		
		++fragCount;
		
		if(stats[2] < 0) {
			if(readCounts) {
				readCounts[bestTemplate]++;
			}
			sfread(stats, sizeof(int), 3, frag_in_raw);
			qseq->len = stats[0];
			header->len = stats[1];
			flag = stats[2];
			sfread(qseq->seq, 1, qseq->len, frag_in_raw);
			sfread(header->seq, 1, header->len, frag_in_raw);
			/* dump frag info */
			alignFrag = smalloc(sizeof(Frag));
			alignFrag->buffer[0] = qseq->len;
			alignFrag->buffer[1] = bestHits;
			alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
			alignFrag->buffer[3] = start;
			alignFrag->buffer[4] = end;
			alignFrag->buffer[5] = header->len;
			alignFrag->buffer[6] = flag;
			alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
			alignFrag->header = ustrdup(header->seq, header->len);
			alignFrag->next = alignFrags[bestTemplate];
			alignFrags[bestTemplate] = alignFrag;
			
			++fragCount;
		}
		
		if(fragCount >= maxFrag) {
			template_fragments[fileCount] = printFrags(alignFrags, DB_size);
			++fileCount;
			fragCount = 0;
			/* control fileamount */
			if(fileCount >= DB_size) {
				template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
				if(!template_fragments) {
					ERROR();
				}
			}
		}
	}
	template_fragments[fileCount] = printFrags(alignFrags, DB_size);
	*Template_fragments = template_fragments;
	
	return ++fileCount;
}

int runConClave_lc(FILE *frag_in_raw, FILE ***Template_fragments, int DB_size, int maxFrag, long unsigned *w_scores, unsigned *fragmentCounts, unsigned *readCounts, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *template_lengths, Qseqs *header, Qseqs *qseq, int *bestTemplates, int *best_start_pos, int *best_end_pos, Frag **alignFrags) {
	
	int i, fileCount, sparse, bestHits, read_score, flag, start, end, bestNum;
	int bestTemplate, best_read_score, fragCount, tmp_start, tmp_end;
	int tmp_template, tmp_tmp_template, stats[5], *qBoundPtr;
	double bestScore, tmp_score;
	FILE **template_fragments;
	Frag *alignFrag;
	
	/* init */
	fileCount = 0;
	fragCount = 0;
	template_fragments = *Template_fragments;
	
	while(fread(stats, sizeof(int), 5, frag_in_raw) && stats[0] != 0) {
		qseq->len = stats[0];
		sparse = stats[1];
		bestHits = abs(sparse);
		read_score = abs(stats[2]);
		header->len = stats[3];
		flag = stats[4];
		
		sfread(qseq->seq, 1, qseq->len, frag_in_raw);
		sfread(header->seq, 1, header->len, frag_in_raw);
		sfread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
		sfread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
		sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
		/* Several mapped templates, choose best */
		if(bestHits > 1) {
			bestTemplate = -1;
			bestScore = 0;
			best_read_score = 0;
			bestNum = 0;
			start = 0;
			end = 0;
			/* iterate hits */
			for(i = 0; i != bestHits; ++i) {
				tmp_tmp_template = bestTemplates[i];
				tmp_start = best_start_pos[i];
				tmp_end = best_end_pos[i];
				if(tmp_tmp_template < 0) {
					tmp_template = -tmp_tmp_template;
				} else {
					tmp_template = tmp_tmp_template;
				}
				tmp_score = 1.0 * alignment_scores[tmp_template] / template_lengths[tmp_template];
				if(tmp_score > bestScore) {
					bestTemplate = tmp_tmp_template;
					best_read_score = alignment_scores[tmp_template];
					bestScore = tmp_score;
					bestNum = uniq_alignment_scores[tmp_template];
					start = tmp_start;
					end = tmp_end;
				} else if(tmp_score == bestScore) {
					if(alignment_scores[tmp_template] > best_read_score) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template];
						bestScore = tmp_score;
						bestNum = uniq_alignment_scores[tmp_template];
						start = tmp_start;
						end = tmp_end;
					} else if(alignment_scores[tmp_template] == best_read_score) {
						if(uniq_alignment_scores[tmp_template] > bestNum) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
							start = tmp_start;
							end = tmp_end;
						} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
							start = tmp_start;
							end = tmp_end;
						}
					}
				}
			}
		} else {
			bestTemplate = *bestTemplates;
			start = *best_start_pos;
			end = *best_end_pos;
		}
		
		/* reverse complement seq */
		if(bestTemplate < 0) {
			bestTemplate = -bestTemplate;
			strrc(qseq->seq, qseq->len);
			flag |= 16;
			
			/* check and invert query bounds */
			if(2 * sizeof(int) + 1 < header->len && header->seq[header->len - 2 * sizeof(int) - 1] == 0) {
				/* get q-bounds */
				qBoundPtr = (int*) (header->seq + (header->len - 2 * sizeof(int)));
				tmp_start = *qBoundPtr;
				tmp_end = *++qBoundPtr;
				/* invert */
				*qBoundPtr = qseq->len - tmp_start;
				*--qBoundPtr = qseq->len - tmp_end;
			}
		}
		w_scores[bestTemplate] += read_score;
		if(fragmentCounts) {
			fragmentCounts[bestTemplate]++;
			readCounts[bestTemplate]++;
		}
		
		/* dump frag info */
		alignFrag = smalloc(sizeof(Frag));
		alignFrag->buffer[0] = qseq->len;
		alignFrag->buffer[1] = bestHits;
		alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
		alignFrag->buffer[3] = start;
		alignFrag->buffer[4] = end;
		alignFrag->buffer[5] = header->len;
		alignFrag->buffer[6] = flag;
		alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
		alignFrag->header = ustrdup(header->seq, header->len);
		alignFrag->next = alignFrags[bestTemplate];
		alignFrags[bestTemplate] = alignFrag;
		
		++fragCount;
		
		if(stats[2] < 0) {
			if(readCounts) {
				readCounts[bestTemplate]++;
			}
			sfread(stats, sizeof(int), 3, frag_in_raw);
			qseq->len = stats[0];
			header->len = stats[1];
			flag = stats[2];
			sfread(qseq->seq, 1, qseq->len, frag_in_raw);
			sfread(header->seq, 1, header->len, frag_in_raw);
			/* dump frag info */
			alignFrag = smalloc(sizeof(Frag));
			alignFrag->buffer[0] = qseq->len;
			alignFrag->buffer[1] = bestHits;
			alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
			alignFrag->buffer[3] = start;
			alignFrag->buffer[4] = end;
			alignFrag->buffer[5] = header->len;
			alignFrag->buffer[6] = flag;
			alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
			alignFrag->header = ustrdup(header->seq, header->len);
			alignFrag->next = alignFrags[bestTemplate];
			alignFrags[bestTemplate] = alignFrag;
			
			++fragCount;
		}
		
		if(fragCount >= maxFrag) {
			template_fragments[fileCount] = printFrags(alignFrags, DB_size);
			++fileCount;
			fragCount = 0;
			/* control fileamount */
			if(fileCount >= DB_size) {
				template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
				if(!template_fragments) {
					ERROR();
				}
			}
		}
	}
	template_fragments[fileCount] = printFrags(alignFrags, DB_size);
	*Template_fragments = template_fragments;
	
	return ++fileCount;
}

int runConClave2(FILE *frag_in_raw, FILE ***Template_fragments, int DB_size, int maxFrag, long unsigned *w_scores, unsigned *fragmentCounts, unsigned *readCounts, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *template_lengths, Qseqs *header, Qseqs *qseq, int *bestTemplates, int *best_start_pos, int *best_end_pos, Frag **alignFrags, long unsigned template_tot_ulen, double scoreT, double evalue) {
	
	int i, j, fileCount, sparse, bestHits, read_score, flag, start, end, rand;
	int template, bestTemplate, best_read_score, fragCount, t_len, tot;
	int tmp_start, tmp_end, tmp_template, tmp_tmp_template;
	int stats[5], *qBoundPtr;
	unsigned randScore;
	long unsigned Nhits, bestNum, score;
	double bestScore, tmp_score, p_value;
	long double expected, q_value;
	FILE **template_fragments;
	Frag *alignFrag;
	
	/* init */
	fileCount = 0;
	fragCount = 0;
	template_fragments = *Template_fragments;
	
	/* find potential template candidates */
	while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
		qseq->len = stats[0];
		sparse = stats[1];
		bestHits = abs(sparse);
		read_score = abs(stats[2]);
		header->len = stats[3];
		
		/* best templates, skip rest */
		sfseek(frag_in_raw, qseq->len + header->len + (2 * bestHits + 1) * sizeof(int), SEEK_CUR);
		sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
		
		/* Several mapped templates, choose best */
		if(bestHits > 1) {
			bestTemplate = -1;
			bestScore = 0;
			best_read_score = 0;
			bestNum = 0;
			/* iterate hits */
			for(i = 0; i != bestHits; ++i) {
				tmp_tmp_template = bestTemplates[i];
				if(tmp_tmp_template < 0) {
					tmp_template = -tmp_tmp_template;
				} else {
					tmp_template = tmp_tmp_template;
				}
				tmp_score = 1.0 * alignment_scores[tmp_template] / template_lengths[tmp_template];
				if(alignment_scores[tmp_template] > best_read_score) {
					bestTemplate = tmp_tmp_template;
					best_read_score = alignment_scores[tmp_template];
					bestScore = tmp_score;
					bestNum = uniq_alignment_scores[tmp_template];
				} else if(alignment_scores[tmp_template] == best_read_score) {
					if(tmp_score > bestScore) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template];
						bestScore = tmp_score;
						bestNum = uniq_alignment_scores[tmp_template];
					} else if(tmp_score == bestScore) {
						if(uniq_alignment_scores[tmp_template] > bestNum) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
						} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
						}
					}
				}
			}
		} else {
			bestTemplate = *bestTemplates;
		}
		w_scores[abs(bestTemplate)] += read_score;
		
		if(stats[2] < 0) {
			sfread(stats, sizeof(int), 2, frag_in_raw);
			sfseek(frag_in_raw, stats[0] + stats[1] + sizeof(int), SEEK_CUR);
		}
	}
	rewind(frag_in_raw);
	
	/* discard insignifiacant templates */
	Nhits = 0;
	template = DB_size;
	while(--template) {
		Nhits += w_scores[template];
	}
	
	template = DB_size;
	while(--template) {
		if((read_score = w_scores[template])) {
			t_len = template_lengths[template];
			expected = t_len;
			expected /= MAX(1, (template_tot_ulen - t_len));
			expected *= (Nhits - read_score);
			q_value = read_score - expected;
			q_value /= (expected + read_score);
			q_value *= read_score - expected;
			p_value  = p_chisqr(q_value);
			if(cmp((p_value <= evalue && read_score > expected), (read_score >= scoreT * t_len)) == 0) {
				w_scores[template] = 0;
			}
		}
	}
	
	/* identify sorting keys */
	while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
		qseq->len = stats[0];
		sparse = stats[1];
		bestHits = abs(sparse);
		read_score = abs(stats[2]);
		header->len = stats[3];
		
		if(bestHits != 1) {
			/* best templates, skip rest */
			sfseek(frag_in_raw, qseq->len + header->len + (2 * bestHits + 1) * sizeof(int), SEEK_CUR);
			sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			bestTemplate = 0;
			i = bestHits;
			while(i--) {
				template = abs(bestTemplates[i]);
				if(w_scores[template]) {
					if(bestTemplate) {
						bestTemplate = 0;
						break;
					} else {
						bestTemplate = template;
					}
				}
			}
			
			if(bestTemplate) {
				uniq_alignment_scores[bestTemplate] += read_score;
			}
		} else {
			/* skip rest */
			sfseek(frag_in_raw, qseq->len + header->len + 4 * sizeof(int), SEEK_CUR);
		}
		
		if(stats[2] < 0) {
			sfread(stats, sizeof(int), 2, frag_in_raw);
			sfseek(frag_in_raw, stats[0] + stats[1] + sizeof(int), SEEK_CUR);
		}
	}
	rewind(frag_in_raw);
	
	/* choose the templates */
	memset(w_scores, 0, DB_size * sizeof(long unsigned));
	while(fread(stats, sizeof(int), 5, frag_in_raw) && stats[0] != 0) {
		qseq->len = stats[0];
		sparse = stats[1];
		bestHits = abs(sparse);
		read_score = abs(stats[2]);
		header->len = stats[3];
		flag = stats[4];
		
		sfread(qseq->seq, 1, qseq->len, frag_in_raw);
		sfread(header->seq, 1, header->len, frag_in_raw);
		sfread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
		sfread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
		sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
		
		/* Several mapped templates, choose best according to sorting keys */
		if(bestHits != 1) {
			
			bestTemplate = 0;
			bestScore = 0;
			start = 0;
			end = 0;
			
			tot = 0;
			i = bestHits;
			while(i--) {
				tot += uniq_alignment_scores[abs(bestTemplates[i])];
			}
			
			if(tot && 16 <= qseq->len) {
				/* get seed */
				rand = qseq->seq[0];
				i = -1;
				j = qseq->len;
				while(++i < 7) {
					rand = (((rand << 2) | qseq->seq[i]) << 2) | qseq->seq[--j];
				}
				/* minimal standard */
				rand = 16807 * (rand % 127773) - 2836 * (rand / 127773);
				if (rand <= 0) {
					rand += 0x7fffffff;
				}
				
				tmp_score = rand;
				tmp_score /= INT_MAX;
				randScore = tmp_score * tot;
				
				score = 0;
				i = 0;
				while(i != bestHits) {
					score += uniq_alignment_scores[abs(bestTemplates[i])];
					if(randScore < score) {
						bestTemplate = bestTemplates[i];
						start = best_start_pos[i];
						end = best_end_pos[i];
						i = bestHits;
					} else {
						++i;
					}
				}
				
				if(bestTemplate == 0) {
					tot = 0;
				}
			} else {
				tot = 0;
			}
			
			if(tot == 0) {
				bestTemplate = 0;
				best_read_score = 0;
				bestNum = 0;
				
				/* iterate hits */
				for(i = 0; i != bestHits; ++i) {
					tmp_tmp_template = bestTemplates[i];
					tmp_start = best_start_pos[i];
					tmp_end = best_end_pos[i];
					if(tmp_tmp_template < 0) {
						tmp_template = -tmp_tmp_template;
					} else {
						tmp_template = tmp_tmp_template;
					}
					tmp_score = 1.0 * alignment_scores[tmp_template] / template_lengths[tmp_template];
					if(alignment_scores[tmp_template] > best_read_score) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template];
						bestScore = tmp_score;
						bestNum = uniq_alignment_scores[tmp_template];
						start = tmp_start;
						end = tmp_end;
					} else if(alignment_scores[tmp_template] == best_read_score) {
						if(tmp_score > bestScore) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
							start = tmp_start;
							end = tmp_end;
						} else if(tmp_score == bestScore) {
							if(uniq_alignment_scores[tmp_template] > bestNum) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							}
						}
					}
				}
			}
		} else {
			bestTemplate = *bestTemplates;
			start = *best_start_pos;
			end = *best_end_pos;
		}
		
		/* reverse complement seq */
		if(bestTemplate < 0) {
			bestTemplate = -bestTemplate;
			strrc(qseq->seq, qseq->len);
			flag |= 16;
			
			/* check and invert query bounds */
			if(2 * sizeof(int) + 1 < header->len && header->seq[header->len - 2 * sizeof(int) - 1] == 0) {
				/* get q-bounds */
				qBoundPtr = (int*) (header->seq + (header->len - 2 * sizeof(int)));
				tmp_start = *qBoundPtr;
				tmp_end = *++qBoundPtr;
				/* invert */
				*qBoundPtr = qseq->len - tmp_start;
				*--qBoundPtr = qseq->len - tmp_end;
			}
		}
		if(bestTemplate) {
			w_scores[bestTemplate] += read_score;
			if(fragmentCounts) {
				fragmentCounts[bestTemplate]++;
				readCounts[bestTemplate]++;
			}
			
			/* dump frag info */
			alignFrag = smalloc(sizeof(Frag));
			alignFrag->buffer[0] = qseq->len;
			alignFrag->buffer[1] = bestHits;
			alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
			alignFrag->buffer[3] = start;
			alignFrag->buffer[4] = end;
			alignFrag->buffer[5] = header->len;
			alignFrag->buffer[6] = flag;
			alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
			alignFrag->header = ustrdup(header->seq, header->len);
			alignFrag->next = alignFrags[bestTemplate];
			alignFrags[bestTemplate] = alignFrag;
			
			++fragCount;
			
			if(stats[2] < 0) {
				if(readCounts) {
					readCounts[bestTemplate]++;
				}
				sfread(stats, sizeof(int), 3, frag_in_raw);
				qseq->len = stats[0];
				header->len = stats[1];
				flag = stats[2];
				sfread(qseq->seq, 1, qseq->len, frag_in_raw);
				sfread(header->seq, 1, header->len, frag_in_raw);
				/* dump frag info */
				alignFrag = smalloc(sizeof(Frag));
				alignFrag->buffer[0] = qseq->len;
				alignFrag->buffer[1] = bestHits;
				alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
				alignFrag->buffer[3] = start;
				alignFrag->buffer[4] = end;
				alignFrag->buffer[5] = header->len;
				alignFrag->buffer[6] = flag;
				alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
				alignFrag->header = ustrdup(header->seq, header->len);
				alignFrag->next = alignFrags[bestTemplate];
				alignFrags[bestTemplate] = alignFrag;
				
				++fragCount;
			}
		} else if(stats[2] < 0) {
			sfread(stats, sizeof(int), 2, frag_in_raw);
			sfseek(frag_in_raw, stats[0] + stats[1] + sizeof(int), SEEK_CUR);
		}
		
		if(fragCount >= maxFrag) {
			template_fragments[fileCount] = printFrags(alignFrags, DB_size);
			++fileCount;
			fragCount = 0;
			/* control fileamount */
			if(fileCount >= DB_size) {
				template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
				if(!template_fragments) {
					ERROR();
				}
			}
		}
	}
	template_fragments[fileCount] = printFrags(alignFrags, DB_size);
	
	return ++fileCount;
}

int runConClave2_lc(FILE *frag_in_raw, FILE ***Template_fragments, int DB_size, int maxFrag, long unsigned *w_scores, unsigned *fragmentCounts, unsigned *readCounts, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *template_lengths, Qseqs *header, Qseqs *qseq, int *bestTemplates, int *best_start_pos, int *best_end_pos, Frag **alignFrags, long unsigned template_tot_ulen, double scoreT, double evalue) {
	
	int i, j, fileCount, sparse, bestHits, read_score, flag, start, end, rand;
	int template, bestTemplate, best_read_score, fragCount, t_len, tot;
	int tmp_start, tmp_end, tmp_template, tmp_tmp_template;
	int stats[5], *qBoundPtr;
	unsigned randScore;
	long unsigned Nhits, bestNum, score;
	double bestScore, tmp_score, p_value;
	long double expected, q_value;
	FILE **template_fragments;
	Frag *alignFrag;
	
	/* init */
	fileCount = 0;
	fragCount = 0;
	template_fragments = *Template_fragments;
	
	/* find potential template candidates */
	while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
		qseq->len = stats[0];
		sparse = stats[1];
		bestHits = abs(sparse);
		read_score = abs(stats[2]);
		header->len = stats[3];
		
		/* best templates, skip rest */
		sfseek(frag_in_raw, qseq->len + header->len + (2 * bestHits + 1) * sizeof(int), SEEK_CUR);
		sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
		
		/* Several mapped templates, choose best */
		if(bestHits > 1) {
			bestTemplate = -1;
			bestScore = 0;
			best_read_score = 0;
			bestNum = 0;
			/* iterate hits */
			for(i = 0; i != bestHits; ++i) {
				tmp_tmp_template = bestTemplates[i];
				if(tmp_tmp_template < 0) {
					tmp_template = -tmp_tmp_template;
				} else {
					tmp_template = tmp_tmp_template;
				}
				tmp_score = 1.0 * alignment_scores[tmp_template] / template_lengths[tmp_template];
				if(tmp_score > bestScore) {
					bestTemplate = tmp_tmp_template;
					best_read_score = alignment_scores[tmp_template];
					bestScore = tmp_score;
					bestNum = uniq_alignment_scores[tmp_template];
				} else if(tmp_score == bestScore) {
					if(alignment_scores[tmp_template] > best_read_score) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template];
						bestScore = tmp_score;
						bestNum = uniq_alignment_scores[tmp_template];
					} else if(alignment_scores[tmp_template] == best_read_score) {
						if(uniq_alignment_scores[tmp_template] > bestNum) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
						} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
						}
					}
				}
			}
		} else {
			bestTemplate = *bestTemplates;
		}
		w_scores[abs(bestTemplate)] += read_score;
		
		if(stats[2] < 0) {
			sfread(stats, sizeof(int), 2, frag_in_raw);
			sfseek(frag_in_raw, stats[0] + stats[1] + sizeof(int), SEEK_CUR);
		}
	}
	rewind(frag_in_raw);
	
	/* discard insignifiacant templates */
	Nhits = 0;
	template = DB_size;
	while(--template) {
		Nhits += w_scores[template];
	}
	
	template = DB_size;
	while(--template) {
		if((read_score = w_scores[template])) {
			t_len = template_lengths[template];
			expected = t_len;
			expected /= MAX(1, (template_tot_ulen - t_len));
			expected *= (Nhits - read_score);
			q_value = read_score - expected;
			q_value /= (expected + read_score);
			q_value *= read_score - expected;
			p_value  = p_chisqr(q_value);
			if(cmp((p_value <= evalue && read_score > expected), (read_score >= scoreT * t_len)) == 0) {
				w_scores[template] = 0;
			}
		}
	}
	
	/* identify sorting keys */
	while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
		qseq->len = stats[0];
		sparse = stats[1];
		bestHits = abs(sparse);
		read_score = abs(stats[2]);
		header->len = stats[3];
		
		if(bestHits != 1) {
			/* best templates, skip rest */
			sfseek(frag_in_raw, qseq->len + header->len + (2 * bestHits + 1) * sizeof(int), SEEK_CUR);
			sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			bestTemplate = 0;
			i = bestHits;
			while(i--) {
				template = abs(bestTemplates[i]);
				if(w_scores[template]) {
					if(bestTemplate) {
						bestTemplate = 0;
						break;
					} else {
						bestTemplate = template;
					}
				}
			}
			
			if(bestTemplate) {
				uniq_alignment_scores[bestTemplate] += read_score;
			}
		} else {
			/* skip rest */
			sfseek(frag_in_raw, qseq->len + header->len + 4 * sizeof(int), SEEK_CUR);
		}
		
		if(stats[2] < 0) {
			sfread(stats, sizeof(int), 2, frag_in_raw);
			sfseek(frag_in_raw, stats[0] + stats[1] + sizeof(int), SEEK_CUR);
		}
	}
	rewind(frag_in_raw);
	
	/* choose the templates */
	memset(w_scores, 0, DB_size * sizeof(long unsigned));
	while(fread(stats, sizeof(int), 5, frag_in_raw) && stats[0] != 0) {
		qseq->len = stats[0];
		sparse = stats[1];
		bestHits = abs(sparse);
		read_score = abs(stats[2]);
		header->len = stats[3];
		flag = stats[4];
		
		sfread(qseq->seq, 1, qseq->len, frag_in_raw);
		sfread(header->seq, 1, header->len, frag_in_raw);
		sfread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
		sfread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
		sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
		
		/* Several mapped templates, choose best according to sorting keys */
		if(bestHits != 1) {
			
			bestTemplate = 0;
			bestScore = 0;
			start = 0;
			end = 0;
			
			tot = 0;
			i = bestHits;
			while(i--) {
				tot += uniq_alignment_scores[abs(bestTemplates[i])];
			}
			
			if(tot && 16 <= qseq->len) {
				/* get seed */
				rand = qseq->seq[0];
				i = -1;
				j = qseq->len;
				while(++i < 7) {
					rand = (((rand << 2) | qseq->seq[i]) << 2) | qseq->seq[--j];
				}
				/* minimal standard */
				rand = 16807 * (rand % 127773) - 2836 * (rand / 127773);
				if (rand <= 0) {
					rand += 0x7fffffff;
				}
				
				tmp_score = rand;
				tmp_score /= INT_MAX;
				randScore = tmp_score * tot;
				
				score = 0;
				i = 0;
				while(i != bestHits) {
					score += uniq_alignment_scores[abs(bestTemplates[i])];
					if(randScore < score) {
						bestTemplate = bestTemplates[i];
						start = best_start_pos[i];
						end = best_end_pos[i];
						i = bestHits;
					} else {
						++i;
					}
				}
				
				if(bestTemplate == 0) {
					tot = 0;
				}
			} else {
				tot = 0;
			}
			
			if(tot == 0) {
				bestTemplate = 0;
				best_read_score = 0;
				bestNum = 0;
				
				/* iterate hits */
				for(i = 0; i != bestHits; ++i) {
					tmp_tmp_template = bestTemplates[i];
					tmp_start = best_start_pos[i];
					tmp_end = best_end_pos[i];
					if(tmp_tmp_template < 0) {
						tmp_template = -tmp_tmp_template;
					} else {
						tmp_template = tmp_tmp_template;
					}
					tmp_score = 1.0 * alignment_scores[tmp_template] / template_lengths[tmp_template];
					if(tmp_score > bestScore) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template];
						bestScore = tmp_score;
						bestNum = uniq_alignment_scores[tmp_template];
						start = tmp_start;
						end = tmp_end;
					} else if(tmp_score == bestScore) {
						if(alignment_scores[tmp_template] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
							start = tmp_start;
							end = tmp_end;
						} else if(alignment_scores[tmp_template] == best_read_score) {
							if(uniq_alignment_scores[tmp_template] > bestNum) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							}
						}
					}
				}
			}
		} else {
			bestTemplate = *bestTemplates;
			start = *best_start_pos;
			end = *best_end_pos;
		}
		
		/* reverse complement seq */
		if(bestTemplate < 0) {
			bestTemplate = -bestTemplate;
			strrc(qseq->seq, qseq->len);
			flag |= 16;
			
			/* check and invert query bounds */
			if(2 * sizeof(int) + 1 < header->len && header->seq[header->len - 2 * sizeof(int) - 1] == 0) {
				/* get q-bounds */
				qBoundPtr = (int*) (header->seq + (header->len - 2 * sizeof(int)));
				tmp_start = *qBoundPtr;
				tmp_end = *++qBoundPtr;
				/* invert */
				*qBoundPtr = qseq->len - tmp_start;
				*--qBoundPtr = qseq->len - tmp_end;
			}
		}
		if(bestTemplate) {
			w_scores[bestTemplate] += read_score;
			if(fragmentCounts) {
				fragmentCounts[bestTemplate]++;
				readCounts[bestTemplate]++;
			}
			
			/* dump frag info */
			alignFrag = smalloc(sizeof(Frag));
			alignFrag->buffer[0] = qseq->len;
			alignFrag->buffer[1] = bestHits;
			alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
			alignFrag->buffer[3] = start;
			alignFrag->buffer[4] = end;
			alignFrag->buffer[5] = header->len;
			alignFrag->buffer[6] = flag;
			alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
			alignFrag->header = ustrdup(header->seq, header->len);
			alignFrag->next = alignFrags[bestTemplate];
			alignFrags[bestTemplate] = alignFrag;
			
			++fragCount;
			
			if(stats[2] < 0) {
				if(readCounts) {
					readCounts[bestTemplate]++;
				}
				sfread(stats, sizeof(int), 3, frag_in_raw);
				qseq->len = stats[0];
				header->len = stats[1];
				flag = stats[2];
				sfread(qseq->seq, 1, qseq->len, frag_in_raw);
				sfread(header->seq, 1, header->len, frag_in_raw);
				/* dump frag info */
				alignFrag = smalloc(sizeof(Frag));
				alignFrag->buffer[0] = qseq->len;
				alignFrag->buffer[1] = bestHits;
				alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
				alignFrag->buffer[3] = start;
				alignFrag->buffer[4] = end;
				alignFrag->buffer[5] = header->len;
				alignFrag->buffer[6] = flag;
				alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
				alignFrag->header = ustrdup(header->seq, header->len);
				alignFrag->next = alignFrags[bestTemplate];
				alignFrags[bestTemplate] = alignFrag;
				
				++fragCount;
			}
		} else if(stats[2] < 0) {
			sfread(stats, sizeof(int), 2, frag_in_raw);
			sfseek(frag_in_raw, stats[0] + stats[1] + sizeof(int), SEEK_CUR);
		}
		
		if(fragCount >= maxFrag) {
			template_fragments[fileCount] = printFrags(alignFrags, DB_size);
			++fileCount;
			fragCount = 0;
			/* control fileamount */
			if(fileCount >= DB_size) {
				template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
				if(!template_fragments) {
					ERROR();
				}
			}
		}
	}
	template_fragments[fileCount] = printFrags(alignFrags, DB_size);
	
	return ++fileCount;
}
