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

#include <stdio.h>
#include <stdlib.h>
#include "pherror.h"
#include "qseqs.h"
#include "updatescores.h"

void update_Scores_MEM(unsigned char *qseq, int q_len, int counter, int score, int *start, int *end, int *template, Qseqs *header, int flag, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw) {
	
	int i, buffer[5];
	
	/* print frag */
	buffer[0] = q_len;
	buffer[1] = counter;
	buffer[2] = score;
	buffer[3] = header->len;
	buffer[4] = flag;
	counter = abs(counter);
	sfwrite(buffer, sizeof(int), 5, frag_out_raw);
	sfwrite(qseq, 1, q_len, frag_out_raw);
	sfwrite(header->seq, 1, header->len, frag_out_raw);
	sfwrite(start, sizeof(int), counter, frag_out_raw);
	sfwrite(end, sizeof(int), counter, frag_out_raw);
	sfwrite(template, sizeof(int), counter, frag_out_raw);
	
	/* update scores */
	if(counter == 1) { /* Only one best match */
		if(*template < 0) {
			*template = -*template;
		}
		alignment_scores[*template] += score;
		uniq_alignment_scores[*template] += score;
	} else {
		i = counter + 1;
		--template;
		while(--i) {
			if(*++template < 0) {
				alignment_scores[-*template] += score;
			} else {
				alignment_scores[*template] += score;
			}
		}
	}
}

void update_Scores_pe_MEM(unsigned char *qseq, int q_len, unsigned char *qseq_r, int qr_len, int counter, int score, int *start, int *end, int *template, Qseqs *header, Qseqs *header_r, int flag, int flag_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw) {
	
	int i, buffer[5];
	
	/* print frag */
	buffer[0] = q_len;
	buffer[1] = counter;
	buffer[2] = -score;
	buffer[3] = header->len;
	buffer[4] = flag;
	counter = abs(counter);
	sfwrite(buffer, sizeof(int), 5, frag_out_raw);
	sfwrite(qseq, 1, q_len, frag_out_raw);
	sfwrite(header->seq, 1, header->len, frag_out_raw);
	sfwrite(start, sizeof(int), counter, frag_out_raw);
	sfwrite(end, sizeof(int), counter, frag_out_raw);
	sfwrite(template, sizeof(int), counter, frag_out_raw);
	
	buffer[0] = qr_len;
	buffer[1] = header_r->len;
	buffer[2] = flag_r;
	sfwrite(buffer, sizeof(int), 3, frag_out_raw);
	sfwrite(qseq_r, 1, qr_len, frag_out_raw);
	sfwrite(header_r->seq, 1, header_r->len, frag_out_raw);
	
	/* update scores */
	if(counter == 1) { /* Only one best match */
		if(*template < 0) {
			*template = -*template;
		}
		alignment_scores[*template] += score;
		uniq_alignment_scores[*template] += score;
	} else {
		i = counter + 1;
		--template;
		while(--i) {
			if(*++template < 0) {
				alignment_scores[-*template] += score;
			} else {
				alignment_scores[*template] += score;
			}
		}
	}
}

void update_Scores_nanoold(unsigned char *qseq, int q_len, double minFrac, int counter, int bestScore, int bestLen, int *start, int *end, int *templates, int *Scores, int *Lengths, Qseqs *header, int flag, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw) {
	
	int i, score, buffer[5];
	int *bestTemplates, *bestStart, *bestEnd;
	double minScore;
	
	/* get best hits, and update ConClave scores */
	bestTemplates = --templates;
	bestStart = --start;
	bestEnd = --end;
	--Scores;
	--Lengths;
	i = counter + 1;
	counter = 0;
	if(minFrac == 1.0) {
		while(--i) {
			score = *++Scores;
			if(*++Lengths == bestLen && score == bestScore) {
				++counter;
				*++bestTemplates = *++templates;
				*++bestStart = *++start;
				*++bestEnd = *++end;
				
				/* update ConClave scores */
				alignment_scores[abs(*templates)] += score;
			} else {
				++templates;
				++start;
				++end;
			}
		}
	} else if(minFrac < 0) {
		minFrac = (-minFrac) * bestScore;
		minScore = minFrac / bestLen;
		while(--i) {
			score = *++Scores;
			if((*++Lengths * minScore <= score) && minFrac <= score) {
				++counter;
				*++bestTemplates = *++templates;
				*++bestStart = *++start;
				*++bestEnd = *++end;
				/* update ConClave scores */
				alignment_scores[abs(*templates)] += score;
			} else {
				++templates;
				++start;
				++end;
			}
		}
	} else {
		minFrac = minFrac * bestScore;
		minScore = minFrac / bestLen;
		while(--i) {
			score = *++Scores;
			if((*++Lengths * minScore <= score) && minFrac <= score) {
				++counter;
				*++bestTemplates = *++templates;
				*++bestStart = *++start;
				*++bestEnd = *++end;
				
				/* update ConClave scores */
				alignment_scores[abs(*templates)] += bestScore;
			} else {
				++templates;
				++start;
				++end;
			}
		}
	}
	
	/* uniq hit */
	if(counter == 1) {
		uniq_alignment_scores[abs(*bestTemplates)] += bestScore;
	} else {
		bestTemplates -= counter - 1;
		bestStart -= counter - 1;
		bestEnd -= counter - 1;
	}
	
	/* print frag */
	buffer[0] = q_len;
	buffer[1] = counter;
	buffer[2] = bestScore;
	buffer[3] = header->len;
	buffer[4] = flag;
	counter = abs(counter);
	sfwrite(buffer, sizeof(int), 5, frag_out_raw);
	sfwrite(qseq, 1, q_len, frag_out_raw);
	sfwrite(header->seq, 1, header->len, frag_out_raw);
	sfwrite(bestStart, sizeof(int), counter, frag_out_raw);
	sfwrite(bestEnd, sizeof(int), counter, frag_out_raw);
	sfwrite(bestTemplates, sizeof(int), counter, frag_out_raw);
}

int update_Scores(unsigned char *qseq, int q_len, double minFrac, int counter, int bestReadScore, double bestScore, int *start, int *end, int *templates, int *Scores, int *Lengths, Qseqs *header, int flag, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw) {
	
	int i, score, buffer[5];
	int *bestTemplates, *bestStart, *bestEnd;
	double minScore;
	
	/* get best hits, and update ConClave scores */
	bestTemplates = --templates;
	bestStart = --start;
	bestEnd = --end;
	--Scores;
	--Lengths;
	i = counter + 1;
	counter = 0;
	if(minFrac == 1.0) {
		while(--i) {
			score = *++Scores;
			minScore = score / *++Lengths;
			if(minScore == bestScore || score == bestReadScore) {
				++counter;
				*++bestTemplates = *++templates;
				*++bestStart = *++start;
				*++bestEnd = *++end;
				
				/* update ConClave scores */
				alignment_scores[abs(*templates)] += score;
			} else {
				++templates;
				++start;
				++end;
			}
		}
	} else if(minFrac < 0) {
		minScore = (-minFrac) * bestScore;
		minFrac = (-minFrac) * bestReadScore;
		while(--i) {
			score = *++Scores;
			if((*++Lengths * minScore <= score) || minFrac <= score) {
				++counter;
				*++bestTemplates = *++templates;
				*++bestStart = *++start;
				*++bestEnd = *++end;
				/* update ConClave scores */
				alignment_scores[abs(*templates)] += score;
			} else {
				++templates;
				++start;
				++end;
			}
		}
	} else {
		minScore = minFrac * bestScore;
		minFrac = minFrac * bestReadScore;
		while(--i) {
			score = *++Scores;
			if((*++Lengths * minScore <= score) || minFrac <= score) {
				++counter;
				*++bestTemplates = *++templates;
				*++bestStart = *++start;
				*++bestEnd = *++end;
				
				/* update ConClave scores */
				alignment_scores[abs(*templates)] += bestReadScore;
			} else {
				++templates;
				++start;
				++end;
			}
		}
	}
	
	/* uniq hit */
	if(counter == 1) {
		uniq_alignment_scores[abs(*bestTemplates)] += bestReadScore;
	} else {
		bestTemplates -= counter - 1;
		bestStart -= counter - 1;
		bestEnd -= counter - 1;
	}
	
	/* print frag */
	buffer[0] = q_len;
	buffer[1] = counter;
	buffer[2] = bestReadScore;
	buffer[3] = header->len;
	buffer[4] = flag;
	counter = abs(counter);
	sfwrite(buffer, sizeof(int), 5, frag_out_raw);
	sfwrite(qseq, 1, q_len, frag_out_raw);
	sfwrite(header->seq, 1, header->len, frag_out_raw);
	sfwrite(bestStart, sizeof(int), counter, frag_out_raw);
	sfwrite(bestEnd, sizeof(int), counter, frag_out_raw);
	sfwrite(bestTemplates, sizeof(int), counter, frag_out_raw);
	
	return counter;
}

int update_Scores_se(unsigned char *qseq, int q_len, double minFrac, int counter, int bestScore, int *start, int *end, int *templates, int *Scores, Qseqs *header, int flag, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw) {
	
	int i, buffer[5];
	int *bestTemplates, *bestStart, *bestEnd;
	
	/* get best hits, and update ConClave scores */
	bestTemplates = --templates;
	bestStart = --start;
	bestEnd = --end;
	--Scores;
	i = counter + 1;
	counter = 0;
	if(minFrac == 1.0) {
		while(--i) {
			if(*++Scores == bestScore) {
				++counter;
				*++bestTemplates = *++templates;
				*++bestStart = *++start;
				*++bestEnd = *++end;
				
				/* update ConClave scores */
				alignment_scores[abs(*templates)] += bestScore;
			} else {
				++templates;
				++start;
				++end;
			}
		}
	} else if(minFrac < 0) {
		minFrac = (-minFrac) * bestScore;
		while(--i) {
			if(minFrac <= *++Scores) {
				++counter;
				*++bestTemplates = *++templates;
				*++bestStart = *++start;
				*++bestEnd = *++end;
				
				/* update ConClave scores */
				alignment_scores[abs(*templates)] += *Scores;
			} else {
				++templates;
				++start;
				++end;
			}
		}
	} else {
		minFrac = minFrac * bestScore;
		while(--i) {
			if(minFrac <= *++Scores) {
				++counter;
				*++bestTemplates = *++templates;
				*++bestStart = *++start;
				*++bestEnd = *++end;
				
				/* update ConClave scores */
				alignment_scores[abs(*templates)] += bestScore;
			} else {
				++templates;
				++start;
				++end;
			}
		}
	}
	
	/* uniq hit */
	if(counter == 1) {
		uniq_alignment_scores[abs(*bestTemplates)] += bestScore;
	} else {
		bestTemplates -= counter - 1;
		bestStart -= counter - 1;
		bestEnd -= counter - 1;
	}
	
	/* print frag */
	buffer[0] = q_len;
	buffer[1] = counter;
	buffer[2] = bestScore;
	buffer[3] = header->len;
	buffer[4] = flag;
	counter = abs(counter);
	sfwrite(buffer, sizeof(int), 5, frag_out_raw);
	sfwrite(qseq, 1, q_len, frag_out_raw);
	sfwrite(header->seq, 1, header->len, frag_out_raw);
	sfwrite(bestStart, sizeof(int), counter, frag_out_raw);
	sfwrite(bestEnd, sizeof(int), counter, frag_out_raw);
	sfwrite(bestTemplates, sizeof(int), counter, frag_out_raw);
	
	return counter;
}

int update_Scores_pe(unsigned char *qseq, int q_len, unsigned char *qseq_r, int qr_len, double minFrac, int counter, int bestScore, int *start, int *end, int *templates, int *Scores, Qseqs *header, Qseqs *header_r, int flag, int flag_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw) {
	
	int i, score, buffer[5];
	int *bestTemplates, *bestStart, *bestEnd;
	
	/* get best hits, and update ConClave scores */
	bestTemplates = --templates;
	bestStart = --start;
	bestEnd = --end;
	--Scores;
	i = counter + 1;
	counter = 0;
	if(minFrac == 1.0) {
		while(--i) {
			score = *++Scores;
			if(score == bestScore) {
				++counter;
				*++bestTemplates = *++templates;
				*++bestStart = *++start;
				*++bestEnd = *++end;
				
				/* update ConClave scores */
				alignment_scores[abs(*templates)] += score;
			} else {
				++templates;
				++start;
				++end;
			}
		}
	} else if(minFrac < 0) {
		minFrac = (-minFrac) * bestScore;
		while(--i) {
			score = *++Scores;
			if(minFrac <= score) {
				++counter;
				*++bestTemplates = *++templates;
				*++bestStart = *++start;
				*++bestEnd = *++end;
				
				/* update ConClave scores */
				alignment_scores[abs(*templates)] += score;
			} else {
				++templates;
				++start;
				++end;
			}
		}
	} else {
		minFrac = minFrac * bestScore;
		while(--i) {
			score = *++Scores;
			if(minFrac <= score) {
				++counter;
				*++bestTemplates = *++templates;
				*++bestStart = *++start;
				*++bestEnd = *++end;
				
				/* update ConClave scores */
				alignment_scores[abs(*templates)] += bestScore;
			} else {
				++templates;
				++start;
				++end;
			}
		}
	}
	
	/* uniq hit */
	if(counter == 1) {
		uniq_alignment_scores[abs(*bestTemplates)] += bestScore;
	} else {
		bestTemplates -= counter - 1;
		bestStart -= counter - 1;
		bestEnd -= counter - 1;
	}
	
	/* print frag */
	buffer[0] = q_len;
	buffer[1] = counter;
	buffer[2] = -bestScore;
	buffer[3] = header->len;
	buffer[4] = flag;
	counter = abs(counter);
	sfwrite(buffer, sizeof(int), 5, frag_out_raw);
	sfwrite(qseq, 1, q_len, frag_out_raw);
	sfwrite(header->seq, 1, header->len, frag_out_raw);
	sfwrite(bestStart, sizeof(int), counter, frag_out_raw);
	sfwrite(bestEnd, sizeof(int), counter, frag_out_raw);
	sfwrite(bestTemplates, sizeof(int), counter, frag_out_raw);
	
	buffer[0] = qr_len;
	buffer[1] = header_r->len;
	buffer[2] = flag_r;
	sfwrite(buffer, sizeof(int), 3, frag_out_raw);
	sfwrite(qseq_r, 1, qr_len, frag_out_raw);
	sfwrite(header_r->seq, 1, header_r->len, frag_out_raw);
	
	return counter;
}

void update_Scores_old(unsigned char *qseq, int q_len, int counter, int score, int *start, int *end, int *template, Qseqs *header, int flag, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw) {
	
	int i, buffer[5];
	
	/* print frag */
	buffer[0] = q_len;
	buffer[1] = counter;
	buffer[2] = score;
	buffer[3] = header->len;
	buffer[4] = flag;
	counter = abs(counter);
	sfwrite(buffer, sizeof(int), 5, frag_out_raw);
	sfwrite(qseq, 1, q_len, frag_out_raw);
	sfwrite(header->seq, 1, header->len, frag_out_raw);
	sfwrite(start, sizeof(int), counter, frag_out_raw);
	sfwrite(end, sizeof(int), counter, frag_out_raw);
	sfwrite(template, sizeof(int), counter, frag_out_raw);
	
	/* update scores */
	if(counter == 1) { //Only one best match
		if(*template < 0) {
			*template = -*template;
		}
		alignment_scores[*template] += score;
		uniq_alignment_scores[*template] += score;
	} else {
		i = counter + 1;
		--template;
		while(--i) {
			if(*++template < 0) {
				alignment_scores[-*template] += score;
			} else {
				alignment_scores[*template] += score;
			}
		}
	}
}

void update_Scores_pe_old(unsigned char *qseq, int q_len, unsigned char *qseq_r, int qr_len, int counter, int score, int *start, int *end, int *template, Qseqs *header, Qseqs *header_r, int flag, int flag_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw) {
	
	int i, buffer[5];
	
	/* print frag */
	buffer[0] = q_len;
	buffer[1] = counter;
	buffer[2] = -score;
	buffer[3] = header->len;
	buffer[4] = flag;
	counter = abs(counter);
	sfwrite(buffer, sizeof(int), 5, frag_out_raw);
	sfwrite(qseq, 1, q_len, frag_out_raw);
	sfwrite(header->seq, 1, header->len, frag_out_raw);
	sfwrite(start, sizeof(int), counter, frag_out_raw);
	sfwrite(end, sizeof(int), counter, frag_out_raw);
	sfwrite(template, sizeof(int), counter, frag_out_raw);
	
	buffer[0] = qr_len;
	buffer[1] = header_r->len;
	buffer[2] = flag_r;
	sfwrite(buffer, sizeof(int), 3, frag_out_raw);
	sfwrite(qseq_r, 1, qr_len, frag_out_raw);
	sfwrite(header_r->seq, 1, header_r->len, frag_out_raw);
	
	/* update scores */
	if(counter == 1) { /* Only one best match */
		if(*template < 0) {
			*template = -*template;
		}
		alignment_scores[*template] += score;
		uniq_alignment_scores[*template] += score;
	} else {
		i = counter + 1;
		--template;
		while(--i) {
			if(*++template < 0) {
				alignment_scores[-*template] += score;
			} else {
				alignment_scores[*template] += score;
			}
		}
	}
}
