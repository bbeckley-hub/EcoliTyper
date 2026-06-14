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
#include "qseqs.h"

void update_Scores_MEM(unsigned char *qseq, int q_len, int counter, int score, int *start, int *end, int *template, Qseqs *header, int flag, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw);
void update_Scores_pe_MEM(unsigned char *qseq, int q_len, unsigned char *qseq_r, int qr_len, int counter, int score, int *start, int *end, int *template, Qseqs *header, Qseqs *header_r, int flag, int flag_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw);
int update_Scores(unsigned char *qseq, int q_len, double minFrac, int counter, int bestReadScore, double bestScore, int *start, int *end, int *templates, int *Scores, int *Lengths, Qseqs *header, int flag, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw);
int update_Scores_se(unsigned char *qseq, int q_len, double minFrac, int counter, int bestScore, int *start, int *end, int *templates, int *Scores, Qseqs *header, int flag, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw);
int update_Scores_pe(unsigned char *qseq, int q_len, unsigned char *qseq_r, int qr_len, double minFrac, int counter, int bestScore, int *start, int *end, int *templates, int *Scores, Qseqs *header, Qseqs *header_r, int flag, int flag_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw);
void update_Scores_old(unsigned char *qseq, int q_len, int counter, int score, int *start, int *end, int *template, Qseqs *header, int flag, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw);
void update_Scores_pe_old(unsigned char *qseq, int q_len, unsigned char *qseq_r, int qr_len, int counter, int score, int *start, int *end, int *template, Qseqs *header, Qseqs *header_r, int flag, int flag_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw);
