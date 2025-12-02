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
#include "frags.h"
#include "qseqs.h"

extern int (*ConClavePtr)(FILE *, FILE ***, int, int, long unsigned *, unsigned *, unsigned *, long unsigned *, long unsigned *, int *, Qseqs *, Qseqs *, int *, int *, int *, Frag **);
extern int (*ConClave2Ptr)(FILE *, FILE ***, int, int, long unsigned *, unsigned *, unsigned *, long unsigned *, long unsigned *, int *, Qseqs *, Qseqs *, int *, int *, int *, Frag **, long unsigned, double, double);

unsigned char * ustrdup(unsigned char *src, size_t n);
int runConClave(FILE *frag_in_raw, FILE ***Template_fragments, int DB_size, int maxFrag, long unsigned *w_scores, unsigned *fragmentCounts, unsigned *readCounts, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *template_lengths, Qseqs *header, Qseqs *qseq, int *bestTemplates, int *best_start_pos, int *best_end_pos, Frag **alignFrags);
int runConClave_lc(FILE *frag_in_raw, FILE ***Template_fragments, int DB_size, int maxFrag, long unsigned *w_scores, unsigned *fragmentCounts, unsigned *readCounts, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *template_lengths, Qseqs *header, Qseqs *qseq, int *bestTemplates, int *best_start_pos, int *best_end_pos, Frag **alignFrags);
int runConClave2(FILE *frag_in_raw, FILE ***Template_fragments, int DB_size, int maxFrag, long unsigned *w_scores, unsigned *fragmentCounts, unsigned *readCounts, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *template_lengths, Qseqs *header, Qseqs *qseq, int *bestTemplates, int *best_start_pos, int *best_end_pos, Frag **alignFrags, long unsigned template_tot_ulen, double scoreT, double evalue);
int runConClave2_lc(FILE *frag_in_raw, FILE ***Template_fragments, int DB_size, int maxFrag, long unsigned *w_scores, unsigned *fragmentCounts, unsigned *readCounts, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *template_lengths, Qseqs *header, Qseqs *qseq, int *bestTemplates, int *best_start_pos, int *best_end_pos, Frag **alignFrags, long unsigned template_tot_ulen, double scoreT, double evalue);
