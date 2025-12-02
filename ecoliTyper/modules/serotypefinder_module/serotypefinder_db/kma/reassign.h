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

#include "compdna.h"
#include "hashmapkma.h"

int reassign_template(HashMapKMA *templates, FILE *name_file, int seq_in, int *template_lengths, Assem *aligned_assem, AssemInfo *matrix, Qseqs *qseq, CompDNA *tseq);
int consensus2qseq(unsigned char *src, int len, Qseqs *dest, const unsigned char *trans);
int reassign_heapify(int *bests, int *lengths, int index);
int reassign_buildheap(int *bests, int *lengths);
int reassign_popheap(int *bests, int *lengths);
int reassign_kmers(const HashMapKMA *templates, int *bestTemplates, int *Score, int *template_lengths, CompDNA *qseq_fw, CompDNA *qseq_rc);
int reassign_getoffset(CompDNA *qseq, long unsigned kmer, int offset);
int reassign_cmpseqs(long unsigned *qseq, long unsigned *target, int len, int offset);
int reassign_testNs(int *N, int start, int end);
int reassign_matchseqs(CompDNA *consensus, CompDNA *qseq);
int reassign_matchseqs_rc(CompDNA *consensus, CompDNA *consensus_rc, CompDNA *candidate);
void reassign_loadseq(CompDNA *qseq, int len, int file, long offset);
int reassign_matrix_bias(AssemInfo *matrix, int bias);
unsigned reassign_matrix_insertions(AssemInfo *matrix);
void reassign_matrix_offset(AssemInfo *matrix, Assem *aligned_assem, int offset, long unsigned *tseq, int t_len);
