/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#define _XOPEN_SOURCE 600
#include <stdio.h>
#include "compdna.h"
#include "hashmap.h"

extern int (*update_DB)(HashMap *, CompDNA *, unsigned, int, double, double, unsigned *, unsigned *, Qseqs *);
extern void (*updateAnnotsPtr)(CompDNA *, int, int, FILE *, unsigned **, unsigned **, unsigned **);
int updateDBs(HashMap *templates, CompDNA *qseq, unsigned template, int MinKlen, double homQ, double homT, unsigned *template_ulengths, unsigned *template_slengths, Qseqs *header);
int updateDBs_sparse(HashMap *templates, CompDNA *qseq, unsigned template, int MinKlen, double homQ, double homT, unsigned *template_ulengths, unsigned *template_slengths, Qseqs *header);
void updateAnnots(CompDNA *qseq, int DB_size, int kmerindex, FILE *seq_out, unsigned **template_lengths, unsigned **template_ulengths, unsigned **template_slengths);
void updateAnnots_sparse(CompDNA *qseq, int DB_size, int kmerindex, FILE *seq_out, unsigned **template_lengths, unsigned **template_ulengths, unsigned **template_slengths);
void dumpSeq(CompDNA *qseq, int kmerindex, FILE *seq_out);
void makeIndexing(CompDNA *compressor, int kmerindex, FILE *seq_out);
