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
#include "penalties.h"
#include "qseqs.h"

#define nameSkip(infile, c) while((c = fgetc(infile)) != '\n' && c != EOF)

int load_DBs_KMA(char *templatefilename, long unsigned **alignment_scores, long unsigned **uniq_alignment_scores, int **template_lengths, unsigned shm);
char * nameLoad(Qseqs *name, FILE *infile);
int runKMA(char *templatefilename, char *outputfilename, char *exePrev, int ConClave, int kmersize, int minlen, Penalties *rewards, int extendedFeatures, double ID_t, double Depth_t, int mq, double scoreT, double mrc, double minFrac, double evalue, double support, int bcd, int ref_fsa, int print_matrix, int print_all, long unsigned tsv, int vcf, int xml, int sam, int nc, int nf, unsigned shm, int thread_num, int maxFrag, int verbose, unsigned preset);
int runKMA_MEM(char *templatefilename, char *outputfilename, char *exePrev, int ConClave, int kmersize, int minlen, Penalties *rewards, int extendedFeatures, double ID_t, double Depth_t, int mq, double scoreT, double mrc, double minFrac, double evalue, double support, int bcd, int ref_fsa, int print_matrix, int print_all, long unsigned tsv, int vcf, int xml, int sam, int nc, int nf, unsigned shm, int thread_num, int maxFrag, int verbose, unsigned preset);
