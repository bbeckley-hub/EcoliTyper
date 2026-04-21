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
#include "compdna.h"
#include "qseqs.h"

#ifndef SPLTDB
typedef struct spltDBbuff SpltDBbuff;
struct spltDBbuff {
	int num;
	int size;
	unsigned char *buff;
	struct spltDBbuff *next;
	struct spltDBbuff *prev;
};
#define SPLTDB 1
#endif

int print_ankers_spltDB(int *out_Tem, CompDNA *qseq, int rc_flag, const Qseqs *header, const int flag, FILE *out);
int print_ankers_Sparse_spltDB(int *out_Tem, CompDNA *qseq, int rc_flag, const Qseqs *header, const int flag, FILE *out);
unsigned get_ankers_spltDB(int *infoSize, int *out_Tem, CompDNA *qseq, Qseqs *header, int *flag, FILE *inputfile);
int runKMA_spltDB(char **templatefilenames, int targetNum, char *outputfilename, int argc, char **argv, int ConClave, int kmersize, int minlen, Penalties *rewards, int extendedFeatures, double ID_t, double Depth_t, int mq, double scoreT, double mrc, double evalue, double support, int bcd, int ref_fsa, int print_matrix, int print_all, long unsigned tsv, int vcf, int xml, int sam, int nc, int nf, unsigned shm, int thread_num, int maxFrag, int verbose);
