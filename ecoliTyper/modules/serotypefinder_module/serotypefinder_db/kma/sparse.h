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
#include "compkmers.h"
#include "hashmapkma.h"
#include "hashtable.h"
#include "qc.h"

int translateToKmersAndDump(long unsigned *Kmers, int n, int max, unsigned char *qseq, int seqlen, const HashMapKMA *templates, FILE *out);
char ** load_DBs_Sparse(char *templatefilename, unsigned **template_lengths, unsigned **template_ulengths, unsigned shm);
void save_kmers_sparse(const HashMapKMA *templates, HashMap_kmers *foundKmers, CompKmers *compressor);
void run_input_sparse(const HashMapKMA *templates, char **inputfiles, int fileCount, int minPhred, int hardmaskQ, int minQ, int fiveClip, int threeClip, int minlen, int maxlen, int kmersize, char *trans, const double *prob, QCstat *qcreport, FILE *out);
int save_kmers_sparse_batch(char *templatefilename, char *outputfilename, char *exePrev, double ID_t, double Depth_t, double evalue, char ss, unsigned shm);
