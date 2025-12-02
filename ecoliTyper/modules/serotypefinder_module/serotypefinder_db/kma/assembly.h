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
#include <pthread.h>
#include <stdio.h>
#include "chain.h"
#include "filebuff.h"
#include "hashmapcci.h"
#include "nw.h"
#include "qseqs.h"

#ifndef ASSEMBLY
typedef struct assem Assem;
typedef struct assembly Assembly;
typedef struct assemInfo AssemInfo;
typedef struct assemble_thread Assemble_thread;

struct assem {
	unsigned char *t;  /* template */
	char *s;  /* score */
	unsigned char *q;  /* query */
	long unsigned depth;
	long unsigned depthVar;
	long unsigned score;
	long unsigned snpSum;
	long unsigned insertSum;
	long unsigned deletionSum;
	unsigned cover;
	unsigned len;
	unsigned aln_len;
	unsigned size;
	unsigned fragmentCountAln;
	unsigned readCountAln;
	unsigned nucHighVar;
	unsigned maxDepth;
	double var;
};

struct assembly {
	short unsigned counts[6];
	unsigned next;
};

struct assemInfo {
	int len;
	int size;
	struct assembly *assmb;
};

struct assemble_thread {
	pthread_t id;
	int num;
	int template;
	int file_count;
	int spin;
	int mq;
	int minlen;
	int bcd;
	int sam;
	int ef;
	int t_len;
	int seq_in;
	int kmersize;
	int thread_num;
	double scoreT;
	double mrc;
	double evalue;
	char *template_name;
	FILE **files;
	FILE *xml_out;
	FileBuff *frag_out;
	Assem *aligned_assem;
	Aln *aligned, *gap_align;
	Qseqs *qseq, *header;
	AssemInfo *matrix;
	AlnPoints *points;
	NWmat *NWmatrices;
	HashMapCCI *template_index;
	Assemble_thread *next;
};
#define ASSEMBLY 1
#endif

extern void * (*assembly_KMA_Ptr)(void *);
extern int (*significantBase)(int, int, double);
extern unsigned char (*baseCall)(unsigned char, unsigned char, int, int, double, Assembly*);
extern void (*alnToMatPtr)(AssemInfo *, Assem *, Aln *, AlnScore, int, int);
void updateMatrix(FileBuff *dest, char *template_name, long unsigned *template_seq, AssemInfo *matrix, int t_len);
int significantNuc(int X, int Y, double evalue);
int significantAnd90Nuc(int X, int Y, double evalue);
int significantAndSupport(int X, int Y, double evalue);
unsigned char baseCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls);
unsigned char orgBaseCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls);
unsigned char refCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls);
unsigned char nanoCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls);
unsigned char refNanoCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls);
void * assemble_KMA_threaded(void *arg);
void * assemble_KMA_dense_threaded(void *arg);
void * skip_assemble_KMA(void *arg);
void alnToMat(AssemInfo *matrix, Assem *aligned_assem, Aln *aligned, AlnScore alnStat, int t_len, int flag);
void alnToMatDense(AssemInfo *matrix, Assem *aligned_assem, Aln *aligned, AlnScore alnStat, int t_len, int flag);
void callConsensus(AssemInfo *matrix, Assem *aligned_assem, long unsigned *seq, int t_len, int bcd, double evalue, int thread_num);
void fixVarOverflow(Assem *aligned_assem, Assembly *assembly, int t_len, int thread_num);
void * assemble_KMA(void *arg);
void complement_swap(short unsigned *c1, short unsigned *c2);
void complement_pos(short unsigned *counts);
unsigned rc_insertion(Assembly *assembly, const unsigned new, const unsigned org);
void assemble_rc(AssemInfo *matrix, int t_len);
int assem_rc(Assem *aligned_assem, unsigned char *trans);
