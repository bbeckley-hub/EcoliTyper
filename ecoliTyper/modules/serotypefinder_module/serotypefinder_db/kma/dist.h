/* Philip T.L.C. Clausen Jan 2020 plan@dtu.dk */

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
#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hashmapkma.h"
#include "matrix.h"
#include "pherror.h"
#include "runkma.h"

#ifndef DIST
typedef struct distThread DistThread;
struct distThread {
	pthread_t id;
	int flag;
	int format;
	int thread_num;
	int *N;
	long unsigned ltdMat;
	long unsigned covMat;
	char *outfileM;
	FILE *name_file;
	HashMapKMA *DB;
	Matrix *Dist;
	Qseqs *template_name;
	struct distThread *next;
};
#define DIST 1
#endif

HashMapKMA * loadValues(const char *filename);
void destroyValues(HashMapKMA *src);
void kmerSimilarity(HashMapKMA *DB, Matrix *Dist, int *N);
void kmerSimilarity_thread(HashMapKMA *DB, Matrix *Dist, int *N, int thread_num, volatile int *lock);
int kmerDist(int Ni, int Nj, int D);
int kmerShared(int Ni, int Nj, int D);
int chi2dist(int Ni, int Nj, int D);
void printIntLtdPhy(char *outfile, Matrix *Dist, int *N, FILE *name_file, Qseqs *template_name, unsigned format, int thread_num, volatile int *lock, const char *method, int (*distPtr)(int, int, int));
double kmerQuery(int Ni, int Nj, int D);
double kmerTemplate(int Ni, int Nj, int D);
double kmerAvg(int Ni, int Nj, int D);
double kmerInvAvg(int Ni, int Nj, int D);
double kmerJaccardDist(int Ni, int Nj, int D);
double kmerJaccardSim(int Ni, int Nj, int D);
double kmerCosineDist(int Ni, int Nj, int D);
double kmerCosineSim(int Ni, int Nj, int D);
double kmerOverlapCoef(int Ni, int Nj, int D);
double kmerInvOverlapCoef(int Ni, int Nj, int D);
void printDoublePhy(char *outfile, Matrix *Dist, int *N, FILE *name_file, Qseqs *template_name, unsigned format, const char *formatString, int ltd, int thread_num, volatile int *lock, const char *method, double (*distPtr)(int, int, int));
long unsigned getPhySize(int flag, int format, long unsigned n, long unsigned *ltdMat, long unsigned *covMat, FILE *name_file);
char * mfile(FILE *outfile, long unsigned size);
void * threadDist(void *arg);
void runDist(char *templatefilename, char *outputfilename, int flag, int format, int disk, int thread_num);
int dist_main(int argc, char *argv[]);
