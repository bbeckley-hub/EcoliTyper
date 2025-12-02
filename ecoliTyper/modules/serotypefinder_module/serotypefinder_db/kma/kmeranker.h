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
#include "penalties.h"

#ifndef KMERANKER
#define KMERANKER 1;
typedef struct kmerAnker KmerAnker;
struct kmerAnker {
	int score;
	int weight;
	int score_len;
	int len_len;
	unsigned start;
	unsigned end;
	unsigned *values;
	struct kmerAnker *descend; /* descending anker */
};
#endif

extern KmerAnker * (*getChainTemplates)(KmerAnker*, const Penalties*, const int*, const int, const int, const int, int*, int*, int*, char*);
extern int (*kmerAnkerScore)(KmerAnker*);
extern const int (*testExtension)(const int, const int, const int);
extern const int (*proxiTestBest)(const double, const int, const int, const int, const int);
extern KmerAnker * (*getBestAnker)(KmerAnker**, unsigned*, const int*);
extern KmerAnker * (*getTieAnker)(int, KmerAnker*, const KmerAnker*);
int ankerScore(KmerAnker *src);
int ankerScoreLen(KmerAnker *src);
const int testExtensionScore(const int q_len, const int t_len, const int best_len);
const int testExtensionScoreLen(const int q_len, const int t_len, const int best_len);
const int proxiTestBestScore(const double proxiScore, const int score, const int q_len, const int t_len, const int best_len);
const int proxiTestBestScoreLen(const double proxiScore, const int score, const int q_len, const int t_len, const int best_len);
const int mrchain(int *bestTemaples, const int *template_lengths, const int q_len, const int maplen);
KmerAnker * getBestChainTemplates(KmerAnker *src, const Penalties *rewards, const int *template_lengths, const int q_len, const int kmersize, const int mlen, int *bests, int *Score, int *extendScore, char *include);
KmerAnker * getProxiChainTemplates(KmerAnker *src, const Penalties *rewards, const int *template_lengths, const int q_len, const int kmersize, const int mlen, int *bests, int *Score, int *extendScore, char *include);
KmerAnker * pruneAnkers(KmerAnker *V_score, int kmersize);
KmerAnker * getBestAnkerScore(KmerAnker **src, unsigned *ties, const int *template_lengths);
KmerAnker * getBestAnkerScoreLen(KmerAnker **src, unsigned *ties, const int *template_lengths);
KmerAnker * getTieAnkerScore(int stop, KmerAnker *src, const KmerAnker *bestScore);
KmerAnker * getTieAnkerScoreLen(int stop, KmerAnker *src, const KmerAnker *bestScore);
int chooseChain(const KmerAnker *bestScore, const KmerAnker *bestScore_r, int cStart, int cStart_r, int *Start, int *Len);
