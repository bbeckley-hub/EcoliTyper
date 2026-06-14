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

#include "penalties.h"

#ifndef NWLOAD
typedef struct nWmat NWmat;
typedef struct alnScore AlnScore;
typedef struct aln Aln;

struct nWmat {
	long NW_s;
	long NW_q;
	unsigned char *E;
	int *D[2];
	int *P[2];
	int *P2[2];
	Penalties *rewards;
};

struct alnScore {
	int score;
	int len;
	int pos;
	int match;
	int tGaps;
	int qGaps;
};

struct aln {
	unsigned char *t;  /* template */
	char *s;  /* score */
	unsigned char *q;  /* query */
	unsigned pos; /* start of aln, relative to template */
	int len;
	unsigned mapQ; /* mapping quality */
	int score; /* aln score */
	int start;
	int end;
};
#define NWLOAD 1
#endif

AlnScore NW(const long unsigned *template, const unsigned char *queryOrg, int k, int t_s, int t_e, int q_s, int q_e, Aln *aligned, NWmat *matrices, int template_length);
AlnScore NW_band(const long unsigned *template, const unsigned char *queryOrg, int k, int t_s, int t_e, int q_s, int q_e, Aln *aligned, int band, NWmat *matrices, int template_length);
AlnScore NW_score(const long unsigned *template, const unsigned char *queryOrg, int k, int t_s, int t_e, int q_s, int q_e, NWmat *matrices, int template_length);
AlnScore NW_band_score(const long unsigned *template, const unsigned char *queryOrg, int k, int t_s, int t_e, int q_s, int q_e, int band, NWmat *matrices, int template_length);
