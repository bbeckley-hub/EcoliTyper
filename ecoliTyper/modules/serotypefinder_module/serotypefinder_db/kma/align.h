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
#include "chain.h"
#include "compdna.h"
#include "hashmapcci.h"
#include "nw.h"

extern AlnScore (*leadTailAlnPtr)(Aln *, Aln *, const long unsigned*, const unsigned char*, int, int, int, const int, NWmat *);
extern void (*trailTailAlnPtr)(Aln *, Aln *, AlnScore *, const long unsigned *, const unsigned char *, int, int, int, int, const int, NWmat *);

AlnScore skipLeadAln(Aln *aligned, Aln *Frag_align, const long unsigned *tseq, const unsigned char *qseq, int t_e, int t_len, int q_e, const int bandwidth, NWmat *matrices);
AlnScore leadTailAln(Aln *aligned, Aln *Frag_align, const long unsigned *tseq, const unsigned char *qseq, int t_e, int t_len, int q_e, const int bandwidth, NWmat *matrices);
void skipTrailAln(Aln *aligned, Aln *Frag_align, AlnScore *Stat, const long unsigned *tseq, const unsigned char *qseq, int t_s, int t_len, int q_s, int q_len, const int bandwidth, NWmat *matrices);
void trailTailAln(Aln *aligned, Aln *Frag_align, AlnScore *Stat, const long unsigned *tseq, const unsigned char *qseq, int t_s, int t_len, int q_s, int q_len, const int bandwidth, NWmat *matrices);
AlnScore KMA(const HashMapCCI *template_index, const unsigned char *qseq, int q_len, int q_start, int q_end, Aln *aligned, Aln *Frag_align, int min, int max, int mq, double scoreT, AlnPoints *points, NWmat *matrices);
AlnScore KMA_score(const HashMapCCI *template_index, const unsigned char *qseq, int q_len, int q_start, int q_end, const CompDNA *qseq_comp, int mq, double scoreT, AlnPoints *points, NWmat *matrices);
int preseed(const HashMapCCI *template_index, unsigned char *qseq, int q_len);
void intcpy(int *dest, int *src, int size);
int anker_rc(const HashMapCCI *template_index, unsigned char *qseq, int q_len, int q_start, int q_end, AlnPoints *points);
int anker_rc_comp(const HashMapCCI *template_index, unsigned char *qseq, unsigned char *qseq_r, CompDNA *qseq_comp, CompDNA *qseq_r_comp, int q_start, int q_end, AlnPoints *points);
