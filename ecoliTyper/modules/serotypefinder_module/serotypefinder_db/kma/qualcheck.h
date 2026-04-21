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
#include "compdna.h"
#include "hashmap.h"

extern int (*QualCheck)(HashMap *templates, CompDNA *, int, double, double, unsigned *, Qseqs *);
int lengthCheck(HashMap *templates, CompDNA *qseq, int MinKlen, double homQ, double homT, unsigned *template_ulengths, Qseqs *header);
int queryCheck(HashMap *templates, CompDNA *qseq, int MinKlen, double homQ, double homT, unsigned *template_ulengths, Qseqs *header);
int templateCheck(HashMap *templates, CompDNA *qseq, int MinKlen, double homQ, double homT, unsigned *template_ulengths, Qseqs *header);
void updateScoreAndTemplate(unsigned *Scores_tot, unsigned *bestTemplates, unsigned *values);
void updateScoreAndTemplateHU(unsigned *Scores_tot, unsigned *bestTemplates, unsigned *values_org);
void addUscore(unsigned *Scores, unsigned *values);
void addUscoreHU(unsigned *Scores, unsigned *values_org);
