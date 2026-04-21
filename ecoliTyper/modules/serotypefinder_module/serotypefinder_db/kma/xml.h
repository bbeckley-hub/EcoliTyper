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

void initXML(FILE *out, const char *templatefilename, const unsigned totFrags, int argc, char **argv);
void capXML(FILE *out);
FILE * openInitXML(const char *filename, const char *templatefilename, const unsigned totFrags, int argc, char **argv);
void closeCapXML(FILE *out);
void newIterXML(FILE *out, const int template, const int t_len, const char *template_name);
double getEntropy(const unsigned char *aligned_assem_q, const int len);
void capIterXML(FILE *out, const int DB_size, const long unsigned seqsize, const int t_len, const int readCounts, const double p_value, const long read_score, const unsigned char *aligned_assem_q, const int len);
void hitXML(FILE *out, const int template, const unsigned char *template_name, const Aln *aligned, const AlnScore *alnStat, const Penalties *rewards, const int flag);
