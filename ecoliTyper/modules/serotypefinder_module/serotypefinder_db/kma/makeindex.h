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
#include "hashmap.h"

extern int (*biasPrintPtr)(FILE*, char*, unsigned char*, int);
extern int (*qualcheck)(CompDNA *, const int);
int biasPrint(FILE *name_out, char *format, unsigned char *name, int bias);
int biasNoPrint(FILE *name_out, char *format, unsigned char *name, int bias);
int lenCheck(CompDNA *qseq, const int minlen);
int internalStopCheck1(CompDNA *qseq);
int internalStopCheck(CompDNA *qseq, const int minlen);
int qualCheck(CompDNA *qseq, const int minlen);
void makeDB(HashMap *templates, int kmerindex, char **inputfiles, int fileCount, char *outputfilename, int appender, char *trans, int MinLen, int MinKlen, double homQ, double homT, unsigned **template_lengths, unsigned **template_ulengths, unsigned **template_slengths);
