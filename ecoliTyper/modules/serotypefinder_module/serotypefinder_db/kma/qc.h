/* Philip T.L.C. Clausen Jan 2025 plan@dtu.dk */

/*
 * Copyright (c) 2025, Philip Clausen, Technical University of Denmark
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

#include <stdio.h>

#ifndef QC
typedef struct qcstat QCstat;
struct qcstat {
	long unsigned fragcount;
	long unsigned org_fragcount;
	long unsigned count;
	long unsigned org_count;
	long unsigned bpcount;
	long unsigned org_bpcount;
	long unsigned totgc;
	long unsigned totns;
	double Eeq;
	int maxlen;
	int qresolution;
	int phredScale;
	int verbose;
	unsigned *ldist;
	unsigned *qdist;
};
#define QC 1
#endif

QCstat * init_QCstat(int verbose);
void destroy_QCstat(QCstat *src);
int rescale_ldist(unsigned *ldist, const int maskold, const int maxlen);
void update_QCstat(QCstat *src, int len, int gc, int ns, double sp);
int print_QCstat(QCstat *src, int minQ, int minPhred, int minmaskQ, int minlen, int maxlen, int fiveClip, int threeClip, FILE *dest);
