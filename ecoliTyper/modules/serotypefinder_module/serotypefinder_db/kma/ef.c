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
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "assembly.h"
#include "ef.h"
#include "stdnuc.h"
#include "threader.h"
#include "vcf.h"
#include "version.h"

void initExtendedFeatures(FILE *out, char *templatefilename, unsigned totFrags, char *cmd) {
	
	char Date[11];
	time_t t1;
	struct tm *tm;
	
	fprintf(out, "## method\tKMA\n");
	fprintf(out, "## version\t%s\n", KMA_VERSION);
	fprintf(out, "## database\t%s\n", noFolder(templatefilename));
	fprintf(out, "## fragmentCount\t%u\n", totFrags);
	time(&t1);
	tm = localtime(&t1);
	strftime(Date, sizeof(Date), "%Y-%m-%d", tm);
	fprintf(out, "## date\t%s\n", Date);
	fprintf(out, "## command\t%s\n", cmd);
	fprintf(out, "# refSequence\treadCount\tfragmentCount\tmapScoreSum\trefCoveredPositions\trefConsensusSum\tbpTotal\tdepthVariance\tnucHighDepthVariance\tdepthMax\tsnpSum\tinsertSum\tdeletionSum\treadCountAln\tfragmentCountAln\n");
}

void getExtendedFeatures(Assem *aligned_assem, AssemInfo *matrix, long unsigned *seq, int t_len, int thread_num) {
	
	static volatile int Lock = 0, next, thread_wait = 0;
	volatile int *excludeMatrix = &Lock;
	unsigned pos, end, chunk, nucHighVar, maxDepth, depthUpdate;
	long unsigned snpSum, insertSum, deletionSum;
	double highVar;
	Assembly *assembly;
	
	/* init */
	assembly = matrix->assmb;
	maxDepth = 0;
	nucHighVar = 0;
	snpSum = 0;
	insertSum = 0;
	deletionSum = 0;
	chunk = 8112;
	lock(excludeMatrix);
	if(!thread_wait) {
		next = 0;
		thread_wait = thread_num;
		aligned_assem->fragmentCountAln = (((aligned_assem->readCountAln >> 1) + (aligned_assem->readCountAln & 1)) <= aligned_assem->fragmentCountAln) ? (aligned_assem->fragmentCountAln) : ((aligned_assem->readCountAln >> 1) + (aligned_assem->readCountAln & 1));
		aligned_assem->nucHighVar = 0;
		aligned_assem->maxDepth = 0;
		aligned_assem->snpSum = 0;
		aligned_assem->insertSum = 0;
		aligned_assem->deletionSum = 0;
	}
	unlock(excludeMatrix);
	
	/* get var */
	fixVarOverflow(aligned_assem, assembly, t_len, thread_num);
	highVar = (long double)(aligned_assem->depth) / t_len + 3 * sqrt(aligned_assem->var);
	
	/* get nucHighVar */
	while(chunk) {
		lock(excludeMatrix);
		/* get next chunk */
		pos = next;
		if((next += chunk) < 0) {
			next = t_len;
		}
		unlock(excludeMatrix);
		
		/* call chunk */
		if(pos < t_len) {
			end = pos + chunk;
			if(t_len <= end) {
				end = t_len - 1;
				chunk = 0;
			}
			while(pos != end) {
				depthUpdate = assembly[pos].counts[0] + assembly[pos].counts[1] + assembly[pos].counts[2] + assembly[pos].counts[3] + assembly[pos].counts[4];
				
				if(pos < t_len) {
					deletionSum += assembly[pos].counts[5];
					snpSum += (depthUpdate - assembly[pos].counts[getNuc(seq, pos)]);
				} else {
					insertSum += depthUpdate;
				}
				
				depthUpdate += assembly[pos].counts[5];
				
				if(maxDepth < depthUpdate) {
					maxDepth = depthUpdate;
				}
				if(highVar < depthUpdate) {
					++nucHighVar;
				}
				if(!(pos = assembly[pos].next)) {
					pos = end;
				}
			}
		} else {
			chunk = 0;
		}
	}
	
	lock(excludeMatrix);
	aligned_assem->nucHighVar += nucHighVar;
	if(aligned_assem->maxDepth < maxDepth) {
		aligned_assem->maxDepth = maxDepth;
	}
	aligned_assem->snpSum += snpSum;
	aligned_assem->insertSum += insertSum;
	aligned_assem->deletionSum += deletionSum;
	--thread_wait;
	unlock(excludeMatrix);
	wait_atomic(thread_wait);
}

void printExtendedFeatures(char *template_name, Assem *aligned_assem, unsigned fragmentCount, unsigned readCount, FILE *outfile) {
	
	if(aligned_assem) {
		fprintf(outfile, "%s\t%u\t%u\t%lu\t%u\t%u\t%lu\t%f\t%u\t%u\t%lu\t%lu\t%lu\t%u\t%u\n", template_name, readCount, fragmentCount, aligned_assem->score, aligned_assem->aln_len, aligned_assem->cover, aligned_assem->depth, aligned_assem->var, aligned_assem->nucHighVar, aligned_assem->maxDepth, aligned_assem->snpSum, aligned_assem->insertSum, aligned_assem->deletionSum, aligned_assem->readCountAln, aligned_assem->fragmentCountAln);
	} else {
		fprintf(outfile, "%s\t%u\t%u\t%u\t%u\t%u\t%u\t%f\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n", template_name, 0, 0, 0, 0, 0, 0, 0.0, 0, 0, 0, 0, 0, 0, 0);
	}
}
