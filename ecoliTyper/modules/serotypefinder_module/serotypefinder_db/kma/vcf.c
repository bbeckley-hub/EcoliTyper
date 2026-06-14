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
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "assembly.h"
#include "filebuff.h"
#include "stdnuc.h"
#include "stdstat.h"
#include "version.h"
#include "vcf.h"

char * noFolder(const char *src) {
	
	int pos;
	
	pos = strlen(src) - 1;
	while(pos && src[pos] != '/') {
		--pos;
	}
	if(pos) {
		++pos;
	}
	
	return ((char *) src) + pos;
}

void initialiseVcf(FileBuff *fileP, char *templateFilename) {
	
	unsigned check, avail;
	char *update;
	
	update = (char *) fileP->next;
	avail = fileP->bytes;
	
	/* header stuff */
	check = sprintf(update, "##fileformat=VCFv4.2\n");
	update += check; avail -= check;
	check = sprintf(update, "##kmaVersion=%s\n", KMA_VERSION);
	update += check; avail -= check;
	
	/* filter stuff */
	check = sprintf(update, "##FILTER=<ID=LowQual,Description=\"Low quality\">\n");
	update += check; avail -= check;
	
	/* INFO stuff */
	check = sprintf(update, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
	update += check; avail -= check;
	check = sprintf(update, "##INFO=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n");
	update += check; avail -= check;
	check = sprintf(update, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Fraction\">\n");
	update += check; avail -= check;
	check = sprintf(update, "##INFO=<ID=RAF,Number=1,Type=Float,Description=\"Revised Allele Fraction\">\n");
	update += check; avail -= check;
	check = sprintf(update, "##INFO=<ID=DEL,Number=1,Type=Float,Description=\"Fraction of Reads Containing Spanning Deletions\">\n");
	update += check; avail -= check;
	check = sprintf(update, "##INFO=<ID=AD6,Number=6,Type=Integer,Description=\"Count of all alternative alleles: A,C,G,T,N,-\">\n");
	update += check; avail -= check;
	
	/* FORMAT stuff */
	check = sprintf(update, "##FORMAT=<ID=Q,Number=1,Type=Float,Description=\"McNemar quantile\">\n");
	update += check; avail -= check;
	check = sprintf(update, "##FORMAT=<ID=P,Number=1,Type=Float,Description=\"McNemar p-value\">\n");
	update += check; avail -= check;
	check = sprintf(update, "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filter\">\n");
	update += check; avail -= check;
	if(templateFilename) {
		check = sprintf(update, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", noFolder(templateFilename));
	} else {
		check = sprintf(update, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tspltDB\n");
	}
	update += check; avail -= check;
	
	fileP->next = (unsigned char *) update;
	fileP->bytes = avail;
	
}

void updateVcf(char *template_name, unsigned char *template_seq, double evalue, double support, int bcd, int t_len, AssemInfo *matrix, int filter, FileBuff *fileP) {
	
	static const char nuc2num[256] = {
		8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		8, 0, 8, 1, 8, 8, 8, 2, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 8, 8, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		8, 0, 8, 1, 8, 8, 8, 2, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 8, 8, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};
	const char *PASS = "PASS", *FAIL = "FAIL", *LowQual = "LowQual", *UNKNOWN = ".";
	int i, pos, nextPos, bestScore, depthUpdate, bestBaseScore, nucNum;
	int template_name_length, check, avail, DP, AD, DEL, QUAL;
	double AF, RAF, Q, P;
	const double lnConst = -10 / log(10);
	char *FILTER, **FILTER_ptr, *update;
	unsigned char nuc, bestNuc;
	const char bases[6] = "ACGTN-";
	Assembly *assembly;
	
	if(filter == 2) {
		FILTER_ptr = &FILTER;
	} else {
		FILTER_ptr = (char **) &UNKNOWN;
	}
	template_name_length = strlen(template_name);
	update = (char *) fileP->next;
	avail = fileP->bytes;
	assembly = matrix->assmb;
	nextPos = 0;
	do {
		pos = nextPos;
		nextPos = assembly[pos].next;
		
		/* handle insertions now */
		nuc = *template_seq++;
		if(pos < t_len) {
			++i;
		} else if(nuc != '-') {
			--template_seq;
			nuc = '-';
		}
		
		/* call base */
		bestNuc = nuc2num[nuc];
		bestScore = assembly[pos].counts[bestNuc];
		depthUpdate = 0;
		for(i = 0; i < 6; ++i) {
			if(bestScore < assembly[pos].counts[i]) {
				bestScore = assembly[pos].counts[i];
				bestNuc = i;
			}
			depthUpdate += assembly[pos].counts[i];
		}
		nucNum = bestNuc;
		bestNuc = bases[bestNuc];
		
		/* Check for minor base call */
		if(!depthUpdate) {
			nucNum = 5;
			bestNuc = '-';
		} else if((bestScore << 1) < depthUpdate) {
			if(bestNuc == '-') {
				bestBaseScore = 0;
				bestNuc = 4;
				for(i = 0; i < 5; ++i) {
					if(bestBaseScore < assembly[pos].counts[i]) {
						bestBaseScore = assembly[pos].counts[i];
						bestNuc = i;
					}
				}
				nucNum = bestNuc;
				bestNuc = tolower(bases[bestNuc]);
			} else {
				bestNuc = tolower(bestNuc);
			}
			bestScore = depthUpdate - assembly[pos].counts[5];
		} else if(depthUpdate < bcd) {
			/* too low depth */
			bestNuc = tolower(bestNuc);
		}
		
		if(bestScore) {
			/* determine base at current position */
			bestNuc = baseCall(bestNuc, nuc, bestScore, depthUpdate, evalue, &assembly[pos]);
			nucNum = nuc2num[bestNuc];
			
			/* INFO */
			DP = depthUpdate;
			AD = assembly[pos].counts[nucNum];
			AF = (double) AD / DP;
			RAF = (double) bestScore / DP;
			DEL = assembly[pos].counts[5];
			Q = pow(depthUpdate - (bestScore << 1), 2) / depthUpdate;
			P = p_chisqr(Q);
			
			/* discard unimportant changes */
			if(nuc != bestNuc || (t_len <= nextPos && *template_seq == '-') || DP < bcd || evalue < P || AD < support * DP) {
				/* QUAL */
				//QUAL = lnConst * log(P);
				QUAL = lnConst * log(binP(DP, AD, 0.25));
				QUAL = (QUAL < 0 || 3079 < QUAL) ? 3079 : QUAL;
				
				/* FILTER */
				if(bcd <= DP && P <= evalue && support * DP <= AD) {
					FILTER = (char *) PASS;
				} else if(bcd <= DP || P <= evalue || support * DP <= AD) {
					FILTER = (char *) LowQual;
				} else {
					FILTER = (char *) FAIL;
				}
				
				if(avail < template_name_length + 167) {
					fileP->bytes = avail;
					writeGzFileBuff(fileP);
					avail = fileP->bytes;
					update = (char *) fileP->next;
				}
				
				strcpy(update, template_name);
				update += template_name_length; avail -= template_name_length;
				*update++ = '\t'; --avail;
				
				if(pos < t_len) {
					check = sprintf(update, "%d", pos + 1);
					update += check; avail -= check;
				} else {
					*update++ = '0';
					--avail;
				}
				*update++ = '\t';
				*update++ = '.';
				*update++ = '\t';
				avail -= 3;
				
				if(nuc != '-') {
					*update++ = nuc;
					--avail;
				} else {
					*update++ = '<';
					*update++ = nuc;
					*update++ = '>';
					avail -= 3;
				}
				*update++ = '\t';
				--avail;
				if(bestNuc == '-') {
					*update++ = '<';
					*update++ = bestNuc;
					*update++ = '>';
					avail -= 3;
				} else {
					*update++ = bestNuc;
					--avail;
				}
				
				check = sprintf(update, "\t%d\t%s\tDP=%d;AD=%d;AF=%.2f;RAF=%.2f;DEL=%d;", QUAL, *FILTER_ptr, DP, AD, AF, RAF, DEL);
				update += check; avail -= check;
				check = sprintf(update, "AD6=%d,%d,%d,%d,%d,%d\t", assembly[pos].counts[0], assembly[pos].counts[1], assembly[pos].counts[2], assembly[pos].counts[3], assembly[pos].counts[4], assembly[pos].counts[5]);
				update += check; avail -= check;
				check = sprintf(update, "Q:P:FT\t%.2f:%4.1e:%s\n", Q, P, FILTER);
				update += check; avail -= check;
			}
		} else if(pos < t_len) {
			FILTER = (char *) FAIL;
			if(avail < template_name_length + 105) {
				fileP->bytes = avail;
				writeGzFileBuff(fileP);
				avail = fileP->bytes;
				update = (char *) fileP->next;
			}
			check = sprintf(update, "%s\t%d\t.\t%c\t%c\t%d\t%s\t", template_name, pos + 1, nuc, '.', 0, *FILTER_ptr);
			update += check; avail -= check;
			check = sprintf(update, "DP=%d;AD=%d;AF=%.2f;RAF=%.2f;DEL=%d;AD6=%d,%d,%d,%d,%d,%d\t", 0, 0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0);
			update += check; avail -= check;
			check = sprintf(update, "Q:P:FT\t%.2f:%4.1e:%s\n", 0.0, 1.0, FILTER);
			update += check; avail -= check;
		
		}
	} while(nextPos != 0);
	
	fileP->next = (unsigned char *) update;
	fileP->bytes = avail;
	
}
