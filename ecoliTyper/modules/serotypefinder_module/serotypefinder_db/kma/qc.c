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
#include <stdlib.h>
#include <math.h>
#include "pherror.h"
#include "qc.h"

QCstat * init_QCstat(int verbose) {
	
	QCstat *dest;
	
	dest = calloc(sizeof(QCstat) + 256 * sizeof(unsigned), 1);
	if(!dest) {
		return 0;
	}
	dest->verbose = verbose;
	dest->qdist = (unsigned *)(dest + 1);
	dest->ldist = calloc(512, sizeof(unsigned));
	if(!dest->ldist) {
		free(dest);
		return 0;
	}
	
	return dest;
}

void destroy_QCstat(QCstat *src) {
	free(src->ldist);
	free(src);
}

int rescale_ldist(unsigned *ldist, const int maskold, const int maxlen) {
	
	unsigned i, masknew, mask, *ptr;
	
	masknew = maskold;
	while(512 <= (maxlen >> ++masknew));
	mask = masknew - maskold;
	i = 0;
	ptr = ldist;
	while(++i < 512) {
		ldist[i >> mask] += *++ptr;
		*ptr = 0;
	}
	
	return masknew;
}

unsigned * rescale_ldist_v1(unsigned *ldist, const int size, const int maxlen) {
	
	unsigned i, *ptr;
	
	ptr = realloc(ldist, (maxlen + 4) * sizeof(unsigned));
	if(!ptr) {
		ERROR();
	}
	ldist = ptr;
	i = maxlen;
	ptr += size - 1;
	while(i <= size) {
		*++ptr = 0;
	}
	
	return ldist;
}

void update_QCstat(QCstat *src, int len, int gc, int ns, double sp) {
	
	src->count++;
	src->bpcount += len;
	src->totgc += gc;
	src->totns += ns;
	src->Eeq += sp;
	if(src->maxlen < len) {
		if(!src->verbose) {
			if(512 <= (len >> src->qresolution)) {
				src->qresolution = rescale_ldist(src->ldist, src->qresolution, len);
			}
		} else { /* currently v >= 1 */
			src->ldist = rescale_ldist_v1(src->ldist, src->maxlen, len);
		}
		src->maxlen = len;
	}
	src->qdist[(int)(ceil(-10 * log10(sp / len)))]++;
	src->ldist[len >> src->qresolution]++;
}

/* adapter trim */
/*
SE, both ends (unless sr)
PE, 3' and only if overlapping.


Build directed de bruijn graph from the end of the sequence where the 
adaptors are expected, from the reads in the first buffer.
Trim the ends of the graph for low frequency k-mers, to get rid of non-adaptor 
sequences. There shouled be a frequency peak somewhere in the "middle", as 
the start of the adaptors are more frequent. I.e. do not trim nodes that are 
not pointed to.
	a. If known adaptors are given, mark terminating nodes and skip low-freq trimming.

Check for adaptors by searching from the end and towards the center of the read.
	a. Mismatches are handles by substituting the next nuc by the graph nuc.
	b. Insertion in seq, omit next nuc in k-mer
	c. Deletion in seq, add nuc from next in graph

Errors in fist n-k seeds can be handled as in BLAST.

For this we need:
	k: k-mer size.
	len: len of area in ends to build graph over.
	size: buffer size.
	edit: max edit dist.
	freq: low frequency threshold, to avoid removing already remoived adaptors.
	(fsa): known adaptor sequences.
	(sensitive): Add BLAST error correction.
	(sensitive): Divide k-mers in two, and space k-mers. So that each k-mer has four representations, allowing at least one error.

*/

/* primer trim */
/*
SE, both ends have primers
PE, start have primer, end if they overlap

Direction is known, but is rc at the 3'.

Must be given in advance, where graph can be terminated on ends (as with 
adaptors).
Trimming them is the same as for adaptors, but with rc too.

*/

/* merge PE */
/*
need:
	k: k-mer size
	len: minimum overlap.
	edit: max edit dist, expect low distance.

Sample k-mers (only unique?) on forward read according to the pigeon hole principle.
Seed-extend in reverse read, and check with max edit distance.

Adaptor trim if insert is less than merged size.

*/


int print_QCstat(QCstat *src, int minQ, int minPhred, int minmaskQ, int minlen, int maxlen, int fiveClip, int threeClip, FILE *dest) {
	
	int i, n, end, n50, scale;
	unsigned *dist;
	long unsigned tot;
	double p;
	
	n = fprintf(dest, "{\n");
	/* trim parameters */
	n += fprintf(dest, "\t\"Maximum Trim length\": %d,\n", maxlen);
	n += fprintf(dest, "\t\"Minimum Trim length\": %d,\n", minlen);
	n += fprintf(dest, "\t\"5\'-clip\": %d,\n", fiveClip);
	n += fprintf(dest, "\t\"3\'-clip\": %d,\n", threeClip);
	if(src->Eeq) {
		n += fprintf(dest, "\t\"Minimum Q\": %d,\n", minQ);
		n += fprintf(dest, "\t\"End Trim Q\": %d,\n", minPhred);
		n += fprintf(dest, "\t\"Hard Mask Q\": %d,\n", minmaskQ);
		n += fprintf(dest, "\t\"Phred Scale\": %d,\n", src->phredScale);
	}
	
	/* results */
	n += fprintf(dest, "\t\"Fragment Count\": %lu,\n", src->fragcount);
	n += fprintf(dest, "\t\"Org. Fragment Count\": %lu,\n", src->org_fragcount);
	n += fprintf(dest, "\t\"Sequence Count\": %lu,\n", src->count);
	n += fprintf(dest, "\t\"Org. Sequence Count\": %lu,\n", src->org_count);
	n += fprintf(dest, "\t\"Bp Count\": %lu,\n", src->bpcount);
	n += fprintf(dest, "\t\"Org. Bp Count\": %lu,\n", src->org_bpcount);
	
	n += fprintf(dest, "\t\"Mean Read Length\": %f,\n", src->count ? ((double)(src->bpcount) / src->count) : 0);
	n += fprintf(dest, "\t\"Org. Mean Read Length\": %f,\n", src->org_count ? ((double)(src->org_bpcount) / src->org_count) : 0);
	n += fprintf(dest, "\t\"GC Content\": %f,\n", (src->bpcount - src->totns) ? ((double)(src->totgc) / (src->bpcount - src->totns)) : 0);
	n += fprintf(dest, "\t\"Max Sequence Length\": %d,\n", src->maxlen);
	
	/* get N50 */
	dist = src->ldist;
	scale = 1 << src->qresolution;
	if((src->maxlen << 1) < src->bpcount) {
		n50 = 0;
		p = 0;
		tot = 0;
		if(src->qresolution) {
			for(i = 0; i < 511; ++i) {
				if(dist[i]) {
					p = (double)(dist[i+1]) / (dist[i] + dist[i+1]);
					tot += (n50 + p * scale) * dist[i];
					if(src->bpcount < (tot << 1)) {
						n50 += p * scale;
						i = 512;
					} else {
						n50 += scale;
					}
				} else {
					n50 += scale;
				}
			}
		} else {
			end = src->verbose ? (src->maxlen + 1) : 512;
			for(i = 0; i < end; ++i) {
				tot += i * dist[i];
				if(src->bpcount < (tot << 1)) {
					n50 = i;
					i = end;
				}
			}
		}
	} else {
		n50 = src->maxlen;
	}
	n += fprintf(dest, "\t\"N50\": %d,\n", n50);
	
	/* print Q distribution */
	if(src->Eeq) {
		dist = src->qdist;
		n += fprintf(dest, "\t\"E(Q)\": %f,\n", -10 * log10(src->Eeq / src->bpcount));
		n += fprintf(dest, "\t\"Q Distribution\": [%d, %d, %d, %d", dist[0], dist[1], dist[2], dist[3]);
		for(i = 4; i < 256; i += 4) {
			n += fprintf(dest, ", %d, %d, %d, %d", dist[i], dist[i + 1], dist[i + 2], dist[i + 3]);
		}
		n += fprintf(dest, "],\n");
	}
	
	/* print length distribution */
	dist = src->ldist;
	n += fprintf(dest, "\t\"Length Resolution\": %d,\n", scale);
	n += fprintf(dest, "\t\"Length Distribution\": [%d, %d, %d, %d", dist[0], dist[1], dist[2], dist[3]);
	end = src->verbose ? (src->maxlen + 1) : 512;
	for(i = 4; i < end; i += 4) {
		n += fprintf(dest, ", %d, %d, %d, %d", dist[i], dist[i + 1], dist[i + 2], dist[i + 3]);
	}
	n += fprintf(dest, "]\n");
	
	/* end */
	n += fprintf(dest, "}\n");
	
	return n;
}
