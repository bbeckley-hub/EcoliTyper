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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "compdna.h"
#include "filebuff.h"
#include "pherror.h"
#include "runinput.h"
#include "qc.h"
#include "qseqs.h"
#include "seqparse.h"
#include "stdstat.h"

void (*printFsa_ptr)(Qseqs*, Qseqs*, Qseqs*, CompDNA*, FILE*) = &printFsa;
void (*printFsa_pair_ptr)(Qseqs*, Qseqs*, Qseqs*, Qseqs*, Qseqs*, Qseqs*, CompDNA*, FILE*) = &printFsa_pair;

unsigned hardmask(unsigned char *seq, unsigned char *qual, int len, const int phredScale, int minQ) {
	
	unsigned n;
	
	if(!len || !minQ) {
		return 0;
	}
	
	minQ += phredScale;
	n = 0;
	--qual;
	do {
		if(*++qual < minQ) {
			*seq = 4;
			++n;
		}
		++seq;
	} while(--len);
	
	return n;
}

int gcontent(unsigned char *seq, int len) {
	
	int gc;
	
	gc = 0;
	++len;
	while(--len) {
		if(*seq == 1 || *seq == 2) {
			++gc;
		}
		++seq;
	}
	
	return gc;
}

double qcstat(const unsigned char *seq, const unsigned char *qual, const int len, const double *prob, const int hardmaskQ, int *GC, int *NS, int *EQ) {
	
	/*
	static const double prob[128] = {
		1.0000000000000000, 0.7943282347242815, 0.6309573444801932, 0.5011872336272722, 0.3981071705534972, 0.3162277660168379, 0.2511886431509580, 0.1995262314968880,
		0.1584893192461113, 0.1258925411794167, 0.1000000000000000, 0.0794328234724281, 0.0630957344480193, 0.0501187233627272, 0.0398107170553497, 0.0316227766016838,
		0.0251188643150958, 0.0199526231496888, 0.0158489319246111, 0.0125892541179417, 0.0100000000000000, 0.0079432823472428, 0.0063095734448019, 0.0050118723362727,
		0.0039810717055350, 0.0031622776601684, 0.0025118864315096, 0.0019952623149689, 0.0015848931924611, 0.0012589254117942, 0.0010000000000000, 0.0007943282347243,
		0.0006309573444802, 0.0005011872336273, 0.0003981071705535, 0.0003162277660168, 0.0002511886431510, 0.0001995262314969, 0.0001584893192461, 0.0001258925411794,
		0.0001000000000000, 0.0000794328234724, 0.0000630957344480, 0.0000501187233627, 0.0000398107170553, 0.0000316227766017, 0.0000251188643151, 0.0000199526231497,
		0.0000158489319246, 0.0000125892541179, 0.0000100000000000, 0.0000079432823472, 0.0000063095734448, 0.0000050118723363, 0.0000039810717055, 0.0000031622776602,
		0.0000025118864315, 0.0000019952623150, 0.0000015848931925, 0.0000012589254118, 0.0000010000000000, 0.0000007943282347, 0.0000006309573445, 0.0000005011872336,
		0.0000003981071706, 0.0000003162277660, 0.0000002511886432, 0.0000001995262315, 0.0000001584893192, 0.0000001258925412, 0.0000001000000000, 0.0000000794328235,
		0.0000000630957344, 0.0000000501187234, 0.0000000398107171, 0.0000000316227766, 0.0000000251188643, 0.0000000199526231, 0.0000000158489319, 0.0000000125892541,
		0.0000000100000000, 0.0000000079432823, 0.0000000063095734, 0.0000000050118723, 0.0000000039810717, 0.0000000031622777, 0.0000000025118864, 0.0000000019952623,
		0.0000000015848932, 0.0000000012589254, 0.0000000010000000, 0.0000000007943282, 0.0000000006309573, 0.0000000005011872, 0.0000000003981072, 0.0000000003162278,
		0.0000000002511886, 0.0000000001995262, 0.0000000001584893, 0.0000000001258925, 0.0000000001000000, 0.0000000000794328, 0.0000000000630957, 0.0000000000501187,
		0.0000000000398107, 0.0000000000316228, 0.0000000000251189, 0.0000000000199526, 0.0000000000158489, 0.0000000000125893, 0.0000000000100000, 0.0000000000079433,
		0.0000000000063096, 0.0000000000050119, 0.0000000000039811, 0.0000000000031623, 0.0000000000025119, 0.0000000000019953, 0.0000000000015849, 0.0000000000012589,
		0.0000000000010000, 0.0000000000007943, 0.0000000000006310, 0.0000000000005012, 0.0000000000003981, 0.0000000000003162, 0.0000000000002512, 0.0000000000001995};
	*/
	/*
	E(Q) = -10 * log_10(sum(10^(-Q/10)) / |Q|) 
	*/
	unsigned i, ns, gc;
	unsigned char *sptr, *qptr;
	double eq;
	
	/* init */
	ns = 0;
	gc = 0;
	eq = 0;
	
	/* sum quality scores */
	i = len + 1;
	sptr = (unsigned char*)(seq) - 1;
	qptr = (unsigned char*)(qual) - 1;
	while(--i) {
		eq += prob[*++qptr];
		if(*++sptr == 4 || *qptr < hardmaskQ) {
			*sptr = 4;
			++ns;
		} else if(*sptr == 1 || *sptr == 2) {
			++gc;
		}
	}
	
	*GC = gc;
	*NS = ns;
	*EQ = ceil(-10 * log10(eq / len));
	
	/* return average Q */
	return eq;
}

int phredStat(unsigned char *seq, unsigned char *qual, int len, const double *prob, const int minPhred, const int minQ, const int hardmaskQ, const int fiveClip, const int threeClip, const int minlen, const int maxlen, int *START, int *END, QCstat *qcreport) {
	
	unsigned i, ns, gc, ns5, gc5, ns3, gc3, l5, l3, start, end;
	unsigned char *sptr5, *qptr5, *sptr3, *qptr3;
	double minP, sp, sp5, sp3;
	
	if(qcreport) {
		qcreport->org_count++;
		qcreport->org_bpcount += len;
	}
	/* too long */
	if(maxlen < len) {
		*START = 0;
		*END = 0;
		return 0;
	}
	
	/* init */
	qptr5 = qual;
	qptr3 = qptr5 + len - 1;
	start = 0;
	end = len;
	
	/* end trim */
	while(start < end && *qptr5 < minPhred) {
		++start;
		++qptr5;
	}
	sptr5 = seq + start;
	while(start < end && *qptr3 < minPhred) {
		--end;
		--qptr3;
	}
	sptr3 = seq + end - 1;
	len = end - start;
	
	if(!minQ && !hardmaskQ && !qcreport) {
		/* set return values */
		*START = start;
		*END = end;
		return len;
	}
	
	/* minQ trim */
	ns = 0;
	gc = 0;
	sp = 0;
	i = len + 1;
	--qptr5;
	--sptr5;
	while(--i) {
		sp += prob[*++qptr5];
		if(*++sptr5 == 4 || *qptr5 < hardmaskQ) {
			*sptr5 = 4;
			++ns;
		} else if(*sptr5 == 1 || *sptr5 == 2) {
			++gc;
		}
	}
	qptr5 = qual + start;
	sptr5 = seq + start;
	
	/* bi-directional phred trim */
	minP = pow(10, (-0.1) * minQ);
	/* set 3' segment, 5' is set when running the algorithm */
	if(minlen <= (len - ns) && (minP * len) < sp) {
		ns5 = 0;
		ns3 = 0;
		gc5 = 0;
		gc3 = 0;
		l5 = 0;
		l3 = 0;
		sp5 = 0;
		sp3 = 0;
		/* high qual part */
		while(l3 < len && minPhred <= *qptr3) {
			sp3 += prob[*qptr3--];
			++l3;
			if(*sptr3 == 1 || *sptr3 == 2) {
				++gc3;
			} else if(*sptr3 == 4) {
				++ns3;
			}
			--sptr3;
			
		}
		/* low qual part */
		while(l3 < len && *qptr3 < minPhred) {
			sp3 += prob[*qptr3--];
			++l3;
			if(*sptr3 == 1 || *sptr3 == 2) {
				++gc3;
			} else if(*sptr3 == 4) {
				++ns3;
			}
			--sptr3;
			
		}
		/* remove ends */
		while(minlen <= (len - ns) && (minP * len) < sp) {
			/* trim lowest quality end */
			if((sp5 * l3) < (sp3 * l5)) { /* (sp5 / l5) < (sp3 / l3) */
				/* trim 3' */
				end -= l3;
				ns -= ns3;
				gc -= gc3;
				len -= l3;
				sp -= sp3;
				ns3 = 0;
				gc3 = 0;
				l3 = 0;
				sp3 = 0;
				/* get next 3' segment */
				/* high qual part */
				while(l3 < len && minPhred <= *qptr3) {
					sp3 += prob[*qptr3--];
					++l3;
					if(*sptr3 == 1 || *sptr3 == 2) {
						++gc3;
					} else if(*sptr3 == 4) {
						++ns3;
					}
					--sptr3;
					
				}
				/* low qual part */
				while(l3 < len && *qptr3 < minPhred) {
					sp3 += prob[*qptr3--];
					++l3;
					if(*sptr3 == 1 || *sptr3 == 2) {
						++gc3;
					} else if(*sptr3 == 4) {
						++ns3;
					}
					--sptr3;
					
				}
			} else {
				/* trim 5' */
				start += l5;
				len -= l5;
				ns -= ns5;
				gc -= gc5;
				sp -= sp5;
				ns5 = 0;
				gc5 = 0;
				l5 = 0;
				sp5 = 0;
				/* get next 5' segment */
				/* high qual part */
				while(l5 < len && minPhred <= *qptr5) {
					sp5 += prob[*qptr5++];
					++l5;
					if(*sptr5 == 1 || *sptr5 == 2) {
						++gc5;
					} else if(*sptr5 == 4) {
						++ns5;
					}
					++sptr5;
					
				}
				/* low qual part */
				while(l5 < len && *qptr5 < minPhred) {
					sp5 += prob[*qptr5++];
					++l5;
					if(*sptr5 == 1 || *sptr5 == 2) {
						++gc5;
					} else if(*sptr5 == 4) {
						++ns5;
					}
					++sptr5;
				}
			}
		}
	}
	
	/* update QC report */
	if(qcreport && minlen <= (len - ns)) {
		update_QCstat(qcreport, len, gc, ns, sp);
	}
	
	/* set return values */
	*START = start;
	*END = end;
	
	return len - ns;
}

int fsastat(unsigned char *seq, int len, const int minlen, const int maxlen, int *START, int *END, QCstat *qcreport) {
	
	unsigned i, start, end, ns, gc;
	unsigned char *sptr;
	
	if(qcreport) {
		qcreport->org_count++;
		qcreport->org_bpcount += len;
	}
	if(maxlen < len) {
		*START = 0;
		*END = 0;
		return 0;
	}
	
	/* init */
	start = 0;
	end = len;
	
	/* end trim */
	sptr = seq + end - 1;
	while(start <= end && *sptr == 4) {
		--end;
		--sptr;
	}
	sptr = seq;
	while(start < end && *sptr == 4) {
		++start;
		++sptr;
	}
	len = end - start;
	
	/* gc- and n- content */
	ns = 0;
	gc = (*sptr == 1 || *sptr == 2) ? 1 : 0;
	i = len;
	while(--i) {
		if(*++sptr == 4) {
			++ns;
		} else if(*sptr == 1 || *sptr == 2) {
			++gc;
		}
	}
	
	if(qcreport && minlen <= len) {
		update_QCstat(qcreport, len, gc, ns, 0);
	}
	
	/* set return values */
	*START = start;
	*END = end;
	
	return len - ns;
}

long unsigned run_input(char **inputfiles, int fileCount, int minPhred, int hardmaskQ, int minQ, int fiveClip, int threeClip, int minlen, int maxlen, char *trans, const double *prob, QCstat *qcreport, FILE *out) {
	
	int fileCounter, phredScale, start, end, len;
	unsigned FASTQ;
	long unsigned count, org_count;
	char *filename;
	Qseqs *header, *qseq, *qual;
	FileBuff *inputfile;
	CompDNA *compressor;
	
	if(minPhred < minQ) {
		minPhred = minQ;
	}
	compressor = smalloc(sizeof(CompDNA));
	allocComp(compressor, 1024);
	header = setQseqs(256);
	qseq = setQseqs(1024);
	qual = setQseqs(1024);
	inputfile = setFileBuff(CHUNK);
	phredScale = 33;
	count = 0;
	org_count = 0;
	
	for(fileCounter = 0; fileCounter < fileCount; ++fileCounter) {
		filename = (char*)(inputfiles[fileCounter]);
		
		/* determine filetype and open it */
		if((FASTQ = openAndDetermine(inputfile, filename)) & 3) {
			fprintf(stderr, "%s\t%s\n", "# Reading inputfile: ", filename);
		}
		
		/* parse the file */
		if(FASTQ & 1) {
			/* get phred scale */
			phredScale = getPhredFileBuff(inputfile);
			fprintf(stderr, "# Phred scale:\t%d\n", phredScale);
			
			/* parse reads */
			while(FileBuffgetFq(inputfile, header, qseq, qual, trans)) {
				++org_count;
				len = phredStat(qseq->seq, qual->seq, qseq->len, prob - phredScale, phredScale + minPhred, minQ, hardmaskQ, fiveClip, threeClip, minlen, maxlen, &start, &end, qcreport);
				
				/* print */
				if(minlen <= len) {
					/* dump seq */
					qseq->seq += start;
					qual->seq += start;
					qseq->len = end - start;
					qual->len = end - start;
					printFsa_ptr(header, qseq, qual, compressor, out);
					qseq->seq -= start;
					qual->seq -= start;
					++count;
				}
			}
		} else if(FASTQ & 2) {
			while(FileBuffgetFsa(inputfile, header, qseq, trans)) {
				++org_count;
				len = fsastat(qseq->seq, qseq->len, minlen, maxlen, &start, &end, qcreport);
				
				if(minlen <= len) {
					/* dump seq */
					qseq->seq += start;
					qseq->len = end - start;
					printFsa_ptr(header, qseq, 0, compressor, out);
					qseq->seq -= start;
					++count;
				}
			}
		}
		
		if(FASTQ & 4) {
			gzcloseFileBuff(inputfile);
		} else {
			closeFileBuff(inputfile);
		}
	}
	
	if(qcreport) {
		qcreport->fragcount += count;
		qcreport->org_fragcount += org_count;
		qcreport->phredScale = phredScale;
	}
	freeComp(compressor);
	free(compressor);
	destroyQseqs(header);
	destroyQseqs(qseq);
	destroyQseqs(qual);
	destroyFileBuff(inputfile);
	
	return count;
}

long unsigned run_input_PE(char **inputfiles, int fileCount, int minPhred, int hardmaskQ, int minQ, int fiveClip, int threeClip, int minlen, int maxlen, char *trans, const double *prob, QCstat *qcreport, FILE *out) {
	
	int fileCounter, phredScale, start, start2, end, end2, len, len2;
	unsigned FASTQ, FASTQ2;
	long unsigned count, org_count;
	char *filename;
	Qseqs *header, *qseq, *qual, *header2, *qseq2, *qual2;
	FileBuff *inputfile, *inputfile2;
	CompDNA *compressor;
	
	if(minPhred < minQ) {
		minPhred = minQ;
	}
	compressor = smalloc(sizeof(CompDNA));
	allocComp(compressor, 1024);
	header = setQseqs(256);
	qseq = setQseqs(1024);
	qual = setQseqs(1024);
	header2 = setQseqs(256);
	qseq2 = setQseqs(1024);
	qual2 = setQseqs(1024);
	inputfile = setFileBuff(CHUNK);
	inputfile2 = setFileBuff(CHUNK);
	phredScale = 33;
	count = 0;
	org_count = 0;
	
	for(fileCounter = 0; fileCounter < fileCount; ++fileCounter) {
		filename = inputfiles[fileCounter];
		/* determine filetype and open it */
		FASTQ = openAndDetermine(inputfile, filename);
		++fileCounter;
		filename = inputfiles[fileCounter];
		FASTQ2 = openAndDetermine(inputfile2, filename);
		if(FASTQ == FASTQ2) {
			fprintf(stderr, "# Reading inputfile:\t%s %s\n", inputfiles[fileCounter-1], filename);
		} else {
			fprintf(stderr, "Inputfiles:\t%s %s\nAre in different format.\n", inputfiles[fileCounter-1], filename);
			FASTQ = 0;
			errno = 1;
		}
		
		/* parse the file */
		if(FASTQ & 1) {
			/* get phred scale */
			phredScale = getPhredFileBuff(inputfile);
			if(phredScale == 0) {
				phredScale = getPhredFileBuff(inputfile2);
			}
			fprintf(stderr, "# Phred scale:\t%d\n", phredScale);
			
			/* parse reads */
			//while(FileBuffgetFq(inputfile, header, qseq, qual) && FileBuffgetFq(inputfile2, header2, qseq2, qual2)) {
			/* this ensures reading of truncated files */
			while((FileBuffgetFq(inputfile, header, qseq, qual, trans) | FileBuffgetFq(inputfile2, header2, qseq2, qual2, trans))) {
				++org_count;
				len = phredStat(qseq->seq, qual->seq, qseq->len, prob - phredScale, phredScale + minPhred, minQ, hardmaskQ, fiveClip, threeClip, minlen, maxlen, &start, &end, qcreport);
				len2 = phredStat(qseq2->seq, qual2->seq, qseq2->len, prob - phredScale, phredScale + minPhred, minQ, hardmaskQ, fiveClip, threeClip, minlen, maxlen, &start2, &end2, qcreport);
				
				/* print */
				qseq->seq += start;
				qual->seq += start;
				qseq->len = end - start;
				qual->len = end - start;
				qseq2->seq += start2;
				qual2->seq += start2;
				qseq2->len = end2 - start2;
				qual2->len = end2 - start2;
				if(minlen <= len && minlen <= len2) {
					printFsa_pair_ptr(header, qseq, qual, header2, qseq2, qual2, compressor, out);
					++count;
				} else if(minlen <= len) {
					printFsa_ptr(header, qseq, qual, compressor, out);
					++count;
				} else if(minlen <= len2) {
					printFsa_ptr(header2, qseq2, qual2, compressor, out);
					++count;
				}
				qseq->seq -= start;
				qual->seq -= start;
				qseq2->seq -= start2;
				qual2->seq -= start2;
			}
		} else if(FASTQ & 2) {
			while((FileBuffgetFsa(inputfile, header, qseq, trans) | FileBuffgetFsa(inputfile2, header2, qseq2, trans))) {
				++org_count;
				len = fsastat(qseq->seq, qseq->len, minlen, maxlen, &start, &end, qcreport);
				len2 = fsastat(qseq2->seq, qseq2->len, minlen, maxlen, &start2, &end2, qcreport);
				
				/* print */
				qseq->seq += start;
				qseq->len = end - start;
				qseq2->seq += start2;
				qseq2->len = end - start;
				if(minlen <= len && minlen <= len2) {
					printFsa_pair_ptr(header, qseq, 0, header2, qseq2, 0, compressor, out);
					++count;
				} else if(minlen <= len) {
					printFsa_ptr(header, qseq, 0, compressor, out);
					++count;
				} else if(minlen <= len2) {
					printFsa_ptr(header2, qseq2, 0, compressor, out);
					++count;
				}
				qseq->seq -= start;
				qseq2->seq -= start2;
			}
		}
		
		--fileCounter;
		if(FASTQ & 4) {
			gzcloseFileBuff(inputfile);
		} else {
			closeFileBuff(inputfile);
		}
		++fileCounter;
		if(FASTQ2 & 4) {
			gzcloseFileBuff(inputfile2);
		} else {
			closeFileBuff(inputfile2);
		}
		
	}
	
	if(qcreport) {
		qcreport->fragcount += count;
		qcreport->org_fragcount += org_count;
		if(qcreport->Eeq) {
			qcreport->phredScale = phredScale;
		}
	}
	freeComp(compressor);
	free(compressor);
	destroyQseqs(header);
	destroyQseqs(qseq);
	destroyQseqs(qual);
	destroyQseqs(header2);
	destroyQseqs(qseq2);
	destroyQseqs(qual2);
	destroyFileBuff(inputfile);
	destroyFileBuff(inputfile2);
	
	return count;
}

long unsigned run_input_INT(char **inputfiles, int fileCount, int minPhred, int hardmaskQ, int minQ, int fiveClip, int threeClip, int minlen, int maxlen, char *trans, const double *prob, QCstat *qcreport, FILE *out) {
	
	int fileCounter, phredScale, start, start2, end, end2, len, len2;
	unsigned FASTQ;
	long unsigned count, org_count;
	char *filename;
	Qseqs *header, *qseq, *qual, *header2, *qseq2, *qual2;
	FileBuff *inputfile;
	CompDNA *compressor;
	
	if(minPhred < minQ) {
		minPhred = minQ;
	}
	compressor = smalloc(sizeof(CompDNA));
	allocComp(compressor, 1024);
	header = setQseqs(256);
	qseq = setQseqs(1024);
	qual = setQseqs(1024);
	header2 = setQseqs(256);
	qseq2 = setQseqs(1024);
	qual2 = setQseqs(1024);
	inputfile = setFileBuff(CHUNK);
	phredScale = 33;
	count = 0;
	org_count = 0;
	
	for(fileCounter = 0; fileCounter < fileCount; ++fileCounter) {
		filename = (char*)(inputfiles[fileCounter]);
		
		/* determine filetype and open it */
		if((FASTQ = openAndDetermine(inputfile, filename)) & 3) {
			fprintf(stderr, "%s\t%s\n", "# Reading inputfile: ", filename);
		}
		
		/* parse the file */
		if(FASTQ & 1) {
			/* get phred scale */
			phredScale = getPhredFileBuff(inputfile);
			fprintf(stderr, "# Phred scale:\t%d\n", phredScale);
			
			/* parse reads */
			while((FileBuffgetFq(inputfile, header, qseq, qual, trans) | FileBuffgetFq(inputfile, header2, qseq2, qual2, trans))) {
				org_count++;
				len = phredStat(qseq->seq, qual->seq, qseq->len, prob - phredScale, phredScale + minPhred, minQ, hardmaskQ, fiveClip, threeClip, minlen, maxlen, &start, &end, qcreport);
				len2 = phredStat(qseq2->seq, qual2->seq, qseq2->len, prob - phredScale, phredScale + minPhred, minQ, hardmaskQ, fiveClip, threeClip, minlen, maxlen, &start2, &end2, qcreport);
				
				/* print */
				qseq->seq += start;
				qual->seq += start;
				qseq->len = end - start;
				qual->len = end - start;
				qseq2->seq += start2;
				qual2->seq += start2;
				qseq2->len = end2 - start2;
				qual2->len = end2 - start2;
				if(minlen <= len && minlen <= len2) {
					printFsa_pair_ptr(header, qseq, qual, header2, qseq2, qual2, compressor, out);
					++count;
				} else if(minlen <= len) {
					printFsa_ptr(header, qseq, qual, compressor, out);
					++count;
				} else if(minlen <= len2) {
					printFsa_ptr(header2, qseq2, qual2, compressor, out);
					++count;
				}
				qseq->seq -= start;
				qual->seq -= start;
				qseq2->seq -= start2;
				qual2->seq -= start2;
			}
		} else if(FASTQ & 2) {
			while((FileBuffgetFsa(inputfile, header, qseq, trans) | FileBuffgetFsa(inputfile, header2, qseq2, trans))) {
				++org_count;
				len = fsastat(qseq->seq, qseq->len, minlen, maxlen, &start, &end, qcreport);
				len2 = fsastat(qseq2->seq, qseq2->len, minlen, maxlen, &start2, &end2, qcreport);
				
				/* print */
				qseq->seq += start;
				qseq->len = end - start;
				qseq2->seq += start2;
				qseq2->len = end - start;
				if(minlen <= len && minlen <= len2) {
					printFsa_pair_ptr(header, qseq, 0, header2, qseq2, 0, compressor, out);
					++count;
				} else if(minlen <= len) {
					printFsa_ptr(header, qseq, 0, compressor, out);
					++count;
				} else if(minlen <= len2) {
					printFsa_ptr(header2, qseq2, 0, compressor, out);
					++count;
				}
				qseq->seq -= start;
				qseq2->seq -= start2;
			}
		}
		
		if(FASTQ & 4) {
			gzcloseFileBuff(inputfile);
		} else {
			closeFileBuff(inputfile);
		}
		
	}
	
	if(qcreport) {
		qcreport->fragcount += count;
		qcreport->org_fragcount += org_count;
		if(qcreport->Eeq) {
			qcreport->phredScale = phredScale;
		}
	}
	freeComp(compressor);
	free(compressor);
	destroyQseqs(header);
	destroyQseqs(qseq);
	destroyQseqs(qual);
	destroyQseqs(header2);
	destroyQseqs(qseq2);
	destroyQseqs(qual2);
	destroyFileBuff(inputfile);
	
	return count;
}

void bootFsa(Qseqs *header, Qseqs *qseq, Qseqs *qual, CompDNA *compressor, FILE *out) {
	
	int i, end, buffer[4];
	
	/* bootstrap in pieces of 1024 */
	buffer[3] = header->len;
	end = qseq->len - 1024;
	for(i = 0; i < end; i += 512) {
		compDNA(compressor, qseq->seq + i, 1024);
		/* print */
		buffer[0] = compressor->seqlen;
		buffer[1] = compressor->complen;
		buffer[2] = compressor->N[0];
		
		sfwrite(buffer, sizeof(int), 4, out);
		sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, out);
		sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], out);
		sfwrite((header->seq + 1), 1, header->len, out);
		resetComp(compressor);
	}
	
	compDNA(compressor, qseq->seq + i, qseq->len - i);
	/* print */
	buffer[0] = compressor->seqlen;
	buffer[1] = compressor->complen;
	buffer[2] = compressor->N[0];
	sfwrite(buffer, sizeof(int), 4, out);
	sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, out);
	sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], out);
	sfwrite((header->seq + 1), 1, header->len, out);
	resetComp(compressor);
}

void printFsa(Qseqs *header, Qseqs *qseq, Qseqs *qual, CompDNA *compressor, FILE *out) {
	
	int buffer[4];
	
	/* translate to 2bit */
	if(qseq->len >= compressor->size) {
		freeComp(compressor);
		allocComp(compressor, qseq->len);
	}
	compDNA(compressor, qseq->seq, qseq->len);
	
	buffer[0] = compressor->seqlen;
	buffer[1] = compressor->complen;
	buffer[2] = compressor->N[0];
	buffer[3] = header->len;
	
	sfwrite(buffer, sizeof(int), 4, out);
	sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, out);
	sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], out);
	sfwrite((header->seq + 1), 1, header->len, out);
	
	resetComp(compressor);
}

void printFsa_pair(Qseqs *header, Qseqs *qseq, Qseqs *qual, Qseqs *header_r, Qseqs *qseq_r, Qseqs *qual_r, CompDNA *compressor, FILE *out) {
	
	int buffer[4];
	
	/* translate to 2bit */
	if(qseq->len >= compressor->size) {
		freeComp(compressor);
		allocComp(compressor, qseq->len);
	}
	compDNA(compressor, qseq->seq, qseq->len);
	
	buffer[0] = compressor->seqlen;
	buffer[1] = compressor->complen;
	buffer[2] = compressor->N[0];
	buffer[3] = -header->len;
	
	sfwrite(buffer, sizeof(int), 4, out);
	sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, out);
	sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], out);
	sfwrite((header->seq + 1), 1, header->len, out);
	resetComp(compressor);
	
	/* translate to 2bit */
	if(qseq_r->len >= compressor->size) {
		freeComp(compressor);
		allocComp(compressor, qseq_r->len);
	}
	compDNA(compressor, qseq_r->seq, qseq_r->len);
	
	buffer[0] = compressor->seqlen;
	buffer[1] = compressor->complen;
	buffer[2] = compressor->N[0];
	buffer[3] = header_r->len;
	
	sfwrite(buffer, sizeof(int), 4, out);
	sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, out);
	sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], out);
	sfwrite((header_r->seq + 1), 1, header_r->len, out);
	resetComp(compressor);
}
