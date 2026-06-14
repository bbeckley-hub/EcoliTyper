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
#include <stdlib.h>
#include "compdna.h"
#include "decon.h"
#include "filebuff.h"
#include "hashmapkma.h"
#include "pherror.h"
#include "seqparse.h"
#include "stdnuc.h"
#include "qseqs.h"
#include "updateindex.h"

int (*deConNode_ptr)(CompDNA *, HashMapKMA *, unsigned **);
int (*addCont)(HashMapKMA *, long unsigned, int, unsigned **);

int hashMap_addCont(HashMapKMA *dest, long unsigned key, int value, unsigned **Values) {
	
	unsigned pos, kpos, *values;
	long unsigned kmer;
	
	kpos = key & dest->size;
	pos = getExistPtr(dest->exist, kpos);
	
	if(pos != dest->null_index) {
		kmer = getKeyPtr(dest->key_index, pos);
		while(key != kmer) {
			++pos;
			if(kpos != (kmer & dest->size)) {
				return 0;
			}
			kmer = getKeyPtr(dest->key_index, pos);
		}
		values = updateValuePtr(Values[getValueIndexPtr(dest->value_index, pos)], value);
		if(values) {
			Values[getValueIndexPtr(dest->value_index, pos)] = values;
			return 1;
		}
	}
	
	return 0;
}

int megaMap_addCont(HashMapKMA *dest, long unsigned index, int value, unsigned **Values) {
	
	long unsigned pos;
	unsigned *values;
	
	if((pos = getExistPtr(dest->exist, index)) != dest->n) {
		values = updateValuePtr(Values[pos], value);
		if(values) {
			Values[pos] = values;
			return 1;
		}
	}
	
	return 0;
}

int deConNode(CompDNA *qseq, HashMapKMA *finalDB, unsigned **Values) {
	
	int i, j, end, mapped_cont, shifter, DB_size, mPos, hLen;
	unsigned kmersize, mlen, flag, cPos, iPos;
	long unsigned mask, mmask, kmer, cmer, hmer, *seq;
	
	if(qseq->seqlen < finalDB->kmersize) {
		return 0;
	}
	
	/* set parameters */
	DB_size = finalDB->DB_size;
	mapped_cont = 0;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (finalDB->kmersize << 1);
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	seq = qseq->seq;
	kmersize = finalDB->kmersize;
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = finalDB->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = finalDB->flag;
	
	j = 0;
	for(i = 1; i <= qseq->N[0]; ++i) {
		/* init k-mer */
		getKmer_macro(kmer, seq, j, cPos, iPos, (shifter + 2));
		cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
		end = qseq->N[i];
		for(j += kmersize - 1;j < end; ++j) {
			/* update k-mer */
			kmer = updateKmer_macro(kmer, seq, j, mask);
			cmer = flag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
			
			/* update hashMap */
			mapped_cont += addCont(finalDB, cmer, DB_size, Values);
		}
		j = end + 1;
	}
	qseq->N[0]--;
	
	return mapped_cont;
}

int deConNode_sparse(CompDNA *qseq, HashMapKMA *finalDB, unsigned **Values) {
	
	int i, j, end, mapped_cont, DB_size, shifter, prefix_len, prefix_shifter;
	int mlen, flag, mPos, hLen, cPos, iPos;
	long unsigned mmask, prefix, kmer, cmer;
	
	if(qseq->seqlen < finalDB->kmersize) {
		return 0;
	}
	DB_size = finalDB->DB_size;
	prefix = finalDB->prefix;
	prefix_len = finalDB->prefix_len;
	prefix_shifter = sizeof(long unsigned) * sizeof(long unsigned) - (prefix_len << 1);
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (finalDB->kmersize << 1);
	mapped_cont = 0;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	mlen = finalDB->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = finalDB->flag;
	j = 0;
	for(i = 1; i <= qseq->N[0]; ++i) {
		end = qseq->N[i] - prefix_len - finalDB->kmersize + 1;
		for(;j < end; ++j) {
			if(prefix_len == 0 || getKmer(qseq->seq, j, prefix_shifter) == prefix) {
				/* get kmer */
				getKmer_macro(kmer, qseq->seq, (j + prefix_len), cPos, iPos, shifter);
				cmer = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
				
				/* update hashMap */
				mapped_cont += addCont(finalDB, cmer, DB_size, Values);
			}
		}
		j = qseq->N[i] + 1;
	}
	qseq->N[0]--;
	return mapped_cont;
}

unsigned deConDB(HashMapKMA *finalDB, char **inputfiles, int fileCount, char *trans, unsigned **Values) {
	
	int FASTQ;
	unsigned fileCounter, mapped_cont;
	char *filename;
	Qseqs *header, *qseq;
	FileBuff *inputfile;
	CompDNA *compressor;
	
	/* allocate */
	compressor = smalloc(sizeof(CompDNA));
	allocComp(compressor, 1024);
	header = setQseqs(1024);
	qseq = setQseqs(1024);
	inputfile = setFileBuff(1024 * 1024);
	
	/* set variables */
	mapped_cont = 0;
	if(--finalDB->size == finalDB->mask) {
		addCont = &megaMap_addCont;
	}
	
	/* iterate inputfiles */
	for(fileCounter = 0; fileCounter < fileCount; ++fileCounter) {
		/* open file */
		filename = inputfiles[fileCounter];
		/* determine filetype and open it */
		if((FASTQ = openAndDetermine(inputfile, filename)) & 3) {
			fprintf(stderr, "%s\t%s\n", "# Reading inputfile: ", filename);
			/* parse the file */
			while(FileBuffgetFsa(inputfile, header, qseq, trans)) {
				fprintf(stderr, "# Decon:\t%s\n", header->seq + 1);
				if(qseq->len > finalDB->kmersize) {
					/* compress DNA */
					if(qseq->len >= compressor->size) {
						freeComp(compressor);
						allocComp(compressor, qseq->len);
					}
					compDNAref(compressor, qseq->seq, qseq->len);
					
					/* Add contamination */
					mapped_cont += deConNode_ptr(compressor, finalDB, Values);
					/* rc */
					comp_rc(compressor);
					mapped_cont += deConNode_ptr(compressor, finalDB, Values);
				}
			}
			
			/* close file buffer */
			if(FASTQ & 4) {
				gzcloseFileBuff(inputfile);
			} else {
				closeFileBuff(inputfile);
			}
		}
	}
	++finalDB->size;
	
	/* clean */
	freeComp(compressor);
	free(compressor);
	destroyQseqs(header);
	destroyQseqs(qseq);
	destroyFileBuff(inputfile);
	
	return mapped_cont;
}
