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
#include <fcntl.h>
#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "align.h"
#include "alnfrags.h"
#include "assembly.h"
#include "chain.h"
#include "compdna.h"
#include "conclave.h"
#include "ef.h"
#include "filebuff.h"
#include "frags.h"
#include "hashmapcci.h"
#include "kmapipe.h"
#include "nw.h"
#include "pherror.h"
#include "printconsensus.h"
#include "qseqs.h"
#include "runkma.h"
#include "sam.h"
#include "spltdb.h"
#include "stdnuc.h"
#include "stdstat.h"
#include "tmp.h"
#include "tsv.h"
#include "updatescores.h"
#include "version.h"
#include "vcf.h"
#include "xml.h"

int print_ankers_spltDB(int *out_Tem, CompDNA *qseq, int rc_flag, const Qseqs *header, const int flag, FILE *out) {
	
	static unsigned target = 1, allIn = 0, thread_num = 1;
	static SpltDBbuff **buffers = 0;
	int num, size, tNum, infoSize[8];
	unsigned char *buff;
	SpltDBbuff *node, *nodeN, *node0;
	
	if(buffers == 0) {
		thread_num = rc_flag;
		buffers = calloc(thread_num, sizeof(SpltDBbuff *));
		if(!buffers) {
			ERROR();
		}
		return 0;
	}
	if(qseq == 0) {
		/* catch remains after last call */
		while(allIn) {
			/* find and write next target */
			node = *buffers;
			num = thread_num;
			tNum = 0;
			while(node == 0 && --num) {
				node = buffers[(tNum = num)];
			}
			node = node->next;
			target = node->num;
			while(--num) {
				if(buffers[num] && buffers[num]->next->num < target) {
					tNum = num;
					node = buffers[num]->next;
					target = node->num;
				}
			}
			sfwrite(node->buff, 1, node->size, out);
			
			/* update cylinder */
			nodeN = node->prev;
			node0 = node->next;
			if(nodeN == node) {
				--allIn;
				nodeN = 0;
			} else {
				nodeN->next = node0;
				node0->prev = nodeN;
			}
			free(node->buff);
			free(node);
			buffers[tNum] = nodeN;	
		}
		sfwrite(&(unsigned){UINT_MAX}, sizeof(unsigned), 1, out);
		sfwrite(&(int){out_Tem[2] - 1}, sizeof(int), 1, out);
		return 0;
	}
	
	num = out_Tem[-1];
	infoSize[0] = num;
	infoSize[1] = qseq->seqlen;
	infoSize[2] = qseq->complen;
	infoSize[3] = qseq->N[0];
	infoSize[4] = rc_flag;
	infoSize[5] = *out_Tem;
	infoSize[6] = header->len;
	infoSize[7] = flag;
	
	if(num == target || thread_num == 1) {
		sfwrite(infoSize, sizeof(int), 8, out);
		sfwrite(qseq->seq, sizeof(long unsigned), qseq->complen, out);
		if(qseq->N[0]) {
			sfwrite(qseq->N + 1, sizeof(int), qseq->N[0], out);
		}
		sfwrite(out_Tem + 1, sizeof(int), *out_Tem, out);
		sfwrite(header->seq, 1, header->len, out);
	} else {
		node = smalloc(sizeof(SpltDBbuff));
		size = header->len + sizeof(int) * (8 + qseq->N[0] + out_Tem[0]) + sizeof(long unsigned) * qseq->complen;
		buff = smalloc(size);
		node->num = num;
		node->size = size;
		node->buff = buff;
		memcpy(buff, infoSize, (size = sizeof(int) * 8));
		buff += size;
		memcpy(buff, qseq->seq, (size = sizeof(long unsigned) * qseq->complen));
		buff += size;
		memcpy(buff, qseq->N + 1, (size = sizeof(int) * qseq->N[0]));
		buff += size;
		memcpy(buff, out_Tem + 1, (size = sizeof(int) * (*out_Tem)));
		buff += size;
		memcpy(buff, header->seq, header->len);
		
		/* find cylinder piece */
		tNum = out_Tem[-2];
		/* extend cylinder piece */
		if((nodeN = buffers[tNum])) {
			node0 = nodeN->next;
			node->prev = nodeN;
			nodeN->next = node;
			node->next = node0;
			node0->prev = node;
		} else {
			node->next = node;
			node->prev = node;
			++allIn;
		}
		buffers[tNum] = node;
	}
		
	/* collect remains */
	while(allIn == thread_num) {
		/* find and write next target */
		node = (*buffers)->next;
		target = node->num;
		num = allIn;
		tNum = 0;
		while(--num) {
			if(buffers[num]->next->num < target) {
				tNum = num;
				node = buffers[num]->next;
				target = node->num;
			}
		}
		sfwrite(node->buff, 1, node->size, out);
		
		/* update cylinder */
		nodeN = node->prev;
		node0 = node->next;
		if(nodeN == node) {
			--allIn;
			nodeN = 0;
		} else {
			nodeN->next = node0;
			node0->prev = nodeN;
		}
		free(node->buff);
		free(node);
		buffers[tNum] = nodeN;
	}
		
	return 0;
}

int print_ankers_Sparse_spltDB(int *out_Tem, CompDNA *qseq, int rc_flag, const Qseqs *header, const int flag, FILE *out) {
	
	static unsigned target = 1, allIn = 0, thread_num = 1;;
	static SpltDBbuff **buffers = 0;
	int num, size, tNum, infoSize[8];
	unsigned char *buff;
	SpltDBbuff *node, *nodeN, *node0;
	
	if(buffers == 0) {
		thread_num = rc_flag;
		buffers = calloc(thread_num, sizeof(SpltDBbuff *));
		if(!buffers) {
			ERROR();
		}
		return 0;
	}
	if(qseq == 0) {
		/* catch remains after last call */
		while(allIn) {
			/* find and write next target */
			node = *buffers;
			num = thread_num;
			tNum = 0;
			while(node == 0 && --num) {
				node = buffers[(tNum = num)];
			}
			node = node->next;
			target = node->num;
			while(--num) {
				if(buffers[num] && buffers[num]->next->num < target) {
					tNum = num;
					node = buffers[num]->next;
					target = node->num;
				}
			}
			sfwrite(node->buff, 1, node->size, out);
			
			/* update cylinder */
			nodeN = node->prev;
			node0 = node->next;
			if(nodeN == node) {
				--allIn;
				nodeN = 0;
			} else {
				nodeN->next = node0;
				node0->prev = nodeN;
			}
			free(node->buff);
			free(node);
			buffers[tNum] = nodeN;	
		}
		sfwrite(&(unsigned){UINT_MAX}, sizeof(unsigned), 1, out);
		sfwrite(&(int){out_Tem[2] - 1}, sizeof(int), 1, out);
		return 0;
	}
	
	num = out_Tem[-1];
	infoSize[0] = num;
	infoSize[1] = qseq->seqlen;
	infoSize[2] = qseq->complen;
	infoSize[3] = qseq->N[0];
	infoSize[4] = -(abs(rc_flag));
	infoSize[5] = *out_Tem;
	infoSize[6] = header->len;
	infoSize[7] = flag;
	
	if(num == target || thread_num == 1) {
		sfwrite(infoSize, sizeof(int), 8, out);
		sfwrite(qseq->seq, sizeof(long unsigned), qseq->complen, out);
		if(qseq->N[0]) {
			sfwrite(qseq->N + 1, sizeof(int), qseq->N[0], out);
		}
		sfwrite(out_Tem + 1, sizeof(int), *out_Tem, out);
		sfwrite(header->seq, 1, header->len, out);
	} else {
		node = smalloc(sizeof(SpltDBbuff));
		size = header->len + sizeof(int) * (8 + qseq->N[0] + out_Tem[0]) + sizeof(long unsigned) * qseq->complen;
		buff = smalloc(size);
		node->num = num;
		node->size = size;
		node->buff = buff;
		memcpy(buff, infoSize, (size = sizeof(int) * 8));
		buff += size;
		memcpy(buff, qseq->seq, (size = sizeof(long unsigned) * qseq->complen));
		buff += size;
		memcpy(buff, qseq->N + 1, (size = sizeof(int) * qseq->N[0]));
		buff += size;
		memcpy(buff, out_Tem + 1, (size = sizeof(int) * (*out_Tem)));
		buff += size;
		memcpy(buff, header->seq, header->len);
		
		/* find cylinder piece */
		tNum = out_Tem[-2];
		/* extend cylinder piece */
		if((nodeN = buffers[tNum])) {
			node0 = nodeN->next;
			node->prev = nodeN;
			nodeN->next = node;
			node->next = node0;
			node0->prev = node;
		} else {
			node->next = node;
			node->prev = node;
			++allIn;
		}
		buffers[tNum] = node;
	}
	
	/* collect remains */
	while(allIn == thread_num) {
		/* find and write next target */
		node = (*buffers)->next;
		target = node->num;
		num = allIn;
		tNum = 0;
		while(--num) {
			if(buffers[num]->next->num < target) {
				tNum = num;
				node = buffers[num]->next;
				target = node->num;
			}
		}
		sfwrite(node->buff, 1, node->size, out);
		
		/* update cylinder */
		nodeN = node->prev;
		node0 = node->next;
		if(nodeN == node) {
			--allIn;
			nodeN = 0;
		} else {
			nodeN->next = node0;
			node0->prev = nodeN;
		}
		free(node->buff);
		free(node);
		buffers[tNum] = nodeN;
	}
	
	return 0;
}

unsigned get_ankers_spltDB(int *infoSize, int *out_Tem, CompDNA *qseq, Qseqs *header, int *flag, FILE *inputfile) {
	
	unsigned num;
	
	if(qseq) {
		qseq->seqlen = infoSize[0];
		qseq->complen = infoSize[1];
		*(qseq->N) = infoSize[2];
		*out_Tem = infoSize[4];
		header->len = infoSize[5];
		*flag = infoSize[6];
		
		/* reallocate */
		if(qseq->size <= qseq->seqlen) {
			free(qseq->N);
			free(qseq->seq);
			if(qseq->seqlen & 31) {
				qseq->size = (qseq->seqlen >> 5) + 1;
				qseq->size <<= 6;
			} else {
				qseq->size = qseq->seqlen << 1;
			}
			
			qseq->seq = calloc(qseq->size >> 5, sizeof(long unsigned));
			qseq->N = malloc((qseq->size + 1) * sizeof(int));
			if(!qseq->seq || !qseq->N) {
				ERROR();
			}
			*(qseq->N) = infoSize[2];
		}
		
		if(header->size < header->len) {
			free(header->seq);
			header->size = header->len;
			header->seq = malloc(header->size);
			if(!header->seq) {
				ERROR();
			}
		}
		
		cfread(qseq->seq, sizeof(long unsigned), qseq->complen, inputfile);
		cfread(qseq->N + 1, sizeof(int), qseq->N[0], inputfile);
		cfread(out_Tem + 1, sizeof(int), *out_Tem, inputfile);
		cfread(header->seq, 1, header->len, inputfile);
	} else {
		/* in the case of equally well scoring DBs */
		sfseek(inputfile, infoSize[1] * sizeof(long unsigned) + infoSize[2] * sizeof(int), SEEK_CUR);
		*out_Tem = infoSize[4];
		cfread(out_Tem + 1, sizeof(int), *out_Tem, inputfile);
		sfseek(inputfile, infoSize[5], SEEK_CUR);
	}
	
	/* get info for next read */
	cfread(&num, sizeof(int), 1, inputfile);
	if(num != UINT_MAX) {
		cfread(infoSize, sizeof(int), 7, inputfile);
	} else {
		cfread(infoSize, sizeof(int), 1, inputfile);
	}
	
	return num;
}

int runKMA_spltDB(char **templatefilenames, int targetNum, char *outputfilename, int argc, char **argv, int ConClave, int kmersize, int minlen, Penalties *rewards, int extendedFeatures, double ID_t, double Depth_t, int mq, double scoreT, double mrc, double evalue, double support, int bcd, int ref_fsa, int print_matrix, int print_all, long unsigned tsv, int vcf, int xml, int sam, int nc, int nf, unsigned shm, int thread_num, int maxFrag, int verbose) {
	
	/* https://www.youtube.com/watch?v=LtXEMwSG5-8 */
	
	int i, j, k, t_len, file_len, template, bestHits, end, aln_len, DB_size;
	int rc_flag, coverScore, status, delta, seq_in_no, progress, sparse;
	int fileCount, flag, flag_r;
	int *template_lengths, *bestTargets, *matched_templates, *bestTemplates;
	int *best_start_pos, *best_end_pos, (*targetInfo)[7], (*ptrInfo)[7];
	unsigned num, target, targetScore, bias;
	unsigned *fragmentCounts, *readCounts, *nums, *uPtr, *dbBiases;
	long best_read_score, read_score, seq_seeker;
	long unsigned Nhits, template_tot_ulen, counter, seqin_size;
	long unsigned *w_scores, *uniq_alignment_scores, *alignment_scores;
	double id, cover, q_id, q_cover, p_value;
	long double depth, q_value, expected;
	char *templatefilename, Date[11];
	FILE **inputfiles, *inputfile, *frag_in_raw;
	FILE *res_out, *tsv_out, *alignment_out, *consensus_out, *frag_out_raw;
	FILE *extendedFeatures_out, *name_file, **template_fragments, *xml_out;
	time_t t0, t1;
	struct tm *tm;
	FileBuff *frag_out, *frag_out_all, *matrix_out, *vcf_out;
	Aln *aligned, *gap_align;
	Assem *aligned_assem;
	Frag **alignFrags;
	CompDNA *qseq_comp, *qseq_r_comp;
	Qseqs *qseq, *qseq_r, *header, *header_r, *template_name;
	AssemInfo *matrix;
	AlnPoints *points;
	NWmat *NWmatrices;
	Assemble_thread *threads, *thread;
	HashMapCCI *template_index;
	
	if(!outputfilename) {
		fprintf(stderr, " No output file specified!\n");
		exit(1);
	}
	
	/* load databases */
	status = 0;
	inputfiles = smalloc(targetNum * sizeof(FILE *));
	nums = smalloc(targetNum * sizeof(unsigned));
	dbBiases = smalloc((targetNum + 1) * sizeof(unsigned)); /* set these */
	bestTargets = smalloc((targetNum + 1) * sizeof(int));
	targetInfo = smalloc(targetNum * 6 * sizeof(int));
	template_name = setQseqs(256);
	DB_size = 0;
	template_lengths = NULL;
	for(i = 0; i < targetNum; ++i) {
		/* get length, size and bias */
		templatefilename = templatefilenames[i];
		file_len = strlen(templatefilename);
		strcat(templatefilename, ".length.b");
		inputfile = sfopen(templatefilename, "rb");
		sfread(&bias, sizeof(int), 1, inputfile);
		dbBiases[i] = DB_size;
		DB_size += bias;
		template_lengths = realloc(template_lengths, DB_size * sizeof(int));
		if(!template_lengths) {
			ERROR();
		}
		sfread(template_lengths + dbBiases[i], sizeof(int), bias, inputfile);
		templatefilename[file_len] = 0;
		fclose(inputfile);
		if(sam) {
			strcat(templatefilename, ".name");
			name_file = sfopen(templatefilename, "rb");
			templatefilename[file_len] = 0;
			saminit(template_name, name_file, template_lengths + dbBiases[i], bias, 0);
			fclose(name_file);
		}
	}
	dbBiases[i] = DB_size;
	if(!kmersize) {
		kmersize = *template_lengths;
	}
	if(kmersize < 4 || 32 < kmersize) {
		kmersize = 16;
	}
	
	/* get scoring arrays */
	matched_templates = smalloc(((DB_size + 1) << 1) * sizeof(int));
	bestTemplates = (matched_templates + 1);
	best_start_pos = calloc((DB_size << 1), sizeof(int));
	best_end_pos = smalloc((DB_size << 1) * sizeof(int));
	alignment_scores = calloc(DB_size, sizeof(long unsigned));
	uniq_alignment_scores = calloc(DB_size, sizeof(long unsigned));
	if(!best_start_pos || !alignment_scores || !uniq_alignment_scores) {
		ERROR();
	}
	
	/* allocate stuff */
	qseq_comp = malloc(sizeof(CompDNA));
	qseq_r_comp = malloc(sizeof(CompDNA));
	if(!qseq_comp || !qseq_r_comp) {
		ERROR();
	}
	delta = 1024;
	allocComp(qseq_comp, delta);
	allocComp(qseq_r_comp, delta);
	qseq = setQseqs(delta);
	qseq_r = setQseqs(delta);
	header = setQseqs(256);
	header_r = setQseqs(256);
	points = seedPoint_init(delta, rewards);
	
	/* open outputfiles */
	file_len = strlen(outputfilename);
	strcat(outputfilename, ".res");
	res_out = sfopen(outputfilename, "w");
	outputfilename[file_len] = 0;
	if(tsv) {
		strcat(outputfilename, ".tsv");
		tsv_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
	} else {
		tsv_out = 0;
	}
	if(nf == 0) {
		strcat(outputfilename, ".frag.gz");
		frag_out = gzInitFileBuff(CHUNK);
		openFileBuff(frag_out, outputfilename, "wb");
		outputfilename[file_len] = 0;
	} else {
		frag_out = 0;
	}
	alignment_out = 0;
	consensus_out = 0;
	if(nc == 0) {
		strcat(outputfilename, ".aln");
		alignment_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".fsa");
		consensus_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
	} else if((nc & 1) == 0) {
		strcat(outputfilename, ".fsa");
		consensus_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
	}
	frag_out_raw = tmpF(0);
	if(!frag_out_raw) {
		ERROR();
	}
	if(print_matrix) {
		matrix_out = gzInitFileBuff(CHUNK);
		strcat(outputfilename, ".mat.gz");
		openFileBuff(matrix_out, outputfilename, "wb");
		outputfilename[file_len] = 0;
	} else {
		matrix_out = 0;
	}
	if(print_all) {
		strcat(outputfilename, ".frag_raw.gz");
		frag_out_all = gzInitFileBuff(CHUNK);
		openFileBuff(frag_out_all, outputfilename, "wb");
		outputfilename[file_len] = 0;
	} else {
		frag_out_all = 0;
	}
	if(vcf) {
		vcf_out = gzInitFileBuff(CHUNK);
		strcat(outputfilename, ".vcf.gz");
		openFileBuff(vcf_out, outputfilename, "wb");
		outputfilename[file_len] = 0;
	} else {
		vcf_out = 0;
	}
	if(extendedFeatures) {
		strcat(outputfilename, ".mapstat");
		extendedFeatures_out = sfopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
		fprintf(extendedFeatures_out, "## method\tKMA\n");
		fprintf(extendedFeatures_out, "## version\t%s\n", KMA_VERSION);
		
		fprintf(extendedFeatures_out, "## databases\t%s", noFolder(*templatefilenames));
		for(i = 1; i < targetNum; ++i) {
			fprintf(extendedFeatures_out, ", %s", noFolder(templatefilenames[i]));
		}
		fprintf(extendedFeatures_out, "\n");
		
		time(&t1);
		tm = localtime(&t1);
		strftime(Date, sizeof(Date), "%Y-%m-%d", tm);
		fprintf(extendedFeatures_out, "## date\t%s\n", Date);
		fprintf(extendedFeatures_out, "## command\t%s", *argv);
		for(i = 1; i < argc; ++i) {
			fprintf(extendedFeatures_out, " %s", argv[i]);
		}
		fprintf(extendedFeatures_out, "\n");
	} else {
		extendedFeatures_out = 0;
	}
	
	if(xml) {
		if(xml == 2) {
			xml_out = openInitXML("--", noFolder(*templatefilenames), **targetInfo, argc, argv);
		} else {
			strcat(outputfilename, ".xml");
			xml_out = openInitXML(outputfilename, noFolder(*templatefilenames), **targetInfo, argc, argv);
			outputfilename[file_len] = 0;
		}
	} else {
		xml_out = 0;
	}
	
	/* open input streams */
	file_len = strlen(outputfilename);
	for(i = 0; i < targetNum; ++i) {
		j = sprintf(outputfilename + file_len, ".%d", i);
		while(!(inputfile = fopen(outputfilename, "rb"))) {
			usleep(100);
		}
		setvbuf(inputfile, NULL, _IOFBF, CHUNK);
		inputfiles[i] = inputfile;
		outputfilename[file_len] = 0;
	}
	
	fprintf(stderr, "# Collecting k-mer scores.\n");
	t0 = clock();
	
	/* collect from DB mappings */
	uPtr = nums + (i = targetNum);
	ptrInfo = targetInfo + i;
	while(i--) {
		cfread(--uPtr, sizeof(unsigned), 1, inputfiles[i]);
		if(*uPtr != UINT_MAX) {
			cfread(--ptrInfo, sizeof(int), 7, inputfiles[i]);
		} else {
			cfread(--ptrInfo, sizeof(int), 1, inputfiles[i]);
			(*--ptrInfo)[3] = INT_MAX;
		}
	}
	
	Nhits = 0;
	target = 0;
	targetScore = 0;
	rc_flag = 0;
	qseq->len = 0;
	qseq_comp->seqlen = 0;
	*bestTargets = 0;
	while(target != UINT_MAX) {
		/* join best templates */
		read_score = 0;
		*matched_templates = 0;
		qseq->len = 0;
		qseq_r->len = 0;
		for(i = 1; i <= *bestTargets; ++i) {
			num = bestTargets[i];
			nums[num] = get_ankers_spltDB(targetInfo[num], matched_templates + (*matched_templates + 1), qseq_comp, header, &flag, inputfiles[num]);
			qseq->len = qseq_comp->seqlen;
			if(matched_templates[*matched_templates + 1]) {
				read_score = 0;
			} else { // PE
				//target = nums[num];
				nums[num] = get_ankers_spltDB(targetInfo[num], matched_templates + (*matched_templates + 1), qseq_r_comp, header_r, &flag_r, inputfiles[num]);
				qseq_r->len = qseq_r_comp->seqlen;
				read_score = 1;
			}
			
			/* bias the templates */
			bias = dbBiases[num];
			j = *matched_templates + 1;
			*matched_templates += matched_templates[*matched_templates + 1];
			while((k = j++) <= *matched_templates) {
				matched_templates[k] = matched_templates[j] + bias;
			}
		}
		
		if(kmersize <= qseq->len) {
			if(delta <= MAX(qseq->len, qseq_r->len)) {
				delta = MAX(qseq->len, qseq_r->len);
				delta <<= 1;
				qseq->size = delta;
				qseq_r->size = delta;
				free(qseq->seq);
				free(qseq_r->seq);
				qseq->seq = smalloc(delta);
				qseq_r->seq = smalloc(delta);
			}
			unCompDNA(qseq_comp, qseq->seq);
			
			best_read_score = targetScore;
			for(i = 1, bestHits = 0; i <= *matched_templates; ++i, ++bestHits) {
				best_end_pos[bestHits] = template_lengths[abs(matched_templates[i])];
			}
			
			if(rc_flag < 0 && 0 < matched_templates[*matched_templates]) {
				bestHits = -bestHits;
			}
			
			if(read_score && kmersize <= qseq_r->len) {
				unCompDNA(qseq_r_comp, qseq_r->seq);
				update_Scores_pe_MEM(qseq->seq, qseq->len, qseq_r->seq, qseq_r->len, bestHits, best_read_score + read_score, best_start_pos, best_end_pos, bestTemplates, header, header_r, flag, flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
			} else {
				update_Scores_MEM(qseq->seq, qseq->len, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
			}
			
			/* dump seq to all */
			if(frag_out_all) {
				updateAllFrag(qseq->seq, qseq->len, abs(bestHits), best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_all);
				if(read_score) {
					updateAllFrag(qseq_r->seq, qseq_r->len, abs(bestHits), read_score, best_start_pos, best_end_pos, bestTemplates, header_r, frag_out_all);
				}
			}
			
			/* verbose */
			if(verbose && verbose++ == 1024) {
				Nhits += verbose - 1;
				fprintf(stderr, "# Scored %ld query sequences.\n", Nhits);
				verbose = 1;
			}
		}
		
		/* remove inferior read matches */
		if(*matched_templates != 0) {
			if(read_score || (flag & 1) == 0 || (flag & 128)) { /* PE or SE or last in pair */
				/* wipe all from same fragment */
				i = targetNum - 1;
				uPtr = nums + i;
				while(i) {
					if(*uPtr == target) {
						*uPtr = get_ankers_spltDB(targetInfo[i], matched_templates, qseq_comp, header, &flag, inputfiles[i]);
					} else {
						--uPtr;
						--i;
					}
				}
			} else if(flag & 64) { /* first in pair */
				/* conserve second non-paired */
				i = targetNum - 1;
				uPtr = nums + i;
				while(i) {
					if(*uPtr == target) {
						if((targetInfo[i][6] & 128) && targetInfo[i][4] == 0) {
							--uPtr;
							--i;
						} else {
							*uPtr = get_ankers_spltDB(targetInfo[i], matched_templates, qseq_comp, header, &flag, inputfiles[i]);
							--uPtr;
							--i;
						}
					} else {
						--uPtr;
						--i;
					}
				}
			}
		}
		
		/* get best templates for next read */
		uPtr = nums;
		target = UINT_MAX;
		targetScore = UINT_MAX;
		rc_flag = 0;
		*bestTargets = 0;
		for(i = 0; i < targetNum; ++i) {
			if(*uPtr < target) {
				target = *uPtr;
				rc_flag = targetInfo[i][3];
				targetScore = abs(rc_flag);
				*bestTargets = 1;
				bestTargets[1] = i;
			} else if(*uPtr == target) {
				if(targetScore < abs(targetInfo[i][3])) {
					rc_flag = targetInfo[i][3];
					targetScore = abs(rc_flag);
					/* clear inferior matches */
					/*
					j = *bestTargets + 1;
					while(--j) {
						if(targetInfo[j][4] == 0) {
							uPtr[j] = get_ankers_spltDB(targetInfo[j], matched_templates, qseq_comp, header, &flag, inputfiles[j]);
						}
						uPtr[j] = get_ankers_spltDB(targetInfo[j], matched_templates, qseq_comp, header, &flag, inputfiles[j]);
					}
					*/
					*bestTargets = 1;
					bestTargets[1] = i;
				} else if(targetScore == abs(targetInfo[i][3])) {
					bestTargets[++*bestTargets] = i;
					rc_flag = rc_flag < 0 ? rc_flag : targetInfo[i][3];
				} else if(*uPtr == UINT_MAX) {
					targetInfo[i][3] = INT_MAX;
				} else {
					if(targetInfo[i][4] == 0) { /* PE */
						*uPtr = get_ankers_spltDB(targetInfo[i], matched_templates, qseq_comp, header, &flag, inputfiles[i]);
					}
					*uPtr = get_ankers_spltDB(targetInfo[i], matched_templates, qseq_comp, header, &flag, inputfiles[i]);
				}
			}
			++uPtr;
		}
	}
	/* verbose */
	if(verbose) {
		Nhits += verbose - 1;
		fprintf(stderr, "# Scored %ld query sequences.\n", Nhits);
		verbose = 1;
	}
	/* get fragmentCount */
	if(extendedFeatures) {
		fprintf(extendedFeatures_out, "## fragmentCount\t%u\n", **targetInfo);
		fprintf(extendedFeatures_out, "# refSequence\treadCount\tfragmentCount\tmapScoreSum\trefCoveredPositions\trefConsensusSum\tbpTotal\tdepthVariance\tnucHighDepthVariance\tdepthMax\tsnpSum\tinsertSum\tdeletionSum\n");
	}
	
	/* close files */
	i = targetNum;
	while(i--) {
		fclose(inputfiles[i]);
	}
	
	i = 0;
	sfwrite(&i, sizeof(int), 1, frag_out_raw);
	fflush(frag_out_raw);
	freeComp(qseq_comp);
	free(qseq_comp);
	freeComp(qseq_r_comp);
	free(qseq_r_comp);
	if(header->size < header_r->size) {
		destroyQseqs(header);
		header = header_r;
	} else {
		destroyQseqs(header_r);
	}
	header_r = 0;
	if(qseq->size < qseq_r->size) {
		destroyQseqs(qseq);
		qseq = qseq_r;
	} else {
		destroyQseqs(qseq_r);
	}
	qseq_r = 0;
	if(frag_out_all) {
		destroyGzFileBuff(frag_out_all);
	}
	t1 = clock();
	fprintf(stderr, "#\n# Time for score collecting:\t%.2f s.\n", difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Sort, output and select k-mer alignments.\n");
	t0 = clock();
	
	/* Get best template for each mapped deltamer/read */
	/* Best hit chosen as: highest mapping score then higest # unique maps */
	alignFrags = calloc(DB_size, sizeof(Frag*));
	w_scores = calloc(DB_size, sizeof(long unsigned));
	template_fragments = calloc(DB_size, sizeof(FILE*));
	if(!alignFrags || !w_scores || !template_fragments) {
		ERROR();
	}
	frag_in_raw = frag_out_raw;
	rewind(frag_in_raw);
	
	/* Patricks features */
	if(extendedFeatures || xml || tsv) {
		fragmentCounts = calloc(DB_size, sizeof(unsigned));
		readCounts = calloc(DB_size, sizeof(unsigned));
		if(!fragmentCounts || !readCounts) {
			ERROR();
		}
	} else {
		fragmentCounts = 0;
		readCounts = 0;
	}
	
	/* Get expected values */
	sparse = 0;
	template_tot_ulen = 0;
	i = DB_size;
	while(--i) {
		template_tot_ulen += template_lengths[i];
	}
	
	/* ConClave */
	if(ConClave == 1) {
		fileCount = ConClavePtr(frag_in_raw, &template_fragments, DB_size, maxFrag, w_scores, fragmentCounts, readCounts, alignment_scores, uniq_alignment_scores, template_lengths, header, qseq, bestTemplates, best_start_pos, best_end_pos, alignFrags);
	} else if(ConClave == 2) {
		fileCount = ConClave2Ptr(frag_in_raw, &template_fragments, DB_size, maxFrag, w_scores, fragmentCounts, readCounts, alignment_scores, uniq_alignment_scores, template_lengths, header, qseq, bestTemplates, best_start_pos, best_end_pos, alignFrags, template_tot_ulen, scoreT, evalue);	
	} else {
		fileCount = 0;
	}
	
	free(alignFrags);
	free(best_start_pos);
	free(best_end_pos);
	free(matched_templates);
	fclose(frag_out_raw);
	
	/* Get expected values */
	Nhits = 0;
	i = DB_size;
	while(--i) {
		Nhits += w_scores[i];
	}
	Nhits = Nhits ? Nhits : 1;
	
	t1 = clock();
	fprintf(stderr, "# Total time for sorting and outputting KMA alignment\t%.2f s.\n", difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Doing local assemblies of found templates, and output results\n");
	t0 = clock();
	
	/* print heading of resistance file: */
	fprintf(res_out, "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value\n");
	if(tsv) {
		initsv(tsv_out, tsv);
	}
	if(vcf) {
		templatefilename = 0;
		initialiseVcf(vcf_out, templatefilename);
	}
	
	/* preallocate assembly matrices */
	matrix = smalloc(sizeof(AssemInfo));
	aligned_assem = smalloc(sizeof(Assem));
	matrix->size = delta;
	for(i = 0; i < DB_size; ++i) {
		if(matrix->size < template_lengths[i]) {
			matrix->size = template_lengths[i];
		}
	}
	if(alnToMatPtr == &alnToMat) {
		matrix->size <<= 1;
	} else {
		matrix->size++;
	}
	matrix->assmb = smalloc(matrix->size * sizeof(Assembly));
	aligned_assem->size = matrix->size;
	aligned_assem->t = smalloc(aligned_assem->size);
	aligned_assem->s = smalloc(aligned_assem->size);
	aligned_assem->q = smalloc(aligned_assem->size);
	seq_in_no = 0;
	template_index = smalloc(sizeof(HashMapCCI));
	template_index->size = 0;
	hashMapCCI_initialize(template_index, matrix->size, kmersize);
	
	/* allocate matrcies for NW */
	i = 1;
	threads = 0;
	while(i < thread_num) {
		/* allocate matrices */
		NWmatrices = smalloc(sizeof(NWmat));
		NWmatrices->NW_s = 1024 * 1024;
		NWmatrices->NW_q = 1024;
		NWmatrices->E = smalloc(NWmatrices->NW_s);
		NWmatrices->D[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
		NWmatrices->P[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
		NWmatrices->D[1] = NWmatrices->D[0] + NWmatrices->NW_q;
		NWmatrices->P[1] = NWmatrices->P[0] + NWmatrices->NW_q;
		NWmatrices->rewards = rewards;
		
		aligned = smalloc(sizeof(Aln));
		gap_align = smalloc(sizeof(Aln));
		aligned->t = smalloc((delta + 1) << 1);
		aligned->s = smalloc((delta + 1) << 1);
		aligned->q = smalloc((delta + 1) << 1);
		gap_align->t = smalloc((delta + 1) << 1);
		gap_align->s = smalloc((delta + 1) << 1);
		gap_align->q = smalloc((delta + 1) << 1);
		
		/* move it to the thread */
		thread = smalloc(sizeof(Assemble_thread));
		thread->num = i;
		thread->thread_num = thread_num;
		thread->mq = mq;
		thread->minlen = minlen;
		thread->scoreT = scoreT;
		thread->mrc = mrc;
		thread->evalue = evalue;
		thread->bcd = bcd;
		thread->sam = sam;
		thread->ef = extendedFeatures;
		thread->seq_in = seq_in_no;
		thread->kmersize = kmersize;
		thread->template = -2;
		thread->file_count = fileCount;
		thread->files = template_fragments;
		thread->frag_out = frag_out;
		thread->xml_out = xml_out;
		thread->aligned_assem = aligned_assem;
		thread->aligned = aligned;
		thread->gap_align = gap_align;
		thread->NWmatrices = NWmatrices;
		thread->matrix = matrix;
		thread->qseq = setQseqs(qseq->size);
		thread->header = setQseqs(header->size);
		thread->points = seedPoint_init(delta, rewards);
		thread->points->len = 0;
		thread->spin = (sparse < 0) ? 10 : 100;
		thread->template_index = template_index;
		thread->next = threads;
		threads = thread;
		
		/* start thread */
		if((errno = pthread_create(&thread->id, NULL, assembly_KMA_Ptr, thread))) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			fprintf(stderr, "Will continue with %d threads.\n", i);
			threads = thread->next;
			free(thread);
			i = thread_num;
		} else {
			++i;
		}
	}
	
	/* start main thread */
	NWmatrices = smalloc(sizeof(NWmat));
	NWmatrices->NW_s = 1024 * 1024;
	NWmatrices->NW_q = 1024;
	NWmatrices->E = smalloc(NWmatrices->NW_s);
	NWmatrices->D[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
	NWmatrices->P[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
	NWmatrices->D[1] = NWmatrices->D[0] + NWmatrices->NW_q;
	NWmatrices->P[1] = NWmatrices->P[0] + NWmatrices->NW_q;
	NWmatrices->rewards = rewards;
	
	aligned = smalloc(sizeof(Aln));
	gap_align = smalloc(sizeof(Aln));
	aligned->t = smalloc((delta + 1) << 1);
	aligned->s = smalloc((delta + 1) << 1);
	aligned->q = smalloc((delta + 1) << 1);
	gap_align->t = smalloc((delta + 1) << 1);
	gap_align->s = smalloc((delta + 1) << 1);
	gap_align->q = smalloc((delta + 1) << 1);
	
	/* move it to the thread */
	thread = smalloc(sizeof(Assemble_thread));
	thread->num = 0;
	thread->thread_num = thread_num;
	thread->mq = mq;
	thread->minlen = minlen;
	thread->scoreT = scoreT;
	thread->mrc = mrc;
	thread->evalue = evalue;
	thread->bcd = bcd;
	thread->sam = sam;
	thread->ef = extendedFeatures;
	thread->seq_in = seq_in_no;
	thread->kmersize = kmersize;
	thread->template = 0;
	thread->file_count = fileCount;
	thread->files = template_fragments;
	thread->frag_out = frag_out;
	thread->xml_out = xml_out;
	thread->aligned_assem = aligned_assem;
	thread->aligned = aligned;
	thread->gap_align = gap_align;
	thread->NWmatrices = NWmatrices;
	thread->matrix = matrix;
	thread->qseq = qseq;
	thread->header = header;
	thread->points = points;
	thread->points->len = 0;
	thread->template_index = template_index;
	thread->next = 0;
	thread->spin = (sparse < 0) ? 10 : 100;
	
	/* Do local assemblies of fragments mapping to the same template */
	depth = 0;
	q_id = 0;
	cover = 0;
	q_cover = 0;
	seq_seeker = 0;
	progress = 0;
	counter = 0;
	if(progress) {
		fprintf(stderr, "# Progress:\t%3d%%\r", 0);
		fflush(stderr);
	} else if(verbose) {
		fprintf(stderr, "# Template\tScore\tProgress\n");
	}
	
	/* get index, seq and names on the fly */
	bias = *dbBiases++;
	templatefilename = *templatefilenames++;
	file_len = strlen(templatefilename);
	strcat(templatefilename, ".seq.b");
	seq_in_no = open(templatefilename, O_RDONLY);
	if(seq_in_no == -1) {
		ERROR();
	}
	seqin_size = 4 * lseek(seq_in_no, 0, SEEK_END);
	if(lseek(seq_in_no, 0, SEEK_SET) != 0) {
		ERROR();
	}
	thread->seq_in = seq_in_no;
	templatefilename[file_len] = 0;
	strcat(templatefilename, ".name");
	name_file = sfopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	alignLoadPtr = &alignLoad_fly_mem;
	if(kmersize < 4 || 32 < kmersize) {
		kmersize = 16;
	}
	
	if(extendedFeatures == 2) {
		printExtendedFeatures(templatefilename, 0, 0, 0, extendedFeatures_out);
	}
	if(assembly_KMA_Ptr == &skip_assemble_KMA) {
		alignLoadPtr = &alignLoad_skip;
	}
	
	for(template = 1; template < DB_size; ++template) {
		if(template == *dbBiases) {
			/* swap indexes */
			templatefilename = *templatefilenames++;
			bias = *dbBiases++;
			file_len = strlen(templatefilename);
			strcat(templatefilename, ".name");
			fclose(name_file);
			name_file = sfopen(templatefilename, "rb");
			templatefilename[file_len] = 0;
			strcat(templatefilename, ".seq.b");
			close(seq_in_no);
			seq_in_no = open(templatefilename, O_RDONLY);
			if(seq_in_no == -1) {
				ERROR();
			}
			seqin_size = 4 * lseek(seq_in_no, 0, SEEK_END);
			if(lseek(seq_in_no, 0, SEEK_SET) != 0) {
				ERROR();
			}
			thread->seq_in = seq_in_no;
			templatefilename[file_len] = 0;
			seq_seeker = 0;
			if(extendedFeatures == 2) {
				printExtendedFeatures(templatefilename, 0, 0, 0, extendedFeatures_out);
			}
		} else if(w_scores[template] > 0) {
			if(progress) {
				counter += w_scores[template];
				fprintf(stderr, "# Progress:\t%3lu%%\r", 100 * counter / Nhits);
				fflush(stderr);
			} else if(verbose) {
				counter += w_scores[template];
				fprintf(stderr, "# %d / %d\t%lu\t%3lu%%\n", template, DB_size, w_scores[template], 100 * counter / Nhits);
			}
			
			/* make p_value to see whether assembly is feasable */
			read_score = w_scores[template];
			t_len = template_lengths[template];
			expected = t_len;
			expected /= (template_tot_ulen - t_len);
			expected *= (Nhits - read_score);
			if(0 < expected) {
				q_value = read_score - expected;
				q_value /= (expected + read_score);
				q_value *= (read_score - expected);
			} else {
				q_value = read_score;
			}
			p_value  = p_chisqr(q_value);
			
			if(cmp((p_value <= evalue && read_score > expected), (read_score >= scoreT * t_len))) {
				/* load DB */
				thread->template_name = nameLoad(template_name, name_file);
				seq_seeker *= sizeof(long unsigned);
				lseek(seq_in_no, seq_seeker, SEEK_CUR);
				seq_seeker = 0;
				thread->template_index->len = 0;
				if(xml) {
					newIterXML(xml_out, template, t_len, thread->template_name);
				}
				
				/* Do assembly */
				//status |= assemblyPtr(aligned_assem, template, template_fragments, fileCount, frag_out, aligned, gap_align, qseq, header, matrix, points, NWmatrices);
				thread->template = template;
				thread->t_len = t_len;
				assembly_KMA_Ptr(thread);
				
				/* Depth, ID and coverage */
				if(aligned_assem->cover > 0) {
					coverScore = aligned_assem->cover;
					depth = aligned_assem->depth;
					depth /= t_len;
					id = 100.0 * coverScore / t_len;
					aln_len = aligned_assem->aln_len;
					q_id = 100.0 * coverScore / aln_len;
					cover = 100.0 * aln_len / t_len;
					q_cover = 100.0 * t_len / aln_len;
				} else {
					id = 0;
					q_id = 0;
					depth = aligned_assem->depth;
					depth /= t_len;
					aln_len = aligned_assem->aln_len;
					cover = 100.0 * aln_len / t_len;
					q_cover = 0;
				}
				
				if(xml) {
					capIterXML(xml_out, DB_size, seqin_size, t_len, readCounts[template], p_value, read_score, aligned_assem->q, aln_len);
				}
				
				if(ID_t <= id && Depth_t <= depth) {
					/* Output result */
					fprintf(res_out, "%s\t%8ld\t%8u\t%8d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
						thread->template_name, read_score, (unsigned) expected, t_len, id, cover, q_id, q_cover, (double) depth, (double) q_value, p_value);
					if(tsv) {
						printsv(tsv_out, tsv, thread->template_name, aligned_assem, t_len, readCounts[template], read_score, expected, q_value, p_value, alignment_scores[template], template_name->seq);
					}
					if(consensus_out) {
						printConsensus(aligned_assem, thread->template_name, alignment_out, consensus_out, ref_fsa);
					}
					/* print matrix */
					if(matrix_out) {
						updateMatrix(matrix_out, thread->template_name, thread->template_index->seq, matrix, t_len);
					}
					if(extendedFeatures) {
						printExtendedFeatures(thread->template_name, aligned_assem, fragmentCounts[template], readCounts[template], extendedFeatures_out);
					}
					if(vcf) {
						updateVcf(thread->template_name, aligned_assem->t, evalue, support, bcd, t_len, matrix, vcf, vcf_out);
					}
				}
			} else {
				if(sam || ID_t == 0.0) {
					seq_seeker *= sizeof(long unsigned);
					lseek(seq_in_no, seq_seeker, SEEK_CUR);
					seq_seeker = 0;
					thread->template_index = alignLoadPtr(thread->template_index, seq_in_no, template_lengths[template], kmersize, 0);
					thread->template_name = nameLoad(template_name, name_file);
					skip_assemble_KMA(thread);
					if(ID_t == 0.0) {
						depth = aligned_assem->depth;
						depth /= t_len;
						fprintf(res_out, "%s\t%8ld\t%8u\t%8d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
							thread->template_name, read_score, (unsigned) expected, t_len, 0.0, 0.0, 0.0, 0.0, (double) depth, (double) q_value, p_value);
						if(tsv) {
							printsv(tsv_out, tsv, thread->template_name, aligned_assem, t_len, readCounts[template], read_score, expected, q_value, p_value, alignment_scores[template], template_name->seq);
						}
						if(extendedFeatures) {
							printExtendedFeatures(thread->template_name, aligned_assem, fragmentCounts[template], readCounts[template], extendedFeatures_out);
						}
					}
				} else {
					nameSkip(name_file, end);
				}
				seq_seeker += ((template_lengths[template] >> 5) + 1);
			}
		} else {
			nameSkip(name_file, end);
			seq_seeker += ((template_lengths[template] >> 5) + 1);	
		}		
	}
	
	if(progress) {
		fprintf(stderr, "\n");
	}
	
	/* join threads */
	thread->template = -1;
	assembly_KMA_Ptr(thread);
	for(thread = threads; thread != 0; thread = thread->next) {
		/* join thread */
		if((errno = pthread_join(thread->id, NULL))) {
			ERROR();
		}
	}
	
	/* Close files */
	close(seq_in_no);
	fclose(res_out);
	if(tsv) {
		fclose(tsv_out);
	}
	if(alignment_out) {
		fclose(alignment_out);
		fclose(consensus_out);
	}
	fclose(name_file);
	if(frag_out) {
		destroyGzFileBuff(frag_out);
	}
	if(matrix_out) {
		destroyGzFileBuff(matrix_out);
	}
	if(extendedFeatures) {
		fclose(extendedFeatures_out);
	}
	if(vcf) {
		destroyGzFileBuff(vcf_out);
	}
	if(xml) {
		closeCapXML(xml_out);
	}
	
	t1 = clock();
	fprintf(stderr, "# Total time used for local assembly: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	
	return status;
}
