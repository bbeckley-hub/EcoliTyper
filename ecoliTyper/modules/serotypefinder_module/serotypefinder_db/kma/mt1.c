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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <time.h>
#include "assembly.h"
#include "chain.h"
#include "compdna.h"
#include "filebuff.h"
#include "hashmapcci.h"
#include "kmapipe.h"
#include "mt1.h"
#include "nw.h"
#include "penalties.h"
#include "pherror.h"
#include "printconsensus.h"
#include "qseqs.h"
#include "runkma.h"
#include "stdnuc.h"
#include "stdstat.h"
#include "tsv.h"
#include "vcf.h"
#include "version.h"
#include "xml.h"

void printFsaMt1(Qseqs *header, Qseqs *qseq, Qseqs *qual, CompDNA *compressor, FILE *out) {
	
	static int buff[8] = {0, 0, 1, 0, 0, 0, 0, 0};
	
	if(header) {
		buff[1] = qseq->len;
		buff[6] = header->len;
		sfwrite(buff, sizeof(int), 8, out);
		sfwrite(qseq->seq, 1, qseq->len, out);
		sfwrite(header->seq + 1, 1, header->len, out);
	} else {
		buff[5] = qseq->len;
	}
}

void printFsa_pairMt1(Qseqs *header, Qseqs *qseq, Qseqs *qual, Qseqs *header_r, Qseqs *qseq_r, Qseqs *qual_r, CompDNA *compressor, FILE *out) {
	
	static int buff[8] = {0, 0, 1, 0, 0, 0, 0, 0};
	
	if(header) {
		buff[1] = qseq->len;
		buff[6] = header->len;
		buff[7] = 97; /* 64 32 1 FirstInPair MateRC PE */
		sfwrite(buff, sizeof(int), 8, out);
		sfwrite(qseq->seq, 1, qseq->len, out);
		sfwrite(header->seq + 1, 1, header->len, out);
		
		buff[1] = qseq_r->len;
		buff[6] = header_r->len;
		buff[7] = 145; /* 128 16 1 LastInPair RC PE */
		strrc(qseq_r->seq, qseq_r->len);
		sfwrite(buff, sizeof(int), 8, out);
		sfwrite(qseq_r->seq, 1, qseq_r->len, out);
		sfwrite(header_r->seq + 1, 1, header_r->len, out);
	} else {
		buff[5] = qseq->len;
	}
	
}

void runKMA_Mt1(char *templatefilename, char *outputfilename, char *exePrev, int kmersize, int minlen, Penalties *rewards, double ID_t, double Depth_t, int mq, double scoreT, double mrc, double evalue, double support, int bcd, int Mt1, int ref_fsa, int print_matrix, long unsigned tsv, int vcf, int xml, int sam, int nc, int nf, int thread_num) {
	
	int i, j, aln_len, t_len, coverScore, file_len, DB_size, delta, seq_in;
	int *template_lengths;
	long unsigned read_score, seeker;
	double p_value, id, q_id, cover, q_cover;
	long double depth;
	FILE *res_out, *tsv_out, *xml_out, *alignment_out, *consensus_out;
	FILE *template_fragments, *DB_file;
	time_t t0, t1;
	FileBuff *frag_out, *matrix_out, *vcf_out;
	Aln *aligned, *gap_align;
	Assem *aligned_assem;
	Qseqs *qseq, *header, *template_name;
	AssemInfo *matrix;
	AlnPoints *points;
	NWmat *NWmatrices;
	Assemble_thread *threads, *thread;
	HashMapCCI *template_index;
	
	/* open pipe */
	//template_fragments = popen(exePrev, "r");
	template_fragments = kmaPipe("-s1", "rb", 0, 0);
	if(!template_fragments) {
		ERROR();
	} else {
		setvbuf(template_fragments, NULL, _IOFBF, CHUNK);
	}
	
	file_len = strlen(outputfilename);
	delta = 1024;
	header = setQseqs(256);
	qseq = setQseqs(delta);
	points = seedPoint_init(delta, rewards);
	
	/* open outputfiles */
	if(outputfilename) {
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
		if(print_matrix) {
			matrix_out = gzInitFileBuff(CHUNK);
			strcat(outputfilename, ".mat.gz");
			openFileBuff(matrix_out, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			matrix_out = 0;
		}
		if(vcf) {
			vcf_out = gzInitFileBuff(CHUNK);
			strcat(outputfilename, ".vcf.gz");
			openFileBuff(vcf_out, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			vcf_out = 0;
		}
		if(xml) {
			if(xml == 2) {
				xml_out = openInitXML("--", templatefilename, 1, 1, &exePrev);
			} else {
				strcat(outputfilename, ".xml");
				xml_out = openInitXML(outputfilename, templatefilename, 1, 1, &exePrev);
				outputfilename[file_len] = 0;
			}
		} else {
			xml_out = 0;
		}
	} else {
		fprintf(stderr, " No output file specified!\n");
		exit(1);
	}
	
	/* load indexing */
	file_len = strlen(templatefilename);
	strcat(templatefilename, ".length.b");
	DB_file = sfopen(templatefilename, "rb");
	sfread(&DB_size, sizeof(int), 1, DB_file);
	
	/* load lengths */
	template_lengths = smalloc(DB_size * sizeof(int));
	sfread(template_lengths, sizeof(int), DB_size, DB_file);
	/*sfseek(DB_file, (2 * DB_size) * sizeof(int), SEEK_CUR);
	if(fread(template_lengths, sizeof(int), DB_size, DB_file) == 0) {
		sfseek(DB_file, sizeof(int), SEEK_SET);
		sfread(template_lengths, sizeof(int), DB_size, DB_file);
	}*/
	templatefilename[file_len] = 0;
	fclose(DB_file);
	if(kmersize < 4) {
		kmersize = *template_lengths;
		if(32 < kmersize || kmersize < 4) {
			kmersize = 16;
		}
	}
	
	/* get seq / index */
	seeker = 0;
	for(i = 2; i <= Mt1; ++i) {
		seeker += ((template_lengths[i - 1] >> 5) + 1) * sizeof(long unsigned);
	}
	
	/* make index */
	*template_lengths = template_lengths[Mt1];
	template_lengths = realloc(template_lengths, sizeof(int));
	if(!template_lengths) {
		ERROR();
	}
	
	strcat(templatefilename, ".seq.b");
	seq_in = open(templatefilename, O_RDONLY);
	if(seq_in == -1) {
		ERROR();
	}
	templatefilename[file_len] = 0;
	template_index = alignLoad_fly(0, seq_in, *template_lengths, kmersize, seeker);
	close(seq_in);
	
	/* get name */
	strcat(templatefilename, ".name");
	DB_file = sfopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	template_name = setQseqs(256);
	i = 1;
	while(i != Mt1 && (j = fgetc(DB_file)) && j != EOF) {
		if(j == '\n') {
			++i;
		}
	}
	nameLoad(template_name, DB_file);
	fclose(DB_file);
	
	if(sam) {
		fprintf(stdout, "@PG\tID:KMA\tPN:kma\tVN:%s\tCL:%s\n", KMA_VERSION, exePrev);
		fprintf(stdout, "@SQ\tSN:%s\tLN:%d\n", template_name->seq, *template_lengths);
	}
	
	fprintf(stderr, "#\n# Doing local assemblies of found templates, and output results\n");
	t0 = clock();
	
	/* print heading of resistance file: */
	fprintf(res_out, "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value\n");
	if(tsv) {
		initsv(tsv_out, tsv);
	}
	if(vcf) {
		initialiseVcf(vcf_out, templatefilename);
	}
	
	/* preallocate assembly matrices */
	matrix = smalloc(sizeof(AssemInfo));
	aligned_assem = smalloc(sizeof(Assem));
	if(alnToMatPtr == &alnToMat) {
		matrix->size = (*template_lengths) << 1;
	} else {
		matrix->size = (*template_lengths) + 1;
	}
	matrix->assmb = smalloc(matrix->size * sizeof(Assembly));
	aligned_assem->size = matrix->size;
	aligned_assem->t = smalloc(aligned_assem->size);
	aligned_assem->s = smalloc(aligned_assem->size);
	aligned_assem->q = smalloc(aligned_assem->size);
	
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
		thread->ef = 0;
		thread->seq_in = 0;
		thread->kmersize = kmersize;
		thread->template = -2;
		thread->file_count = 1;
		thread->files = &template_fragments;
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
		thread->spin = 10;
		
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
	thread->ef = 0;
	thread->seq_in = 0;
	thread->kmersize = kmersize;
	thread->template = 0;
	thread->file_count = 1;
	thread->files = &template_fragments;
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
	thread->next = 0;
	thread->spin = 10;
	
	/* Do local assemblies of fragments mapping to the same template */
	depth = 0;
	q_id = 0;
	cover = 0;
	q_cover = 0;
	/* Do assembly */
	//assemblyPtr(aligned_assem, 0, &template_fragments, 1, frag_out, aligned, gap_align, qseq, header, matrix, points, NWmatrices);
	thread->template_name = (char *) template_name->seq;
	thread->template_index = template_index;
	if(xml) {
		newIterXML(xml_out, Mt1, *template_lengths, thread->template_name);
	}
	thread->t_len = *template_lengths;
	assembly_KMA_Ptr(thread);
	
	/* make p_value */
	read_score = aligned_assem->score;
	t_len = *template_lengths;
	p_value = p_chisqr(read_score);
	aln_len = 0;
	
	if(cmp((p_value <= evalue && read_score > 0), read_score >= scoreT * t_len)) {
		
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
		}
		if(ID_t <= id && 0 < id && Depth_t <= depth) {
			/* Output result */
			fprintf(res_out, "%s\t%8lu\t%8d\t%8d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
				thread->template_name, read_score, 0, t_len, id, cover, q_id, q_cover, (double) depth, (double) read_score, p_value);
			if(tsv) {
				printsv(tsv_out, tsv, thread->template_name, aligned_assem, t_len, aligned_assem->readCountAln, read_score, 0, (double) read_score, p_value, read_score, template_name->seq);
			}
			if(consensus_out) {
				printConsensus(aligned_assem, thread->template_name, alignment_out, consensus_out, ref_fsa);
			}
			/* print matrix */
			if(matrix_out) {
				updateMatrix(matrix_out, thread->template_name, template_index->seq, matrix, t_len);
			}
			if(vcf) {
				updateVcf(thread->template_name, aligned_assem->t, evalue, support, bcd, t_len, matrix, vcf, vcf_out);
			}
		}
		/* destroy this DB index */
		hashMapCCI_destroy(template_index);
	} else if(ID_t == 0.0) {
		fprintf(res_out, "%s\t%8ld\t%8u\t%8d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
				thread->template_name, read_score, 0, t_len, 0.0, 0.0, 0.0, 0.0, (double) depth, (double) read_score, p_value);
		if(tsv) {
			printsv(tsv_out, tsv, thread->template_name, aligned_assem, t_len, aligned_assem->readCountAln, read_score, 0, (double) read_score, p_value, read_score, template_name->seq);
		}
	}
	
	if(xml) {
		capIterXML(xml_out, 1, t_len, t_len, read_score, p_value, read_score, aligned_assem->q, aln_len);
		closeCapXML(xml_out);
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
	fclose(res_out);
	if(tsv) {
		fclose(tsv_out);
	}
	if(consensus_out) {
		fclose(consensus_out);
	}
	if(alignment_out) {
		fclose(alignment_out);
	}
	if(frag_out) {
		destroyGzFileBuff(frag_out);
	}
	if(matrix_out) {
		destroyGzFileBuff(matrix_out);
	}
	if(vcf) {
		destroyGzFileBuff(vcf_out);
	}
	
	t1 = clock();
	fprintf(stderr, "# Total time used for local assembly: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
}
