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
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ankers.h"
#include "compkmers.h"
#include "filebuff.h"
#include "hashmapkma.h"
#include "hashtable.h"
#include "kmapipe.h"
#include "kmmap.h"
#include "pherror.h"
#include "runinput.h"
#include "qc.h"
#include "qseqs.h"
#include "savekmers.h"
#include "seqparse.h"
#include "sparse.h"
#include "stdnuc.h"
#include "stdstat.h"
#ifdef _WIN32
typedef int key_t;
#define ftok(charPtr, integer) (0)
#define shmget(key, size, permission) ((size != 0) ? (-1) : (-key))
#define shmat(shmid, NULL_Ptr, integer) (NULL)
#else
#include <sys/ipc.h>
#include <sys/shm.h>
#endif

int translateToKmersAndDump(long unsigned *Kmers, int n, int max, unsigned char *qseq, int seqlen, const HashMapKMA *templates, FILE *out) {
	
	int i, end, rc, prefix_len, kmersize, shifter, mlen, flag, hLen, mPos;
	long unsigned prefix, pmask, mask, mmask, pmer, kmer, hmer, cmer;
	
	/* init */
	prefix = templates->prefix;
	prefix_len = templates->prefix_len;
	pmask = 0xFFFFFFFFFFFFFFFF >> (64 - (prefix_len << 1));
	kmersize = templates->kmersize;
	shifter = 64 - (kmersize << 1);
	mask = 0xFFFFFFFFFFFFFFFF >> shifter;
	mlen = templates->mlen;
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	flag = templates->flag;
	hLen = kmersize;
	pmer = 0;
	
	if(prefix_len) {
		for(rc = 0; rc < 2; ++rc) {
			i = 0;
			while(i < seqlen) {
				end = charpos(qseq, 4, i, seqlen);
				if(end == -1) {
					end = seqlen;
				}
				if(i < end - kmersize - prefix_len) {
					pmer = makeKmer(qseq, i, prefix_len - 1);
					i += (prefix_len - 1);
					end -= kmersize;
				} else {
					i = end + 1;
				}
				
				while(i < end) {
					pmer = ((pmer << 2) | qseq[i]) & pmask;
					++i;
					if(pmer == prefix) {
						kmer = makeKmer(qseq, i, kmersize);
						Kmers[n] = flag ? getCmer(kmer, &mPos, &hLen, shifter, mlen, mmask) : kmer;
						if(++n == max) {
							sfwrite(Kmers, sizeof(long unsigned), n, out);
							n = 0;
						}
					}
				}
				i = end + kmersize + 1;
			}
			
			/* rc */
			strrc(qseq, seqlen);
		}
	} else {
		for(rc = 0; rc < 2; ++rc) {
			i = 0;
			while(i < seqlen) {
				end = charpos(qseq, 4, i, seqlen);
				if(end == -1) {
					end = seqlen;
				}
				kmer = makeKmer(qseq, i, kmersize - 1);
				cmer = flag ? initCmer(kmer, &mPos, &hmer, &hLen, shifter + 2, kmersize, mlen, mmask) : kmer;
				
				i += kmersize - 1;
				while(i < end) {
					kmer = ((kmer << 2) | qseq[i]) & mask;
					cmer = flag ? updateCmer(cmer, &mPos, &hmer, &hLen, kmer, kmersize, mlen, mmask) : kmer;
					Kmers[n] = cmer;
					if(++n == max) {
						sfwrite(Kmers, sizeof(long unsigned), n, out);
						n = 0;
					}
					++i;
				}
				i = end + kmersize + 1;
			}
			/* rc */
			strrc(qseq, seqlen);
		}
	}
	return n;
}

char ** load_DBs_Sparse(char *templatefilename, unsigned **template_lengths, unsigned **template_ulengths, unsigned shm) {
	
	/* load DBs needed for KMA */
	int file_len, shmid, DB_size;
	long unsigned i, j, file_size;
	char **template_names;
	FILE *DB_file;
	key_t key;
	
	/* Open DB */
	file_len = strlen(templatefilename);
	strcat(templatefilename, ".length.b");
	templatefilename[file_len + 9] = 0;
	DB_file = sfopen(templatefilename, "rb");
	
	/* allocate DBs */
	sfread(&DB_size, sizeof(int), 1, DB_file);
	if(shm & 4) {
		sfseek(DB_file, 0, SEEK_END);
		file_size = ftell(DB_file) - sizeof(int);
		key = ftok(templatefilename, 'l');
		shmid = shmget(key, file_size, 0666);
		if(shmid < 0) {
			fprintf(stderr, "No shared length\n");
			exit(errno);
		} else {
			*template_lengths = shmat(shmid, NULL, 0);
		}
		*template_lengths += DB_size;
		*template_ulengths = *template_lengths + DB_size;
	} else {
		*template_lengths = smalloc(DB_size * sizeof(int));
		*template_ulengths = smalloc(DB_size * sizeof(int));
		
		/* load lengths */
		sfseek(DB_file, DB_size * sizeof(int), SEEK_CUR);
		i = fread(*template_lengths, sizeof(int), DB_size, DB_file);
		i += fread(*template_ulengths, sizeof(int), DB_size, DB_file);	
		if(!i) {
			fprintf(stderr, "DB needs to sparse indexed, to run a sparse mapping.\n");
			exit(1);
		}
	}
	templatefilename[file_len] = 0;
	fclose(DB_file);
	
	/* load names */
	strcat(templatefilename, ".name");
	templatefilename[file_len + 5] = 0;
	DB_file = sfopen(templatefilename, "rb");
	
	/* get size of file */
	sfseek(DB_file, 0, SEEK_END);
	file_size = ftell(DB_file);
	template_names = malloc(DB_size * sizeof(char*));
	if(!template_names) {
		ERROR();
	}
	
	/* load file */
	if(shm & 16) {
		key = ftok(templatefilename, 'n');
		shmid = shmget(key, file_size, 0666);
		if(shmid < 0) {
			fprintf(stderr, "No shared name\n");
			ERROR();
		} else {
			template_names[0] = shmat(shmid, NULL, 0);
		}
		for(i = 0, j = 2; j < DB_size; ++i) {
			if(template_names[0][i] == 0) {
				template_names[j] = template_names[0] + i + 1;
				++j;
			}
		}
	} else {
		rewind(DB_file);
		template_names[0] = malloc(file_size);
		if(!template_names[0]) {
			ERROR();
		}
		sfread(template_names[0], 1, file_size, DB_file);
		template_names[0][file_size - 1] = 0;
		
		template_names[1] = template_names[0];
		for(i = 0, j = 2; j < DB_size; ++i) {
			if(template_names[0][i] == '\n') {
				template_names[0][i] = 0;
				template_names[j] = template_names[0] + i + 1;
				++j;
			}
		}
	}
	templatefilename[file_len] = 0;
	fclose(DB_file);
	
	return template_names;
}

void save_kmers_sparse(const HashMapKMA *templates, HashMap_kmers *foundKmers, CompKmers *compressor) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	and is the time determining step */
	
	int i;
	
	for(i = 0; i < compressor->n; ++i) {
		if(hashMap_get(templates, compressor->kmers[i])) {
			hashMap_kmers_CountIndex(foundKmers, compressor->kmers[i]);
		}
	}
}

void run_input_sparse(const HashMapKMA *templates, char **inputfiles, int fileCount, int minPhred, int hardmaskQ, int minQ, int fiveClip, int threeClip, int minlen, int maxlen, int kmersize, char *trans, const double *prob, QCstat *qcreport, FILE *out) {
	
	int fileCounter, phredScale, start, end, len;
	unsigned FASTQ;
	long unsigned count, org_count;
	char *filename;
	Qseqs *qseq, *qual;
	FileBuff *inputfile;
	CompKmers *Kmers;
	
	if(minPhred < minQ) {
		minPhred = minQ;
	}
	Kmers = smalloc(sizeof(CompKmers));
	allocCompKmers(Kmers, 1024);
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
		Kmers->n = 0;
		if(FASTQ & 1) {
			/* get phred scale */
			phredScale = getPhredFileBuff(inputfile);
			fprintf(stderr, "# Phred scale:\t%d\n", phredScale);
			
			/* parse reads */
			while(FileBuffgetFqSeq(inputfile, qseq, qual, trans)) {
				++org_count;
				len = phredStat(qseq->seq, qual->seq, qseq->len, prob - phredScale, phredScale + minPhred, minQ, hardmaskQ, fiveClip, threeClip, minlen, maxlen, &start, &end, qcreport);
				
				/* print */
				if(minlen <= len) {
					/* translate to kmers */
					qseq->len = end - start;
					Kmers->n = translateToKmersAndDump(Kmers->kmers, Kmers->n, Kmers->size, qseq->seq + start, qseq->len, templates, out);
					++count;
					
				}
			}
			if(Kmers->n) {
				sfwrite(Kmers->kmers, sizeof(long unsigned), Kmers->n, out);
				Kmers->n = 0;
			}
		} else if(FASTQ & 2) {
			while(FileBuffgetFsaSeq(inputfile, qseq, trans)) {
				++org_count;
				len = fsastat(qseq->seq, qseq->len, minlen, maxlen, &start, &end, qcreport);
				
				/* translate to kmers */
				if(minlen <= len) {
					Kmers->n = translateToKmersAndDump(Kmers->kmers, Kmers->n, Kmers->size, qseq->seq + start, end - start, templates, out);
					++count;
				}
			}
			if(Kmers->n) {
				sfwrite(Kmers->kmers, sizeof(long unsigned), Kmers->n, out);
				Kmers->n = 0;
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
	free(Kmers->kmers);
	free(Kmers);
	destroyQseqs(qseq);
	destroyQseqs(qual);
	destroyFileBuff(inputfile);
}

int save_kmers_sparse_batch(char *templatefilename, char *outputfilename, char *exePrev, double ID_t, double Depth_t, double evalue, char ss, unsigned shm) {
	
	int i, file_len, stop, template, status, contamination, score_add, deCon;
	unsigned *Scores, *w_Scores, *SearchList;
	unsigned *template_lengths, *template_ulengths;
	long unsigned Ntot, score, tmp_score, score_tot_add;
	long unsigned *Scores_tot, *w_Scores_tot;
	char **template_names;
	double expected, q_value, p_value, etta, depth, cover, query_cover;
	double tot_depth, tot_cover, tot_query_cover, tmp_depth, tmp_cover;
	double tmp_expected, tmp_q, tmp_p;
	FILE *inputfile, *templatefile, *sparse_out;
	time_t t0, t1;
	HashMapKMA *templates;
	HashMap_kmers *foundKmers;
	HashTable *kmerList, *deConTable, *node, *prev, **Collecter;
	Hit Nhits, w_Nhits;
	CompKmers *Kmers;
	
	/* here */
	/* split input up */
	
	/* open pipe */
	//inputfile = popen(exePrev, "r");
	status = 0;
	inputfile = kmaPipe(exePrev, "rb", 0, 0);
	if(!inputfile) {
		ERROR();
	}
	
	/* Load hashMap */
	t0 = clock();
	file_len = strlen(templatefilename);
	if((deCon = deConPrintPtr == deConPrint)) {
		strcat(templatefilename, ".decon.comp.b");
	} else {
		strcat(templatefilename, ".comp.b");
	}
	templatefile = sfopen(templatefilename, "rb" );
	templates = smalloc(sizeof(HashMapKMA));
	hashMap_get = &hashMap_getGlobal;
	if((shm & 1) || (deCon && (shm & 2))) {
		hashMapKMA_load_shm(templates, templatefile, templatefilename);
	} else if(shm & 32) {
		hashMapKMAmmap(templates, templatefile);
	} else {
		if(hashMapKMA_load(templates, templatefile, templatefilename)) {
			fprintf(stderr, "Wrong format of DB.\n");
			exit(1);
		}
	}
	fclose(templatefile);
	templatefilename[file_len] = 0;
	
	/* load template attributes */
	template_names = load_DBs_Sparse(templatefilename, &template_lengths, &template_ulengths, shm);
	
	/* open output file */
	if(strcmp(outputfilename, "--") == 0) {
		sparse_out = stdout;
	} else {
		file_len = strlen(outputfilename);
		strcat(outputfilename, ".spa");
		sparse_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
	}
	if(!sparse_out) {
		ERROR();
	}
	
	/* set hashMap for found kmers */
	foundKmers = smalloc(sizeof(HashMap_kmers));
	/* foundKmers->size = templates->n; */
	foundKmers->flag = templates->flag;
	foundKmers->size = *template_lengths;
	i = templates->DB_size;
	while(--i) {
		if(foundKmers->size < template_lengths[i]) {
			foundKmers->size = template_lengths[i];
		}
	}
	foundKmers->size <<= 1;
	foundKmers->size = MIN(1048576, foundKmers->size);
	hashMap_kmers_initialize(foundKmers, foundKmers->size);
	
	t1 = clock();
	fprintf(stderr, "#\n# Total time used for DB loading: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	t0 = clock();
	fprintf(stderr, "# Finding k-mers\n");
	
	Kmers = smalloc(sizeof(CompKmers));
	allocCompKmers(Kmers, 1024);
	Ntot = 0;
	/* count kmers */
	while((Kmers->n = fread(Kmers->kmers, sizeof(long unsigned), Kmers->size, inputfile))) {
		Ntot += Kmers->n;
		save_kmers_sparse(templates, foundKmers, Kmers);
	}
	kmaPipe(0, 0, inputfile, &status);
	if(kmaPipe == &kmaPipeFork) {
		t1 = clock();
		fprintf(stderr, "#\n# Total time used to identify k-mers in query: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	} else {
		fprintf(stderr, "# k-mers in query identified\n#\n");
	}
	t0 = clock();
	fprintf(stderr, "# Finding best matches and output results.\n");
	
	Scores = calloc(templates->DB_size, sizeof(unsigned));
	Scores_tot = calloc(templates->DB_size, sizeof(long unsigned));
	SearchList = malloc(templates->DB_size * sizeof(unsigned));
	if(!Scores || !Scores_tot || !SearchList) {
		ERROR();
	}
	etta = 1.0e-6;
	expected = 0;
	q_value = 0;
	p_value = 0;
	
	fprintf(sparse_out, "#Template\tNum\tScore\tExpected\tTemplate_length\tQuery_Coverage\tTemplate_Coverage\tDepth\ttot_query_Coverage\ttot_template_Coverage\ttot_depth\tq_value\tp_value\n");
	stop = 0;
	if(ss == 'n') {
		/* collect scores */
		kmerList = collect_Kmers(templates, Scores, Scores_tot, foundKmers, &Nhits);
		
		fprintf(stderr, "# Total number of matches: %lu of %lu kmers\n", Nhits.tot, Ntot);
		
		/* get hits */
		for(template = 0; template < templates->DB_size; ++template) {
			score = Scores_tot[template];
			cover = 100.0 * Scores[template] / template_ulengths[template];
			depth = 1.0 * score / template_lengths[template];
			
			expected = 1.0 * (Nhits.tot - score) * template_ulengths[template] / (templates->n - template_ulengths[template] + etta);
			q_value = (score - expected) * (score - expected) / (score + expected);
			p_value = p_chisqr(q_value);
			/* validate best match */
			if(cover && ID_t <= cover && Depth_t <= depth) {
				/* output results */
				query_cover = 100.0 * Scores_tot[template] / Ntot;
				/* calc tot values */
				tot_cover = 100.0 * Scores[template] / template_ulengths[template];
				tot_depth = 1.0 * Scores_tot[template] / template_lengths[template];
				tot_query_cover = 100.0 * Scores_tot[template] / Ntot;
				fprintf(sparse_out, "%s\t%d\t%lu\t%d\t%d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n", 
					template_names[template], template, score, (int) expected, template_ulengths[template], query_cover, cover, depth, tot_query_cover, tot_cover, tot_depth, q_value, p_value);
			}
		}
	} else if(deCon) {
		/* start by removing contamination and collect scores */
		contamination = templates->DB_size;
		Collecter = collect_Kmers_deCon(templates, Scores, Scores_tot, foundKmers, &Nhits, contamination);
		kmerList = Collecter[0];
		deConTable = Collecter[1];
		free(Collecter);
		
		fprintf(stderr, "# Total number of matches: %lu of %lu kmers\n", Nhits.tot, Ntot);
		/* copy scores */
		w_Scores = smalloc(templates->DB_size * sizeof(unsigned));
		w_Scores_tot = smalloc(templates->DB_size * sizeof(long unsigned));
		
		for(i = 0; i < templates->DB_size; ++i) {
			w_Scores[i] = Scores[i];
			w_Scores_tot[i] = Scores_tot[i];
			if(Scores[i] == 0) {
				SearchList[i] = 0;
			} else {
				SearchList[i] = 1;
			}
		}
		w_Nhits.n = Nhits.n;
		w_Nhits.tot = Nhits.tot;
		
		if(w_Scores[contamination] > 0 || w_Scores_tot[contamination] > 0) {
			fprintf(stderr, " Failed at removing contamination\n");
			exit(1);
		}
		SearchList[contamination] = 0;
		
		/* get best matches */
		while(! stop) {
			/* get best match, depth first then coverage */
			depth = 0;
			cover = 0;
			score = 0;
			template = 0;
			if(ss == 'q') {
				for(i = 0; i < templates->DB_size; ++i) {
					if(SearchList[i] && w_Scores_tot[i] >= score) {
						tmp_cover = 100.0 * w_Scores[i] / template_ulengths[i];
						tmp_score = w_Scores_tot[i];
						tmp_depth = 1.0 * tmp_score / template_lengths[i];
						if(ID_t <= tmp_cover && Depth_t <= tmp_depth) {
							if(tmp_score > score || (tmp_cover > cover || (tmp_cover == cover && (tmp_depth > depth || (tmp_depth == depth && template_ulengths[i] > template_ulengths[template]))))) {
								/* calculate p_value */
								tmp_expected = 1.0 * (Nhits.tot - w_Scores_tot[i]) * template_ulengths[i] / (templates->n - template_ulengths[i] + etta);
								tmp_q = (tmp_score - tmp_expected) * (tmp_score - tmp_expected) / (tmp_score + tmp_expected);
								tmp_p = p_chisqr(tmp_q);
								if(tmp_p <= evalue && tmp_score > tmp_expected) {
									score = tmp_score;
									cover = tmp_cover;
									depth = tmp_depth;
									template = i;
									expected = tmp_expected;
									p_value = tmp_p;
									q_value = tmp_q;
								} else {
									SearchList[i] = 0;
								}
							}
						} else {
							SearchList[i] = 0;
						}
					}
				}
			} else if (ss == 'd') {
				for(i = 0; i < templates->DB_size; ++i) {
					if(SearchList[i]) {
						tmp_cover = 100.0 * w_Scores[i] / template_ulengths[i];
						tmp_score = w_Scores_tot[i];
						tmp_depth = 1.0 * tmp_score / template_lengths[i];
						if(ID_t <= tmp_cover && Depth_t <= tmp_depth) {
							if(tmp_depth > depth || (tmp_depth == depth && (tmp_cover > cover || (tmp_cover == cover && (tmp_score > score || (tmp_score == score && template_ulengths[i] > template_ulengths[template])))))) {
								/* calculate p_value */
								tmp_expected = 1.0 * (Nhits.tot - w_Scores_tot[i]) * template_ulengths[i] / (templates->n - template_ulengths[i] + etta);
								tmp_q = (tmp_score - tmp_expected) * (tmp_score - tmp_expected) / (tmp_score + tmp_expected);
								tmp_p = p_chisqr(tmp_q);
								if(tmp_p <= evalue && tmp_score > tmp_expected) {
									score = tmp_score;
									cover = tmp_cover;
									depth = tmp_depth;
									template = i;
									expected = tmp_expected;
									p_value = tmp_p;
									q_value = tmp_q;
								} else {
									SearchList[i] = 0;
								}
							}
						} else {
							SearchList[i] = 0;
						}
					}
				}
			} else {
				for(i = 0; i < templates->DB_size; ++i) {
					if(SearchList[i]) {
						tmp_cover = 100.0 * w_Scores[i] / template_ulengths[i];
						tmp_score = w_Scores_tot[i];
						tmp_depth = 1.0 * tmp_score / template_lengths[i];
						if(ID_t <= tmp_cover && Depth_t <= tmp_depth) {
							if(tmp_cover > cover || (tmp_cover == cover && (tmp_depth > depth || (tmp_depth == depth && (tmp_score > score || (tmp_score == score && template_ulengths[i] > template_ulengths[template])))))) {
								/* calculate p_value */
								tmp_expected = 1.0 * (Nhits.tot - w_Scores_tot[i]) * template_ulengths[i] / (templates->n - template_ulengths[i] + etta);
								tmp_q = (tmp_score - tmp_expected) * (tmp_score - tmp_expected) / (tmp_score + tmp_expected);
								tmp_p = p_chisqr(tmp_q);
								if(tmp_p <= evalue && tmp_score > tmp_expected) {
									score = tmp_score;
									cover = tmp_cover;
									depth = tmp_depth;
									template = i;
									expected = tmp_expected;
									p_value = tmp_p;
									q_value = tmp_q;
								} else {
									SearchList[i] = 0;
								}
							}
						} else {
							SearchList[i] = 0;
						}
					}
				}
			}
			
			/* validate best match */
			if(cover && ID_t <= cover && Depth_t <= depth) {
				/* with draw contamination k-mers matching this template */
				score_add = 0;
				score_tot_add = 0;
				if(deConTable != 0) {
					prev = 0;
					node = deConTable;
					while(node != 0) {
						if(intpos_bin(node->value, template) != -1) {
							++score_add;
							score_tot_add += node->key;
							if(prev == 0) {
								deConTable = node->next;
								free(node->value);
								free(node);
								node = deConTable;
							} else {
								prev->next = node->next;
								free(node->value);
								free(node);
								node = prev->next;
							}
						} else {
							prev = node;
							node = node->next;
						}
					}
				}
				
				/* Calculate new attributes */
				query_cover = 100.0 * (w_Scores_tot[template] + score_tot_add) / Ntot;
				cover = 100.0 * (w_Scores[template] + score_add) / template_ulengths[template];
				depth = 1.0 * (w_Scores_tot[template] + score_tot_add) / template_lengths[template];
				/* calc tot values */
				tot_cover = 100.0 * (Scores[template] + score_add) / template_ulengths[template];
				tot_depth = 1.0 * (Scores_tot[template] + score_tot_add) / template_lengths[template];
				tot_query_cover = 100.0 * (Scores_tot[template] + score_tot_add) / Ntot;
				
				/* output results */
				fprintf(sparse_out, "%s\t%d\t%lu\t%d\t%d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%f\t%e\n", 
					template_names[template], template, score, (int) expected, template_lengths[template], query_cover, cover, depth, tot_query_cover, tot_cover, tot_depth, q_value, p_value);
				
				/* update scores */
				kmerList = withDraw_Kmers(w_Scores, w_Scores_tot, kmerList, template, &w_Nhits);
				
				if(w_Scores[template] != 0 || w_Scores_tot[template] != 0) {
					fprintf(stderr, "# Failed updating the scores\n");
					SearchList[template] = 0;
				} else {
					SearchList[template] = 0;
				}
				if(kmerList == 0) {
					stop = 1;
				}
			} else {
				stop = 1;
			}
		}
	} else {
		/* collect scores */
		kmerList = collect_Kmers(templates, Scores, Scores_tot, foundKmers, &Nhits);
		
		fprintf(stderr, "# Total number of matches: %lu of %lu kmers\n", Nhits.tot, Ntot);
		/* copy scores */
		w_Scores = smalloc(templates->DB_size * sizeof(int));
		w_Scores_tot = smalloc(templates->DB_size * sizeof(long unsigned));
		
		for(i = 0; i < templates->DB_size; ++i) {
			w_Scores[i] = Scores[i];
			w_Scores_tot[i] = Scores_tot[i];
			if(Scores[i] == 0) {
				SearchList[i] = 0;
			} else {
				SearchList[i] = 1;
			}
		}
		w_Nhits.n = Nhits.n;
		w_Nhits.tot = Nhits.tot;
		
		if(kmerList == 0) {
			stop = 1;
		}
		
		while(!stop) {
			/* get best match, depth first then coverage */
			depth = 0;
			cover = 0;
			score = 0;
			template = 0;
			if(ss == 'q') {
				for(i = 0; i < templates->DB_size; ++i) {
					if(SearchList[i] && w_Scores_tot[i] >= score) {
						tmp_cover = 100.0 * w_Scores[i] / template_ulengths[i];
						tmp_score = w_Scores_tot[i];
						tmp_depth = 1.0 * tmp_score / template_lengths[i];
						if(ID_t <= tmp_cover && Depth_t <= tmp_depth) {
							if(tmp_score > score || (tmp_cover > cover || (tmp_cover == cover && (tmp_depth > depth || (tmp_depth == depth && template_ulengths[i] > template_ulengths[template]))))) {
								/* calculate p_value */
								tmp_expected = 1.0 * (Nhits.tot - w_Scores_tot[i]) * template_ulengths[i] / (templates->n - template_ulengths[i] + etta);
								tmp_q = (tmp_score - tmp_expected) * (tmp_score - tmp_expected) / (tmp_score + tmp_expected);
								tmp_p = p_chisqr(tmp_q);
								if(tmp_p <= evalue && tmp_score > tmp_expected) {
									score = tmp_score;
									cover = tmp_cover;
									depth = tmp_depth;
									template = i;
									expected = tmp_expected;
									p_value = tmp_p;
									q_value = tmp_q;
								} else {
									SearchList[i] = 0;
								}
							}
						} else {
							SearchList[i] = 0;
						}
					}
				}
			} else if (ss == 'd') {
				for(i = 0; i < templates->DB_size; ++i) {
					if(SearchList[i]) {
						tmp_cover = 100.0 * w_Scores[i] / template_ulengths[i];
						tmp_score = w_Scores_tot[i];
						tmp_depth = 1.0 * tmp_score / template_lengths[i];
						if(ID_t <= tmp_cover && Depth_t <= tmp_depth) {
							if(tmp_depth > depth || (tmp_depth == depth && (tmp_cover > cover || (tmp_cover == cover && (tmp_score > score || (tmp_score == score && template_ulengths[i] > template_ulengths[template])))))) {
								/* calculate p_value */
								tmp_expected = 1.0 * (Nhits.tot - w_Scores_tot[i]) * template_ulengths[i] / (templates->n - template_ulengths[i] + etta);
								tmp_q = (tmp_score - tmp_expected) * (tmp_score - tmp_expected) / (tmp_score + tmp_expected);
								tmp_p = p_chisqr(tmp_q);
								if(tmp_p <= evalue && tmp_score > tmp_expected) {
									score = tmp_score;
									cover = tmp_cover;
									depth = tmp_depth;
									template = i;
									expected = tmp_expected;
									p_value = tmp_p;
									q_value = tmp_q;
								} else {
									SearchList[i] = 0;
								}
							}
						} else {
							SearchList[i] = 0;
						}
					}
				}
			} else {
				for(i = 0; i < templates->DB_size; ++i) {
					if(SearchList[i]) {
						tmp_cover = 100.0 * w_Scores[i] / template_ulengths[i];
						tmp_score = w_Scores_tot[i];
						tmp_depth = 1.0 * tmp_score / template_lengths[i];
						if(ID_t <= tmp_cover && Depth_t <= tmp_depth) {
							if(tmp_cover > cover || (tmp_cover == cover && (tmp_depth > depth || (tmp_depth == depth && (tmp_score > score || (tmp_score == score && template_ulengths[i] > template_ulengths[template])))))) {
								/* calculate p_value */
								tmp_expected = 1.0 * (Nhits.tot - w_Scores_tot[i]) * template_ulengths[i] / (templates->n - template_ulengths[i] + etta);
								tmp_q = (tmp_score - tmp_expected) * (tmp_score - tmp_expected) / (tmp_score + tmp_expected);
								tmp_p = p_chisqr(tmp_q);
								if(tmp_p <= evalue && tmp_score > tmp_expected) {
									score = tmp_score;
									cover = tmp_cover;
									depth = tmp_depth;
									template = i;
									expected = tmp_expected;
									p_value = tmp_p;
									q_value = tmp_q;
								} else {
									SearchList[i] = 0;
								}
							}
						} else {
							SearchList[i] = 0;
						}
					}
				}
			}
			
			/* validate best match */
			if(cover && ID_t <= cover && Depth_t <= depth) {
				/* output results */
				query_cover = 100.0 * w_Scores_tot[template] / Ntot;
				/* calc tot values */
				tot_cover = 100.0 * Scores[template] / template_ulengths[template];
				tot_depth = 1.0 * Scores_tot[template] / template_lengths[template];
				tot_query_cover = 100.0 * Scores_tot[template] / Ntot;
				fprintf(sparse_out, "%s\t%d\t%lu\t%d\t%d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n", 
					template_names[template], template, score, (int) expected, template_ulengths[template], query_cover, cover, depth, tot_query_cover, tot_cover, tot_depth, q_value, p_value);
				
				/* update scores */
				kmerList = withDraw_Kmers(w_Scores, w_Scores_tot, kmerList, template, &w_Nhits);
				if(w_Scores[template] != 0 || w_Scores_tot[template] != 0) {
					fprintf(stderr, "# Failed updating the scores\n");
					SearchList[template] = 0;
				} else {
					SearchList[template] = 0;
				}
				if(kmerList == 0) {
					stop = 1;
				}
			} else {
				stop = 1;
			}
		}
	}
	fclose(sparse_out);
	t1 = clock();
	fprintf(stderr, "# Total for finding and outputting best matches: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	
	return status;
}
