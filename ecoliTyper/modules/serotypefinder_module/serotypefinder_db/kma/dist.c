/* Philip T.L.C. Clausen Jan 2020 plan@dtu.dk */

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
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/param.h>
#include "dist.h"
#include "hashmapkma.h"
#include "matrix.h"
#include "pherror.h"
#include "runkma.h"
#include "threader.h"
#include "tmp.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.\n", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

HashMapKMA * loadValues(const char *filename) {
	
	long unsigned check, size;
	FILE *file;
	HashMapKMA *dest;
	
	/* init */
	file = sfopen(filename, "rb");
	dest = smalloc(sizeof(HashMapKMA));
	
	/* load sizes */
	sfread(&dest->DB_size, sizeof(unsigned), 1, file);
	sfread(&dest->mlen, sizeof(unsigned), 1, file);
	sfread(&dest->prefix_len, sizeof(unsigned), 1, file);
	sfread(&dest->prefix, sizeof(long unsigned), 1, file);
	sfread(&dest->size, sizeof(long unsigned), 1, file);
	sfread(&dest->n, sizeof(long unsigned), 1, file);
	sfread(&dest->v_index, sizeof(long unsigned), 1, file);
	sfread(&dest->null_index, sizeof(long unsigned), 1, file);
	
	dest->mask = 0;
	dest->mask = (~dest->mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (dest->mlen << 1));
	dest->shmFlag = 0;
	
	/* simple check for old indexing */
	if(dest->size < dest->n) {
		free(dest);
		fclose(file);
		return 0;
	}
	
	/* exist */
	size = dest->size;
	if((dest->size - 1) == dest->mask) {
		/* mega */
		if(dest->v_index <= UINT_MAX) {
			size *= sizeof(unsigned);
		} else {
			size *= sizeof(long unsigned);
		}
		
		/* load */
		dest->exist = smalloc(size);
		check = fread(dest->exist, 1, size, file);
		if(check != size) {
			free(dest->exist);
			free(dest);
			fclose(file);
			return 0;
		}
		dest->exist_l = (long unsigned *)(dest->exist);
		dest->shmFlag |= 1;
	} else {
		if(dest->n <= UINT_MAX) {
			size *= sizeof(unsigned);
		} else {
			size *= sizeof(long unsigned);
		}
		
		/* skip */
		dest->exist = 0;
		dest->exist_l = 0;
		sfseek(file, size, SEEK_CUR);
	}
	
	/* values */
	size = dest->v_index;
	if(dest->DB_size < USHRT_MAX) {
		size *= sizeof(short unsigned);
	} else {
		size *= sizeof(unsigned);
	}
	dest->values = smalloc(size);
	check = fread(dest->values, 1, size, file);
	if(check != size) {
		free(dest->exist);
		free(dest->values);
		free(dest);
		fclose(file);
		return 0;
	}
	dest->values_s = (short unsigned *)(dest->values);
	
	/* check for megaMap */
	if(!dest->exist) {
		/* skip kmers */
		size = dest->n + 1;
		if(dest->mlen <= 16) {
			size *= sizeof(unsigned);
		} else {
			size *= sizeof(long unsigned);
		}
		dest->key_index = 0;
		dest->key_index_l = 0;
		sfseek(file, size, SEEK_CUR);
		
		/* value indexes */
		size = dest->n;
		if(dest->v_index < UINT_MAX) {
			size *= sizeof(unsigned);
		} else {
			size *= sizeof(long unsigned);
		}
		dest->exist = smalloc(size);
		check = fread(dest->exist, 1, size, file);
		if(check != size) {
			free(dest->exist);
			free(dest->values);
			free(dest);
			fclose(file);
			return 0;
		}
		dest->exist_l = (long unsigned *)(dest->exist);
	}
	
	if(fread(&dest->kmersize, sizeof(unsigned), 1, file)) {
		sfread(&dest->flag, sizeof(unsigned), 1, file);
	} else {
		dest->kmersize = dest->mlen;
		dest->flag = 0;
	}
	
	fclose(file);
	return dest;
}

void destroyValues(HashMapKMA *src) {
	
	free(src->exist);
	free(src->values);
	free(src);
}

void kmerSimilarity(HashMapKMA *DB, Matrix *Dist, int *N) {
	
	int i, j, el, vs, v_i, **D, *Di;
	unsigned *exist, *values_i, *values_j;
	long unsigned n, pos, *exist_l;
	short unsigned *values_si, *values_sj;
	
	/* init */
	el = DB->v_index < UINT_MAX;
	vs = DB->DB_size < USHRT_MAX;
	exist = DB->exist - 1;
	exist_l = DB->exist_l - 1;
	D = Dist->mat;
	
	/* get values */
	n = DB->n;
	while(n) {
		pos = el ? *++exist : *++exist_l;
		if(pos != 1) {
			if(vs) {
				values_si = DB->values_s + pos;
				i = *values_si + 1;
				values_si += i;
				while(--i) {
					j = i;
					values_sj = --values_si;
					v_i = *values_si - 1;
					Di = D[v_i];
					while(--j) {
						++Di[*--values_sj - 1];
					}
					++N[v_i];
				}
			} else {
				values_i = DB->values + pos;
				i = *values_i + 1;
				values_i += i;
				while(--i) {
					j = i;
					values_j = --values_i;
					v_i = *values_i - 1;
					Di = D[v_i];
					while(--j) {
						++Di[*--values_j - 1];
					}
					++N[v_i];
				}
			}
			--n;
		}
	}
	
	Dist->n = DB->DB_size - 1;
	destroyValues(DB);
}

void kmerSimilarity_thread(HashMapKMA *DB, Matrix *Dist, int *N, int thread_num, volatile int *lock) {
	
	static volatile int thread_wait = 0;
	static volatile long unsigned next;
	static long unsigned n;
	int i, j, el, vs, v_i, chunk, chunkster, **D, *Di, *Dptr;
	unsigned *exist, *values_i, *values_j;
	long unsigned pos, size, *exist_l;
	short unsigned *values_si, *values_sj;
	
	/* init */
	el = DB->v_index < UINT_MAX;
	vs = DB->DB_size < USHRT_MAX;
	D = Dist->mat;
	chunk = 131072;
	size = ((DB->size - 1) == DB->mask) ? DB->size : DB->n;
	
	/* static init */
	lockTime(lock, 10);
	if(!thread_wait) {
		thread_wait = thread_num;
		next = 0;
		n = DB->n;
		Dist->n = DB->DB_size - 1;
	}
	unlock(lock);
	
	/* get values */
	while(n && chunk) {
		/* get next chunk */
		lockTime(lock, 10);
		pos = next;
		if(size < (next += chunk)) {
			next = size;
		}
		unlock(lock);
		
		/* update chunk */
		exist = DB->exist + pos - 1;
		exist_l = DB->exist_l + pos - 1;
		if(size < pos + chunk) {
			chunkster = size - pos + 1;
			chunk = 0;
		} else {
			chunkster = (chunk + 1);
		}
		while(n && --chunkster) {
			pos = el ? *++exist : *++exist_l;
			if(pos != 1) {
				if(vs) {
					values_si = DB->values_s + pos;
					i = *values_si + 1;
					values_si += i;
					while(--i) {
						j = i;
						values_sj = --values_si;
						v_i = *values_si - 1;
						Di = D[v_i];
						while(--j) {
							Dptr = Di + (*--values_sj - 1);
							__sync_add_and_fetch(Dptr, 1);
						}
						Dptr = N + v_i;
						__sync_add_and_fetch(Dptr, 1);
					}
				} else {
					values_i = DB->values + pos;
					i = *values_i + 1;
					values_i += i;
					while(--i) {
						j = i;
						values_j = --values_i;
						v_i = *values_i - 1;
						Di = D[v_i];
						while(--j) {
							Dptr = Di + (*--values_j - 1);
							__sync_add_and_fetch(Dptr, 1);
						}
						Dptr = N + v_i;
						__sync_add_and_fetch(Dptr, 1);
					}
				}
				__sync_sub_and_fetch(&n, 1);
			}
		}
	}
	
	if(__sync_sub_and_fetch(&thread_wait, 1)) {
		wait_atomic(thread_wait);
	} else {
		destroyValues(DB);
	}
}

int kmerDist(int Ni, int Nj, int D) {
	D = Ni + Nj - (D << 1);
	return (D < 0) ? 0 : D;
}

int kmerShared(int Ni, int Nj, int D) {
	return D;
}

int chi2dist(int Ni, int Nj, int D) {
	long d = (Ni + Nj - (D << 1));
	d = d * d / (Ni + Nj);
	return (d < 0) ? 0 : d;
}

void printIntLtdPhy(char *outfile, Matrix *Dist, int *N, FILE *name_file, Qseqs *template_name, unsigned format, int thread_num, volatile int *lock, const char *method, int (*distPtr)(int, int, int)) {
	
	static volatile int thread_wait1 = 0, thread_wait2 = 0;
	static int next_i, next_j, entrance = 0;
	static long unsigned row_bias = 0;
	volatile int *thread_wait;
	int i, j, j_end, d, chunk, N_i, *Nj, **D, *Di;
	char *outfile_chunk, *name;
	
	/* init */
	D = Dist->mat;
	chunk = 65536;
	lockTime(lock, 10);
	if(!row_bias) {
		if(format & 4) {
			row_bias += sprintf(outfile, "# %-35s\n", method);
		}
		row_bias += sprintf(outfile + row_bias, "%10d", Dist->n);
		next_i = -1;
		next_j = 0;
		if(++entrance & 1) {
			thread_wait1 = thread_num;
			thread_wait = &thread_wait1;
		} else {
			thread_wait2 = thread_num;
			thread_wait = &thread_wait2;
		}
	} else if(entrance & 1) {
		thread_wait = &thread_wait1;
	} else {
		thread_wait = &thread_wait2;
	}
	unlock(lock);
	
	while(next_i < Dist->n) {
		lockTime(lock, 10);
		/* check for new row */
		if(Dist->n <= next_i) {
			chunk = 0;
		} else if(next_i <= next_j) {
			i = ++next_i;
			if(i) {
				row_bias += (i - 1) * 11;
			}
			j = 0;
			next_j = chunk;
			if(next_i < Dist->n) {
				name = nameLoad(template_name, name_file);
				if(format & 1) {
					row_bias += sprintf(outfile + row_bias, "\n%s", name);
				} else {
					row_bias += sprintf(outfile + row_bias, "\n%-10.10s", name);
				}
			} else {
				outfile[row_bias] = '\n';
				++row_bias;
				chunk = 0;
			}
		} else {
			i = next_i;
			j = next_j;
			next_j += chunk;
		}
		outfile_chunk = outfile + row_bias + j * 11;
		unlock(lock);
		
		if(chunk) {
			N_i = N[i];
			j_end = (i < j + chunk) ? i : (j + chunk);
			Di = D[i] + --j;
			Nj = N + j;
			while(++j < j_end) {
				d = distPtr(N_i, *++Nj, *++Di);
				outfile_chunk += sprintf(outfile_chunk, "\t%10d", (d < 0 ? 0 : d));
			}
			/* fix null character */
			*outfile_chunk = (j == i) ? '\n' : '\t';
		}
	}
	
	/* wait for matrix to finish */
	lockTime(lock, 10);
	if(--*thread_wait) {
		unlock(lock);
		wait_atomic(*thread_wait);
	} else {
		sfseek(name_file, 0, SEEK_SET);
		row_bias = 0;
		unlock(lock);
	}
}

double kmerQuery(int Ni, int Nj, int D) {
	double d = 100.0 * D / Ni;
	return (d < 0) ? 0.0 : ((100 < d) ? 100.0 : d);
}

double kmerTemplate(int Ni, int Nj, int D) {
	double d = 100.0 * D / Nj;
	return (d < 0) ? 0.0 : ((100 < d) ? 100.0 : d);
}

double kmerAvg(int Ni, int Nj, int D) {
	double d = 200.0 * D / (Ni + Nj);
	return (d < 0) ? 0.0 : ((100 < d) ? 100.0 : d);
}

double kmerInvAvg(int Ni, int Nj, int D) {
	double d = 100.0 - 200.0 * D / (Ni + Nj);
	return (d < 0) ? 0.0 : ((100 < d) ? 100.0 : d);
}

double kmerJaccardDist(int Ni, int Nj, int D) {
	double d = 1.0 - (double)(D) / (Ni + Nj - D);
	return (d < 0) ? 0.0 : ((1 < d) ? 1.0 : d);
}

double kmerJaccardSim(int Ni, int Nj, int D) {
	double d = (double)(D) / (Ni + Nj - D);
	return (d < 0) ? 0.0 : ((1 < d) ? 1.0 : d);
}

double kmerCosineDist(int Ni, int Nj, int D) {
	double d = 1.0 - (double)(D) / (sqrt(Ni) * sqrt(Nj));
	return (d < 0) ? 0.0 : ((1 < d) ? 1.0 : d);
}

double kmerCosineSim(int Ni, int Nj, int D) {
	double d = (double)(D) / (sqrt(Ni) * sqrt(Nj));
	return (d < 0) ? 0.0 : ((1 < d) ? 1.0 : d);
}

double kmerOverlapCoef(int Ni, int Nj, int D) {
	double d = (double)(D) / (Ni < Nj ? Ni : Nj);
	return (d < 0) ? 0.0 : ((1 < d) ? 1.0 : d);
}

double kmerInvOverlapCoef(int Ni, int Nj, int D) {
	double d = 1.0 - (double)(D) / (Ni < Nj ? Ni : Nj);
	return (d < 0) ? 0.0 : ((1 < d) ? 1.0 : d);
}

void printDoublePhy(char *outfile, Matrix *Dist, int *N, FILE *name_file, Qseqs *template_name, unsigned format, const char *formatString, int ltd, int thread_num, volatile int *lock, const char *method, double (*distPtr)(int, int, int)) {
	
	static volatile int thread_wait1 = 0, thread_wait2 = 0;
	static int next_i, next_j, entrance = 0;
	static long unsigned row_bias = 0;
	volatile int *thread_wait;
	int i, j, j_end, chunk, N_i, *Nj, **D, *Di;
	char endChar, *outfile_chunk, *name;
	double d;
	
	/* init */
	D = Dist->mat;
	chunk = 65536;
	lockTime(lock, 10);
	if(!row_bias) {
		if(format & 4) {
			row_bias += sprintf(outfile, "# %-35s\n", method);
		}
		row_bias += sprintf(outfile + row_bias, "%10d", Dist->n);
		next_i = -1;
		next_j = Dist->n;
		if(++entrance & 1) {
			thread_wait1 = thread_num;
			thread_wait = &thread_wait1;
		} else {
			thread_wait2 = thread_num;
			thread_wait = &thread_wait2;
		}
	} else if(entrance & 1) {
		thread_wait = &thread_wait1;
	} else {
		thread_wait = &thread_wait2;
	}
	unlock(lock);
	
	while(next_i < Dist->n) {
		lockTime(lock, 10);
		/* check for new row */
		if(Dist->n <= next_i) {
			chunk = 0;
		} else if((ltd && next_i <= next_j) || (Dist->n <= next_j)) {
			i = ++next_i;
			if(i) {
				row_bias += (ltd ? (i - 1) : Dist->n) * 11;
			}
			j = 0;
			next_j = chunk;
			if(next_i < Dist->n) {
				name = nameLoad(template_name, name_file);
				if(format & 1) {
					row_bias += sprintf(outfile + row_bias, "\n%s", name);
				} else {
					row_bias += sprintf(outfile + row_bias, "\n%-10.10s", name);
				}
			} else {
				outfile[row_bias] = '\n';
				++row_bias;
				chunk = 0;
			}
		} else {
			i = next_i;
			j = next_j;
			next_j += chunk;
		}
		outfile_chunk = outfile + row_bias + j * 11;
		unlock(lock);
		
		if(chunk) {
			if(ltd) {
				j_end = (i < j + chunk) ? i : (j + chunk);
				endChar = (j_end == i) ? '\n' : '\t';
			} else {
				j_end = (Dist->n < j + chunk) ? Dist->n : (j + chunk);
				endChar = (j_end == Dist->n) ? '\n' : '\t';
			}
			N_i = N[i];
			Di = D[i] + --j;
			Nj = N + j;
			
			if(j_end < i) {
				while(++j < j_end) {
					d = distPtr(N_i, *++Nj, *++Di);
					outfile_chunk += sprintf(outfile_chunk, formatString, (d < 0 ? 0 : d));
				}
			} else if(i < j) {
				while(++j < j_end) {
					d = distPtr(N_i, *++Nj, D[j][i]);
					outfile_chunk += sprintf(outfile_chunk, formatString, (d < 0 ? 0 : d));
				}
			} else {
				while(++j < j_end) {
					if(j < i) {
						d = distPtr(N_i, *++Nj, *++Di);
						outfile_chunk += sprintf(outfile_chunk, formatString, (d < 0 ? 0 : d));
					} else if(i < j) {
						d = distPtr(N_i, *++Nj, D[j][i]);
						outfile_chunk += sprintf(outfile_chunk, formatString, (d < 0 ? 0 : d));
					} else {
						outfile_chunk += sprintf(outfile_chunk, formatString, 100.0);
						++Nj;
					}
				}
			}
			/* fix null character */
			*outfile_chunk = endChar;
		}
	}
	
	/* wait for matrix to finish */
	lockTime(lock, 10);
	if(--*thread_wait) {
		unlock(lock);
		wait_atomic(*thread_wait);
	} else {
		sfseek(name_file, 0, SEEK_SET);
		row_bias = 0;
		unlock(lock);
	}
}

long unsigned getPhySize(int flag, int format, long unsigned n, long unsigned *ltdMat, long unsigned *covMat, FILE *name_file) {
	
	long unsigned size;
	
	/* get name name */
	if(format & 1) {
		sfseek(name_file, 0, SEEK_END);
		size = ftell(name_file);
		rewind(name_file);
	} else {
		size = n * 11;
	}
	if(format & 4) {
		size += 38;
	}
	
	/* phy size */
	size += 11;
	
	/* ltd cell size */
	*ltdMat = size + (((n - 1) * (n - 2)) >> 1) * 11;
	
	/* cov cell size */
	*covMat = size + (n - 1) * (n - 1) * 11;
	
	/* get number of matrices to make */
	size = 0;
	if(flag & 4) {
		size += *covMat;
		flag ^= 4;
	}
	if(flag & 8) {
		size += *covMat;
		flag ^= 8;
	}
	n = 0;
	while(flag) {
		n += flag & 1;
		flag >>= 1;
	}
	
	return size + n * *ltdMat;
}

char * mfile(FILE *outfile, long unsigned size) {
	
	char *outfileM;
	
	if(fseek(outfile, size - 1, SEEK_SET) || putc(0, outfile) == EOF) {
		ERROR();
	}
	rewind(outfile);
	outfileM = mmap(0, size, PROT_WRITE, MAP_SHARED, fileno(outfile), 0);
	if(outfileM == MAP_FAILED) {
		ERROR();
	}
	posix_madvise(outfileM, size, POSIX_MADV_SEQUENTIAL);
	
	return outfileM;
}

void * threadDist(void *arg) {
	
	static volatile int Lock = 0;
	volatile int *lock = &Lock;
	DistThread *thread = arg;
	int flag, format, thread_num, *N;
	long unsigned ltdMat, covMat;
	char *outfileM;
	FILE *name_file;
	HashMapKMA *DB;
	Matrix *Dist;
	Qseqs *template_name;
	
	/* init */
	flag = thread->flag;
	format = thread->format;
	thread_num = thread->thread_num;
	N = thread->N;
	ltdMat = thread->ltdMat;
	covMat = thread->covMat;
	outfileM = thread->outfileM;
	name_file = thread->name_file;
	DB = thread->DB;
	Dist = thread->Dist;
	template_name = thread->template_name;
	
	/* get kmer similarities and lengths */
	if(thread_num != 1) {
		kmerSimilarity_thread(DB, Dist, N, thread_num, lock);
	} else {
		kmerSimilarity(DB, Dist, N);
	}
	
	/* k-mer dist, lt */
	if(flag & 1) {
		printIntLtdPhy(outfileM, Dist, N, name_file, template_name, format, thread_num, lock, "k-mer distance", &kmerDist);
		outfileM += ltdMat;
	}
	
	/* k-mer shared, lt */
	if(flag & 2) {
		printIntLtdPhy(outfileM, Dist, N, name_file, template_name, format, thread_num, lock, "shared k-mers", &kmerShared);
		outfileM += ltdMat;
	}
	
	/* query cov, asym */
	if(flag & 4) {
		printDoublePhy(outfileM, Dist, N, name_file, template_name, format, "\t%10.6f", 0, thread_num, lock, "Query k-mer coverage [%]", &kmerQuery);
		outfileM += covMat;
	}
	
	/* template cov, asym */
	if(flag & 8) {
		printDoublePhy(outfileM, Dist, N, name_file, template_name, format, "\t%10.6f", 0, thread_num, lock, "Template k-mer coverage [%]", &kmerTemplate);
		outfileM += covMat;
	}
	
	/* avg. cov, lt */
	if(flag & 16) {
		printDoublePhy(outfileM, Dist, N, name_file, template_name, format, "\t%10.6f", 1, thread_num, lock, "Avg. k-mer coverage [%]", &kmerAvg);
		outfileM += ltdMat;
	}
	
	/* inv. avg. cov, lt */
	if(flag & 32) {
		printDoublePhy(outfileM, Dist, N, name_file, template_name, format, "\t%10.6f", 1, thread_num, lock, "Inverse Avg. k-mer coverage", &kmerInvAvg);
		outfileM += ltdMat;
	}
	
	/* Jaccard dist */
	if(flag & 64) {
		printDoublePhy(outfileM, Dist, N, name_file, template_name, format, "\t%.8f", 1, thread_num, lock, "Jaccard Distance", &kmerJaccardDist);
		outfileM += ltdMat;
	}
	
	/* Jaccard similarity */
	if(flag & 128) {
		printDoublePhy(outfileM, Dist, N, name_file, template_name, format, "\t%.8f", 1, thread_num, lock, "Jaccard Similarity", &kmerJaccardSim);
		outfileM += ltdMat;
	}
	
	/* Cosine dist */
	if(flag & 256) {
		printDoublePhy(outfileM, Dist, N, name_file, template_name, format, "\t%.8f", 1, thread_num, lock, "Cosine distance", &kmerCosineDist);
		outfileM += ltdMat;
	}
	
	/* Cosine similarity */
	if(flag & 512) {
		printDoublePhy(outfileM, Dist, N, name_file, template_name, format, "\t%.8f", 1, thread_num, lock, "Cosine similarity", &kmerCosineSim);
		outfileM += ltdMat;
	}
	
	/* Szymkiewicz–Simpson similarity */
	if(flag & 1024) {
		printDoublePhy(outfileM, Dist, N, name_file, template_name, format, "\t%.8f", 1, thread_num, lock, "Szymkiewicz–Simpson similarity", &kmerOverlapCoef);
		outfileM += ltdMat;
	}
	
	/* Szymkiewicz–Simpson dissimilarity */
	if(flag & 2048) {
		printDoublePhy(outfileM, Dist, N, name_file, template_name, format, "\t%.8f", 1, thread_num, lock, "Szymkiewicz–Simpson dissimilarity", &kmerInvOverlapCoef);
		outfileM += ltdMat;
	}
	
	/* Chi-square distance */
	if(flag & 4096) {
		printIntLtdPhy(outfileM, Dist, N, name_file, template_name, format, thread_num, lock, "Chi-square distance", &chi2dist);
		outfileM += ltdMat;
	}
	
	return NULL;
}

void runDist(char *templatefilename, char *outputfilename, int flag, int format, int disk, int thread_num) {
	
	int i, file_len, *N;
	unsigned DB_size;
	long unsigned out_size, ltdMat, covMat;
	char *outfileM;
	FILE *outfile, *name_file;
	DistThread *thread, *threads;
	HashMapKMA *DB;
	Matrix *Dist;
	Qseqs *template_name;
	
	/* init */
	outfile = sfopen(outputfilename, "wb+");
	file_len = strlen(templatefilename);
	
	/* load k-mer links from KMA DB */
	strcpy(templatefilename + file_len, ".comp.b");
	if(!(DB = loadValues(templatefilename))) {
		fprintf(stderr, "Wrong format of DB.\n");
		exit(1);
	}
	templatefilename[file_len] = 0;
	
	/* load names */
	strcpy(templatefilename + file_len, ".name");
	name_file = sfopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	template_name = setQseqs(1024);
	
	/* allocate output file */
	DB_size = DB->DB_size;
	out_size = getPhySize(flag, format, DB_size, &ltdMat, &covMat, name_file);
	outfileM = mfile(outfile, out_size);
	
	/* allocate matrices */
	if(!(N = calloc(DB_size, sizeof(int)))) {
		ERROR();
	}
	if(disk) {
		Dist = ltdMatrix_minit(DB_size);
	} else {
		Dist = ltdMatrix_init(DB_size);
	}
	
	/* thread out */
	thread = 0;
	threads = 0;
	i = thread_num;
	while(i--) {
		/* init */
		thread = smalloc(sizeof(DistThread));
		thread->flag = flag;
		thread->format = format;
		thread->thread_num = thread_num;
		thread->N = N;
		thread->ltdMat = ltdMat;
		thread->covMat = covMat;
		thread->outfileM = outfileM;
		thread->name_file = name_file;
		thread->DB = DB;
		thread->Dist = Dist;
		thread->template_name = template_name;
		thread->next = threads;
		threads = thread;
		/* start */
		if(i && (errno = pthread_create(&thread->id, NULL, &threadDist, thread))) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			fprintf(stderr, "Will continue with %d threads.\n", thread_num - i);
			i = 0;
		}
	}
	/* start main thread */
	thread->id = 0;
	threadDist(thread);
	
	/* join threads */
	while((thread = thread->next)) {
		/* join thread */
		if((errno = pthread_join(thread->id, NULL))) {
				ERROR();
		}
	}
	
	/* clean */
	while(threads) {
		thread = threads->next;
		free(threads);
		threads = thread;
	}
	outfileM -= out_size;
	msync(outfileM, out_size, MS_SYNC);
	munmap(outfileM, out_size);
	fclose(outfile);
	fclose(name_file);
	free(N);
	if(disk) {
		Matrix_mdestroy(Dist);
	} else {
		Matrix_destroy(Dist);
	}
	destroyQseqs(template_name);
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#kma dist calculates distances between templates from a kma index\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-t_db", "Template DB", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "DB");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-f", "Output flags", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fh", "Help on option \"-f\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-d", "Distance method", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-dh", "Help on option \"-d\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-m", "Allocate matrix on the disk", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-tmp", "Set directory for temporary files", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-t", "Number of threads", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

/* main */
int dist_main(int argc, char *argv[]) {
	
	int args, flag, format, mmap, thread_num, file_len;
	char *arg, *errorMsg, *templatefilename, *outputfilename;
	
	/* init */
	flag = 1;
	format = 1;
	mmap = 0;
	thread_num = 1;
	file_len = 0;
	templatefilename = 0;
	outputfilename = 0;
	
	/* parse cmd-line */
	args = 1;
	while(args < argc) {
		arg = argv[args];
		if(*arg++ == '-') {
			 if(strcmp(arg, "t_db") == 0) {
				if(++args < argc) {
					file_len = strlen(argv[args]);
					templatefilename = smalloc(file_len + 64);
					strcpy(templatefilename, argv[args]);
				} else {
					missArg("\"-o\"");
				}
			} else if(strcmp(arg, "o") == 0) {
				if(++args < argc) {
					outputfilename = argv[args];
				} else {
					missArg("\"-o\"");
				}
			} else if(strcmp(arg, "f") == 0) {
				if(++args < argc) {
					format = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-f\"");
					}
				} else {
					missArg("\"-f\"");
				}
			} else if(strcmp(arg, "fh") == 0) {
				fprintf(stdout, "# Format flags output, add them to combine them.\n");
				fprintf(stdout, "#\n");
				fprintf(stdout, "#%9d\t%s\n", 1, "Relaxed Phylip");
				fprintf(stdout, "#%9d\t%s\n", 4, "Include distance method(s) in phylip file");
				return 0;
			} else if(strcmp(arg, "d") == 0) {
				if(++args < argc) {
					flag = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-d\"");
					}
				} else {
					missArg("\"-d\"");
				}
			} else if(strcmp(arg, "dh") == 0) {
				fprintf(stdout, "# Distance / Similarity calculation methods, add them to combine them:\n");
				fprintf(stdout, "#\n");
				fprintf(stdout, "#%9d\t%s\n", 1, "k-mer hamming distance");
				fprintf(stdout, "#%9d\t%s\n", 2, "Shared k-mers");
				fprintf(stdout, "#%9d\t%s\n", 4, "k-mer query coverage");
				fprintf(stdout, "#%9d\t%s\n", 8, "k-mer template coverage");
				fprintf(stdout, "#%9d\t%s\n", 16, "k-mer avg. coverage");
				fprintf(stdout, "#%9d\t%s\n", 32, "k-mer inv. avg. coverage");
				fprintf(stdout, "#%9d\t%s\n", 64, "Jaccard distance");
				fprintf(stdout, "#%9d\t%s\n", 128, "Jaccard similarity");
				fprintf(stdout, "#%9d\t%s\n", 256, "Cosine distance");
				fprintf(stdout, "#%9d\t%s\n", 512, "Cosine similarity");
				fprintf(stdout, "#%9d\t%s\n", 1024, "Szymkiewicz–Simpson similarity");
				fprintf(stdout, "#%9d\t%s\n", 2048, "Szymkiewicz–Simpson dissimilarity");
				fprintf(stdout, "#%9d\t%s\n", 4096, "Chi-square distance");
				fprintf(stdout, "#\n");
				return 0;
			} else if(strcmp(arg, "m") == 0) {
				mmap = 1;
			} else if(strcmp(argv[args], "-tmp") == 0) {
				if(++args < argc) {
					if(argv[args][0] != '-') {
						tmpF(argv[args]);
					} else {
						invaArg("\"-tmp\"");
					}
				}
			} else if(strcmp(argv[args], "-t") == 0) {
				if(++args < argc) {
					thread_num = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-t\"");
					}
				} else {
					missArg("\"-t\"");
				}
			} else if(strcmp(arg, "h") == 0) {
				return helpMessage(stdout);
			} else {
				fprintf(stderr, "Unknown option:%s\n", arg - 1);
				return helpMessage(stderr);
			}
		} else {
			fprintf(stderr, "Unknown argument:%s\n", arg - 1);
			return helpMessage(stderr);
		}
		++args;
	}
	
	if(!templatefilename) {
		fprintf(stderr, "Too few arguments handed.\n");
		exit(1);
	}
	if(!outputfilename) {
		outputfilename = smalloc(file_len + 64);
		file_len = sprintf(outputfilename, "%s.phy", templatefilename);
	}
	
	runDist(templatefilename, outputfilename, flag, format, mmap, thread_num);
	
	return 0;
}
