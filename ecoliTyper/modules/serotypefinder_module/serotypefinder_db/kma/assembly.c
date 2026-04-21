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
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "align.h"
#include "assembly.h"
#include "chain.h"
#include "ef.h"
#include "filebuff.h"
#include "hashmapcci.h"
#include "kmapipe.h"
#include "nw.h"
#include "pherror.h"
#include "qseqs.h"
#include "sam.h"
#include "stdnuc.h"
#include "stdstat.h"
#include "threader.h"
#include "xml.h"
#define mrcheck(mrc, Stat, q_len, t_len) ((mrc * q_len <= Stat.len - Stat.qGaps) || (mrc * t_len <= Stat.len - Stat.tGaps))

void * (*assembly_KMA_Ptr)(void *) = &assemble_KMA;
int (*significantBase)(int, int, double) = &significantNuc;
unsigned char (*baseCall)(unsigned char, unsigned char, int, int, double, Assembly*) = &baseCaller;
void (*alnToMatPtr)(AssemInfo *, Assem *, Aln *, AlnScore, int, int) = &alnToMat;

void updateFrags(FileBuff *dest, Qseqs *qseq, Qseqs *header, char *template_name, int *stats) {
	
	int check, avail;
	char *update;
	
	avail = dest->bytes;
	check = 47 + qseq->len + header->len + strlen(template_name);
	
	/* flush buffer */
	if(avail < check) {
		writeGzFileBuff(dest);
		
		/* seq is too big, reallocate buffer */
		if(dest->bytes < check) {
			resetGzFileBuff(dest, check << 1);
		}
		avail = dest->bytes;
	}
	
	/* update buffer with fragment */
	memcpy(dest->next, qseq->seq, qseq->len);
	dest->next += qseq->len;
	avail -= qseq->len;
	
	/* stats */
	update = (char *) dest->next;
	check = sprintf(update, "\t%d\t%d\t%d\t%d\t%s\t%s\n", stats[0], stats[1], stats[2], stats[3], template_name, header->seq);
	dest->next += check;
	avail -= check;
	dest->bytes = avail;
	
	/* equivalent with:
	fprintf(frag_out, "%s\t%d\t%d\t%d\t%d\t%s\t%s\n", qseq->seq, stats[0], stats[1], stats[2], stats[3], template_names[template], header->seq);
	*/
}

void updateMatrix(FileBuff *dest, char *template_name, long unsigned *template_seq, AssemInfo *matrix, int t_len) {
	
	unsigned i, pos, check, avail, asm_len;
	char *update;
	const char bases[6] = "ACGTN-";
	Assembly *assembly;
	
	/* check buffer capacity */
	check = strlen(template_name) + 2;
	if(dest->bytes < check) {
		writeGzFileBuff(dest);
	}
	update = (char *) dest->next;
	avail = dest->bytes - check;
	
	/* fill in header */
	check -= 2;
	*update++ = '#';
	memcpy(update, template_name, check);
	update += check;
	*update++ = '\n';
	
	/* fill in rows */
	asm_len = matrix->len;
	assembly = matrix->assmb;
	i = 0;
	for(pos = 0; asm_len != 0; --asm_len, pos = assembly[pos].next) {
		/* check buffer capacity */
		if(avail < 38) {
			dest->bytes = avail;
			writeGzFileBuff(dest);
			avail = dest->bytes;
			update = (char *) dest->next;
		}
		
		/* update with row */
		if(pos < t_len) {
			check = sprintf(update, "%c\t%hu\t%hu\t%hu\t%hu\t%hu\t%hu\n", bases[getNuc(template_seq, i)], assembly[pos].counts[0], assembly[pos].counts[1], assembly[pos].counts[2], assembly[pos].counts[3], assembly[pos].counts[4], assembly[pos].counts[5]);
			++i;
		} else {
			check = sprintf(update, "-\t%hu\t%hu\t%hu\t%hu\t%hu\t%hu\n", assembly[pos].counts[0], assembly[pos].counts[1], assembly[pos].counts[2], assembly[pos].counts[3], assembly[pos].counts[4], assembly[pos].counts[5]);
		}
		avail -= check;
		update += check;
	}
	
	/* update with last newline */
	if(avail == 0) {
		writeGzFileBuff(dest);
		avail = dest->bytes;
		update = (char *) dest->next;
	}
	*update++ = '\n';
	dest->next = (unsigned char *) update;
	dest->bytes = avail - 1;
}


int significantNuc(int X, int Y, double evalue) {
	return (Y < X && p_chisqr(pow(X - Y, 2) / (X + Y)) <= evalue);
}

int significantAnd90Nuc(int X, int Y, double evalue) {
	return (Y < X && (9 * (X + Y) <= 10 * X) && p_chisqr(pow(X - Y, 2) / (X + Y)) <= evalue);
}

int significantAndSupport(int X, int Y, double evalue) {
	
	static double support = 0;
	
	if(support == 0) {
		support = evalue;
	}
	
	return (Y < X && (support * (X + Y) <= X) && p_chisqr(pow(X - Y, 2) / (X + Y)) <= evalue);
}

unsigned char baseCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls) {
	
	/* determine base at current position */
	if(depthUpdate == 0) {
		bestNuc = '-';
	} else {
		/* Use MnNemars test to test significance of the base call */
		if(significantBase(bestScore, depthUpdate - bestScore, evalue) == 0) {
			if(bestNuc == '-' && tNuc != '-' && bestScore != depthUpdate) {
				bestNuc = 'n';
			} else {
				bestNuc = tolower(bestNuc);
			}
		}
	}
	
	return bestNuc;
}

unsigned char orgBaseCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls) {
	
	/* determine base at current position */
	if(depthUpdate == 0 || bestNuc == '-') {
		bestNuc = '-';
	} else if(significantBase(bestScore, depthUpdate - bestScore, evalue) == 0) { /* McNemars test */
		bestNuc = tolower(bestNuc);
	}
	
	return bestNuc;
}

unsigned char refCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls) {
	
	/* determine base at current position */
	if(depthUpdate == 0 || (bestNuc == '-' && tNuc != '-')) {
		bestNuc = 'n';
	} else if(significantBase(bestScore, depthUpdate - bestScore, evalue) == 0) {
		bestNuc = tolower(bestNuc);
	}
	
	return bestNuc;
}

unsigned char nanoCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls) {
	
	int j, bestBaseScore;
	const char bases[6] = "ACGTN-";
	
	/* determine base at current position */
	if(depthUpdate == 0) {
		bestNuc = '-';
	} else {
		/* Use MC Neymars test to test significance of the base call */
		if(significantBase(bestScore, depthUpdate - bestScore, evalue) == 0) {
			if(bestNuc == '-' && tNuc != '-' && bestScore != depthUpdate) {
				bestBaseScore = 0;
				for(j = 0; j < 5; ++j) {
					if(bestBaseScore < calls->counts[j]) {
						bestBaseScore = calls->counts[j];
						bestNuc = j;
					}
				}
				if(bestBaseScore == 0) {
					bestNuc = '-';
				} else {
					bestNuc = tolower(bases[bestNuc]);
				}
			} else {
				bestNuc = tolower(bestNuc);
			}
		}
	}
	
	return bestNuc;
}

unsigned char refNanoCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls) {
	
	int j, bestBaseScore;
	const char bases[6] = "ACGTN-";
	
	/* determine base at current position */
	if(depthUpdate == 0) {
		bestNuc = 'n';
	} else {
		/* Use MC Neymars test to test significance of the base call */
		if(significantBase(bestScore, depthUpdate - bestScore, evalue) == 0) {
			if(bestNuc == '-') {
				bestBaseScore = 0;
				for(j = 0; j < 5; ++j) {
					if(bestBaseScore < calls->counts[j]) {
						bestBaseScore = calls->counts[j];
						bestNuc = j;
					}
				}
				if(bestBaseScore == 0) {
					bestNuc = 'n';
				} else {
					bestNuc = tolower(bases[bestNuc]);
				}
			} else {
				bestNuc = tolower(bestNuc);
			}
		} else if(bestNuc == '-') {
			bestNuc = 'n';
		}
	}
	
	return bestNuc;
}

void * assemble_KMA_threaded(void *arg) {
	
	static volatile int Lock[3] = {0, 0, 0}, mainTemplate = -2, thread_wait = 0;
	static char *template_name;
	static HashMapCCI *template_index;
	volatile int *excludeIn = &Lock[0], *excludeOut = &Lock[1], *excludeMatrix = &Lock[2];
	Assemble_thread *thread = arg;
	int i, j, t_len, aln_len, start, end, bias, myBias, gaps, pos, spin, sam;
	int read_score, depthUpdate, bestBaseScore, bestScore, template, asm_len;
	int nextTemplate, file_i, file_count, delta, thread_num, mq, status, bcd;
	int minlen, q_start, q_end, stats[5], buffer[8], *qBoundPtr;
	unsigned coverScore;
	long unsigned depth, depthVar;
	short unsigned *counts;
	const char bases[6] = "ACGTN-";
	double score, scoreT, mrc, evalue;
	unsigned char bestNuc;
	FILE **files, *file, *xml_out;
	AlnScore alnStat;
	Assembly *assembly;
	FileBuff *frag_out;
	Assem *aligned_assem;
	Aln *aligned, *gap_align;
	Qseqs *qseq, *header;
	AssemInfo *matrix;
	AlnPoints *points;
	NWmat *NWmatrices;
	
	/* get input */
	template = thread->template;
	file_count = thread->file_count;
	files = thread->files;
	xml_out = thread->xml_out;
	frag_out = thread->frag_out;
	aligned_assem = thread->aligned_assem;
	aligned = thread->aligned;
	gap_align = thread->gap_align;
	qseq = thread->qseq;
	header = thread->header;
	matrix = thread->matrix;
	points = thread->points;
	NWmatrices = thread->NWmatrices;
	delta = qseq->size;
	mq = thread->mq;
	minlen = thread->minlen;
	scoreT = thread->scoreT;
	mrc = thread->mrc;
	evalue = thread->evalue;
	bcd = thread->bcd;
	sam = thread->sam;
	spin = thread->spin;
	thread_num = thread->thread_num;
	
	if(template != -2) {
		/* all assemblies done, 
		signal threads to return */
		if(template == -1) {
			lock(excludeMatrix);
			mainTemplate = template;
			unlock(excludeMatrix);
			return NULL;
		}
		
		/* Allocate assembly arrays */
		lock(excludeMatrix);
		template_name = thread->template_name;
		template_index = thread->template_index;
		t_len = template_index->len;
		matrix->len = 0;
		
		/* start threads */
		aligned_assem->score = 0;
		aligned_assem->fragmentCountAln = 0;
		aligned_assem->readCountAln = 0;
		mainTemplate = template;
		thread_wait = thread_num;
		unlock(excludeMatrix);
		template = -2;
	}
	
	do {
		while(template == mainTemplate) {
			usleep(100);
		}
		lock(excludeMatrix);
		template = mainTemplate;
		if(template != -1) {
			t_len = template_index->len;
		}
		unlock(excludeMatrix);
		if(template == -1) {
			return NULL;
		}
		
		/* load reads of this template */
		file_i = 0;
		while(file_i < file_count) {
			lockTime(excludeIn, spin);
			file = files[file_i];
			if(file != 0) {
				read_score = fread(buffer, sizeof(int), 8, file);
				if((nextTemplate = buffer[0]) == template) {
					/* load frag */
					qseq->len = buffer[1];
					stats[0] = buffer[2];
					read_score = buffer[3];
					stats[2] = buffer[4];
					stats[3] = buffer[5];
					header->len = buffer[6];
					stats[4] = buffer[7];
					
					if(qseq->size < qseq->len) {
						free(qseq->seq);
						qseq->size = qseq->len << 1;
						qseq->seq = smalloc(qseq->size);
					}
					if(header->size < header->len) {
						header->size = header->len + 1;
						free(header->seq);
						header->seq = smalloc(header->size);
					}
					sfread(qseq->seq, 1, qseq->len, file);
					sfread(header->seq, 1, header->len, file);
					unlock(excludeIn);
					
					if(delta < qseq->len) {
						delta = qseq->len << 1;
						free(aligned->t);
						free(aligned->s);
						free(aligned->q);
						free(gap_align->t);
						free(gap_align->s);
						free(gap_align->q);
						aligned->t = smalloc((delta + 1) << 1);
						aligned->s = smalloc((delta + 1) << 1);
						aligned->q = smalloc((delta + 1) << 1);
						gap_align->t = smalloc((delta + 1) << 1);
						gap_align->s = smalloc((delta + 1) << 1);
						gap_align->q = smalloc((delta + 1) << 1);
					}
					
					/* q-bound */
					if(2 * sizeof(int) + 1 < header->len && header->seq[header->len - 2 * sizeof(int) - 1] == 0) {
						qBoundPtr = (int*) (header->seq + (header->len - 2 * sizeof(int)));
						q_start = *qBoundPtr;
						q_end = *++qBoundPtr;
					} else {
						q_start = 0;
						q_end = qseq->len;
					}
					
					/* Update assembly with read */
					if(read_score || anker_rc(template_index, qseq->seq, qseq->len, q_start, q_end, points)) {
						/* Start with alignment */
						if(stats[3] <= stats[2]) {
							stats[2] = 0;
							stats[3] = t_len;
						}
						
						alnStat = KMA(template_index, qseq->seq, qseq->len, q_start, q_end, aligned, gap_align, stats[2], MIN(t_len, stats[3]), mq, scoreT, points, NWmatrices);
						
						/* get read score */
						aln_len = alnStat.len;
						start = alnStat.pos;
						end = start + aln_len - alnStat.tGaps;
						
						/* Get normed score check read coverage */
						read_score = alnStat.score;
						if(minlen <= aln_len && mrcheck(mrc, alnStat, qseq->len, t_len)) {
							score = 1.0 * read_score / aln_len;
						} else {
							read_score = 0;
							score = 0;
						}
						
						if(0 < read_score && scoreT <= score) {
							stats[1] = read_score;
							stats[2] = start;
							stats[3] = end;
							if(t_len < end) {
								stats[3] -= t_len;
							}
							/* Update backbone and counts */
							//lock(excludeMatrix);
							lockTime(excludeMatrix, 10)
							aligned_assem->score += read_score;
							if(!(stats[4] & 2) || (stats[4] & 64)) {
								++aligned_assem->fragmentCountAln;
							}
							++aligned_assem->readCountAln;
							
							/* init matrix */
							if(!matrix->len) {
								matrix->len = t_len;
								if(matrix->size < (t_len << 1)) {
									matrix->size = (t_len << 1);
									free(matrix->assmb);
									matrix->assmb = smalloc(matrix->size * sizeof(Assembly));
								}
								
								/* cpy template seq */
								assembly = matrix->assmb - 1;
								i = 0;
								while(i != t_len) {
									counts = (++assembly)->counts;
									*counts = 0;
									*++counts = 0;
									*++counts = 0;
									*++counts = 0;
									*++counts = 0;
									*++counts = 0;
									assembly->next = ++i;
								}
								/* circularize */
								assembly->next = 0;
								/*
								for(i = 0, j = 1; i < t_len; ++i, ++j) {
									assembly[i].counts[0] = 0;
									assembly[i].counts[1] = 0;
									assembly[i].counts[2] = 0;
									assembly[i].counts[3] = 0;
									assembly[i].counts[4] = 0;
									assembly[i].counts[5] = 0;
									assembly[i].next = j;
								}
								assembly[t_len - 1].next = 0;
								*/
							}
							
							/* diff */
							i = 0;
							pos = start;
							assembly = matrix->assmb;
							while(i < aln_len) {
								if(aligned->t[i] == 5) { // Template gap, insertion
									if(t_len <= pos) {
										if(!++assembly[pos].counts[aligned->q[i]]) {
											assembly[pos].counts[aligned->q[i]] = USHRT_MAX;
										}
										++i;
										pos = assembly[pos].next;
									} else {
										/* get estimate for non insertions */
										myBias = 0;
										for(j = 0; j < 6; ++j) {
											myBias += assembly[pos].counts[j];
										}
										
										/* find position of insertion */
										gaps = pos;
										if(pos != 0) {
											--pos;
										} else {
											pos = t_len - 1;
										}
										while(assembly[pos].next != gaps) {
											pos = assembly[pos].next;
										}
										while(i < aln_len && aligned->t[i] == 5) {
											assembly[pos].next = matrix->len++;
											if(matrix->len == matrix->size) {
												matrix->size <<= 1;
												matrix->assmb = realloc(assembly, matrix->size * sizeof(Assembly));
												if(!matrix->assmb) {
													matrix->size >>= 1;
													matrix->size += 1024;
													matrix->assmb = realloc(assembly, matrix->size * sizeof(Assembly));
													if(!matrix->assmb) {
														ERROR();
													}
												}
												assembly = matrix->assmb;
											}
											pos = assembly[pos].next;
											assembly[pos].next = gaps;
											assembly[pos].counts[0] = 0;
											assembly[pos].counts[1] = 0;
											assembly[pos].counts[2] = 0;
											assembly[pos].counts[3] = 0;
											assembly[pos].counts[4] = 0;
											assembly[pos].counts[5] = myBias < USHRT_MAX ? myBias : USHRT_MAX;
											assembly[pos].counts[aligned->q[i]] = 1;
											
											++i;
										}
										pos = assembly[pos].next;
									}
								} else if(t_len <= pos) { // Old template gap, not present in this read
									if(!++assembly[pos].counts[5]) {
										assembly[pos].counts[5] = USHRT_MAX;
									}
									pos = assembly[pos].next;
								} else {
									if(!++assembly[pos].counts[aligned->q[i]]) {
										assembly[pos].counts[aligned->q[i]] = USHRT_MAX;
									}
									++i;
									pos = assembly[pos].next;
								}
							}
							unlock(excludeMatrix);
							
							/* Convert fragment */
							for(i = 0; i < qseq->len; ++i) {
								 qseq->seq[i] = bases[qseq->seq[i]];
							}
							qseq->seq[qseq->len] = 0;
							
							/* Save fragment */
							if(frag_out) {
								lockTime(excludeOut, 10);
								updateFrags(frag_out, qseq, header, template_name, stats);
								unlock(excludeOut);
							}
							
							if(sam) {
								header->seq[header->len - 1] = 0;
								samwrite(qseq, header, 0, template_name, aligned, stats);
							}
							
							if(xml_out) {
								hitXML(xml_out, template, header->seq, aligned, &alnStat, NWmatrices->rewards, stats[4]);
							}
						} else if(sam && !(sam & 2096)) {
							stats[1] = read_score;
							stats[2] = start;
							stats[3] = end;
							header->seq[header->len - 1] = 0;
							nibble2base(qseq->seq, qseq->len);
							if(read_score) {
								samwrite(qseq, header, 0, template_name, aligned, stats);
							} else {
								stats[4] |= 4;
								stats[1] = stats[4];
								samwrite(qseq, header, 0, template_name, 0, stats);
							}
						}
					} else if(sam && !(sam & 2096)) {
						stats[4] |= 4;
						stats[1] = stats[4];
						header->seq[header->len - 1] = 0;
						nibble2base(qseq->seq, qseq->len);
						samwrite(qseq, header, 0, template_name, 0, stats);
					}
				} else if(nextTemplate == -1) {
					if(template) {
						fclose(file);
					} else {
						kmaPipe(0, 0, file, &status);
						errno |= status;
					}
					files[file_i] = 0;
					unlock(excludeIn);
					++file_i;
				} else if(nextTemplate < template) {
					/* Move pointer forward */
					sfseek(file, buffer[1] + buffer[6], SEEK_CUR);
					unlock(excludeIn);
				} else {
					/* Move pointer back */
					sfseek(file, (-8) * sizeof(int), SEEK_CUR);
					unlock(excludeIn);
					++file_i;
				}
			} else {
				unlock(excludeIn);
				++file_i;
			}
		}
		lock(excludeIn);
		--thread_wait;
		unlock(excludeIn);
	} while(thread->num != 0);
	
	wait_atomic(thread_wait);
	
	if(aligned_assem->score == 0) {
		aligned_assem->cover = 0;
		aligned_assem->depth = 0;
		aligned_assem->depthVar = 0;
		aligned_assem->t[0] = 0;
		aligned_assem->s[0] = 0;
		aligned_assem->q[0] = 0;
		aligned_assem->len = 0;
		aligned_assem->aln_len = 0;
		
		return NULL;
	}
	
	/* diff */
	/* pre on dense */
	/* Pepare and make alignment on consensus */
	asm_len = matrix->len;
	assembly = matrix->assmb;
	if(aligned_assem->size <= asm_len) {
		aligned_assem->size = (asm_len + 1) << 1;
		free(aligned_assem->t);
		free(aligned_assem->s);
		free(aligned_assem->q);
		aligned_assem->t = smalloc(aligned_assem->size);
		aligned_assem->s = smalloc(aligned_assem->size);
		aligned_assem->q = smalloc(aligned_assem->size);
	}
	
	/* Call nucleotides for the consensus */
	/* diff */
	i = 0;
	pos = 0;
	depth = 0;
	depthVar = 0;
	aln_len = 0;
	while(i < asm_len) {
		/* call template */
		if(pos < t_len) {
			bestNuc = getNuc(template_index->seq, pos);
		} else {
			bestNuc = 5;
		}
		aligned_assem->t[i] = bases[bestNuc];
		
		/* call query */
		bestScore = assembly[pos].counts[bestNuc];
		depthUpdate = 0;
		for(j = 0; j < 6; ++j) {
			if(bestScore < assembly[pos].counts[j]) {
				bestScore = assembly[pos].counts[j];
				bestNuc = j;
			}
			depthUpdate += assembly[pos].counts[j];
		}
		bestNuc = bases[bestNuc];
		
		/* check for minor base call */
		if(!depthUpdate) {
			bestNuc = '-';
		} else if((bestScore << 1) < depthUpdate) {
			if(bestNuc == '-') {
				bestBaseScore = assembly[pos].counts[4];
				bestNuc = 4;
				for(j = 0; j < 4; ++j) {
					if(bestBaseScore < assembly[pos].counts[j]) {
						bestBaseScore = assembly[pos].counts[j];
						bestNuc = j;
					}
				}
				bestNuc = tolower(bases[bestNuc]);
			} else {
				bestNuc = tolower(bestNuc);
			}
			bestScore = depthUpdate - assembly[pos].counts[5];
		} else if(depthUpdate < bcd) {
			/* too low depth */
			bestNuc = tolower(bestNuc);
		}
		
		/* determine base at current position */
		/*
		if(bcd <= depthUpdate) {
			bestNuc = baseCall(bestNuc, aligned_assem->t[i], bestScore, depthUpdate, evalue, &assembly[pos]);
		} else {
			bestNuc = baseCall('-', aligned_assem->t[i], 0, 0, evalue, &assembly[pos]);
		}
		*/
		bestNuc = baseCall(bestNuc, aligned_assem->t[i], bestScore, depthUpdate, evalue, &assembly[pos]);
		aligned_assem->q[i] = bestNuc;
		
		if(bestNuc != '-') {
			depth += depthUpdate;
			depthVar += (depthUpdate * depthUpdate);
			++aln_len;
		}
		
		++i;
		pos = assembly[pos].next;
	}
	
	/* Trim alignment on consensus */
	coverScore = 0;
	bias = 0;
	for(i = 0; i < asm_len; ++i) {
		if(aligned_assem->t[i] == '-' && aligned_assem->q[i] == '-') {
			++bias;
		} else {
			aligned_assem->t[i - bias] = aligned_assem->t[i];
			aligned_assem->q[i - bias] = aligned_assem->q[i];
			if(tolower(aligned_assem->t[i]) == tolower(aligned_assem->q[i])) {
				aligned_assem->s[i - bias] = '|';
				++coverScore;
			} else {
				aligned_assem->s[i - bias] = '_';
			}
		}
	}
	asm_len -= bias;
	aligned_assem->t[asm_len] = 0;
	aligned_assem->s[asm_len] = 0;
	aligned_assem->q[asm_len] = 0;
	aligned_assem->cover = coverScore;
	aligned_assem->depth = depth;
	aligned_assem->depthVar = depthVar;
	aligned_assem->len = asm_len;
	aligned_assem->aln_len = aln_len;
	
	return NULL;
}

void * assemble_KMA_dense_threaded(void *arg) {
	
	static volatile int Lock[3] = {0, 0, 0}, mainTemplate = -2, thread_wait = 0;
	static char *template_name;
	static HashMapCCI *template_index;
	volatile int *excludeIn = &Lock[0], *excludeOut = &Lock[1], *excludeMatrix = &Lock[2];
	Assemble_thread *thread = arg;
	int i, j, t_len, aln_len, start, end, file_i, file_count, template, spin;
	int pos, read_score, bestScore, depthUpdate, bestBaseScore, nextTemplate;
	int sam, thread_num, mq, status, bcd, minlen, q_start, q_end, *qBoundPtr;
	int stats[5], buffer[8];
	unsigned coverScore, delta;
	long unsigned depth, depthVar;
	short unsigned *counts;
	const char bases[6] = "ACGTN-";
	double score, scoreT, mrc, evalue;
	unsigned char bestNuc;
	FILE **files, *file, *xml_out;
	AlnScore alnStat;
	Assembly *assembly;
	FileBuff *frag_out;
	Assem *aligned_assem;
	Aln *aligned, *gap_align;
	Qseqs *qseq, *header;
	AssemInfo *matrix;
	AlnPoints *points;
	NWmat *NWmatrices;
	
	/* get input */
	template = thread->template;
	file_count = thread->file_count;
	files = thread->files;
	xml_out = thread->xml_out;
	frag_out = thread->frag_out;
	aligned_assem = thread->aligned_assem;
	aligned = thread->aligned;
	gap_align = thread->gap_align;
	qseq = thread->qseq;
	header = thread->header;
	matrix = thread->matrix;
	points = thread->points;
	NWmatrices = thread->NWmatrices;
	delta = qseq->size;
	mq = thread->mq;
	minlen = thread->minlen;
	scoreT = thread->scoreT;
	mrc = thread->mrc;
	evalue = thread->evalue;
	bcd = thread->bcd;
	sam = thread->sam;
	spin = thread->spin;
	thread_num = thread->thread_num;
	
	if(template != -2) {
		/* all assemblies done, 
		signal threads to return */
		if(template == -1) {
			lock(excludeOut);
			mainTemplate = template;
			unlock(excludeOut);
			return NULL;
		}
		
		/* Allocate assembly arrays */
		lock(excludeOut);
		template_name = thread->template_name;
		template_index = thread->template_index;
		t_len = template_index->len;
		matrix->len = t_len;
		
		/* diff */
		if(aligned_assem->size <= t_len) {
			aligned_assem->size = t_len + 1;
			free(aligned_assem->t);
			free(aligned_assem->s);
			free(aligned_assem->q);
			aligned_assem->t = smalloc(t_len + 1);
			aligned_assem->s = smalloc(t_len + 1);
			aligned_assem->q = smalloc(t_len + 1);
		}
		if(matrix->size <= t_len) {
			matrix->size = t_len + 1;
			free(matrix->assmb);
			matrix->assmb = smalloc(matrix->size * sizeof(Assembly));
		}
		
		/* cpy template seq */
		assembly = matrix->assmb - 1;
		i = 0;
		while(i != t_len) {
			aligned_assem->t[i] = getNuc(template_index->seq, i);
			counts = (++assembly)->counts;
			*counts = 0;
			*++counts = 0;
			*++counts = 0;
			*++counts = 0;
			*++counts = 0;
			*++counts = 0;
			assembly->next = ++i;
		}
		/* circularize */
		assembly->next = 0;
		assembly = matrix->assmb;
		/*
		for(i = 0, j = 1; i < t_len; ++i, ++j) {
			assembly[i].counts[0] = 0;
			assembly[i].counts[1] = 0;
			assembly[i].counts[2] = 0;
			assembly[i].counts[3] = 0;
			assembly[i].counts[4] = 0;
			assembly[i].counts[5] = 0;
			assembly[i].next = j;
		}
		assembly[t_len - 1].next = 0;
		*/
		
		/* start threads */
		aligned_assem->score = 0;
		aligned_assem->fragmentCountAln = 0;
		aligned_assem->readCountAln = 0;
		mainTemplate = template;
		thread_wait = thread_num;
		unlock(excludeOut);
		template = -2;
	}
	
	do {
		while(template == mainTemplate) {
			usleep(100);
		}
		lock(excludeOut);
		template = mainTemplate;
		if(template != -1) {
			t_len = template_index->len;
			assembly = matrix->assmb;
		}
		unlock(excludeOut);
		if(template == -1) {
			return NULL;
		}
		
		/* load reads of this template */
		file_i = 0;
		while(file_i < file_count) {
			lockTime(excludeIn, spin);
			file = files[file_i];
			if(file != 0) {
				read_score = fread(buffer, sizeof(int), 8, file);
				if((nextTemplate = buffer[0]) == template) {
					/* load frag */
					qseq->len = buffer[1];
					stats[0] = buffer[2];
					read_score = buffer[3];
					stats[2] = buffer[4];
					stats[3] = buffer[5];
					header->len = buffer[6];
					stats[4] = buffer[7];
					
					if(qseq->size < qseq->len) {
						free(qseq->seq);
						qseq->size = qseq->len << 1;
						qseq->seq = malloc(qseq->size);
						if(!qseq->seq) {
							ERROR();
						}
					}
					if(header->size < header->len) {
						header->size = header->len + 1;
						free(header->seq);
						header->seq = malloc(header->size);
						if(!header->seq) {
							ERROR();
						}
					}
					sfread(qseq->seq, 1, qseq->len, file);
					sfread(header->seq, 1, header->len, file);
					unlock(excludeIn);
					
					if(delta < qseq->size) {
						delta = qseq->size;
						free(aligned->t);
						free(aligned->s);
						free(aligned->q);
						free(gap_align->t);
						free(gap_align->s);
						free(gap_align->q);
						aligned->t = malloc((delta + 1) << 1);
						aligned->s = malloc((delta + 1) << 1);
						aligned->q = malloc((delta + 1) << 1);
						gap_align->t = malloc((delta + 1) << 1);
						gap_align->s = malloc((delta + 1) << 1);
						gap_align->q = malloc((delta + 1) << 1);
						if(!aligned->t || !aligned->s || !aligned->q || !gap_align->t || !gap_align->s || !gap_align->q) {
							ERROR();
						}
					}
					
					/* q-bound */
					if(2 * sizeof(int) + 1 < header->len && header->seq[header->len - 2 * sizeof(int) - 1] == 0) {
						qBoundPtr = (int*) (header->seq + (header->len - 2 * sizeof(int)));
						q_start = *qBoundPtr;
						q_end = *++qBoundPtr;
					} else {
						q_start = 0;
						q_end = qseq->len;
					}
					
					/* Update assembly with read */
					if(read_score || anker_rc(template_index, qseq->seq, qseq->len, q_start, q_end, points)) {
						if(stats[3] <= stats[2]) {
							stats[2] = 0;
							stats[3] = t_len;
						}
						/* Start with alignment */
						alnStat = KMA(template_index, qseq->seq, qseq->len, q_start, q_end, aligned, gap_align, stats[2], MIN(t_len, stats[3]), mq, scoreT, points, NWmatrices);
						
						/* get read score */
						aln_len = alnStat.len;
						start = alnStat.pos;
						end = start + aln_len - alnStat.tGaps;
						
						/* Get normed score check read coverage */
						read_score = alnStat.score;
						if(minlen <= aln_len && mrcheck(mrc, alnStat, qseq->len, t_len)) {
							score = 1.0 * read_score / aln_len;
						} else {
							read_score = 0;
							score = 0;
						}
						
						if(0 < read_score && scoreT <= score) {
							
							stats[1] = read_score;
							stats[2] = start;
							stats[3] = end;
							if(t_len < end) {
								stats[3] -= t_len;
							}
							/* Update backbone and counts */
							//lock(excludeMatrix);
							lockTime(excludeMatrix, 10)
							aligned_assem->score += read_score;
							if(!(stats[4] & 2) || (stats[4] & 64)) {
								++aligned_assem->fragmentCountAln;
							}
							++aligned_assem->readCountAln;
							
							/* diff */
							for(i = 0, pos = start; i < aln_len; ++i) {
								if(aligned->t[i] == aligned_assem->t[pos]) {
									if(!++assembly[pos].counts[aligned->q[i]]) {
										assembly[pos].counts[aligned->q[i]] = USHRT_MAX;	
									}
									pos = assembly[pos].next;
								}
							}
							unlock(excludeMatrix);
							
							/* Convert fragment */
							for(i = 0; i < qseq->len; ++i) {
								 qseq->seq[i] = bases[qseq->seq[i]];
							}
							qseq->seq[qseq->len] = 0;
							
							/* Save fragment */
							if(frag_out) {
								lockTime(excludeOut, 10);
								updateFrags(frag_out, qseq, header, template_name, stats);
								unlock(excludeOut);
							}
							
							if(sam) {
								header->seq[header->len - 1] = 0;
								samwrite(qseq, header, 0, template_name, aligned, stats);
							}
							
							if(xml_out) {
								hitXML(xml_out, template, header->seq, aligned, &alnStat, NWmatrices->rewards, stats[4]);
							}
						} else if(sam && !(sam & 2096)) {
							stats[1] = read_score;
							stats[2] = start;
							stats[3] = end;
							/* flag */
							header->seq[header->len - 1] = 0;
							nibble2base(qseq->seq, qseq->len);
							if(read_score) {
								samwrite(qseq, header, 0, template_name, aligned, stats);
							} else {
								stats[4] |= 4;
								stats[1] = stats[4];
								samwrite(qseq, header, 0, template_name, 0, stats);
							}
						}
					} else if(sam && !(sam & 2096)) {
						stats[4] |= 4;
						stats[1] = stats[4];
						header->seq[header->len - 1] = 0;
						nibble2base(qseq->seq, qseq->len);
						samwrite(qseq, header, 0, template_name, 0, stats);
					}
				} else if (nextTemplate == -1) {
					if(template) {
						fclose(file);
					} else {
						kmaPipe(0, 0, file, &status);
						errno |= status;
					}
					files[file_i] = 0;
					unlock(excludeIn);
					++file_i;
				} else if(nextTemplate < template) {
					/* Move pointer forward */
					sfseek(file, buffer[1] + buffer[6], SEEK_CUR);
					unlock(excludeIn);
				} else {
					/* Move pointer back */
					sfseek(file, (-8) * sizeof(int), SEEK_CUR);
					unlock(excludeIn);
					++file_i;
				}
			} else {
				unlock(excludeIn);
				++file_i;
			}
		}
		
		lock(excludeIn);
		--thread_wait;
		unlock(excludeIn);
		
	} while(thread->num != 0);
	
	wait_atomic(thread_wait);
	
	if(aligned_assem->score == 0) {
		aligned_assem->cover = 0;
		aligned_assem->depth = 0;
		aligned_assem->depthVar = 0;
		aligned_assem->t[0] = 0;
		aligned_assem->s[0] = 0;
		aligned_assem->q[0] = 0;
		aligned_assem->len = 0;
		aligned_assem->aln_len = 0;
		
		return NULL;
	}
	
	/* Make consensus assembly by majority voting */
	/* diff */
	depth = 0;
	depthVar = 0;
	coverScore = 0;
	aln_len = 0;
	for(i = 0; i < t_len; ++i) {
		/* call template */
		bestNuc = getNuc(template_index->seq, i);
		aligned_assem->t[i] = bases[bestNuc];
		
		/* call query */
		bestScore = assembly[i].counts[bestNuc];
		depthUpdate = 0;
		for(j = 0; j < 6; ++j) {
			if(bestScore < assembly[i].counts[j]) {
				bestScore = assembly[i].counts[j];
				bestNuc = j;
			}
			depthUpdate += assembly[i].counts[j];
		}
		bestNuc = bases[bestNuc];
		
		/* Check for minor base call */
		if(!depthUpdate) {
			bestNuc = '-';
		} else if((bestScore << 1) < depthUpdate) {
			if(bestNuc == '-') {
				bestBaseScore = 0;
				bestNuc = 4;
				for(j = 0; j < 5; ++j) {
					if(bestBaseScore < assembly[i].counts[j]) {
						bestBaseScore = assembly[i].counts[j];
						bestNuc = j;
					}
				}
				bestNuc = tolower(bases[bestNuc]);
			} else {
				bestNuc = tolower(bestNuc);
			}
			bestScore = depthUpdate - assembly[i].counts[5];
		} else if(depthUpdate < bcd) {
			/* too low depth */
			bestNuc = tolower(bestNuc);
		}
		
		/* determine base at current position */
		/*
		if(bcd <= depthUpdate) {
			bestNuc = baseCall(bestNuc, aligned_assem->t[i], bestScore, depthUpdate, evalue, &assembly[i]);
		} else {
			bestNuc = baseCall('-', aligned_assem->t[i], 0, 0, evalue, &assembly[i]);
		}
		*/
		bestNuc = baseCall(bestNuc, aligned_assem->t[i], bestScore, depthUpdate, evalue, &assembly[i]);
		aligned_assem->q[i] = bestNuc;
		
		if(bestNuc != '-') {
			depth += depthUpdate;
			depthVar += (depthUpdate * depthUpdate);
			++aln_len;
		}
		
		if(tolower(aligned_assem->q[i]) == tolower(aligned_assem->t[i])) {
			aligned_assem->s[i] = '|';
			++coverScore;
		} else {
			aligned_assem->s[i] = '_';
		}
	}
	aligned_assem->t[t_len] = 0;
	aligned_assem->s[t_len] = 0;
	aligned_assem->q[t_len] = 0;
	aligned_assem->cover = coverScore;
	aligned_assem->depth = depth;
	aligned_assem->depthVar = depthVar;
	aligned_assem->len = t_len;
	aligned_assem->aln_len = aln_len;
	
	return NULL;
}

void * skip_assemble_KMA(void *arg) {
	
	Assemble_thread *thread = arg;
	int template, t_len, sam, nextTemplate, file_count, file_i, status;
	int stats[5], buffer[8];
	char *template_name;
	FILE *file, **files;
	Assem *aligned_assem;
	FileBuff *frag_out;
	Qseqs *qseq, *header;
	
	if((template = thread->template) < 0) {
		return NULL;
	}
	
	/* get input */
	sam = thread->sam;
	t_len = thread->template_index->len;
	template_name = thread->template_name;
	file_count = thread->file_count;
	files = thread->files;
	frag_out = thread->frag_out;
	aligned_assem = thread->aligned_assem;
	qseq = thread->qseq;
	header = thread->header;
	thread->matrix->len = 0;
	
	aligned_assem->var = 0;
	aligned_assem->nucHighVar = 0;
	aligned_assem->maxDepth = 0;
	aligned_assem->snpSum = 0;
	aligned_assem->insertSum = 0;
	aligned_assem->deletionSum = 0;
	aligned_assem->score = 0;
	aligned_assem->cover = 0;
	aligned_assem->depth = 0;
	aligned_assem->depthVar = 0;
	aligned_assem->t[0] = 0;
	aligned_assem->s[0] = 0;
	aligned_assem->q[0] = 0;
	aligned_assem->len = t_len;
	aligned_assem->aln_len = t_len;
	aligned_assem->fragmentCountAln = 0;
	aligned_assem->readCountAln = 0;
	
	/* load reads of this template */
	file_i = 0;
	while(file_i < file_count) {
		file = files[file_i];
		if(file != 0) {
			nextTemplate = fread(buffer, sizeof(int), 8, file);
			if((nextTemplate = buffer[0]) == template) {
				/* load frag */
				qseq->len = buffer[1];
				stats[0] = buffer[2];
				stats[2] = buffer[4];
				stats[3] = buffer[5];
				header->len = buffer[6];
				stats[4] = buffer[7];
				
				if(qseq->size < qseq->len) {
					free(qseq->seq);
					qseq->size = qseq->len << 1;
					qseq->seq = smalloc(qseq->size);
				}
				if(header->size < header->len) {
					header->size = header->len + 1;
					free(header->seq);
					header->seq = smalloc(header->size);
				}
				sfread(qseq->seq, 1, qseq->len, file);
				sfread(header->seq, 1, header->len, file);
				
				/* Update with read */
				aligned_assem->depth += qseq->len;
				if(sam) {
					stats[4] |= 4;
					stats[1] = stats[4];
					header->seq[header->len - 1] = 0;
					nibble2base(qseq->seq, qseq->len);
					samwrite(qseq, header, 0, template_name, 0, stats);
				}
				/* Save fragment */
				if(frag_out) {
					updateFrags(frag_out, qseq, header, template_name, stats);
				}
			} else if (nextTemplate == -1) {
				if(template) {
					fclose(file);
				} else {
					kmaPipe(0, 0, file, &status);
					errno |= status;
				}
				files[file_i] = 0;
				++file_i;
			} else if(nextTemplate < template) {
				/* Move pointer forward */
				sfseek(file, buffer[1] + buffer[6], SEEK_CUR);
			} else {
				/* Move pointer back */
				sfseek(file, (-8) * sizeof(int), SEEK_CUR);
				++file_i;
			}
		} else {
			++file_i;
		}
	}
	aligned_assem->cover = 0; 
	aligned_assem->aln_len = 0;//(1 - exp((-1.0) * aligned_assem->depth / t_len)) * t_len; // expected coverage from depth
	
	return NULL;
}

void alnToMat(AssemInfo *matrix, Assem *aligned_assem, Aln *aligned, AlnScore alnStat, int t_len, int flag) {
	
	static volatile int Lock = 0;
	volatile int *excludeMatrix = &Lock;
	int i, pos, aln_len, start, read_score, myBias, tmp, gaps;
	short unsigned *counts;
	Assembly *assembly;
	
	/* init */
	aln_len = alnStat.len;
	start = alnStat.pos;
	read_score = alnStat.score;
	
	
	/* Update backbone and counts */
	//lockTime(excludeMatrix, 10)
	lock(excludeMatrix);
	aligned_assem->score += read_score;
	if(!(flag & 2) || (flag & 64)) {
		++aligned_assem->fragmentCountAln;
	}
	++aligned_assem->readCountAln;
	
	/* trim trailing gaps */
	i = aln_len - 1;
	while(i && (aligned->t[i] == 5 || aligned->q[i] == 5)) {
		--i;
	}
	aln_len = i + 1;
	
	/* trim leading gaps */
	i = 0;
	while(i < aln_len && (aligned->t[i] == 5 || aligned->q[i] == 5)) {
		if(aligned->q[i] == 5) {
			++start;
		}
		++i;
	}
	
	/* diff */
	pos = start;
	assembly = matrix->assmb;
	while(i < aln_len) {
		if(aligned->t[i] == 5) { // Template gap, insertion
			if(t_len <= pos) {
				if(!++assembly[pos].counts[aligned->q[i]]) {
					assembly[pos].counts[aligned->q[i]] = USHRT_MAX;
				}
				++i;
				pos = assembly[pos].next;
			} else {
				/* check if insertion already present at pos */
				gaps = pos;
				pos = pos ? (pos - 1) : (t_len - 1);
				
				/* find position of insertion, if it already exists */
				while(assembly[pos].next != gaps) {
					pos = assembly[pos].next;
				}
				
				/* get number of bases not supporting the insertion */
				counts = assembly[pos].counts;
				myBias = *counts;
				myBias += *++counts;
				myBias += *++counts;
				myBias += *++counts;
				myBias += *++counts;
				myBias += *++counts;
				
				/* get minimum of bases sorrounding new insertion */
				counts = assembly[gaps].counts;
				tmp = *counts;
				tmp += *++counts;
				tmp += *++counts;
				tmp += *++counts;
				tmp += *++counts;
				tmp += *++counts;
				myBias = (tmp < myBias) ? tmp : (myBias - 1);
				if(USHRT_MAX < myBias) {
					myBias = USHRT_MAX;
				}
				
				/* update new insertion */
				while(i < aln_len && aligned->t[i] == 5) {
					assembly[pos].next = matrix->len++;
					if(matrix->len == matrix->size) {
						matrix->size <<= 1;
						matrix->assmb = realloc(assembly, matrix->size * sizeof(Assembly));
						if(!matrix->assmb) {
							matrix->size >>= 1;
							matrix->size += 1024;
							matrix->assmb = realloc(assembly, matrix->size * sizeof(Assembly));
							if(!matrix->assmb) {
								ERROR();
							}
						}
						assembly = matrix->assmb;
					}
					pos = assembly[pos].next;
					assembly[pos].next = gaps;
					counts = assembly[pos].counts;
					counts[0] = 0;
					counts[1] = 0;
					counts[2] = 0;
					counts[3] = 0;
					counts[4] = 0;
					counts[5] = myBias;
					counts[aligned->q[i]] = 1;
					
					++i;
				}
				pos = assembly[pos].next;
			}
		} else if(t_len <= pos) { // Old template gap, not present in this read
			if(!++assembly[pos].counts[5]) {
				assembly[pos].counts[5] = USHRT_MAX;
			}
			pos = assembly[pos].next;
		} else {
			if(!++assembly[pos].counts[aligned->q[i]]) {
				assembly[pos].counts[aligned->q[i]] = USHRT_MAX;
			}
			++i;
			pos = assembly[pos].next;
		}
	}
	unlock(excludeMatrix);
}

void alnToMatDense(AssemInfo *matrix, Assem *aligned_assem, Aln *aligned, AlnScore alnStat, int t_len, int flag) {
	
	static volatile int Lock = 0;
	volatile int *excludeMatrix = &Lock;
	int i, pos, aln_len, start, read_score;
	Assembly *assembly;
	
	/* init */
	aln_len = alnStat.len;
	start = alnStat.pos;
	read_score = alnStat.score;
	
	
	/* Update backbone and counts */
	//lockTime(excludeMatrix, 10)
	lock(excludeMatrix);
	aligned_assem->score += read_score;
	if(!(flag & 2) || (flag & 64)) {
		++aligned_assem->fragmentCountAln;
	}
	++aligned_assem->readCountAln;
	
	/* trim trailing gaps */
	i = aln_len - 1;
	while(i < aln_len && (aligned->t[i] == 5 || aligned->q[i] == 5)) {
		--i;
	}
	aln_len = i + 1;
	
	/* trim leading gaps */
	i = 0;
	while(i && (aligned->t[i] == 5 || aligned->q[i] == 5)) {
		if(aligned->q[i] == 5) {
			++start;
		}
		++i;
	}
	
	/* diff */
	pos = start;
	assembly = matrix->assmb;
	while(i < aln_len) {
		if(aligned->t[i] != 5) {
			if(!++assembly[pos].counts[aligned->q[i]]) {
				assembly[pos].counts[aligned->q[i]] = USHRT_MAX;	
			}
			pos = assembly[pos].next;
		}
		++i;
	}
	unlock(excludeMatrix);
}

void callConsensus(AssemInfo *matrix, Assem *aligned_assem, long unsigned *seq, int t_len, int bcd, double evalue, int thread_num) {
	
	const char bases[6] = "ACGTN-";
	static volatile int Lock = 0, next, thread_wait = 0;
	volatile int *excludeMatrix = &Lock;
	int i, j, pos, end ,asm_len, aln_len, bestScore, bestBaseScore, chunk;
	int coverScore;
	long unsigned depth, depthVar, depthUpdate;
	unsigned char bestNuc;
	Assembly *assembly;
	
	/* init */
	lock(excludeMatrix);
	if(!thread_wait) {
		next = 0;
		thread_wait = thread_num;
	}
	unlock(excludeMatrix);
	asm_len = matrix->len;
	assembly = matrix->assmb;
	depth = 0;
	depthVar = 0;
	aln_len = 0;
	coverScore = 0;
	
	/* call in chunk of 8112 ~= 1 MB */
	chunk = alnToMatPtr != &alnToMatDense ? asm_len : 8112;
	while(chunk) {
		lock(excludeMatrix);
		/* get next chunk */
		i = next;
		if((next += chunk) < 0) {
			next = asm_len;
		}
		unlock(excludeMatrix);
		
		/* call chunk */
		if(i < asm_len) {
			end = i + chunk;
			if(asm_len < end) {
				end = asm_len;
				chunk = 0;
			}
			pos = i;
			while(i < end) {
				/* call template */
				if(pos < t_len) {
					bestNuc = getNuc(seq, pos);
				} else {
					bestNuc = 5;
				}
				aligned_assem->t[i] = bases[bestNuc];
				
				/* call query */
				bestScore = assembly[pos].counts[bestNuc];
				depthUpdate = 0;
				for(j = 0; j < 6; ++j) {
					if(bestScore < assembly[pos].counts[j]) {
						bestScore = assembly[pos].counts[j];
						bestNuc = j;
					}
					depthUpdate += assembly[pos].counts[j];
				}
				bestNuc = bases[bestNuc];
				
				/* check for minor base call */
				if(!depthUpdate) {
					bestNuc = '-';
				} else if((bestScore << 1) < depthUpdate) {
					if(bestNuc == '-') {
						bestBaseScore = assembly[pos].counts[4];
						bestNuc = 4;
						for(j = 0; j < 4; ++j) {
							if(bestBaseScore < assembly[pos].counts[j]) {
								bestBaseScore = assembly[pos].counts[j];
								bestNuc = j;
							}
						}
						bestNuc = tolower(bases[bestNuc]);
					} else {
						bestNuc = tolower(bestNuc);
					}
					bestScore = depthUpdate - assembly[pos].counts[5];
				} else if(depthUpdate < bcd) {
					/* too low depth */
					bestNuc = tolower(bestNuc);
				}
				
				/* determine base at current position */
				/*
				if(bcd <= depthUpdate) {
					bestNuc = baseCall(bestNuc, aligned_assem->t[i], bestScore, depthUpdate, evalue, &assembly[pos]);
				} else {
					bestNuc = baseCall('-', aligned_assem->t[i], 0, 0, evalue, &assembly[pos]);
				}
				*/
				bestNuc = baseCall(bestNuc, aligned_assem->t[i], bestScore, depthUpdate, evalue, &assembly[pos]);
				aligned_assem->q[i] = bestNuc;
				
				if(bestNuc != '-') {
					depth += depthUpdate;
					depthVar += (depthUpdate * depthUpdate);
					++aln_len;
					
					if(pos < t_len && aligned_assem->t[i] == toupper(bestNuc)) {
						++coverScore;
						aligned_assem->s[i] = '|';
					} else {
						aligned_assem->s[i] = '_';
					}
				} else {
					aligned_assem->s[i] = '_';
				}
				
				++i;
				pos = assembly[pos].next;
			}
		} else {
			chunk = 0;
		}
	}
	
	/* update aligned_assem */
	lock(excludeMatrix);
	aligned_assem->depth += depth;
	aligned_assem->depthVar += depthVar;
	aligned_assem->len = asm_len;
	aligned_assem->aln_len += aln_len;
	aligned_assem->cover += coverScore;
	--thread_wait;
	unlock(excludeMatrix);
	wait_atomic(thread_wait);
}

void fixVarOverflow(Assem *aligned_assem, Assembly *assembly, int t_len, int thread_num) {
	
	static volatile int Lock = 0, next, thread_wait = 0;
	volatile int *excludeMatrix = &Lock;
	int pos, end, chunk, depthUpdate;
	double var, depth, tmp;
	
	/*
	var(x) = E(x^2) - E(x)^2 = sum((x_i - E(x))^2) / n
	*/
	
	tmp = aligned_assem->depth;
	tmp /= t_len;
	var = aligned_assem->depthVar;
	var /= t_len;
	var -= (tmp * tmp);
	if(0 <= var) {
		aligned_assem->var = var;
		return;
	}
	
	/* init */
	var = 0;
	depth = (long double)(aligned_assem->depth) / t_len;
	chunk = 8112;
	lock(excludeMatrix);
	if(!thread_wait) {
		next = 0;
		thread_wait = thread_num;
	}
	unlock(excludeMatrix);
	
	/* get variance */
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
			if(t_len < end) {
				end = t_len;
				chunk = 0;
			}
			while(pos < end) {
				depthUpdate = assembly[pos].counts[0] + assembly[pos].counts[1] + assembly[pos].counts[2] + assembly[pos].counts[3] + assembly[pos].counts[4] + assembly[pos].counts[5];
				tmp = (depthUpdate - depth);
				var += tmp * tmp / t_len;
				++pos;
			}
		} else {
			chunk = 0;
		}
	}
	
	lock(excludeMatrix);
	aligned_assem->var += var;
	--thread_wait;
	unlock(excludeMatrix);
	wait_atomic(thread_wait);
}

void * assemble_KMA(void *arg) {
	
	const char bases[6] = "ACGTN-";
	static volatile int thread_wait = 0, thread_init = 0, thread_begin = 0;
	static volatile int mainTemplate = -2, next, Lock[3] = {0, 0, 0};
	static int t_len, load, seq_in;
	static char *template_name;
	static HashMapCCI *template_index;
	volatile int *excludeIn = &Lock[0], *excludeOut = &Lock[1], *excludeMatrix = &Lock[2];
	Assemble_thread *thread = arg;
	int i, minlen, aln_len, kmersize, sam, chunk, ef, template;
	int read_score, asm_len, nextTemplate, file_i, file_count, delta, status;
	int thread_num, mq, bcd, start, end, q_start, q_end, Wl;
	int stats[5], buffer[8], *qBoundPtr;
	short unsigned *counts;
	double score, scoreT, mrc, evalue;
	unsigned char *q;
	FILE **files, *file, *xml_out;
	AlnScore alnStat;
	Assembly *assembly;
	FileBuff *frag_out;
	Assem *aligned_assem;
	Aln *aligned, *gap_align;
	Qseqs *qseq, *header;
	AssemInfo *matrix;
	AlnPoints *points;
	NWmat *NWmatrices;
	
	/* get input */
	template = thread->template;
	file_count = thread->file_count;
	files = thread->files;
	xml_out = thread->xml_out;
	frag_out = thread->frag_out;
	aligned_assem = thread->aligned_assem;
	aligned = thread->aligned;
	gap_align = thread->gap_align;
	qseq = thread->qseq;
	header = thread->header;
	matrix = thread->matrix;
	points = thread->points;
	NWmatrices = thread->NWmatrices;
	delta = qseq->size;
	mq = thread->mq;
	minlen = thread->minlen;
	scoreT = thread->scoreT;
	mrc = thread->mrc;
	evalue = thread->evalue;
	bcd = thread->bcd;
	sam = thread->sam;
	ef = thread->ef;
	//spin = thread->spin;
	thread_num = thread->thread_num;
	seq_in = thread->seq_in;
	kmersize = thread->kmersize;
	Wl = points->rewards->Wl;
	
	if(template != -2) {
		wait_atomic(thread_begin);
		
		/* all assemblies done, 
		signal threads to return */
		if(template == -1) {
			lock(excludeMatrix);
			mainTemplate = template;
			unlock(excludeMatrix);
			return NULL;
		}
		
		/* Allocate assembly arrays */
		lock(excludeMatrix);
		template_name = thread->template_name;
		template_index = thread->template_index;
		t_len = thread->t_len;
		load = (template_index->len != 0) ? 0 : 1;
		matrix->len = 0;
		seq_in = thread->seq_in;
		
		/* start threads */
		aligned_assem->score = 0;
		aligned_assem->fragmentCountAln = 0;
		aligned_assem->readCountAln = 0;
		aligned_assem->cover = 0;
		aligned_assem->depth = 0;
		aligned_assem->depthVar = 0;
		aligned_assem->var = 0;
		aligned_assem->nucHighVar = 0;
		aligned_assem->maxDepth = 0;
		aligned_assem->snpSum = 0;
		aligned_assem->insertSum = 0;
		aligned_assem->deletionSum = 0;
		aligned_assem->t[0] = 0;
		aligned_assem->s[0] = 0;
		aligned_assem->q[0] = 0;
		aligned_assem->len = 0;
		aligned_assem->aln_len = 0;
		mainTemplate = template;
		thread_wait = thread_num;
		thread_init = thread_num;
		thread_begin = thread_num;
		unlock(excludeMatrix);
		template = -2;
	}
	
	do {
		while(template == mainTemplate) {
			usleep(100);
		}
		lock(excludeMatrix);
		template = mainTemplate;
		unlock(excludeMatrix);
		if(template == -1) {
			return NULL;
		}
		
		/* load index */
		if(load) {
			hashMapCCI_load_thread(template_index, seq_in, t_len, kmersize, thread_num);
		}
		
		/* init matrix */
		lock(excludeMatrix);
		if(matrix->len == 0) {
			matrix->len = t_len;
			if(matrix->size < (t_len << 1)) {
				matrix->size = (t_len << 1);
				free(matrix->assmb);
				matrix->assmb = smalloc(matrix->size * sizeof(Assembly));
			}
			next = 0;
		}
		unlock(excludeMatrix);
		
		/* init in chunks of 8112 ~= 1 MB*/
		chunk = 8112;
		while(chunk) {
			/* get next chunk */
			lock(excludeMatrix);
			i = next;
			if((next += chunk) < 0) {
				next = t_len;
			}
			unlock(excludeMatrix);
			
			/* fill in chunk */
			if(i < t_len) {
				end = i + chunk;
				if(t_len < end) {
					end = t_len;
				}
				assembly = matrix->assmb + i - 1;
				while(i != end) {
					counts = (++assembly)->counts;
					*counts = 0;
					*++counts = 0;
					*++counts = 0;
					*++counts = 0;
					*++counts = 0;
					*++counts = 0;
					assembly->next = ++i;
				}
				if(t_len <= end) {
					/* circularize */
					assembly->next = 0;
					chunk = 0;
				}
			} else {
				chunk = 0;
			}
		}
		
		/* wait for init to finish */
		lock(excludeIn);
		--thread_init;
		unlock(excludeIn);
		wait_atomic(thread_init);
		
		/* load reads of this template */
		file_i = 0;
		while(file_i < file_count) {
			//lockTime(excludeIn, spin);
			lock(excludeIn);
			file = files[file_i];
			if(file != 0) {
				read_score = fread(buffer, sizeof(int), 8, file);
				if((nextTemplate = buffer[0]) == template) {
					/* load frag */
					qseq->len = buffer[1];
					stats[0] = buffer[2];
					read_score = buffer[3];
					stats[2] = buffer[4];
					stats[3] = buffer[5];
					header->len = buffer[6];
					stats[4] = buffer[7];
					
					if(qseq->size < qseq->len) {
						free(qseq->seq);
						qseq->size = qseq->len << 1;
						qseq->seq = smalloc(qseq->size);
					}
					if(header->size < header->len) {
						header->size = header->len + 1;
						free(header->seq);
						header->seq = smalloc(header->size);
					}
					sfread(qseq->seq, 1, qseq->len, file);
					sfread(header->seq, 1, header->len, file);
					unlock(excludeIn);
					
					if(delta < qseq->len) {
						delta = qseq->len << 1;
						free(aligned->t);
						free(aligned->s);
						free(aligned->q);
						free(gap_align->t);
						free(gap_align->s);
						free(gap_align->q);
						aligned->t = smalloc((delta + 1) << 1);
						aligned->s = smalloc((delta + 1) << 1);
						aligned->q = smalloc((delta + 1) << 1);
						gap_align->t = smalloc((delta + 1) << 1);
						gap_align->s = smalloc((delta + 1) << 1);
						gap_align->q = smalloc((delta + 1) << 1);
					}
					
					/* q-bound */
					if(2 * sizeof(int) + 1 < header->len && header->seq[header->len - 2 * sizeof(int) - 1] == 0) {
						qBoundPtr = (int*) (header->seq + (header->len - 2 * sizeof(int)));
						q_start = *qBoundPtr;
						q_end = *++qBoundPtr;
					} else {
						q_start = 0;
						q_end = qseq->len;
					}
					
					/* Update assembly with read */
					if(read_score || anker_rc(template_index, qseq->seq, qseq->len, q_start, q_end, points)) {
						/* Start with alignment */
						if(stats[3] <= stats[2]) {
							stats[2] = 0;
							stats[3] = t_len;
						}
						
						alnStat = KMA(template_index, qseq->seq, qseq->len, q_start, q_end, aligned, gap_align, stats[2], MIN(t_len, stats[3]), mq, scoreT, points, NWmatrices);
						
						/* get read score */
						aln_len = alnStat.len;
						start = alnStat.pos;
						end = start + aln_len - alnStat.tGaps;
						if(t_len < end) {
							/* circular aln */
							end -= t_len;
						}
						
						/* Get normed score check read coverage */
						read_score = alnStat.score;
						/* full gene award */
						if(qseq->len <= aln_len || t_len <= aln_len) {
							read_score -= Wl;
						} else {
							if(start == 0) {
								read_score -= Wl;
							}
							if(end == t_len) {
								read_score -= Wl;
							}
						}
						alnStat.score = read_score;
						
						/* Get normed score check read coverage */
						if(minlen <= aln_len && mrcheck(mrc, alnStat, qseq->len, t_len)) {
							score = 1.0 * read_score / aln_len;
						} else {
							read_score = 0;
							score = 0;
						}
						
						if(0 < read_score && scoreT <= score) {
							stats[1] = read_score;
							stats[2] = start;
							stats[3] = end;
							if(t_len < end) {
								stats[3] -= t_len;
							}
							/* Update backbone and counts */
							alnToMatPtr(matrix, aligned_assem, aligned, alnStat, t_len, stats[4]);
							
							/* Convert fragment */
							q = qseq->seq;
							i = qseq->len + 1;
							while(--i) {
								*q = bases[*q];
								++q;
							}
							*q = 0;
							
							/* Save fragment */
							if(frag_out) {
								lock(excludeOut);
								updateFrags(frag_out, qseq, header, template_name, stats);
								unlock(excludeOut);
							}
							
							if(sam) {
								header->seq[header->len - 1] = 0;
								samwrite(qseq, header, 0, template_name, aligned, stats);
							}
							
							if(xml_out) {
								hitXML(xml_out, template, header->seq, aligned, &alnStat, NWmatrices->rewards, stats[4]);
							}
						} else if(sam && !(sam & 2096)) {
							stats[1] = read_score;
							stats[2] = start;
							stats[3] = end;
							header->seq[header->len - 1] = 0;
							nibble2base(qseq->seq, qseq->len);
							if(read_score) {
								samwrite(qseq, header, 0, template_name, aligned, stats);
							} else {
								stats[4] |= 4;
								stats[1] = stats[4];
								samwrite(qseq, header, 0, template_name, 0, stats);
							}
						}
					} else if(sam && !(sam & 2096)) {
						stats[4] |= 4;
						stats[1] = stats[4];
						header->seq[header->len - 1] = 0;
						nibble2base(qseq->seq, qseq->len);
						samwrite(qseq, header, 0, template_name, 0, stats);
					}
				} else if(nextTemplate == -1) {
					if(template) {
						fclose(file);
					} else {
						kmaPipe(0, 0, file, &status);
						errno |= status;
					}
					files[file_i] = 0;
					unlock(excludeIn);
					++file_i;
				} else if(nextTemplate < template) {
					/* Move pointer forward */
					sfseek(file, buffer[1] + buffer[6], SEEK_CUR);
					unlock(excludeIn);
				} else {
					/* Move pointer back */
					sfseek(file, (-8) * sizeof(int), SEEK_CUR);
					unlock(excludeIn);
					++file_i;
				}
			} else {
				unlock(excludeIn);
				++file_i;
			}
		}
		/* wait for last reads to align and update */
		lock(excludeIn);
		--thread_wait;
		unlock(excludeIn);
		wait_atomic(thread_wait);
		
		/* collect consensus */
		if(aligned_assem->score) {
			/* Pepare and make alignment on consensus */
			asm_len = matrix->len;
			assembly = matrix->assmb;
			lock(excludeMatrix);
			if(aligned_assem->size <= asm_len) {
				aligned_assem->size = (asm_len + 1) << 1;
				free(aligned_assem->t);
				free(aligned_assem->s);
				free(aligned_assem->q);
				aligned_assem->t = smalloc(aligned_assem->size);
				aligned_assem->s = smalloc(aligned_assem->size);
				aligned_assem->q = smalloc(aligned_assem->size);
			}
			unlock(excludeMatrix);
			
			/* Call nucleotides for the consensus */
			callConsensus(matrix, aligned_assem, template_index->seq, t_len, bcd, evalue, thread_num);
			
			if(ef) {
				getExtendedFeatures(aligned_assem, matrix, template_index->seq, t_len, thread_num);
			}
		} else {
			asm_len = 0;
		}
		lock(excludeMatrix);
		--thread_begin;
		unlock(excludeMatrix);
	} while(thread->num != 0);
	
	/* Terminate alignment on consensus */
	aligned_assem->t[asm_len] = 0;
	aligned_assem->s[asm_len] = 0;
	aligned_assem->q[asm_len] = 0;
	
	return NULL;
}

void complement_swap(short unsigned *c1, short unsigned *c2) {
	
	/* swap counts between two positions and complement them */
	short unsigned tmp;
	
	/* a1 -> t2 */
	tmp = c1[0];
	c1[0] = c2[3];
	c2[3] = tmp;
	
	/* c1 -> g2 */
	tmp = c1[1];
	c1[1] = c2[2];
	c2[2] = tmp;
	
	/* g1 -> c2 */
	tmp = c1[2];
	c1[2] = c2[1];
	c2[1] = tmp;
	
	/* t1 -> a2 */
	tmp = c1[3];
	c1[3] = c2[0];
	c2[0] = tmp;
	
	/* n1 -> n2 */
	tmp = c1[4];
	c1[4] = c2[4];
	c2[4] = tmp;
	
	/* -1 -> -2 */
	tmp = c1[5];
	c1[5] = c2[5];
	c2[5] = tmp;
}

void complement_pos(short unsigned *counts) {
	
	short unsigned tmp;
	
	/* at -> ta */
	tmp = counts[0];
	counts[0] = counts[3];
	counts[3] = tmp;
	
	/* cg -> gc */
	tmp = counts[1];
	counts[1] = counts[2];
	counts[2] = tmp;
}

unsigned rc_insertion(Assembly *assembly, const unsigned new, const unsigned org) {
	
	/* rc insertion ending at 'stop', 
	nodes sorounding insertion will not be updated.
	
	assembly	assembly matrix
	new			new end position of insertion
	org			org end position of insertion
	*/
	unsigned prev, current, next;
	Assembly *ptr;
	
	/* get start of insertion */
	ptr = assembly + (org - 1);
	prev = new;
	current = ptr->next;
	ptr = assembly + current;
	
	/* go over inserted positions */
	while(current != org) {
		/* complent counts */
		complement_pos(ptr->counts);
		
		/* reverse direction */
		next = ptr->next;
		ptr->next = prev;
		
		/* update prev and current, and get next insertion */
		prev = current;
		current = next;
		ptr = assembly + next;
	}
	
	/* return coordinate of last inserted count */
	return prev;
}

void assemble_rc(AssemInfo *matrix, int t_len) {
	
	/* reverse complement assembly matrix */
	int i, j;
	unsigned insert;
	Assembly *assembly;
	
	insert = 0;
	assembly = matrix->assmb;
	for(i = 0, j = t_len - 1; i < j; ++i, --j) {
		/* swap and complement counts */
		complement_swap(assembly[i].counts, assembly[j].counts);
		
		/* check insertions on last position */
		if(i) {
			assembly[i - 1].next = (assembly[j].next < t_len) ? i : rc_insertion(assembly, i, j + 1);
		}
		
		/* update with insertion of last position in previous round */
		assembly[j].next = insert;
		
		/* check first position and save for next round */
		insert = (assembly[i].next < t_len) ? j : rc_insertion(assembly, j, i + 1);
	}
	
	/* check if number of bases are uneven */
	if(i == j) {
		/* complement position */
		complement_pos(assembly[i].counts);
		
		/* check insertion */
		if(i) {
			assembly[i - 1].next = (assembly[j].next < t_len) ? i : rc_insertion(assembly, j, j + 1);
		}
	}
	
	/* update with insertion of last position in previous round */
	assembly[j].next = insert;
}

int assem_rc(Assem *aligned_assem, unsigned char *trans) {
	
	int len, t_len;
	char *s, *s_r;
	unsigned char nuc, *t, *t_r, *q, *q_r;
	
	/* init */
	t_len = 0;
	len = aligned_assem->len;
	t = aligned_assem->t;
	t_r = t-- + len;
	s = aligned_assem->s;
	s_r = s-- + len;
	q = aligned_assem->q;
	q_r = q-- + len;
	len = (len >> 1) + 1;
	
	while(--len) {
		/* t */
		nuc = trans[*++t];
		*t = trans[*--t_r];
		*t_r = nuc;
		t_len += (*t != '-') + (*t_r != '-');
		
		/* s */
		nuc = *++s;
		*s = *--s_r;
		*s_r = nuc;
		
		/* q */
		nuc = trans[*++q];
		*q = trans[*--q_r];
		*q_r = nuc;
	}
	
	if(aligned_assem->len & 1) {
		t_len += (*++t != '-');
		*t = trans[*t];
		++q;
		*q = trans[*q];
	}
	
	return t_len;
}
