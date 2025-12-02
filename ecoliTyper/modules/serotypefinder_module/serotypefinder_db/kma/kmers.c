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
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include "ankers.h"
#include "compdna.h"
#include "hashmapkma.h"
#include "kmapipe.h"
#include "kmeranker.h"
#include "kmers.h"
#include "kmmap.h"
#include "penalties.h"
#include "pherror.h"
#include "qseqs.h"
#include "savekmers.h"
#include "spltdb.h"
#ifndef _WIN32
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/types.h>
#else
typedef int key_t;
#define ftok(charPtr, integer) (0)
#define shmget(key, size, permission) ((size != 0) ? (-1) : (-key))
#define shmat(shmid, NULL_Ptr, integer) (NULL)
#define shmdt(dest) fprintf(stderr, "sysV not available on Windows.\n")
#define shmctl(shmid, cmd, buf) fprintf(stderr, "sysV not available on Windows.\n")
#endif

int save_kmers_batch(char *templatefilename, char *exePrev, unsigned shm, int thread_num, const int exhaustive, Penalties *rewards, FILE *out, int sam, int minlen, double mrs, double coverT, double minFrac) {
	
	int i, file_len, shmid, deCon, *bestTemplates, *template_lengths;
	long unsigned *softProxi;
	FILE *inputfile, *templatefile;
	time_t t0, t1;
	key_t key;
	HashMapKMA *templates;
	KmerScan_thread *threads, *thread;
	
	/* open pipe */
	//inputfile = popen(exePrev, "r");
	inputfile = kmaPipe(exePrev, "rb", 0, 0);
	if(!inputfile) {
		ERROR();
	}
	t0 = clock();
	
	/* do not output unmapped sam reads */
	if(sam != 1 || out == stdout) {
		sam = 0;
	}
	
	/* load hashMap */
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
		if(hashMapKMA_load(templates, templatefile, templatefilename) == 1) {
			fprintf(stderr, "Wrong format of DB.\n");
			exit(1);
		}
	}
	templatefilename[file_len] = 0;
	fclose(templatefile);
	
	/* check if DB is sparse */
	if(templates->prefix_len != 0 || templates->prefix != 0) {
		/* set pointers to sparse detection */
		if(printPtr != &print_ankers_spltDB) {
			printPtr = &print_ankers_Sparse;
		} else {
			printPtr = &print_ankers_Sparse_spltDB;
		}
		if(deConPrintPtr != &deConPrint) {
			deConPrintPtr = printPtr;
		}
		if(templates->prefix_len == 0 && get_kmers_for_pair_ptr != &get_kmers_for_pair_count) {
			if(kmerScan == &save_kmers) {
				kmerScan = &save_kmers_pseuodeSparse;
			} else {
				kmerScan = &save_kmers_sparse_chain;
			}
			//kmerScan = &save_kmers_pseuodeSparse;
			get_kmers_for_pair_ptr = &get_kmers_for_pair_pseoudoSparse;
		} else {
			if(kmerScan == &save_kmers) {
				kmerScan = &save_kmers_Sparse;
			} else {
				kmerScan = &save_kmers_sparse_chain;
			}
			//kmerScan = &save_kmers_Sparse;
			get_kmers_for_pair_ptr = &get_kmers_for_pair_Sparse;
		}
	}
	
	/* allocate scoring arrays */
	if(printPtr == &print_ankers_spltDB || printPtr == &print_ankers_Sparse_spltDB) {
		printPtr(0, 0, thread_num, 0, 0, 0);
	}
	/* proxi */
	softProxi = 0;
	if(minFrac < 1.0) {
		if(minFrac < 0) {
			minFrac = -minFrac;
			softProxi = calloc(templates->DB_size + 3, sizeof(long unsigned));
			if(!softProxi) {
				ERROR();
			}
		}
		getMatch((int *)(&minFrac), 0);
		getMatch((int *)(softProxi), 0);
		getMatchSparse((int *)(&minFrac), 0, 0, 0, 0, 0);
		getMatchSparse((int *)(softProxi), 0, 0, 0, 0, 0);
		getSecondForce((int *)(&minFrac), (int *)(softProxi), 0, 0, 0, 0);
		getSecondPen((int *)(&minFrac), (int *)(softProxi), 0, 0, 0, 0, 0, 0);
		getF((int *)(&minFrac), (int *)(softProxi), 0, 0, 0);
		ankerAndClean((int *)(&minFrac), (int *)(softProxi), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
		ankerAndClean_MEM((int *)(&minFrac), (int *)(softProxi), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
		getProxiChainTemplates(0, (const Penalties *)(&minFrac), (int *)(softProxi), 0, 0, 0, 0, 0, 0, 0);
	}
	template_lengths = 0;
	if(kmerScan == &save_kmers_HMM || kmerScan == &save_kmers_chain || kmerScan == &save_kmers_sparse_chain) {
		/* load lengths */
		strcat(templatefilename, ".length.b");
		templatefile = sfopen(templatefilename, "rb");
		sfread(&templates->DB_size, sizeof(int), 1, templatefile);
		if(shm & 4) {
			key = ftok(templatefilename, 'l');
			shmid = shmget(key, templates->DB_size * sizeof(int), 0666);
			if(shmid < 0) {
				fprintf(stderr, "No shared length\n");
				exit(1);
			} else {
				template_lengths = shmat(shmid, NULL, 0);
			}
		} else {
			template_lengths = smalloc(templates->DB_size * sizeof(int));
			sfread(template_lengths, sizeof(int), templates->DB_size, templatefile);
		}
		templatefilename[file_len] = 0;
		fclose(templatefile);
		
		if(kmerScan == &save_kmers_HMM) {
			save_kmers_HMM(templates, 0, &(int){thread_num}, template_lengths, 0, 0, 0, 0, 0, 0, minlen, 0, 0);
		} else {
			kmerScan(0, 0, &(int){thread_num}, (int *)(&coverT), (int *)(&mrs), template_lengths, 0, 0, 0, 0, minlen, 0, 0);
		}
	}
	
	t1 = clock();
	fprintf(stderr, "#\n# Total time used for DB loading: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	t0 = clock();
	fprintf(stderr, "# Finding k-mer ankers\n");
	
	/* initialize threads */
	i = 1;
	threads = 0;
	while(i < thread_num) {
		
		thread = smalloc(sizeof(KmerScan_thread));
		thread->num = i;
		thread->exhaustive = exhaustive;
		thread->sam = sam;
		thread->bestScore = 0;
		thread->bestScore_r = 0;
		thread->bestTemplates = calloc((templates->DB_size << 1) + 4, sizeof(int));
		thread->bestTemplates_r = calloc((templates->DB_size << 1) + 4, sizeof(int));
		if(!thread->bestTemplates || !thread->bestTemplates_r) {
			ERROR();
		}
		thread->templates = templates;
		thread->inputfile = inputfile;
		thread->rewards = rewards;
		thread->out = out;
		thread->next = threads;
		threads = thread;
		
		/* start thread */
		if((errno = pthread_create(&thread->id, NULL, &save_kmers_threaded, thread))) {
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
	thread = smalloc(sizeof(KmerScan_thread));
	thread->num = 0;
	thread->bestScore = 0;
	thread->bestScore_r = 0;
	bestTemplates = calloc((templates->DB_size << 1) + 4, sizeof(int));
	thread->bestTemplates = bestTemplates;
	thread->bestTemplates_r = calloc(templates->DB_size + 4, sizeof(int));
	if(!bestTemplates || !thread->bestTemplates_r) {
		ERROR();
	}
	thread->templates = templates;
	thread->inputfile = inputfile;
	thread->rewards = rewards;
	thread->exhaustive = exhaustive;
	thread->sam = sam;
	thread->out = out;
	
	/* start k-mer search */
	save_kmers_threaded(thread);
	
	/* join threads */
	for(thread = threads; thread; thread = thread->next) {
		/* join thread */
		if((errno = pthread_join(thread->id, NULL))) {
			ERROR();
		} else if(bestTemplates[2] < thread->bestTemplates[2]) {
			bestTemplates[2] = thread->bestTemplates[2];
		}
	}
	
	/* print remaining buffer */
	if(printPtr == &print_ankers_spltDB || printPtr == &print_ankers_Sparse_spltDB) {
		printPtr(bestTemplates, 0, 0, 0, 0, out);
	} else {
		/* print number of fragments and send terminating signal*/
		sfwrite(&(int){-bestTemplates[2]}, sizeof(int), 1, out);
	}
	if(softProxi) {
		sfwrite(softProxi, sizeof(int), 6, out);
		sfwrite(softProxi, sizeof(long unsigned), templates->DB_size, out);
	}
	
	kmaPipe(0, 0, inputfile, &i);
	if(kmaPipe == &kmaPipeFork) {
		t1 = clock();
		fprintf(stderr, "#\n# Total time used ankering query: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	} else {
		fprintf(stderr, "# Query ankered\n#\n");
	}
	
	/* clean up */
	if(!((shm & 1) || (deCon && (shm & 2)))) {
		hashMapKMA_destroy(templates);
	}
	if(template_lengths && !(shm & 4)) {
		free(template_lengths);
	}
	if(kmerScan == &save_kmers_chain || kmerScan == &save_kmers_sparse_chain) {
		kmerScan(0, 0, &(int){thread_num}, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	}
	for(thread = threads; thread; thread = threads) {
		threads = thread->next;
		free(thread->bestTemplates);
		free(thread->bestTemplates_r);
		free(thread);
	}
	
	return i;
}
