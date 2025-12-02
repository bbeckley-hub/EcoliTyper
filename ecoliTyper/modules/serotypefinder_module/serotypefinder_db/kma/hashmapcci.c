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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hashmapcci.h"
#include "pherror.h"
#include "stdnuc.h"
#include "stdstat.h"
#include "threader.h"
#ifndef _WIN32
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/types.h>
#else
typedef int key_t;
#define ftok(charPtr, integer) (0)
#define shmget(key, size, permission) ((size != 0) ? (-1) : (-key))
#define shmat(shmid, NULL_Ptr, integer) (NULL)
#endif

HashMapCCI * (*alignLoadPtr)(HashMapCCI *, int, int, int, long unsigned) = &alignLoad_fly;

long unsigned hashMapCCI_initialize(HashMapCCI *dest, int len, int kmerindex) {
	
	long unsigned size;
	
	dest->len = len;
	dest->cci_next = 1;
	dest->kmerindex = kmerindex;
	dest->bandwidth = 64;
	dest->ndup = 0;
	/* convert to 2 times nearest power of 2 */
	size = (len << 1) - 1;
	size |= size >> 1;
	size |= size >> 2;
	size |= size >> 4;
	size |= size >> 8;
	size |= size >> 16;
	size |= size >> 32;
	dest->mask = size++;
	if(dest->size < size) {
		if(dest->size) {
			free(dest->index);
			free(dest->seq);
		}
		dest->size = size;
		
		/* add closed chains */
		dest->cci_size = size + 1;//(size >> 1) + (size >> 2) + 2;
		dest->index = calloc(size + dest->cci_size, sizeof(int));
		dest->seq = malloc(((len >> 5) + 1) * sizeof(long unsigned));
		if(!dest->index || !dest->seq) {
			ERROR();
		}
		size = 0;
	} else {
		dest->cci_size = size + 1;//(size >> 1) + (size >> 2) + 2;
		size += dest->cci_size;
	}
	dest->cci_avail = dest->cci_size - 1;
	dest->chain = dest->index + dest->mask; /* 1-indexed */
	
	return size;
}

void hashMapCCI_destroy(HashMapCCI *dest) {
	
	if(dest) {
		if(dest->index) {
			free(dest->index);
		}
		if(dest->seq) {
			free(dest->seq);
		}
		free(dest);
	}
}

int hashMapCCI_get(const HashMapCCI *dest, long unsigned key, unsigned shifter) {
	
	unsigned cPos, iPos;
	int pos, *chain;
	long unsigned index, kmer;
	
	/* get hash */
	murmur(index, key);
	index &= dest->mask;
	
	/* check index */
	if((pos = dest->index[index]) == 0) {
		return 0;
	} else if(0 < pos) {
		getKmer_macro(kmer, dest->seq, (pos - 1), cPos, iPos, shifter);
		return (kmer == key) ? pos : 0;
	}
	
	/* check chain */
	chain = dest->chain - pos - 1;
	while((pos = *++chain)) {
		pos = abs(pos) - 1;
		getKmer_macro(kmer, dest->seq, pos, cPos, iPos, shifter);
		if(kmer == key) {
			return *chain;
		}
	}
	
	return 0;
}

int hashMapCCI_get_bound(const HashMapCCI *dest, long unsigned key, int min, int max, unsigned shifter) {
	
	long unsigned index;
	int pos, apos, *chain;
	
	/* get hash */
	murmur(index, key);
	index &= dest->mask;
	
	/* check index */
	if((pos = dest->index[index]) == 0) {
		return 0;
	} else if(0 < pos) {
		return (min < pos && pos < max && getKmer(dest->seq, pos - 1, shifter) == key) ? pos : 0;
	}
	
	/* check chain */
	chain = dest->chain - pos - 1;
	while((pos = *++chain)) {
		apos = abs(pos);
		if(min < apos && apos < max && getKmer(dest->seq, apos - 1, shifter) == key) {
			return pos;
		}
	}
	
	return 0;
}

int * hashMapCCI_getDubPos(const HashMapCCI *dest, long unsigned key, int value, unsigned shifter) {
	
	long unsigned index;
	int pos, *chain;
	
	/* check input */
	if(0 < value) {
		return 0;
	}
	
	/* get hash */
	murmur(index, key);
	index &= dest->mask;
	
	/* check index */
	if(0 <= (pos = dest->index[index])) {
		return 0;
	}
	
	/* check chain */
	chain = dest->chain - pos - 1;
	while((pos = *++chain)) {
		if(pos == value) {
			return chain;
		}
	}
	
	return 0;
}

int * hashMapCCI_getNextDubPos(const HashMapCCI *dest, int *chain, long unsigned key, int min, int max, unsigned shifter) {
	
	int pos;
	
	min = -min;
	max = -max;
	
	/* check chain */
	while((pos = *++chain)) {
		if(pos < 0 && max < pos && pos < min && getKmer(dest->seq, -pos - 1, shifter) == key) {
			return chain;
		}
	}
	
	return 0;
}

int defragChain(HashMapCCI *dest, int size, int shifter) {
	
	/* defragmentize the chains to make space for new chain, return 0 on failure */
	int newsize, cci_size, pos, apos, newpos, index, fulldefrag;
	int *newchain, *chain, *next;
	unsigned cPos, iPos;
	long unsigned ipos, kmer;
	
	/* check if chain is available */
	if(size <= dest->cci_avail) {
		dest->cci_avail -= size;
		pos = dest->cci_next;
		dest->cci_next += size;
		return pos;
	}
	dest->cci_avail = 0;
	fulldefrag = dest->cci_next == 1;
	
	/* sync cci_next / move to next available slot */
	pos = dest->cci_next;
	chain = dest->chain + pos;
	next = chain--;
	cci_size = dest->cci_size - 1;
	while(pos < cci_size && *chain != *next) {
		++pos;
		chain = next++;
	}
	/* end of chains, start over */
	if(cci_size <= pos) {
		if(fulldefrag) {
			return 0;
		}
		dest->cci_next = 1;
		return defragChain(dest, size, shifter);
	}
	dest->cci_next = pos;
	
	/* defragmentize to make space */
	newsize = 0;
	newpos = pos; /* index of newchain */
	index = pos; /* index of chain */
	newchain = chain++; /* newchain is one behind */
	next = chain + 1;
	/* find new fragment of min size */
	while(index < cci_size && newsize < size) {
		next = chain + 1;
		while(index < cci_size && *chain == *next && newsize < size) {
			++newsize;
			++index;
			chain = next++;
		}
		if(index < cci_size && newsize < size) {
			apos = abs(*next) - 1;
			getKmer_macro(kmer, dest->seq, apos, cPos, iPos, shifter);
			murmur(ipos, kmer);
			ipos &= dest->mask;
			dest->index[ipos] = -newpos;
			chain = next;
			++index;
			while(*chain) {
				*++newchain = *chain;
				*chain++ = 0;
				++newpos;
				++index;
			}
			*++newchain = 0;
			++newpos;
		}
	}
	++newchain;
	dest->cci_next = (index < cci_size) ? index : 1;
	
	/* check if defrag was succesfull */
	if(newsize < size && !fulldefrag) {
		dest->cci_next = 1;
		return defragChain(dest, size, shifter);
	}
	
	return (newsize < size) ? 0 : newpos;
}

int newChain(HashMapCCI *dest, int pos, int newpos, long unsigned kmer, int shifter) {
	
	/* add new chain to hashmap */
	int cPos, iPos, *chain;
	long unsigned index, key;
	
	if(4 <= dest->cci_avail) {
		index = dest->cci_next;
		dest->cci_next += 4;
		dest->cci_avail -= 4;
	} else if(3 <= dest->cci_avail) {
		index = dest->cci_next;
		dest->cci_next += 3;
		dest->cci_avail -= 3;
	} else {
		/* no space left, defrag chains */
		index = defragChain(dest, 3, shifter);
	}
	chain = dest->chain + index;
	
	/* check duplication */
	getKmer_macro(key, dest->seq, (pos - 1), cPos, iPos, shifter);
	if(key == kmer) {
		*chain = -pos;
		*++chain = -newpos;
	} else {
		*chain = pos;
		*++chain = newpos;
	}
	
	return -index;
}

int extendChain(HashMapCCI *dest, int chainpos, int newpos, long unsigned kmer, int shifter) {
	
	int pos, apos, dup, *chain, *newchain;
	unsigned cPos, iPos;
	long unsigned index, key;
	
	pos = chainpos;
	chain = dest->chain + pos;
	dup = 1;
	while(*chain && dup) {
		apos = abs(*chain) - 1;
		getKmer_macro(key, dest->seq, apos, cPos, iPos, shifter);
		if(key == kmer) {
			dup = 0;
			if(0 < *chain) {
				*chain = -*chain;
				dest->ndup++;
			}
		}
		++pos;
		++chain;
	}
	
	/* shuffle duplicates */
	/* count duplicates, and select new dup with uniform possibility */
	if(dup == 0) {
		/* mark as duplicate */
		newpos = -newpos;
		
		/* select duplicate with uniform distribution */
		dest->ndup++;
		murmur(index, dest->ndup);
		index %= dest->ndup;
		
		/* find position of swap */
		while(*chain) {
			/* check k-mer */
			apos = abs(*chain) - 1;
			getKmer_macro(key, dest->seq, apos, cPos, iPos, shifter);
			if(key == kmer && !index--) { /* swap */
				apos = *chain;
				*chain = newpos;
				newpos = apos;
			}
			++pos;
			++chain;
		}
	} else {
		index = 0;
	}
	
	/* check if chain is sufficent */
	if(index) { /* newpos is duplicated and falls outside bandwidth */
		dest->ndup--;
		return chainpos;
	} else if(*++chain == 0) { /* already space enough */
		*--chain = newpos;
		/* move next avail chain */
		if(dest->cci_avail && pos + 1 == dest->cci_next) {
			++dest->cci_next;
			--dest->cci_avail;
		}
	} else if(chainpos != 1 && (dest->chain[chainpos - 2] == 0 || chainpos == 2)) { /* check previous chain */
		/* move chain one back */
		newchain = dest->chain + chainpos;
		chain = newchain--;
		while(*chain) {
			*newchain++ = *chain++;
		}
		*newchain = newpos;
		--chainpos;
		/* threading -> unordered */
		//dest->chain[--chainpos] = newpos;
	} else if((pos = defragChain(dest, (pos - chainpos + 2), shifter))) { /* run defragmentation of chains */
		/* mv to new pos */
		murmur(index, kmer);
		index &= dest->mask;
		chainpos = -(dest->index[index]);
		
		chain = dest->chain + chainpos;
		newchain = dest->chain + pos;
		while(*chain) {
			*newchain++ = *chain;
			*chain++ = 0;
		}
		*newchain = newpos;
		chainpos = pos;
	} else { /* defrag failed -> bias chains by one to the end of chains */
		/* get new chainpos */
		murmur(index, kmer);
		index &= dest->mask;
		chainpos = -(dest->index[index]);
		
		/* go to end of chain */
		chain = dest->chain + chainpos;
		while(*++chain);
		
		/* bias by one to end of chains */
		newchain = chain++;
		while(*newchain != *chain) {
			if((dup = *newchain) == 0) {
				/* add new chain index */
				apos = abs(newpos) - 1;
				getKmer_macro(kmer, dest->seq, apos, cPos, iPos, shifter);
				murmur(index, kmer);
				index &= dest->mask;
				--(dest->index[index]);
			}
			*newchain = newpos;
			newpos = dup;
			newchain = chain++;
		}
	}
	
	return -chainpos;
}

void hashMapCCI_add(HashMapCCI *dest, long unsigned key, int newpos, unsigned shifter) {
	
	int pos, *index_ptr;
	long unsigned index;
	
	if(key == 0) {
		/* likely undefined region */
		return;
	}
	
	/* get hash */
	murmur(index, key);
	index_ptr = dest->index + (index & dest->mask);
	
	/* check index */
	if((pos = *index_ptr) == 0) {
		*index_ptr = newpos;
	} else if(0 < pos) {
		/* new chain */
		*index_ptr = newChain(dest, pos, newpos, key, shifter);
	} else {
		/* extend chain */
		*index_ptr = extendChain(dest, -pos, newpos, key, shifter);
	}
}

void hashMapCCI_add_thread(HashMapCCI *dest, long unsigned key, int newpos, unsigned shifter) {
	
	static volatile int Lock[2] = {0, 0};
	volatile int *indexLock = &Lock[0], *cciLock = &Lock[1];
	int pos, *index_ptr;
	long unsigned index;
	
	if(key == 0) {
		/* likely undefined region */
		return;
	}
	
	/* get hash */
	murmur(index, key);
	index_ptr = dest->index + (index & dest->mask);
	
	/* check index */
	lockTime(indexLock, 1);
	if((pos = *index_ptr) == 0) {
		*index_ptr = newpos;
		unlock(indexLock);
	} else {
		unlock(indexLock);
		/* check chain */
		lockTime(cciLock, 5);
		if((pos = *index_ptr) < 0) {
			/* extend chain */
			*index_ptr = extendChain(dest, -pos, newpos, key, shifter);
		} else {
			*index_ptr = newChain(dest, pos, newpos, key, shifter);
		}
		unlock(cciLock);
	}
}

HashMapCCI * hashMapCCI_load(HashMapCCI *src, int seq, int len, int kmersize) {
	
	int i, end, shifter, cPos, iPos;
	long size;
	long unsigned kmer;
	
	/* init */
	if(src == 0) {
		src = smalloc(sizeof(HashMapCCI));
		src->size = 0;
	}
	if((size = hashMapCCI_initialize(src, len, kmersize))) {
		memset(src->index, 0, size * sizeof(int));
	}
	
	/* get seq */
	size = ((src->len >> 5) + 1) * sizeof(long unsigned);
	if(size != read(seq, src->seq, size)) {
		if(0 < size) {
			fprintf(stderr, "Corrupted *.seq.b\n");
			exit(1);
		} else {
			ERROR();
		}
	}
	
	/* add k-mers */
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (src->kmerindex << 1);
	end = len - kmersize + 1;
	for(i = 0; i < end; ++i) {
		getKmer_macro(kmer, src->seq, i, cPos, iPos, shifter);
		hashMapCCI_add(src, kmer, i + 1, shifter);
	}
	
	return src;
}

HashMapCCI * hashMapCCI_load_thread(HashMapCCI *src, int seq, int len, int kmersize, int thread_num) {
	
	static volatile int Lock = 0, next = 1, thread_wait = 0;
	static long unsigned size;
	volatile int *lock = &Lock;
	int i, end, shifter, chunk, cPos, iPos;
	long check;
	long unsigned kmer;
	
	/* init */
	lock(lock);
	if(src->len == 0) {
		size = hashMapCCI_initialize(src, len, kmersize);
		thread_wait = thread_num;
		next = 0;
		unlock(lock);
		
		/* get seq */
		check = ((src->len >> 5) + 1) * sizeof(long unsigned);
		if(check != read(seq, src->seq, check)) {
			if(check < 0) {
				ERROR();
			} else {
				fprintf(stderr, "Corrupted *.seq.b\n");
				exit(1);
			}
		}
	} else {
		unlock(lock);
	}
	
	/* set in chunks of 16224 */
	chunk = 16224;
	while(chunk) {
		lock(lock);
		i = next;
		if((next += chunk) < 0) {
			next = size;
		}
		unlock(lock);
		
		if(i < size) {
			if(size < i + chunk) {
				chunk = size - i;
			}
			memset(src->index + i, 0, chunk * sizeof(int));
		} else {
			chunk = 0;
		}
	}
	
	/* this version slurps cpu in the locks of hashMapCCI_add_thread */
	/* wait for init and seq to finish */
	/*
	int stop;
	lock(lock);
	if(--thread_wait == 0) {
		next = 0;
	}
	unlock(lock);
	wait_atomic(thread_wait);
	*/
	/* add k-mers */
	/*
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	end = len - kmersize + 1;
	chunk = 16224;
	while(chunk) {
		lock(lock);
		i = next;
		if((next += chunk) < 0) {
			next = end;
		}
		unlock(lock);
		
		if(i < end) {
			stop = (end <= i + chunk) ? end : (i + chunk);
			--i;
			while(++i < stop) {
				hashMapCCI_add_thread(src, getKmer(src->seq, i, shifter), i + 1, shifter);
			}
		} else {
			chunk = 0;
		}
	}
	*/
	
	lock(lock);
	if(--thread_wait == 0) {
		/* add k-mers */
		shifter = sizeof(long unsigned) * sizeof(long unsigned) - (src->kmerindex << 1);
		end = len - kmersize + 1;
		i = -1;
		while(++i < end) {
			getKmer_macro(kmer, src->seq, i, cPos, iPos, shifter);
			hashMapCCI_add(src, kmer, i + 1, shifter);
		}
	}
	unlock(lock);
	
	/* might not be completely loaded here, 
	thus external thread_wait is required */
	return src;
}

void hashMapCCI_dump(HashMapCCI *src, FILE *seq) {
	cfwrite(src->seq, sizeof(long unsigned), (src->len >> 5) + 1, seq);
}

HashMapCCI * alignLoad_fly(HashMapCCI *dest, int seq_in, int len, int kmersize, long unsigned seq_index) {
	
	/* move file pointer */
	lseek(seq_in, seq_index, SEEK_SET);
	
	return hashMapCCI_load(dest, seq_in, len, kmersize);
}

HashMapCCI * alignLoad_fly_mem(HashMapCCI *dest, int seq_in, int len, int kmersize, long unsigned seq_index) {
	
	return hashMapCCI_load(dest, seq_in, len, kmersize);
}

HashMapCCI * alignLoad_skip(HashMapCCI *dest, int seq_in, int len, int kmersize, long unsigned seq_index) {
	
	dest->len = len;
	dest->kmerindex = kmersize;
	
	return dest;
}
