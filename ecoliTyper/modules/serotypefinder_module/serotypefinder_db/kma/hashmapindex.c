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
#include "hashmapindex.h"
#include "pherror.h"
#include "stdnuc.h"
#include "stdstat.h"
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

void (*destroyPtr)(HashMap_index *) = &alignClean;
HashMap_index * (*alignLoadPtr)(HashMap_index *, int, int, int, int, long unsigned, long unsigned) = &alignLoad_fly;

int hashMap_index_initialize(HashMap_index *dest, int len, int kmerindex) {
	
	long unsigned size;
	
	dest->len = len;
	dest->kmerindex = kmerindex;
	/* convert to 2 times nearest power of 2 */
	size = (len << 1) - 1;
	size |= size >> 1;
	size |= size >> 2;
	size |= size >> 4;
	size |= size >> 8;
	size |= size >> 16;
	size |= size >> 32;
	if(dest->size < size) {
		if(dest->size) {
			free(dest->index);
			free(dest->seq);
		}
		dest->size = size++;
		dest->index = calloc(size, sizeof(int));
		dest->seq = malloc(((len >> 5) + 1) * sizeof(long unsigned));
		if(!dest->index || !dest->seq) {
			ERROR();
		}
		
		return 0;
	}
	
	return 1;
}

void hashMap_index_set(HashMap_index *dest) {
	
	memset(dest->index, -1, (dest->size + 1) * sizeof(int));
	memset(dest->seq, 0, ((dest->len >> 5) + 1) * sizeof(long unsigned));
}

void hashMap_index_destroy(HashMap_index *dest) {
	
	if(dest->index) {
		free(dest->index);
	}
	if(dest->seq) {
		free(dest->seq);
	}
}

int hashMap_index_get(const HashMap_index *dest, long unsigned key, unsigned shifter) {
	
	unsigned index;
	int pos;
	
	murmur(index, key);
	index &= dest->size;
	while(index < dest->size && (pos = dest->index[index]) != 0) {
		if(getKmer(dest->seq, abs(pos) - 1, shifter) == key) {
			return pos;
		}
		++index;
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; ++index) {
			if(getKmer(dest->seq, abs(pos) - 1, shifter) == key) {
				return pos;
			}
		}
	}
	
	return 0;
}

int hashMap_index_get_bound(const HashMap_index *dest, long unsigned key, int min, int max, unsigned shifter) {
	
	unsigned index;
	int pos;
	
	murmur(index, key);
	index &= dest->size;
	
	while(index < dest->size && (pos = dest->index[index]) != 0) {
		if(min < abs(pos) && abs(pos) < max && getKmer(dest->seq, abs(pos) - 1, shifter) == key) {
			return pos;
		}
		++index;
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; ++index) {
			if(min < abs(pos) && abs(pos) < max && getKmer(dest->seq, abs(pos) - 1, shifter) == key) {
				return pos;
			}
		}
	}
	
	return 0;
}

int hashMap_index_getDubPos(const HashMap_index *dest, long unsigned key, int value, unsigned shifter) {
	
	unsigned index;
	int pos;
	
	murmur(index, key);
	index &= dest->size;
	while(index < dest->size && (pos = dest->index[index]) != 0) {
		if(pos == value) {
			return index;
		}
		++index;
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; ++index) {
			if(pos == value) {
				return index;
			}
		}
	}
	
	return -1;
}

int hashMap_index_getNextDubPos(const HashMap_index *dest, long unsigned key, int min, int max, unsigned index, unsigned shifter) {
	
	int pos;
	
	min = -min;
	max = -max;
	
	while(++index < dest->size && (pos = dest->index[index]) != 0) {
		if(max < pos && pos < min && getKmer(dest->seq, -pos - 1, shifter) == key) {
			return index;
		}
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; ++index) {
			if(max < pos && pos < min && getKmer(dest->seq, -pos - 1, shifter) == key) {
				return index;
			}
		}
	}
	
	return -1;
}

void hashMap_index_add(HashMap_index *dest, long unsigned key, int newpos, unsigned shifter) {
	
	int pos, neg;
	unsigned index;
	
	if(key == 0) {
		/* likely undefined region */
		return;
	}
	
	neg = 1;
	++newpos;
	murmur(index, key);
	index &= dest->size;
	while(index < dest->size && (pos = dest->index[index]) != 0) {
		if(pos > 0) {
			if(getKmer(dest->seq, pos - 1, shifter) == key) {
				dest->index[index] = -pos;
				neg = -1;
			}
		} else {
			if(getKmer(dest->seq, -pos - 1, shifter) == key) {
				neg = -1;
			}
		}
		++index;
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; ++index) {
			if(pos > 0) {
				if(getKmer(dest->seq, pos - 1, shifter) == key) {
					dest->index[index] = -pos;
					neg = -1;
				}
			} else {
				if(getKmer(dest->seq, -pos - 1, shifter) == key) {
					neg = -1;
				}
			}
		}
	}
	
	if(index < dest->size) {
		dest->index[index] = neg * newpos;
	}
	
}

void hashMapIndex_add(HashMap_index *dest, long unsigned key, int newpos) {
	
	int index, pos, neg, shifter;
	
	neg = 1;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (dest->kmerindex << 1);
	++newpos;
	murmur(index, key);
	index &= dest->size;
	while(index < dest->size && (pos = dest->index[index]) != 0) {
		if(pos > 0) {
			if(getKmer(dest->seq, pos - 1, shifter) == key) {
				dest->index[index] = -pos;
				neg = -1;
			}
		} else {
			if(getKmer(dest->seq, -pos - 1, shifter) == key) {
				neg = -1;
			}
		}
		++index;
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; ++index) {
			if(pos > 0) {
				if(getKmer(dest->seq, pos - 1, shifter) == key) {
					dest->index[index] = -pos;
					neg = -1;
				}
			} else {
				if(getKmer(dest->seq, -pos - 1, shifter) == key) {
					neg = -1;
				}
			}
		}
	}
	
	if(index < dest->size) {
		dest->index[index] = neg * newpos;
	}
	
}

HashMap_index * hashMap_index_load(HashMap_index *src, int seq, int index, int len, int kmersize) {
	
	if(src == 0) {
		src = smalloc(sizeof(HashMap_index));
		src->size = 0;
	}
	if(hashMap_index_initialize(src, len, kmersize)) {
		src->size = len << 1;
	}
	
	read(seq, src->seq, ((src->len >> 5) + 1) * sizeof(long unsigned));
	read(index, src->index, src->size * sizeof(int));
	
	return src;
}

void hashMap_index_dump(HashMap_index *src, FILE *seq, FILE *index) {
	cfwrite(src->seq, sizeof(long unsigned), (src->len >> 5) + 1, seq);
	cfwrite(src->index, sizeof(int), src->size, index);
}

HashMap_index * hashMap_index_build(HashMap_index *src, int seq, int len, int kmersize) {
	
	int i, end, shifter;
	
	if(src == 0) {
		src = smalloc(sizeof(HashMap_index));
		src->size = 0;
	}
	
	if(hashMap_index_initialize(src, len, kmersize)) {
		memset(src->index, 0, src->size * sizeof(int));
	}
	
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (src->kmerindex << 1);
	
	read(seq, src->seq, ((src->len >> 5) + 1) * sizeof(long unsigned));
	
	end = len - kmersize + 1;
	for(i = 0; i < end; ++i) {
		hashMap_index_add(src, getKmer(src->seq, i, shifter), i, shifter);
	}
	
	return src;
}

HashMap_index * alignLoad_fly(HashMap_index *dest, int seq_in, int index_in, int len, int kmersize, long unsigned seq_index, long unsigned index_index) {
	
	/* move file pointer */
	lseek(index_in, index_index, SEEK_SET);
	lseek(seq_in, seq_index, SEEK_SET);
	
	return hashMap_index_load(dest, seq_in, index_in, len, kmersize);
}

HashMap_index * alignLoad_fly_mem(HashMap_index *dest, int seq_in, int index_in, int len, int kmersize, long unsigned seq_index, long unsigned index_index) {
	
	return hashMap_index_load(dest, seq_in, index_in, len, kmersize);
}

HashMap_index * alignLoad_fly_build(HashMap_index *dest, int seq_in, int index_in, int len, int kmersize, long unsigned seq_index, long unsigned index_index) {
	
	/* move file pointer */
	lseek(seq_in, seq_index, SEEK_SET);
	
	return hashMap_index_build(dest, seq_in, len, kmersize);
}

HashMap_index * alignLoad_fly_build_mem(HashMap_index *dest, int seq_in, int index_in, int len, int kmersize, long unsigned seq_index, long unsigned index_index) {
	
	return hashMap_index_build(dest, seq_in, len, kmersize);
}

HashMap_index * alignLoad_skip(HashMap_index *dest, int seq_in, int index_in, int len, int kmersize, long unsigned seq_index, long unsigned index_index) {
	
	dest->len = len;
	dest->kmerindex = kmersize;
	
	return dest;
}

HashMap_index * alignLoad_fly_shm(HashMap_index *dest, int seq_in, int index_in, int len, int kmersize, long unsigned seq_index, long unsigned index_index) {
	
	static HashMap_index *src = 0;
	
	if(kmersize == 0) {
		if(src == 0) {
			src = smalloc(sizeof(HashMap_index));
		}
		src->seq = dest->seq;
		src->index = dest->index;
	} else {
		dest->len = len;
		dest->size = len << 1;
		dest->kmerindex = kmersize;
		dest->seq = src->seq + (seq_index / sizeof(long unsigned));
		dest->index = src->index + ((index_index - sizeof(int)) / sizeof(int));
	}
	
	return dest;
}

HashMap_index * alignLoad_shm_initial(char *templatefilename, int file_len, int seq_in, int index_in, int kmersize) {
	
	key_t key;
	int shmid;
	HashMap_index *dest;
	long unsigned size;
	
	dest = smalloc(sizeof(HashMap_index));
	dest->len = 0;
	dest->size = 0;
	dest->kmerindex = kmersize;
	
	templatefilename[file_len] = 0;
	strcat(templatefilename, ".index.b");
	key = ftok(templatefilename, 'i');
	read(index_in, &kmersize, sizeof(int));
	size = lseek(index_in, 0, SEEK_END) - sizeof(int);
	shmid = shmget(key, size, 0666);
	if(shmid < 0) {
		fprintf(stderr, "SHM error i.\n");
		exit(1);
	} else {
		dest->index = shmat(shmid, NULL, 0);
	}
	
	templatefilename[file_len] = 0;
	strcat(templatefilename, ".seq.b");
	key = ftok(templatefilename, 's');
	size = lseek(seq_in, 0, SEEK_END);
	shmid = shmget(key, size, 0666);
	if(shmid < 0) {
		fprintf(stderr, "SHM error s.\n");
		exit(1);
	} else {
		dest->seq = shmat(shmid, NULL, 0);
	}
	lseek(seq_in, 0, SEEK_SET);
	templatefilename[file_len] = 0;
	
	return alignLoad_fly_shm(dest, 0, 0, 0, 0, 0, 0);
}

void alignClean(HashMap_index *template_index) {
	/* Overwrite later */
	/*
	if(template_index) {
		hashMap_index_destroy(template_index);
	}
	*/
}

void alignClean_shm(HashMap_index *template_index) {
	
	if(template_index) {
		template_index->seq = 0;
		template_index->index = 0;
	}
}
