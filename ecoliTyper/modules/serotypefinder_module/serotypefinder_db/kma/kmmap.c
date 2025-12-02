/* Philip T.L.C. Clausen May 2019 plan@dtu.dk */

/*
 * Copyright (c) 2019, Philip Clausen, Technical University of Denmark
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
#include <stdio.h>
#include <stdlib.h>
#include "hashmapkma.h"
#include "pherror.h"
#include "stdnuc.h"
#ifdef _WIN32
#define mmap(addr, len, prot, flags, fd, offset) (0; fprintf(stderr, "mmap not available on windows.\n"); exit(1););
#define munmap(addr, len) (-1);
#else
#include <sys/mman.h>
#endif

int hashMapKMAmmap(HashMapKMA *dest, FILE *file) {
	
	int fd;
	unsigned *uptr;
	long unsigned Size, size, *luptr;
	unsigned char *data;
	
	/* mmap data */
	fd = fileno(file);
	sfseek(file, 0, SEEK_END);
	Size = ftell(file);
	data = mmap(0, Size, PROT_READ, MAP_SHARED, fd, 0);
	if(data == MAP_FAILED) {
		ERROR();
	}
	
	/* get data */
	uptr = (unsigned *) data;
	dest->DB_size = *uptr++;
	dest->mlen = *uptr++;
	dest->prefix_len = *uptr++;
	luptr = (long unsigned *) uptr;
	dest->prefix = *luptr++;
	dest->size = *luptr++;
	dest->n = *luptr++;
	dest->v_index = *luptr++;
	dest->null_index = *luptr++;
	data = (unsigned char *) luptr;
	Size -= (3 * sizeof(unsigned) + 5 * sizeof(long unsigned));
	dest->mask = 0;
	dest->mask = (~dest->mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (dest->mlen << 1));
	dest->shmFlag = 16;
	
	/* exist */
	size = dest->size;
	if((dest->size - 1) == dest->mask) {
		if(dest->v_index <= UINT_MAX) {
			size *= sizeof(unsigned);
			getExistPtr = &getExist;
		} else {
			size *= sizeof(long unsigned);
			getExistPtr = &getExistL;
		}
	} else {
		if(dest->n <= UINT_MAX) {
			size *= sizeof(unsigned);
			getExistPtr = &getExist;
		} else {
			size *= sizeof(long unsigned);
			getExistPtr = &getExistL;
		}
	}
	dest->exist = (unsigned *) data;
	dest->exist_l = (long unsigned *) data;
	data += size;
	Size -= size;
	
	/* values */
	size = dest->v_index;
	if(dest->DB_size < USHRT_MAX) {
		size *= sizeof(short unsigned);
		getValuePtr = &getValueS;
		intpos_bin_contaminationPtr = &intpos_bin_contamination_s;
	} else {
		size *= sizeof(unsigned);
		getValuePtr = &getValue;
		intpos_bin_contaminationPtr = &intpos_bin_contamination;
	}
	
	dest->values = (unsigned *) data;
	dest->values_s = (short unsigned *) data;
	dest->shmFlag |= 2;
	data += size;
	Size -= size;
	
	/* check for megaMap */
	if((dest->size - 1) == dest->mask) {
		dest->key_index = 0;
		dest->key_index_l = 0;
		dest->value_index = 0;
		dest->value_index_l = 0;
		hashMap_get = &megaMap_getGlobal;
	} else {
		hashMap_get = &hashMap_getGlobal;
		/* kmers */
		size = dest->n + 1;
		if(dest->mlen <= 16) {
			size *= sizeof(unsigned);
			getKeyPtr = &getKey;
		} else {
			size *= sizeof(long unsigned);
			getKeyPtr = &getKeyL;
		}
		
		dest->key_index = (unsigned *) data;
		dest->key_index_l = (long unsigned *) data;
		dest->shmFlag |= 4;
		data += size;
		Size -= size;
		
		/* value indexes */
		size = dest->n;
		if(dest->v_index < UINT_MAX) {
			size *= sizeof(unsigned);
			getValueIndexPtr = &getValueIndex;
		} else {
			size *= sizeof(long unsigned);
			getValueIndexPtr = &getValueIndexL;
		}
		dest->value_index = (unsigned *) data;
		dest->value_index_l = (long unsigned *) data;
		dest->shmFlag |= 8;
		data += size;
		Size -= size;
		
		/* make indexing a masking problem */
		--dest->size;
	}
	
	if(Size) {
		uptr = (unsigned *) data;
		dest->kmersize = *uptr++;
		dest->flag = *uptr++;
	} else {
		dest->kmersize = dest->mlen;
		/* mark as invalid for later unmapping */
		++dest->mlen;
		dest->flag = 0;
	}
	
	return 0;
}

void hashMapKMA_munmap(HashMapKMA *dest) {
	
	int unit;
	long unsigned size;
	unsigned char *data;
	
	if(dest && dest->shmFlag & 16) {
		data = (unsigned char *) dest->exist;
		size = 3 * sizeof(unsigned) + 5 * sizeof(long unsigned);
		data -= size;
		if((dest->size - 1) == dest->mask) {
			if(dest->v_index <= UINT_MAX) {
				unit = sizeof(unsigned);
			} else {
				unit = sizeof(long unsigned);
			}
		} else {
			if(dest->n <= UINT_MAX) {
				unit = sizeof(unsigned);
			} else {
				unit = sizeof(long unsigned);
			}
		}
		size += dest->size * unit;
		unit = (dest->DB_size < USHRT_MAX) ? sizeof(short unsigned) : sizeof(unsigned);
		size += dest->v_index * unit;
		if((dest->size - 1) != dest->mask) {
			unit = (dest->mlen <= 16) ? sizeof(unsigned) : sizeof(long unsigned);
			size += (dest->n + 1) * unit;
			unit = (dest->v_index < UINT_MAX) ? sizeof(unsigned) : sizeof(long unsigned);
			size += dest->n * unit;
		}
		if(dest->mlen <= dest->kmersize) {
			size += 2 * sizeof(unsigned);
		}
		
		if(data && dest->shmFlag & 1 && munmap(data, size) < 0) {
			ERROR();
		}
		free(dest);
	}
}
