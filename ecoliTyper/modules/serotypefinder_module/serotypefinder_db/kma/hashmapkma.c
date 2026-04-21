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
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "hashmapkma.h"
#include "pherror.h"
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
#include <sys/types.h>
#endif

void (*hashMapKMA_destroy)(HashMapKMA *) = &hashMapKMA_free;
long unsigned (*getExistPtr)(const unsigned *, const long unsigned);
long unsigned (*getKeyPtr)(const unsigned *, const long unsigned);
long unsigned (*getValueIndexPtr)(const unsigned *, const long unsigned);
unsigned * (*getValuePtr)(const HashMapKMA *, const long unsigned);
unsigned * (*hashMap_get)(const HashMapKMA *, const long unsigned) = &hashMap_getGlobal;
int (*intpos_bin_contaminationPtr)(const unsigned *, const int);
int (*getSizePtr)(const unsigned *);
void (*hashMapKMA_addKey_ptr)(HashMapKMA *, long unsigned, long unsigned);
void (*hashMapKMA_addValue_ptr)(HashMapKMA *, long unsigned, long unsigned);
void (*hashMapKMA_addExist_ptr)(HashMapKMA *, long unsigned, long unsigned);

long unsigned getExist(const unsigned *exist, const long unsigned pos) {
	return exist[pos];
}

long unsigned getExistL(const unsigned *exist, const long unsigned pos) {
	return ((long unsigned *) exist)[pos];
}

long unsigned getKey(const unsigned *key_index, const long unsigned pos) {
	return key_index[pos];
}

long unsigned getKeyL(const unsigned *key_index, const long unsigned pos) {
	return ((long unsigned *) key_index)[pos];
}

long unsigned getValueIndex(const unsigned *value_index, const long unsigned pos) {
	return value_index[pos];
}

long unsigned getValueIndexL(const unsigned *value_index, const long unsigned pos) {
	return ((long unsigned *) value_index)[pos];
}

unsigned * getValue(const HashMapKMA *dest, const long unsigned pos) {
	return (dest->values + pos);
}

unsigned * getValueS(const HashMapKMA *dest, const long unsigned pos) {
	return (unsigned *)(dest->values_s + pos);
}

int getSize(const unsigned *values) {
	return *values * sizeof(unsigned) + sizeof(unsigned);
}

int getSizeS(const unsigned *values) {
	return *((short unsigned *)(values)) * sizeof(short unsigned) + sizeof(short unsigned);
}

int intpos_bin_contamination(const unsigned *str1, const int str2) {
	
	int pos, upLim, downLim, template;
	
	upLim = *str1;
	if(upLim == 0) {
		return -1;
	}
	
	downLim = 1;
	pos = (upLim + downLim) >> 1;
	while(0 < (upLim - downLim)) {
		template = str1[pos];
		if(template == str2) {
			return pos;
		} else if(template < str2) {
			downLim = pos + 1;
		} else {
			upLim = pos - 1;
		}
		pos = (upLim + downLim) >> 1;
	}
	if(str1[pos] == str2) {
		return pos;
	}
	
	return -1;
}

int intpos_bin_contamination_s(const unsigned *Str1, const int str2) {
	
	int pos, upLim, downLim, template;
	short unsigned *str1 = (short unsigned *) Str1;
	
	upLim = *str1;
	if(upLim == 0) {
		return -1;
	}
	
	downLim = 1;
	pos = (upLim + downLim) >> 1;
	while(0 < (upLim - downLim)) {
		template = str1[pos];
		if(template == str2) {
			return pos;
		} else if(template < str2) {
			downLim = pos + 1;
		} else {
			upLim = pos - 1;
		}
		pos = (upLim + downLim) >> 1;
	}
	if(str1[pos] == str2) {
		return pos;
	}
	
	return -1;
}

unsigned * hashMap_getGlobal(const HashMapKMA *templates, const long unsigned key) {
	
	unsigned flag;
	long unsigned pos, kpos, kmer;
	
	if((flag = templates->flag)) {
		murmur(kpos, key);
		kpos &= templates->size;
	} else {
		kpos = key & templates->size;
	}
	
	pos = getExistPtr(templates->exist, kpos);
	
	if(pos != templates->null_index) {
		kmer = getKeyPtr(templates->key_index, pos);
		while(key != kmer) {
			if(flag) {
				murmur(kmer, kmer);
			}
			if(kpos != (kmer & templates->size)) {
				return 0;
			}
			kmer = getKeyPtr(templates->key_index, ++pos);
		}
		return getValuePtr(templates, getValueIndexPtr(templates->value_index, pos));
	}
	
	return 0;
}

void loadPrefix(HashMapKMA *dest, FILE *file) {
	
	long unsigned Size, size;
	
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
	dest->exist = 0;
	dest->exist_l = 0;
	dest->values = 0;
	dest->values_s = 0;
	dest->key_index = 0;
	dest->key_index_l = 0;
	dest->value_index = 0;
	dest->value_index_l = 0;
	
	/* get index of kmersize and flag */
	Size = 3 * sizeof(unsigned) + 5 * sizeof(long unsigned);
	/* exist */
	size = dest->size;
	if((dest->size - 1) == dest->mask) {
		if(dest->v_index <= UINT_MAX) {
			size *= sizeof(unsigned);
		} else {
			size *= sizeof(long unsigned);
		}
	} else {
		if(dest->n <= UINT_MAX) {
			size *= sizeof(unsigned);
		} else {
			size *= sizeof(long unsigned);
		}
	}
	Size += size;
	
	/* values */
	size = dest->v_index;
	if(dest->DB_size < USHRT_MAX) {
		size *= sizeof(short unsigned);
	} else {
		size *= sizeof(unsigned);
	}
	Size += size;
	
	if((dest->size - 1) != dest->mask) {
		/* kmers */
		size = (dest->n + 1);
		if(dest->mlen <= 16) {
			size *= sizeof(unsigned);
		} else {
			size *= sizeof(long unsigned);
		}
		Size += size;
		/* value indexes */
		size = dest->n;
		if(dest->v_index < UINT_MAX) {
			size *= sizeof(unsigned);
		} else {
			size *= sizeof(long unsigned);
		}
		Size += size;
	}
	sfseek(file, Size, SEEK_SET);
	
	if(fread(&dest->kmersize, sizeof(unsigned), 1, file)) {
		sfread(&dest->flag, sizeof(unsigned), 1, file);
		setCmerPointers(dest->flag);
	} else {
		dest->kmersize = dest->mlen;
		dest->flag = 0;
	}
}

unsigned * megaMap_getGlobal(const HashMapKMA *templates, const long unsigned key) {
	
	long unsigned pos;
	
	if((pos = getExistPtr(templates->exist, key & templates->mask)) != 1) {
		return getValuePtr(templates, pos);
	}
	
	return 0;
}

int hashMapKMA_load(HashMapKMA *dest, FILE *file, const char *filename) {
	
	key_t key;
	int shmid;
	long unsigned check, size, seekSize;
	
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
		return 1;
	}
	
	/* check shared memory, else load */
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
	key = ftok(filename, 'e');
	shmid = shmget(key, size, 0666);
	if(shmid < 0) {
		/* not shared, load */
		errno = 0;
		dest->exist = smalloc(size);
		check = fread(dest->exist, 1, size, file);
		if(check != size) {
			return 1;
		}
		seekSize = 0;
		dest->shmFlag |= 1;
	} else {
		/* found */
		dest->exist = shmat(shmid, NULL, 0);
		seekSize = size;
	}
	dest->exist_l = (long unsigned *)(dest->exist);
	
	/* values */
	size = dest->v_index;
	if(dest->DB_size < USHRT_MAX) {
		size *= sizeof(short unsigned);
		getValuePtr =&getValueS;
		intpos_bin_contaminationPtr = &intpos_bin_contamination_s;
	} else {
		size *= sizeof(unsigned);
		getValuePtr =&getValue;
		intpos_bin_contaminationPtr = &intpos_bin_contamination;
	}
	key = ftok(filename, 'v');
	shmid = shmget(key, size, 0666);
	if(shmid < 0) {
		/* not shared, load */
		errno = 0;
		dest->values = smalloc(size);
		if(seekSize) {
			sfseek(file, seekSize, SEEK_CUR);
		}
		seekSize = 0;
		check = fread(dest->values, 1, size, file);
		if(check != size) {
			return 1;
		}
		dest->shmFlag |= 2;
	} else {
		/* found */
		dest->values = shmat(shmid, NULL, 0);
		seekSize += size;
	}
	dest->values_s = (short unsigned *)(dest->values);
	
	/* check for megaMap */
	if((dest->size - 1) == dest->mask) {
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
		key = ftok(filename, 'k');
		shmid = shmget(key, size, 0666);
		if(shmid < 0) {
			/* not shared, load */
			errno = 0;
			dest->key_index = smalloc(size);
			if(seekSize) {
				sfseek(file, seekSize, SEEK_CUR);
			}
			seekSize = 0;
			check = fread(dest->key_index, 1, size, file);
			if(check != size) {
				return 1;
			}
			dest->shmFlag |= 4;
		} else {
			/* found */
			dest->key_index = shmat(shmid, NULL, 0);
			seekSize += size;
		}
		dest->key_index_l = (long unsigned *)(dest->key_index);
		
		/* value indexes */
		size = dest->n;
		if(dest->v_index < UINT_MAX) {
			size *= sizeof(unsigned);
			getValueIndexPtr = &getValueIndex;
		} else {
			size *= sizeof(long unsigned);
			getValueIndexPtr = &getValueIndexL;
		}
		key = ftok(filename, 'i');
		shmid = shmget(key, size, 0666);
		if(shmid < 0) {
			/* not shared, load */
			errno = 0;
			dest->value_index = smalloc(size);
			if(seekSize) {
				sfseek(file, seekSize, SEEK_CUR);
			}
			seekSize = 0;
			check = fread(dest->value_index, 1, size, file);
			if(check != size) {
				return 1;
			}
			dest->shmFlag |= 8;
		} else {
			/* found */
			dest->value_index = shmat(shmid, NULL, 0);
			seekSize += size;
		}
		dest->value_index_l = (long unsigned *)(dest->value_index);
		
		/* make indexing a masking problem */
		--dest->size;
	}
	
	if(seekSize) {
		sfseek(file, seekSize, SEEK_CUR);
		seekSize = 0;
	}
	if(fread(&dest->kmersize, sizeof(unsigned), 1, file)) {
		sfread(&dest->flag, sizeof(unsigned), 1, file);
		setCmerPointers(dest->flag);
	} else {
		dest->kmersize = dest->mlen;
		dest->flag = 0;
	}
	
	return 0;
}

void hashMapKMA_load_shm(HashMapKMA *dest, FILE *file, const char *filename) {
	
	key_t key;
	int shmid;
	long unsigned size, seekSize;
	
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
	
	/* check shared memory */
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
	key = ftok(filename, 'e');
	shmid = shmget(key, size, 0666);
	seekSize = size;
	if(shmid < 0) {
		/* not shared */
		fprintf(stderr, "DB e not shared, see kma_shm\n");
		exit(1);
	} else {
		/* found */
		dest->exist = shmat(shmid, NULL, 0);
	}
	dest->exist_l = (long unsigned *)(dest->exist);
	
	/* values */
	size = dest->v_index;
	if(dest->DB_size < USHRT_MAX) {
		size *= sizeof(short unsigned);
		getValuePtr =&getValueS;
		intpos_bin_contaminationPtr = &intpos_bin_contamination_s;
	} else {
		size *= sizeof(unsigned);
		getValuePtr =&getValue;
		intpos_bin_contaminationPtr = &intpos_bin_contamination;
	}
	key = ftok(filename, 'v');
	shmid = shmget(key, size, 0666);
	seekSize += size;
	if(shmid < 0) {
		/* not shared */
		fprintf(stderr, "DB v not shared, see kma_shm\n");
		exit(1);
	} else {
		/* found */
		dest->values = shmat(shmid, NULL, 0);
	}
	dest->values_s = (short unsigned *)(dest->values);
	
	/* check for megaMap */
	if((dest->size - 1) == dest->mask) {
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
		key = ftok(filename, 'k');
		shmid = shmget(key, size, 0666);
		seekSize += size;
		if(shmid < 0) {
			/* not shared */
			fprintf(stderr, "DB k not shared, see kma_shm\n");
			exit(1);
		} else {
			/* found */
			dest->key_index = shmat(shmid, NULL, 0);
		}
		dest->key_index_l = (long unsigned *)(dest->key_index);
		
		/* value indexes */
		size = dest->n;
		if(dest->v_index < UINT_MAX) {
			size *= sizeof(unsigned);
			getValueIndexPtr = &getValueIndex;
		} else {
			size *= sizeof(long unsigned);
			getValueIndexPtr = &getValueIndexL;
		}
		key = ftok(filename, 'i');
		shmid = shmget(key, size, 0666);
		seekSize += size;
		if(shmid < 0) {
			/* not shared */
			fprintf(stderr, "DB i not shared, see kma_shm\n");
			exit(1);
		} else {
			/* found */
			dest->value_index = shmat(shmid, NULL, 0);
		}
		dest->value_index_l = (long unsigned *)(dest->value_index);
		
		/* make indexing a masking problem */
		--dest->size;
	}
	
	if(seekSize) {
		sfseek(file, seekSize, SEEK_CUR);
		seekSize = 0;
	}
	if(fread(&dest->kmersize, sizeof(unsigned), 1, file)) {
		sfread(&dest->flag, sizeof(unsigned), 1, file);
		setCmerPointers(dest->flag);
	} else {
		dest->kmersize = dest->mlen;
		dest->flag = 0;
	}
}

int hashMapKMAload(HashMapKMA *dest, FILE *file) {
	
	long unsigned check, size;
	
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
	
	/* exist */
	size = dest->size;
	if((dest->size - 1) == dest->mask) {
		if(dest->v_index <= UINT_MAX) {
			size *= sizeof(unsigned);
			getExistPtr = &getExist;
			hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
		} else {
			size *= sizeof(long unsigned);
			getExistPtr = &getExistL;
			hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
		}
	} else {
		if(dest->n <= UINT_MAX) {
			size *= sizeof(unsigned);
			getExistPtr = &getExist;
			hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
		} else {
			size *= sizeof(long unsigned);
			getExistPtr = &getExistL;
			hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
		}
	}
	dest->exist = smalloc(size);
	check = fread(dest->exist, 1, size, file);
	if(check != size) {
		return 1;
	}
	dest->exist_l = (long unsigned *)(dest->exist);
	dest->shmFlag |= 1;
	
	/* values */
	size = dest->v_index;
	if(dest->DB_size < USHRT_MAX) {
		size *= sizeof(short unsigned);
		getValuePtr = &getValueS;
		getSizePtr = &getSizeS;
	} else {
		size *= sizeof(unsigned);
		getValuePtr = &getValue;
		getSizePtr = &getSize;
	}
	dest->values = smalloc(size);
	check = fread(dest->values, 1, size, file);
	if(check != size) {
		return 1;
	}
	dest->values_s = (short unsigned *)(dest->values);
	dest->shmFlag |= 2;
	
	/* check for megaMap */
	if((dest->size - 1) == dest->mask) {
		dest->key_index = 0;
		dest->key_index_l = 0;
		dest->value_index = 0;
		dest->value_index_l = 0;
	} else {
		/* kmers */
		size = dest->n + 1;
		if(dest->mlen <= 16) {
			size *= sizeof(unsigned);
			getKeyPtr = &getKey;
		} else {
			size *= sizeof(long unsigned);
			getKeyPtr = &getKeyL;
		}
		dest->key_index = smalloc(size);
		check = fread(dest->key_index, 1, size, file);
		if(check != size) {
			return 1;
		}
		dest->key_index_l = (long unsigned *)(dest->key_index);
		dest->shmFlag |= 4;
		
		/* value indexes */
		size = dest->n;
		if(dest->v_index < UINT_MAX) {
			size *= sizeof(unsigned);
			hashMapKMA_addValue_ptr = &hashMapKMA_addValue;
			getValueIndexPtr = &getValueIndex;
		} else {
			size *= sizeof(long unsigned);
			hashMapKMA_addValue_ptr = &hashMapKMA_addValueL;
			getValueIndexPtr = &getValueIndexL;
		}
		dest->value_index = smalloc(size);
		check = fread(dest->value_index, 1, size, file);
		if(check != size) {
			return 1;
		}
		dest->value_index_l = (long unsigned *)(dest->value_index);
		dest->shmFlag |= 8;
	}
	
	if(fread(&dest->kmersize, sizeof(unsigned), 1, file)) {
		sfread(&dest->flag, sizeof(unsigned), 1, file);
		setCmerPointers(dest->flag);
	} else {
		dest->kmersize = dest->mlen;
		dest->flag = 0;
	}
	
	return 0;
}

void hashMapKMA_dump(HashMapKMA *dest, FILE *out) {
	
	/* dump sizes */
	cfwrite(&dest->DB_size, sizeof(unsigned), 1, out);
	cfwrite(&dest->mlen, sizeof(unsigned), 1, out);
	cfwrite(&dest->prefix_len, sizeof(unsigned), 1, out);
	cfwrite(&dest->prefix, sizeof(long unsigned), 1, out);
	cfwrite(&dest->size, sizeof(long unsigned), 1, out);
	cfwrite(&dest->n, sizeof(long unsigned), 1, out);
	cfwrite(&dest->v_index, sizeof(long unsigned), 1, out);
	cfwrite(&dest->null_index, sizeof(long unsigned), 1, out);
	
	/* dump arrays */
	if(dest->n <= UINT_MAX) {
		cfwrite(dest->exist, sizeof(unsigned), dest->size, out);
		getExistPtr = &getExist;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
	} else {
		cfwrite(dest->exist_l, sizeof(long unsigned), dest->size, out);
		getExistPtr = &getExistL;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
	}
	
	if(dest->DB_size < USHRT_MAX) {
		cfwrite(dest->values_s, sizeof(short unsigned), dest->v_index, out);
		getValuePtr = &getValueS;
		getSizePtr = &getSizeS;
	} else {
		cfwrite(dest->values, sizeof(unsigned), dest->v_index, out);
		getValuePtr = &getValue;
		getSizePtr = &getSize;
	}
	
	if(dest->mlen <= 16) {
		cfwrite(dest->key_index, sizeof(unsigned), dest->n + 1, out);
		getKeyPtr = &getKey;
	} else {
		cfwrite(dest->key_index_l, sizeof(long unsigned), dest->n + 1, out);
		getKeyPtr = &getKeyL;
	}
	
	if(dest->v_index < UINT_MAX) {
		cfwrite(dest->value_index, sizeof(unsigned), dest->n, out);
		hashMapKMA_addValue_ptr = &hashMapKMA_addValue;
		getValueIndexPtr = &getValueIndex;
	} else {
		cfwrite(dest->value_index_l, sizeof(long unsigned), dest->n, out);
		hashMapKMA_addValue_ptr = &hashMapKMA_addValueL;
		getValueIndexPtr = &getValueIndexL;
	}
	
	cfwrite(&dest->kmersize, sizeof(unsigned), 1, out);
	cfwrite(&dest->flag, sizeof(unsigned), 1, out);
}

void megaMapKMA_dump(HashMapKMA *dest, FILE *out) {
	
	/* dump sizes */
	cfwrite(&dest->DB_size, sizeof(unsigned), 1, out);
	cfwrite(&dest->mlen, sizeof(unsigned), 1, out);
	cfwrite(&dest->prefix_len, sizeof(unsigned), 1, out);
	cfwrite(&dest->prefix, sizeof(long unsigned), 1, out);
	cfwrite(&dest->size, sizeof(long unsigned), 1, out);
	cfwrite(&dest->n, sizeof(long unsigned), 1, out);
	cfwrite(&dest->v_index, sizeof(long unsigned), 1, out);
	cfwrite(&dest->null_index, sizeof(long unsigned), 1, out);
	
	/* dump arrays */
	if(dest->v_index <= UINT_MAX) {
		cfwrite(dest->exist, sizeof(unsigned), dest->size, out);
		getExistPtr = &getExist;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
	} else {
		cfwrite(dest->exist_l, sizeof(long unsigned), dest->size, out);
		getExistPtr = &getExistL;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
	}
	
	if(dest->DB_size < USHRT_MAX) {
		cfwrite(dest->values_s, sizeof(short unsigned), dest->v_index, out);
		getValuePtr = &getValueS;
		getSizePtr = &getSizeS;
	} else {
		cfwrite(dest->values, sizeof(unsigned), dest->v_index, out);
		getValuePtr = &getValue;
		getSizePtr = &getSize;
	}
	
	cfwrite(&dest->kmersize, sizeof(unsigned), 1, out);
	cfwrite(&dest->flag, sizeof(unsigned), 1, out);
}

void hashMapKMA_addKey(HashMapKMA *dest, long unsigned index, long unsigned key) {
	dest->key_index[index] = key;
}

void hashMapKMA_addKeyL(HashMapKMA *dest, long unsigned index, long unsigned key) {
	dest->key_index_l[index] = key;
}

void hashMapKMA_addValue(HashMapKMA *dest, long unsigned index, long unsigned v_index) {
	dest->value_index[index] = v_index;
}

void hashMapKMA_addValueL(HashMapKMA *dest, long unsigned index, long unsigned v_index) {
	dest->value_index_l[index] = v_index;
}

void hashMapKMA_addExist(HashMapKMA *dest, long unsigned index, long unsigned relative) {
	dest->exist[index] = relative;
}

void hashMapKMA_addExistL(HashMapKMA *dest, long unsigned index, long unsigned relative) {
	dest->exist_l[index] = relative;
}

void hashMapKMA_free(HashMapKMA *dest) {
	
	if(dest) {
		if(dest->exist && dest->shmFlag & 1) {
			free(dest->exist);
		}
		if(dest->values && dest->shmFlag & 2) {
			free(dest->values);
		}
		if(dest->key_index && dest->shmFlag & 4) {
			free(dest->key_index);
		}
		if(dest->value_index && dest->shmFlag & 8) {
			free(dest->value_index);
		}
		free(dest);
	}
}
