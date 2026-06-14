/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
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
//#include <sys/param.h>
#include "compress.h"
#include "hashmap.h"
#include "hashmapkma.h"
#include "pherror.h"
#include "valueshash.h"
#include "stdstat.h"
#ifdef _WIN32
#define mmap(addr, len, prot, flags, fd, offset) (0; fprintf(stderr, "mmap not available on windows.\n"); exit(1););
#define munmap(addr, len) (-1);
#else
#include <sys/mman.h>
#endif

unsigned (*valuesSize)(unsigned *);

static void mmapinit(HashMapKMA *finalDB, long unsigned size, FILE *out) {
	
	unsigned char *data;
	unsigned *udata;
	
	/* allocate file */
	sfseek(out, size - 1, SEEK_SET);
	putc(0, out);
	rewind(out);
	
	/* map file */
	data = mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(out), 0);
	finalDB->exist = (unsigned *) data;
	if(finalDB->exist == MAP_FAILED) {
		ERROR();
	}
	posix_madvise(finalDB->exist, size, POSIX_MADV_SEQUENTIAL);
	*finalDB->exist++ = finalDB->DB_size;
	*finalDB->exist++ = finalDB->mlen;
	*finalDB->exist++ = finalDB->prefix_len;
	finalDB->exist_l = (long unsigned *) finalDB->exist;
	*finalDB->exist_l++ = finalDB->prefix;
	*finalDB->exist_l++ = finalDB->size;
	*finalDB->exist_l++ = finalDB->n;
	*finalDB->exist_l++ = finalDB->v_index;
	*finalDB->exist_l++ = finalDB->null_index;
	finalDB->exist = (unsigned *) finalDB->exist_l;
	
	data += size - 2 * sizeof(unsigned);
	udata = (unsigned *) data;
	*udata++ = finalDB->kmersize;
	*udata++ = finalDB->flag;
}

static void rmemcpy(unsigned char *dst, unsigned char *src, size_t n) {
	
	dst += n;
	src += n;
	++n;
	while(--n) {
		*--dst = *--src;
	}
}

HashMapKMA * compressKMA_DB(HashMap *templates, FILE *out) {
	
	unsigned char *data;
	long unsigned i, j, check, size, v_size;
	long unsigned index, t_index, v_index, new_index, null_index;
	unsigned flag, swap, *values, *finalV;
	short unsigned *values_s, *finalV_s;
	HashMapKMA *finalDB;
	ValuesHash *shmValues;
	ValuesTable *node, *next, *table;
	HashTable *node_t, *next_t, *table_t;
	
	/* convert templates to linked list */
	table_t = 0;
	i = templates->size + 1;
	while(i--) {
		for(node_t = templates->table[i]; node_t != 0; node_t = next_t) {
			next_t = node_t->next;
			node_t->next = table_t;
			table_t = node_t;
		}
	}
	free(templates->table);
	templates->table = 0;
	
	/* prepare final DB */
	swap = 0;
	v_size = 0;
	size = 5 * sizeof(unsigned) + 5 * sizeof(long unsigned);
	check = 0;
	check = ~check;
	check >>= 32;
	fprintf(stderr, "# Preparing compressed DB.\n");
	finalDB = smalloc(sizeof(HashMapKMA));
	/* Fill in known values */
	finalDB->size = templates->size + 1;
	finalDB->n = templates->n;
	finalDB->mask = templates->mask;
	finalDB->prefix_len = templates->prefix_len;
	finalDB->prefix = templates->prefix;
	finalDB->kmersize = templates->kmersize;
	finalDB->DB_size = templates->DB_size;
	finalDB->mlen = templates->mlen;
	flag = (finalDB->flag = templates->flag);
	
	/* allocate existence */
	if(finalDB->n <= check) {
		finalDB->exist = malloc(finalDB->size * sizeof(unsigned));
		finalDB->exist_l = 0;
		size += finalDB->size * sizeof(unsigned);
		hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
	} else {
		finalDB->exist = 0;
		finalDB->exist_l = malloc(finalDB->size * sizeof(long unsigned));
		size += finalDB->size * sizeof(long unsigned);
		hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
	}
	
	if(finalDB->mlen <= 16) {
		finalDB->key_index = malloc((finalDB->n + 1) * sizeof(unsigned));
		finalDB->key_index_l = 0;
		size += (finalDB->n + 1) * sizeof(unsigned);
		hashMapKMA_addKey_ptr = &hashMapKMA_addKey;
	} else {
		finalDB->key_index = 0;
		finalDB->key_index_l = malloc((finalDB->n + 1) * sizeof(long unsigned));
		size += (finalDB->n + 1) * sizeof(long unsigned);
		hashMapKMA_addKey_ptr = &hashMapKMA_addKeyL;
	}
	finalDB->value_index = malloc(finalDB->n * sizeof(unsigned));
	size += finalDB->n * sizeof(unsigned);
	
	if((!finalDB->exist && !finalDB->exist_l) || (!finalDB->key_index && !finalDB->key_index_l) || !finalDB->value_index) {
		free(finalDB->exist);
		free(finalDB->exist_l);
		free(finalDB->key_index);
		free(finalDB->key_index_l);
		free(finalDB->value_index);
		
		/* fail over to swap */
		fprintf(stderr, "# Fail over to swap.\n");
		swap = 1;
		
		/* map file */
		mmapinit(finalDB, size, out);
		if(finalDB->n <= UINT_MAX) {
			finalDB->exist_l = 0;
			finalDB->key_index = (finalDB->exist + finalDB->size);
		} else {
			finalDB->exist = 0;
			finalDB->key_index = (unsigned *) (finalDB->exist_l + finalDB->size);
		}
		if(16 < finalDB->mlen) {
			finalDB->key_index_l = (long unsigned *) finalDB->key_index;
			finalDB->key_index = 0;
			finalDB->value_index = (unsigned *) (finalDB->key_index_l + (finalDB->n + 1));
		} else {
			finalDB->value_index = finalDB->key_index + (finalDB->n + 1);
		}
	}
	
	null_index = finalDB->n;
	finalDB->null_index = null_index;
	/* fill with null_indexes */
	i = finalDB->size;
	while(i--) {
		hashMapKMA_addExist_ptr(finalDB, i, null_index);
	}
	
	/* get relative indexes */
	fprintf(stderr, "# Calculating relative indexes.\n");
	hashMapKMA_addValue_ptr = &hashMapKMA_addValue;
	node_t = table_t;
	--finalDB->size;
	shmValues = initialize_hashValues(null_index, finalDB->DB_size);
	t_index = 0;
	v_index = 0;
	while(node_t != 0) {
		/* get index */
		if(flag) {
			murmur(index, node_t->key);
			index &= finalDB->size;
		} else {
			index = (node_t->key & finalDB->size);
		}
		hashMapKMA_addExist_ptr(finalDB, index, t_index);
		
		/* mv chain */
		do {
			next_t = node_t->next;
			
			/* add kmer */
			hashMapKMA_addKey_ptr(finalDB, t_index, node_t->key);
			
			/* the actual value index */
			new_index = valuesHash_add(shmValues, node_t->value, v_index);
			
			if(new_index == v_index) {
				v_index += valuesSize(node_t->value);
				if(check <= v_index) {
					fprintf(stderr, "# Compression overflow.\n");
					check = 0;
					check = ~check;
					hashMapKMA_addValue_ptr = &hashMapKMA_addValueL;
					getValueIndexPtr = &getValueIndexL;
					++finalDB->size;
					if(swap) {
						/* sync and unmap */
						if(finalDB->exist) {
							finalDB->exist_l = (long unsigned *) finalDB->exist;
						}
						finalDB->exist_l -= 5;
						finalDB->exist = (unsigned *) finalDB->exist_l;
						finalDB->exist -= 3;
						msync(finalDB->exist, size, MS_ASYNC);
						munmap(finalDB->exist, size);
						
						/* prep new size */
						size += finalDB->n * (sizeof(long unsigned) - sizeof(unsigned));
						sfseek(out, size, SEEK_SET);
						fputc(0, out);
						rewind(out);
						
						/* map file */
						finalDB->exist = mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(out), 0);
						if(finalDB->exist == MAP_FAILED) {
							ERROR();
						}
						posix_madvise(finalDB->exist, size, POSIX_MADV_SEQUENTIAL);
						finalDB->exist_l = (long unsigned *) (finalDB->exist + 3);
						finalDB->exist_l += 5;
						
						/* set pointers */
						if(finalDB->n <= check) {
							finalDB->exist_l = 0;
							finalDB->key_index = (finalDB->exist + finalDB->size);
						} else {
							finalDB->exist = 0;
							finalDB->key_index = (unsigned *) (finalDB->exist_l + finalDB->size);
						}
						finalDB->value_index = finalDB->key_index + (finalDB->n + 1);
						finalDB->value_index_l = (long unsigned *) finalDB->value_index;
						if(16 < finalDB->mlen) {
							finalDB->key_index_l = (long unsigned *) finalDB->key_index;
							finalDB->key_index = 0;
						}
					} else {
						finalDB->value_index_l = realloc(finalDB->value_index, finalDB->n * sizeof(long unsigned));
						if(!finalDB->value_index_l) {
							fprintf(stderr, "# Fail over to swap.\n");
							/* map file */
							size += finalDB->n * (sizeof(long unsigned) - sizeof(unsigned));
							values = finalDB->exist ? finalDB->exist : (unsigned *) finalDB->exist_l;
							mmapinit(finalDB, size, out);
							
							if(finalDB->n <= check) {
								memcpy(finalDB->exist, values, finalDB->size * sizeof(unsigned));
								free(values);
								values = finalDB->key_index;
								finalDB->exist_l = 0;
								finalDB->key_index = (finalDB->exist + finalDB->size);
							} else {
								memcpy(finalDB->exist_l, values, finalDB->size * sizeof(long unsigned));
								free(values);
								values = finalDB->key_index;
								finalDB->exist = 0;
								finalDB->key_index = (unsigned *) (finalDB->exist_l + finalDB->size);
							}
							
							if(16 < finalDB->mlen) {
								finalDB->key_index_l = (long unsigned *) finalDB->key_index;
								finalDB->key_index = 0;
								memcpy(finalDB->key_index_l, values, (finalDB->n + 1) * sizeof(long unsigned));
								finalDB->value_index_l = finalDB->key_index_l + (finalDB->n + 1);
							} else {
								memcpy(finalDB->key_index, values, (finalDB->n + 1) * sizeof(unsigned));
								finalDB->value_index_l = (long unsigned *) (finalDB->key_index + (finalDB->n + 1));
							}
							free(values);
							swap = 2;
						}
					}
					--finalDB->size;
					finalDB->value_index = (unsigned *)(finalDB->value_index_l);
					j = finalDB->n;
					while(j--) {
						finalDB->value_index_l[j] = finalDB->value_index[j];
					}
					if(swap == 2) {
						free(finalDB->value_index);
						swap = 1;
					}
					finalDB->value_index = 0;
				}
			} else {
				/* values were duplicated, clean up */
				free(node_t->value);
			}
			
			hashMapKMA_addValue_ptr(finalDB, t_index, new_index);
			++t_index;
			
			/* clean */
			free(node_t);
			if((node_t = next_t)) {
				if(flag) {
					murmur(new_index, node_t->key);
					new_index &= finalDB->size;
				} else {
					new_index = node_t->key & finalDB->size;
				}
			} else {
				/* break chain */
				new_index = ~index;
			}
		} while(new_index == index);
	}
	/* convert valuesHash to a linked list */
	table = 0;
	i = shmValues->size;
	while(i--) {
		for(node = shmValues->table[i]; node != 0; node = next) {
			next = node->next;
			node->next = table;
			table = node;
		}
	}
	free(shmValues->table);
	
	/* make compressed values */
	fprintf(stderr, "# Finalizing indexes.\n");
	finalDB->v_index = v_index;
	if(!swap) {
		if(finalDB->DB_size < USHRT_MAX) {
			finalDB->values = 0;
			finalDB->values_s = calloc(v_index, sizeof(short unsigned));
			size += v_index * sizeof(short unsigned);
			if(!finalDB->values_s) {
				swap = 1;
			}
		} else {
			finalDB->values = calloc(v_index, sizeof(unsigned));
			finalDB->values_s = 0;
			size += v_index * sizeof(unsigned);
			if(!finalDB->values) {
				swap = 1;
			}
		}
		
		/* dump part of DB */
		if(swap) {
			++finalDB->size;
			fprintf(stderr, "# Fail over to swap.\n");
			if(finalDB->exist) {
				values = finalDB->exist;
			} else {
				values = (unsigned *) finalDB->exist_l;
			}
			
			/* map file */
			mmapinit(finalDB, size, out);
			
			/* dump arrays */
			if(finalDB->n <= UINT_MAX) {
				memcpy(finalDB->exist, values, finalDB->size * sizeof(unsigned));
				data = (unsigned char *)(finalDB->exist + finalDB->size);
				getExistPtr = &getExist;
				hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
			} else {
				memcpy(finalDB->exist_l, values, finalDB->size * sizeof(long unsigned));
				data = (unsigned char *)(finalDB->exist_l + finalDB->size);
				getExistPtr = &getExistL;
				hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
			}
			free(values);
			
			/* skip values */
			if(finalDB->DB_size < USHRT_MAX) {
				data += finalDB->v_index * sizeof(short unsigned);
				getValuePtr = &getValueS;
				getSizePtr = &getSizeS;
			} else {
				data += finalDB->v_index * sizeof(unsigned);
				getValuePtr = &getValue;
				getSizePtr = &getSize;
			}
			
			if(finalDB->mlen <= 16) {
				memcpy(data, finalDB->key_index, (finalDB->n + 1) * sizeof(unsigned));
				data += (finalDB->n + 1) * sizeof(unsigned);
				free(finalDB->key_index);
				getKeyPtr = &getKey;
			} else {
				memcpy(data, finalDB->key_index_l, (finalDB->n + 1) * sizeof(long unsigned));
				data += (finalDB->n + 1) * sizeof(long unsigned);
				free(finalDB->key_index_l);
				getKeyPtr = &getKeyL;
			}
			
			if(finalDB->v_index < UINT_MAX) {
				memcpy(data, finalDB->value_index, finalDB->n * sizeof(unsigned));
				free(finalDB->value_index);
				hashMapKMA_addValue_ptr = &hashMapKMA_addValue;
				getValueIndexPtr = &getValueIndex;
			} else {
				memcpy(data, finalDB->value_index_l, finalDB->n * sizeof(long unsigned));
				free(finalDB->value_index_l);
				hashMapKMA_addValue_ptr = &hashMapKMA_addValueL;
				getValueIndexPtr = &getValueIndexL;
			}
			--finalDB->size;
		}
	} else {
		/* sync and unmap */
		if(finalDB->exist) {
			finalDB->exist_l = (long unsigned *) finalDB->exist;
		}
		finalDB->exist_l -= 5;
		finalDB->exist = (unsigned *) finalDB->exist_l;
		finalDB->exist -= 3;
		msync(finalDB->exist, size, MS_ASYNC);
		munmap(finalDB->exist, size);
		
		/* make new mapping */
		if(finalDB->DB_size < USHRT_MAX) {
			size += v_index * sizeof(short unsigned);
		} else {
			size += v_index * sizeof(unsigned);
		}
		++finalDB->size;
		mmapinit(finalDB, size, out);
		--finalDB->size;
		swap = 2;
	}
	
	/* set pointers */
	if(swap) {
		++finalDB->size;
		if(finalDB->n <= UINT_MAX) {
			finalDB->values = finalDB->exist + finalDB->size;
		} else {
			finalDB->values = (unsigned *) (finalDB->exist_l + finalDB->size);
		}
		if(finalDB->DB_size < USHRT_MAX) {
			finalDB->values_s = (short unsigned *) finalDB->values;
			finalDB->key_index = (unsigned *) (finalDB->values_s + finalDB->v_index);
		} else {
			finalDB->key_index = finalDB->values + finalDB->v_index;
		}
		if(finalDB->mlen <= 16) {
			finalDB->value_index = finalDB->key_index + (finalDB->n + 1);
		} else {
			finalDB->key_index_l = (long unsigned *) finalDB->key_index;
			finalDB->value_index = (unsigned *) (finalDB->key_index_l + (finalDB->n + 1));
		}
		if(UINT_MAX <= finalDB->v_index) {
			finalDB->value_index_l = (long unsigned *) finalDB->value_index;
		}
		--finalDB->size;
	}
	
	/* move pieces to correct location, rmemcpy */
	if(swap == 2) {
		/* value_index */
		if(finalDB->mlen <= 16) {
			values = finalDB->values + (finalDB->n + 1);
		} else {
			values = (unsigned *) (((long unsigned *) finalDB->values) + (finalDB->n + 1));
		}
		v_size = finalDB->n * (UINT_MAX <= finalDB->v_index ? sizeof(long unsigned) : sizeof(unsigned));
		rmemcpy((unsigned char *) finalDB->value_index, (unsigned char *) values, v_size);
		
		/* key_index */
		v_size = (finalDB->n + 1) * (finalDB->mlen <= 16 ? sizeof(unsigned) : sizeof(long unsigned));
		rmemcpy((unsigned char *) finalDB->key_index, (unsigned char *) finalDB->values, v_size);
		
		swap = 1;
	}
	
	/* move values */
	if(finalDB->DB_size < USHRT_MAX) {
		for(node = table; node != 0; node = next) {
			next = node->next;
			values_s = (short unsigned *)(node->values);
			finalV_s = finalDB->values_s + node->v_index - 1;
			i = 2 + *values_s--;
			while(--i) {
				*++finalV_s = *++values_s;
			}
			free(node->values);
			free(node);
			
			/*
			next = node->next;
			values_s = (short unsigned *)(node->values);
			for(i = node->v_index, j = 0; j <= *values_s; ++i, ++j) {
				finalDB->values_s[i] = values_s[j];
			}
			free(values_s);
			free(node);
			*/
		}
	} else {
		for(node = table; node != 0; node = next) {
			next = node->next;
			values = node->values;
			finalV = finalDB->values + node->v_index - 1;
			i = 2 + *values--;
			while(--i) {
				*++finalV = *++values;
			}
			free(node->values);
			free(node);
			
			/*
			next = node->next;
			values = node->values;
			for(i = node->v_index, j = 0; j <= *values; ++i, ++j) {
				finalDB->values[i] = values[j];
			}
			free(values);
			free(node);
			*/
		}
	}
	
	/* add terminating key */
	i = 0;
	if(finalDB->flag) {
		if(finalDB->mlen <= 16) { 
			murmur(j, finalDB->key_index[finalDB->n - 1]);
			j &= finalDB->size;
			do {
				murmur(index, finalDB->key_index[i]);
				index &= finalDB->size;
				++i;
			} while(j == index);
			finalDB->key_index[finalDB->n] = finalDB->key_index[i];
		} else {
			murmur(j, finalDB->key_index_l[finalDB->n - 1]);
			j &= finalDB->size;
			do {
				murmur(index, finalDB->key_index_l[i]);
				index &= finalDB->size;
				++i;
			} while(j == index);
			finalDB->key_index_l[finalDB->n] = finalDB->key_index_l[i];
		}
	} else {
		if(finalDB->mlen <= 16) { 
			j = finalDB->key_index[finalDB->n - 1] & finalDB->size;
			while(j == (finalDB->key_index[i] & finalDB->size)) {
				++i;
			}
			finalDB->key_index[finalDB->n] = finalDB->key_index[i];
		} else {
			j = finalDB->key_index_l[finalDB->n - 1] & finalDB->size;
			while(j == (finalDB->key_index_l[i] & finalDB->size)) {
				++i;
			}
			finalDB->key_index_l[finalDB->n] = finalDB->key_index_l[i];
		}
	}
	
	/* dump final DB */
	fprintf(stderr, "# Dumping compressed DB\n");	
	++finalDB->size;
	
	if(swap) {
		/* save prev mapping */
		if(finalDB->exist) {
			finalDB->exist_l = (long unsigned *) finalDB->exist;
		}
		*--finalDB->exist_l = finalDB->null_index;
		*--finalDB->exist_l = finalDB->v_index;
		*--finalDB->exist_l = finalDB->n;
		*--finalDB->exist_l = finalDB->size;
		*--finalDB->exist_l = finalDB->prefix;
		finalDB->exist = (unsigned *) finalDB->exist_l;
		*--finalDB->exist = finalDB->prefix_len;
		*--finalDB->exist = finalDB->mlen;
		*--finalDB->exist = finalDB->DB_size;
		msync(finalDB->exist, size, MS_ASYNC);
		munmap(finalDB->exist, size);
		free(finalDB);
		finalDB = 0;
	} else {
		hashMapKMA_dump(finalDB, out);
	}
	
	return finalDB;
}

HashMapKMA * compressKMA_megaDB(HashMap *templates, FILE *out) {
	
	long unsigned i, j, v_index, new_index, null_index, size, v_size;
	unsigned check, swap, *values, *finalV;
	short unsigned *values_s, *finalV_s;
	HashMapKMA *finalDB;
	ValuesHash *shmValues;
	ValuesTable *node, *next, *table;
	
	/* Fill in known values */
	hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
	swap = 0;
	size = 0;
	v_size = 0;
	check = 0;
	check = ~check;
	finalDB = smalloc(sizeof(HashMapKMA));
	finalDB->size = templates->size + 1;
	finalDB->n = templates->n;
	finalDB->mask = templates->mask;
	finalDB->prefix_len = templates->prefix_len;
	finalDB->prefix = templates->prefix;
	finalDB->kmersize = templates->kmersize;
	finalDB->DB_size = templates->DB_size;
	finalDB->mlen = templates->mlen;
	finalDB->flag = templates->flag;
	
	/* allocate existence */
	finalDB->exist = malloc(finalDB->size * sizeof(unsigned));
	if(!finalDB->exist) {
		/* fail over to swap */
		fprintf(stderr, "# Fail over to swap.\n");
		swap = 1;
		
		/* map file */
		size = 5 * sizeof(unsigned) + 5 * sizeof(long unsigned) + finalDB->size * sizeof(unsigned);
		mmapinit(finalDB, size, out);
	}
	finalDB->exist_l = 0;
	finalDB->key_index = 0;
	finalDB->value_index = 0;
	
	/* get relative indexes */
	fprintf(stderr, "# Calculating relative indexes.\n");
	null_index = finalDB->n;
	v_index = 0;
	while(v_index != finalDB->size && templates->values[v_index] != 0) {
		finalDB->exist[v_index] = v_index;
		++v_index;
	}
	for(i = v_index; i != finalDB->size; ++i) {
		if(templates->values[i]) {
			finalDB->exist[i] = v_index;
			templates->values[v_index] = templates->values[i];
			templates->values[i] = 0;
			++v_index;
		} else {
			finalDB->exist[i] = null_index;
		}
	}
	
	/* decrease size of values to what is actually used */
	templates->values = realloc(templates->values, templates->n * sizeof(unsigned *));
	if(!templates->values) {
		ERROR();
	}
	
	/* get compressed indexes */
	fprintf(stderr, "# Compressing indexes.\n");
	v_index = 0;
	shmValues = initialize_hashValues(null_index, finalDB->DB_size);
	i = finalDB->size;
	j = 0;
	while(i--) {
		if(finalDB->exist[i] != null_index) {
			values = templates->values[finalDB->exist[i]];
			
			/* the actual index */
			new_index = valuesHash_add(shmValues, values, v_index);
			
			
			if(new_index == v_index) {
				v_index += valuesSize(values);
				if(check < v_index) {
					fprintf(stderr, "# Compression overflow.\n");
					j = 1;
					hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
					break;
				}
			} else {
				/* values were duplicated, clean up */
				free(values);
				templates->values[finalDB->exist[i]] = 0;
			}
			/* update to new index */
			finalDB->exist[i] = new_index;
		} else {
			finalDB->exist[i] = 1;
		}
	}
	if(j) {
		fprintf(stderr, "# Bypassing overflow.\n");
		if(swap & 1) {
			/* save prev mapping */
			finalDB->exist_l = (long unsigned *) finalDB->exist;
			finalDB->exist_l -= 5;
			finalDB->exist = (unsigned *) finalDB->exist_l;
			finalDB->exist -= 3;
			msync(finalDB->exist, size, MS_ASYNC);
			munmap(finalDB->exist, size);
			
			/* allocate file */
			size = 5 * sizeof(unsigned) + 5 * sizeof(long unsigned) + finalDB->size * sizeof(long unsigned);
			sfseek(out, size, SEEK_SET);
			fputc(0, out);
			rewind(out);
			
			/* map file */
			finalDB->exist = mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(out), 0);
			if(finalDB->exist == MAP_FAILED) {
				ERROR();
			}
			posix_madvise(finalDB->exist, size, POSIX_MADV_SEQUENTIAL);
			finalDB->exist_l = (long unsigned *) (finalDB->exist + 3);
			finalDB->exist_l += 5;
		} else {
			finalDB->exist_l = realloc(finalDB->exist, finalDB->size * sizeof(long unsigned));
			if(!finalDB->exist_l) {
				/* swap */
				fprintf(stderr, "Fail over to swap.\n");
				swap |= 1;
				/* map file */
				size = 5 * sizeof(unsigned) + 5 * sizeof(long unsigned) + finalDB->size * sizeof(long unsigned);
				values = finalDB->exist;
				mmapinit(finalDB, size, out);
				memcpy(finalDB->exist_l, values, finalDB->size * sizeof(unsigned));
				free(finalDB->exist);
			}
		}
		finalDB->exist = (unsigned *)(finalDB->exist_l);
		j = finalDB->size;
		while(j--) {
			finalDB->exist_l[j] = finalDB->exist[j];
		}
		finalDB->exist = 0;
		finalDB->exist_l[i] = new_index;
		
		while(i--) {
			if(finalDB->exist_l[i] != null_index) {
				values = templates->values[finalDB->exist_l[i]];
				
				/* the actual index */
				new_index = valuesHash_add(shmValues, values, v_index);
				
				
				if(new_index == v_index) {
					v_index += valuesSize(values);
				} else {
					/* values were duplicated, clean up */
					free(values);
					templates->values[finalDB->exist_l[i]] = 0;
				}
				/* update to new index */
				finalDB->exist_l[i] = new_index;
				
			} else {
				finalDB->exist_l[i] = 1;
			}
		}
		fprintf(stderr, "# Overflow bypassed.\n");
	}
	free(templates->values);
	
	/* convert valuesHash to a linked list */
	table = 0;
	i = shmValues->size;
	while(i--) {
		for(node = shmValues->table[i]; node != 0; node = next) {
			next = node->next;
			node->next = table;
			table = node;
		}
	}
	free(shmValues->table);
	
	/* make compressed values */
	fprintf(stderr, "# Finalizing indexes.\n");
	finalDB->v_index = v_index;
	finalDB->null_index = 1;
	
	/* try to allocate first, else gradually fail over over to swap */
	finalDB->values = 0;
	finalDB->values_s = 0;
	while(finalDB->values == 0 && finalDB->values_s == 0 && (swap & 4) == 0) {
		if(finalDB->DB_size < USHRT_MAX) {
			finalDB->values = 0;
			finalDB->values_s = calloc(v_index, sizeof(short unsigned));
			if(!finalDB->values_s) {
				swap += 2;
			}
		} else {
			finalDB->values = calloc(v_index, sizeof(unsigned));
			finalDB->values_s = 0;
			if(!finalDB->values) {
				swap += 2;
			}
		}
		
		if(swap & 2) {
			/* free exist */
			size = 5 * sizeof(unsigned) + 5 * sizeof(long unsigned);
			if(swap & 1) {
				/* save prev mapping */
				finalDB->exist_l = (long unsigned *) finalDB->exist;
				*--finalDB->exist_l = finalDB->null_index;
				*--finalDB->exist_l = finalDB->v_index;
				*--finalDB->exist_l = finalDB->n;
				*--finalDB->exist_l = finalDB->size;
				*--finalDB->exist_l = finalDB->prefix;
				finalDB->exist = (unsigned *) finalDB->exist_l;
				*--finalDB->exist = finalDB->prefix_len;
				*--finalDB->exist = finalDB->mlen;
				*--finalDB->exist = finalDB->DB_size;
				msync(finalDB->exist, size, MS_ASYNC);
				munmap(finalDB->exist, size);
				if(finalDB->v_index <= UINT_MAX) {
					munmap(finalDB->exist, finalDB->size * sizeof(unsigned));
				} else {
					munmap(finalDB->exist_l, finalDB->size * sizeof(long unsigned));
				}
			} else {
				/* dump piece of DB */
				cfwrite(&finalDB->DB_size, sizeof(unsigned), 1, out);
				cfwrite(&finalDB->mlen, sizeof(unsigned), 1, out);
				cfwrite(&finalDB->prefix_len, sizeof(unsigned), 1, out);
				cfwrite(&finalDB->prefix, sizeof(long unsigned), 1, out);
				cfwrite(&finalDB->size, sizeof(long unsigned), 1, out);
				cfwrite(&finalDB->n, sizeof(long unsigned), 1, out);
				cfwrite(&finalDB->v_index, sizeof(long unsigned), 1, out);
				cfwrite(&finalDB->null_index, sizeof(long unsigned), 1, out);
				if(finalDB->v_index <= UINT_MAX) {
					cfwrite(finalDB->exist, sizeof(unsigned), finalDB->size, out);
					size += finalDB->size * sizeof(unsigned);
				} else {
					cfwrite(finalDB->exist_l, sizeof(long unsigned), finalDB->size, out);
					size += finalDB->size * sizeof(long unsigned);
				}
				free(finalDB->exist);
			}
			finalDB->exist = 0;
			
			/* set poiters */
			if(finalDB->v_index <= UINT_MAX) {
				getExistPtr = &getExist;
				hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
				size += finalDB->size * sizeof(unsigned);
			} else {
				getExistPtr = &getExistL;
				hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
				size += finalDB->size * sizeof(long unsigned);
			}
		}
	}
	
	/* both allocations failed */
	if(swap & 4) {
		fprintf(stderr, "Fail over to swap.\n");
		v_size = size;
		if(finalDB->DB_size < USHRT_MAX) {
			size += v_index * sizeof(short unsigned);
		} else {
			size += v_index * sizeof(unsigned);
		}
		sfseek(out, size, SEEK_SET);
		fputc(0, out);
		rewind(out);
		finalDB->values = mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(out), 0);
		if(finalDB->values == MAP_FAILED) {
			ERROR();
		}
		posix_madvise(finalDB->values, size, POSIX_MADV_SEQUENTIAL);
		finalDB->values = (unsigned *)(((unsigned char *) finalDB->values) + v_size);
		if(finalDB->DB_size < USHRT_MAX) {
			finalDB->values_s = (short unsigned *) finalDB->values;
		}
	}
	
	/* move values */
	if(finalDB->DB_size < USHRT_MAX) {
		for(node = table; node != 0; node = next) {
			next = node->next;
			values_s = (short unsigned *)(node->values);
			finalV_s = finalDB->values_s + node->v_index - 1;
			i = 2 + *values_s--;
			while(--i) {
				*++finalV_s = *++values_s;
			}
			free(node->values);
			free(node);
			
			/*
			next = node->next;
			values_s = (short unsigned *)(node->values);
			for(i = node->v_index, j = 0; j <= *values_s; ++i, ++j) {
				finalDB->values_s[i] = values_s[j];
			}
			free(values_s);
			free(node);
			*/
		}
	} else {
		for(node = table; node != 0; node = next) {
			next = node->next;
			values = node->values;
			finalV = finalDB->values + node->v_index - 1;
			i = 2 + *values--;
			while(--i) {
				*++finalV = *++values;
			}
			free(node->values);
			free(node);
			
			/*
			next = node->next;
			values = node->values;
			for(i = node->v_index, j = 0; j <= *values; ++i, ++j) {
				finalDB->values[i] = values[j];
			}
			free(values);
			free(node);
			*/
		}
	}
	/* dump final DB */
	fprintf(stderr, "# Dumping compressed DB\n");
	if(swap < 2) {
		megaMapKMA_dump(finalDB, out);
	} else {
		/* dump values */
		if(swap & 4) {
			finalDB->values = (unsigned *)(((unsigned char *) finalDB->values) - v_size);
			finalDB->exist_l = (long unsigned *)(finalDB->values + 3);
			finalDB->exist_l[3] = finalDB->v_index;
			msync(finalDB->values, size, MS_ASYNC);
			munmap(finalDB->values, size);
		}
		
		/* set pointers */
		if(finalDB->DB_size < USHRT_MAX) {
			getValuePtr = &getValueS;
			getSizePtr = &getSizeS;
		} else {
			getValuePtr = &getValue;
			getSizePtr = &getSize;
		}
		
		if(swap & 2) {
			free(finalDB->values);
		}
		free(finalDB);
		finalDB = 0;
	}
	
	return finalDB;
}

void compressKMA_deconDB(HashMapKMA *finalDB, unsigned **Values) {
	
	long unsigned i, j, v_index, new_index, check;
	unsigned *values;
	short unsigned *values_s;
	ValuesHash *shmValues;
	ValuesTable *node, *next, *table;
	
	/* prepare final DB */
	check = 0;
	check = ~check;
	if(finalDB->v_index < UINT_MAX) {
		check >>= 32;
		hashMapKMA_addValue_ptr = &hashMapKMA_addValue;
		getValueIndexPtr = &getValueIndex;
	} else {
		hashMapKMA_addValue_ptr = &hashMapKMA_addValueL;
		getValueIndexPtr = &getValueIndexL;
	}
	i = finalDB->n;
	shmValues = initialize_hashValues(finalDB->n, finalDB->DB_size);
	v_index = 0;
	while(i--) {
		/* the actual value index */
		values = Values[i];
		new_index = valuesHash_add(shmValues, values, v_index);
		
		if(new_index == v_index) {
			v_index += valuesSize(values);
			if(check <= v_index) {
				fprintf(stderr, "# Compression overflow.\n");
				check = 0;
				check = ~check;
				hashMapKMA_addValue_ptr = &hashMapKMA_addValueL;
				getValueIndexPtr = &getValueIndexL;
				finalDB->value_index_l = realloc(finalDB->value_index, finalDB->n * sizeof(long unsigned));
				if(!finalDB->value_index_l) {
					ERROR();
				}
				finalDB->value_index = (unsigned *)(finalDB->value_index_l);
				j = finalDB->n;
				while(j--) {
					finalDB->value_index_l[j] = finalDB->value_index[j];
				}
				finalDB->value_index = 0;
			}
		} else {
			free(values);
		}
		
		hashMapKMA_addValue_ptr(finalDB, i, new_index);
	}
	free(Values);
	/* convert valuesHash to a linked list */
	table = 0;
	i = shmValues->size;
	while(i--) {
		for(node = shmValues->table[i]; node != 0; node = next) {
			next = node->next;
			node->next = table;
			table = node;
		}
	}
	free(shmValues->table);
	
	/* make compressed values */
	fprintf(stderr, "# Finalizing indexes.\n");
	finalDB->v_index = v_index;
	if(finalDB->DB_size < USHRT_MAX) {
		finalDB->values = 0;
		finalDB->values_s = calloc(v_index, sizeof(short unsigned));
		if(!finalDB->values_s) {
			ERROR();
		}
		/* move values */
		for(node = table; node != 0; node = next) {
			next = node->next;
			values_s = (short unsigned *)(node->values);
			for(i = node->v_index, j = 0; j <= *values_s; ++i, ++j) {
				finalDB->values_s[i] = values_s[j];
			}
			free(values_s);
			free(node);
		}
	} else {
		finalDB->values = calloc(v_index, sizeof(unsigned));
		finalDB->values_s = 0;
		if(!finalDB->values) {
			ERROR();
		}
		/* move values */
		for(node = table; node != 0; node = next) {
			next = node->next;
			values = node->values;
			for(i = node->v_index, j = 0; j <= *values; ++i, ++j) {
				finalDB->values[i] = values[j];
			}
			free(values);
			free(node);
		}
	}
}

void compressKMA_deconMegaDB(HashMapKMA *finalDB, unsigned **Values) {
	
	long unsigned i, j, v_index, new_index, pos, check;
	unsigned *values;
	short unsigned *values_s;
	ValuesHash *shmValues;
	ValuesTable *node, *next, *table;
	
	fprintf(stderr, "# Compressing indexes.\n");
	check = 0;
	check = ~check;
	if(finalDB->v_index < UINT_MAX) {
		check >>= 32;
		getExistPtr = &getExist;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
	} else {
		getExistPtr = &getExistL;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
	}
	i = finalDB->size;
	shmValues = initialize_hashValues(finalDB->n, finalDB->DB_size);
	v_index = 0;
	while(i--) {
		if((pos = getExistPtr(finalDB->exist, i)) != finalDB->n) {
			values = Values[pos];
			new_index = valuesHash_add(shmValues, values, v_index);
			
			if(new_index == v_index) {
				v_index += valuesSize(values);
				if(check <= v_index) {
					fprintf(stderr, "# Compression overflow.\n");
					check = 0;
					check = ~check;
					getExistPtr = &getExistL;
					hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
					finalDB->exist_l = realloc(finalDB->exist, finalDB->size * sizeof(long unsigned));
					if(!finalDB->value_index_l) {
						ERROR();
					}
					finalDB->exist = (unsigned *)(finalDB->exist_l);
					j = finalDB->size;
					while(j--) {
						finalDB->exist_l[j] = finalDB->exist[j];
					}
					finalDB->exist = 0;
				}
			} else {
				free(values);
			}
			hashMapKMA_addExist_ptr(finalDB, i, new_index);
		} else {
			hashMapKMA_addExist_ptr(finalDB, i, 1);
		}
	}
	free(Values);
	/* convert valuesHash to a linked list */
	table = 0;
	i = shmValues->size;
	while(i--) {
		for(node = shmValues->table[i]; node != 0; node = next) {
			next = node->next;
			node->next = table;
			table = node;
		}
	}
	free(shmValues->table);
	
	/* make compressed values */
	fprintf(stderr, "# Finalizing indexes.\n");
	finalDB->v_index = v_index;
	finalDB->null_index = 1;
	if(finalDB->DB_size < USHRT_MAX) {
		finalDB->values = 0;
		finalDB->values_s = calloc(v_index, sizeof(short unsigned));
		if(!finalDB->values_s) {
			ERROR();
		}
		/* move values */
		for(node = table; node != 0; node = next) {
			next = node->next;
			values_s = (short unsigned *)(node->values);
			for(i = node->v_index, j = 0; j <= *values_s; ++i, ++j) {
				finalDB->values_s[i] = values_s[j];
			}
			free(values_s);
			free(node);
		}
	} else {
		finalDB->values = calloc(v_index, sizeof(unsigned));
		finalDB->values_s = 0;
		if(!finalDB->values) {
			ERROR();
		}
		/* move values */
		for(node = table; node != 0; node = next) {
			next = node->next;
			values = node->values;
			for(i = node->v_index, j = 0; j <= *values; ++i, ++j) {
				finalDB->values[i] = values[j];
			}
			free(values);
			free(node);
		}
	}
}
