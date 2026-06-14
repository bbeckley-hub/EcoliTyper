/* Philip T.L.C. Clausen Mar 2025 plan@dtu.dk */

/*
 * Copyright (c) 2025, Philip Clausen, Technical University of Denmark
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
#include "merge.h"
#include "hashmapkma.h"
#include "kmmap.h"
#include "middlelayer.h"
#include "pherror.h"
#include "stdstat.h"
#include "tmp.h"
#include "version.h"
#ifdef _WIN32
#define mmap(addr, len, prot, flags, fd, offset) (0; fprintf(stderr, "mmap not available on windows.\n"); exit(1););
#define munmap(addr, len) (-1);
#else
#include <sys/mman.h>
#endif

long unsigned bucket_insertionsort(unsigned *key_index, unsigned *value_index, long unsigned *value_index_l, long unsigned mask_new, long unsigned mask_org, unsigned flag) {
	
	unsigned *start, *end;
	long unsigned n, min, mink, kmer, ks, ke, pos, index;
	
	/* find end of bucket */
	n = 0;
	pos = 0;
	start = key_index;
	end = start;
	mink = *start;
	ks = mink;
	if(flag) {
		murmur(ks, ks);
	}
	min = ks & mask_new;
	ks &= mask_org;
	ke = ks;
	while(ks == ke) {
		kmer = *++end;
		ke = kmer;
		if(flag) {
			murmur(ke, ke);
		}
		index = ke & mask_new;
		ke &= mask_org;
		if((ks != ke) || (min <= index && (min < index || mink < kmer))) { /* check if sorted */
			if(ks == ke) {
				min = index;
				mink = kmer;
			}
			++pos;
		}
		++n;
	}
	
	/* already sorted */
	if(pos == n) {
		return n;
	}
	
	/* run insertion sort */
	while(start < end) {
		/* find min */
		kmer = *start;
		ks = kmer;
		if(flag) {
			murmur(ks, ks);
		}
		min = ks & mask_new;
		mink = kmer;
		pos = 0;
		index = 1;
		while(++start < end) {
			kmer = *start;
			ks = kmer;
			if(flag) {
				murmur(ks, ks);
			}
			ks &= mask_new;
			if(ks <= min && (ks < min || kmer < mink)) { /* secondary sort on kmer */
				min = ks;
				mink = kmer;
				pos = index;
			}
			++index;
		}
		
		/* swap */
		if(value_index) {
			ks = *value_index;
			*value_index = value_index[pos];
			value_index[pos] = ks;
			++value_index;
		} else {
			ks = *value_index_l;
			*value_index_l = value_index_l[pos];
			value_index_l[pos] = ks;
			++value_index_l;
		}
		ks = *key_index;
		*key_index = key_index[pos];
		key_index[pos] = ks;
		start = ++key_index;
	}
	
	return n;
}

long unsigned bucket_insertionsort_l(long unsigned *key_index, unsigned *value_index, long unsigned *value_index_l, long unsigned mask_new, long unsigned mask_org, unsigned flag) {
	
	long unsigned n, min, mink, kmer, ks, ke, pos, index, *start, *end;
	
	/* find end of bucket */
	n = 0;
	pos = 0;
	start = key_index;
	end = start;
	mink = *start;
	ks = mink;
	if(flag) {
		murmur(ks, ks);
	}
	min = ks & mask_new;
	ks &= mask_org;
	ke = ks;
	while(ks == ke) {
		kmer = *++end;
		ke = kmer;
		if(flag) {
			murmur(ke, ke);
		}
		index = ke & mask_new;
		ke &= mask_org;
		if((ks != ke) || (min <= index && (min < index || mink < kmer))) { /* check if sorted */
			if(ks == ke) {
				min = index;
				mink = kmer;
			}
			++pos;
		}
		++n;
	}
	
	/* already sorted */
	if(pos == n) {
		return n;
	}
	
	/* run insertion sort */
	while(start < end) {
		/* find min */
		kmer = *start;
		ks = kmer;
		if(flag) {
			murmur(ks, ks);
		}
		min = ks & mask_new;
		mink = kmer;
		pos = 0;
		index = 1;
		while(++start < end) {
			kmer = *start;
			ks = kmer;
			if(flag) {
				murmur(ks, ks);
			}
			ks &= mask_new;
			if(ks <= min && (ks < min || kmer < mink)) { /* secondary sort on kmer */
				min = ks;
				mink = kmer;
				pos = index;
			}
			++index;
		}
		
		/* swap */
		if(value_index) {
			ks = *value_index;
			*value_index = value_index[pos];
			value_index[pos] = ks;
			++value_index;
		} else {
			ks = *value_index_l;
			*value_index_l = value_index_l[pos];
			value_index_l[pos] = ks;
			++value_index_l;
		}
		ks = *key_index;
		*key_index = key_index[pos];
		key_index[pos] = ks;
		start = ++key_index;
	}
	
	return n;
}

void hashMapKMA_sortbuckets(HashMapKMA *dest, HashMapKMA *src) {
	
	unsigned flag, *exist, *key_index, *value_index;
	long unsigned n, null_index, size, mask, mask_src, index, kmer;
	long unsigned *exist_l, *key_index_l, *value_index_l;
	
	/* load value_indexes */
	if(src->v_index < UINT_MAX) {
		memcpy(dest->value_index, src->value_index, src->n * sizeof(unsigned));
		value_index = dest->value_index;
		value_index_l = 0;
	} else {
		memcpy(dest->value_index_l, src->value_index_l, src->n * sizeof(long unsigned));
		value_index_l = dest->value_index_l;
		value_index = 0;
	}
	if(dest->mlen <= 16) {
		key_index = dest->key_index;
		key_index_l = 0;
	} else {
		key_index = 0;
		key_index_l = dest->key_index_l;
	}
	
	/* sort buckets */
	flag = dest->flag;
	mask = dest->size;
	mask_src = src->size;
	size = src->n;
	while(size) {
		if(key_index) {
			n = bucket_insertionsort(key_index, value_index, value_index_l, mask, mask_src, flag);
			key_index += n;
		} else {
			n = bucket_insertionsort_l(key_index_l, value_index, value_index_l, mask, mask_src, flag);
			key_index_l += n;
		}
		if(value_index) {
			value_index += n;
		} else {
			value_index_l += n;
		}
		size -= n;
	}
	
	/* reset exist */
	if(dest->exist) {
		exist = dest->exist - 1;
		exist_l = 0;
	} else {
		exist = 0;
		exist_l = dest->exist_l - 1;
	}
	null_index = src->n;
	dest->n = null_index;
	size = dest->size + 2;
	if(exist) {
		while(--size) {
			*++exist = null_index;
		}
		exist = dest->exist;
	} else {
		while(--size) {
			*++exist_l = null_index;
		}
		exist_l = dest->exist_l;
	}
	
	/* add links to exist */
	if(dest->mlen <= 16) {
		key_index = dest-> key_index - 1;
		key_index_l = 0;
	} else {
		key_index = 0;
		key_index_l = dest->key_index_l - 1;
	}
	index = 0;
	size = null_index;
	while(size) {
		kmer = key_index ? *++key_index : *++key_index_l;
		if(flag) {
			murmur(kmer, kmer);
		}
		kmer &= mask;
		if(exist && exist[kmer] == null_index) {
			exist[kmer] = index;
		} else if(exist_l && exist_l[kmer] == null_index) {
			exist_l[kmer] = index;
		}
		++index;
		--size;
	}
	
}

FILE * hashMapKMA_dumpbuckets(HashMapKMA *src) {
	
	unsigned flag, *exist, *key_index, *value_index;
	long unsigned n, index, null_index, mask, kmer, key;
	long unsigned buff[2], *exist_l, *key_index_l, *value_index_l;
	FILE *tmp;
	
	/* init */
	tmp = tmpF(0);
	flag = src->flag;
	null_index = src->null_index;
	mask = src->size;
	if(src->exist) {
		exist = src->exist - 1;
		exist_l = 0;
	} else {
		exist = 0;
		exist_l = src->exist_l - 1;
	}
	if(src->mlen <= 16) {
		key_index = src->key_index;
		key_index_l = 0;
	} else {
		key_index = 0;
		key_index_l = src->key_index_l;
	}
	if(src->v_index < UINT_MAX) {
		value_index = src->value_index;
		value_index_l = 0;
	} else {
		value_index = 0;
		value_index_l = src->value_index_l;
	}
	
	/* get k-mers ordered according to dest */
	n = src->n;
	key = 0;
	while(n) {
		index = exist ? *++exist : *++exist_l;
		if(index != null_index) {
			buff[0] = key_index ? key_index[index] : key_index_l[index];
			do {
				/* get pair */
				buff[1] = value_index ? value_index[index] : value_index_l[index];
				
				/* dump key-value */
				fwrite(buff, sizeof(long unsigned), 2, tmp);
				--n;
				
				/* get next kmer in bucket */
				++index;
				kmer = key_index ? key_index[index] : key_index_l[index];
				buff[0] = kmer;
				if(flag) {
					murmur(kmer, kmer);
				}
			} while(key == (kmer & mask));
		}
		++key;
	}
	
	/* flush and rewind */
	fflush(tmp);
	rewind(tmp);
	
	return tmp;
}

long unsigned getV_index(long unsigned *exist, long unsigned size, long unsigned null_index, long unsigned v_indexes, MiddleLayer *middle, MiddleLayer *alternative) {
	
	long unsigned n, index, v_index, *layer, *altlayer;
	
	/* get v_index */
	n = 0;
	v_index = 0;
	--exist;
	++size;
	while(--size) {
		index = *++exist;
		if(index != null_index) {
			if(index < middle->n) {
				layer = middle->layer + (index << 1);
				if(layer[1] != null_index) {
					/* get signature */
					if(index != (middle->n - 1)) {
						v_index += layer[2] - *layer;
					} else {
						v_index += v_indexes - *layer;
					}
					/* mark as visited */
					layer[1] = null_index;
				}
			} else {
				index -= middle->n;
				altlayer = alternative->layer + ((index << 1) + index);
				if(altlayer[2] != null_index) {
					/* get first signature */
					layer = middle->layer + (*altlayer << 1);
					v_index += layer[2] - *layer - 1;
					*altlayer = *layer;
					
					/* get second signature */
					layer = middle->layer + (altlayer[1] << 1);
					if(altlayer[1] != (middle->n - 1)) {
						v_index += layer[2] - *layer;
					} else {
						v_index += v_indexes - *layer;
					}
					altlayer[1] = *layer;
					
					/* mark as visited */
					altlayer[2] = null_index;
				}
			}
			++n;
		}
	}
	
	return v_index;
}

unsigned * adjustV_index(long unsigned *exist_l, long unsigned size, long unsigned v_index) {
	
	unsigned *exist;
	long unsigned n, *eptr;
	
	exist = (unsigned *)(exist_l);
	if(v_index <= UINT_MAX) {
		--exist;
		eptr = exist_l - 1;
		n = size + 1;
		while(--n) {
			*++exist = *++eptr;
		}
		exist = realloc(exist_l, size * sizeof(unsigned));
		if(!exist) {
			exist = (unsigned *)(exist_l);
		} else {
			exist_l = (long unsigned *)(exist);
		}
		/*
		worst case new indexes:
		n1 + 2 * n2
		
		worst case values is one on each org, so new will be:
		v_index1 + v_index2 - n1
		
		Thus following is satified (guarantying no overflows):
		n1 + 2 * n2 <= v_index1 + v_index2 - n1 <->
		2*n1 <= v_index1 && 2*n2 <= v_index2
		*/
	}
	
	return exist;
}

long unsigned add_pairs(long unsigned *exist, long unsigned size, long unsigned null_index, long unsigned v_indexes, MiddleLayer *middle, MiddleLayer *alternative, HashMapKMA *t2) {
	
	unsigned *exist_t, *value_index;
	long unsigned n, index, v_index;
	long unsigned *layer, *exist_lt, *value_index_l;
	
	if(t2->kmersize <= 16) {
		exist_t = t2->key_index - 1;
		exist_lt = 0;
	} else {
		exist_t = 0;
		exist_lt = t2->key_index_l - 1;
	}
	if(t2->v_index < UINT_MAX) {
		value_index = t2->value_index - 1;
		value_index_l = 0;
	} else {
		value_index = 0;
		value_index_l = t2->value_index_l - 1;
	}
	n = 0;
	++size;
	while(--size) {
		index = exist_t ? *++exist_t : *++exist_lt;
		v_index = value_index ? *++value_index : *++value_index_l;
		if(exist[index] == null_index) { /* new k-mer */
			exist[index] = MiddleLayer_search(middle, v_indexes + v_index);
			++n;
		} else { /* get combo of signatures */ 
			/* search (new) combo */
			layer = middle->layer + (exist[index] << 1);
			layer[1] = layer[1] ? layer[1] : (middle->n + alternative->n);
			exist[index] = middle->n + AlternativeLayer_add(alternative, layer[1] - middle->n, exist[index], MiddleLayer_search(middle, v_indexes + v_index));
		}
	}
	
	return n;
}

void hashMapKMA_merge(HashMapKMA *dest, MiddleLayer *middle, MiddleLayer *alternative, FILE *tmp_1, FILE *tmp_2, long unsigned v_indexes) {
	
	unsigned flag, *exist, *key_index;
	long unsigned n, size, v_index, null_index, index, value;
	long unsigned kmer, kmer1, kmer2, buff1[2], buff2[2];
	long unsigned *layer, *exist_l, *key_index_l, *value_index_l;
	
	/* init */
	n = 0;
	null_index = dest->null_index;
	v_index = dest->v_index;
	size = dest->size;
	flag = dest->flag;
	if(dest->n <= UINT_MAX) {
		exist = dest->exist - 1;
		exist_l = 0;
	} else {
		exist = 0;
		exist_l = dest->exist_l - 1;
	}
	if(dest->mlen <= 16) {
		key_index = dest->key_index - 1;
		key_index_l = 0;
	} else {
		key_index = 0;
		key_index_l = dest->key_index_l - 1;
	}
	value_index_l = dest->value_index_l - 1;
	
	index = 0;
	if(fread(buff1, sizeof(long unsigned), 2, tmp_1)) {
		if(flag) {
			murmur(kmer1, buff1[0]);
			kmer1 &= size;
		} else {
			kmer1 = buff1[0] & size;
		}
	} else {
		/* no more k-mers */
		buff1[0] = index - 1;
		kmer1 = index - 1;
	}
	if(fread(buff2, sizeof(long unsigned), 2, tmp_2)) {
		if(flag) {
			murmur(kmer2, buff2[0]);
			kmer2 &= size;
		} else {
			kmer2 = buff2[0] & size;
		}
	} else {
		/* no more k-mers */
		buff2[0] = index - 1;
		kmer2 = index - 1;
	}
	kmer = 0;
	
	/* iterate pairs */
	index = 0;
	while(index <= size) {
		if(index == kmer1 || index == kmer2) { /* occupied */
			/* update hashmap */
			if(exist) {
				*++exist = n;
			} else {
				*++exist_l = n;
			}
			
			/* add key-value pairs of bucket */
			while(index == kmer1 || index == kmer2) {
				if(buff1[0] == buff2[0]) { /* new combo */
					/* get position in middle layer */
					value = MiddleLayer_search(middle, buff1[1]);
					layer = middle->layer + (value << 1);
					/* get combo */
					kmer = buff1[0];
					layer[1] = layer[1] ? layer[1] : (middle->n + alternative->n);
					value = middle->n + AlternativeLayer_add(alternative, layer[1] - middle->n, value, MiddleLayer_search(middle, v_index + buff2[1]));
					
					/* load next pairs */
					if(fread(buff1, sizeof(long unsigned), 2, tmp_1)) {
						if(flag) {
							murmur(kmer1, buff1[0]);
							kmer1 &= size;
						} else {
							kmer1 = buff1[0] & size;
						}
					} else {
						/* no more k-mers */
						buff1[0] = kmer1 - 1;
						kmer1 = index - 1;
					}
					if(fread(buff2, sizeof(long unsigned), 2, tmp_2)) {
						if(flag) {
							murmur(kmer2, buff2[0]);
							kmer2 &= size;
						} else {
							kmer2 = buff2[0] & size;
						}
					} else {
						/* no more k-mers */
						buff2[0] = kmer2 - 1;
						kmer2 = index - 1;
					}
				} else if(kmer1 == kmer2) { /* sync k-mers */
					if(buff1[0] < buff2[0]) { /* add k-mer from t1 */
						kmer = buff1[0];
						value = MiddleLayer_search(middle, buff1[1]);
						/* load next pair */
						if(fread(buff1, sizeof(long unsigned), 2, tmp_1)) {
							if(flag) {
								murmur(kmer1, buff1[0]);
								kmer1 &= size;
							} else {
								kmer1 = buff1[0] & size;
							}
						} else {
							/* no more k-mers */
							buff1[0] = kmer1 - 1;
							kmer1 = index - 1;
						}
					} else { /* add k-mer from t2 */
						kmer = buff2[0];
						value = MiddleLayer_search(middle, v_index + buff2[1]);
						/* load next pair */
						if(fread(buff2, sizeof(long unsigned), 2, tmp_2)) {
							if(flag) {
								murmur(kmer2, buff2[0]);
								kmer2 &= size;
							} else {
								kmer2 = buff2[0] & size;
							}
						} else {
							/* no more k-mers */
							buff2[0] = kmer2 - 1;
							kmer2 = index - 1;
						}
					}
				} else if(kmer1 == index) { /* add k-mer from t1 */
					kmer = buff1[0];
					value = MiddleLayer_search(middle, buff1[1]);
					/* load next pair */
					if(fread(buff1, sizeof(long unsigned), 2, tmp_1)) {
						if(flag) {
							murmur(kmer1, buff1[0]);
							kmer1 &= size;
						} else {
							kmer1 = buff1[0] & size;
						}
					} else {
						/* no more k-mers */
						buff1[0] = kmer1 - 1;
						kmer1 = index - 1;
					}
				} else { /* add k-mer from t2 */
					kmer = buff2[0];
					value = MiddleLayer_search(middle, v_index + buff2[1]);
					/* load next pair */
					if(fread(buff2, sizeof(long unsigned), 2, tmp_2)) {
						if(flag) {
							murmur(kmer2, buff2[0]);
							kmer2 &= size;
						} else {
							kmer2 = buff2[0] & size;
						}
					} else {
						/* no more k-mers */
						buff2[0] = kmer2 - 1;
						kmer2 = index - 1;
					}
				}
				
				/* update key-value indexes */
				if(key_index) {
					*++key_index = kmer;
				} else {
					*++key_index_l = kmer;
				}
				*++value_index_l = value;
				++n;
			}
		} else {
			if(exist) {
				*++exist = dest->n;
			} else {
				*++exist_l = dest->n;
			}
		}
		++index;
	}
	
	/* terminate pairs */
	if(flag) {
		murmur(kmer1, kmer);
		kmer1 &= size;
		++kmer;
		murmur(kmer2, kmer);
		kmer2 &= size;
		while(kmer1 == kmer2) {
			++kmer;
			murmur(kmer2, kmer);
			kmer2 &= size;
		}
	} else {
		kmer1 = kmer & size;
		kmer2 = ++kmer & size;
		while(kmer1 == kmer2) {
			kmer2 = ++kmer & size;
		}
	}
	if(key_index) {
		*++key_index = kmer;
	} else {
		*++key_index_l = kmer;
	}
	
	/* test if everything was read */
	if(dest->n != n) {
		fprintf(stderr, "Did not get the expected pairs., %lu, %lu\n", dest->n, n);
		exit(1);
	}
	
	/* adjust new alternative layer */
	AlternativeLayer_readjust(alternative);
	
	/* get v_index */
	fprintf(stderr, "# Getting new v_index.\n");
	v_index = getV_index(dest->value_index_l, n, null_index, v_indexes, middle, alternative);
	dest->v_index = v_index;
	
	/* adjust size of value_index */
	dest->value_index = adjustV_index(dest->value_index_l, dest->n, v_index);
	dest->value_index_l = (long unsigned *)(dest->value_index);
}

unsigned loadValues1(unsigned *values, short unsigned *values_s, unsigned *values1, short unsigned *values1_s, long unsigned index) {
	
	unsigned n, size;
	
	/* get signature */
	if(values1) {
		values1 += index;
		n = *values1;
	} else {
		values1_s += index;
		n = *values1_s;
	}
	
	/* cp signature */
	size = n + 1;
	if(values1) {
		*values = n;
		while(--size) {
			*++values = *++values1;
		}
	} else if(values) {
		*values = n;
		while(--size) {
			*++values = *++values1_s;
		}
	} else {
		*values_s = n;
		while(--size) {
			*++values_s = *++values1_s;
		}
	}
	
	return n + 1;
}

unsigned loadValues2(unsigned *values, short unsigned *values_s, unsigned *values2, short unsigned *values2_s, long unsigned index, unsigned offset) {
	
	unsigned n, size;
	
	/* get signature */
	if(values2) {
		values2 += index;
		n = *values2;
	} else {
		values2_s += index;
		n = *values2_s;
	}
	
	/* cp signature */
	size = n + 1;
	if(values2) {
		*values = n;
		while(--size) {
			*++values = *++values2 + offset;
		}
	} else if(values) {
		*values = n;
		while(--size) {
			*++values = *++values2_s + offset;
		}
	} else {
		*values_s = n;
		while(--size) {
			*++values_s = *++values2_s + offset;
		}
	}
	
	return n + 1;
}

unsigned loadValues12(unsigned *values, short unsigned *values_s, unsigned *values1, short unsigned *values1_s, unsigned *values2, short unsigned *values2_s, long unsigned index1, long unsigned index2, unsigned offset) {
	
	unsigned n, size, *values_ptr;
	short unsigned *values_sptr;
	
	/* get signature from t1 */
	if(values1) {
		values1 += index1;
		n = *values1;
	} else {
		values1_s += index1;
		n = *values1_s;
	}
	
	/* cp signature from t1 */
	size = n + 1;
	values_ptr = values;
	values_sptr = values_s;
	if(values1) {
		*values = n;
		while(--size) {
			*++values_ptr = *++values1;
		}
	} else if(values) {
		*values = n;
		while(--size) {
			*++values_ptr = *++values1_s;
		}
	} else {
		*values_s = n;
		while(--size) {
			*++values_sptr = *++values1_s;
		}
	}
	
	/* get signature from t2 */
	if(values2) {
		values2 += index2;
		n = *values2;
	} else {
		values2_s += index2;
		n = *values2_s;
	}
	
	/* append signature from t2 */
	size = n + 1;
	if(values2) {
		n += *values;
		*values = n;
		while(--size) {
			*++values_ptr = *++values2 + offset;
		}
	} else if(values) {
		n += *values;
		*values = n;
		while(--size) {
			*++values_ptr = *++values2_s + offset;
		}
	} else {
		n += *values_s;
		*values_s = n;
		while(--size) {
			*++values_sptr = *++values2_s + offset;
		}
	}
	
	return n + 1;
}

void hashMapKMA_dumpmerge(HashMapKMA *src, HashMapKMA *t1, HashMapKMA *t2, MiddleLayer *middle, MiddleLayer *alternative, FILE *out) {
	
	unsigned offset, *values, *values1, *values2, *exist;
	long unsigned size, n, mn, index, null_index, v_index, v_index1;
	long unsigned *layer, *exist_l;
	short unsigned *values_s, *values1_s, *values2_s;
	
	/* init */
	offset = t1->DB_size - 1;
	n = 0;
	mn = middle->n;
	v_index1 = t1->v_index;
	null_index = src->null_index;
	src->null_index = (src->size != (src->mask + 1)) ? src->n : 1;
	
	/* dump sizes */
	cfwrite(&src->DB_size, sizeof(unsigned), 1, out);
	cfwrite(&src->mlen, sizeof(unsigned), 1, out);
	cfwrite(&src->prefix_len, sizeof(unsigned), 1, out);
	cfwrite(&src->prefix, sizeof(long unsigned), 1, out);
	cfwrite(&src->size, sizeof(long unsigned), 1, out);
	cfwrite(&src->n, sizeof(long unsigned), 1, out);
	cfwrite(&src->v_index, sizeof(long unsigned), 1, out);
	cfwrite(&src->null_index, sizeof(long unsigned), 1, out);
	
	if(src->size != (src->mask + 1)) { /* hashmap */
		size = (src->n <= UINT_MAX) ? sizeof(unsigned) : sizeof(long unsigned);
		cfwrite(src->exist, size, src->size, out);
		free(src->exist);
		src->exist = 0;
		src->exist_l = 0;
	} else { /* direct */
		/* make space for exist for when it is done */
		size = (src->v_index <= UINT_MAX) ? sizeof(unsigned) : sizeof(long unsigned);
		sfseek(out, size * src->size, SEEK_CUR);
	}
	
	/* set pointers to values */
	if(src->DB_size < USHRT_MAX) {
		values_s = smalloc((src->DB_size + 1) * sizeof(short unsigned));
		values = 0;
	} else {
		values_s = 0;
		values = smalloc((src->DB_size + 1) * sizeof(unsigned));
	
	}
	if(t1->DB_size < USHRT_MAX) {
		values1_s = t1->values_s;
		values1 = 0;
	} else {
		values1_s = 0;
		values1 = t1->values;
	}
	if(t2->DB_size < USHRT_MAX) {
		values2_s = t2->values_s;
		values2 = 0;
	} else {
		values2_s = 0;
		values2 = t2->values;
	}
	
	/* dump new values */
	fprintf(stderr, "# Creating new signatures.\n");
	v_index = 0;
	if(src->exist) { /* direct */
		size = src->size + 1;
		src->null_index = 1;
		if(src->v_index <= UINT_MAX) {
			exist = src->exist - 1;
			exist_l = 0;
		} else {
			exist = 0;
			exist_l = src->exist_l - 1;
		}
	} else { /* hashmap */
		size = src->n + 1;
		if(src->v_index < UINT_MAX) {
			exist = src->value_index - 1;
			exist_l = 0;
		} else {
			exist = 0;
			exist_l = src->value_index_l - 1;
		}
	}
	while(--size) {
		index = exist ? *++exist : *++exist_l;
		if(index != null_index) {
			/* get signature */
			if(index < mn) { /* old signature */
				layer = middle->layer + (index << 1);
				if(layer[1]) { /* not dumped */
					/* get index */
					index = v_index;
					
					/* load signature */
					if(*layer < v_index1) {
						v_index += loadValues1(values, values_s, values1, values1_s, *layer);
					} else {
						v_index += loadValues2(values, values_s, values2, values2_s, *layer - v_index1, offset);
					}
					
					/* dump signature */
					if(values) {
						cfwrite(values, sizeof(unsigned), *values + 1, out);
					} else {
						cfwrite(values_s, sizeof(short unsigned), *values_s + 1, out);
					}
					
					/* mark layer as dumped */
					*layer = index;
					layer[1] = 0;
				} else {
					/* get index */
					index = *layer;
				}
			} else { /* new signature */
				index -= mn;
				layer = alternative->layer + ((index << 1) + index);
				if(layer[2]) { /* not dumped */
					/* get index */
					index = v_index;
					
					/* load signatures */
					v_index += loadValues12(values, values_s, values1, values1_s, values2, values2_s, *layer, layer[1] - v_index1, offset);
					
					/* dump signature */
					if(values) {
						cfwrite(values, sizeof(unsigned), *values + 1, out);
					} else {
						cfwrite(values_s, sizeof(short unsigned), *values_s + 1, out);
					}
					
					/* mark layer as dumped */
					*layer = index;
					layer[1] = 0;
					layer[2] = 0;
				} else {
					/* get index */
					index = *layer;
				}
			}
			if(exist) {
				*exist = index;
			} else {
				*exist_l = index;
			}
			++n;
		} else { /* not entered if hashmap */
			if(exist) {
				*exist = 1;
			} else {
				*exist_l = 1;
			}
		}
	}
	
	if(v_index != src->v_index) {
		fprintf(stderr, "New signatures does not match expected v_index\n");
		
	}
	
	/* clean */
	if(values) {
		free(values);
	} else {
		free(values_s);
	}
	
	/* dump remaining indexes */
	if(!src->exist) {
		if(src->mlen <= 16) {
			cfwrite(src->key_index, sizeof(unsigned), src->n + 1, out);
			free(src->key_index);
		} else {
			cfwrite(src->key_index_l, sizeof(long unsigned), src->n + 1, out);
			free(src->key_index_l);
		}
		if(src->v_index < UINT_MAX) {
			cfwrite(src->value_index, sizeof(unsigned), src->n, out);
			free(src->value_index);
		} else {
			cfwrite(src->value_index_l, sizeof(long unsigned), src->n, out);
			free(src->value_index_l);
		}
	}
	cfwrite(&src->kmersize, sizeof(unsigned), 1, out);
	cfwrite(&src->flag, sizeof(unsigned), 1, out);
	
	/* dump exist */
	if(src->exist) {
		sfseek(out, 3 * sizeof(unsigned) + 5 * sizeof(long unsigned), SEEK_SET);
		if(src->v_index <= UINT_MAX) {
			cfwrite(src->exist, sizeof(unsigned), src->size, out);
			free(src->exist);
		} else {
			cfwrite(src->exist_l, sizeof(long unsigned), src->size, out);
			free(src->exist_l);
		}
		sfseek(out, 0, SEEK_END);
	}
	fflush(out);
	free(src);
}

HashMapKMA * merge_kmersignatures(HashMapKMA *t1, HashMapKMA *t2, MiddleLayer *middle, MiddleLayer *alternative) {
	
	/* output hash with exist, key_index, value_index */
	unsigned flag, *exist, *exist_t, *value_index, *key_index, *key_index_t;
	long unsigned index, null_index, null_index_t, size, n, v_index, mask;
	long unsigned key, kmer, kpos;
	long unsigned *layer, *exist_l, *exist_lt, *value_index_l;
	long unsigned *key_index_l, *key_index_tl;
	FILE *tmp_1, *tmp_2;
	HashMapKMA *dest;
	
	/* init */
	dest = smalloc(sizeof(HashMapKMA));
	dest->size = t1->size;
	dest->n = t1->n;
	dest->mask = t1->mask;
	null_index = middle->n + t2->n; /* unreachable index in middle-/alternative-layer */
	dest->null_index = null_index;
	dest->v_index = t1->v_index;
	dest->kmersize = t1->kmersize;
	dest->prefix_len = t1->prefix_len;
	dest->prefix = t1->prefix;
	dest->DB_size = t1->DB_size + t2->DB_size - 1;
	dest->shmFlag = t1->shmFlag;
	dest->mlen = t1->mlen;
	dest->flag = t1->flag;
	dest->exist = 0;
	dest->exist_l = 0;
	dest->values = 0;
	dest->values_s = 0;
	dest->key_index = 0;
	dest->key_index_l = 0;
	dest->value_index = 0;
	dest->value_index_l = 0;
	flag = dest->flag;
	
	/* direct on t1 */
	if(t1->size == (t1->mask + 1)) {
		exist_l = smalloc(t1->size * sizeof(long unsigned));
		dest->exist_l = exist_l--;
		
		/* populate new exist with t1 */
		fprintf(stderr, "# Getting middlelayer signatures from first index.\n");
		if(t1->v_index < UINT_MAX) {
			exist_t = t1->exist - 1;
			exist_lt = 0;
		} else {
			exist_t = 0;
			exist_lt = t1->exist_l - 1;
		}
		size = dest->size + 1;
		null_index_t = t1->null_index;
		while(--size) {
			index = exist_t ? *++exist_t : *++exist_lt;
			if(index == null_index_t) {
				/* empty */
				*++exist_l = null_index;
			} else {
				/* get index to middlelayer */
				*++exist_l = MiddleLayer_search(middle, index);
			}
		}
		
		/* merge with t2 */
		fprintf(stderr, "# Merging middlelayer signatures with second index.\n");
		if(t2->size == (t2->mask + 1)) { /* direct on t2 */
			if(t2->v_index <= UINT_MAX) {
				exist_t = t2->exist - 1;
				exist_lt = 0;
			} else {
				exist_t = 0;
				exist_lt = t2->exist_l - 1;
			}
			size = dest->size + 1;
			null_index_t = t2->null_index;
			exist_l = dest->exist_l;
			while(--size) {
				index = exist_t ? *++exist_t : *++exist_lt;
				if(index != null_index_t) {
					if(*exist_l == null_index) { /* new k-mer */
						*exist_l = MiddleLayer_search(middle, t1->v_index + index);
						dest->n++;
					} else { /* get combo of signatures */
						/* search (new) combo */
						layer = middle->layer + (*exist_l << 1);
						layer[1] = layer[1] ? layer[1] : (middle->n + alternative->n);
						*exist_l = middle->n + AlternativeLayer_add(alternative, layer[1] - middle->n, *exist_l, MiddleLayer_search(middle, t1->v_index + index));
					}
				}
				++exist_l;
			}
		} else { /* hashmap on t2*/
			dest->n += add_pairs(dest->exist_l, t2->n, null_index, t1->v_index, middle, alternative, t2);
		}
		
		/* adjust new alternative layer */
		AlternativeLayer_readjust(alternative);
		
		/* get v_index */
		fprintf(stderr, "# Getting new v_index.\n");
		v_index = getV_index(dest->exist_l, dest->size, null_index, t1->v_index + t2->v_index, middle, alternative);
		dest->v_index = v_index;
		
		/* adjust size of exist */
		dest->exist = adjustV_index(dest->exist_l, dest->size, v_index);
		dest->exist_l = (long unsigned *)(dest->exist_l);
	} else { /* hashmap on both */
		/* load exist and k-mers */
		if(t1->n <= UINT_MAX) {
			dest->exist = smalloc((t1->size + 1) * sizeof(unsigned));
			memcpy(dest->exist, t1->exist, (t1->size + 1) * sizeof(unsigned));
			exist = dest->exist;
			exist_l = 0;
		} else {
			dest->exist_l = smalloc((t1->size + 1) * sizeof(long unsigned));
			memcpy(dest->exist_l, t1->exist_l, (t1->size + 1) * sizeof(long unsigned));
			exist = 0;
			exist_l = dest->exist_l;
		}
		if(t1->mlen <= 16) {
			dest->key_index = smalloc((t1->n + 1) * sizeof(unsigned));
			memcpy(dest->key_index, t1->key_index, (t1->n + 1) * sizeof(unsigned));
			key_index = dest->key_index;
			key_index_l = 0;
			key_index_t = t2->key_index - 1;
			key_index_tl = 0;
		} else {
			dest->key_index_l = smalloc((t1->n + 1) * sizeof(long unsigned));
			memcpy(dest->key_index_l, t1->key_index_l, (t1->n + 1) * sizeof(long unsigned));
			key_index = 0;
			key_index_l = dest->key_index_l;
			key_index_t = 0;
			key_index_tl = t2->key_index_l - 1;
		}
		
		/* get new n */
		fprintf(stderr, "# Identifying unique k-mers.\n");
		size = t2->n + 1;
		mask = dest->size;
		n = dest->n;
		while(--size) {
			key = key_index_t ? *++key_index_t : *++key_index_tl;
			if(flag) {
				murmur(kpos, key);
				kpos &= mask;
			} else {
				kpos = key & mask;
			}
			
			index = exist ? exist[kpos] : exist_l[kpos];
			if(index != t1->null_index) {
				kmer = key_index ? key_index[index] : key_index_l[index];
				while(key != kmer) {
					if(flag) {
						murmur(kmer, kmer);
					}
					if(kpos != (kmer & mask)) {
						/* new k-mer */
						++n;
						/* break loop */
						kmer = key;
					} else {
						kmer = key_index ? key_index[++index] : key_index_l[++index];
					}
				}
			} else {
				/* new k-mer */
				++n;
			}
		}
		
		/* adjust dest to new size */
		fprintf(stderr, "# Adjusting size of new index.\n");
		dest->size++;
		if(dest->size <= n && (dest->mask + 1) <= (dest->size << 2)) { /* direct */
			dest->size = dest->mask + 1;
			if(dest->mlen <= 16) {
				free(dest->key_index);
				dest->key_index = 0;
			} else {
				free(dest->key_index_l);
				dest->key_index_l = 0;
			}
			if(exist) {
				free(dest->exist);
				dest->exist = 0;
			} else {
				free(dest->exist_l);
				dest->exist_l = 0;
			}
			exist_l = smalloc(dest->size * sizeof(long unsigned));
			dest->exist_l = exist_l--;
			
			/* init new exist */
			size = dest->size + 1;
			while(--size) {
				*++exist_l = null_index;
			}
			
			/* populate new exist with t1 */
			fprintf(stderr, "# Getting middlelayer signatures from first index.\n");
			if(t1->kmersize <= 16) {
				exist_t = t1->key_index - 1;
				exist_lt = 0;
			} else {
				exist_t = 0;
				exist_lt = t1->key_index_l - 1;
			}
			if(t1->v_index < UINT_MAX) {
				value_index = t1->value_index - 1;
				value_index_l = 0;
			} else {
				value_index = 0;
				value_index_l = t1->value_index_l - 1;
			}
			size = t1->n + 1;
			exist_l = dest->exist_l;
			while(--size) {
				index = exist_t ? *++exist_t : *++exist_lt;
				v_index = value_index ? *++value_index : *++value_index_l;
				exist_l[index] = MiddleLayer_search(middle, v_index);
			}
			
			/* merge with t2 */
			fprintf(stderr, "# Merging middlelayer signatures with second index.\n");
			dest->n += add_pairs(exist_l, t2->n, null_index, t1->v_index, middle, alternative, t2);
			
			/* adjust new alternative layer */
			AlternativeLayer_readjust(alternative);
			
			/* get v_index */
			fprintf(stderr, "# Getting new v_index.\n");
			v_index = getV_index(dest->exist_l, dest->size, null_index, t1->v_index + t2->v_index, middle, alternative);
			dest->v_index = v_index;
			
			/* adjust size of exist */
			dest->exist = adjustV_index(dest->exist_l, dest->size, v_index);
			dest->exist_l = (long unsigned *)(dest->exist_l);
		} else { /* hashmap */
			/* adjust size of hashmap */
			if(dest->size <= n) {
				dest->size <<= 1;
				free(dest->exist);
				free(dest->exist_l);
				if(n <= UINT_MAX) {
					dest->exist = malloc(dest->size * sizeof(unsigned));
					dest->exist_l = 0;
				} else {
					dest->exist = 0;
					dest->exist_l = malloc(dest->size * sizeof(long unsigned));
				}
			} else if(t1->n <= UINT_MAX && UINT_MAX < n) {
				dest->exist_l = realloc(dest->exist, dest->size * sizeof(long unsigned));
				dest->exist = 0;
			}
			
			/* adjust key_index */
			if(t1->n != n) {
				if(dest->key_index) {
					dest->key_index = realloc(dest->key_index, (n + 1) * sizeof(unsigned));
				} else {
					dest->key_index_l = realloc(dest->key_index_l, (n + 1) * sizeof(long unsigned));
				}
			}
			
			/* adjust value_index */
			dest->value_index_l = malloc(n * sizeof(long unsigned));
			dest->value_index = (unsigned *)(dest->value_index_l);
			
			/* check that everything was allocated properly */
			if((!dest->exist && !dest->exist_l) || (!dest->key_index && !dest->key_index_l) || (!dest->value_index && !dest->value_index)) {
				ERROR();
			}
			
			/* sort k-mer buckets from t1 based on dest, and dump on tmp */
			fprintf(stderr, "# Sorting buckets of first index w.r.t. new index.\n");
			dest->null_index = t1->null_index;
			dest->size--;
			hashMapKMA_sortbuckets(dest, t1);
			tmp_1 = hashMapKMA_dumpbuckets(dest);
			
			/* load k-mer from t2 */
			fprintf(stderr, "# Sorting buckets of second index w.r.t. new index.\n");
			if(dest->key_index) {
				memcpy(dest->key_index, t2->key_index, (t2->n + 1) * sizeof(unsigned));
			} else {
				memcpy(dest->key_index_l, t2->key_index_l, (t2->n + 1) * sizeof(long unsigned));
			}
			/* sort k-mer buckets from t2 based on dest, and dump on tmp */
			dest->null_index = t2->null_index;
			hashMapKMA_sortbuckets(dest, t2);
			tmp_2 = hashMapKMA_dumpbuckets(dest);
			
			/* merge in key- and value_indexes */
			fprintf(stderr, "# Merging signatures between indexes.\n");
			dest->n = n;
			dest->v_index = t1->v_index;
			dest->null_index = null_index;
			hashMapKMA_merge(dest, middle, alternative, tmp_1, tmp_2, t1->v_index + t2->v_index);
			dest->size++;
			fclose(tmp_1);
			fclose(tmp_2);
		}
	}
	
	return dest;
}

int merge(char *templatefilename, char *templatefilename1, char *templatefilename2) {
	
	int order;
	long unsigned org_split, alt_split;
	FILE *templatefile, *templatefile1, *templatefile2;
	HashMapKMA *t, *t1, *t2;
	MiddleLayer *middle, *alternative;
	
	/* init */
	order = 0;
	templatefile = sfopen(templatefilename, "wb");
	templatefile1 = sfopen(templatefilename1, "rb");
	templatefile2 = sfopen(templatefilename2, "rb");
	
	/* check compatability */
	t1 = smalloc(sizeof(HashMapKMA));
	t2 = smalloc(sizeof(HashMapKMA));
	loadPrefix(t1, templatefile1);
	loadPrefix(t2, templatefile2);
	if(hashMapKMA_compatible(t1, t2)) {
		fprintf(stderr, "# Hashmaps are compatible.\n");
		rewind(templatefile1);
		rewind(templatefile2);
	} else {
		fprintf(stderr, "Hashmaps are not compatible.\n");
		fclose(templatefile);
		fclose(templatefile1);
		fclose(templatefile2);
		free(t1);
		free(t2);
		return order;
	}
	
	/* load hashmaps on disk, and get largets hashmap on t1 */
	if(t1->size < t2->size) {
		order = 2;
		hashMapKMAmmap(t2, templatefile1);
		hashMapKMAmmap(t1, templatefile2);
	} else {
		order = 1;
		hashMapKMAmmap(t1, templatefile1);
		hashMapKMAmmap(t2, templatefile2);
	}
	
	/* create middlelayer with unique ids for each signature */
	fprintf(stderr, "# Adding middle layer.\n");
	middle = MiddleLayer_init(1048576);
	org_split = MiddleLayer_populate(middle, t1, 0);
	if(org_split != t1->v_index) {
		fprintf(stderr, "middle t1 (%lu) does not match v_index1 (%lu)\n", t1->v_index, org_split);
	}
	alt_split = MiddleLayer_populate(middle, t2, t1->v_index);
	if(alt_split != (t1->v_index + t2->v_index)) {
		fprintf(stderr, "middle t2 (%lu) does not match v_index2 (%lu)\n", t1->v_index + t2->v_index, alt_split);
	}
	MiddleLayer_readjust(middle); /* free up unused space */
	alternative = AlternativeLayer_init(1024);
	
	/* merge hashmaps */
	fprintf(stderr, "# Merging signatures.\n");
	t = merge_kmersignatures(t1, t2, middle, alternative);
	
	/* dump and free new hashmap */
	fprintf(stderr, "# Dumping new index.\n");
	hashMapKMA_dumpmerge(t, t1, t2, middle, alternative, templatefile);
	
	/* clean up */
	fprintf(stderr, "# Clean up new index.\n");
	hashMapKMA_munmap(t1);
	hashMapKMA_munmap(t2);
	MiddleLayer_free(middle);
	MiddleLayer_free(alternative);
	fclose(templatefile);
	fclose(templatefile1);
	fclose(templatefile2);
	
	return order;
}

int merge_lengths(char *outname, char *inname1, char *inname2) {
	
	int check;
	unsigned DB_size, n1, n2, *lengths;
	FILE *out, *in1, *in2;
	size_t size;
	
	/* init */
	out = sfopen(outname, "wb");
	in1 = sfopen(inname1, "rb");
	in2 = sfopen(inname2, "rb");
	
	/* get DB_size */
	sfread(&n1, sizeof(unsigned), 1, in1);
	sfread(&n2, sizeof(unsigned), 1, in2);
	DB_size = n1 + --n2;
	sfwrite(&DB_size, sizeof(unsigned), 1, out);
	
	/* merge lengths */
	lengths = smalloc(DB_size * sizeof(unsigned));
	sfread(lengths, sizeof(unsigned), n1, in1);
	sfseek(in2, sizeof(unsigned), SEEK_CUR); /* skip template zero */
	sfread(lengths + n1, sizeof(unsigned), n2, in2);
	sfwrite(lengths, sizeof(unsigned), DB_size, out);
	
	/* merge slengths */
	size = fread(lengths, sizeof(unsigned), n1, in1);
	check = fseek(in2, sizeof(unsigned), SEEK_CUR); /* skip template zero */
	size += fread(lengths + n1, sizeof(unsigned), n2, in2);
	if(!check && size == DB_size) {
		sfwrite(lengths, sizeof(unsigned), DB_size, out);
	}
	
	/* merge ulengths */
	size = fread(lengths, sizeof(unsigned), n1, in1);
	check = fseek(in2, sizeof(unsigned), SEEK_CUR); /* skip template zero */
	size += fread(lengths + n1, sizeof(unsigned), n2, in2);
	if(!check && size == DB_size) {
		sfwrite(lengths, sizeof(unsigned), DB_size, out);
	}
	
	/* clean up */
	free(lengths);
	fclose(in1);
	fclose(in2);
	return fclose(out);
}

int cat(char *outname, char *inname1, char *inname2) {
	
	char *buff;
	FILE *out, *in;
	size_t size;
	
	/* init */
	buff = smalloc(1048576);
	out = sfopen(outname, "wb");
	
	/* get first name */
	in = sfopen(inname1, "rb");
	while((size = fread(buff, 1, 1048576, in))) {
		cfwrite(buff, 1, size, out);
	}
	fclose(in);
	
	/* get second name */
	in = sfopen(inname2, "rb");
	while((size = fread(buff, 1, 1048576, in))) {
		cfwrite(buff, 1, size, out);
	}
	fclose(in);
	free(buff);
	
	return fclose(out);
}

static void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# kma merge merges two indexed kma databases.\n");
	fprintf(helpOut, "# %16s\t%-32s\t%s\n", "Options:", "Desc:", "Default:");
	fprintf(helpOut, "# %16s\t%-32s\t%s\n", "-o", "Output prefix", "");
	fprintf(helpOut, "# %16s\t%-32s\t%s\n", "-t_db", "Add to DB", "");
	fprintf(helpOut, "# %16s\t%-32s\t%s\n", "-s_db", "DB to merge", "");
	fprintf(helpOut, "# %16s\t%-32s\t%s\n", "-tmp", "Set directory for temporary files", "");
	fprintf(helpOut, "# %16s\t%-32s\t%s\n", "-v", "Version", "");
	fprintf(helpOut, "# %16s\t%-32s\t%s\n", "-h", "Shows this help message", "");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int merge_main(int argc, char *argv[]) {
	
	int args, order, o_len, t_len, s_len;
	char *outfilename, *templatefilename, *secondfilename, *tmpfilename;
	
	/* init */
	o_len = 0;
	t_len = 0;
	s_len = 0;
	outfilename = 0;
	templatefilename = 0;
	secondfilename = 0;
	
	/* PARSE COMMAND LINE OPTIONS */
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-o") == 0) {
			++args;
			if(args < argc) {
				o_len = strlen(argv[args]);
				outfilename = smalloc(o_len + 64);
				strcpy(outfilename, argv[args]);
			}
		} else if(strcmp(argv[args], "-t_db") == 0) {
			++args;
			if(args < argc) {
				t_len = strlen(argv[args]);
				templatefilename = smalloc(t_len + 64);
				strcpy(templatefilename, argv[args]);
			}
		} else if(strcmp(argv[args], "-s_db") == 0) {
			++args;
			if(args < argc) {
				s_len = strlen(argv[args]);
				secondfilename = smalloc(s_len + 64);
				strcpy(secondfilename, argv[args]);
			}
		} else if(strcmp(argv[args], "-tmp") == 0) {
			if(++args < argc) {
				if(argv[args][0] != '-') {
					if(argv[args][strlen(argv[args]) - 1] != '/') {
						fprintf(stderr, "Invalid output directory specified.\n");
						exit(1);
					}
					tmpF(argv[args]);
				} else {
					--args;
				}
			}
		} else if(strcmp(argv[args], "-v") == 0) {
			fprintf(stdout, "KMA_index-%s\n", KMA_VERSION);
			exit(0);
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else {
			fprintf(stderr, "# Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, "# Printing help message:\n");
			helpMessage(1);
		}
		++args; 
	}
	
	/* check if all input was parsed */
	if(!o_len || !t_len || !s_len) {
		fprintf(stderr, "Insufficient number of agruments parsed.\n");
		helpMessage(1);
	}
	if(strcmp(templatefilename, secondfilename) == 0) {
		fprintf(stderr, "Indexes to merge cannot be the same.\n");
		return 1;
	}
	
	/* merge *.comp.b */
	fprintf(stderr, "# Merging *.comp.b\n");
	strcpy(outfilename + o_len, ".comp.b");
	strcpy(templatefilename + t_len, ".comp.b");
	strcpy(secondfilename + s_len, ".comp.b");
	order = merge(outfilename, templatefilename, secondfilename);
	if(order == 2) { /* switch order of merged DBs */
		tmpfilename = templatefilename;
		templatefilename = secondfilename;
		secondfilename = tmpfilename;
		order = t_len;
		t_len = s_len;
		s_len = order;
	}
	
	/* merge *.length.b */
	fprintf(stderr, "# Merging *.length.b\n");
	strcpy(outfilename + o_len, ".length.b");
	strcpy(templatefilename + t_len, ".length.b");
	strcpy(secondfilename + s_len, ".length.b");
	merge_lengths(outfilename, templatefilename, secondfilename);
	
	/* merge *.seq.b */
	fprintf(stderr, "# Merging *.seq.b\n");
	strcpy(outfilename + o_len, ".seq.b");
	strcpy(templatefilename + t_len, ".seq.b");
	strcpy(secondfilename + s_len, ".seq.b");
	cat(outfilename, templatefilename, secondfilename);
	
	/* merge *.name */
	fprintf(stderr, "# Merging *.name\n");
	strcpy(outfilename + o_len, ".name");
	strcpy(templatefilename + t_len, ".name");
	strcpy(secondfilename + s_len, ".name");
	cat(outfilename, templatefilename, secondfilename);
	
	return 0;
}
