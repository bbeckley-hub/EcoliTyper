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
#include <stdlib.h>
#include <stdio.h>
#include <sys/mman.h>
#include <sys/param.h>
#include "matrix.h"
#include "pherror.h"
#include "tmp.h"

Matrix * matrix_init(unsigned size) {
	
	int i, **ptr, *src;
	long unsigned Size;
	Matrix *dest;
	
	dest = smalloc(sizeof(Matrix));
	dest->n = 0;
	dest->size = size;
	dest->mat = smalloc(size * sizeof(int *));
	Size = size;
	Size *= size;
	Size *= sizeof(int);
	*(dest->mat) = smalloc(Size);
	
	/* set matrix rows */
	ptr = dest->mat;
	src = *ptr;
	i = size;
	*ptr = src;
	while(--i) {
		src += size;
		*++ptr = src;
	}
	
	return dest;
}

Matrix * ltdMatrix_init(unsigned size) {
	
	int i, **ptr, *src;
	long unsigned Size;
	Matrix *dest;
	
	dest = smalloc(sizeof(Matrix));
	dest->n = 0;
	dest->size = size;
	dest->mat = smalloc(size * sizeof(int *));
	Size = size;
	Size *= (size - 1);
	Size *= (sizeof(double) / 2);
	*(dest->mat) = calloc(Size, 1);
	if(!*(dest->mat)) {
		ERROR();
	}
	
	/* set matrix rows */
	ptr = dest->mat;
	src = *ptr;
	i = 0;
	*ptr++ = src;
	while(--size) {
		*ptr++ = src + i;
		src += i++;
	}
	
	return dest;
}

Matrix * ltdMatrix_minit(long unsigned size) {
	
	int i, n, **ptr, *src;
	FILE *tmp;
	Matrix *dest;
	
	dest = smalloc(sizeof(Matrix));
	dest->n = 0;
	dest->size = size;
	dest->mat = smalloc(size * sizeof(int *));
	n = size;
	size = size * (size - 1) * sizeof(int) / 2;
	
	/* get matrix */
	tmp = tmpF(0);
	if(fseek(tmp, size - 1, SEEK_SET) || putc(0, tmp) == EOF) {
		ERROR();
	}
	fflush(tmp);
	sfseek(tmp, 0, SEEK_SET);
	*(dest->mat) = mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
	if(*(dest->mat) == MAP_FAILED) {
			ERROR();
	}
	posix_madvise(*(dest->mat), size, POSIX_MADV_SEQUENTIAL);
	
	/* set matrix rows */
	ptr = dest->mat;
	src = *ptr;
	i = 0;
	*ptr++ = src;
	while(--n) {
		*ptr++ = src + i;
		src += i++;
	}
	
	return dest;
}

void ltdMatrix_realloc(Matrix *src, unsigned size) {
	
	int i, **ptr, *mat;
	long unsigned Size;
	
	Size = size;
	Size *= (size - 1);
	Size *= (sizeof(int) / 2);
	*(src->mat) = realloc(*(src->mat), Size);
	src->mat = realloc(src->mat, size * sizeof(int *));
	if(!src->mat || !*(src->mat)) {
		ERROR();
	}
	src->size = size;
	
	/* set matrix rows */
	ptr = src->mat;
	mat = *ptr;
	i = 0;
	*ptr++ = mat;
	while(--size) {
		*ptr++ = mat + i;
		mat += i++;
	}
}

void Matrix_destroy(Matrix *src) {
	
	free(*(src->mat));
	free(src->mat);
	free(src);
}

void Matrix_mdestroy(Matrix *src) {
	
	long unsigned size;
	
	size = src->size;
	size *= (src->size - 1);
	size *= sizeof(int) / 2;
	msync(*(src->mat), size, MS_SYNC);
	munmap(*(src->mat), size);
	free(src->mat);
	free(src);
}

void ltdMatrix_popArrange(Matrix *mat, unsigned pos) {
	
	int i, n, *dest, *src;
	
	n = --mat->n;
	if(pos != n) {
		/* row to be emptied */
		dest = mat->mat[pos];
		/* row to be moved up */
		src = mat->mat[n];
		
		/* copy last row into "first" row */
		i = pos + 1;
		while(--i) {
			*dest++ = *src++;
		}
		
		/* tilt remaining part of last row into column "pos" */
		i = pos;
		++src;
		while(++i < n) {
			mat->mat[i][pos] = *src++;
		}
	}
}

int ltdMatrix_add(Matrix *src) {
	
	int i, *ptr;
	
	/* realloc */
	if(++src->n == src->size) {
		ltdMatrix_realloc(src, src->size << 1);
	}
	
	/* init new row */
	i = src->n;
	ptr = src->mat[i - 1] - 1;
	while(--i) {
		*++ptr = 0;
	}
	
	return src->n - 1;
}
