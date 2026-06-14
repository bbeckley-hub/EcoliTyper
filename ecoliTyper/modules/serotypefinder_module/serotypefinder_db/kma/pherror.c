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
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "pherror.h"

void * smalloc(const size_t size) {
	
	void *dest = malloc(size);
	if(!dest) {
		ERROR();
	}
	
	return dest;
}

FILE * sfopen(const char *filename, const char *mode) {
	
	FILE *file = fopen(filename, mode);
	if(!file) {
		fprintf(stderr, "Filename:\t%s\n", filename);
		ERROR();
	}
	
	return file;
}

void cfread(void *src, size_t size, size_t nmemb, FILE *stream) {
	
	unsigned char *ptr = src;
	nmemb = fread(ptr, 1, (size *= nmemb), stream);
	if((size -= nmemb)) {
		ptr += nmemb;
		while(size) {
			if(0 < (nmemb = read(fileno(stream), ptr, size))) {
				size -= nmemb;
				ptr += nmemb;
			} else if(nmemb == 0 || (nmemb == -1 && (errno & EAGAIN))) {
				usleep(1000);
			} else {
				ERROR();
			}
		}
	}
}

void cfwrite(const void *src, size_t size, size_t nmemb, FILE *stream) {
	
	unsigned char *ptr;
	
	size *= nmemb;
	ptr = (unsigned char *)(src);
	while(size) {
		nmemb = fwrite(ptr, 1, size, stream);
		if(nmemb == 0) {
			ERROR();
		}
		size -= nmemb;
		ptr += nmemb;
	}
	
}
