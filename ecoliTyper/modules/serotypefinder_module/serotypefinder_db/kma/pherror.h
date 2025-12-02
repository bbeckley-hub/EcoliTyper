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
#include <stdio.h>
#include <string.h>
#include <unistd.h>

/* errno error message */
#define ERROR() fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno)); exit(errno);

/* succes or exit */
#define sfwrite(ptr, size, nmemb, stream) if(fwrite(ptr, size, nmemb, stream) != nmemb) {if(errno) {ERROR();} else {fprintf(stderr, "Writing error.\n"); exit(1);}}
#define sfread(ptr, size, nmemb, stream) if(fread(ptr, size, nmemb, stream) != nmemb) {if(errno) {ERROR();} else {fprintf(stderr, "Reading error.\n"); exit(1);}}
#define sfseek(stream, offset, whence) if(fseek(stream, offset, whence)) {if(errno) {ERROR();} else {fprintf(stderr, "fseek error.\n"); exit(1);}}
void * smalloc(const size_t size);
FILE * sfopen(const char *filename, const char *mode);

/* cyclic */
void cfread(void *src, size_t size, size_t nmemb, FILE *stream);
void cfwrite(const void *src, size_t size, size_t nmemb, FILE *stream);
