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

#include <stdio.h>
#include <zlib.h>

#ifndef FILEBUFF
typedef struct fileBuff FileBuff;
struct fileBuff {
	int bytes;
	int buffSize;
	unsigned char *buffer;
	unsigned char *inBuffer;
	unsigned char *next;
	FILE *file;
	z_stream *strm;
	int z_err;
};
#define FILEBUFF 1
#define CHUNK 1048576
#define GZIP_ENCODING 16
#define ENABLE_ZLIB_GZIP 32
#endif

/* pointer to load buffer from a regular or gz file stream */
extern int (*buffFileBuff)(FileBuff *);
int BuffgzFileBuff(FileBuff *dest);
void init_gzFile(FileBuff *inputfile);
FileBuff * setFileBuff(int buffSize);
void openFileBuff(FileBuff *dest, char *filename, char *mode);
void closeFileBuff(FileBuff *dest);
void gzcloseFileBuff(FileBuff *dest);
void destroyFileBuff(FileBuff *dest);
int buff_FileBuff(FileBuff *dest);
z_stream * strm_init();
FileBuff * gzInitFileBuff(int size);
void resetGzFileBuff(FileBuff *dest, int size);
void writeGzFileBuff(FileBuff *dest);
void closeGzFileBuff(FileBuff *dest);
void destroyGzFileBuff(FileBuff *dest);
