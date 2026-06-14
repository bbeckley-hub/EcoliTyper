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
#include <stdlib.h>
#include <zlib.h>
#include "filebuff.h"
#include "pherror.h"
#include "qseqs.h"

int (*buffFileBuff)(FileBuff *) = &BuffgzFileBuff;

int BuffgzFileBuff(FileBuff *dest) {
	
	int status;
	z_stream *strm;
	
	/* check compressed buffer, and load it */
	strm = dest->strm;
	if(strm->avail_in == 0) {
		strm->avail_in = fread(dest->inBuffer, 1, dest->buffSize, dest->file);
		strm->next_in = (unsigned char*) dest->inBuffer;
		if(strm->avail_in == 0) {
			dest->bytes = 0;
			dest->next = dest->buffer;
			return 0;
		}
	}
	
	/* reset uncompressed buffer */
	strm->avail_out = dest->buffSize;
	strm->next_out = (unsigned char*) dest->buffer;
	
	/* uncompress buffer */
	status = inflate(strm, Z_NO_FLUSH);
	dest->z_err = status;
	
	/* concatenated file */
	if(status == Z_STREAM_END && strm->avail_out == dest->buffSize) {
		inflateReset(strm);
		return BuffgzFileBuff(dest);
	}
	
	if(status == Z_OK || status == Z_STREAM_END) {
		dest->bytes = dest->buffSize - strm->avail_out;
		dest->next = dest->buffer;
		if(status == Z_OK && dest->bytes == 0) {
			return BuffgzFileBuff(dest);
		}
	} else {
		dest->bytes = 0;
		dest->next = dest->buffer;
		fprintf(stderr, "Gzip error %d\n", status);
		errno |= status;
	}
	
	return dest->bytes;
}

void init_gzFile(FileBuff *inputfile) {
	
	int status;
	unsigned char *tmp;
	z_stream *strm;
	
	/* set inBuffer, for compressed format */
	if(inputfile->inBuffer) {
		tmp = inputfile->buffer;
		inputfile->buffer = inputfile->inBuffer;
		inputfile->inBuffer = tmp;
	} else {
		inputfile->inBuffer = inputfile->buffer;
		inputfile->buffer = smalloc(inputfile->buffSize + 1);
		inputfile->buffer[inputfile->buffSize] = 0;
	}
	inputfile->next = inputfile->buffer;
	
	/* set the compressed stream */
	strm = inputfile->strm;
	if(!strm && !(strm = malloc(sizeof(z_stream)))) {
		ERROR();
	}
	strm->zalloc = Z_NULL;
	strm->zfree  = Z_NULL;
	strm->opaque = Z_NULL;
	status = inflateInit2(strm, 15 | ENABLE_ZLIB_GZIP);
	if(status < 0) {
		fprintf(stderr, "Gzip error %d\n", status);
		exit(status);
	}
	strm->next_in = inputfile->inBuffer;
	strm->avail_in = inputfile->bytes;
	strm->avail_out = 0;
	inputfile->strm = strm;
	inputfile->z_err = Z_OK;
	
	inputfile->bytes = BuffgzFileBuff(inputfile);
}

FileBuff * setFileBuff(int buffSize) {
	
	FileBuff *dest;
	
	dest = smalloc(sizeof(FileBuff));
	dest->file = 0;
	dest->inBuffer = 0;
	dest->strm = 0;
	dest->buffSize = buffSize;
	dest->buffer = smalloc(buffSize);
	dest->next = dest->buffer;
	
	return dest;
}

void openFileBuff(FileBuff *dest, char *filename, char *mode) {
	dest->file = sfopen(filename, mode);
}

void closeFileBuff(FileBuff *dest) {
	fclose(dest->file);
	dest->file = 0;
}


void gzcloseFileBuff(FileBuff *dest) {
	
	int status;
	if((status = inflateEnd(dest->strm)) != Z_OK) {
		fprintf(stderr, "Gzip error %d\n", status);
		errno |= status;
	}
	if(dest->z_err != Z_STREAM_END) {
		fprintf(stderr, "Unexpected end of file\n");
		errno |= status;
	}
	fclose(dest->file);
	dest->file = 0;
	dest->strm->avail_out = 0;
}

void destroyFileBuff(FileBuff *dest) {
	free(dest->buffer);
	free(dest->inBuffer);
	free(dest->strm);
	free(dest);
}

int buff_FileBuff(FileBuff *dest) {
	dest->bytes = fread(dest->buffer, 1, dest->buffSize, dest->file);
	dest->next = dest->buffer;
	return dest->bytes;
}

z_stream * strm_init() {
	
	z_stream *strm;
	int status;
	
	strm = smalloc(sizeof(z_stream));
	strm->zalloc = Z_NULL;
	strm->zfree  = Z_NULL;
	strm->opaque = Z_NULL;
	
	status = deflateInit2(strm, 1, Z_DEFLATED, 31 | GZIP_ENCODING, 9, Z_FILTERED);
	if(status < 0) {
		fprintf(stderr, "Gzip error %d\n", status);
		exit(status);
	}
	
	return strm;
}

FileBuff * gzInitFileBuff(int size) {
	
	FileBuff *dest;
	
	dest = smalloc(sizeof(FileBuff));
	dest->file = 0;
	dest->bytes = size;
	dest->buffSize = size;
	dest->strm = strm_init();
	dest->buffer = smalloc(size);
	dest->inBuffer = smalloc(size);
	dest->next = dest->buffer;
	
	return dest;
}

void resetGzFileBuff(FileBuff *dest, int size) {
	
	dest->bytes = size;
	dest->buffSize = size;
	free(dest->buffer);
	free(dest->inBuffer);
	dest->buffer = smalloc(size);
	dest->inBuffer = smalloc(size);
	dest->next = dest->buffer;
}

void writeGzFileBuff(FileBuff *dest) {
	
	int check = Z_OK;
	z_stream *strm = dest->strm;
	strm->avail_in = dest->buffSize - dest->bytes;
	strm->next_in = dest->buffer;
	strm->avail_out = 0;
	
	while(strm->avail_out == 0 && check != Z_STREAM_END) {
		strm->avail_out = dest->buffSize;
		strm->next_out = dest->inBuffer;
		check = deflate(strm, Z_NO_FLUSH);
		sfwrite(dest->inBuffer, 1, dest->buffSize - strm->avail_out, dest->file);
	}
	dest->bytes = dest->buffSize;
	dest->next = dest->buffer;
}

void closeGzFileBuff(FileBuff *dest) {
	
	int check = Z_OK;
	z_stream *strm = dest->strm;
	strm->avail_in = dest->buffSize - dest->bytes;
	strm->next_in = dest->buffer;
	strm->avail_out = 0;
	
	while(strm->avail_out == 0 && check != Z_STREAM_END) {
		strm->avail_out = dest->buffSize;
		strm->next_out = dest->inBuffer;
		check = deflate(strm, Z_FINISH);
		sfwrite(dest->inBuffer, 1, dest->buffSize - strm->avail_out, dest->file);
	}
	deflateEnd(strm);
	fclose(dest->file);
}

void destroyGzFileBuff(FileBuff *dest) {
	
	closeGzFileBuff(dest);
	free(dest->strm);
	free(dest->buffer);
	free(dest->inBuffer);
	free(dest);
}
