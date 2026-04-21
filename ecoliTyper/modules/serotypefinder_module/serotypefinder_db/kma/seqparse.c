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

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "filebuff.h"
#include "pherror.h"
#include "qseqs.h"
#include "seqparse.h"

int openAndDetermine(FileBuff *inputfile, char *filename) {
	
	unsigned FASTQ;
	short unsigned *check;
	
	/* determine filetype and open it */
	FASTQ = 0;
	if(*filename == '-' && strcmp(filename + 1, "-") == 0) {
		inputfile->file = stdin;
	} else {
		openFileBuff(inputfile, filename, "rb");
	}
	inputfile->buffer[0] = 0;
	if(buff_FileBuff(inputfile)) {
		check = (short unsigned *) inputfile->buffer;
		if(*check == 35615) {
			FASTQ = 4;
			init_gzFile(inputfile);
			buffFileBuff = &BuffgzFileBuff;
		} else {
			buffFileBuff = &buff_FileBuff;
		}
	}
	if(inputfile->buffer[0] == '@') { //FASTQ
		FASTQ |= 1;
	} else if(inputfile->buffer[0] == '>') { //FASTA
		FASTQ |= 2;
	} else if(inputfile->buffer[0] == 0) {
		fprintf(stderr, "Empty file:\t%s\n", filename);
	} else {
		fprintf(stderr, "Cannot determine format of file:\t%s\n", filename);
		errno |= 1;
	}
	
	return FASTQ;
}

int FileBuffgetFsa(FileBuff *src, Qseqs *header, Qseqs *qseq, char *trans) {
	
	unsigned char *buff, *seq;
	int size, avail;
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	qseq->len = 0;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get header */
	seq = header->seq;
	size = header->size;
	while((*seq++ = *buff++) != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		if(--size == 0) {
			size = header->size;
			header->size <<= 1;
			header->seq = realloc(header->seq, header->size);
			if(!header->seq) {
				ERROR();
			}
			seq = header->seq + size;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	/* chomp header */
	while(isspace(*--seq)) {
		++size;
	}
	*++seq = 0;
	header->len = header->size - size + 1;
	/* get qseq */
	seq = qseq->seq;
	size = qseq->size;
	while(*buff != '>') {
		*seq = trans[*buff++];
		if(((*seq) >> 3) == 0) {
			if(--size == 0) {
				size = qseq->size;
				qseq->size <<= 1;
				qseq->seq = realloc(qseq->seq, qseq->size);
				if(!qseq->seq) {
					ERROR();
				}
				seq = qseq->seq + size;
			} else {
				++seq;
			}
		}
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				/* chomp header */
				while(*--seq == 8) {
					++size;
				}
				*++seq = 0;
				qseq->len = qseq->size - size;
				
				src->bytes = 0;
				src->next = buff;
				return 1;
			}
			buff = src->buffer;
		}
	}
	/* chomp header */
	while(*--seq == 8) {
		++size;
	}
	*++seq = 0;
	qseq->len = qseq->size - size;
	
	src->bytes = avail;
	src->next = buff;
	
	return 1;
}

int FileBuffgetFsaSeq(FileBuff *src, Qseqs *qseq, char *trans) {
	
	unsigned char *buff, *seq;
	int size, avail;
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	qseq->len = 0;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* skip header */
	while(*buff++ != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get qseq */
	seq = qseq->seq;
	size = qseq->size;
	while(*buff != '>') {
		*seq = trans[*buff++];
		if(((*seq) >> 3) == 0) {
			if(--size == 0) {
				size = qseq->size;
				qseq->size <<= 1;
				qseq->seq = realloc(qseq->seq, qseq->size);
				if(!qseq->seq) {
					ERROR();
				}
				seq = qseq->seq + size;
			} else {
				++seq;
			}
		}
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				/* chomp header */
				while(*--seq == 8) {
					++size;
				}
				*++seq = 0;
				qseq->len = qseq->size - size;
				
				src->bytes = 0;
				src->next = buff;
				return 1;
			}
			buff = src->buffer;
		}
	}
	
	/* chomp qseq */
	while(*--seq == 8) {
		++size;
	}
	*++seq = 0;
	qseq->len = qseq->size - size;
	
	src->bytes = avail;
	src->next = buff;
	
	return 1;
}

int FileBuffgetFq(FileBuff *src, Qseqs *header, Qseqs *qseq, Qseqs *qual, char *trans) {
	
	unsigned char *buff, *seq;
	int size, avail;
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	qseq->len = 0;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	if(*buff != '@') {
		fprintf(stderr, "Malformed input.\n");
		errno |= 1;
		return 0;
	}
	
	/* get header */
	seq = header->seq;
	size = header->size;
	while((*seq++ = *buff++) != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		if(--size == 0) {
			size = header->size;
			header->size <<= 1;
			header->seq = realloc(header->seq, header->size);
			if(!header->seq) {
				ERROR();
			}
			seq = header->seq + size;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	/* chomp header */
	while(isspace(*--seq)) {
		++size;
	}
	*++seq = 0;
	header->len = header->size - size + 1;
	
	/* get qseq */
	seq = qseq->seq;
	size = qseq->size;
	while((*seq++ = trans[*buff++]) != 16) {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		if(--size == 0) {
			size = qseq->size;
			qseq->size <<= 1;
			qseq->seq = realloc(qseq->seq, qseq->size);
			if(!qseq->seq) {
				ERROR();
			}
			seq = qseq->seq + size;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	/* chomp header */
	while(*--seq == 8) {
		++size;
	}
	*++seq = 0;
	qseq->len = qseq->size - size;
	
	/* skip info */
	while(*buff++ != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get quality */
	if(qual->size != qseq->size) {
		qual->size = qseq->size;
		free(qual->seq);
		qual->seq = malloc(qual->size);
		if(!qual->seq) {
			ERROR();
		}
	}
	qual->len = qseq->len;
	if(qual->len < avail) {
		memcpy(qual->seq, buff, qual->len);
		avail -= qual->len;
		buff += qual->len;
	} else {
		seq = qual->seq;
		size = qual->len;
		while(size) {
			if(size < avail) {
				memcpy(seq, buff, size);
				size = 0;
				avail -= size;
				buff += size;
			} else {
				memcpy(seq, buff, avail);
				size -= avail;
				seq += avail;
				if((avail = buffFileBuff(src)) == 0) {
					return 0;
				}
				buff = src->buffer;
			}
		}
	}
	qual->seq[qual->len] = 0;
	
	/* skip newline */
	while(*buff++ != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				/* warning */
				fprintf(stderr, "Truncated file.\n");
				return 0;
			}
			buff = src->buffer;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 1;
		}
		src->bytes = avail;
		buff = src->buffer;
	}
	
	src->bytes = avail;
	src->next = buff;
	
	return 1;
}

int FileBuffgetFqSeq(FileBuff *src, Qseqs *qseq, Qseqs *qual, char *trans) {
	
	unsigned char *buff, *seq;
	int size, avail;
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	qseq->len = 0;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	if(*buff != '@') {
		fprintf(stderr, "Malformed input.\n");
		errno |= 1;
		return 0;
	}
	
	/* skip header */
	while(*buff++ != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get qseq */
	seq = qseq->seq;
	size = qseq->size;
	while((*seq++ = trans[*buff++]) != 16) {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		if(--size == 0) {
			size = qseq->size;
			qseq->size <<= 1;
			qseq->seq = realloc(qseq->seq, qseq->size);
			if(!qseq->seq) {
				ERROR();
			}
			seq = qseq->seq + size;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	/* chomp header */
	while(*--seq == 8) {
		++size;
	}
	*++seq = 0;
	qseq->len = qseq->size - size;
	
	/* skip info */
	while(*buff++ != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get quality */
	if(qual->size != qseq->size) {
		qual->size = qseq->size;
		free(qual->seq);
		qual->seq = malloc(qual->size);
		if(!qual->seq) {
			ERROR();
		}
	}
	qual->len = qseq->len;
	if(qual->len < avail) {
		memcpy(qual->seq, buff, qual->len);
		avail -= qual->len;
		buff += qual->len;
	} else {
		seq = qual->seq;
		size = qual->len;
		while(size) {
			if(size < avail) {
				memcpy(seq, buff, size);
				size = 0;
				avail -= size;
				buff += size;
			} else {
				memcpy(seq, buff, avail);
				size -= avail;
				seq += avail;
				if((avail = buffFileBuff(src)) == 0) {
					return 0;
				}
				buff = src->buffer;
			}
		}
	}
	qual->seq[qual->len] = 0;
	
	/* skip newline */
	while(*buff++ != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 1;
		}
		src->bytes = avail;
		buff = src->buffer;
	}
	
	src->bytes = avail;
	src->next = buff;
	
	return 1;
}

int getPhredFileBuff(FileBuff *dest) {
	
	int seek, avail, scale, len, maxlen;
	unsigned char *buff;
	
	avail = dest->bytes;
	scale = 33;
	maxlen = 0;
	buff = dest->next;
	
	while(avail) {
		seek = 3;
		while(seek && --avail) {
			if(*++buff == '\n') {
				--seek;
			}
		}
		
		len = 0;
		seek = avail ? 1 : 0;
		while(seek && --avail) {
			if(*++buff == '\n') {
				seek = 0;
			} else if(*buff < 33) {
				return 0;
			} else if(53 < *buff && *buff < 59) {
				return 33;
			} else if(94 < *buff) {
				scale = 64;
			}
			++len;
		}
		if(maxlen < len) {
			maxlen = len;
		}
	}
	
	return maxlen <= 301 ? scale : 33;
}
