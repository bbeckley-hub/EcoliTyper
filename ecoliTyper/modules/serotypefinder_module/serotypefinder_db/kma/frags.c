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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "filebuff.h"
#include "frags.h"
#include "pherror.h"
#include "qseqs.h"
#include "threader.h"
#include "tmp.h"

FILE * printFrags(Frag **alignFrags, int DB_size) {
	
	int i;
	FILE *OUT;
	Frag *alignFrag, *next;
	
	if(!(OUT = tmpF(0))) {
		fprintf(stderr, "Could not create tmp files.\n");
		ERROR();
	}
	for(i = 0; i < DB_size; ++i) {
		if(alignFrags[i]) {
			for(alignFrag = alignFrags[i]; alignFrag != 0; alignFrag = next) {
				next = alignFrag->next;
				
				sfwrite(&i, sizeof(int), 1, OUT);
				sfwrite(alignFrag->buffer, sizeof(int), 7, OUT);
				sfwrite(alignFrag->qseq, 1, alignFrag->buffer[0], OUT);
				sfwrite(alignFrag->header, 1, alignFrag->buffer[5], OUT);
				
				free(alignFrag->qseq);
				free(alignFrag->header);
				free(alignFrag);
			}
			alignFrags[i] = 0;
		}
	}
	sfwrite(&(int){-1}, sizeof(int), 1, OUT);
	fflush(OUT);
	rewind(OUT);
	
	return OUT;
}

void updateAllFrag(unsigned char *qseq, int q_len, int bestHits, int best_read_score, int *best_start_pos, int *best_end_pos, int *bestTemplates, Qseqs *header, FileBuff *dest) {
	
	static volatile int Lock = 0;
	volatile int *lock = &Lock;
	int i, check, avail;
	char *update;
	const char bases[6] = "ACGTN-";
	
	if(!q_len) {
		return;
	}
	
	lock(lock);
	check = q_len;
	avail = dest->bytes;
	
	if(avail < check) {
		writeGzFileBuff(dest);
		
		/* seq is too big, reallocate buffer */
		if(dest->bytes < check) {
			resetGzFileBuff(dest, check << 1);
		}
		avail = dest->bytes;
	}
	
	/* copy seq */
	update = (char *) dest->next;
	++q_len;
	--qseq;
	while(--q_len) {
		*update++ = bases[*++qseq];
	}
	avail -= check;
	
	check = 33;
	if(avail < check) {
		dest->bytes = avail;
		writeGzFileBuff(dest);
		avail = dest->bytes;
		update = (char *) dest->next;
	}
	check = sprintf(update, "\t%d\t%d\t%d", bestHits, best_read_score, *best_start_pos);
	avail -= check;
	update += check;
	
	for(i = 1; i < bestHits; ++i) {
		if(avail < 11) {
			dest->bytes = avail;
			writeGzFileBuff(dest);
			avail = dest->bytes;
			update = (char *) dest->next;
		}
		check = sprintf(update, ",%d", best_start_pos[i]);
		avail -= check;
		update += check;
	}
	
	if(avail < 11) {
		dest->bytes = avail;
		writeGzFileBuff(dest);
		avail = dest->bytes;
		update = (char *) dest->next;
	}
	check = sprintf(update, "\t%d", *best_end_pos);
	avail -= check;
	update += check;
	for(i = 1; i < bestHits; ++i) {
		if(avail < 11) {
			dest->bytes = avail;
			writeGzFileBuff(dest);
			avail = dest->bytes;
			update = (char *) dest->next;
		}
		check = sprintf(update, ",%d", best_end_pos[i]);
		avail -= check;
		update += check;
	}
	
	if(avail < 11) {
		dest->bytes = avail;
		writeGzFileBuff(dest);
		avail = dest->bytes;
		update = (char *) dest->next;
	}
	check = sprintf(update, "\t%d", *bestTemplates);
	avail -= check;
	update += check;
	for(i = 1; i < bestHits; ++i) {
		if(avail < 12) {
			dest->bytes = avail;
			writeGzFileBuff(dest);
			avail = dest->bytes;
			update = (char *) dest->next;
		}
		check = sprintf(update, ",%d", bestTemplates[i]);
		avail -= check;
		update += check;
	}
	
	/* q-bound */
	if(2 * sizeof(int) + 1 < header->len && header->seq[header->len - 2 * sizeof(int) - 1] == 0) {
		i = header->len - 2 * sizeof(int) - 1;
	} else {
		i = header->len;
	}
	check = i + 1;
	if(avail < check) {
		dest->bytes = avail;
		writeGzFileBuff(dest);
		avail = dest->bytes;
		update = (char *) dest->next;
	}
	*update++ = '\t';
	memcpy(update, header->seq, i--);
	update += i;
	*update++ = '\n';
	
	dest->bytes = avail - check;
	dest->next = (unsigned char *) update;
	unlock(lock);
}
