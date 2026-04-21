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


#include <stdlib.h>
#include "compkmers.h"
#include "pherror.h"

void allocCompKmers(CompKmers *compressor, int size) {
	
	compressor->n = 0;
	compressor->size = size;
	compressor->kmers = smalloc(size * sizeof(long unsigned));
}

void reallocCompKmers(CompKmers *compressor, int size) {
	
	compressor->kmers = realloc(compressor->kmers, size * sizeof(long unsigned));
	if(!compressor->kmers) {
		ERROR();
	}
	compressor->size = size;
	
}

int pushCompKmers(CompKmers *compressor, char *qseq, int kmersize) {
	
	int i;
	long unsigned key = 0;
	
	for(i = 0; i < kmersize; ++i) {
		if(qseq[i] == 4) {
			return 0;
		} else {
			key = (key << 2) | qseq[i];
		}
	}
	
	compressor->kmers[compressor->n] = key;
	return 1;
}
