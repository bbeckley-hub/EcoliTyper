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
#include <stdlib.h>
#include "pherror.h"
#include "seqmenttree.h"

SeqmentTree * initializeSeqmentTree(const long unsigned size) {
	
	SeqmentTree *dest;
	
	dest = smalloc(sizeof(SeqmentTree));
	dest->n = 0;
	dest->root = smalloc((dest->size = size) * sizeof(SeqmentTrees));
	
	return dest;
}

SeqmentTree * initSeqmentTree(SeqmentTree *src, const unsigned start, const unsigned end) {
	
	SeqmentTrees *node;
	
	if(!src) {
		src = initializeSeqmentTree(64);
	}
	
	src->n = 1;
	node = src->root;
	node->start = start;
	node->end = end;
	node->covered = end - start;
	node->branch[0] = 0;
	node->branch[1] = 0;
	
	return src;
}

SeqmentTrees * rcpSeqmentTree(SeqmentTree *src, SeqmentTrees *branch) {
	
	SeqmentTrees *dest;
	
	/* cp branch */
	dest = (src->root + src->n++);
	*dest = *branch;
	
	/* add subtree of the branch */
	if(*(branch->branch)) {
		branch->branch[0] = rcpSeqmentTree(src, branch->branch[0]);
		branch->branch[1] = rcpSeqmentTree(src, branch->branch[1]);
	}
	
	return dest;
}

void reallocSeqmentTree(SeqmentTree *src, const long unsigned size) {
	
	SeqmentTrees *org;
	
	/* save old tree */
	org = src->root;
	
	/* allocate new tree */
	src->size = size;
	src->root = smalloc(src->size * sizeof(SeqmentTrees));
	
	/* cp old tree into new */
	src->n = 0;
	src->root = rcpSeqmentTree(src, org);
	
	/* clean up */
	free(org);
}

void resizeSeqmentTree(SeqmentTree *src) {
	
	SeqmentTrees *org;
	
	/* save old tree */
	org = src->root;
	
	/* allocate new tree */
	src->size <<= 1;
	src->root = smalloc(src->size * sizeof(SeqmentTrees));
	
	/* cp old tree into new */
	src->n = 0;
	src->root = rcpSeqmentTree(src, org);
	
	/* clean up */
	free(org);
}

unsigned addSeqmentTrees(SeqmentTrees *root, SeqmentTrees *node) {
	
	unsigned pos, covered;
	SeqmentTrees *bud;
	
	if(*(root->branch)) { /* search */
		/* adjust limits of root */
		if(node->start < root->start && root->end < node->end) {
			root->start = node->start;
			root->end = node->end;
			root->covered = node->covered;
			node->covered = 0;
			*(root->branch) = 0;
			return root->covered;
		} else if(root->end < node->end) {
			root->end = node->end;
		} else if(node->start < root->start) {
			root->start = node->start;
		}
		
		/* search tree */
		pos = root->branch[1]->start;
		if(node->end < pos) { /* left */
			root->covered = root->branch[1]->covered + addSeqmentTrees(root->branch[0], node);
		} else if(pos <= node->start) { /* right */
			root->covered = root->branch[0]->covered + addSeqmentTrees(root->branch[1], node);
		} else { /* split */
			/* calculate right side */
			pos = node->start;
			node->start = root->branch[0]->end + 1;
			node->covered = node->end - node->start;
			covered = addSeqmentTrees(root->branch[1], node);
			
			/* calculate left side */
			node->start = pos;
			node->end = root->branch[0]->end;
			node->covered = node->end - node->start;
			root->covered = covered + addSeqmentTrees(root->branch[0], node);
		}
	} else if(node->end < root->start || root->end < node->start) { /* new leaf */
		/* create and grow bud */
		bud = node + 1;
		bud->start = root->start;
		bud->end = root->end;
		bud->covered = root->covered;
		*(bud->branch) = 0;
		
		if(node->end < root->start) {
			/* right bud */
			root->start = node->start;
			root->branch[0] = node;
			root->branch[1] = bud;
		} else {
			/* left bud */
			root->end = node->end;
			root->branch[0] = bud;
			root->branch[1] = node;
		}
		root->covered += node->covered;
	} else { /* extend leaf */
		if(node->start < root->start) { /* extend left */
			root->start = node->start;
		}
		if(root->end < node->end) { /* extend right */
			root->end = node->end;
		}
		node->covered = 0;
		root->covered = root->end - root->start;
	}
	
	
	return root->covered;
}

int growSeqmentTree(SeqmentTree *src, const unsigned start, const unsigned end) {
	
	SeqmentTrees *node;
	
	/* make room for new anker */
	if(src->size <= src->n + 2) {
		resizeSeqmentTree(src);
	} else if(src->n == 0) {
		initSeqmentTree(src, start, end);
		return end - start;
	}
	
	/* make new leaf */
	node = src->root + src->n;
	node->start = start;
	node->end = end;
	node->covered = end - start;
	*(node->branch) = 0;
	
	src->root->covered = addSeqmentTrees(src->root, node);
	
	if(node->covered) {
		src->n += 2;
	}
	
	return src->root->covered;
}

unsigned queSeqmentTree(SeqmentTrees *src, const unsigned start, const unsigned end) {
	
	if(end < src->start || src->end < start) {
		/* miss */
		return 0;
	} else if(start <= src->start && src->end <= end) {
		/* leaf covered by query */
		return src->covered;
	} else if(*(src->branch)) {
		/* check next */
		return queSeqmentTree(src->branch[0], start, end) + queSeqmentTree(src->branch[1], start, end);
	} else if(src->start <= start && end <= src->end) {
		return end - start;
	} else if(src->start <= start && start < src->end) {
		/* left side */
		return src->end - start;
	} else if(src->start < end && end <= src->end) {
		/* right side */
		return end - src->start;
	}
	
	return 0;
}
