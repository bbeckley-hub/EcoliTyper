/* fuzzy k-mer */
/*

T NNNNNNNN
Q NNNN*NNN

T NNNNNNNN -> NNNNNNNN
Q NNNN-NNN -> NNNN***-

T NNNN-NNNN -> NNNNNNNN-
Q NNNN*NNN- -> NNNN****-

Could we do a 1-shift when both halfs miss.


*/

/* adapter trim */
/*
SE, both ends (unless sr)
PE, 3' and only if overlapping.


Build directed de bruijn graph from the end of the sequence where the 
adaptors are expected, from the reads in the first buffer.
Trim the ends of the graph for low frequency k-mers, to get rid of non-adaptor 
sequences. There shouled be a frequency peak somewhere in the "middle", as 
the start of the adaptors are more frequent. I.e. do not trim nodes that are 
not pointed to.
	a. If known adaptors are given, mark terminating nodes and skip low-freq trimming.

Check for adaptors by searching from the end and towards the center of the read.
	a. Mismatches are handles by substituting the next nuc by the graph nuc.
	b. Insertion in seq, omit next nuc in k-mer
	c. Deletion in seq, add nuc from next in graph

Errors in fist n-k seeds can be handled as in BLAST.

For this we need:
	k: k-mer size.
	len: len of area in ends to build graph over.
	size: buffer size.
	edit: max edit dist.
	freq: low frequency threshold, to avoid removing already remoived adaptors.
	(fsa): known adaptor sequences.
	(sensitive): Add BLAST error correction.
	(sensitive): Divide k-mers in two, and space k-mers. So that each k-mer has four representations, allowing at least one error.

*/

/* primer trim */
/*
SE, both ends have primers
PE, start have primer, end if they overlap

Direction is known, but is rc at the 3'.

Must be given in advance, where graph can be terminated on ends (as with 
adaptors).
Trimming them is the same as for adaptors, but with rc too.

*/

/* merge PE */
/*
need:
	k: k-mer size
	len: minimum overlap.
	edit: max edit dist, expect low distance.

Sample k-mers (only unique?) on forward read according to the pigeon hole principle.
Seed-extend in reverse read, and check with max edit distance.

Adaptor trim if insert is less than merged size.

*/

#include <stdint.h>

int fuzzykmercmp(uint64_t k1, uint64_t k2) {
	
	/* fuzzy comparison of k-mers, return 0 for no mismatches, 1 for 1 
	mismatch and 2 for for >=2 mismatches */
	
	if(k1 == k2) { /* zero mismatches */
		return 0;
	}
	
	/* check each half the k-mers recursively */
	/*
	uint32_t shift = 32;
	uint32_t mask = UINT32_MAX;
	while(2 < shift) {
		if((k1 >> shift) == (k2 >> shift)) { // first half match
			k1 &= mask;
			k2 &= mask;
		} else if((k1 & mask) == (k2 & mask)) { // second half match
			k1 >>= shift;
			k2 >>= shift;
		} else { // at least two errors
			return 2;
		}
		shift >>= 1;
		mask >>= shift;
	}
	
	// cmp last two nucleotides
	return ((k1 >> shift) != (k2 >> shift)) + ((k1 & mask) != (k2 & mask));
	*/
	
	/* unroll */
	/* 16-mer */
	if((k1 & 4294967295) == (k2 & 4294967295)) { /* second half match */
		k1 >>= 32;
		k2 >>= 32;
	} else if((k1 >> 32) != (k2 >> 32)) { /* more than one error */
		return 2;
	}
	/* 8-mer */
	if((k1 & 65535) == (k2 & 65535)) { /* second half match */
		k1 >>= 16;
		k2 >>= 16;
	} else if((k1 >> 16) != (k2 >> 16)) { /* more than one error */
		return 2;
	}
	/* 4-mer */
	if((k1 & 255) == (k2 & 255)) { /* second half match */
		k1 >>= 8;
		k2 >>= 8;
	} else if((k1 >> 8) != (k2 >> 8)) { /* more than one error */
		return 2;
	}
	/* 2-mer */
	if((k1 & 15) == (k2 & 15)) { /* second half match */
		return ((k1 >> 6) != (k2 >> 6)) + ((k1 & 48) != (k2 & 48));
	} else if((k1 >> 4) != (k2 >> 4)) { /* more than one error */
		return 2;
	}
	
	/* last two nucleotides */
	return ((k1 >> 2) != (k2 >> 2)) + ((k1 & 3) != (k2 & 3));
}

int fuzzyindelcmp(uint64_t k1, uint64_t k2, uint32_t shift, uint64_t mask) {
	
	/* fuzzy comparison of k-mer indels, return 0 for no mismatches/indels, 
	1 for 1 indel and 2 for for >=2 indels */
	
	if(k1 == k2) { /* zero mismatches */
		return 0;
	}
	
	/* check each half the k-mers recursively */
	int indel = 0; /* shifted left or right */
	while(1 < shift) {
		if(indel <= 0 && (k1 >> 2) == (k2 & (mask | (mask << (shift - 2))))) { // test right, full match
			return 1;
		} else if(0 <= indel && (k2 >> 2) == (k1 & (mask | (mask << (shift - 2))))) { // test left, full match
			return 1;
		} else if(indel <= 0 && ((k1 >> 2) & mask) == (k2 & mask)) { // test right, last half match
			++indel;
			k1 >>= shift + 2;
			k2 >>= shift;
		} else if(indel <= 0 && ((k1 >> (shift - 2)) & mask) == (k2 >> shift)) { // test right, first half match
			++indel;
			k1 &= mask >> 2;
			k2 &= mask;
		} else if(0 <= indel && ((k2 >> 2) & mask) == (k1 & mask)) { // test left, last half match
			--indel;
			k2 >>= shift + 2;
			k1 >>= shift;
		} else if(0 <= indel && ((k2 >> (shift - 2)) & mask) == (k1 >> shift)) { // test left, first half match
			--indel;
			k2 &= mask >> 2;
			k1 &= mask;
		} else {
			return 2;
		}
		shift >>= 1;
		mask >>= shift;
	}
	
	return 1;
}

int fuzzierkmercmp(uint64_t k1, uint64_t k2) {
	
	/* fuzzy comparison of k-mers, return 0 for no mismatches/indels, 1 for 1 
	mismatch/indel and 2 for for >=2 mismatches/indels */
	
	if(k1 == k2) { /* zero mismatches */
		return 0;
	}
	
	/* check each half the k-mers recursively */
	/*
	uint32_t shift = 32;
	uint32_t mask = UINT32_MAX;
	while(1 < shift) {
		if((k1 >> shift) == (k2 >> shift)) { // first half match
			k1 &= mask;
			k2 &= mask;
		} else if((k1 & mask) == (k2 & mask)) { // second half match
			k1 >>= shift;
			k2 >>= shift;
		} else { // check indel
			return fuzzyindelcmp(k1, k2, shift, mask);
		}
		shift >>= 1;
		mask >>= shift;
	}
	
	// cmp last two nucleotides
	return 1;
	*/
	
	/* unroll */
	/* 16-mer */
	if((k1 & 4294967295) == (k2 & 4294967295)) { /* second half match */
		k1 >>= 32;
		k2 >>= 32;
	} else if((k1 >> 32) != (k2 >> 32)) { /* more than one mismatch */
		return fuzzyindelcmp(k1, k2, 32, 4294967295);
	}
	/* 8-mer */
	if((k1 & 65535) == (k2 & 65535)) { /* second half match */
		k1 >>= 16;
		k2 >>= 16;
	} else if((k1 >> 16) != (k2 >> 16)) { /* more than one mismatch */
		return fuzzyindelcmp(k1 & 4294967295, k2 & 4294967295, 16, 65535);
	}
	/* 4-mer */
	if((k1 & 255) == (k2 & 255)) { /* second half match */
		k1 >>= 8;
		k2 >>= 8;
	} else if((k1 >> 8) != (k2 >> 8)) { /* more than one mismatch */
		return fuzzyindelcmp(k1 & 65535, k2 & 65535, 8, 255);
	}
	/* 2-mer */
	if((k1 & 15) == (k2 & 15)) { /* second half match */
		return ((k1 >> 6) != (k2 >> 6)) + ((k1 & 48) != (k2 & 48));
	} else if((k1 >> 4) != (k2 >> 4)) { /* more than one mismatch */
		return fuzzyindelcmp(k1 & 255, k2 & 255, 4, 15);
	}
	
	/* last two nucleotides */
	return (((k1 >> 2) == (k2 >> 2)) || ((k1 & 3) == (k2 & 3)) || ((k1 >> 2) == (k2 & 3)) || ((k1 & 3) == (k2 >> 2))) ? 1 : 2;
}
