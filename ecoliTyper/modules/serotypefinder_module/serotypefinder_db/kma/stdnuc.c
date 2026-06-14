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

#include "pherror.h"
#include "stdnuc.h"

long unsigned (*initCmer)(long unsigned, int*, long unsigned*, int*, const unsigned, const int, const int, const long unsigned) = &initHmer;
long unsigned (*initCmerR)(long unsigned, int*, long unsigned*, int*, const unsigned, const int, const int, const long unsigned) = &initHmer;
long unsigned (*updateCmer)(long unsigned, int*, long unsigned*, int*, long unsigned, const int, const int, const long unsigned) = &updateKmerHom;
long unsigned (*updateCmerR)(long unsigned, int*, long unsigned*, int*, long unsigned, const int, const int, const long unsigned) = &updateKmerHomR;
long unsigned (*getCmer)(long unsigned, int*, int*, const unsigned, const int, const long unsigned) = &getHmer;
long unsigned (*getCmerR)(long unsigned, int*, int*, const unsigned, const int, const long unsigned) = &getHmer;

long unsigned updateKmer(const long unsigned pKmer, const long unsigned nuc, const long unsigned mask) {
	return (((pKmer << 2) | nuc) & mask);
}

long unsigned updateKmerHom(long unsigned pHmer, int *Pos, long unsigned *PHmer, int *H_len, long unsigned nKmer, const int kmersize, const int mlen, const long unsigned mmask) {
	
	int hLen;
	
	/* init */
	pHmer = *PHmer;
	
	/* update trailing base */
	if((pHmer & 3) != (nKmer & 3)) {
		/*
		pHmer <<= 2;
		pHmer |= (kmer & 3);
		*/
		pHmer = (pHmer << 2) | (nKmer & 3);
		++*H_len;
	}
	
	/* update leading base */
	hLen = *H_len << 1;
	nKmer >>= ((kmersize << 1) - 2);
	if((pHmer >> hLen) != nKmer) {
		/* first base is duplicated to acoid leading 0/A matches */
		/*
		mask = 0xFFFFFFFFFFFFFFFF >> (64 - (kmersize << 1)); normal mask
		mask = 0xFFFFFFFFFFFFFFFF >> (64 - (hLen - 2)); remove first base (+ redundant one)
		mask = 0xFFFFFFFFFFFFFFFF >> (66 - hLen); remove first base (+ redundant one)
		pHmer &= mask;
		pHmer |= (fnuc << (hLen - 2));
		*/
		pHmer = (pHmer & (0xFFFFFFFFFFFFFFFF >> (66 - hLen))) | (nKmer << (hLen - 2));
		--*H_len;
	}
	*PHmer = pHmer;
	
	return (kmersize != 16) ? pHmer : (pHmer & 0xFFFFFFFF);
}

long unsigned updateKmerHomR(long unsigned pHmer, int *Pos, long unsigned *PHmer, int *H_len, long unsigned nKmer, const int kmersize, const int mlen, const long unsigned mmask) {
	
	int hLen;
	
	/* init */
	pHmer = *PHmer;
	
	/* update trailing base */
	if((pHmer & 3) != (nKmer & 3)) {
		pHmer >>= 2;
		--*H_len;
	}
	
	/* update leading base */
	hLen = *H_len << 1;
	nKmer >>= ((kmersize << 1) - 2);
	if((pHmer >> hLen) != nKmer) {
		/* first base is duplicated to acoid leading 0/A matches */
		/*
		mask = 0xFFFFFFFFFFFFFFFF >> (64 - hLen); remove redundant base
		pHmer &= mask;
		nKmer = (nKmer << 2) | nKmer;
		nKmer <<= hLen;
		pHmer |= nKmer;
		*/
		
		pHmer = (pHmer & (0xFFFFFFFFFFFFFFFF >> (64 - hLen))) | (((nKmer << 2) | nKmer) << hLen);
		++*H_len;
	}
	*PHmer = pHmer;
	
	return (kmersize != 16) ? pHmer : (pHmer & 0xFFFFFFFF);
}

long unsigned updateKmerMin(long unsigned pKmer, int *Pos, long unsigned *PHmer, int *H_len, long unsigned nKmer, const int kmersize, const int mlen, const long unsigned mmask) {
	
	int i;
	/*
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	*/
	
	/* check new */
	if((nKmer & mmask) <= pKmer) {
		pKmer = (nKmer & mmask);
		*Pos = kmersize - mlen + 1;
	} else if(!--*Pos) { /* find new */
		pKmer = (nKmer & mmask);
		nKmer >>= 2;
		*Pos = (i = (kmersize - mlen + 1));
		while(--i) {
			if((nKmer & mmask) <= pKmer) {
				pKmer = (nKmer & mmask);
				*Pos = i;
			}
			nKmer >>= 2;
		}
	}
	
	return pKmer;
}

long unsigned updateKmerMinR(long unsigned pKmer, int *Pos, long unsigned *PHmer, int *H_len, long unsigned nKmer, const int kmersize, const int mlen, const long unsigned mmask) {
	
	int i;
	/*
	mmask = 0xFFFFFFFFFFFFFFFF >> (64 - (mlen << 1));
	*/
	
	/* check new */
	i = (kmersize - mlen) << 1;
	if((nKmer >> i) <= pKmer) {
		pKmer = nKmer >> i;
		*Pos = (i >> 1) + 1;
	} else if(!--*Pos) {
		/* get minimizer */
		pKmer = nKmer >> i;
		*Pos = (i >> 1) + 1;
		while(i) {
			if(((nKmer >> (i -= 2)) & mmask) < pKmer) {
				pKmer = ((nKmer >> i) & mmask);
				*Pos = (i >> 1) + 1;
			}
		}
	}
	
	return pKmer;
}

long unsigned updateKmerHomMin(long unsigned pKmer, int *Pos, long unsigned *PHmer, int *H_len, long unsigned nKmer, const int kmersize, const int mlen, const long unsigned mmask) {
	
	int hLen;
	long unsigned pHmer;
	
	/* homopolymer compression */
	pHmer = *PHmer;
	
	/* update trailing base */
	if((pHmer & 3) != (nKmer & 3)) {
		/*
		pHmer <<= 2;
		pHmer |= (kmer & 3);
		*/
		pHmer = (pHmer << 2) | (nKmer & 3);
		++*H_len;
	}
	
	/* update leading base */
	hLen = *H_len << 1;
	nKmer >>= ((kmersize << 1) - 2);
	if((pHmer >> hLen) != nKmer) {
		/* first base is duplicated to acoid leading 0/A matches */
		pHmer = (pHmer & (0xFFFFFFFFFFFFFFFF >> (66 - hLen))) | (nKmer << (hLen - 2));
		--*H_len;
	} else {
		/* maintain minimizer existence */
		++*Pos;
	}
	*PHmer = pHmer;
	
	/* minimizer update */
	if(*H_len <= mlen) {
		pHmer &= mmask;
		*Pos = 1;
	} else {
		pHmer = updateKmerMin(pKmer, Pos, PHmer, H_len, pHmer, *H_len, mlen, mmask);
	}
	
	return pHmer;
}

long unsigned updateKmerHomMinR(long unsigned pKmer, int *Pos, long unsigned *PHmer, int *H_len, long unsigned nKmer, const int kmersize, const int mlen, const long unsigned mmask) {
	
	int hLen;
	long unsigned pHmer, mask;
	
	/* homopolymer compression */
	pHmer = *PHmer;
	
	/* update trailing base */
	if((pHmer & 3) != (nKmer & 3)) {
		pHmer >>= 2;
		--*H_len;
	} else {
		/* maintain minimizer existence */
		++*Pos;
	}
	
	/* update leading base */
	hLen = *H_len << 1;
	nKmer >>= ((kmersize << 1) - 2);
	mask = (0xFFFFFFFFFFFFFFFF >> (64 - hLen));
	if((pHmer >> hLen) != nKmer) {
		/* first base is duplicated to acoid leading 0/A matches */
		pHmer = (pHmer & mask) | (((nKmer << 2) | nKmer) << hLen);
		++*H_len;
		mask = ((mask << 2) | 3);
	}
	*PHmer = pHmer;
	
	/* minimizer update */
	if(*H_len <= mlen) {
		pHmer &= mmask;
		*Pos = 1;
	} else {
		pHmer = updateKmerMinR(pKmer, Pos, PHmer, H_len, (pHmer & mask), *H_len, mlen, mmask);
	}
	
	return pHmer;
}

long unsigned initHmer(long unsigned kmer, int *Pos, long unsigned *hukmer, int *H_len, const unsigned shifter, const int kmersize, const int mlen, const long unsigned mmask) {
	
	/* init homopolymer compression for updates */
	*Pos = 1;
	*hukmer = getHmer(kmer, Pos, H_len, shifter, mlen, mmask);
	*Pos = kmersize; // Mitigate unsigned overflow
	
	return *hukmer;
}

long unsigned initMmer(long unsigned kmer, int *Pos, long unsigned *hukmer, int *H_len, const unsigned shifter, const int kmersize, const int mlen, const long unsigned mmask) {
	
	/* init minimizer compression for updates */
	*Pos = 1;
	*hukmer = 0;
	*H_len = kmersize;
	
	return 0;
}

long unsigned initHMmer(long unsigned kmer, int *Pos, long unsigned *hukmer, int *H_len, const unsigned shifter, const int kmersize, const int mlen, const long unsigned mmask) {
	
	/* init homopolymer minimizer compression for updates */
	*Pos = 1;
	*hukmer = getHmer(kmer, Pos, H_len, shifter, mlen, mmask);
	
	return getMmer(*hukmer, Pos, H_len, 0, mlen, mmask);
}

long unsigned initHMmerR(long unsigned kmer, int *Pos, long unsigned *hukmer, int *H_len, const unsigned shifter, const int kmersize, const int mlen, const long unsigned mmask) {
	
	/* init homopolymer minimizer compression for rc updates */
	*Pos = 1;
	*hukmer = getHmer(kmer, Pos, H_len, shifter, mlen, mmask);
	
	return getMmerR(*hukmer, Pos, H_len, 0, mlen, mmask);;
}

long unsigned getKmer(long unsigned *compressor, unsigned cPos, const unsigned shifter) {
	
	/* see getKmer_macro(kmer, Comp, cPos, iPos, shifter) */
	unsigned iPos = (cPos & 31) << 1;
	cPos >>= 5;
	
	return (iPos <= shifter) ? ((compressor[cPos] << iPos) >> shifter) : (((compressor[cPos] << iPos) | (compressor[cPos + 1] >> (64-iPos))) >> shifter);
}

long unsigned getHmer(long unsigned kmer, int *Pos, int *H_len, const unsigned shifter, const int mlen, const long unsigned mmask) {
	
	int i, hLen, nuc;
	long unsigned hmer;
	
	/* transform to homopolymer compressed kmer */
	/*
	shifter = 64 - (kmersize << 1);
	iPos = (kmersize << 1) = 64 - (64 - (kmersize << 1)) = 64 - shifter;
	*/
	i = 62 - shifter;
	hmer = kmer >> i;
	hmer = (hmer << 2) | hmer;
	hLen = 1;
	while(0 <= (i -= 2)) {
		if((nuc = ((kmer >> i) & 3)) != (hmer & 3)) {
			hmer = (hmer << 2) | nuc;
			++hLen;
		}
	}
	
	*H_len = hLen;
	return (*Pos != 16) ? hmer : (hmer & 0xFFFFFFFF);
}

long unsigned getMmer(long unsigned kmer, int *Pos, int *kmersize, const unsigned shifter, const int mlen, const long unsigned mmask) {
	
	int i;
	long unsigned mmer;
	
	/* get minimizer */
	if((i = (*kmersize - mlen + 1)) <= 0) {
		*Pos = 1;
		return (kmer & mmask);
	}
	mmer = (kmer & mmask);
	kmer >>= 2;
	*Pos = i;
	while(--i) {
		if((kmer & mmask) < mmer) {
			mmer = (kmer & mmask);
			*Pos = i;
		}
		kmer >>= 2;
	}
	
	return mmer;
}

long unsigned getMmerR(long unsigned kmer, int *Pos, int *kmersize, const unsigned shifter, const int mlen, const long unsigned mmask) {
	
	int i;
	long unsigned mmer;
	
	/* get minimizer */
	if((i = ((*kmersize - mlen) << 1)) <= 0) {
		*Pos = 1;
		return (kmer & mmask);
	}
	mmer = (kmer >> i) & mmask;
	*Pos = (i >> 1) + 1;
	while(i) {
		if(((kmer >> (i -= 2)) & mmask) < mmer) {
			mmer = ((kmer >> i) & mmask);
			*Pos = (i >> 1) + 1;
		}
	}
	
	return mmer;
}

long unsigned getHMmer(long unsigned kmer, int *Pos, int *H_len, const unsigned shifter, const int mlen, const long unsigned mmask) {
	
	int i, hLen;
	long unsigned hmer, mmer;
	
	/* transform to homopolymer compressed kmer */
	i = 62 - shifter;
	hmer = kmer >> i;
	hmer = (hmer << 2) | hmer;
	hLen = 1;
	while(0 <= (i -= 2)) {
		if((mmer = ((kmer >> i) & 3)) != (hmer & 3)) {
			hmer = (hmer << 2) | mmer;
			++hLen;
		}
	}
	*H_len = hLen;
	
	/* get minimizer */
	if(mlen < hLen) {
		mmer = (hmer & mmask);
		hmer >>= 2;
		i = (hLen - mlen + 1);
		*Pos = i;
		while(--i) {
			if((hmer & mmask) < mmer) {
				mmer = (hmer & mmask);
				*Pos = i;
			}
			hmer >>= 2;
		}
	} else {
		mmer = (hmer & mmask);
		*Pos = 1;
	}
	return mmer;
}

void setCmerPointers(unsigned flag) {
	
	if((flag &= 3)) {
		if(flag == 1) {
			initCmer = &initHmer;
			initCmerR = &initHmer;
			updateCmer = &updateKmerHom;
			updateCmerR = &updateKmerHomR;
			getCmer = &getHmer;
			getCmerR = &getHmer;
		} else if(flag == 2) {
			initCmer = &initMmer;
			initCmerR = &initMmer;
			updateCmer = &updateKmerMin;
			updateCmerR = &updateKmerMinR;
			getCmer = &getMmer;
			getCmerR = &getMmerR;
		} else if(flag == 3) {
			initCmer = &initHMmer;
			initCmerR = &initHMmerR;
			updateCmer = &updateKmerHomMin;
			updateCmerR = &updateKmerHomMinR;
			getCmer = &getHMmer;
			getCmerR = &getHMmer;
		}
	}
}

long unsigned makeKmer(const unsigned char *qseq, unsigned pos, unsigned size) {
	
	long unsigned key = qseq[pos];
	
	size += pos;
	for(++pos; pos < size; ++pos) {
		key = (key << 2) | qseq[pos];
	}
	
	return key;
}

int charpos(const unsigned char *src, unsigned char target, int start, int len) {
	
	unsigned char *ptr;
	
	ptr = (unsigned char *) src + --start;
	while(++start < len) {
		if(*++ptr == target) {
			return start;
		}
	}
	
	return -1;
}

void strrc(unsigned char *qseq, int q_len) {
	
	int i, j, seqlen;
	unsigned char carry, comp[6] = {3, 2, 1, 0, 4, 5};
	
	seqlen = q_len >> 1;
	
	for(i = 0, j = q_len - 1; i < seqlen; ++i, --j) {
		carry = comp[qseq[i]];
		qseq[i] = comp[qseq[j]];
		qseq[j] = carry;
	}
	if(q_len & 1) {
		qseq[seqlen] = comp[qseq[seqlen]];
	}
	
}

void strtranslate(unsigned char *strp, char *trans) {
	--strp;
	while(*++strp) {
		*strp = trans[*strp];
	}
}

void nibble2base(unsigned char *seq, int len) {
	
	const char bases[6] = "ACGTN-";
	
	seq += len;
	*seq-- = 0;
	++len;
	while(--len) {
		*seq = bases[*seq];
		--seq;
	}	
}

unsigned char * base2nibble_table() {
	
	int i;
	unsigned char *table;
	
	table = smalloc(256);
	for(i = 0; i < 256; ++i) {
		table[i] = 5;
	}
	table['A'] = 0;
	table['C'] = 1;
	table['G'] = 2;
	table['T'] = 3;
	table['N'] = 4;
	table['a'] = 0;
	table['c'] = 1;
	table['g'] = 2;
	table['t'] = 3;
	table['n'] = 4;
	table['R'] = 0;
	table['Y'] = 1;
	table['S'] = 2;
	table['W'] = 3;
	table['K'] = 2;
	table['M'] = 0;
	table['B'] = 1;
	table['D'] = 0;
	table['H'] = 3;
	table['V'] = 2;
	table['X'] = 4;
	table['r'] = 0;
	table['y'] = 1;
	table['s'] = 2;
	table['w'] = 3;
	table['k'] = 2;
	table['m'] = 0;
	table['b'] = 1;
	table['d'] = 0;
	table['h'] = 3;
	table['v'] = 2;
	table['x'] = 4;
	table['U'] = 3;
	table['u'] = 3;
	table['-'] = 5;
	
	return table;
}

unsigned char * complement_table() {
	
	int i;
	unsigned char *table;
	
	table = smalloc(256);
	for(i = 0; i < 256; ++i) {
		table[i] = i;
	}
	table['A'] = 'T';
	table['C'] = 'G';
	table['G'] = 'C';
	table['T'] = 'A';
	table['a'] = 't';
	table['c'] = 'g';
	table['g'] = 'c';
	table['t'] = 'a';
	
	table['R'] = 'Y';
	table['Y'] = 'R';
	table['S'] = 'W';
	table['W'] = 'S';
	table['K'] = 'M';
	table['M'] = 'K';
	table['B'] = 'D';
	table['D'] = 'B';
	table['H'] = 'V';
	table['V'] = 'H';
	table['U'] = 'A';
	
	table['r'] = 'y';
	table['y'] = 'r';
	table['s'] = 'w';
	table['w'] = 's';
	table['k'] = 'm';
	table['m'] = 'k';
	table['b'] = 'd';
	table['d'] = 'b';
	table['h'] = 'v';
	table['v'] = 'h';
	table['u'] = 'a';
	
	return table;
}
