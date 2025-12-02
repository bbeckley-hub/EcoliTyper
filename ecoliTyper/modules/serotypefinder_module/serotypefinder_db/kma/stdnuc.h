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

#define getNuc(Comp,pos) ((Comp[pos >> 5] << ((pos & 31) << 1)) >> 62)
#define updateKmer_macro(kmer, Comp, pos, mask) ((((kmer << 2) | (getNuc(Comp,pos))) & mask))
//rShifter = 64 - (kmersize << 1)
#define updateKmerR_macro(kmer, Comp, pos, rShifter) ((kmer >> 2) | ((Comp[pos >> 5] << ((pos & 31) << 1)) >> rShifter))
#define setEx(src, pos)(src[pos >> 3] |= (1 << (pos & 7)))
#define unsetEx(src, pos)(src[pos >> 3] ^= (1 << (pos & 7)))
#define getEx(src, pos)((src[pos >> 3] >> (pos & 7)) & 1)
#define getKmer_macro(kmer, Comp, pos, cPos, iPos, shifter) \
		iPos = (pos & 31) << 1;\
		cPos = pos >> 5;\
		kmer = (iPos <= shifter) ? ((Comp[cPos] << iPos) >> shifter) : (((Comp[cPos] << iPos) | (Comp[cPos + 1] >> (64-iPos))) >> shifter);

extern long unsigned (*initCmer)(long unsigned, int*, long unsigned*, int*, const unsigned, const int, const int, const long unsigned);
extern long unsigned (*initCmerR)(long unsigned, int*, long unsigned*, int*, const unsigned, const int, const int, const long unsigned);
extern long unsigned (*updateCmer)(long unsigned, int*, long unsigned*, int*, long unsigned, const int, const int, const long unsigned);
extern long unsigned (*updateCmerR)(long unsigned, int*, long unsigned*, int*, long unsigned, const int, const int, const long unsigned);
extern long unsigned (*getCmer)(long unsigned, int*, int*, const unsigned, const int, const long unsigned);
extern long unsigned (*getCmerR)(long unsigned, int*, int*, const unsigned, const int, const long unsigned);
long unsigned updateKmer(const long unsigned pKmer, const long unsigned nuc, const long unsigned mask);
long unsigned updateKmerHom(long unsigned pHmer, int *Pos, long unsigned *PHmer, int *H_len, long unsigned nKmer, const int kmersize, const int mlen, const long unsigned mmask);
long unsigned updateKmerHomR(long unsigned pHmer, int *Pos, long unsigned *PHmer, int *H_len, long unsigned nKmer, const int kmersize, const int mlen, const long unsigned mmask);
long unsigned updateKmerMin(long unsigned pKmer, int *Pos, long unsigned *PHmer, int *H_len, long unsigned nKmer, const int kmersize, const int mlen, const long unsigned mmask);
long unsigned updateKmerMinR(long unsigned pKmer, int *Pos, long unsigned *PHmer, int *H_len, long unsigned nKmer, const int kmersize, const int mlen, const long unsigned mmask);
long unsigned updateKmerHomMin(long unsigned pKmer, int *Pos, long unsigned *PHmer, int *H_len, long unsigned nKmer, const int kmersize, const int mlen, const long unsigned mmask);
long unsigned updateKmerHomMinR(long unsigned pKmer, int *Pos, long unsigned *PHmer, int *H_len, long unsigned nKmer, const int kmersize, const int mlen, const long unsigned mmask);
long unsigned initHmer(long unsigned kmer, int *Pos, long unsigned *hukmer, int *H_len, const unsigned shifter, const int kmersize, const int mlen, const long unsigned mmask);
long unsigned initMmer(long unsigned kmer, int *Pos, long unsigned *hukmer, int *H_len, const unsigned shifter, const int kmersize, const int mlen, const long unsigned mmask);
long unsigned initHMmer(long unsigned kmer, int *Pos, long unsigned *hukmer, int *H_len, const unsigned shifter, const int kmersize, const int mlen, const long unsigned mmask);
long unsigned initKmerR(long unsigned kmer, int *Pos, long unsigned *hukmer, int *H_len, const unsigned shifter, const int kmersize, const int mlen, const long unsigned mmask);
long unsigned initHMmerR(long unsigned kmer, int *Pos, long unsigned *hukmer, int *H_len, const unsigned shifter, const int kmersize, const int mlen, const long unsigned mmask);
long unsigned getKmer(long unsigned *compressor, unsigned cPos, const unsigned shifter);
long unsigned getHmer(long unsigned kmer, int *Pos, int *H_len, const unsigned shifter, const int mlen, const long unsigned mmask);
long unsigned getMmer(long unsigned kmer, int *Pos, int *kmersize, const unsigned shifter, const int mlen, const long unsigned mmask);
long unsigned getMmerR(long unsigned kmer, int *Pos, int *kmersize, const unsigned shifter, const int mlen, const long unsigned mmask);
long unsigned getHMmer(long unsigned kmer, int *Pos, int *H_len, const unsigned shifter, const int mlen, const long unsigned mmask);
void setCmerPointers(unsigned flag);
long unsigned makeKmer(const unsigned char *qseq, unsigned pos, unsigned size);
int charpos(const unsigned char *src, unsigned char target, int start, int len);
void strrc(unsigned char *qseq, int q_len);
void strtranslate(unsigned char *strp, char *trans);
void nibble2base(unsigned char *seq, int len);
unsigned char * base2nibble_table();
unsigned char * complement_table();
