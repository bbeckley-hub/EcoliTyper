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
#include <string.h>
#include "nw.h"
#include "penalties.h"
#include "pherror.h"
#include "stdnuc.h"

AlnScore NW2P(const long unsigned *template, const unsigned char *queryOrg, int k, int t_s, int t_e, int q_s, int q_e, Aln *aligned, NWmat *matrices, int template_length) {
	
	int m, n, t_len, q_len, thisScore, nuc_pos, W1, U, W2, U2, MM;
	int pos[2], *tmp, **d, *D_ptr, *D_prev, Q, Q_prev, *P_ptr, *P_prev;
	int Q2, Q2_prev, *P2_ptr, *P2_prev;
	unsigned char *query, t_nuc, *E, *E_ptr, e;
	AlnScore Stat;
	Penalties *rewards;
	
	rewards = matrices->rewards;
	W1 = rewards->W1;
	U = rewards->U;
	W2 = rewards->W2;
	U2 = rewards->U2;
	MM = rewards->MM;
	d = rewards->d;
	q_len = q_e - q_s;
	t_len = t_e - t_s;
	aligned->start = 0;
	aligned->end = 0;
	if(t_len < 0) {
		t_len += template_length;
	}
	query = (unsigned char*)(queryOrg + q_s);
	Stat.pos = 0;
	
	if(t_len == 0 || q_len == 0) {
		if(t_len == q_len) {
			Stat.len = 0;
			Stat.match = 0;
			Stat.tGaps = 0;
			Stat.qGaps = 0;
			Stat.score = 0;
			aligned->s[0] = 0;
		} else if(t_len == 0) {
			Stat.len = q_len;
			Stat.match = 0;
			Stat.tGaps = q_len;
			Stat.qGaps = 0;
			Stat.score = W1 + (q_len - 1) * U;
			if((W2 + (q_len - 1) * U2) < Stat.score) {
				Stat.score = (W2 + (q_len - 1) * U2);
			}
			memset(aligned->s, '_', q_len);
			aligned->s[q_len] = 0;
			memset(aligned->t, 5, q_len);
			memcpy(aligned->q, query, q_len);
		} else {
			Stat.len = t_len;
			Stat.match = 0;
			Stat.tGaps = 0;
			Stat.qGaps = t_len;
			Stat.score = W1 + (t_len - 1) * U;
			if((W2 + (t_len - 1) * U2) < Stat.score) {
				Stat.score = (W2 + (t_len - 1) * U2);
			}
			memset(aligned->s, '_', t_len);
			aligned->s[t_len] = 0;
			memset(aligned->q, 5, t_len);
			m = t_len;
			nuc_pos = (t_e ? t_e : aligned->pos) - 1;
			while(m--) {
				aligned->t[m] = getNuc(template, nuc_pos);
				if(--nuc_pos < 0) {
					nuc_pos = aligned->pos - 1;
				}
			}
		}
		return Stat;
	}
	
	/* check matrix size */
	if(matrices->NW_q <= q_len) {
		matrices->NW_q = q_len << 1;
		free(matrices->D[0]);
		free(matrices->P[0]);
		free(matrices->P2[0]);
		matrices->D[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->P[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->P2[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->D[1] = matrices->D[0] + matrices->NW_q;
		matrices->P[1] = matrices->P[0] + matrices->NW_q;
		matrices->P2[1] = matrices->P[0] + matrices->NW_q;
	}
	if(matrices->NW_s <= ((q_len + 1) * (t_len + 1))) {
		matrices->NW_s = ((q_len + 2) * (t_len + 2));
		free(matrices->E);
		matrices->E = smalloc(matrices->NW_s);
	}
	
	/* fill in start penalties */
	D_ptr = matrices->D[0];
	D_prev = matrices->D[1];
	P_ptr = matrices->P[0];
	P_prev = matrices->P[1];
	P2_ptr = matrices->P2[0];
	P2_prev = matrices->P2[1];
	E = matrices->E;
	thisScore = (t_len + q_len) * (MM + U + W1);
	Stat.score = thisScore;
	if(0 < k) {
		E_ptr = E;
		for(m = 0; m < t_len; ++m) {
			E_ptr[q_len] = 0;
			E_ptr += (q_len + 1);
		}
		if(k == 1) {
			for(n = q_len - 1; n >= 0; --n) {
				Q = W1 + (q_len - 1 - n) * U;
				Q2 = W2 + (q_len - 1 - n) * U2;
				if(Q2 < Q) {
					D_prev[n] = W1 + (q_len - 1 - n) * U;
					E_ptr[n] = 2;
				} else {
					D_prev[n] = W2 + (q_len - 1 - n) * U2;
					E_ptr[n] = 4;
				}
				P_prev[n] = thisScore;
				P2_prev[n] = thisScore;
			}
			E_ptr[q_len - 1] = 86;
			E_ptr[q_len] = 0;
			D_prev[q_len] = 0;
			P_prev[q_len] = 0;
			P2_prev[q_len] = 0;
		} else {
			for(n = q_len; n >= 0; --n) {
				D_prev[n] = 0;
				P_prev[n] = thisScore;
				P2_prev[n] = thisScore;
				E_ptr[n] = 0;
			}
		}
	} else {
		E_ptr = E;
		for(m = 0; m < t_len; ++m) {
			Q = W1 + (t_len - 1 - m) * U;
			Q2 = W2 + (t_len - 1 - m) * U2;
			if(Q2 < Q) {
				E_ptr[q_len] = 3;
			} else {
				E_ptr[q_len] = 5;
			}
			E_ptr += (q_len + 1);
		}
		E_ptr[-1] = 168;
		
		for(n = q_len - 1; n >= 0; --n) {
			Q = W1 + (q_len - 1 - n) * U;
			Q2 = W2 + (q_len - 1 - n) * U2;
			if(Q2 < Q) {
				D_prev[n] = Q;
				E_ptr[n] = 2;
			} else {
				D_prev[n] = Q2;
				E_ptr[n] = 4;
			}
			P_prev[n] = thisScore;
			P2_prev[n] = thisScore;
		}
		E_ptr[q_len - 1] = 86;
		E_ptr[q_len] = 0;
		D_prev[q_len] = 0;
		P_prev[q_len] = 0;
		P2_prev[q_len] = 0;
	}
	E_ptr -= (q_len + 1);
	
	/*
	E  0
	D  1
	Q  2, 16
	P  3, 32
	Q2 4, 64
	P2 5, 128
	*/
	
	/* Perform NW */
	pos[0] = 0;
	for(m = t_len - 1, nuc_pos = t_e - 1; m >= 0; --m, --nuc_pos) {
		
		if(nuc_pos < 0) {
			nuc_pos = template_length - 1;
		}
		
		D_ptr[q_len] = (0 < k) ? 0 : (W1 + (t_len - 1 - m) * U);
		Q_prev = (t_len + q_len) * (MM + U + W1);
		Q2_prev = (t_len + q_len) * (MM + U2 + W2);
		t_nuc = getNuc(template, nuc_pos);
		for(n = q_len - 1; n >= 0; --n) {
			E_ptr[n] = 0;
			
			/* update Q and P, gap openings */
			Q = D_ptr[n + 1] + W1;
			P_ptr[n] = D_prev[n] + W1;
			if(Q < P_ptr[n]) {
				D_ptr[n] = P_ptr[n];
				e = 3;
			} else {
				D_ptr[n] = Q;
				e = 2;
			}
			
			/* update Q and P, gap extensions */
			/* mark bit 4 and 5 as possible gap-opennings, if necesarry */
			thisScore = Q_prev + U;
			if(Q < thisScore) {
				Q = thisScore;
				if(D_ptr[n] <= thisScore) {
					D_ptr[n] = thisScore;
					e = 2;
				}
			} else {
				E_ptr[n] |= 16;
			}
			thisScore = P_prev[n] + U;
			if(P_ptr[n] < thisScore) {
				P_ptr[n] = thisScore;
				if(D_ptr[n] <= thisScore) {
					D_ptr[n] = thisScore;
					e = 3;
				}
			} else {
				E_ptr[n] |= 32;
			}
			
			/* update Q2 and P2, gap openings */
			Q2 = D_ptr[n + 1] + W2;
			P2_ptr[n] = D_prev[n] + W2;
			if(Q2 < P2_ptr[n]) {
				D_ptr[n] = P2_ptr[n];
				e = 5;
			} else {
				D_ptr[n] = Q2;
				e = 4;
			}
			
			/* update Q2 and P2, gap extensions */
			/* mark bit 6 and 7 as possible gap-opennings, if necesarry */
			thisScore = Q2_prev + U2;
			if(Q2 < thisScore) {
				Q2 = thisScore;
				if(D_ptr[n] <= thisScore) {
					D_ptr[n] = thisScore;
					e = 4;
				}
			} else {
				E_ptr[n] |= 64;
			}
			thisScore = P2_prev[n] + U2;
			if(P2_ptr[n] < thisScore) {
				P2_ptr[n] = thisScore;
				if(D_ptr[n] <= thisScore) {
					D_ptr[n] = thisScore;
					e = 5;
				}
			} else {
				E_ptr[n] |= 128;
			}
			
			/* Update D, match */
			thisScore = D_prev[n + 1] + d[t_nuc][query[n]];
			if(D_ptr[n] <= thisScore) {
				D_ptr[n] = thisScore;
				E_ptr[n] |= 1;
			} else {
				E_ptr[n] |= e;
			}
			
			Q_prev = Q;
		}
		
		E_ptr -= (q_len + 1);
		
		if(k < 0 && Stat.score < *D_ptr) {
			Stat.score = *D_ptr;
			pos[0] = m;
		}
		
		tmp = D_ptr;
		D_ptr = D_prev;
		D_prev = tmp;
		
		tmp = P_ptr;
		P_ptr = P_prev;
		P_prev = tmp;
		
		tmp = P2_ptr;
		P2_ptr = P2_prev;
		P2_prev = tmp;
	}
	E_ptr = E;
	
	/* get start position of alignment */
	if(k < 0) {
		pos[1] = 0;
		if(k == -2) {
			for(n = 0; n < q_len; ++n) {
				if(Stat.score <= D_prev[n]) {
					Stat.score = D_prev[n];
					pos[0] = 0;
					pos[1] = (aligned->start = n);
				}
			}
		}
	} else {
		Stat.score = *D_prev;
		pos[0] = 0;
		pos[1] = 0;
	}
	
	/* make back tracking */
	m = pos[0];
	n = pos[1];
	E_ptr = E + (m * (q_len + 1));
	nuc_pos = m + t_s;
	Stat.len = 0;
	Stat.match = 0;
	Stat.tGaps = 0;
	Stat.qGaps = 0;
	while(E_ptr[n] != 0) {
		if(nuc_pos == template_length) {
			nuc_pos = 0;
		}
		e = E_ptr[n] & 15;
		if(e == 1) {
			aligned->t[Stat.len] = getNuc(template, nuc_pos);
			aligned->q[Stat.len] = query[n];
			aligned->s[Stat.len] = (aligned->t[Stat.len] == aligned->q[Stat.len]) ? '|' : '_';
			++Stat.match;
			++nuc_pos;
			E_ptr += (q_len + 1);
			++n;
		} else if(e == 2) {
			/* Q */
			while(!(E_ptr[n] & 16)) {
				aligned->t[Stat.len] = 5;
				aligned->q[Stat.len] = query[n];
				aligned->s[Stat.len] = '_';
				++n;
				++Stat.len;
				++Stat.tGaps;
			}
			aligned->t[Stat.len] = 5;
			aligned->q[Stat.len] = query[n];
			aligned->s[Stat.len] = '_';
			++n;
			++Stat.tGaps;
		} else if(e == 3) {
			/* P */
			while(!(E_ptr[n] & 32)) {
				aligned->t[Stat.len] = getNuc(template, nuc_pos);
				aligned->q[Stat.len] = 5;
				aligned->s[Stat.len] = '_';	
				++nuc_pos;
				E_ptr += (q_len + 1);
				++Stat.len;
				++Stat.qGaps;
			}
			aligned->t[Stat.len] = getNuc(template, nuc_pos);
			aligned->q[Stat.len] = 5;
			aligned->s[Stat.len] = '_';		
			++nuc_pos;
			E_ptr += (q_len + 1);
			++Stat.qGaps;
		} else if(e == 4) {
			/* Q2 */
			while(!(E_ptr[n] & 64)) {
				aligned->t[Stat.len] = 5;
				aligned->q[Stat.len] = query[n];
				aligned->s[Stat.len] = '_';
				++n;
				++Stat.len;
				++Stat.tGaps;
			}
			aligned->t[Stat.len] = 5;
			aligned->q[Stat.len] = query[n];
			aligned->s[Stat.len] = '_';
			++n;
			++Stat.tGaps;
		} else /*if(e == 5)*/ {
			/* P2 */
			while(!(E_ptr[n] & 128)) {
				aligned->t[Stat.len] = getNuc(template, nuc_pos);
				aligned->q[Stat.len] = 5;
				aligned->s[Stat.len] = '_';	
				++nuc_pos;
				E_ptr += (q_len + 1);
				++Stat.len;
				++Stat.qGaps;
			}
			aligned->t[Stat.len] = getNuc(template, nuc_pos);
			aligned->q[Stat.len] = 5;
			aligned->s[Stat.len] = '_';		
			++nuc_pos;
			E_ptr += (q_len + 1);
			++Stat.qGaps;
		}
		++Stat.len;
	}
	aligned->s[Stat.len] = 0;
	aligned->end = q_len - n;
	
	return Stat;
}

AlnScore NW_band(const long unsigned *template, const unsigned char *queryOrg, int k, int t_s, int t_e, int q_s, int q_e, Aln *aligned, int band, NWmat *matrices, int template_length) {
	
	int m, n, t_len, q_len, thisScore, nuc_pos, pos[2];
	int bq_len, halfBand, sn, en, sq, eq, q_pos, c_pos, W1, U, MM;
	int *D_ptr, *D_prev, Q, Q_prev, *P_ptr, *P_prev, *tmp, **d;
	unsigned char *query, t_nuc, *E, *E_ptr, e;
	AlnScore Stat;
	Penalties *rewards;
	
	rewards = matrices->rewards;
	W1 = rewards->W1;
	U = rewards->U;
	MM = rewards->MM;
	d = rewards->d;
	q_len = q_e - q_s;
	t_len = t_e - t_s;
	aligned->start = 0;
	aligned->end = 0;
	if(t_len < 0) {
		t_len += template_length;
	}
	query = (unsigned char*)(queryOrg + q_s);
	Stat.pos = 0;
	
	if(t_len == 0 || q_len == 0) {
		if(t_len == q_len) {
			Stat.score = 0;
			Stat.len = 0;
			Stat.match = 0;
			Stat.tGaps = 0;
			Stat.qGaps = 0;
			aligned->s[0] = 0;
		} else if(t_len == 0) {
			Stat.len = q_len;
			Stat.match = 0;
			Stat.tGaps = q_len;
			Stat.qGaps = 0;
			Stat.score = W1 + (q_len - 1) * U;
			memset(aligned->s, '_', q_len);
			aligned->s[q_len] = 0;
			memset(aligned->t, 5, q_len);
			memcpy(aligned->q, query, q_len);
		} else {
			Stat.len = t_len;
			Stat.match = 0;
			Stat.tGaps = 0;
			Stat.qGaps = t_len;
			Stat.score = W1 + (t_len - 1) * U;
			memset(aligned->s, '_', t_len);
			aligned->s[t_len] = 0;
			memset(aligned->q, 5, t_len);
			m = t_len;
			nuc_pos = (t_e ? t_e : aligned->pos) - 1;
			while(m--) {
				aligned->t[m] = getNuc(template, nuc_pos);
				if(--nuc_pos < 0) {
					nuc_pos = aligned->pos - 1;
				}
			}
		}
		return Stat;
	}
	
	/* ensure that band is equal */
	if(band & 1) {
		++band;
	}
	halfBand = band >> 1;
	
	/* check matrix size */
	if(matrices->NW_q <= band) {
		matrices->NW_q = band << 1;
		free(matrices->D[0]);
		free(matrices->P[0]);
		matrices->D[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->P[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->D[1] = matrices->D[0] + matrices->NW_q;
		matrices->P[1] = matrices->P[0] + matrices->NW_q;
	}
	if(matrices->NW_s <= ((band + 2) * (t_len + 1))) {
		matrices->NW_s = ((band + 3) * (t_len + 2));
		free(matrices->E);
		matrices->E = smalloc(matrices->NW_s);
	}
	
	/* fill in start penalties */
	bq_len = band + 1; /* (band + 1) ~ q_len */
	D_ptr = matrices->D[0];
	D_prev = matrices->D[1];
	P_ptr = matrices->P[0];
	P_prev = matrices->P[1];
	E = matrices->E;
	thisScore = (t_len + q_len) * (MM + U + W1);
	Stat.score = thisScore;
	E_ptr = E + (t_len * (bq_len + 1));
	c_pos = (t_len + q_len) >> 1;
	
	sn = q_len - 1 - (c_pos - halfBand);
	if(k != 2) {
		for(n = sn - 1; n >= 0; --n) {
			D_prev[n] = W1 + (sn - n - 1) * U;
			P_prev[n] = thisScore;
			E_ptr[n] = 3;
		}
		E_ptr[sn - 1] = 18;
		E_ptr[sn] = 0;
		D_prev[sn] = 0;
		P_prev[sn] = 0;
	} else {
		for(n = sn; n >= 0; --n) {
			D_prev[n] = 0;
			P_prev[n] = thisScore;
			E_ptr[n] = 0;
		}
	}
	E_ptr -= (bq_len + 1);
	
	/* Perform banded NW */
	pos[0] = 0;
	pos[1] = 0;
	en = 0;
	c_pos = (t_len + q_len) >> 1;
	for(m = t_len - 1, nuc_pos = t_e - 1; m >= 0; --m, --nuc_pos, --c_pos) {
		
		if(nuc_pos < 0) {
			nuc_pos = template_length - 1; 
		}
		
		/* get banded boundaries, w.r.t. query */
		sq = c_pos + halfBand;
		eq = c_pos - halfBand;
		
		if(eq < 0) {
			eq = 0;
			++en;
		} else {
			en = 0;
		}
		
		/* get start penalties */
		Q_prev = (t_len + q_len) * (MM + U + W1);
		/* check boundaries */
		if(sq < (q_len - 1)) {
			sn = bq_len - 1;
			D_ptr[bq_len] = (t_len + q_len) * (MM + U + W1);
			E_ptr[bq_len] = 37;
		} else {
			sq = q_len - 1;
			sn = en + (q_len - eq);
			D_ptr[sn] = (0 < k) ? 0 : (W1 + (t_len - 1 - m) * U);
			E_ptr[sn] = (0 < k) ? 0 : 37;
			--sn;
		}
		
		t_nuc = getNuc(template, nuc_pos);
		for(n = sn, q_pos = sq; n > en; --q_pos, --n) {
			E_ptr[n] = 0;
			
			/* update Q and P, gap openings */
			Q = D_ptr[n + 1] + W1;
			P_ptr[n] = D_prev[n - 1] + W1;
			if(Q < P_ptr[n]) {
				D_ptr[n] = P_ptr[n];
				e = 4;
			} else {
				D_ptr[n] = Q;
				e = 2;
			}
			
			/* update Q and P, gap extensions */
			/* mark bit 4 and 5 as possible gap-opennings, if necesarry */
			thisScore = Q_prev + U;
			if(Q < thisScore) {
				Q = thisScore;
				if(D_ptr[n] <= thisScore) {
					D_ptr[n] = thisScore;
					e = 3;
				}
			} else {
				E_ptr[n] |= 16;
			}
			thisScore = P_prev[n - 1] + U;
			if(P_ptr[n] < thisScore) {
				P_ptr[n] = thisScore;
				if(D_ptr[n] <= thisScore) {
					D_ptr[n] = thisScore;
					e = 5;
				}
			} else {
				E_ptr[n] |= 32;
			}
			
			/* Update D, match */
			thisScore = D_prev[n] + d[t_nuc][query[q_pos]];
			if(D_ptr[n] <= thisScore) {
				D_ptr[n] = thisScore;
				E_ptr[n] |= 1;
			} else {
				E_ptr[n] |= e;
			}
			
			Q_prev = Q;
		}
		
		/* handle banded boundary */
		E_ptr[n] = 0;
		
		/* update Q gap */
		Q = D_ptr[n + 1] + W1;
		thisScore = Q_prev + U;
		if(Q < thisScore) {
			Q = thisScore;
			e = 3;
		} else {
			e = 2;
			E_ptr[n] |= 16;
		}
		
		/* update unavailable P gap */
		P_ptr[n] = (t_len + q_len) * (MM + U + W1);
		
		/* Update D */
		D_ptr[n] = D_prev[n] + d[t_nuc][query[q_pos]];
		
		/* set D to max, and set E */
		if(Q <= D_ptr[n]) {
			E_ptr[n] |= 1;
		} else {
			D_ptr[n] = Q;
			E_ptr[n] |= e;
		}
		
		/* continue as usual */
		E_ptr -= (bq_len + 1);
		
		if(eq == 0 && k < 0 && Stat.score < D_ptr[n]) {
			Stat.score = D_ptr[n];
			pos[0] = m;
			pos[1] = n;
		}
		
		tmp = D_ptr;
		D_ptr = D_prev;
		D_prev = tmp;
		
		tmp = P_ptr;
		P_ptr = P_prev;
		P_prev = tmp;
	}
	E_ptr = E;
	
	/* get start position of alignment */
	q_pos = 0;
	if(pos[0] == 0) {
		pos[1] = en;
		Stat.score = D_prev[en];
	}
	if(k == -2) {
		for(n = en; n < bq_len; ++n) {
			if(Stat.score <= D_prev[n]) {
				Stat.score = D_prev[n];
				pos[0] = 0;
				pos[1] = n;
				q_pos = n - en;
			}
		}
	}
	aligned->start = q_pos;
	
	/* back tracking */
	m = pos[0];
	n = pos[1];
	E_ptr = E + (m * (bq_len + 1));
	nuc_pos = m + t_s;
	Stat.len = 0;
	Stat.match = 0;
	Stat.tGaps = 0;
	Stat.qGaps = 0;
	while(E_ptr[n] != 0) {
		if(nuc_pos == template_length) {
			nuc_pos = 0;
		}
		if((E_ptr[n] & 7) == 1) {
			aligned->t[Stat.len] = getNuc(template, nuc_pos);
			aligned->q[Stat.len] = query[q_pos];
			aligned->s[Stat.len] = (aligned->t[Stat.len] == aligned->q[Stat.len]) ? '|' : '_';
			++Stat.match;
			++nuc_pos;
			E_ptr += (bq_len + 1);
			++q_pos;
		} else if((E_ptr[n] & 7) >= 4) {
			while(!(E_ptr[n] >> 4)) {
				aligned->t[Stat.len] = getNuc(template, nuc_pos);
				aligned->q[Stat.len] = 5;
				aligned->s[Stat.len] = '_';	
				++nuc_pos;
				E_ptr += (bq_len + 1);
				--n;
				++Stat.len;
				++Stat.qGaps;
			}
			aligned->t[Stat.len] = getNuc(template, nuc_pos);
			aligned->q[Stat.len] = 5;
			aligned->s[Stat.len] = '_';		
			++nuc_pos;
			E_ptr += (bq_len + 1);
			--n;
			++Stat.qGaps;
		} else {
			while(!(E_ptr[n] >> 3)) {
				aligned->t[Stat.len] = 5;
				aligned->q[Stat.len] = query[q_pos];
				aligned->s[Stat.len] = '_';
				++n;
				++q_pos;
				++Stat.len;
				++Stat.tGaps;
			}
			aligned->t[Stat.len] = 5;
			aligned->q[Stat.len] = query[q_pos];
			aligned->s[Stat.len] = '_';
			++n;
			++q_pos;
			++Stat.tGaps;
		}
		++Stat.len;
	}
	aligned->s[Stat.len] = 0;
	aligned->end = q_len - q_pos;
	
	return Stat;
}

AlnScore NW_score(const long unsigned *template, const unsigned char *queryOrg, int k, int t_s, int t_e, int q_s, int q_e, NWmat *matrices, int template_length) {
	
	int m, n, t_len, q_len, thisScore, nuc_pos, W1, U, MM, pos[2];
	int *D_ptr, *D_prev, Q, Q_prev, *P_ptr, *P_prev, *tmp, **d;
	unsigned char *query, t_nuc, *E, *E_ptr, e;
	AlnScore Stat;
	Penalties *rewards;
	
	rewards = matrices->rewards;
	W1 = rewards->W1;
	U = rewards->U;
	MM = rewards->MM;
	d = rewards->d;
	t_len = t_e - t_s;
	q_len = q_e - q_s;
	if(t_len < 0) {
		t_len += template_length;
	}
	query = (unsigned char*)(queryOrg + q_s);
	Stat.pos = 0;
	
	if(t_len == 0 || q_len == 0) {
		if(t_len == q_len) {
			Stat.len = 0;
			Stat.match = 0;
			Stat.tGaps = 0;
			Stat.qGaps = 0;
			Stat.score = 0;
		} else if(t_len == 0) {
			Stat.len = q_len;
			Stat.match = 0;
			Stat.tGaps = q_len;
			Stat.qGaps = 0;
			Stat.score = W1 + (q_len - 1) * U;
		} else {
			Stat.len = t_len;
			Stat.match = 0;
			Stat.tGaps = 0;
			Stat.qGaps = t_len;
			Stat.score = W1 + (t_len - 1) * U;
		}
		return Stat;
	}
	
	/* check matrix size */
	if(matrices->NW_q <= q_len) {
		matrices->NW_q = q_len << 1;
		free(matrices->D[0]);
		free(matrices->P[0]);
		matrices->D[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->P[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->D[1] = matrices->D[0] + matrices->NW_q;
		matrices->P[1] = matrices->P[0] + matrices->NW_q;
	}
	if(matrices->NW_s <= ((q_len + 1) * (t_len + 1))) {
		matrices->NW_s = ((q_len + 2) * (t_len + 2));
		free(matrices->E);
		matrices->E = smalloc(matrices->NW_s);
	}
	
	/* fill in start penalties */
	D_ptr = matrices->D[0];
	D_prev = matrices->D[1];
	P_ptr = matrices->P[0];
	P_prev = matrices->P[1];
	E = matrices->E;
	thisScore = (t_len + q_len) * (MM + U + W1);
	Stat.score = thisScore;
	if(0 < k) {
		E_ptr = E;
		for(m = 0; m < t_len; ++m) {
			E_ptr[q_len] = 0;
			E_ptr += (q_len + 1);
		}
		if(k == 1) {
			for(n = q_len - 1; n >= 0; --n) {
				D_prev[n] = W1 + (q_len - 1 - n) * U;
				P_prev[n] = thisScore;
				E_ptr[n] = 3;
			}
			E_ptr[q_len - 1] = 18;
			E_ptr[q_len] = 0;
			D_prev[q_len] = 0;
			P_prev[q_len] = 0;
		} else {
			for(n = q_len; n >= 0; --n) {
				D_prev[n] = 0;
				P_prev[n] = thisScore;
				E_ptr[n] = 0;
			}
		}
	} else {
		E_ptr = E;
		for(m = 0; m < t_len; ++m) {
			E_ptr[q_len] = 5;
			E_ptr += (q_len + 1);
		}
		E_ptr[-1] = 36;
		
		for(n = q_len - 1; n >= 0; --n) {
			D_prev[n] = W1 + (q_len - 1 - n) * U;
			P_prev[n] = thisScore;
			E_ptr[n] = 3;
		}
		E_ptr[q_len - 1] = 18;
		E_ptr[q_len] = 0;
		D_prev[q_len] = 0;
		P_prev[q_len] = 0;
	}
	E_ptr -= (q_len + 1);
	
	/* Perform NW */
	pos[0] = 0;
	for(m = t_len - 1, nuc_pos = t_e - 1; m >= 0; --m, --nuc_pos) {
		
		if(nuc_pos < 0) {
			nuc_pos = template_length - 1;
		}
		
		D_ptr[q_len] = (0 < k) ? 0 : (W1 + (t_len - 1 - m) * U);
		Q_prev = (t_len + q_len) * (MM + U + W1);
		
		t_nuc = getNuc(template, nuc_pos);
		for(n = q_len - 1; n >= 0; --n) {
			E_ptr[n] = 0;
			
			/* update Q and P, gap openings */
			Q = D_ptr[n + 1] + W1;
			P_ptr[n] = D_prev[n] + W1;
			if(Q < P_ptr[n]) {
				D_ptr[n] = P_ptr[n];
				e = 4;
			} else {
				D_ptr[n] = Q;
				e = 2;
			}
			
			/* update Q and P, gap extensions */
			/* mark bit 4 and 5 as possible gap-opennings, if necesarry */
			thisScore = Q_prev + U;
			if(Q < thisScore) {
				Q = thisScore;
				if(D_ptr[n] <= thisScore) {
					D_ptr[n] = thisScore;
					e = 3;
				}
			} else {
				E_ptr[n] |= 16;
			}
			thisScore = P_prev[n] + U;
			if(P_ptr[n] < thisScore) {
				P_ptr[n] = thisScore;
				if(D_ptr[n] <= thisScore) {
					D_ptr[n] = thisScore;
					e = 5;
				}
			} else {
				E_ptr[n] |= 32;
			}
			
			/* Update D, match */
			thisScore = D_prev[n + 1] + d[t_nuc][query[n]];
			if(D_ptr[n] <= thisScore) {
				D_ptr[n] = thisScore;
				E_ptr[n] |= 1;
			} else {
				E_ptr[n] |= e;
			}
			
			Q_prev = Q;
		}
		
		E_ptr -= (q_len + 1);
		
		if(k < 0 && Stat.score < *D_ptr) {
			Stat.score = *D_ptr;
			pos[0] = m;
		}
		
		tmp = D_ptr;
		D_ptr = D_prev;
		D_prev = tmp;
		
		tmp = P_ptr;
		P_ptr = P_prev;
		P_prev = tmp;
	}
	E_ptr = E;
	
	/* get start position of alignment */
	if(k < 0) {
		pos[1] = 0;
		if(k == -2) {
			for(n = 0; n < q_len; ++n) {
				if(Stat.score <= D_prev[n]) {
					Stat.score = D_prev[n];
					pos[0] = 0;
					pos[1] = n;
				}
			}
		}
	} else {
		Stat.score = *D_prev;
		pos[0] = 0;
		pos[1] = 0;
	}
	
	/* make back tracking */
	m = pos[0];
	E_ptr = E + (m * (q_len + 1));
	n = pos[1];
	nuc_pos = m + t_s;
	Stat.len = 0;
	Stat.match = 0;
	Stat.tGaps = 0;
	Stat.qGaps = 0;
	while(E_ptr[n] != 0) {
		if(nuc_pos == template_length) {
			nuc_pos = 0;
		}
		if((E_ptr[n] & 7) == 1) {
			++Stat.match;
			++nuc_pos;
			E_ptr += (q_len + 1);
			++n;
		} else if((E_ptr[n] & 7) >= 4) {
			while(!(E_ptr[n] >> 4)) {
				++nuc_pos;
				E_ptr += (q_len + 1);
				++Stat.len;
				++Stat.qGaps;
			}
			++Stat.qGaps;
			++nuc_pos;
			E_ptr += (q_len + 1);
		} else {
			while(!(E_ptr[n] >> 3)) {
				++n;
				++Stat.len;
				++Stat.tGaps;
			}
			++Stat.tGaps;
			++n;
		}
		++Stat.len;
	}
	
	return Stat;
}

AlnScore NW_band_score(const long unsigned *template, const unsigned char *queryOrg, int k, int t_s, int t_e, int q_s, int q_e, int band, NWmat *matrices, int template_length) {
	
	int m, n, t_len, q_len, thisScore, nuc_pos, pos[2];
	int bq_len, halfBand, sn, en, sq, eq, q_pos, c_pos, W1, U, MM;
	int *D_ptr, *D_prev, Q, Q_prev, *P_ptr, *P_prev, *tmp, **d;
	unsigned char *query, t_nuc, *E, *E_ptr, e;
	AlnScore Stat;
	Penalties *rewards;
	
	rewards = matrices->rewards;
	W1 = rewards->W1;
	U = rewards->U;
	MM = rewards->MM;
	d = rewards->d;
	t_len = t_e - t_s;
	q_len = q_e - q_s;
	if(t_len < 0) {
		t_len += template_length;
	}
	query = (unsigned char*)(queryOrg + q_s);
	Stat.pos = 0;
	
	if(t_len == 0 || q_len == 0) {
		if(t_len == q_len) {
			Stat.len = 0;
			Stat.match = 0;
			Stat.tGaps = 0;
			Stat.qGaps = 0;
			Stat.score = 0;
		} else if(t_len == 0) {
			Stat.len = q_len;
			Stat.match = 0;
			Stat.tGaps = q_len;
			Stat.qGaps = 0;
			Stat.score = W1 + (q_len - 1) * U;
		} else {
			Stat.len = t_len;
			Stat.match = 0;
			Stat.tGaps = 0;
			Stat.qGaps = t_len;
			Stat.score = W1 + (t_len - 1) * U;
		}
		return Stat;
	}
	
	/* ensure that band is equal */
	if(band & 1) {
		++band;
	}
	halfBand = band >> 1;
	
	/* check matrix size */
	if(matrices->NW_q <= band) {
		matrices->NW_q = band << 1;
		free(matrices->D[0]);
		free(matrices->P[0]);
		matrices->D[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->P[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->D[1] = matrices->D[0] + matrices->NW_q;
		matrices->P[1] = matrices->P[0] + matrices->NW_q;
	}
	if(matrices->NW_s <= ((band + 2) * (t_len + 1))) {
		matrices->NW_s = ((band + 3) * (t_len + 2));
		free(matrices->E);
		matrices->E = smalloc(matrices->NW_s);
	}
	
	/* fill in start penalties */
	bq_len = band + 1; /* (band + 1) ~ q_len */
	D_ptr = matrices->D[0];
	D_prev = matrices->D[1];
	P_ptr = matrices->P[0];
	P_prev = matrices->P[1];
	E = matrices->E;
	thisScore = (t_len + q_len) * (MM + U + W1);
	Stat.score = thisScore;
	E_ptr = E + (t_len * (bq_len + 1));
	c_pos = (t_len + q_len) >> 1;
	
	sn = q_len - 1 - (c_pos - halfBand);
	if(k != 2) {
		for(n = sn - 1; n >= 0; --n) {
			D_prev[n] = W1 + (sn - n - 1) * U;
			P_prev[n] = thisScore;
			E_ptr[n] = 3;
		}
		E_ptr[sn - 1] = 18;
		E_ptr[sn] = 0;
		D_prev[sn] = 0;
		P_prev[sn] = 0;
	} else {
		for(n = sn; n >= 0; --n) {
			D_prev[n] = 0;
			P_prev[n] = thisScore;
			E_ptr[n] = 0;
		}
	}
	E_ptr -= (bq_len + 1);
	
	
	/* Perform banded NW */
	pos[0] = 0;
	pos[1] = 0;
	en = 0;
	c_pos = (t_len + q_len) >> 1;
	for(m = t_len - 1, nuc_pos = t_e - 1; m >= 0; --m, --nuc_pos, --c_pos) {
		
		if(nuc_pos < 0) {
			nuc_pos = template_length - 1; 
		}
		
		/* get banded boundaries, w.r.t. query */
		sq = c_pos + halfBand;
		eq = c_pos - halfBand;
		
		if(eq < 0) {
			eq = 0;
			++en;
		} else {
			en = 0;
		}
		
		/* get start penalties */
		Q_prev = (t_len + q_len) * (MM + U + W1);
		/* check boundaries */
		if(sq < (q_len - 1)) {
			sn = bq_len - 1;
			D_ptr[bq_len] = (t_len + q_len) * (MM + U + W1);
			E_ptr[bq_len] = 37;
		} else {
			sq = q_len - 1;
			sn = en + (q_len - eq);
			D_ptr[sn] = (0 < k) ? 0 : (W1 + (t_len - 1 - m) * U);
			E_ptr[sn] = (0 < k) ? 0 : 37;
			--sn;
		}
		
		t_nuc = getNuc(template, nuc_pos);
		for(n = sn, q_pos = sq; n > en; --q_pos, --n) {
			E_ptr[n] = 0;
			
			/* update Q and P, gap openings */
			Q = D_ptr[n + 1] + W1;
			P_ptr[n] = D_prev[n - 1] + W1;
			if(Q < P_ptr[n]) {
				D_ptr[n] = P_ptr[n];
				e = 4;
			} else {
				D_ptr[n] = Q;
				e = 2;
			}
			
			/* update Q and P, gap extensions */
			/* mark bit 4 and 5 as possible gap-opennings, if necesarry */
			thisScore = Q_prev + U;
			if(Q < thisScore) {
				Q = thisScore;
				if(D_ptr[n] <= thisScore) {
					D_ptr[n] = thisScore;
					e = 3;
				}
			} else {
				E_ptr[n] |= 16;
			}
			thisScore = P_prev[n - 1] + U;
			if(P_ptr[n] < thisScore) {
				P_ptr[n] = thisScore;
				if(D_ptr[n] <= thisScore) {
					D_ptr[n] = thisScore;
					e = 5;
				}
			} else {
				E_ptr[n] |= 32;
			}
			
			/* Update D, match */
			thisScore = D_prev[n] + d[t_nuc][query[q_pos]];
			if(D_ptr[n] <= thisScore) {
				D_ptr[n] = thisScore;
				E_ptr[n] |= 1;
			} else {
				E_ptr[n] |= e;
			}
			
			Q_prev = Q;
		}
		
		/* handle banded boundary */
		E_ptr[n] = 0;
		
		/* update Q gap */
		Q = D_ptr[n + 1] + W1;
		thisScore = Q_prev + U;
		if(Q < thisScore) {
			Q = thisScore;
			e = 3;
		} else {
			e = 2;
			E_ptr[n] |= 16;
		}
		
		/* update unavailable P gap */
		P_ptr[n] = (t_len + q_len) * (MM + U + W1);
		
		/* Update D */
		D_ptr[n] = D_prev[n] + d[t_nuc][query[q_pos]];
		
		/* set D to max, and set E */
		if(Q <= D_ptr[n]) {
			E_ptr[n] |= 1;
		} else {
			D_ptr[n] = Q;
			E_ptr[n] |= e;
		}
		
		/* continue as usual */
		E_ptr -= (bq_len + 1);
		
		if(eq == 0 && k < 0 && Stat.score < D_ptr[n]) {
			Stat.score = D_ptr[n];
			pos[0] = m;
			pos[1] = n;
		}
		
		tmp = D_ptr;
		D_ptr = D_prev;
		D_prev = tmp;
		
		tmp = P_ptr;
		P_ptr = P_prev;
		P_prev = tmp;
	}
	E_ptr = E;
	
	/* get start position of alignment */
	q_pos = 0;
	if(pos[0] == 0) {
		pos[1] = en;
		Stat.score = D_prev[en];
	}
	if(k == -2) {
		for(n = en; n < bq_len; ++n) {
			if(Stat.score <= D_prev[n]) {
				Stat.score = D_prev[n];
				pos[0] = 0;
				pos[1] = n;
				q_pos = n;
			}
		}
	}
	
	/* back tracking */
	m = pos[0];
	n = pos[1];
	E_ptr = E + (m * (bq_len + 1));
	nuc_pos = m + t_s;
	Stat.len = 0;
	Stat.match = 0;
	Stat.tGaps = 0;
	Stat.qGaps = 0;
	while(E_ptr[n] != 0) {
		if(nuc_pos == template_length) {
			nuc_pos = 0;
		}
		if((E_ptr[n] & 7) == 1) {
			++Stat.match;
			++nuc_pos;
			E_ptr += (bq_len + 1);
			++q_pos;
		} else if((E_ptr[n] & 7) >= 4) {
			while(!(E_ptr[n] >> 4)) {
				++nuc_pos;
				E_ptr += (bq_len + 1);
				--n;
				++Stat.len;
				++Stat.qGaps;
			}
			++Stat.qGaps;
			++nuc_pos;
			E_ptr += (bq_len + 1);
			--n;
		} else {
			while(!(E_ptr[n] >> 3)) {
				++n;
				++q_pos;
				++Stat.len;
				++Stat.tGaps;
			}
			++Stat.tGaps;
			++n;
			++q_pos;
		}
		++Stat.len;
	}
	
	return Stat;
}
