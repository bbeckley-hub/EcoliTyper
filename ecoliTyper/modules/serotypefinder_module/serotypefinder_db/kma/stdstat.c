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

#include <math.h>
#include "stdstat.h"

int (*cmp)(int, int) = &cmp_or;

int cmp_or(int t, int q) {
	return (t || q);
}

int cmp_and(int t, int q) {
	return (t && q);
}

int cmp_true(int t, int q) {
	return 1;
}

double fastp(long double q) {
	/* P-value from quantile in a chi-square distribution */
	double p = 1.0;
	if(q > 114.5242) {
		p = 1e-26;
	} else if(q > 109.9604) {
		p = 1e-25;
	} else if(q > 105.3969) {
		p = 1e-24;
	} else if(q > 100.8337) {
		p = 1e-23;
	} else if(q > 96.27476) {
		p = 1e-22;
	} else if(q > 91.71701) {
		p = 1e-21;
	} else if(q > 87.16164) {
		p = 1e-20;
	} else if(q > 82.60901) {
		p = 1e-19;
	} else if(q > 78.05917) {
		p = 1e-18;
	} else if(q > 73.51245) {
		p = 1e-17;
	} else if(q > 68.96954) {
		p = 1e-16;
	} else if(q > 64.43048) {
		p = 1e-15;
	} else if(q > 59.89615) {
		p = 1e-14;
	} else if(q > 55.36699) {
		p = 1e-13;
	} else if(q > 50.84417) {
		p = 1e-12;
	} else if(q > 46.32844) {
		p = 1e-11;
	} else if(q > 41.82144) {
		p = 1e-10;
	} else if(q > 37.32489) {
		p = 1e-9;
	} else if(q > 32.84127) {
		p = 1e-8;
	} else if(q > 28.37395) {
		p = 1e-7;
	} else if(q > 23.92814) {
		p = 1e-6;
	} else if(q > 19.51139) {
		p = 1e-5;
	} else if(q > 15.13671) {
		p = 1e-4;
	} else if(q > 10.82759) {
		p = 1e-3;
	} else if(q > 6.634897) {
		p = 0.01;
	} else if(q > 3.841443) {
		p = 0.05;
	} else if(q > 2.705532) {
		p = 0.1;
	} else if(q > 2.072251) {
		p = 0.15;
	} else if(q > 1.642374) {
		p = 0.2;
	} else if(q > 1.323304) {
		p = 0.25;
	} else if(q > 1.074194) {
		p = 0.3;
	} else if(q > 0.8734571) {
		p = 0.35;
	} else if(q > 0.7083263) {
		p = 0.4;
	} else if(q > 0.5706519) {
		p = 0.45;
	} else if(q > 0.4549364) {
		p = 0.5;
	} else if(q > 0.3573172) {
		p = 0.55;
	} else if(q > 0.2749959) {
		p = 0.6;
	} else if(q > 0.2059001) {
		p = 0.65;
	} else if(q > 0.1484719) {
		p = 0.7;
	} else if(q > 0.1015310) {
		p = 0.75;
	} else if(q > 0.06418475) {
		p = 0.8;
	} else if(q > 0.03576578) {
		p = 0.85;
	} else if(q > 0.01579077) {
		p = 0.9;
	} else if(q > 0.00393214) {
		p = 0.95;
	} else if(q >= 0.0) {
		p = 1.0;
	} else {
		p = 1.00 - fastp(-1 * q);
	}
	return p;
}

double p_chisqr(long double q) {
	
	if(q < 0) {
		/* Handle negative overflow */
		return 1e-26;
	} else if(q > 49) {
		/* Get p-val from table, to avoid overflow */
		return fastp(q);
	}
	/* claculate p-value */
	return 1 - 1.772453850 * erf(sqrt(0.5 * q)) / tgamma(0.5);
}

double power(double x, unsigned n) {
	
	double y;
	
	if(n) {
		y = power(x, n >> 1);
		return (n & 1) ? y * y * x : y * y;
	}
	
	return 1.0;
}

double binP(int n, int k, double p) {
	
	int i, j, nk;
	double P, q, pq;
	
	/*
		P = n! / (k! (n-k)!) * p^k * q^(n - k)
	*/
	
	q = 1 - p;
	if(k == 0) {
		P = power(q, n);
		return P != 0.0 ? P : 1.0e-308;
	} else if(n == k) {
		P = power(p, n);
		return P != 0.0 ? P : 1.0e-308;
	} else if(p == 0 || q == 0) {
		return 0.0;
	}
	
	P = 1.0;
	nk = n - k;
	pq = p * q;
	i = n + 1;
	if(k < nk) {
		j = k + 1;
	} else {
		j = nk + 1;
	}
	
	while(--j) {
		P *= (--i * pq / j);
	}
	
	if(nk < k) {
		P *= power(p, k - nk);
	} else if(k < nk) {
		P *= power(q, nk - k);
	}
	
	return P != 0.0 ? P : 1.0e-308;
}

unsigned minimum(unsigned *src, unsigned n) {
	
	unsigned min;
	
	min = *src;
	while(--n) {
		if(*++src < min) {
			min = *src;
		}
	}
	
	return min;
}

/*
double eQual(unsigned char *qual, int len, int phredScale, int minQ) {
	
	//E(Q) = -10 * log_10(sum(10^(-Q/10)) / |Q|) 
	unsigned i;
	double sum;
	
	if(!len || !minQ) {
		return 0;
	}
	
	// sum quality scores
	i = len;
	sum = pow(10, (phredScale - *qual) / 10.0);
	while(--i) {
		sum += pow(10, (phredScale - *++qual) / 10.0);
	}
	
	// return average
	return -10 * log10(((double)(sum)) / len);
}
*/

unsigned eQual(unsigned char *qual, const int len, const int minQ, const double *prob) {
	
	/*
	static const double prob[128] = {
		1.0000000000000000, 0.7943282347242815, 0.6309573444801932, 0.5011872336272722, 0.3981071705534972, 0.3162277660168379, 0.2511886431509580, 0.1995262314968880,
		0.1584893192461113, 0.1258925411794167, 0.1000000000000000, 0.0794328234724281, 0.0630957344480193, 0.0501187233627272, 0.0398107170553497, 0.0316227766016838,
		0.0251188643150958, 0.0199526231496888, 0.0158489319246111, 0.0125892541179417, 0.0100000000000000, 0.0079432823472428, 0.0063095734448019, 0.0050118723362727,
		0.0039810717055350, 0.0031622776601684, 0.0025118864315096, 0.0019952623149689, 0.0015848931924611, 0.0012589254117942, 0.0010000000000000, 0.0007943282347243,
		0.0006309573444802, 0.0005011872336273, 0.0003981071705535, 0.0003162277660168, 0.0002511886431510, 0.0001995262314969, 0.0001584893192461, 0.0001258925411794,
		0.0001000000000000, 0.0000794328234724, 0.0000630957344480, 0.0000501187233627, 0.0000398107170553, 0.0000316227766017, 0.0000251188643151, 0.0000199526231497,
		0.0000158489319246, 0.0000125892541179, 0.0000100000000000, 0.0000079432823472, 0.0000063095734448, 0.0000050118723363, 0.0000039810717055, 0.0000031622776602,
		0.0000025118864315, 0.0000019952623150, 0.0000015848931925, 0.0000012589254118, 0.0000010000000000, 0.0000007943282347, 0.0000006309573445, 0.0000005011872336,
		0.0000003981071706, 0.0000003162277660, 0.0000002511886432, 0.0000001995262315, 0.0000001584893192, 0.0000001258925412, 0.0000001000000000, 0.0000000794328235,
		0.0000000630957344, 0.0000000501187234, 0.0000000398107171, 0.0000000316227766, 0.0000000251188643, 0.0000000199526231, 0.0000000158489319, 0.0000000125892541,
		0.0000000100000000, 0.0000000079432823, 0.0000000063095734, 0.0000000050118723, 0.0000000039810717, 0.0000000031622777, 0.0000000025118864, 0.0000000019952623,
		0.0000000015848932, 0.0000000012589254, 0.0000000010000000, 0.0000000007943282, 0.0000000006309573, 0.0000000005011872, 0.0000000003981072, 0.0000000003162278,
		0.0000000002511886, 0.0000000001995262, 0.0000000001584893, 0.0000000001258925, 0.0000000001000000, 0.0000000000794328, 0.0000000000630957, 0.0000000000501187,
		0.0000000000398107, 0.0000000000316228, 0.0000000000251189, 0.0000000000199526, 0.0000000000158489, 0.0000000000125893, 0.0000000000100000, 0.0000000000079433,
		0.0000000000063096, 0.0000000000050119, 0.0000000000039811, 0.0000000000031623, 0.0000000000025119, 0.0000000000019953, 0.0000000000015849, 0.0000000000012589,
		0.0000000000010000, 0.0000000000007943, 0.0000000000006310, 0.0000000000005012, 0.0000000000003981, 0.0000000000003162, 0.0000000000002512, 0.0000000000001995};
	*/
	/*
	E(Q) = -10 * log_10(sum(10^(-Q/10)) / |Q|) 
	*/
	unsigned i;
	double sum;
	
	if(!len || !minQ) {
		return 0;
	}
	
	/* sum quality scores */
	i = len;
	sum = prob[*qual];
	while(--i) {
		sum += prob[*++qual];
	}
	
	/* return average */
	return ceil(-10 * log10(sum / len));
}
