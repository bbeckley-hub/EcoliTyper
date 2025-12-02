/* Philip T.L.C. Clausen May 2021 plan@dtu.dk */

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
#include <stdio.h>
#include "compdna.h"
#include "pherror.h"
#include "runinput.h"
#include "qc.h"
#include "qseqs.h"
#include "trim.h"

void printTrimFsa(Qseqs *header, Qseqs *qseq, Qseqs *qual, CompDNA *compressor, FILE *out) {
	
	const char bases[6] = "ACGTN-";
	static FILE *OUT = 0;
	int i, error;
	unsigned char *ptr;
	
	/* init */
	if(!OUT) {
		OUT = out;
		return;
	}
	error = 0;
	
	/* translate to human readable bases */
	ptr = qseq->seq;
	i = qseq->len + 1;
	while(--i) {
		*ptr = bases[*ptr];
		++ptr;
	}
	
	/* print */
	if(qual) {
		sfwrite(header->seq, sizeof(unsigned char), header->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(qseq->seq, sizeof(unsigned char), qseq->len, OUT);
		sfwrite("\n+\n", sizeof(unsigned char), 3, OUT);
		sfwrite(qual->seq, sizeof(unsigned char), qual->len, OUT);
		error += (putc('\n', OUT) == EOF);
	} else {
		sfwrite(header->seq, sizeof(unsigned char), header->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(qseq->seq, sizeof(unsigned char), qseq->len, OUT);
		error += (putc('\n', OUT) == EOF);
	}
	
	if(error) {
		ERROR();
	}
}

void printTrimFsa_pair(Qseqs *header, Qseqs *qseq, Qseqs *qual, Qseqs *header_r, Qseqs *qseq_r, Qseqs *qual_r, CompDNA *compressor, FILE *out) {
	
	const char bases[6] = "ACGTN-";
	static FILE *OUT = 0;
	int i, error;
	unsigned char *ptr;
	
	/* init */
	if(!OUT) {
		OUT = out;
		return;
	}
	error = 0;
	
	/* translate to human readable bases */
	ptr = qseq->seq;
	i = qseq->len + 1;
	while(--i) {
		*ptr = bases[*ptr];
		++ptr;
	}
	ptr = qseq_r->seq;
	i = qseq_r->len + 1;
	while(--i) {
		*ptr = bases[*ptr];
		++ptr;
	}
	
	/* print */
	if(qual) {
		sfwrite(header->seq, sizeof(unsigned char), header->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(qseq->seq, sizeof(unsigned char), qseq->len, OUT);
		sfwrite("\n+\n", sizeof(unsigned char), 3, OUT);
		sfwrite(qual->seq, sizeof(unsigned char), qual->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(header_r->seq, sizeof(unsigned char), header_r->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(qseq_r->seq, sizeof(unsigned char), qseq_r->len, OUT);
		sfwrite("\n+\n", sizeof(unsigned char), 3, OUT);
		sfwrite(qual_r->seq, sizeof(unsigned char), qual_r->len, OUT);
		error += (putc('\n', OUT) == EOF);
	} else {
		sfwrite(header->seq, sizeof(unsigned char), header->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(qseq->seq, sizeof(unsigned char), qseq->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(header_r->seq, sizeof(unsigned char), header_r->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(qseq_r->seq, sizeof(unsigned char), qseq_r->len, OUT);
		error += (putc('\n', OUT) == EOF);
	}
	
	if(error) {
		ERROR();
	}
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#kma trim trims sequences\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file(s)", "STDIN");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ipe", "Paired input files", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-int", "Inerleaved input file(s)", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "STDOUT");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-qc", "Report QC, repeat for verbose", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ml", "Minimum length", "16");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-xl", "Maximum length", "2147483647");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mp", "Minimum phred", "20");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mi", "Minimum internal phred score", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-eq", "Minimum average quality", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-5p", "Trim 5 prime", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-3p", "Trim 3 prime", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	
	return (out == stderr);
}

int trim_main(int argc, char *argv[]) {
	
	static const double prob[256] = {1.00000000000000000000000000000000, 0.79432823472428149003121689020190, 0.63095734448019324958067954867147, 0.50118723362727224390766878059367, 0.39810717055349720272516833574628, 0.31622776601683794117647607890831, 0.25118864315095795758381314044527, 0.19952623149688791803768594945723,
		0.15848931924611134314240246112604, 0.12589254117941672816982645599637, 0.10000000000000000555111512312578, 0.07943282347242813790089144276862, 0.06309573444801930275360746236402, 0.05011872336272722022743053571503, 0.03981071705534971333362292966740, 0.03162277660168379134209004632794,
		0.02511886431509579437060253326308, 0.01995262314968878833432164299211, 0.01584893192461113431424024611260, 0.01258925411794166795975691286458, 0.01000000000000000020816681711722, 0.00794328234724281379008914427686, 0.00630957344480193027536074623640, 0.00501187233627271942065783960629,
		0.00398107170553496907822177419689, 0.00316227766016837939441752602932, 0.00251188643150957943706025332631, 0.00199526231496887892016833809805, 0.00158489319246111256406228662286, 0.00125892541179416618882247469458, 0.00100000000000000002081668171172, 0.00079432823472428131395678407856,
		0.00063095734448019298416798772422, 0.00050118723362727198543387086005, 0.00039810717055349691866419914454, 0.00031622776601683793944175260293, 0.00025118864315095795454804705749, 0.00019952623149688787575380122252, 0.00015848931924611125640622866229, 0.00012589254117941660804022574460,
		0.00010000000000000000479217360239, 0.00007943282347242805821203176508, 0.00006309573444801929299578790999, 0.00005011872336272725139824299467, 0.00003981071705534969457692534567, 0.00003162277660168379529942797590, 0.00002511886431509577106025582482, 0.00001995262314968878622012740665,
		0.00001584893192461110734471120554, 0.00001258925411794166114283575336, 0.00001000000000000000081803053914, 0.00000794328234724280480476363980, 0.00000630957344480192963839196990, 0.00000501187233627271480602234296, 0.00000398107170553496911887935567, 0.00000316227766016837919112961869,
		0.00000251188643150957719072887721, 0.00000199526231496887870671603539, 0.00000158489319246111090387771000, 0.00000125892541179416611428357534, 0.00000099999999999999995474811183, 0.00000079432823472428050165218766, 0.00000063095734448019296383919699, 0.00000050118723362727143825058693,
		0.00000039810717055349692247584741, 0.00000031622776601683791911296187, 0.00000025118864315095769789706404, 0.00000019952623149688787067160354, 0.00000015848931924611109038777100, 0.00000012589254117941661672231345, 0.00000009999999999999999547481118, 0.00000007943282347242804751824081,
		0.00000006309573444801929638391970, 0.00000005011872336272714382505869, 0.00000003981071705534968960060678, 0.00000003162277660168379191129619, 0.00000002511886431509577111319538, 0.00000001995262314968878640541586, 0.00000001584893192461110970052159, 0.00000001258925411794166101048686,
		0.00000001000000000000000020922561, 0.00000000794328234724282195718082, 0.00000000630957344480191723068278, 0.00000000501187233627271471337811, 0.00000000398107170553496896006068, 0.00000000316227766016837952200186, 0.00000000251188643150958199168515, 0.00000000199526231496887450463852,
		0.00000000158489319246111092869313, 0.00000000125892541179416626648481, 0.00000000100000000000000006228159, 0.00000000079432823472428217503857, 0.00000000063095734448019166102973, 0.00000000050118723362727142997878, 0.00000000039810717055349693736510, 0.00000000031622776601683795220019,
		0.00000000025118864315095716519275, 0.00000000019952623149688748148313, 0.00000000015848931924611108769943, 0.00000000012589254117941661630872, 0.00000000010000000000000000364322, 0.00000000007943282347242789438643, 0.00000000006309573444801916610297, 0.00000000005011872336272714816776,
		0.00000000003981071705534969502898, 0.00000000003162277660168379392755, 0.00000000002511886431509571975045, 0.00000000001995262314968874750208, 0.00000000001584893192461110747747, 0.00000000001258925411794166163087, 0.00000000000999999999999999939497, 0.00000000000794328234724278911553,
		0.00000000000630957344480191677186, 0.00000000000501187233627271465522, 0.00000000000398107170553496950290, 0.00000000000316227766016837939275, 0.00000000000251188643150957181349, 0.00000000000199526231496887466943, 0.00000000000158489319246111099009, 0.00000000000125892541179416608231,
		0.00000000000099999999999999997989, 0.00000000000079432823472427893175, 0.00000000000063095734448019171758, 0.00000000000050118723362727146552, 0.00000000000039810717055349692000, 0.00000000000031622776601683791908, 0.00000000000025118864315095719145, 0.00000000000019952623149688746694,
		0.00000000000015848931924611108891, 0.00000000000012589254117941662843, 0.00000000000010000000000000000304, 0.00000000000007943282347242789317, 0.00000000000006309573444801916418, 0.00000000000005011872336272714403, 0.00000000000003981071705534969326, 0.00000000000003162277660168379569,
		0.00000000000002511886431509571851, 0.00000000000001995262314968874606, 0.00000000000001584893192461110952, 0.00000000000001258925411794166158, 0.00000000000000999999999999999999, 0.00000000000000794328234724278869, 0.00000000000000630957344480191689, 0.00000000000000501187233627271456,
		0.00000000000000398107170553496948, 0.00000000000000316227766016837941, 0.00000000000000251188643150957185, 0.00000000000000199526231496887469, 0.00000000000000158489319246111095, 0.00000000000000125892541179416628, 0.00000000000000100000000000000008, 0.00000000000000079432823472427887,
		0.00000000000000063095734448019173, 0.00000000000000050118723362727146, 0.00000000000000039810717055349695, 0.00000000000000031622776601683793, 0.00000000000000025118864315095717, 0.00000000000000019952623149688748, 0.00000000000000015848931924611109, 0.00000000000000012589254117941662,
		0.00000000000000010000000000000000, 0.00000000000000007943282347242789, 0.00000000000000006309573444801943, 0.00000000000000005011872336272714, 0.00000000000000003981071705534953, 0.00000000000000003162277660168380, 0.00000000000000002511886431509572, 0.00000000000000001995262314968883,
		0.00000000000000001584893192461111, 0.00000000000000001258925411794161, 0.00000000000000001000000000000000, 0.00000000000000000794328234724279, 0.00000000000000000630957344480194, 0.00000000000000000501187233627271, 0.00000000000000000398107170553495, 0.00000000000000000316227766016838,
		0.00000000000000000251188643150957, 0.00000000000000000199526231496888, 0.00000000000000000158489319246111, 0.00000000000000000125892541179416, 0.00000000000000000100000000000000, 0.00000000000000000079432823472428, 0.00000000000000000063095734448019, 0.00000000000000000050118723362727,
		0.00000000000000000039810717055350, 0.00000000000000000031622776601684, 0.00000000000000000025118864315096, 0.00000000000000000019952623149689, 0.00000000000000000015848931924611, 0.00000000000000000012589254117942, 0.00000000000000000010000000000000, 0.00000000000000000007943282347243,
		0.00000000000000000006309573444802, 0.00000000000000000005011872336273, 0.00000000000000000003981071705535, 0.00000000000000000003162277660168, 0.00000000000000000002511886431510, 0.00000000000000000001995262314969, 0.00000000000000000001584893192461, 0.00000000000000000001258925411794,
		0.00000000000000000001000000000000, 0.00000000000000000000794328234724, 0.00000000000000000000630957344480, 0.00000000000000000000501187233627, 0.00000000000000000000398107170553, 0.00000000000000000000316227766017, 0.00000000000000000000251188643151, 0.00000000000000000000199526231497,
		0.00000000000000000000158489319246, 0.00000000000000000000125892541179, 0.00000000000000000000100000000000, 0.00000000000000000000079432823472, 0.00000000000000000000063095734448, 0.00000000000000000000050118723363, 0.00000000000000000000039810717055, 0.00000000000000000000031622776602,
		0.00000000000000000000025118864315, 0.00000000000000000000019952623150, 0.00000000000000000000015848931925, 0.00000000000000000000012589254118, 0.00000000000000000000010000000000, 0.00000000000000000000007943282347, 0.00000000000000000000006309573445, 0.00000000000000000000005011872336,
		0.00000000000000000000003981071706, 0.00000000000000000000003162277660, 0.00000000000000000000002511886432, 0.00000000000000000000001995262315, 0.00000000000000000000001584893192, 0.00000000000000000000001258925412, 0.00000000000000000000001000000000, 0.00000000000000000000000794328235,
		0.00000000000000000000000630957344, 0.00000000000000000000000501187234, 0.00000000000000000000000398107171, 0.00000000000000000000000316227766, 0.00000000000000000000000251188643, 0.00000000000000000000000199526231, 0.00000000000000000000000158489319, 0.00000000000000000000000125892541,
		0.00000000000000000000000100000000, 0.00000000000000000000000079432823, 0.00000000000000000000000063095734, 0.00000000000000000000000050118723, 0.00000000000000000000000039810717, 0.00000000000000000000000031622777, 0.00000000000000000000000025118864, 0.00000000000000000000000019952623,
		0.00000000000000000000000015848932, 0.00000000000000000000000012589254, 0.00000000000000000000000010000000, 0.00000000000000000000000007943282, 0.00000000000000000000000006309573, 0.00000000000000000000000005011872, 0.00000000000000000000000003981072, 0.00000000000000000000000003162278};
	int i, args, fileCounter, fileCounter_PE, fileCounter_INT, fileCount;
	int minPhred, minmaskQ, minQ, fiveClip, threeClip, minlen, maxlen;
	long unsigned totFrags;
	char **inputfiles, **inputfiles_PE, **inputfiles_INT;
	char *to2Bit, *outputfilename, *exeBasic;
	FILE *out, *out_int, *out_json;
	QCstat *qcreport; 
	
	/* set defaults */
	minPhred = 20;
	minmaskQ = 0;
	minQ = 0;
	fiveClip = 0;
	threeClip = 0;
	minlen = 16;
	maxlen = 2147483647;
	fileCounter = 0;
	fileCounter_PE = 0;
	fileCounter_INT = 0;
	inputfiles = 0;
	inputfiles_PE = 0;
	inputfiles_INT = 0;
	outputfilename = 0;
	out = stdout;
	out_int = stdout;
	out_json = stderr;
	qcreport = 0;
	
	/* PARSE COMMAND LINE OPTIONS */
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-i") == 0) {
			++args;
			fileCount = fileCounter;
			for(i = args; i < argc && (strncmp(argv[i], "-", 1) != 0 || strcmp(argv[i], "--") == 0); ++i) {
				++fileCounter;
			}
			if(fileCounter == 0) {
				fprintf(stderr, "No files were specified.\n");
				exit(1);
			} else {
				inputfiles = realloc(inputfiles, fileCounter * sizeof(char *));
				if(!inputfiles) {
					ERROR();
				}
			}
			
			for(i = fileCount; i < fileCounter; ++i, ++args) {
				inputfiles[i] = argv[args];
			}
			--args;
		} else if(strcmp(argv[args], "-ipe") == 0) {
			++args;
			fileCount = fileCounter_PE;
			for(i = args; i < argc && strncmp(argv[i], "-", 1) != 0; ++i) {
				++fileCounter_PE;
			}
			if(fileCounter_PE % 2) {
				fprintf(stderr, "Uneven number of paired end files.\n");
				exit(1);
			} else if(fileCounter_PE == 0) {
				fprintf(stderr, "No paired end files were specified.\n");
				exit(1);
			} else {
				inputfiles_PE = realloc(inputfiles_PE, fileCounter_PE * sizeof(char *));
				if(!inputfiles_PE) {
					ERROR();
				}
			}
			
			for(i = fileCount; i < fileCounter_PE; ++i, ++args) {
				inputfiles_PE[i] = argv[args];
			}
			--args;
		} else if(strcmp(argv[args], "-int") == 0) {
			++args;
			fileCount = fileCounter_INT;
			for(i = args; i < argc && (strncmp(argv[i], "-", 1) != 0 || strcmp(argv[i], "--") == 0); ++i) {
				++fileCounter_INT;
			}
			if(fileCounter_INT == 0) {
				fprintf(stderr, "No interleaved files were specified.\n");
				exit(1);
			}
			inputfiles_INT = realloc(inputfiles_INT, fileCounter_INT * sizeof(char *));
			if(!inputfiles_INT) {
				ERROR();
			}
			for(i = fileCount; i < fileCounter_INT; ++i, ++args) {
				inputfiles_INT[i] = argv[args];
			}
			--args;
		} else if(strcmp(argv[args], "-o") == 0) {
			++args;
			if(args < argc) {
				outputfilename = smalloc(strlen(argv[args]) + 64);
				strcpy(outputfilename, argv[args]);
			}
		} else if(strcmp(argv[args], "-qc") == 0) {
			if(qcreport) {
				qcreport->verbose++;
			} else {
				qcreport = init_QCstat(0);
				if(!qcreport) {
					ERROR();
				}
			}
		} else if(strcmp(argv[args], "-ml") == 0) {
			++args;
			if(args < argc) {
				minlen = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid minimum length parsed\n");
					exit(1);
				}
			}
		} else if(strcmp(argv[args], "-xl") == 0) {
			++args;
			if(args < argc) {
				maxlen = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid minimum length parsed\n");
					exit(1);
				}
			}
		} else if(strcmp(argv[args], "-mp") == 0) {
			++args;
			if(args < argc) {
				minPhred = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid minimum phred score parsed\n");
					exit(1);
				}
			}
		} else if(strcmp(argv[args], "-mi") == 0) {
			++args;
			if(args < argc) {
				minmaskQ = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid internal minimum phred score parsed\n");
					exit(1);
				}
			}
		} else if(strcmp(argv[args], "-eq") == 0) {
			++args;
			if(args < argc) {
				minQ = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid average quality score parsed\n");
					exit(1);
				}
			}
		} else if(strcmp(argv[args], "-5p") == 0) {
			++args;
			if(args < argc) {
				fiveClip = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-5p\".\n");
					exit(1);
				}
			}
		} else if(strcmp(argv[args], "-3p") == 0) {
			++args;
			if(args < argc) {
				threeClip = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-3p\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-h") == 0) {
			return helpMessage(stdout);
		} else {
			fprintf(stderr, " Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, " Printing help message:\n");
			return helpMessage(stderr);
		}
		++args;
	}
	
	/* set ptrs */
	if(outputfilename) {
		i = strlen(outputfilename);
		sprintf(outputfilename + i, ".fq");
		out = sfopen(outputfilename, "wb");
		if(fileCounter_PE + fileCounter_INT) {
			sprintf(outputfilename + i, "_int.fq");
			out_int = sfopen(outputfilename, "wb");
		}
		if(qcreport) {
			sprintf(outputfilename + i, ".json");
			out_json = sfopen(outputfilename, "wb");
		}
		outputfilename[i] = 0;
	}
	printFsa_ptr = &printTrimFsa;
	printFsa_ptr(0, 0, 0, 0, out);
	printFsa_pair_ptr = &printTrimFsa_pair;
	printFsa_pair_ptr(0, 0, 0, 0, 0, 0, 0, out_int);
	
	/* set to2Bit conversion */
	to2Bit = smalloc(384); /* 128 * 3 = 384 -> OS independent */
	for(i = 0; i < 384; ++i) {
		to2Bit[i] = 8;
	}
	to2Bit += 128;
	to2Bit['\n'] = 16;
	to2Bit['A'] = 0;
	to2Bit['C'] = 1;
	to2Bit['G'] = 2;
	to2Bit['T'] = 3;
	to2Bit['N'] = 4;
	to2Bit['a'] = 0;
	to2Bit['c'] = 1;
	to2Bit['g'] = 2;
	to2Bit['t'] = 3;
	to2Bit['n'] = 4;
	to2Bit['R'] = 0;
	to2Bit['Y'] = 1;
	to2Bit['S'] = 2;
	to2Bit['W'] = 3;
	to2Bit['K'] = 2;
	to2Bit['M'] = 0;
	to2Bit['B'] = 1;
	to2Bit['D'] = 0;
	to2Bit['H'] = 3;
	to2Bit['V'] = 2;
	to2Bit['X'] = 4;
	to2Bit['r'] = 0;
	to2Bit['y'] = 1;
	to2Bit['s'] = 2;
	to2Bit['w'] = 3;
	to2Bit['k'] = 2;
	to2Bit['m'] = 0;
	to2Bit['b'] = 1;
	to2Bit['d'] = 0;
	to2Bit['h'] = 3;
	to2Bit['v'] = 2;
	to2Bit['x'] = 4;
	to2Bit['U'] = 3;
	to2Bit['u'] = 3;
	
	/* trim sequences */
	totFrags = 0;
	if(!fileCounter && !fileCounter_PE && !fileCounter_INT) {
		fileCounter = 1;
		inputfiles = smalloc(sizeof(char *));
		*inputfiles = "--";
	}
	if(minPhred < minmaskQ) {
		minPhred = minmaskQ;
	}
	
	/* SE */
	if(fileCounter > 0) {
		totFrags += run_input(inputfiles, fileCounter, minPhred, minmaskQ, minQ, fiveClip, threeClip, minlen, maxlen, to2Bit, prob, qcreport, out);
	}
	
	/* PE */
	if(fileCounter_PE > 0) {
		totFrags += run_input_PE(inputfiles_PE, fileCounter_PE, minPhred, minmaskQ, minQ, fiveClip, threeClip, minlen, maxlen, to2Bit, prob, qcreport, out);
	}
	
	/* INT */
	if(fileCounter_INT > 0) {
		totFrags += run_input_INT(inputfiles_INT, fileCounter_INT, minPhred, minmaskQ, minQ, fiveClip, threeClip, minlen, maxlen, to2Bit, prob, qcreport, out);
	}
	
	/* print QC */
	if(qcreport) {
		print_QCstat(qcreport, minQ, minPhred, minmaskQ, minlen, maxlen, fiveClip, threeClip, out_json);
		destroy_QCstat(qcreport);
	}
	
	/* clean up */
	free(outputfilename);
	if(out != stdout) {
		fclose(out);
	}
	if(out_int != stdout) {
		fclose(out_int);
	}
	if(out_json != stderr) {
		fclose(out_json);
	}
	
	fprintf(stderr, "#\n# Total number of query fragment after trimming:\t%lu\n", totFrags);
	
	return 0;
}
