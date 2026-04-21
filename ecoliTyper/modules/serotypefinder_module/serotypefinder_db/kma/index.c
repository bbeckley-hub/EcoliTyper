/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#define _XOPEN_SOURCE 600
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "compress.h"
#include "decon.h"
#include "hashmap.h"
#include "hashmapkma.h"
#include "index.h"
#include "loadupdate.h"
#include "makeindex.h"
#include "pherror.h"
#include "qualcheck.h"
#include "stdnuc.h"
#include "stdstat.h"
#include "updateindex.h"
#include "valueshash.h"
#include "version.h"

static void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# kma_index creates the databases needed to run KMA, from a list of fasta files given.\n");
	fprintf(helpOut, "# Options are:\t\tDesc:\t\t\t\t\tDefault:\n");
	fprintf(helpOut, "#\n");
	fprintf(helpOut, "#\t-i\t\tInput/query file name (STDIN: \"--\")\tNone\n");
	fprintf(helpOut, "#\t-o\t\tOutput file\t\t\t\tInput/template file\n");
	fprintf(helpOut, "#\t-batch\t\tBatch input file\n");
	fprintf(helpOut, "#\t-deCon\t\tFile with contamination (STDIN: \"--\")\tNone/False\n");
	fprintf(helpOut, "#\t-batchD\t\tBatch decon file\n");
	fprintf(helpOut, "#\t-t_db\t\tAdd to existing DB\t\t\tNone/False\n");
	fprintf(helpOut, "#\t-k\t\tKmersize\t\t\t\t16\n");
	fprintf(helpOut, "#\t-m\t\tMinimizer size\t\t\t\t16/False\n");
	fprintf(helpOut, "#\t-hc\t\tHomopolymer compression\t\t\t\tFalse\n");
	fprintf(helpOut, "#\t-c\t\tTemplates are open reading frames\n");
	fprintf(helpOut, "#\t-C\t\tTemplates are CDSs\n");
	fprintf(helpOut, "#\t-ML\t\tMinimum length of templates\t\tkmersize (16)\n");
	fprintf(helpOut, "#\t-CS\t\tStart Chain size\t\t\t1 M\n");
	fprintf(helpOut, "#\t-ME\t\tMega DB\t\t\t\t\tFalse\n");
	fprintf(helpOut, "#\t-NI\t\tDo not dump *.index.b\t\t\tFalse\n");
	fprintf(helpOut, "#\t-Sparse\t\tMake Sparse DB ('-' for no prefix)\tNone/False\n");
	fprintf(helpOut, "#\t-ht\t\tHomology template\t\t\t1.0\n");
	fprintf(helpOut, "#\t-hq\t\tHomology query\t\t\t\t1.0\n");
	fprintf(helpOut, "#\t-and\t\tBoth homolgy thresholds\n#\t\t\thas to be reached\t\t\tor\n");
	fprintf(helpOut, "#\t-nbp\t\tNo bias print\t\t\t\tFalse\n");
	fprintf(helpOut, "#\t-v\t\tVersion\n");
	fprintf(helpOut, "#\t-h\t\tShows this help message\n");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int index_main(int argc, char *argv[]) {
	
	int i, args, stop, filecount, deconcount, sparse_run, size, mapped_cont;
	int file_len, appender, prefix_len, MinLen, MinKlen;
	unsigned kmersize, mlen, flag, kmerindex, megaDB, **Values;
	unsigned *template_lengths, *template_slengths, *template_ulengths;
	long unsigned initialSize, prefix, mask;
	double homQ, homT;
	char **inputfiles, *outputfilename, *templatefilename, **deconfiles;
	char *to2Bit, *line, *exeBasic;
	unsigned char *update;
	FILE *inputfile, *out;
	time_t t0, t1;
	HashMap *templates;
	HashMapKMA *finalDB;
	
	if(argc == 1) {
		fprintf(stderr, "# Too few arguments handed.\n");
		helpMessage(-1);
	} else if(sizeof(long unsigned) != 8) {
		fprintf(stderr, "Need a 64-bit system.\n");
		exit(1);
	}
	
	/* set defaults */
	initialSize = 1048576;
	templates = 0;
	kmersize = 16;
	kmerindex = 16;
	mlen = 0;
	flag = 0;
	sparse_run = 0;
	appender = 0;
	MinLen = 0;
	MinKlen = 1;
	prefix_len = 0;
	prefix = 0;
	homQ = 1;
	homT = 1;
	cmp = &cmp_or;
	template_ulengths = 0;
	template_slengths = 0;
	filecount = 0;
	deconcount = 0;
	outputfilename = 0;
	templatefilename = 0;
	megaDB = 0;
	inputfiles = smalloc(sizeof(char*));
	deconfiles = smalloc(sizeof(char*));
	to2Bit = smalloc(384);
	
	/* set to2Bit */
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
	
	/* Future RNA encoding */
	to2Bit['U'] = 3;
	to2Bit['u'] = 3;
	
	
	/* PARSE COMMAND LINE OPTIONS */
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-i") == 0) {
			stop = 0;
			++args;
			while(stop == 0 && args < argc) {
				if(strncmp(argv[args], "-", 1) != 0 || strcmp(argv[args], "--") == 0) {
					++filecount;
					inputfiles = realloc(inputfiles, filecount * sizeof(char*));
					if(inputfiles == NULL) {
						ERROR();
					}
					inputfiles[filecount - 1] = argv[args];
					++args;
				} else {
					stop = 1;
				}
			}
			--args;
		} else if(strcmp(argv[args], "-o") == 0) {
			++args;
			if(args < argc) {
				outputfilename = smalloc(strlen(argv[args]) + 64);
				strcpy(outputfilename, argv[args]);
			}
		} else if(strcmp(argv[args], "-deCon") == 0) {
			stop = 0;
			++args;
			while(stop == 0 && args < argc) {
				if(strncmp(argv[args], "-", 1) != 0 || strcmp(argv[args], "--") == 0) {
					deconfiles = realloc(deconfiles, (deconcount + 1) * sizeof(char*));
					if(deconfiles == NULL) {
						ERROR();
					}
					deconfiles[deconcount] = argv[args];
					++deconcount;
					++args;
				} else {
					stop = 1;
				}
			}
			if(deconcount == 0) {
				fprintf(stderr, "No deCon file specified.\n");
				exit(1);
			}
			--args;
		} else if(strcmp(argv[args], "-t_db") == 0) {
			++args;
			if(args < argc) {
				templatefilename = malloc(strlen(argv[args]) + 64);
				if(!templatefilename) {
					ERROR();
				}
				strcpy(templatefilename, argv[args]);
			}
		} else if(strcmp(argv[args], "-k") == 0) {
			++args;
			if(args < argc) {
				kmersize = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid kmersize parsed\n");
					exit(1);
				} else if(kmersize == 0) {
					fprintf(stderr, "# Invalid kmersize parsed, using default\n");
					kmersize = 16;
				} else if(kmersize > 31) {
					fprintf(stderr, "# Invalid kmersize parsed, max size is 31\n");
					exit(1);
				}
				kmerindex = kmersize;
			}
		} else if(strcmp(argv[args], "-m") == 0) {
			++args;
			if(args < argc) {
				mlen = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0 || mlen == 0) {
					fprintf(stderr, "# Invalid minimizer size parsed\n");
					exit(1);
				} else if(mlen > 30) {
					mlen = 30;
				}
			} else {
				fprintf(stderr, "# Missing minimizer size\n");
				exit(1);
			}
			flag |= 2;
		} else if(strcmp(argv[args], "-hc") == 0) {
			flag |= 1;
		} else if(strcmp(argv[args], "-k_t") == 0) {
			++args;
			if(args < argc) {
				kmersize = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid kmersize parsed\n");
					exit(1);
				} else if(kmersize == 0) {
					fprintf(stderr, "# Invalid kmersize parsed, using default\n");
					kmersize = 16;
				} else if(kmersize > 31) {
					fprintf(stderr, "# Invalid kmersize parsed, max size is 31\n");
					exit(1);
				}
			}
		} else if(strcmp(argv[args], "-k_i") == 0) {
			++args;
			if(args < argc) {
				kmerindex = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid kmersize parsed\n");
					exit(1);
				} else if(kmerindex == 0) {
					fprintf(stderr, "# Invalid kmersize parsed, using default\n");
					kmerindex = 16;
				} else if(kmerindex > 31) {
					fprintf(stderr, "# Invalid kmersize parsed, max size is 31\n");
					exit(1);
				}
			}
		} else if(strcmp(argv[args], "-CS") == 0) {
			++args;
			if(args < argc) {
				
				size = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid start size parsed\n");
					exit(1);
				}
				initialSize = pow(2, ceil(log(size)/log(2))) + 0.5;
				initialSize *= 1048576;
				if(initialSize == 0) {
					fprintf(stderr, "# Invalid Chain Size parsed, using default\n");
					initialSize = 1048576;
				}
			}
		} else if(strcmp(argv[args], "-and") == 0) {
			cmp = &cmp_and;
		} else if(strcmp(argv[args], "-c") == 0) {
			qualcheck = &qualCheck;
		} else if(strcmp(argv[args], "-C") == 0) {
			qualcheck = &internalStopCheck;
		} else if(strcmp(argv[args], "-ML") == 0) {
			++args;
			if(args < argc) {
				MinLen = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid minimum length parsed\n");
					exit(1);
				} else if(MinLen <= 0) {
					fprintf(stderr, "# Invalid minimum length parsed, using default\n");
					MinLen = 0;
				}
			}
		} else if(strcmp(argv[args], "-hq") == 0) {
			++args;
			if(args < argc) {
				homQ = strtod(argv[args], &exeBasic);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-hq\".\n");
					exit(1);
				} else if(homQ < 0) {
					fprintf(stderr, "Invalid -hq\n");
					homQ = 1.0;
				}
			}
		} else if(strcmp(argv[args], "-ht") == 0) {
			++args;
			if(args < argc) {
				homT = strtod(argv[args], &exeBasic);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-ht\".\n");
					exit(1);
				} else if(homT < 0) {
					fprintf(stderr, "Invalid -hq\n");
					homT = 1.0;
				}
			}
		} else if(strcmp(argv[args], "-batch") == 0) {
			++args;
			if(args < argc) {
				inputfile = sfopen(argv[args], "rb");
				sfseek(inputfile, 0, SEEK_END);
				size = ftell(inputfile) + 1;
				rewind(inputfile);
				
				++filecount;
				inputfiles = realloc(inputfiles, filecount * sizeof(char*));
				if(!inputfiles) {
					ERROR();
				}
				inputfiles[filecount - 1] = malloc(size);
				if(!inputfiles[filecount - 1]) {
					ERROR();
				}
				
				/* get number of file */
				sfread(inputfiles[filecount - 1], 1, size - 1, inputfile);
				fclose(inputfile);
				
				i = size;
				size = 0;
				line = inputfiles[filecount - 1];
				while(--i) {
					if(line[i] == '\n') {
						size++;
						while(isspace(line[i])) {
							line[i] = 0;
							--i;
						}
					}
				}
				--size;
				inputfiles = realloc(inputfiles, (filecount + size) * sizeof(char*));
				if(!inputfiles) {
					ERROR();
				}
				for(i = size; i; --i) {
					while(*line != 0) {
						++line;
					}
					while(*line == 0) {
						++line;
					}
					inputfiles[filecount] = line;
					++filecount;
				}
			}
		} else if(strcmp(argv[args], "-batchD") == 0) {
			++args;
			if(args < argc) {
				inputfile = sfopen(argv[args], "rb");
				sfseek(inputfile, 0, SEEK_END);
				size = ftell(inputfile) + 1;
				rewind(inputfile);
				
				++deconcount;
				deconfiles = realloc(deconfiles, deconcount * sizeof(char*));
				if(!inputfiles) {
					ERROR();
				}
				deconfiles[deconcount - 1] = malloc(size);
				if(!deconfiles[deconcount - 1]) {
					ERROR();
				}
				
				/* get number of file */
				sfread(deconfiles[deconcount - 1], 1, size - 1, inputfile);
				fclose(inputfile);
				
				i = size;
				size = 0;
				line = deconfiles[deconcount - 1];
				while(--i) {
					if(line[i] == '\n') {
						size++;
						while(isspace(line[i])) {
							line[i] = 0;
							--i;
						}
					}
				}
				--size;
				deconfiles = realloc(deconfiles, (deconcount + size) * sizeof(char*));
				if(!inputfiles) {
					ERROR();
				}
				for(i = size; i; --i) {
					while(*line != 0) {
						++line;
					}
					while(*line == 0) {
						++line;
					}
					deconfiles[deconcount] = line;
					++deconcount;
				}
			}
		} else if(strcmp(argv[args], "-Sparse") == 0) {
			sparse_run = 1;
			++args;
			if(args < argc) {
				if(strcmp(argv[args], "-") == 0) {
					prefix_len = 0;
					prefix = 1;
				} else {
					prefix_len = strlen(argv[args]);
					prefix = 0;
					update = (unsigned char *) argv[args];
					for(i = 0; i < prefix_len; ++i) {
						prefix = (prefix << 2) | to2Bit[update[i]];
						if(to2Bit[update[i]] > 3) {
							fprintf(stderr, "Invalid prefix.\n");
							exit(1);
						}
					}
					if(prefix_len == 0) {
						fprintf(stderr, "Invalid prefix.\n");
						exit(1);
					}
				}
			}
		} else if(strcmp(argv[args], "-ME") == 0) {
			megaDB = 1;
		} else if(strcmp(argv[args], "-NI") == 0) {
			
		} else if(strcmp(argv[args], "-nbp") == 0) {
			biasPrintPtr = &biasNoPrint;
		} else if(strcmp(argv[args], "-v") == 0) {
			fprintf(stdout, "KMA_index-%s\n", KMA_VERSION);
			exit(0);
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else {
			fprintf(stderr, "# Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, "# Printing help message:\n");
			helpMessage(1);
		}
		++args;
	}
	
	/* check for sufficient input */
	if(filecount == 0 && deconcount == 0) {
		fprintf(stderr, "No inputfiles defined.\n");
		helpMessage(-1);
	} else if(filecount == 0 && deconcount != 0 && templatefilename == 0) {
		fprintf(stderr, "Nothing to update.\n");
		exit(0);
	} else if(outputfilename == 0 && templatefilename != 0) {
		outputfilename = smalloc((strlen(templatefilename) + 64));
		strcpy(outputfilename, templatefilename);
	} else if(outputfilename == 0 && filecount != 0) {
		outputfilename = smalloc((strlen(*inputfiles) + 64));
		strcpy(outputfilename, *inputfiles);
	} else if(outputfilename == 0) {
		fprintf(stderr, "Output destination not defined.\n");
		helpMessage(-1);
	} else if(kmersize <= mlen) {
		fprintf(stderr, "Minimizer size cannot be smaller than kmersize.\n");
		exit(1);
	}
	file_len = strlen(outputfilename);
	
	if(!mlen) {
		mlen = kmersize;
	}
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (mlen << 1));
	if(megaDB) {
		initialSize = mask + 1;
	} else if(initialSize >= (mask + 1)) {
		initialSize = (mask + 1);
		megaDB = 1;
	}
	setCmerPointers(flag);
	
	/* load DB */
	if(templatefilename != 0) {
		/* load */
		fprintf(stderr, "# Loading database: %s\n", outputfilename);
		finalDB = smalloc(sizeof(HashMapKMA));
		kmerindex = load_DBs(templatefilename, outputfilename, &template_lengths, &template_ulengths, &template_slengths, finalDB);
		kmersize = finalDB->kmersize;
		mask = 0;
		mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
		prefix = finalDB->prefix;
		prefix_len = finalDB->prefix_len;
		
		/* determine params based on loaded DB */
		if(prefix_len == 0 && prefix == 0) {
			sparse_run = 0;
		} else {
			sparse_run = 1;
		}
		if(finalDB->mask == finalDB->size) {
			megaDB = 1;
			initialSize = finalDB->mask + 1;
		} else if(megaDB == 0) {
			initialSize = finalDB->size + 1;
		}
		appender = 1;
	} else {
		finalDB = 0;
		appender = 0;
	}
	
	/* function pointers */
	hashMap_add = &hashMap_addKMA;
	hashMapGet = &hashMapGetValue;
	addCont = &hashMap_addCont;
	
	if(finalDB == 0 || finalDB->DB_size < USHRT_MAX) {
		updateValuePtr = &updateShortValue;
		valuesKeyPtr = &huValuesKey;
		cmpValuesPtr = &cmpHuValues;
		valuesSize = &huSize;
	} else {
		updateValuePtr = &updateValue;
		valuesKeyPtr = &valuesKey;
		cmpValuesPtr = &cmpValues;
		valuesSize = &uSize;
	}
	
	if(sparse_run) {
		update_DB = &updateDBs_sparse;
		updateAnnotsPtr = &updateAnnots_sparse;
		if(prefix_len == 0) {
			prefix = 1;
		}
	} else {
		update_DB = &updateDBs;
		updateAnnotsPtr = &updateAnnots;
	}
	
	if(prefix_len != 0) {
		deConNode_ptr = &deConNode_sparse;
	} else {
		deConNode_ptr = &deConNode;
	}
	if(megaDB) {
		hashMap_add = &megaMap_addKMA;
		hashMapGet = &megaMap_getValue;
		addCont = &megaMap_addCont;
	}
	
	/* set homology check */
	if(MinLen > (kmersize + prefix_len + 1)) {
		MinKlen = 2 * (MinLen - kmersize - prefix_len + 1);
		for(i = 0; i < prefix_len; ++i) {
			MinKlen /= 4;
		}
	} else {
		MinLen = MAX(kmersize, kmerindex);
	}
	if(homT < 1) {
		QualCheck = &templateCheck;
	} else if(homQ < 1) {
		QualCheck = &queryCheck;
	} else {
		QualCheck = &lengthCheck;
	}
	
	/* update DBs */
	if(filecount != 0) {
		if(finalDB) {
			/* convert */
			templates = hashMapKMA_openChains(finalDB);
		} else {
			/* create */
			templates = hashMap_initialize(initialSize, kmersize, mlen, flag);
			template_lengths = smalloc(1024 * sizeof(unsigned));;
			if(sparse_run) {
				templates->prefix = prefix;
				templates->prefix_len = prefix_len;
				template_slengths = smalloc(1024 * sizeof(unsigned));
				template_ulengths = smalloc(1024 * sizeof(unsigned));
				*template_lengths = kmerindex;
				*template_slengths = 1024;
				*template_ulengths = 1024;
			} else {
				template_slengths = 0;
				template_ulengths = 0;
				*template_lengths = 1024;
			}
		}
		fprintf(stderr, "# Indexing databases.\n");
		t0 = clock();
		makeDB(templates, kmerindex, inputfiles, filecount, outputfilename, appender, to2Bit, MinLen, MinKlen, homQ, homT, &template_lengths, &template_ulengths, &template_slengths);
		t1 = clock();
		fprintf(stderr, "#\n# Total time used for DB indexing: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
		free(template_lengths);
		free(template_slengths);
		free(template_ulengths);
		
		/* compress db */
		fprintf(stderr, "# Compressing templates\n");
		t0 = clock();
		if(templates->DB_size < USHRT_MAX) {
			valuesKeyPtr = &huValuesKey;
			cmpValuesPtr = &cmpHuValues;
			valuesSize = &huSize;
		} else {
			valuesKeyPtr = &valuesKey;
			cmpValuesPtr = &cmpValues;
			valuesSize = &uSize;
		}
		strcat(outputfilename, ".comp.b");
		out = sfopen(outputfilename, "wb+");
		if(templates->table != 0) {
			finalDB = compressKMA_DB(templates, out);
		} else {
			finalDB = compressKMA_megaDB(templates, out);
		}
		fclose(out);
		outputfilename[file_len] = 0;
		free(templates);
		fprintf(stderr, "# Template database created.\n");
		t1 = clock();
		fprintf(stderr, "#\n# Total time used for DB compression: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	} else {
		++finalDB->size;
	}
	
	/* decontaminate */
	if(deconcount != 0) {
		/* open values */
		fprintf(stderr, "# Opening values\n");
		if(!finalDB) {
			strcat(outputfilename, ".comp.b");
			out = sfopen(outputfilename, "rb");
			finalDB = smalloc(sizeof(HashMapKMA));
			if(hashMapKMAload(finalDB, out)) {
				fprintf(stderr, "Wrong format of DB\n");
				exit(1);
			}
			fclose(out);
			outputfilename[file_len] = 0;
		}
		Values = hashMapKMA_openValues(finalDB);
		
		/* get decontamination info */
		fprintf(stderr, "# Adding decontamination information\n");
		t0 = clock();
		mapped_cont = deConDB(finalDB, deconfiles, deconcount, to2Bit, Values);
		fprintf(stderr, "# Contamination information added.\n");
		fprintf(stderr, "# %d kmers mapped to the DB.\n", mapped_cont);
		fprintf(stderr, "# Contamination mapped to %f %% of the DB.\n", 100.0 * mapped_cont / finalDB->n);
		
		/* compress DB */
		fprintf(stderr, "# Compressing templates\n");
		if(finalDB->DB_size < USHRT_MAX) {
			valuesKeyPtr = &huValuesKey;
			cmpValuesPtr = &cmpHuValues;
			valuesSize = &huSize;
		} else {
			valuesKeyPtr = &valuesKey;
			cmpValuesPtr = &cmpValues;
			valuesSize = &uSize;
		}
		if((finalDB->size - 1) != finalDB->mask) {
			compressKMA_deconDB(finalDB, Values);
		} else {
			compressKMA_deconMegaDB(finalDB, Values);
		}
		t1 = clock();
		fprintf(stderr, "#\n# Total time used for DB decontamination: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
		
		/* dump DB */
		fprintf(stderr, "# Dumping DB.\n");
		strcat(outputfilename, ".decon.comp.b");
		out = sfopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
		if((finalDB->size - 1) != mask) {
			hashMapKMA_dump(finalDB, out);
		} else {
			megaMapKMA_dump(finalDB, out);
		}
		fclose(out);
		outputfilename[file_len] = 0;
	}
	
	return 0;
}
