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
#ifndef _WIN32
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/types.h>
#else
typedef int key_t;
#define ftok(charPtr, integer) (0)
#define shmget(key, size, permission) ((size != 0) ? (-1) : (-key))
#define shmat(shmid, NULL_Ptr, integer) (NULL)
#define shmdt(dest) fprintf(stderr, "sysV not available on Windows.\n")
#define shmctl(shmid, cmd, buf) fprintf(stderr, "sysV not available on Windows.\n")
#endif
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "pherror.h"
#include "hashmapkma.h"
#include "shm.h"
#include "stdnuc.h"
#include "version.h"

void hashMap_shm_detach(HashMapKMA *dest) {
	shmdt(dest->exist);
	shmdt(dest->values);
	shmdt(dest->key_index);
	shmdt(dest->value_index);
}

int hashMapKMA_setupSHM(HashMapKMA *dest, FILE *file, const char *filename) {
	
	int shmid, mlen, status;
	unsigned DB_size;
	long unsigned mask, size;
	key_t key;
	
	/* load sizes */
	sfread(&DB_size, sizeof(unsigned), 1, file);
	sfread(&dest->mlen, sizeof(unsigned), 1, file);
	sfread(&dest->prefix_len, sizeof(unsigned), 1, file);
	sfread(&dest->prefix, sizeof(long unsigned), 1, file);
	sfread(&dest->size, sizeof(long unsigned), 1, file);
	sfread(&dest->n, sizeof(long unsigned), 1, file);
	sfread(&dest->v_index, sizeof(long unsigned), 1, file);
	sfread(&dest->null_index, sizeof(long unsigned), 1, file);
	mlen = dest->mlen;
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (mlen << 1));
	status = 0;
	
	/* check shared memory, else load */
	size = dest->size;
	if((dest->size - 1) == mask) {
		if(dest->v_index <= UINT_MAX) {
			size *= sizeof(unsigned);
		} else {
			size *= sizeof(long unsigned);
		}
	} else {
		if(dest->n <= UINT_MAX) {
			size *= sizeof(unsigned);
		} else {
			size *= sizeof(long unsigned);
		}
	}
	key = ftok(filename, 'e');
	shmid = shmget(key, size, IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap e\n");
		sfseek(file, size, SEEK_CUR);
		dest->exist = 0;
		status = 1;
	} else {
		dest->exist = shmat(shmid, NULL, 0);
		cfread(dest->exist, 1, size, file);
	}
	
	/* values */
	size = dest->v_index;
	if(DB_size < USHRT_MAX) {
		size *= sizeof(short unsigned);
	} else {
		size *= sizeof(unsigned);
	}
	key = ftok(filename, 'v');
	shmid = shmget(key, size, IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap v\n");
		sfseek(file, size, SEEK_CUR);
		dest->values = 0;
		status = 1;
	} else {
		/* found */
		dest->values = shmat(shmid, NULL, 0);
		cfread(dest->values, 1, size, file);
	}
	if((dest->size - 1) == mask) {
		return status;
	}
	
	/* kmers */
	size = dest->n + 1;
	if(dest->mlen <= 16) {
		size *= sizeof(unsigned);
	} else {
		size *= sizeof(long unsigned);
	}
	key = ftok(filename, 'k');
	shmid = shmget(key, size, IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap k\n");
		sfseek(file, size, SEEK_CUR);
		dest->values = 0;
		status = 1;
	} else {
		/* found */
		dest->key_index = shmat(shmid, NULL, 0);
		cfread(dest->key_index, 1, size, file);
	}
	
	/* value indexes */
	size = dest->n;
	if(dest->v_index < UINT_MAX) {
		size *= sizeof(unsigned);
	} else {
		size *= sizeof(long unsigned);
	}
	key = ftok(filename, 'i');
	shmid = shmget(key, size, IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap i\n");
		sfseek(file, size, SEEK_CUR);
		dest->value_index = 0;
		status = 1;
	} else {
		/* found */
		dest->value_index = shmat(shmid, NULL, 0);
		cfread(dest->value_index, 1, size, file);
	}
	
	if(fread(&dest->kmersize, sizeof(unsigned), 1, file)) {
		sfread(&dest->flag, sizeof(unsigned), 1, file);
		setCmerPointers(dest->flag);
	} else {
		dest->kmersize = dest->mlen;
		dest->flag = 0;
	}
	
	return status;
}

void hashMapKMA_destroySHM(HashMapKMA *dest, FILE *file, const char *filename) {
	
	int shmid, mlen;
	unsigned DB_size;
	long unsigned mask, size;
	key_t key;
	
	/* load sizes */
	sfread(&DB_size, sizeof(unsigned), 1, file);
	sfread(&dest->mlen, sizeof(unsigned), 1, file);
	sfread(&dest->prefix_len, sizeof(unsigned), 1, file);
	sfread(&dest->prefix, sizeof(long unsigned), 1, file);
	sfread(&dest->size, sizeof(long unsigned), 1, file);
	sfread(&dest->n, sizeof(long unsigned), 1, file);
	sfread(&dest->v_index, sizeof(long unsigned), 1, file);
	sfread(&dest->null_index, sizeof(long unsigned), 1, file);
	mlen = dest->mlen;
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (mlen << 1));
	
	/* check shared memory, and destroy */
	size = dest->size;
	if((dest->size - 1) == mask) {
		if(dest->v_index <= UINT_MAX) {
			size *= sizeof(unsigned);
		} else {
			size *= sizeof(long unsigned);
		}
	} else {
		if(dest->n <= UINT_MAX) {
			size *= sizeof(unsigned);
		} else {
			size *= sizeof(long unsigned);
		}
	}
	key = ftok(filename, 'e');
	shmid = shmget(key, size, 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
		sfseek(file, size, SEEK_CUR);
	}
	
	/* values */
	size = dest->v_index;
	if(DB_size < USHRT_MAX) {
		size *= sizeof(short unsigned);
	} else {
		size *= sizeof(unsigned);
	}
	key = ftok(filename, 'v');
	shmid = shmget(key, size, 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
		sfseek(file, size, SEEK_CUR);
	}
	
	/* kmers */
	size = dest->n + 1;
	if(dest->mlen <= 16) {
		size *= sizeof(unsigned);
	} else {
		size *= sizeof(long unsigned);
	}
	key = ftok(filename, 'k');
	shmid = shmget(key, size, 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
		sfseek(file, size, SEEK_CUR);
	}
	
	/* value indexes */
	size = dest->n;
	if(dest->v_index < UINT_MAX) {
		size *= sizeof(unsigned);
	} else {
		size *= sizeof(long unsigned);
	}
	key = ftok(filename, 'i');
	shmid = shmget(key, size, 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
		sfseek(file, size, SEEK_CUR);
	}
	
	if(fread(&dest->kmersize, sizeof(unsigned), 1, file)) {
		sfread(&dest->flag, sizeof(unsigned), 1, file);
		setCmerPointers(dest->flag);
	} else {
		dest->kmersize = dest->mlen;
		dest->flag = 0;
	}
}

int * length_setupSHM(FILE *file, const char *filename) {
	
	int shmid, *template_lengths;
	long unsigned size;
	key_t key;
	
	/* load size */
	sfseek(file, 0, SEEK_END);
	size = ftell(file) - sizeof(int);
	sfseek(file, sizeof(int), SEEK_SET);
	
	key = ftok(filename, 'l');
	shmid = shmget(key, size, IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared length\n");
		template_lengths = 0;
	} else {
		template_lengths = shmat(shmid, NULL, 0);
		cfread(template_lengths, sizeof(unsigned), size / sizeof(unsigned), file);
	}
	
	return template_lengths;
}

void length_destroySHM(FILE *file, const char *filename) {
	
	int shmid;
	long unsigned size;
	key_t key;
	
	/* load size */
	sfseek(file, 0, SEEK_END);
	size = ftell(file) - sizeof(int);
	
	key = ftok(filename, 'l');
	shmid = shmget(key, size, 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
	}
	
}

long unsigned * seq_setupSHM(FILE *file, const char *filename) {
	
	int shmid;
	long unsigned size, *seq;
	key_t key;
	
	/* load size */
	sfseek(file, 0, SEEK_END);
	size = ftell(file);
	rewind(file);
	
	key = ftok(filename, 's');
	shmid = shmget(key, size, IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared length\n");
		seq = 0;
	} else {
		seq = shmat(shmid, NULL, 0);
		cfread(seq, sizeof(long unsigned), size / sizeof(long unsigned), file);
	}
	
	return seq;
}

void seq_destroySHM(FILE *file, const char *filename) {
	
	int shmid;
	long unsigned size;
	key_t key;
	
	/* load size */
	sfseek(file, 0, SEEK_END);
	size = ftell(file);
	
	key = ftok(filename, 's');
	shmid = shmget(key, size, 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
	}
}

char * name_setupSHM(FILE *file, const char *filename) {
	
	int i, shmid;
	long unsigned size;
	char *template_names;
	key_t key;
	
	/* load size */
	sfseek(file, 0, SEEK_END);
	size = ftell(file);
	rewind(file);
	
	key = ftok(filename, 'n');
	shmid = shmget(key, size, IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared length\n");
		template_names = 0;
	} else {
		template_names = shmat(shmid, NULL, 0);
		cfread(template_names, 1, size, file);
		for(i = 0; i < size; ++i) {
			if(template_names[i] == '\n') {
				template_names[i] = 0;
			}
		}
	}
	
	return template_names;
}

void name_destroySHM(FILE *file, const char *filename) {
	
	int shmid;
	long unsigned size;
	key_t key;
	
	/* load size */
	sfseek(file, 0, SEEK_END);
	size = ftell(file);
	
	key = ftok(filename, 'n');
	shmid = shmget(key, size, 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
	}
}

static void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# kma_shm sets up a shared database (sysV) for mapping with KMA.\n");
	fprintf(helpOut, "# Options are:\t\tDesc:\t\t\t\tDefault:\tRequirements:\n");
	fprintf(helpOut, "#\n");
	fprintf(helpOut, "#\t-t_db\t\tTemplate DB\t\t\tNone\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-destroy\tDestroy shared DB\t\tFalse\n");
	fprintf(helpOut, "#\t-shmLvl\t\tLevel of shared memory\t\t1\n");
	fprintf(helpOut, "#\t-shm-h\t\tExplain shm levels\n");
	fprintf(helpOut, "#\t-v\t\tVersion\n");
	fprintf(helpOut, "#\t-h\t\tShows this help message\n");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int shm_main(int argc, char *argv[]) {
	
	int args, file_len, destroy, status, *template_lengths;
	unsigned shmLvl;
	long unsigned *seq;
	char *templatefilename, *template_names;
	HashMapKMA *templates;
	FILE *file;
	
	#ifdef _WIN32
	fprintf(stderr, "sysV not available on Windows.\n");
	return 1;
	#endif
	
	/* SET DEFAULTS */
	templatefilename = 0;
	destroy = 0;
	shmLvl = 1;
	status = 0;
	
	/* PARSE COMMAND LINE OPTIONS */
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-t_db") == 0) {
			++args;
			if(args < argc) {
				templatefilename = malloc(strlen(argv[args]) + 64);
				if(!templatefilename) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				strcpy(templatefilename, argv[args]);
			}
		} else if(strcmp(argv[args], "-destroy") == 0) {
			destroy = 1;
		} else if(strcmp(argv[args], "-shmLvl") == 0) {
			++args;
			if(args < argc) {
				shmLvl = atoi(argv[args]);
				if(!shmLvl) {
					fprintf(stderr, "Invalid shmLvl\n");
					exit(0);
				}
			}
		} else if(strcmp(argv[args], "-v") == 0) {
			fprintf(stdout, "KMA_SHM-%s\n", KMA_VERSION);
			exit(0);
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else if(strcmp(argv[args], "-shm-h") == 0) {
			fprintf(stderr, "# Flags for shared memory, add them to combine them\n");
			fprintf(stderr, "# After shm is setup, the DB-files should not be changed.\n");
			fprintf(stderr, "#\n");
			fprintf(stderr, "#\tDB piece\t\tFlag\n");
			fprintf(stderr, "#\t*.comp.b\t\t1\n");
			fprintf(stderr, "#\t*.decon.comp.b\t\t2\n");
			fprintf(stderr, "#\t*.length.b\t\t4\n");
			fprintf(stderr, "#\t*.seq.b\t\t\t8\n");
			fprintf(stderr, "#\t*.name\t\t\t16\n");
			//fprintf(stderr, "#\tall\t\t\t31\n");
			exit(0);
		} else {
			fprintf(stderr, "# Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, "# Printing help message:\n");
			helpMessage(-1);
		}
		++args;
	}
	if(templatefilename == 0) {
		fprintf(stderr, "# Too few arguments handed\n");
		fprintf(stderr, "# Printing help message:\n");
		helpMessage(-1);
	}
	
	file_len = strlen(templatefilename);
	templates = malloc(sizeof(HashMapKMA));
	if(!templates) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	/* setup or destroy shm */
	if(destroy) {
		/* *comp.b */
		if(shmLvl & 1) {
			strcat(templatefilename, ".comp.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				status |= errno;
			} else {
				hashMapKMA_destroySHM(templates, file, templatefilename);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *decon.comp.b */
		if(shmLvl & 2) {
			strcat(templatefilename, ".decon.comp.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				status |= errno;
			} else {
				hashMapKMA_destroySHM(templates, file, templatefilename);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *.length.b */
		if(shmLvl & 4) {
			strcat(templatefilename, ".length.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				status |= errno;
			} else {
				length_destroySHM(file, templatefilename);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *.seq.b */
		if(shmLvl & 8) {
			/* *.seq.b */
			strcat(templatefilename, ".seq.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				status |= errno;
			} else {
				seq_destroySHM(file, templatefilename);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *.name */
		if(shmLvl & 16) {
			strcat(templatefilename, ".name");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				status |= errno;
			} else {
				name_destroySHM(file, templatefilename);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
	} else {
		/* *.comp.b */
		if(shmLvl & 1) {
			strcat(templatefilename, ".comp.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				status |= errno;
			} else {
				status |= hashMapKMA_setupSHM(templates, file, templatefilename);
				hashMap_shm_detach(templates);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *.decon.comp.b */
		if(shmLvl & 2) {
			strcat(templatefilename, ".decon.comp.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				status |= errno;
			} else {
				status |=  hashMapKMA_setupSHM(templates, file, templatefilename);
				hashMap_shm_detach(templates);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *.length.b */
		if(shmLvl & 4) {
			strcat(templatefilename, ".length.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				status |= errno;
			} else {
				template_lengths = length_setupSHM(file, templatefilename);
				if(template_lengths) {
					shmdt(template_lengths);
				} else {
					status |= 1;
				}
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *.seq.b*/
		if(shmLvl & 8) {
			/* *.seq.b */
			strcat(templatefilename, ".seq.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				status |= errno;
			} else {
				seq = seq_setupSHM(file, templatefilename);
				if(seq) {
					shmdt(seq);
				} else {
					status |= 1;
				}
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *.name */
		if(shmLvl & 16) {
			strcat(templatefilename, ".name");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				status |= errno;
			} else {
				template_names = name_setupSHM(file, templatefilename);
				if(template_names) {
					shmdt(template_names);
				} else {
					status |= 1;
				}
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
	}
	
	/* Set sys shm to 1 GB
	* sudo sysctl -w kern.sysv.shmmax=1073741824
	* sudo sysctl -w kern.sysv.shmall=1073741824
	*
	* Set sys shm to 2 GB
	* sudo sysctl -w kern.sysv.shmmax=2147483648
	* sudo sysctl -w kern.sysv.shmall=2147483648
	*
	* Set sys shm to 3 GB
	* sudo sysctl -w kern.sysv.shmmax=3221225472
	* sudo sysctl -w kern.sysv.shmall=3221225472
	*
	* check status:
	* ipcs -a
	*/
	
	return status;
}
