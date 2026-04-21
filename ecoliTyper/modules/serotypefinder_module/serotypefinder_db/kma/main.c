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
#include <string.h>
#include "cmp.h"
#include "db.h"
#include "dist.h"
#include "kma.h"
#include "index.h"
#include "merge.h"
#include "shm.h"
#include "seq2fasta.h"
#include "trim.h"
#include "update.h"

static int helpmessage(FILE *out) {
	
	fprintf(out, "# KMA enables alignment towards databases, using two k-mer mapping steps and one alignment step.\n");
	fprintf(out, "# %16s\t%-32s\n", "Options are:", "Desc:");
	fprintf(out, "# %16s\t%-32s\n", "", "Alignment and mapping");
	fprintf(out, "# %16s\t%-32s\n", "index", "Indexing of databases");
	fprintf(out, "# %16s\t%-32s\n", "shm", "Shared memory");
	fprintf(out, "# %16s\t%-32s\n", "seq2fasta", "Conversion of database to fasta");
	fprintf(out, "# %16s\t%-32s\n", "dist", "Calculate distance measures between templates");
	fprintf(out, "# %16s\t%-32s\n", "db", "Make statistics on KMA db");
	fprintf(out, "# %16s\t%-32s\n", "merge", "Merge two indexed kma databases");
	fprintf(out, "# %16s\t%-32s\n", "cmp", "Compare two indexed kma databases");
	fprintf(out, "# %16s\t%-32s\n", "update", "Update database to current version");
	fprintf(out, "# %16s\t%-32s\n", "trim", "trim sequences");
	fprintf(out, "# %16s\t%-32s\n", "-c", "Citation");
	fprintf(out, "# %16s\t%-32s\n", "-v", "Version");
	fprintf(out, "# %16s\t%-32s\n", "-h", "Help on alignment and mapping");
	return 1;
}

int main(int argc, char *argv[]) {
	
	int status;
	
	if(argc != 1) {
		if(**++argv == '-') {
			status = kma_main(argc, --argv);
		} else if(--argc && strcmp(*argv, "index") == 0) {
			status = index_main(argc, argv);
		} else if(strcmp(*argv, "shm") == 0) {
			status = shm_main(argc, argv);
		} else if(strcmp(*argv, "seq2fasta") == 0) {
			status = seq2fasta_main(argc, argv);
		} else if(strcmp(*argv, "dist") == 0) {
			status = dist_main(argc, argv);
		} else if(strcmp(*argv, "merge") == 0) {
			status = merge_main(argc, argv);
		} else if(strcmp(*argv, "cmp") == 0) {
			status = cmp_main(argc, argv);
		} else if(strcmp(*argv, "update") == 0) {
			status = update_main(argc, argv);
		} else if(strcmp(*argv, "db") == 0) {
			status = db_main(argc, argv);
		} else if(strcmp(*argv, "trim") == 0) {
			status = trim_main(argc, argv);
		} else {
			fprintf(stderr, "Invalid option:\t%s\n", *argv);
			status = helpmessage(stderr);
		}
	} else {
		status = helpmessage(stdout);
	}
	
	return status;
}
