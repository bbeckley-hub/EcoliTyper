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
#include <stdio.h>
#include <string.h>
#include "pherror.h"
#include "threader.h"
#include "tmp.h"

FILE * tmpF(const char *location) {
	
	static int tmpNum = 0;
	static char *tmpname = 0, *dirname = 0, *filename = 0;
	static volatile int Lock = 0;
	volatile int *lock = &Lock;
	int fd;
	FILE *file;
	
	lock(lock);
	if(location) {
		/* set location for tmpfiles */
		if(filename) {
			free(filename);
			filename = 0;
		} else if(dirname) {
			free(dirname);
			dirname = 0;
		}
		if(*location) {
			if(location[strlen(location) - 1] != '/') {
				filename = smalloc(strlen(location) + 15);
				tmpname = filename + sprintf(filename, "%s.tmp", location);
				dirname = 0;
			} else {
				dirname = smalloc(strlen(location) + 12);
				tmpname = dirname + sprintf(dirname, "%s.kma-", location);
				strcpy(tmpname, "XXXXXX");
				filename = 0;
			}
		}
		file = 0;
	} else if(dirname) {
		if((fd = mkstemp(dirname)) == -1) {
			ERROR();
		}
		if((file = fdopen(fd, "wb+"))) {
			unlink(dirname);
		}
		strcpy(tmpname, "XXXXXX");
	} else if(filename) {
		/* open tmpfile on previous location */
		fd = sprintf(tmpname, "%d", tmpNum++);
		if((file = fopen(filename, "wb+"))) {
			unlink(filename);
		}
		*tmpname = 0;
	} else if(!tmpNum) {
		/* set default according to TMPDIR */
		if((tmpname = getenv("TMPDIR"))) {
			dirname = smalloc(strlen(tmpname) + 13);
			tmpname = dirname + sprintf(dirname, "%s/.kma-", tmpname);
			strcpy(tmpname, "XXXXXX");
		}
		tmpNum = 1;
		unlock(lock);
		
		return tmpF(0);
	} else {
		/* open a "normal" tmp file */
		file = tmpfile();
	}
	unlock(lock);
	
	return file;
}
