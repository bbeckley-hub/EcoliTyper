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
#include <pthread.h>
#include <unistd.h>
#include "kma.h"
#include "kmapipe.h"
#include "pherror.h"
#include "threader.h"
#ifdef _WIN32
#include <windows.h>
#else
#include <paths.h>
#include <sys/types.h>
#include <sys/wait.h>
#endif

FILE * (*kmaPipe)(const char*, const char*, FILE*, int*) = &kmaPipeThread;

void * pipeThreader(void *arg) {
	
	char *cmd[2];
	Pid *pidPtr = arg;
	
	cmd[0] = pidPtr->cmd;
	cmd[1] = (char *) pidPtr->ioStream;
	
	/* start child work */
	kma_main(0, cmd);
	
	/* close stream */
	fclose(pidPtr->ioStream);
	
	return NULL;
}

FILE * kmaPipeThread(const char *cmd, const char *type, FILE *ioStream, int *status) {
	
	/* kmaPipe is a combination of popen and pclose, but allows for binary mode */
	static volatile int Lock = 0;
	static Pid *pidlist = 0;
	volatile int *lock = &Lock;
	int pdes[2];
	pthread_t id;
	Pid *src, *last, *dest;
	
	if(cmd && type) {
		/* check mode */
		if((*type != 'r' && *type != 'w') || (type[1] != 0 && type[2] != 0)) {
			errno = EINVAL;
			ERROR();
		}
		
		/* create pipe */
		dest = malloc(sizeof(Pid));
		if(!dest) {
			ERROR();
		} else if(pipe(pdes) != 0) {
			ERROR();
		}
		dest->cmd = (char *) cmd;
		
		/* convert pipe */
		if (*type == 'r') {
			dest->fp = fdopen(pdes[0], type);
			dest->ioStream = fdopen(pdes[1], "wb");
		} else {
			dest->fp = fdopen(pdes[1], type);
			dest->ioStream = fdopen(pdes[0], "rb");
		}
		if(!dest->fp || !dest->ioStream) {
			ERROR();
		}
		
		/* spawn process */
		if((errno = pthread_create(&dest->id, NULL, &pipeThreader, dest))) {
			ERROR();
		}
		
		/* Link into list of file descriptors. */
		lock(lock);
		dest->next = pidlist;
		pidlist = dest;
		unlock(lock);
		
		return dest->fp;
	} else {
		*status = 0;
		/* Get file pointer. */
		lock(lock);
		for (last = 0, src = pidlist; src->fp != ioStream; last = src, src = src->next) {
			if(!src) {
				*status = 1;
				return 0;
			}
		}
		unlock(lock);
		
		/* close stream and get exit status */
		*status = 0;
		id = src->id;
		//if((errno = pthread_join(id, NULL))) {
		if((errno = pthread_join(id, (void *) &status))) {
			ERROR();
		}
		fclose(ioStream);
		
		/* Remove the entry from the linked list. */
		lock(lock);
		for (last = 0, src = pidlist; src->id != id; last = src, src = src->next) {
			if(!src) {
				*status = 1;
				return 0;
			}
		}
		if (!last) {
			pidlist = src->next;
		} else {
			last->next = src->next;
		}
		unlock(lock);
		free(src);
		
		return 0;
	}
}

FILE * kmaPipeFork(const char *cmd, const char *type, FILE *ioStream, int *status) {
	
	/* kmaPipe is a combination of popen and pclose, but allows for binary mode */
	static volatile int Lock = 0;
	static Pid *pidlist = 0;
	volatile int *lock = &Lock;
	int exit_status, pdes[2];
	char *argv[2];
	pid_t id;
	volatile pid_t pid;
	Pid *src, *last, *dest;
	
	if(cmd && type) {
		/* check mode */
		if((*type != 'r' && *type != 'w') || (type[1] != 0 && type[2] != 0)) {
			errno = EINVAL;
			ERROR();
		}
		
		/* create pipe */
		dest = malloc(sizeof(Pid));
		if(!dest) {
			ERROR();
		} else if(pipe(pdes) != 0) {
			ERROR();
		}
		dest->cmd = (char *) cmd;
		
		/* spawn process */
		pid = fork();
		if(pid < 0) {
			ERROR();
		} else if(pid == 0) {
			errno = 0;
			lock(lock);
			for(src = pidlist; src; src = src->next) {
				close(fileno(src->fp));
			}
			unlock(lock);
			
			if(*type == 'r') {
				close(pdes[0]);
				dest->ioStream = fdopen(pdes[1], "wb");
			} else {
				close(pdes[1]);
				dest->ioStream = fdopen(pdes[0], "rb");
			}
			if(!dest->ioStream) {
				ERROR();
			}
			
			/* start child work */
			argv[0] = dest->cmd;
			argv[1] = (char *) dest->ioStream;
			exit_status = kma_main(0, argv);
			
			/* close stream */
			fclose(dest->ioStream);
			
			/* kill child */
			_exit(exit_status);
		} else if(*type == 'r') {
			close(pdes[1]);
			dest->fp = fdopen(pdes[0], type);
		} else {
			close(pdes[0]);
			dest->fp = fdopen(pdes[1], type);
		}
		if(!dest->fp) {
			ERROR();
		}
		
		/* Link into list of file descriptors. */
		lock(lock);
		dest->pid = pid;
		dest->next = pidlist;
		pidlist = dest;
		unlock(lock);
		
		return dest->fp;
	} else {
		*status = 0;
		/* Get file pointer. */
		lock(lock);
		for (last = 0, src = pidlist; src->fp != ioStream; last = src, src = src->next) {
			if(!src) {
				*status = 1;
				unlock(lock);
				return 0;
			}
		}
		unlock(lock);
		
		/* close stream and get exit status */
		*status = 0;
		#ifndef _WIN32
		while ((pid = waitpid(src->pid, status, 0)) == -1 && errno == EINTR) {
			usleep(100);
		}
		if(WIFEXITED(*status)) {
			exit_status = WEXITSTATUS(*status);
			*status = exit_status;
		} else {
			*status = 1;
		}
		#else
		WaitForSingleObject(src->pid, INFINITE);
		#endif
		fclose(ioStream);
		
		/* Remove the entry from the linked list. */
		lock(lock);
		id = src->pid;
		for (last = 0, src = pidlist; src->pid != id; last = src, src = src->next) {
			if(!src) {
				*status = 1;
				return 0;
			}
		}
		if (!last) {
			pidlist = src->next;
		} else {
			last->next = src->next;
		}
		unlock(lock);
		free(src);
		
		return 0;
	}
}
