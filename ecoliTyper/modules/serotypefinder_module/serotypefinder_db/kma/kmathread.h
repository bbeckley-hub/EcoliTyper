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
#include <pthread.h>
#include <unistd.h>

int usleep(unsigned usec);

#define lock(exclude) while(__sync_lock_test_and_set(exclude, 1)) {while(*exclude) {usleep(100);}}
#define lockTime(exclude, spin) while(__sync_lock_test_and_set(exclude, 1)) {while(*exclude) {usleep(spin);}}
#define unlock(exclude) (__sync_lock_release(exclude))
#define wait_atomic(src) while(src) {usleep(100);}
