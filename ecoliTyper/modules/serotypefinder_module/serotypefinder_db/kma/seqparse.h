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

#include "filebuff.h"
#include "qseqs.h"

/* determine format */
int openAndDetermine(FileBuff *inputfile, char *filename);
/* get entry from fastafile */
int FileBuffgetFsa(FileBuff *src, Qseqs *header, Qseqs *qseq, char *trans);
int FileBuffgetFsaSeq(FileBuff *src, Qseqs *qseq, char *trans);
/* get entry from fastq file */
int FileBuffgetFq(FileBuff *src, Qseqs *header, Qseqs *qseq, Qseqs *qual, char *trans);
int FileBuffgetFqSeq(FileBuff *src, Qseqs *qseq, Qseqs *qual, char *trans);
/* get phred scale */
int getPhredFileBuff(FileBuff *dest);
