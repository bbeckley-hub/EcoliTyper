/* Philip T.L.C. Clausen Sep 2021 plan@dtu.dk */

/*
 * Copyright (c) 2021, Philip Clausen, Technical University of Denmark
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

#include <stdio.h>
#include <math.h>
#include "assembly.h"
#include "tsv.h"
#define printsvField(flag, outfile, format, feature) if(flag & 1) {fprintf(outfile, format, feature, (flag >>= 1) ? '\t' : '\n');} else {flag >>= 1;}

void initsv(FILE *outfile, long unsigned flag) {
	
	/* mask off invalid flags */
	flag &= 131071;
	
	printsvField(flag, outfile, "%s%c", "Template_Name");
	printsvField(flag, outfile, "%s%c", "Template_Length");
	printsvField(flag, outfile, "%s%c", "Template_Identity");
	printsvField(flag, outfile, "%s%c", "Template_Coverage");
	printsvField(flag, outfile, "%s%c", "Template_Depth");
	printsvField(flag, outfile, "%s%c", "Query_Identity");
	printsvField(flag, outfile, "%s%c", "Query_Coverage");
	printsvField(flag, outfile, "%s%c", "Query_Depth");
	printsvField(flag, outfile, "%s%c", "Read_Count_Map");
	printsvField(flag, outfile, "%s%c", "Read_Count_Aln");
	printsvField(flag, outfile, "%s%c", "Score");
	printsvField(flag, outfile, "%s%c", "Expected");
	printsvField(flag, outfile, "%s%c", "q_value");
	printsvField(flag, outfile, "%s%c", "p_value");
	printsvField(flag, outfile, "%s%c", "ConClave_Score");
	printsvField(flag, outfile, "%s%c", "ConClave_Quality");
	printsvField(flag, outfile, "%s%c", "Reassign_template");
}

void printsv(FILE *outfile, long unsigned flag, char *template_name, Assem *aligned_assem, int t_len, unsigned readCount, long unsigned read_score, double expected, double q_value, double p_value, long unsigned ConClave_Score, unsigned char *reassign) {
	
	//printsv(tsv_out, tsv, template_name, aligned_assem, t_len, readCount, read_score, expected, q_value, p_value, alignment_scores[template]);
	
	/* mask off invalid flags */
	flag &= 131071;
	
	/* 0 no file */
	if(!flag) {
		return;
	}
	
	/* 1 Template_Name */
	printsvField(flag, outfile, "%s%c", template_name);
	
	/* 2	Template_Length */
	printsvField(flag, outfile, "%d%c", t_len);
	
	/* 4	Template_Identity */
	printsvField(flag, outfile, "%f%c", 100.0 * aligned_assem->cover / t_len);
	
	/* 8	Template_Coverage */
	printsvField(flag, outfile, "%f%c", 100.0 * aligned_assem->aln_len / t_len);
	
	/* 16	Template_Depth */
	printsvField(flag, outfile, "%f%c", (double)(aligned_assem->depth) / t_len);
	
	/* 32	Query_Identity */
	printsvField(flag, outfile, "%f%c", 100.0 * aligned_assem->cover / aligned_assem->aln_len);
	
	/* 64	Query_Coverage */
	printsvField(flag, outfile, "%f%c", 100.0 * t_len / aligned_assem->aln_len);
	
	/* 128	Query_Depth */
	printsvField(flag, outfile, "%f%c", (double)(aligned_assem->depth) / aligned_assem->aln_len);
	
	/* 256	Read_Count_Map */
	printsvField(flag, outfile, "%d%c", readCount);
	
	/* 512	Read_Count_Aln */
	printsvField(flag, outfile, "%d%c", aligned_assem->readCountAln);
	
	/* 1024	Score */
	printsvField(flag, outfile, "%lu%c", read_score);
	
	/* 2048	Expected */
	printsvField(flag, outfile, "%f%c", expected);
	
	/* 4096	q_value */
	printsvField(flag, outfile, "%f%c", q_value);
	
	/* 8192	p_value */
	printsvField(flag, outfile, "%e%c", p_value);
	
	/* 16384	ConClave_Score */
	printsvField(flag, outfile, "%lu%c", ConClave_Score);
	
	/* 32768	ConClave_Quality */
	printsvField(flag, outfile, "%f%c", 40.0 * read_score / ConClave_Score * log(read_score));
	
	/* 65536	ConClave_Quality */
	printsvField(flag, outfile, "%s%c", reassign);
}
