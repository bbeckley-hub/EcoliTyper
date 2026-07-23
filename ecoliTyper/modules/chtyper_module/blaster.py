#!/home/data1/tools/bin/anaconda/bin/python
from __future__ import division
import sys, os, time, random, re, subprocess
sys.path.append(os.path.abspath("submodules/blaster/submodules/biopython"))
sys.path.append(os.path.abspath("submodules/biopython"))
sys.path.append(os.path.abspath("submodules/"))
from Bio.Blast import NCBIXML
from Bio import SeqIO
#from string import maketrans
import collections

def Blaster(inputfile, databases, db_path, out_path='', min_cov=0.6, threshold=0.9, blast='blastn', cut_off=True):
    min_cov = 100 * float(min_cov)
    threshold = 100 * float(threshold)
    
    gene_align_query = dict()
    gene_align_homo = dict()
    gene_align_sbjct = dict()
    results = dict()
    
    for db in databases:
        db_file = "%s/%s.fsa" % (db_path, db)
        os.system("mkdir -p %s/tmp" % (out_path))
        os.system("chmod 775 %s/tmp" % (out_path))
        out_file = "%s/tmp/out_%s.xml" % (out_path, db)
        
        cmd = "%s -subject %s -query %s -out %s -outfmt '5' -perc_identity %s -dust 'no'" % (
            blast, db_file, inputfile, out_file, threshold)
        sys.stderr.write('LOG: executing - %s\n' % cmd)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        
        result_handle = open(out_file)
        blast_records = NCBIXML.parse(result_handle)
        gene_results = dict()
        best_hsp = dict()
        gene_split = collections.defaultdict(dict)
        gene_align_query[db] = dict()
        gene_align_homo[db] = dict()
        gene_align_sbjct[db] = dict()
        
        for blast_record in blast_records:
            query = blast_record.query
            blast_record.alignments.sort(key=lambda align: -max((len(hsp.query) * (int(hsp.identities) / float(len(hsp.query))) for hsp in align.hsps)))
            for alignment in blast_record.alignments:
                best_e_value = 1
                best_bit = 0
                for hsp in alignment.hsps:
                    if hsp.expect < best_e_value or hsp.bits > best_bit:
                        best_e_value = hsp.expect
                        best_bit = hsp.bits
                        tmp = alignment.title.split(" ")
                        sbjct_header = tmp[1]
                        bit = hsp.bits
                        sbjct_length = alignment.length
                        sbjct_start = hsp.sbjct_start
                        sbjct_end = hsp.sbjct_end
                        gaps = hsp.gaps
                        query_string = str(hsp.query)
                        homo_string = str(hsp.match)
                        sbjct_string = str(hsp.sbjct)
                        contig_name = query.replace(">", "")
                        query_start = hsp.query_start
                        query_end = hsp.query_end
                        HSP_length = len(query_string)
                        perc_ident = int(hsp.identities) / float(HSP_length) * 100
                        strand = 0
                        coverage = ((int(HSP_length) - int(gaps)) / float(sbjct_length))
                        perc_coverage = ((int(HSP_length) - int(gaps)) / float(sbjct_length)) * 100
                        if int(HSP_length) == int(sbjct_length):
                            cal_score = perc_ident * coverage * 100
                        else:
                            cal_score = perc_ident * coverage
                        hit_id = "%s:%s..%s:%s:%f" % (contig_name, query_start, query_end, sbjct_header, cal_score)
                        
                        if sbjct_start > sbjct_end:
                            tmp = sbjct_start
                            sbjct_start = sbjct_end
                            sbjct_end = tmp
                            query_string = reversecomplement(query_string)
                            homo_string = homo_string[::-1]
                            sbjct_string = reversecomplement(sbjct_string)
                            strand = 1
                        
                        if cut_off == True:
                            if perc_coverage > 20:
                                best_hsp = {
                                    'evalue': hsp.expect, 'sbjct_header': sbjct_header, 'bit': bit,
                                    'perc_ident': perc_ident, 'sbjct_length': sbjct_length,
                                    'sbjct_start': sbjct_start, 'sbjct_end': sbjct_end,
                                    'gaps': gaps, 'query_string': query_string,
                                    'homo_string': homo_string, 'sbjct_string': sbjct_string,
                                    'contig_name': contig_name, 'query_start': query_start,
                                    'query_end': query_end, 'HSP_length': HSP_length, 'coverage': coverage,
                                    'cal_score': cal_score, 'hit_id': hit_id, 'strand': strand,
                                    'perc_coverage': perc_coverage
                                }
                        else:
                            best_hsp = {
                                'evalue': hsp.expect, 'sbjct_header': sbjct_header, 'bit': bit,
                                'perc_ident': perc_ident, 'sbjct_length': sbjct_length,
                                'sbjct_start': sbjct_start, 'sbjct_end': sbjct_end,
                                'gaps': gaps, 'query_string': query_string,
                                'homo_string': homo_string, 'sbjct_string': sbjct_string,
                                'contig_name': contig_name, 'query_start': query_start,
                                'query_end': query_end, 'HSP_length': HSP_length, 'coverage': coverage,
                                'cal_score': cal_score, 'hit_id': hit_id, 'strand': strand,
                                'perc_coverage': perc_coverage
                            }
                
                if best_hsp:
                    save = 1
                    if gene_results:
                        tmp_gene_split = gene_split
                        tmp_results = gene_results
                        save, gene_split, gene_results = compare_results(save, best_hsp, tmp_results, tmp_gene_split)
                    if save == 1:
                        gene_results[hit_id] = best_hsp
        
        # --- FIXED: iterate over a copy of keys ---
        for hit_id in list(gene_results.keys()):
            hit = gene_results[hit_id]
            perc_coverage = hit['perc_coverage']
            
            if hit['sbjct_header'] in gene_split and len(gene_split[hit['sbjct_header']]) > 1:
                new_length = calculate_new_length(gene_split, gene_results, hit)
                hit['split_length'] = new_length
                perc_coverage = new_length / float(hit['sbjct_length']) * 100
            
            if perc_coverage >= min_cov:
                if hit['coverage'] == 1:
                    gene_align_query[db][hit_id] = hit['query_string']
                    gene_align_homo[db][hit_id] = hit['homo_string']
                    gene_align_sbjct[db][hit_id] = hit['sbjct_string']
                else:
                    for seq_record in SeqIO.parse(db_file, "fasta"):
                        if seq_record.description == hit['sbjct_header']:
                            gene_align_sbjct[db][hit_id] = str(seq_record.seq)
                            break
                    contig = ''
                    for seq_record in SeqIO.parse(inputfile, "fasta"):
                        if seq_record.description == hit['contig_name']:
                            contig = str(seq_record.seq)
                            break
                    query_seq, homo_seq = get_query_align(hit, contig)
                    gene_align_query[db][hit_id] = query_seq
                    gene_align_homo[db][hit_id] = homo_seq
            else:
                del gene_results[hit_id]
                if hit['sbjct_header'] in gene_split:
                    del gene_split[hit['sbjct_header']]
        
        if gene_results:
            results[db] = gene_results
        else:
            results[db] = "No hit found"
    
    return (results, gene_align_query, gene_align_homo, gene_align_sbjct)

def reversecomplement(seq):
   # Make reverse complement strand
   trans = str.maketrans("ATGC", "TACG")
   return seq.translate(trans)[::-1]

# Function for comparing hits and saving only the best hit
def compare_results(save, best_hsp, tmp_results, tmp_gene_split):
    hit_id = best_hsp['hit_id']
    new_start_query = best_hsp['query_start']
    new_end_query = best_hsp['query_end']
    new_start_sbjct = int(best_hsp['sbjct_start'])
    new_end_sbjct = int(best_hsp['sbjct_end'])
    new_score = best_hsp['cal_score']
    new_db_hit = best_hsp['sbjct_header']
    new_contig = best_hsp['contig_name']
    new_HSP = best_hsp['HSP_length']
    
    # --- FIXED: iterate over a copy of keys ---
    for hit in list(tmp_results.keys()):
        hit_data = tmp_results[hit]
        old_start_query = hit_data['query_start']
        old_end_query = hit_data['query_end']
        old_start_sbjct = int(hit_data['sbjct_start'])
        old_end_sbjct = int(hit_data['sbjct_end'])
        old_score = hit_data['cal_score']
        old_db_hit = hit_data['sbjct_header']
        old_contig = hit_data['contig_name']
        old_HSP = hit_data['HSP_length']
        
        remove_old = 0
        
        if new_db_hit == old_db_hit:
            if new_start_sbjct < old_start_sbjct or new_end_sbjct > old_end_sbjct:
                tmp_gene_split[old_db_hit][hit_id] = 1
                if hit not in tmp_gene_split[old_db_hit]:
                    tmp_gene_split[old_db_hit][hit] = 1
            else:
                if new_score > old_score:
                    remove_old = 1
                    if new_db_hit in tmp_gene_split and hit_id not in tmp_gene_split[new_db_hit]:
                        tmp_gene_split[new_db_hit][hit_id] = 1
                else:
                    save = 0
                    if hit_id != hit:
                        if new_db_hit in tmp_gene_split and hit_id in tmp_gene_split[new_db_hit]:
                            del tmp_gene_split[new_db_hit][hit_id]
                    break
        
        if new_contig == old_contig:
            if old_start_query == new_start_query and old_end_query == new_end_query:
                if best_hsp['perc_ident'] > hit_data['perc_ident']:
                    remove_old = 1
                    if new_db_hit in tmp_gene_split and hit_id not in tmp_gene_split[new_db_hit]:
                        tmp_gene_split[new_db_hit][hit_id] = 1
                elif best_hsp['perc_ident'] == hit_data['perc_ident']:
                    if new_db_hit in tmp_gene_split and hit_id not in tmp_gene_split[new_db_hit]:
                        tmp_gene_split[new_db_hit][hit_id] = 1
                else:
                    save = 0
                    if new_db_hit in tmp_gene_split and hit_id in tmp_gene_split[new_db_hit]:
                        del tmp_gene_split[new_db_hit][hit_id]
                    break
            elif (max(old_end_query, new_end_query) - min(old_start_query, new_start_query)) <= ((old_end_query - old_start_query) + (new_end_query - new_start_query)):
                if new_score > old_score:
                    remove_old = 1
                    if new_db_hit in tmp_gene_split and hit_id not in tmp_gene_split[new_db_hit]:
                        tmp_gene_split[new_db_hit][hit_id] = 1
                elif new_score == old_score:
                    if int(best_hsp['perc_coverage']) == 100 and int(hit_data['perc_coverage']) == 100 and new_HSP > old_HSP:
                        remove_old = 1
                    if new_db_hit in tmp_gene_split and hit_id not in tmp_gene_split[new_db_hit]:
                        tmp_gene_split[new_db_hit][hit_id] = 1
                else:
                    if new_db_hit in tmp_gene_split and hit_id in tmp_gene_split[new_db_hit]:
                        del tmp_gene_split[new_db_hit][hit_id]
                    save = 0
                    break
        
        if remove_old == 1:
            del tmp_results[hit]
            if old_db_hit in tmp_gene_split and hit in tmp_gene_split[old_db_hit]:
                del tmp_gene_split[old_db_hit][hit]
    
    return save, tmp_gene_split, tmp_results

# Function for calcualting new length if the gene is split on several contigs
def calculate_new_length(gene_split, gene_results, hit):
   # Looping over splitted hits and calculate new length
   first = 1
   for split in gene_split[hit['sbjct_header']]:
      #print(split)
      #print("Datebase: %s, cov: %f"%(db, perc_coverage))
      
      new_start = int(gene_results[split]['sbjct_start'])
      new_end = int(gene_results[split]['sbjct_end'])
      
      # Get the frist HSP
      if first == 1:
         new_length = int(gene_results[split]['HSP_length'])
         old_start = new_start
         old_end = new_end
         first = 0
         continue
      if new_start < old_start:
         new_length = new_length + (old_start - new_start)
         old_start = new_start

      if new_end > old_end:
         new_length = new_length + (new_end - old_end)
         old_end = new_end
   
   return(new_length)

# Function for extracting extra seqeunce data to the query alignment if the full reference length
# are not covered
def get_query_align(hit, contig):
   
   # Getting data needed to extract sequences
   query_seq = hit['query_string']
   homo_seq = hit['homo_string']
   sbjct_start = int(hit['sbjct_start'])
   sbjct_end = int(hit['sbjct_end'])
   query_start = int(hit['query_start'])
   query_end = int(hit['query_end'])
   length = int(hit['sbjct_length'])
   
   # If the alignment doesn't start at the first position data is added to the begnning
   if sbjct_start!= 1:
      missing = sbjct_start - 1
      
      if query_start >= missing and hit['strand'] != 1 or hit['strand'] == 1 and missing <= (len(contig) - query_end):
         # Getting the query sequence
         # If the the hit is on the other strand the characters are reversed
         if hit['strand'] == 1:
            start_pos = query_end
            end_pos = query_end + missing
            chars = contig[start_pos:end_pos]
            chars = reversecomplement(chars)
         else:
            start_pos = query_start - missing - 1
            end_pos = query_start - 1
            chars = contig[start_pos:end_pos]
         
         query_seq = chars + str(query_seq)
      else:
         # Getting the query sequence
         # If the the hit is on the other strand the characters are reversed
         if hit['strand'] == 1:
            if query_end == len(contig):
               query_seq = "-" * missing + str(query_seq)
            else:
               start_pos = query_end
               chars = contig[start_pos:]
               chars = reversecomplement(chars)
               
               query_seq = "-" * (missing - len(chars)) + chars + str(query_seq)
         elif query_start < 3:
            query_seq = "-" * missing + str(query_seq)
         else:
            end_pos = query_start - 2
            chars = contig[0:end_pos]
            
            query_seq = "-" * (missing - len(chars)) + chars + str(query_seq)
         
      # Adding to the homo sequence
      spaces = " " * missing
      homo_seq = str(spaces) + str(homo_seq)
   
   # If the alignment dosen't end and the last position data is added to the end
   if sbjct_end < length:
      missing = length - sbjct_end
      
      if missing <= (len(contig) - query_end) and hit['strand'] != 1 or hit['strand'] == 1 and query_start >= missing:
         # Getting the query sequence
         # If the the hit is on the other strand the characters are reversed
         if hit['strand'] == 1:
            start_pos = query_start - missing - 1
            end_pos = query_start - 1
            chars = contig[start_pos:end_pos]
            chars = reversecomplement(chars)
         else:
            start_pos = query_end
            end_pos = query_end + missing
            chars = contig[start_pos:end_pos]
         
         query_seq = query_seq + chars
      else:
         # If the hit is on the other strand the characters are reversed
         if hit['strand'] == 1:
            if query_start < 3:
               query_seq = query_seq + "-" * missing
            else:
               end_pos = query_start - 2
               chars = contig[0:end_pos]
               chars = reversecomplement(chars)
               
               query_seq = query_seq + chars + "-" * (missing - len(chars))
         elif query_end == len(contig):
            query_seq = query_seq + "-" * missing
         else:
            start_pos = query_end
            chars = contig[start_pos:]
            
            query_seq = query_seq + chars + "-" * (missing - len(chars))
            
      # Adding to the homo sequence
      spaces = " " * int(missing)
      homo_seq = str(homo_seq) + str(spaces)
   
   return query_seq, homo_seq
