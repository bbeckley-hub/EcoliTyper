#!/home/data1/tools/bin/anaconda/bin/python
from __future__ import division
import sys, os, time, random, re, subprocess
from argparse import ArgumentParser
from tabulate import tabulate
from blaster import Blaster

##########################################################################
# FUNCTIONS
##########################################################################

def text_table(title, headers, rows, table_format='psql'):
   ''' Create text table

   USAGE:
      >>> from tabulate import tabulate
      >>> title = 'My Title'
      >>> headers = ['A','B']
      >>> rows = [[1,2],[3,4]]
      >>> print(text_table(title, headers, rows))
      +-----------+
      | My Title  |
      +-----+-----+
      |   A |   B |
      +=====+=====+
      |   1 |   2 |
      |   3 |   4 |
      +-----+-----+
   '''
   # Create table
   table = tabulate(rows, headers, tablefmt=table_format)
   # Prepare title injection
   width = len(table.split('\n')[0])
   tlen = len(title)
   if tlen + 4 > width:
      # Truncate oversized titles
      tlen = width - 4
      title = title[:tlen]
   spaces = width - 2 - tlen
   left_spacer = ' '*int(spaces / 2)
   right_spacer = ' '*(spaces - len(left_spacer))
   # Update table with title
   table = '\n'.join(['+%s+'%('-'*(width-2)),
                      '|%s%s%s|'%(left_spacer, title, right_spacer),
                      table, '\n'])
   return table

##########################################################################
#	DEFINE GLOBAL VARIABLES
##########################################################################

global database_path, databases, min_length, threshold

##########################################################################
# PARSE COMMAND LINE OPTIONS
##########################################################################

parser = ArgumentParser()
parser.add_argument("-i", "--inputfile", dest="inputfile",help="Input file", default='')
parser.add_argument("-o", "--outputPath", dest="out_path",help="Path to blast output", default='')
parser.add_argument("-b", "--blastPath", dest="blast_path",help="Path to blast", default='blastn')
parser.add_argument("-p", "--databasePath", dest="db_path",help="Path to the databases", default='')
parser.add_argument("-d", "--databases", dest="databases",help="Databases chosen to search in - if non is specified all is used", default=None)
parser.add_argument("-l", "--min_cov", dest="min_cov",help="Minimum coverage", default=0.60)
parser.add_argument("-t", "--threshold", dest="threshold",help="Blast threshold for identity", default=0.90)
args = parser.parse_args()


##########################################################################
# MAIN
##########################################################################

# Defining varibales

min_cov = args.min_cov
threshold = args.threshold

# Check if valid database is provided
if args.db_path is None:
      sys.exit("Input Error: No database directory was provided!\n")
elif not os.path.exists(args.db_path):
   sys.exit("Input Error: The specified database directory does not"
                       " exist!\n")
else:
   # Check existence of config file
   db_config_file = '%s/config'%(args.db_path)
   if not os.path.exists(db_config_file):
      sys.exit("Input Error: The database config file could not be "
                          "found!")
   # Save path
   db_path = args.db_path

# Check if valid input file is provided
if args.inputfile is None:
   sys.exit("Input Error: No Input were provided!\n")
elif not os.path.exists(args.inputfile):
   sys.exit("Input Error: Input file does not exist!\n")
else:
    inputfile = args.inputfile

# Check if valid output directory is provided
if not os.path.exists(args.out_path):
   sys.exit("Input Error: Output dirctory does not exists!\n")
else:
   out_path = args.out_path

# Check if valid path to BLAST is provided
if not os.path.exists(args.blast_path):
   sys.exit("Input Error: The path to BLAST does not exists!\n")
else:
   blast = args.blast_path

# Check if databases and config file are correct/correponds
if args.databases is '':
      sys.exit("Input Error: No database was specified!\n")
else:
   dbs = dict()
   extensions = []
   databases = []
   with open(db_config_file) as f:
      for l in f:
         l = l.strip()
         if l == '': continue
         if l[0] == '#':
            if 'extensions:' in l:
               extensions = [s.strip() for s in l.split('extensions:')[-1].split(',')]
            continue
         tmp = l.split('\t')
         if len(tmp) != 3:
            sys.exit(("Input Error: Invalid line in the database "
                      "config file!\nA proper entry requires 3 tab "
                      "separated columns!\n%s")%(l))
         db_prefix = tmp[0].strip()
         name = tmp[1].split('#')[0].strip()
         description = tmp[2]
         # Check if all db files are present
         for ext in extensions:
            db = "%s/%s.%s"%(db_path, db_prefix, ext)
            if not os.path.exists(db):
               sys.exit(("Input Error: The database file (%s) "
                         "could not be found!")%(db_path))
         if db_prefix not in dbs:
            dbs[db_prefix] = []
            databases.append(db_prefix)
         dbs[db_prefix].append(name)
   if len(dbs) == 0:
      sys.exit("Input Error: No databases were found in the "
                          "database config file!")

   if args.databases is None:
      # Use all available databases from the config file
      pass
   else:
      # Handle multiple databases
      args.databases = args.databases.split(',')
      # Check that the ResFinder DBs are valid
      databases = []
      for db_prefix in args.databases:
         if db_prefix in dbs:
            databases.append(db_prefix)
         else:
            sys.exit("Input Error: Provided database was not "
                              "recognised! (%s)\n"%db_prefix)

# Calling blast and parsing output
results, query_align, homo_align, sbjct_align = Blaster(inputfile, databases,
                                                        db_path, out_path,
                                                        min_cov, threshold, blast)
# Making output files
tab_file = open(out_path+"/results_tab.txt", 'w')
table_file = open(out_path+"/results_table.txt", 'w')
txt_file = open(out_path+"/results.txt", 'w')
ref_file = open(out_path+"/Reference_gene_seq.fsa", 'w')
hit_file = open(out_path+"/Hit_in_genome_seq.fsa", 'w')

# Write the header for the tab file
header_line =  ["Allele type", "Identity", "Alignment Length/Gene Length", "Position in reference", "Contig", "Position in contig", "Accession no."]
tab_file.write("%s\n"%("\t".join(header_line)))

# Getting and writing out the results
titles = dict()
rows = dict()
headers = dict()
txt_file_seq_text = dict()
split_print = list()

for db in databases:
   profile = ''.join(dbs[db])
   if results[db] == "No hit found":
      #titles[db] = profile
      #headers[db] = results[db]
      #rows[db] = ""
      table_file.write("%s\n%s\n\n"%(profile, results[db]))
   else:
      titles[db] = "%s"%(profile)
      headers[db] = header_line
      table_file.write("%s\n"%(profile))
      table_file.write("%s\n"%("\t".join(header_line)))
      rows[db] = list()
      txt_file_seq_text[db] = list()
      for hit in results[db]:
         header = results[db][hit]["sbjct_header"]
         if header in split_print:
            next
         else:
            split_print.append(header)
         tmp = header.split("_")
         if len(tmp) == 2:
            gene = tmp[0]
            acc = tmp[1]
         else:
            gene = header
            acc = "Not available"
         ID = results[db][hit]["perc_ident"]
         sbjt_length = results[db][hit]["sbjct_length"]
         #
         # If else clause contain duplicate code. This could be avoided using
         # the get method of dict:
         # length = results[db][hit].get("split_length", "HSP_length")
         #
         if "split_length" in results[db][hit]:
            HSP = results[db][hit]["split_length"]
            positions_contig = "%s..%s"%(results[db][hit]["query_start"], results[db][hit]["query_end"])
            positions_ref = "%s..%s"%(results[db][hit]["sbjct_start"], results[db][hit]["sbjct_end"])
            contig_name = results[db][hit]["contig_name"]
         else:
            HSP = results[db][hit]["HSP_length"]
            positions_contig = "%s..%s"%(results[db][hit]["query_start"], results[db][hit]["query_end"])
            positions_ref = "%s..%s"%(results[db][hit]["sbjct_start"], results[db][hit]["sbjct_end"])
            contig_name = results[db][hit]["contig_name"]

         # Write tabels
         table_file.write("%s\t%.2f\t%s/%s\t%s\t%s\t%s\t%s\n"%(gene, ID, HSP, sbjt_length, positions_ref, contig_name, positions_contig, acc))
         tab_file.write("%s\t%.2f\t%s/%s\t%s\t%s\t%s\t%s\n"%(gene, ID, HSP, sbjt_length, positions_ref, contig_name, positions_contig, acc))

         # Saving the output to write the txt result table
         hsp_length = "%s/%s"%(HSP, sbjt_length)
         rows[db].append([gene, ID, hsp_length, positions_ref, contig_name, positions_contig, acc])

         # Writing subjet/ref sequence
         ref_seq = sbjct_align[db][hit]
         ref_file.write(">%s_%s\n"%(gene, acc))
         for i in range(0, len(ref_seq), 60):
            ref_file.write("%s\n"%(ref_seq[i:i + 60]))

         # Getting the header and text for the txt file output
         #>aac(2')-Ic: PERFECT MATCH, ID: 100.00%, HSP/Length: 546/546, Positions in reference: 1..546, Contig name: gi|375294201|ref|NC_016768.1|, Position: 314249..314794
         sbjct_start = results[db][hit]["sbjct_start"]
         sbjct_end = results[db][hit]["sbjct_end"]
         text = "%s, ID: %.2f %%, Alignment Length/Gene Length: %s/%s, Positions in reference: %s..%s, Contig name: %s, Position: %s"%(gene, ID, HSP, sbjt_length, sbjct_start, sbjct_end, contig_name, positions_contig)
         hit_file.write(">%s\n"%text)

         # Writing query/hit sequence
         hit_seq = query_align[db][hit]
         for i in range(0, len(hit_seq), 60):
            hit_file.write("%s\n"%(hit_seq[i:i + 60]))

         # Saving the output to print the txt result file allignemts
         txt_file_seq_text[db].append((text, ref_seq, homo_align[db][hit], hit_seq))

      table_file.write("\n")

# Writing the txt file

# Writing table txt for all hits
for db in titles:
   # Txt file table
   table = text_table(titles[db], headers[db], rows[db])
   txt_file.write(table)

# Writing alignment txt for all hits
for db in titles:
   # Txt file alignments
   txt_file.write("##################### %s #####################\n"%(db))
   for text in txt_file_seq_text[db]:
      txt_file.write("%s\n\n"%(text[0]))
      for i in range(0, len(text[1]), 60):
         txt_file.write("%s\n"%(text[1][i:i + 60]))
         txt_file.write("%s\n"%(text[2][i:i + 60]))
         txt_file.write("%s\n\n"%(text[3][i:i + 60]))
      txt_file.write("\n")
