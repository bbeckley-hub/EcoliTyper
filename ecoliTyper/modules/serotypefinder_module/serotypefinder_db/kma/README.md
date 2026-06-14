[![install with bioconda](https://img.shields.io/badge/install%20with-genomicepidemiology-brightgreen.svg?style=flat)](https://anaconda.org/genomicepidemiology/kma)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/kma/README.html)
# Getting Started #

```
git clone https://bitbucket.org/genomicepidemiology/kma.git
cd kma && make

./kma index -i templates.fsa.gz -o templates
./kma -i reads_se.fq.gz -o output/name -t_db templates
./kma -ipe reads_1.fq.gz reads_2.fq.gz -o output/name -t_db templates
```

# Introduction #
KMA is a mapping method designed to map raw reads directly against redundant databases, in an 
ultra-fast manner using seed and extend. KMA is particulary good at aligning 
high quality reads against highly redundant databases, where unique matches often does 
not exist. It works for long low quality reads as well, such as those from Nanopore. 
Non-unique matches are resolved using the "ConClave" sorting scheme, and a consensus sequence are outputtet
in addition to other common attributes, based on our users demands.

If you use KMA for your published research, then please cite:
Philip T.L.C. Clausen, Frank M. Aarestrup & Ole Lund, 
"Rapid and precise alignment of raw reads against redundant databases with KMA", 
BMC Bioinformatics, 2018;19:307.


# Usage #
A more detailed guide on KMA and its options can be found in the pdf "KMAspecification.pdf".
For practical reasons you might want to add kma to your path, this is usually done with:

```
mv kma ~/bin/
```

## Indexing ##
In order to use KMA for mapping, the databases need to indexed. 
This is done with kma index, the most important options are described below:

```
-i Input fasta file(s), space separated. By default kma index reads from stdin.
-o Output name, the name given to the database.
-k kmersize used for indexing the database.
-k_t kmersize used to identify template candidates when running KMA.
-k_i kmersize used when performing alignments between two sequences.
```

Example of use:

```
kma index -i templates.fsa.gz -o database/name
```

## Mapping ##
The default settings of KMA is set so that it should work on most cases, 
but it is evident to use the right options for the right problems.
Some of the most important options:

```
-i Inputfile(s), default is read from stdin. All input options takes as many files as you wish in fastq or fasta format, space separated.
-ipe Inputfile(s), paired end. The file pairs should be placed right after each other.
-int Inputfile(s), interleaved.
-o Output destination.
-t_db Database from indexing.
-mem_mode *.index and *.seq are not loaded into memory, which enables one to map against larger databases. Templates are chosen using k-mer counting.
-dense Skip insertions when making the consensus sequence.
-ref_fsa - will be substituted with n in the consensus sequence.
-matrix Gives the counts all all called bases at each position in each mapped template. Columns are: reference base, A count, C count, G count, T count, N count, - count.
-mp Minimum phred-score.
-Mt1 Match to only one template in the database.
-ID Minimum identity to output template match.
-apm Paired end method, p: Reward if pairing the reads, u: unite best template matches in each read if possible, f force paired reads to pair.
-1t1 One read to one template, no splicing performed. Well suited for short reads and whole genome mapping.
-bc90 Basecalls should be significantly overrepresented, and have at least 90% agreement.
-bcNano Basecalls optimized for nanopore sequencing.
-mrs minimum alignment score normalized to alignment length.
```

Examples of running KMA:

Short read mapping, one read maps only to one template:
```
kma -i singleEndReads.fq.gz -ipe pairedEnd_*.fq.gz -o output/name -t_db database/name -1t1
```

Long read mapping against a database of genes:
```
kma -i someLongReads.fq.gz -o output/name -t_db database/name -bcNano -bc 0.7
```

Whole genome mapping with nanopore reads:
```
kma -i nanoporeReads.fq.gz -o output/name -t_db database/name -mem_mode -mp 20 -mrs 0.0 -bcNano -bc 0.7
```

Whole genome mapping against a single genome in the database (still nanopore), specified by template number (here 2).
This is the same as the line number of the wanted sequences in the \*.name file from the indexing.
```
kma -i nanoporeReads.fq.gz -o output/name -t_db database/name -mp 20 -bcNano -bc 0.7 -Mt1 2
```


# Result Explanation #
When the mapping is done KMA will produce the following files:

1. \*.res A result overview giving the most common statistics for each mapped template.
2. \*.fsa The consensus sequences drawn from the alignments.
3. \*.aln The consensus alignment of the reads against their template.
4. \*.frag.gz Mapping information on each mapped read, columns are: read, number of equally well mapping templates, mapping score, start position, end position (w.r.t. template), the choosen template.
5. \*.mat.gz Base counts on each position in each template, (only if -matrix is enabled)

# Shared memory #
The databases of KMA can be put into shared memory, this enables you to align several 
samples at ones while only having the database loaded in one place. 
KMA_SHM can be used to put the databases into shared memory, and most importantly take them 
down afterwards. It is of the highest concern that these shared databases are used correctly, 
as they will take up a part of memory and they will exist there until they are taken down again, 
the computer is restarted or computer breaks down.

Example of setting up a database in shared memory.
```
kma shm -t_db database/name -shmLvl 1
```

-shmLvl specifies how much of the database there should be stored in shared memory, use.
-shm-h to get an overview of the different levels.

Example of taking it down again, always remember to do this then it is no longer needed:
```
kma shm -t_db database/name -shmLvl 1 -destroy
```

# Installation Requirements #
In order to install KMA, you need to have a C-compiler and zlib development files installed.
Zlib development files can be installed on unix systems with:
```
sudo apt-get install libz-dev
```

# Help #
Usage and options are available with the "-h" option on all three programs.
If in doubt, please mail any concerns or problems to: *plan@dtu.dk*.

# Citation #
1. Philip T.L.C. Clausen, Frank M. Aarestrup & Ole Lund, "Rapid and precise alignment of raw reads against redundant databases with KMA", BMC Bioinformatics, 2018;19:307.

# License #
Copyright (c) 2017, Philip Clausen, Technical University of Denmark
All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
