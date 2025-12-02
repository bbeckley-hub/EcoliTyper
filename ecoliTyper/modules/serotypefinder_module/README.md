SerotypeFinder
===================

This project documents SerotypeFinder service


Documentation
=============

## What is it?

The SerotypeFinder service contains one python script *serotypefinder.py* which is the script of the latest
version of the SerotypeFinder service. SerotypeFinder identifies the serotype in total or partial sequenced
isolates of E. coli.

## Content of the repository
1. serotypefinder.py     - the program
2. README.md
3. Dockerfile   - dockerfile for building the serotypefinder docker container
4. test.fsa     - test fasta file


## Installation

Setting up SerotypeFinder program   
**Warning:** Due to bugs in the BioPython 1.74, if not using the Docker container, do not use that version if not using Python 3.7.
```bash
# Go to wanted location for serotypefinder
cd /path/to/some/dir
# Clone and enter the serotypefinder directory
git clone https://bitbucket.org/genomicepidemiology/serotypefinder.git
cd serotypefinder
```

Build Docker container
```bash
# Build container
docker build -t serotypefinder .
```

#Download and install SerotypeFinder database
```bash
# Go to the directory where you want to store the serotypefinder database
cd /path/to/some/dir
# Clone database from git repository (develop branch)
git clone https://bitbucket.org/genomicepidemiology/serotypefinder_db.git
cd serotypefinder_db
STFinder_DB=$(pwd)
# Install SerotypeFinder database with executable kma_index program
python3 INSTALL.py kma_index
```

If kma_index has no bin install please install kma_index from the kma repository:
https://bitbucket.org/genomicepidemiology/kma

## Dependencies
In order to run the program without using docker, Python 3.5 (or newer) should be installed along with the following versions of the modules (or newer).

#### Modules
- cgecore 1.5.5
- tabulate 0.7.7

Modules can be installed using the following command. Here, the installation of the module cgecore is used as an example:
```bash
pip3 install cgecore
```
#### KMA and BLAST
Additionally KMA and BLAST version 2.8.1 or newer should be installed.
The newest version of KMA and BLAST can be installed from here:
```url
https://bitbucket.org/genomicepidemiology/kma
```

```url
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

## Usage

The program can be invoked with the -h option to get help and more information of the service.
Run Docker container


```bash
# Run serotypefinder container
docker run --rm -it \
       -v $STFinder_DB:/database \
       -v $(pwd):/workdir \
       serotypefinder -i [INPUTFILE] -o . [-x] [-mp] [-t] [-l] [-p] [-d] [-tmp]
```

When running the docker file you have to mount 2 directory: 
 1. serotypefinder_db (SerotypeFinder database) downloaded from bitbucket
 2. An output/input folder from where the input file can be reached and an output files can be saved. 
Here we mount the current working directory (using $pwd) and use this as the output directory, 
the input file should be reachable from this directory as well.
 
-i INPUTFILE        input file (fasta or fastq) relative to pwd 

-o OUTDIR        output directory relative to pwd

-x	extended output    Will create an extented output

-mp METHOD_PATH        Path to executable of the method to be used (kma or blast)

-d DATABASE        Choose specific database

-p DB_PATH        Path to database directory

-tmp TMP_DIR        Temporary directory for storage of results from external software

-l MINCOV        Set threshold for minimum coverage

-t THRESHOLD        Set threshold for minimum identity


## Web-server

A webserver implementing the methods is available at the [CGE website](http://www.genomicepidemiology.org/) and can be found here:
https://cge.cbs.dtu.dk/services/SerotypeFinder/


## The Latest Version


The latest version can be found at
https://bitbucket.org/genomicepidemiology/serotypefinder/overview

## Documentation


The documentation available as of the date of this release can be found at
https://bitbucket.org/genomicepidemiology/serotypefinder/overview.


Citation
=======

When using the method please cite:

Flemming Scheutz, SSI
[Epub ahead of print]

References
=======

1. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture and applications. BMC Bioinformatics 2009; 10:421. 
2. Clausen PTLC, Aarestrup FM, Lund O. Rapid and precise alignment of raw reads against redundant databases with KMA. BMC Bioinformatics 2018; 19:307. 

License
=======

Copyright (c) 2014, Ole Lund, Technical University of Denmark
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
