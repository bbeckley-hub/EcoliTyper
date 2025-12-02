### Example Command

```bash
python3 CHTyper-1.0.py -p [DB_PATH] -i [INPUT_FASTA] -o [OUTPUT_PATH] -b [path/to/blast/binary/]
```

### Database

```bash
git clone https://bitbucket.org/genomicepidemiology/chtyper_db/src/master/
```

### BLAST

You need to have BLAST installed - which can be obtained from:

```url
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

### Running CHTyper


```bash
Command specifications: 

"-i", Input file
"-o", Path to blast output
"-b", Path to blast
"-p", Path to the databases
"-d", Databases chosen to search in - if non is specified all is used
"-l", Minimum coverage, default=0.60
"-t", Blast threshold for identity, default=0.90
```

## Citation
For publication of results, please cite:

CHTyper, a Web Tool for Subtyping of Extraintestinal Pathogenic Escherichia coli Based on the fumC and fimH Alleles.
Roer L, Johannesen TB, Hansen F, Stegger M, Tchesnokova V, Sokurenko E, et al.
J. Clin. Microbiol. 2018;56:e00063-18.
https://doi.org/10.1128/JCM.00063-18
