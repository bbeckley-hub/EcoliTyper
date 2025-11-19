# EcoliTyper 🧬

**Unified MLST + Serotyping + Clermont Phylotyping for *Escherichia coli***

[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17226894.svg)](https://doi.org/10.5281/zenodo.17226894)

## 📖 Overview

**EcoliTyper** is a comprehensive, high-performance genotyping tool for *Escherichia coli* that unifies three essential typing methods into a single, streamlined workflow:

- 🔬 **MLST (Multi-Locus Sequence Typing)** - For precise strain classification
- 🦠 **Serotyping** - For O and H antigen identification using CGE's highly curated database
- 🧬 **Clermont Phylotyping** - For rapid phylogenetic group determination

### 🚀 Key Features

- **⚡ Blazing Fast**: Parallel processing with configurable threads for high-throughput analysis
- **🎯 Unified Workflow**: Single command for all three typing methods
- **📁 Wildcard Support**: Process multiple genomes with glob patterns (`*.fna`, `*.fasta`)
- **🔧 Self-Contained**: Bundled with CGE's highly curated SerotypeFinder database
- **📊 Comprehensive Output**: Multiple formats (TSV, JSON) for easy downstream analysis
- **🎨 User-Friendly**: Beautiful ASCII art interface with science quotes

## 🛠️ Installation (RECOMMENDED)

```bash
conda create -n ecolityper -c bbeckley-hub -c conda-forge ecolityper python=3.9 -y
conda activate ecolityper
```

### Prerequisites

- Python 3.6 or higher
- BLAST+ tools (`blastn`, `makeblastdb`)
- Perl (for MLST)

### Quick Install

```bash
# Clone the repository
git clone https://github.com/bbeckley-hub/ecoliTyper.git
cd ecoliTyper

# Install the package
pip install -e .

# Install EzClermont (required for phylotyping)
pip install ezclermont

# Environment Check
ecolityper --check
```

## 🚀 Usage

```bash
# Basic Usage
# Single genome analysis
ecolityper -i genome.fna -o results

# Multiple genomes with wildcards
ecolityper -i "*.fna" -o results --threads 8

# Multiple specific files
ecolityper -i genome1.fna genome2.fna genome3.fna -o results

# Advanced Options
# High-performance with 16 threads
ecolityper -i "data/*.fasta" -o analysis_results --threads 16

# Check environment and exit
ecolityper --check

# Show version information
ecolityper --version
```

### Command Line Options

- `-i, --inputs`      Input genome FASTA files (supports globs, e.g. '*.fasta')
- `-o, --outdir`      Output directory (default: ecolityper_results)
- `--threads`         Number of parallel workers (default: CPU count)
- `--check`           Check environment and exit
- `--version`         Print version banner and exit

## 📊 Output Files

EcoliTyper generates comprehensive output in multiple formats:

### Main Output Files

- `ecolityper_summary.tsv` - Combined results for all samples
- `mlst_results.tsv` - Detailed MLST results
- `serotype_results.tsv` - O and H antigen typing results
- `clermont_results.tsv` - Phylotyping results with method
- `ecolityper_run_meta.json` - Run metadata and tool versions

### Per-Sample Files

- `{sample}.ecolityper.json` - Complete results in JSON format
- `{sample}_serotype.json` - Raw SerotypeFinder JSON output

### Output Columns

| Sample    | MLST Scheme     | ST  | O-type | H-type | Clermont Phylotype | Method |
|-----------|-----------------|-----|--------|--------|-------------------|--------|
| ecoli.fna | ecoli_achtman_4 | 156 | O173   | H28    | B1                | PCR    |

## 🏗️ Technical Details

### Performance

- **Parallel Processing**: Utilizes Python's ThreadPoolExecutor for efficient multi-core usage
- **Memory Efficient**: Processes samples sequentially with minimal memory footprint
- **Fast Execution**: Typical analysis time: 10-15 seconds per genome

### Typing Methods

- **MLST**: Uses the standard Achtman 7-gene scheme for E. coli
- **Serotyping**: Leverages CGE's highly curated SerotypeFinder database with BLAST-based identification
- **Clermont Typing**: Implements the EzClermont tool for rapid phylogroup determination

### Database Integration

- **MLST**: Bundled MLST tool with comprehensive allele database
- **Serotyping**: Integrated CGE SerotypeFinder with regularly updated O and H antigen databases
- **Clermont**: EzClermont Python package for consistent phylogroup calling

## 🎯 Use Cases

### 🏥 Clinical Microbiology

- Rapid strain characterization for outbreak investigations
- Surveillance of antimicrobial resistant clones
- Virulence factor association studies

### 🔬 Research Applications

- Population genetics studies
- Evolutionary analysis of E. coli lineages
- Comparative genomics projects

### 🐄 Veterinary and Food Safety

- Source tracking in foodborne outbreaks
- Zoonotic transmission studies
- Agricultural surveillance

## 📝 Citation

If you use EcoliTyper in your research, please cite:

**Brown, B. (2025). EcoliTyper: Unified MLST + Serotyping + Clermont typing for Escherichia coli. Zenodo. https://doi.org/10.5281/zenodo.17226894**

```bibtex
@software{brown2025ecolityper,
  title = {EcoliTyper: Unified MLST + Serotyping + Clermont typing for Escherichia coli},
  author = {Brown, Beckley},
  year = {2025},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.17226894},
  url = {https://doi.org/10.5281/zenodo.17226894}
}
```

## 🤝 Contributing

We welcome contributions! Please feel free to submit pull requests, report bugs, or suggest new features.

### Development Setup

```bash
git clone https://github.com/bbeckley-hub/ecoliTyper.git
cd ecoliTyper
pip install -e .
```

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🙏 Acknowledgments

- CGE (Center for Genomic Epidemiology) for the excellent SerotypeFinder tool and database
- EzClermont developers for the reliable phylotyping implementation
- MLST developers for maintaining the comprehensive typing scheme
- The open-source bioinformatics community for invaluable tools and resources

## 📞 Contact

Beckley Brown  
📧 brownbeckley94@gmail.com  
🌐 GitHub Profile  
🔗 EcoliTyper Repository
