
```bash
███████╗ ██████╗ ██████╗ ██╗     ██╗████████╗██╗   ██╗██████╗ ███████╗██████╗ 
██╔════╝██╔════╝██╔═══██╗██║     ██║╚══██╔══╝╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
█████╗  ██║     ██║   ██║██║     ██║   ██║    ╚████╔╝ ██████╔╝█████╗  ██████╔╝
██╔══╝  ██║     ██║   ██║██║     ██║   ██║     ╚██╔╝  ██╔═══╝ ██╔══╝  ██╔══██╗
███████╗╚██████╗╚██████╔╝███████╗██║   ██║      ██║   ██║     ███████╗██║  ██║
╚══════╝ ╚═════╝ ╚═════╝ ╚══════╝╚═╝   ╚═╝      ╚═╝   ╚═╝     ╚══════╝╚═╝  ╚═╝
```

# 🧬 EcoliTyper

**A species-optimized computational pipeline for comprehensive genotyping and surveillance of *Escherichia coli***

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8-3.11](https://img.shields.io/badge/Python-3.8--3.11-blue.svg)](https://www.python.org/)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](https://github.com/bbeckley-hub/EcoliTyper)

## 🌟 Overview

**EcoliTyper** is a revolutionary bioinformatics pipeline that transforms fragmented *E. coli* genomic analysis into a unified, intelligent workflow. Say goodbye to juggling multiple disconnected tools! 🎯

> *"We shall not cease from exploration, and the end of all our exploring will be to arrive where we started and know the place for the first time."* - T.S. Eliot

### 🚀 Why EcoliTyper?

| Before EcoliTyper 😫 | With EcoliTyper 🎉 |
|---------------------|-------------------|
| 7+ different tools to run manually | **Single command** does it all |
| Manual data synthesis & cross-referencing | **Automated pattern discovery** & risk assessment |
| Installation conflicts & dependency hell | **Self-contained Conda package** |
| Hours of manual analysis | **Intelligent insights** in minutes |
| Fragmented, disconnected results | **Coherent biological narratives** |

## ✨ Key Features

### 🧩 Comprehensive Integration
- **🧬 MLST Typing** - Multi-Locus Sequence Typing
- **🔍 Serotyping** - O and H antigen determination
- **🎯 CH Typing** - High-resolution *fumC/fimH* typing
- **🌳 Phylogrouping** - Clermont 2013 scheme
- **💊 AMR Profiling** - Comprehensive resistance gene detection
- **🦠 Virulence Analysis** - Pathogenicity assessment
- **📊 Lineage Context** - Curated database of high-risk clones

### 🧠 Intelligent Analytics
- **🔬 Cross-genome pattern discovery**
- **⚠️ Automated clinical risk assessment**
- **📈 Gene frequency analysis**
- **🎚️ Tiered alert system** (CRITICAL, WARNING)
- **🌍 Epidemiological contextualization**

### ⚡ Performance Optimized
- **🚀 Parallel execution architecture**
- **🎛️ Dynamic resource allocation**
- **📦 Self-contained dependency management**
- **⚖️ Memory-aware processing**

## 🛠️ Installation

### Prerequisites
- **Conda** (Miniconda or Anaconda)

### Quick Install (Recommended)
```bash
conda create -n ecolityper -c bbeckley-hub -c conda-forge ecolityper -y
conda activate ecolityper
```

### Development Install
```bash
git clone https://github.com/bbeckley-hub/EcoliTyper.git
cd EcoliTyper
pip install -e .
```

### Database Setup
```bash
# Update AMRfinderPlus databases
amrfinder -u

# Setup ABRicate databases
abricate --setupdb
```

## 🎯 Usage

### Basic Single Genome Analysis
```bash
ecolityper -i genome.fna -o results/
```

### High-Throughput Batch Processing
```bash
ecolityper -i "*.fna" -o batch_results --threads 8
```

### Customized Analysis (Skip Specific Modules)
```bash
ecolityper -i "*.fasta" -o analysis --threads 16 --skip-lineage
```

### Complete Command Reference
```bash
usage: ecolityper [-h] -i INPUT -o OUTPUT [-t THREADS] [--skip-amrfinder]
                  [--skip-abricate] [--skip-mlst] [--skip-serotyping]
                  [--skip-chtyper] [--skip-phylogrouping] [--skip-lineage]

EcoliTyper: Complete E. coli Typing Pipeline

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTA file(s) - can use glob patterns like
                        "*.fna" or "*.fasta"
  -o OUTPUT, --output OUTPUT
                        Output directory for all results
  -t THREADS, --threads THREADS
                        Number of threads (default: 2)
  --skip-amrfinder      Skip AMRfinderPlus analysis
  --skip-abricate       Skip ABRicate analysis
  --skip-mlst           Skip MLST analysis
  --skip-serotyping     Skip serotyping analysis
  --skip-chtyper        Skip CH typing analysis
  --skip-phylogrouping  Skip phylogrouping analysis
  --skip-lineage        Skip lineage reference generation
```

## 📊 Output Structure

```
results/
├── 📄 individual_reports/          # Per-genome HTML reports
├── 📊 cross_genome_analysis/       # Population-level insights
├── 🔬 mlst_results/                # Sequence typing results
├── 🧪 serotype_results/            # O:H antigen typing
├️── 🎯 chtyper_results/            # fumC/fimH typing
├️── 🌳 phylogroup_results/         # Clermont phylogrouping
├️── 💊 amrfinder_results/          # NCBI AMR detection
├️── 🦠 abricate_results/           # Multi-database screening
├️── 📈 lineage_context/            # Epidemiological context
└️── 📋 summary_files/              # Consolidated TSV/JSON files
```

## 🎨 Example Reports

### 🔍 Cross-Genome Pattern Discovery
![Cross-genome Analysis](https://via.placeholder.com/800x400/4A90E2/FFFFFF?text=Interactive+Gene+Frequency+Analysis)

### 💊 Antimicrobial Resistance Profiling
![AMR Profile](https://via.placeholder.com/800x400/50E3C2/FFFFFF?text=Resistance+Gene+Distribution+Map)

### 🦠 Virulence Factor Analysis
![Virulence Profile](https://via.placeholder.com/800x400/B8E986/FFFFFF?text=Virulence+Factor+Prevalence)

## ⚡ Performance Benchmarks

| Scenario | Genomes | Time | Hardware |
|----------|---------|------|----------|
| 🐢 Standard Workstation | 30 genomes | ~80-150 min | 2 CPU cores, 8GB RAM |
| 🚀 High-Performance | 30 genomes | **~41 min** | 16 CPU cores, 16GB RAM |
| ⚡ Single Genome | 1 genome | 1-6 min | Variable |

## 📋 Supported Input Formats

- `.fna` ✅
- `.fasta` ✅  
- `.fa` ✅
- `.fsa` ✅

## 🆚 Feature Comparison

| Feature | EcoliTyper | ECTyper | Bactopia | Mykrobe |
|---------|------------|---------|----------|---------|
| **Species Focus** | 🎯 *E. coli* specific | *E. coli* specific | Multi-species | Multi-species |
| **MLST** | ✅ | ❌ | ✅ | ❌ |
| **Serotyping** | ✅ | ✅ | ✅ | ❌ |
| **CH Typing** | ✅ | ❌ | ❌ | ❌ |
| **Phylogrouping** | ✅ | ❌ | ✅ | ❌ |
| **Cross-genome Analysis** | ✅ | ❌ | ❌ | ❌ |
| **Lineage Database** | ✅ | ❌ | ❌ | ❌ |
| **Installation** | ⚡ Fast (minutes) | Fast | Slow (complex) | Fast |

## 🎓 Citation

If you use EcoliTyper in your research, please cite:

```bibtex
@software{beckley2024ecolityper,
  title = {EcoliTyper: A species-optimized computational pipeline for comprehensive genotyping and surveillance of Escherichia coli},
  author = {Beckley, Brown and Amarh, Vincent},
  year = {2025},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/bbeckley-hub/EcoliTyper}},
  doi = {10.5281/zenodo.17226894}
}
```

### 📚 Integrated Tool Citations

Please also cite these foundational tools:

```bibtex
% SerotypeFinder
@article{joensen2015rapid,
  title={Rapid and easy in silico serotyping of Escherichia coli using whole genome sequencing (WGS) data},
  author={Joensen, KG and others},
  journal={Journal of Clinical Microbiology},
  year={2015}
}

% MLST
@article{larsen2012multilocus,
  title={Multilocus sequence typing of total genome sequenced bacteria},
  author={Larsen, M and others},
  journal={Journal of Clinical Microbiology},
  year={2012}
}

% KMA
@article{clausen2018rapid,
  title={Rapid and precise alignment of raw reads against redundant databases with KMA},
  author={Clausen, P and others},
  journal={BMC Bioinformatics},
  year={2018}
}

% CH Typing
@article{roer2018chtyper,
  title={CHTyper, a web tool for subtyping of extraintestinal pathogenic Escherichia coli},
  author={Roer, L and others},
  journal={Journal of Clinical Microbiology},
  year={2018}
}

% EzClermont
@article{waters2020easy,
  title={Easy phylotyping of Escherichia coli via the EzClermont web app and command-line tool},
  author={Waters, NR and others},
  journal={Access Microbiology},
  year={2020}
}
```

## 🤝 Contributing

We welcome contributions! Please feel free to submit issues, fork the repository, and create pull requests.

1. 🍴 Fork the repository
2. 🌿 Create a feature branch (`git checkout -b amazing-feature`)
3. 💾 Commit your changes (`git commit -m 'Add amazing feature'`)
4. 🚀 Push to the branch (`git push origin amazing-feature`)
5. 🔔 Open a Pull Request

## 🐛 Issue Reporting

Found a bug? Have a feature request? Please let us know by [creating an issue](https://github.com/bbeckley-hub/EcoliTyper/issues).

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👥 Authors

- **Brown Beckley** - *Creator & Lead Developer* - [GitHub](https://github.com/bbeckley-hub)
- **Dr. Vincent Amarh** - *Supervisor & Advisor*

## 🙏 Acknowledgments

- The open-source bioinformatics community 👐
- Developers of integrated foundational tools 🔧
- Center for Genomic Epidemiology & PubMLST 🗄️
- Colleagues who provided invaluable feedback during testing 🧪

## 📞 Contact

**Beckley Brown**  
📧 [brownbeckley94@gmail.com](mailto:brownbeckley94@gmail.com)  
🌐 [GitHub Profile](https://github.com/bbeckley-hub)  
🔗 [EcoliTyper Repository](https://github.com/bbeckley-hub/EcoliTyper)

---

<div align="center">

**⭐ Star us on GitHub if you find this tool useful!**

*Transforming fragmented genomic data into coherent biological narratives* 🧬✨

</div>
