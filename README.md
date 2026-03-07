---
```markdown
███████╗ ██████╗ ██████╗ ██╗     ██╗████████╗██╗   ██╗██████╗ ███████╗██████╗ 
██╔════╝██╔════╝██╔═══██╗██║     ██║╚══██╔══╝╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
█████╗  ██║     ██║   ██║██║     ██║   ██║    ╚████╔╝ ██████╔╝█████╗  ██████╔╝
██╔══╝  ██║     ██║   ██║██║     ██║   ██║     ╚██╔╝  ██╔═══╝ ██╔══╝  ██╔══██╗
███████╗╚██████╗╚██████╔╝███████╗██║   ██║      ██║   ██║     ███████╗██║  ██║
╚══════╝ ╚═════╝ ╚═════╝ ╚══════╝╚═╝   ╚═╝      ╚═╝   ╚═╝     ╚══════╝╚═╝  ╚═╝
```

<div align="center">

# 🧬 **EcoliTyper v1.0.0**

### A species-optimized computational pipeline for comprehensive genotyping and surveillance of **_Escherichia coli_**

**Complete *E. coli* genomic analysis in minutes — not hours**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8-3.14](https://img.shields.io/badge/Python-3.8--3.14-blue.svg)](https://www.python.org/)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](https://github.com/bbeckley-hub/EcoliTyper)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17761775.svg)](https://doi.org/10.5281/zenodo.17761775)
[![GitHub stars](https://img.shields.io/github/stars/bbeckley-hub/EcoliTyper)](https://github.com/bbeckley-hub/EcoliTyper/stargazers)

**Perfect for clinical microbiology, outbreak investigations, and genomic research.**

</div>

---

## 📋 **Table of Contents**

- [🌟 Overview](#-overview)
- [✨ Core Features](#-core-features)
- [🛠️ Installation](#️-installation)
- [🎯 Usage Examples](#-usage-examples)
- [📊 Output Structure](#-output-structure)
- [🎨 Interactive Report Features](#-interactive-report-features)
- [🔗 Integrated External Tools & Dependencies](#-integrated-external-tools--dependencies)
- [🌍 EcoliDB Lineage Database](#-ecolityper-ecolidb-lineage-database)
- [⚡ Performance Benchmarks](#-performance-benchmarks)
- [🆚 Competitive Comparison](#-competitive-comparison)
- [📚 Citation](#-citation)
- [❓ Frequently Asked Questions](#-frequently-asked-questions)
- [🤝 Contributing](#-contributing)
- [🐛 Issue Reporting](#-issue-reporting)
- [⚠️ Limitations & Considerations](#️-limitations--considerations)
- [📜 License & Third-Party Components](#-license--third-party-components)
- [👥 Authors & Affiliations](#-authors--affiliations)
- [🙏 Acknowledgements](#-acknowledgements)
- [🔮 Future Development Roadmap](#-future-development-roadmap)
- [📞 Support & Community](#-support--community)

---

## 🌟 **Overview**

**EcoliTyper** is a revolutionary bioinformatics pipeline designed to eliminate workflow fragmentation in *E. coli* genomic surveillance. By integrating **seven core genotyping analyses** into a single automated workflow, EcoliTyper transforms disconnected genomic data into coherent biological narratives with actionable public health intelligence.

> *"From fragmented analysis to integrated insight in one command"*

### 🚀 **The EcoliTyper Advantage**

| Traditional Workflow 😫 | EcoliTyper Solution 🎉 |
|------------------------|-----------------------|
| 7+ independent tools required | **Single unified pipeline** |
| Manual data integration & synthesis | **Automated cross-genome pattern discovery** |
| Hours of manual curation | **Intelligent risk assessment & alerting** |
| Disconnected epidemiological context | **Integrated lineage database of high-risk clones** |
| Multiple output formats to reconcile | **Consolidated HTML report + structured data (TSV/JSON)** |
| Complex installation & dependencies | **Self-contained Conda package** |

**Key Achievement:** Processes 30 *E. coli* genomes in **~41 minutes** on 16 CPU cores with **perfect concordance** against reference tools.

---

## ✨ **Core Features**

### 🧩 **Comprehensive *E. coli* Typing Suite**
- **🧬 Multi-Locus Sequence Typing (MLST)** – Achtman scheme with PubMedST database
- **🔍 In silico Serotyping** – O and H antigen determination via SerotypeFinder (≥90% coverage/identity)
- **🎯 CH Typing** – High-resolution *fumC/fimH* typing for fine-scale discrimination
- **🌳 Clermont Phylogrouping** – Evolutionary context with 2013 scheme (8 phylogroups)
- **💊 Antimicrobial Resistance Profiling** – Dual screening via ABRicate (9 databases) & NCBI-AMRFinderPlus
- **🦠 Virulence Factor Detection** – Comprehensive pathogenicity assessment
- **📊 Plasmid Replicon Typing** – Mobile genetic element characterization

### 🧠 **Intelligent Analytics Layer**
- **🔬 Cross-genome pattern discovery** – Automated gene frequency analysis & distribution mapping
- **⚠️ Rule-based clinical risk assessment** – Hierarchical alerting (CARBAPENEMASE > ESBL > COLISTIN-RES)
- **🌍 Integrated lineage database** – Manually curated reference of high-risk clones (ST131, ST1193, etc.)
- **📈 Population-level insights** – Immediate epidemiological overview of resistance cassettes & virulence profiles

### ⚡ **Performance Optimized Architecture**
- **🚀 Hybrid parallel execution** – Inter-module & intra-module parallelization
- **🎛️ Dynamic resource allocation** – Automatic scaling with genome complexity
- **⚖️ Memory-aware processing** – Strategic sequential execution for resource-intensive operations
- **🔄 Robust error handling** – Graceful recovery with checkpointing & automated cleanup

---

## 🛠️ **Installation**

### Quick Install (Recommended)
```bash
# Create and activate environment
conda create -n ecolityper -c conda-forge -c bioconda -c bbeckley-hub ecolityper -y
conda activate ecolityper
```

### From Source
```bash
git clone https://github.com/bbeckley-hub/EcoliTyper.git
cd EcoliTyper
conda env create -f environment.yml
conda activate ecolityper
pip install -e .
```

### System Requirements
- **Minimum:** 2 CPU cores, 8 GB RAM
- **Recommended:** 8+ CPU cores, 16+ GB RAM for batch processing
- **OS:** Linux, macOS, or Windows (WSL2 recommended for Windows)

---

## 🎯 **Usage Examples**

### Basic Single Genome Analysis
```bash
ecolityper -i genome.fasta -o results_directory/
```

### High-Throughput Batch Processing
```bash
# Process all FASTA files in current directory
ecolityper -i "*.fasta" -o batch_results --threads 8

# Process specific pattern
ecolityper -i "GCF_*.fna" -o surveillance_run --threads 16
```

### Customized Analysis Workflows
```bash
# Skip specific modules for faster processing
ecolityper -i isolates/ -o quick_typing --skip-amrfinder --skip-visualization

# Minimum typing only
ecolityper -i sample.fna -o basic_results --skip-lineage --skip-summary
```

### Complete Command Reference
```
usage: ecolityper [-h] -i INPUT -o OUTPUT [-t THREADS] [--skip-amrfinder]
                  [--skip-abricate] [--skip-mlst] [--skip-serotyping]
                  [--skip-chtyper] [--skip-phylogrouping] [--skip-lineage]
                  [--skip-summary] [--skip-visualization]

EcoliTyper: Complete E. coli Typing Pipeline

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTA file(s) - can use glob patterns like "*.fna" or "*.fasta"
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
  --skip-summary        Skip summary report generation
  --skip-visualization  Skip visualization generation

Examples:
  ecolityper -i genome.fna -o results/
  ecolityper -i "*.fna" -o batch_results --threads 8
  ecolityper -i "*.fasta" -o analysis --threads 16 --skip-lineage
  ecolityper -i "genome*.fa" -o results/ --threads 4

Supported FASTA formats: .fna, .fasta, .fa, .fsa

Analysis Modules:
  • MLST (Multi-Locus Sequence Typing)
  • Serotyping (O and H antigen determination)
  • CH Typing (FumC and FimH typing)
  • Phylogrouping (Clermont algorithm)
  • ABRicate (Resistance/Virulence/Plasmid screening)
  • AMRfinderPlus (NCBI AMR gene detection)
  • Lineage reference database
  • Summary Reports (HTML summary reports)
  • Visualizations (Charts and visualizations)

Output: Comprehensive results for all analyses in organized directories
```

---

## 📊 **Output Structure**

```
results_directory/
├── 📄 abricate_results/              # Multi-database screening (CARD, ResFinder, VFDB, etc.)
│   ├── ecoli_*_summary.json         # Consolidated JSON summaries
│   ├── ecoli_*_summary_report.html  # Interactive HTML reports
│   └── per_sample/                  # Individual genome results
├── 🔬 amrfinder_results/             # NCBI AMRFinderPlus outputs
│   ├── ecoli_amrfinder_summary.tsv
│   ├── ecoli_amrfinder_summary_report.html
│   └── per_sample/
├── 🎯 chtyper_results/               # High-resolution CH typing
│   ├── chtyper_results.tsv
│   ├── chtyper_results.html
│   └── per_sample/
├── 🧬 mlst_results/                  # Multi-Locus Sequence Typing
│   ├── mlst_summary.tsv
│   ├── mlst_summary.html
│   └── per_sample/
├── 🌳 phylogrouping_results/         # Clermont phylogrouping
│   ├── phylogrouping_results.tsv
│   ├── phylogrouping_results.html
│   └── per_sample/
├── 🔍 serotyping_results/            # O:H antigen typing
│   ├── serotype_analysis_report.tsv
│   ├── serotype_analysis_report.html
│   └── per_sample/
├── 🌍 lineage_results/               # Epidemiological context
│   └── ecoli_comprehensive_reference.html
├── 📈 summary_results/               # Consolidated reports
│   └── GENIUS_ULTIMATE_REPORTS/
│       ├── genius_ultimate_report.html     # Main interactive report
│       ├── genius_ultimate_report.json
│       ├── amr_genes.csv
│       ├── virulence_genes.csv
│       └── pattern_discovery.csv
└── 🎨 visualization_results/         # Publication-ready figures
    └── ECOLI_VISUALIZATIONS/
        ├── PDF/     # Vector graphics
        ├── PNG/     # Raster images
        ├── SVG/     # Scalable vector graphics
        └── DATA/    # Source data for figures
```

---

## 🎨 **Interactive Report Features**

### **Main Dashboard**
- **Sample Overview**: Quick glance at typing results across all genomes
- **Risk Alert Panel**: Automatic flagging of high-priority resistance markers
- **Epidemiological Context**: Lineage information for identified clones

### **Cross-Genome Analysis**
- **Gene Frequency Tables**: Prevalence of AMR/virulence genes across population
- **Pattern Discovery**: Identification of common resistance cassettes
- **Distribution Maps**: Visual representation of gene carriage

### **Visualization Gallery**
- **Stacked Bar Charts**: MLST, serotype, and phylogroup distributions
- **Violin Plots**: Quantitative metrics distribution
- **Pie Charts**: Phylogroup and serotype proportions
- **Heatmaps**: Gene presence/absence patterns

---

## 🔗 **Integrated External Tools & Dependencies**

EcoliTyper integrates several powerful open-source tools and databases. These are **not bundled directly in this repository**. Instead, they are automatically installed as **dependencies via Conda** (as defined in `environment.yml`). The MIT license that applies to the EcoliTyper pipeline code does not cover these external tools. Each tool is used under the terms of its own license, and we gratefully acknowledge their authors.

| Tool/Database | Purpose | Source | License |
|---------------|---------|--------|---------|
| **MLST** | Multi-locus sequence typing | [tseemann/mlst](https://github.com/tseemann/mlst) | GPL v2 |
| **ABRicate** | Mass screening for resistance/virulence | [tseemann/abricate](https://github.com/tseemann/abricate) | GPL v2 |
| **AMRFinderPlus** | AMR gene detection | [ncbi/amr](https://github.com/ncbi/amr) | Public Domain |
| **SerotypeFinder** | O:H antigen typing | [CGE](https://bitbucket.org/genomicepidemiology/serotypefinder_db/) | Apache 2.0 |
| **CHTyper DB** | *fumC/fimH* typing | [CGE](https://bitbucket.org/genomicepidemiology/chtyper_db/) | Free for research |
| **ezClermont** | Phylogrouping | [https://github.com/nickp60/ezClermont](https://github.com/nickp60/ezClermont) | MIT |

### **AMR & Virulence Databases (via ABRicate)**
| Database | Purpose | License |
|----------|---------|---------|
| **CARD** | Comprehensive antibiotic resistance | ODbL |
| **ResFinder** | Acquired antimicrobial resistance | Free for research |
| **NCBI** | NCBI bacterial AMR reference | Public Domain |
| **ARG-ANNOT** | Antibiotic resistance gene annotation | Free for research |
| **MEGARES** | Comprehensive resistance database | Free for research |
| **VFDB** | Virulence factors | Free for research |
| **EcoH** | *E. coli* hemolysins | Free for research |
| **Ecoli_VF** | *E. coli* virulence factors | Free for research |
| **PlasmidFinder** | Plasmid replicons | Free for research |

---

## 🌍 **EcoliDB Lineage Database** (EcoliTyper)

### **Overview**
EcoliTyper includes **EcoliDB**, a manually curated comprehensive reference database for rapid *E. coli* lineage contextualization. This database associates sequence types with clinical pathotypes, serotypes, and risk profiles to inform public health analysis.

### **Database Statistics**
- **12 Sequence Types** with detailed epidemiological profiles
- **13 Pathotypes** categorized (Diarrheagenic, Extraintestinal, Hybrid, Animal, Mucosal)
- **13 Serotypes** with clinical associations
- **8 Phylogroups** according to Clermont scheme
- **4 Carbapenemase Types** for resistance profiling
- **79 Scientific References** supporting the data

### **Included High-Risk Clones**
| Sequence Type | Risk Level | Primary Pathotype | Key Features |
|--------------|------------|-------------------|--------------|
| **ST131** | VERY HIGH | UPEC/ExPEC | Global MDR pandemic clone, CTX-M-15, fluoroquinolone resistance |
| **ST1193** | HIGH | UPEC/ExPEC | Emerging fluoroquinolone-resistant, community-associated UTIs |
| **ST95** | VERY HIGH | NMEC/ExPEC | Neonatal meningitis, high virulence, O18:H7 serotype |
| **ST405** | VERY HIGH | ExPEC | Global MDR, carbapenemase producers (OXA-48, NDM) |
| **ST410** | VERY HIGH | ExPEC | Emerging MDR, OXA-181/NDM-5 carbapenemases |
| **ST648** | VERY HIGH | Zoonotic MDR | Pan-drug resistance emerging, significant One Health concern |
| **ST11** | VERY HIGH | EHEC | O157:H7, hemorrhagic colitis, HUS risk |
| **ST10** | LOW-MODERATE | Commensal/Pathogenic | Diverse genetic background for horizontal gene transfer |
| **ST117** | MODERATE | APEC | Avian pathogenic, poultry industry concern |
| **ST69** | HIGH | Hybrid UPEC/EAEC | Uropathogenic/diarrheagenic hybrid |
| **ST73** | HIGH | Classic UPEC | Community-associated UTIs, high virulence |
| **ST88** | HIGH | NMEC/ExPEC | Meningitis-associated, less common than ST95 |

### **Accessing the Lineage Database**
The lineage database is automatically generated during analysis and can be found at:
```
lineage_results/ecoli_comprehensive_reference.html
```

This interactive HTML file provides:
- **Search functionality** by sequence type, serotype, or resistance profile
- **Risk categorization** (HIGH, MODERATE, LOW)
- **Geographical distribution** maps
- **Treatment recommendations** based on resistance profiles
- **Key references** for each lineage

---

## ⚡ **Performance Benchmarks**

| Scenario | Genomes | Time | Hardware | Speed per Genome |
|----------|---------|------|----------|------------------|
| Standard Workstation | 30 genomes | 80-150 min | 2 CPU cores, 8GB RAM | 3-6 min |
| High-Performance Server | 30 genomes | **41 min** | 16 CPU cores, 16GB RAM | **1.2 min** |
| Single Genome | 1 genome | 1-6 min | Variable | - |

### **Validation Accuracy**
- **100% concordance** with standalone reference tools (mlst, SerotypeFinder, ezClermont)
- **Perfect typing** of reference strains (K-12 MG1655, O157:H7, O18ac:H7)
- **Robust performance** across diverse clinical and reference isolates

---

## 🆚 **Competitive Comparison**

| Feature | EcoliTyper | ECTyper | Bactopia | Mykrobe |
|---------|------------|---------|----------|---------|
| **Primary Focus** | *E. coli* integrated genotyping | *E. coli* serotyping | Multi-species generalist | AMR prediction |
| **MLST** | ✅ Achtman scheme | ❌ | ✅ | ❌ |
| **Serotyping** | ✅ O:H (SerotypeFinder) | ✅ | Limited | ❌ |
| **CH Typing** | ✅ *fumC/fimH* | ❌ | ❌ | ❌ |
| **Clermont Phylogrouping** | ✅ 2013 scheme | ❌ | ✅ | ❌ |
| **AMR Profiling** | ✅ ABRicate + AMRFinderPlus | Limited | ✅ AMRFinder | ✅ Core function |
| **Virulence Screening** | ✅ 9 databases | Shiga toxins only | Limited | ❌ |
| **Cross-genome Analysis** | ✅ Automated pattern discovery | ❌ | ❌ | ❌ |
| **Lineage Database** | ✅ Curated high-risk clones | ❌ | ❌ | ❌ |
| **Output Formats** | HTML, TSV, JSON, text | Various | Various | Various |
| **Installation** | ⚡ Single Conda package | Moderate | Complex (Nextflow) | Simple |
| **Typing Speed (30 genomes)** | **41 minutes** | N/A | ~120 minutes | N/A |

**Reference Tools:**
- **Mykrobe:** [https://github.com/Mykrobe-tools/mykrobe](https://github.com/Mykrobe-tools/mykrobe)
- **Bactopia:** [https://github.com/bactopia/bactopia](https://github.com/bactopia/bactopia)
- **ECTyper:** [https://github.com/phac-nml/irida-plugin-ectyper](https://github.com/phac-nml/irida-plugin-ectyper)

---

## 📚 **Citation**

If you use EcoliTyper in your research, please cite:

```bibtex
@software{beckley2025ecolityper,
  title = {EcoliTyper: A species-optimized computational pipeline for comprehensive genotyping and surveillance of Escherichia coli},
  author = {Beckley, B. and Amarh, V.},
  year = {2025},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/bbeckley-hub/EcoliTyper}},
  doi = {10.5281/zenodo.17761775}
}
```

### **Third-Party Tool Citations**
EcoliTyper integrates several third-party tools. Please cite them when using corresponding modules:

**Serotyping & CH Typing**
```bibtex
@article{joensen2015rapid,
  author = {Joensen, K. G. et al.},
  title = {Rapid and easy in silico serotyping of Escherichia coli using whole genome sequencing data},
  journal = {Journal of Clinical Microbiology},
  year = {2015}
}

@article{roer2018chtyper,
  author = {Roer, L. et al.},
  title = {CHTyper, a web tool for subtyping of extraintestinal pathogenic Escherichia coli},
  journal = {Journal of Clinical Microbiology},
  year = {2018}
}

**ABRicate (Torsten Seemann)**
```bibtex

@software{seemann_abricate_2018,
  author = {Seemann, T.},
  title = {ABRicate: Mass screening of contigs for antimicrobial resistance and virulence genes},
  year = {2028},
  publisher = {GitHub},
  url = {https://github.com/tseemann/abricate}
}
**AMR**
```bibtex
@article{feldgarden2019validating,
  author = {Feldgarden, M. et al.},
  title = {Validating the AMRFinder Tool and Resistance Gene Database},
  journal = {Antimicrobial Agents and Chemotherapy},
  year = {2019}
}
```

**Phylogrouping**
```bibtex
@article{waters2020easy,
  author = {Waters, N. R. et al.},
  title = {Easy phylotyping of Escherichia coli via the EzClermont web app},
  journal = {Access Microbiology},
  year = {2020}
}
```

---

## ❓ **Frequently Asked Questions**

### **General Questions**

**Q: What makes EcoliTyper different from other typing tools?**
A: EcoliTyper is specifically optimized for *E. coli* and integrates 7 complementary typing methods into a single pipeline with automated cross-genome pattern discovery and a curated lineage database for epidemiological context.

**Q: Can I use EcoliTyper for other bacterial species?**
A: No, EcoliTyper is specifically optimized for *Escherichia coli*. The algorithms, thresholds, and databases are tailored for this species.

### **Installation & Setup**

**Q: How much disk space is required?**
A: Approximately 5-10 GB for the Conda environment and databases. Additional space is needed for input genomes and output files.

### **Analysis & Results**

**Q: How accurate is EcoliTyper compared to standalone tools?**
A: EcoliTyper shows **100% concordance** with standalone reference tools (mlst, SerotypeFinder, ezClermont) for standard typing methods on validated reference strains.

**Q: What should I do if I find a novel sequence type not in the database?**
A: Please report it as a GitHub issue with supporting references. We actively maintain and expand the lineage database.

---

## 🤝 **Contributing**

We welcome contributions from the community! Here's how you can help:

1. 🍴 Fork the repository
2. 🌿 Create a feature branch (`git checkout -b feature/amazing-feature`)
3. 💾 Commit your changes (`git commit -m 'Add amazing feature'`)
4. 🚀 Push to the branch (`git push origin feature/amazing-feature`)
5. 🔔 Open a Pull Request

**Areas for Contribution:**
- Database expansion and curation
- Additional typing schemes
- Performance optimizations
- Visualization enhancements
- Documentation improvements

---

## 📜 **License & Third-Party Components**

### **EcoliTyper Core Code**
The EcoliTyper pipeline code (the workflow engine, report generation, HTML templates, and Python modules written by the authors) is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.

### **Third-Party Tool Licenses**
EcoliTyper executes several external bioinformatics tools, which are installed as Conda dependencies. Each tool is the property of its respective developers and is used under its own license:

| Tool | License |
|------|---------|
| **MLST** (Torsten Seemann) | GPL v2 |
| **ABRicate** (Torsten Seemann) | GPL v2 |
| **AMRFinderPlus** (NCBI) | Public Domain |
| **SerotypeFinder** (CGE) | Apache 2.0 |
| **CH Typing databases** (CGE) | Free for research |
| **ezClermont** | MIT |

By using EcoliTyper, you agree to comply with the licenses of these third-party tools and databases.

---

## 👥 **Authors & Affiliations**

### **Primary Authors**
- **Brown Beckley** – *Creator & Lead Developer*
  Department of Medical Biochemistry, University of Ghana Medical School, Accra, Ghana
  Department of Biochemistry and Biotechnology, Kwame Nkrumah University of Science and Technology, Kumasi, Ghana
  📧 [brownbeckley94@gmail.com](mailto:brownbeckley94@gmail.com)

- **Dr. Vincent Amarh** – *Supervisor & Advisor*
  Department of Medical Biochemistry, University of Ghana Medical School, Accra, Ghana

---

## 🔮 **Future Development Roadmap**

### **Short-term Goals (2025)**
- Regular database updates
- Enhanced visualization capabilities
- Improved documentation and tutorials

### **Medium-term Goals (2026)**
- Integration with raw read analysis pipelines
- Real-time database update mechanisms
- Cloud deployment options (Docker, Singularity)

### **Long-term Vision**
- AI/ML models for predictive analytics
- Web interface for non-command-line users
- Expanded lineage database with global collaborations
- Integration with public health surveillance systems

---

<div align="center">

## **⭐ Star us on GitHub if you find EcoliTyper useful!**

*Transforming fragmented genomic surveillance into integrated public health intelligence* 🧬✨

**"From sequences to surveillance in one command"**

---

**Join the Fight Against Antimicrobial Resistance**

Antimicrobial resistance (AMR) represents one of the most significant global health threats of our time. We invite researchers, clinicians, and public health professionals to collaborate with us in expanding and validating our *E. coli* database, sharing regional epidemiological data, and advancing AMR surveillance.

**Together, we can enhance global AMR monitoring and develop more effective treatment strategies.**

</div>

