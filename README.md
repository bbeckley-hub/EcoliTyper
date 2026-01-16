
```
███████╗ ██████╗ ██████╗ ██╗     ██╗████████╗██╗   ██╗██████╗ ███████╗██████╗ 
██╔════╝██╔════╝██╔═══██╗██║     ██║╚══██╔══╝╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
█████╗  ██║     ██║   ██║██║     ██║   ██║    ╚████╔╝ ██████╔╝█████╗  ██████╔╝
██╔══╝  ██║     ██║   ██║██║     ██║   ██║     ╚██╔╝  ██╔═══╝ ██╔══╝  ██╔══██╗
███████╗╚██████╗╚██████╔╝███████╗██║   ██║      ██║   ██║     ███████╗██║  ██║
╚══════╝ ╚═════╝ ╚═════╝ ╚══════╝╚═╝   ╚═╝      ╚═╝   ╚═╝     ╚══════╝╚═╝  ╚═╝
```
</div>
# 🧬 EcoliTyper v1.0.0

**A species-optimized computational pipeline for comprehensive genotyping and surveillance of *Escherichia coli***

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8-3.14](https://img.shields.io/badge/Python-3.8--3.14-blue.svg)](https://www.python.org/)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](https://github.com/bbeckley-hub/EcoliTyper)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17761775.svg)](https://doi.org/10.5281/zenodo.17761775)
[![GitHub stars](https://img.shields.io/github/stars/bbeckley-hub/EcoliTyper)](https://github.com/bbeckley-hub/EcoliTyper/stargazers)

## 🌟 Overview

**EcoliTyper** is a revolutionary bioinformatics pipeline designed to eliminate workflow fragmentation in *E. coli* genomic surveillance. By integrating **seven core genotyping analyses** into a single automated workflow, EcoliTyper transforms disconnected genomic data into coherent biological narratives with actionable public health intelligence.

> *"From fragmented analysis to integrated insight in one command"*

### 🚀 The EcoliTyper Advantage

| Traditional Workflow 😫 | EcoliTyper Solution 🎉 |
|------------------------|-----------------------|
| 7+ independent tools required | **Single unified pipeline** |
| Manual data integration & synthesis | **Automated cross-genome pattern discovery** |
| Hours of manual curation | **Intelligent risk assessment & alerting** |
| Disconnected epidemiological context | **Integrated lineage database of high-risk clones** |
| Multiple output formats to reconcile | **Consolidated HTML report + structured data (TSV/JSON)** |
| Complex installation & dependencies | **Self-contained Conda package** |

**Key Achievement:** Processes 30 *E. coli* genomes in ~41 minutes on 16 CPU cores with perfect concordance against reference tools.

## ✨ Core Features

### 🧩 **Comprehensive *E. coli* Typing Suite**
- **🧬 Multi-Locus Sequence Typing (MLST)** - Achtman scheme with PubMedST database
- **🔍 In silico Serotyping** - O and H antigen determination via SerotypeFinder (≥90% coverage/identity)
- **🎯 CH Typing** - High-resolution *fumC/fimH* typing for fine-scale discrimination
- **🌳 Clermont Phylogrouping** - Evolutionary context with 2013 scheme (8 phylogroups)
- **💊 Antimicrobial Resistance Profiling** - Dual screening via ABRicate (9 databases) & NCBI-AMRFinderPlus
- **🦠 Virulence Factor Detection** - Comprehensive pathogenicity assessment
- **📊 Plasmid Replicon Typing** - Mobile genetic element characterization

### 🧠 **Intelligent Analytics Layer**
- **🔬 Cross-genome pattern discovery** - Automated gene frequency analysis & distribution mapping
- **⚠️ Rule-based clinical risk assessment** - Hierarchical alerting (CARBAPENEMASE > ESBL > COLISTIN-RES)
- **🌍 Integrated lineage database** - Manually curated reference of high-risk clones (ST131, ST1193, etc.)
- **📈 Population-level insights** - Immediate epidemiological overview of resistance cassettes & virulence profiles

### ⚡ **Performance Optimized Architecture**
- **🚀 Hybrid parallel execution** - Inter-module & intra-module parallelization
- **🎛️ Dynamic resource allocation** - Automatic scaling with genome complexity
- **⚖️ Memory-aware processing** - Strategic sequential execution for resource-intensive operations
- **🔄 Robust error handling** - Graceful recovery with checkpointing & automated cleanup

## 🛠️ Installation

### Quick Install (Recommended)
```bash
# Create and activate environment
conda create -n ecolityper-c conda-forge -c bioconda  -c bbeckley-hub ecolityper -y
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

## 🎯 Usage Examples

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
```bash
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

## 📊 Output Structure

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

## 🎨 Interactive Report Features

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

## 🔗 Integrated Databases & Resources

### **Core Typing Databases**
| Database | Purpose | Source | Version |
|----------|---------|--------|---------|
| **PubMedST** | MLST typing | https://pubmlst.org/ | Latest |
| **SerotypeFinder DB** | O:H antigen typing | https://bitbucket.org/genomicepidemiology/serotypefinder_db/ | 2.0.1 |
| **CH Typing DB** | *fumC/fimH* typing | https://bitbucket.org/genomicepidemiology/chtyper_db/ | Latest |
| **Clermont Scheme** | Phylogrouping | Integrated via ezClermont | 2013 |

### **AMR & Virulence Databases (via ABRicate)**
| Database | Purpose | Coverage |
|----------|---------|----------|
| **CARD** | Comprehensive antibiotic resistance | 5,000+ genes |
| **ResFinder** | Acquired antimicrobial resistance | 3,000+ genes |
| **NCBI** | NCBI bacterial AMR reference | 2,500+ genes |
| **ARG-ANNOT** | Antibiotic resistance gene annotation | 2,000+ genes |
| **MEGARES** | Comprehensive resistance database | 8,000+ genes |
| **VFDB** | Virulence factors | 2,500+ genes |
| **EcoH** | *E. coli* hemolysins | 100+ genes |
| **Ecoli_VF** | *E. coli* virulence factors | 500+ genes |
| **PlasmidFinder** | Plasmid replicons | 500+ types |

### **NCBI AMRFinderPlus Database**
- **Coverage**: 6,000+ resistance genes and variants
- **Update Frequency**: Weekly
- **Scope**: Both acquired genes and chromosomal mutations


## 🌍 EcoliDB Lineage Database

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

### **Future Development: AI-Powered Prediction**
We are developing machine learning and AI approaches to:
- Predict complete pattern combinations for rapid *E. coli* characterization
- Integrate results from EcoliTyper for enhanced predictive analytics
- Develop models for emerging resistance pattern prediction

**Follow our GitHub repository for upcoming releases:** https://github.com/bbeckley-hub/EcoliTyper

## ⚡ Performance Benchmarks

### **Processing Times**
| Scenario | Genomes | Time | Hardware | Speed per Genome |
|----------|---------|------|----------|------------------|
| Standard Workstation | 30 genomes | 80-150 min | 2 CPU cores, 8GB RAM | 3-6 min |
| High-Performance Server | 30 genomes | **41 min** | 16 CPU cores, 16GB RAM | **1.2 min** |
| Single Genome | 1 genome | 1-6 min | Variable | - |

### **Validation Accuracy**
- **100% concordance** with standalone reference tools (mlst, SerotypeFinder, ezClermont)
- **Perfect typing** of reference strains (K-12 MG1655, O157:H7, O18ac:H7)
- **Robust performance** across diverse clinical and reference isolates

## 🆚 Competitive Comparison

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
- **Mykrobe:** https://github.com/Mykrobe-tools/mykrobe
- **Bactopia:** https://github.com/bactopia/bactopia
- **ECTyper:** https://github.com/phac-nml/irida-plugin-ectyper

## 📚 Citation

If you use EcoliTyper in your research, please cite:

```bibtex
@software{beckley2025ecolityper,
  title = {EcoliTyper: A species-optimized computational pipeline for comprehensive genotyping and surveillance of Escherichia coli},
  author = {Beckley, Brown and Amarh, Vincent},
  year = {2025},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/bbeckley-hub/EcoliTyper}},
  doi = {10.5281/zenodo.17761775}
}
```

### **Third-Party Tool Citations**
EcoliTyper integrates several third-party tools. Please cite them when using corresponding modules:

```bibtex
% SerotypeFinder
@article{joensen2015rapid,
  title={Rapid and easy in silico serotyping of Escherichia coli using whole genome sequencing (WGS) data},
  author={Joensen, KG and others},
  journal={Journal of Clinical Microbiology},
  year={2015}
}

% MLST (mlst tool)
@article{larsen2012multilocus,
  title={Multilocus sequence typing of total genome sequenced bacteria},
  author={Larsen, M and others},
  journal={Journal of Clinical Microbiology},
  year={2012}
}

% ABRicate
@software{seemann2020abricate,
  title = {ABRicate: Mass screening of contigs for antimicrobial resistance or virulence genes},
  author = {Seemann, Torsten},
  year = {2020},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/tseemann/abricate}}
}

% AMRFinderPlus
@article{feldgarden2019validating,
  title={Validating the AMRFinder Tool and Resistance Gene Database by Using Antimicrobial Resistance Genotype-Phenotype Correlations in a Collection of Isolates},
  author={Feldgarden, M and others},
  journal={Antimicrobial Agents and Chemotherapy},
  year={2019}
}

% EzClermont
@article{waters2020easy,
  title={Easy phylotyping of Escherichia coli via the EzClermont web app and command-line tool},
  author={Waters, NR and others},
  journal={Access Microbiology},
  year={2020}
}

% CH Typing
@article{roer2018chtyper,
  title={CHTyper, a web tool for subtyping of extraintestinal pathogenic Escherichia coli},
  author={Roer, L and others},
  journal={Journal of Clinical Microbiology},
  year={2018}
}
```

## ❓ Frequently Asked Questions (FAQ)

### **General Questions**

#### Q: What makes EcoliTyper different from other typing tools?
**A:** EcoliTyper is specifically optimized for *E. coli* and integrates 7 complementary typing methods into a single pipeline with automated cross-genome pattern discovery and a curated lineage database for epidemiological context.

#### Q: Can I use EcoliTyper for other bacterial species?
**A:** No, EcoliTyper is specifically optimized for *Escherichia coli*. The algorithms, thresholds, and databases are tailored for this species. For other bacteria, consider generalist pipelines like Bactopia.

#### Q: What input formats are supported?
**A:** EcoliTyper requires assembled genomes in FASTA format (.fna, .fasta, .fa, .fsa). It does not process raw reads directly.

### **Installation & Setup**

#### Q: Can I install EcoliTyper without Conda?
**A:** While Conda is recommended for managing complex dependencies, you can install from source using pip. However, this requires manual installation of all dependencies. 

#### Q: How much disk space is required?
**A:** Approximately 5-10 GB for the Conda environment and databases. Additional space is needed for input genomes and output files.

### **Analysis & Results**

#### Q: How accurate is EcoliTyper compared to standalone tools?
**A:** EcoliTyper shows 100% concordance with standalone reference tools (mlst, SerotypeFinder, ezClermont) for standard typing methods on validated reference strains.

#### Q: Can I customize the analysis thresholds?
**A:** Currently, thresholds are optimized based on published recommendations (e.g., ≥90% coverage/identity for serotyping, ≥80% for gene detection). 

#### Q: How does the lineage database get updated?
**A:** The EcoliDB lineage database is manually curated and included in the package. Updates will be released with new versions of EcoliTyper. Users can contribute new lineages via GitHub issues or direct email.

#### Q: What should I do if I find a novel sequence type not in the database?
**A:** Please report it as a GitHub issue with supporting references. We actively maintain and expand the lineage database.

### **Performance & Troubleshooting**

#### Q: Why is my analysis taking longer than expected?
**A:** Processing time depends on genome size, complexity, and available resources. Large accessory genomes or many contigs increase processing time. 

#### Q: How do I interpret the risk assessment alerts?
**A:** Alerts follow a hierarchical system:
- **CARBAPENEMASE**: Highest priority (e.g., blaKPC, blaNDM)
- **ESBL**: Extended-spectrum β-lactamases (e.g., blaCTX-M)
- **COLISTIN-RES**: Plasmid-mediated colistin resistance (mcr genes)
- **OTHER CRITICAL**: Other high-priority resistance mechanisms

### **Collaboration & Contribution**

#### Q: How can I contribute to the lineage database?
**A:** Submit new lineages or updates via GitHub issues with:
1. Sequence type information
2. Epidemiological data (geography, prevalence)
3. Published references (PMID/DOI)
4. Resistance and virulence profiles

#### Q: Can I use EcoliTyper in clinical diagnostics?
**A:** EcoliTyper provides genotypic predictions. For clinical decision-making, confirmatory phenotypic testing is essential. Always follow local regulations and guidelines.


## 🤝 Contributing

We welcome contributions from the community! Here's how you can help:

### **Development Workflow**
1. 🍴 Fork the repository
2. 🌿 Create a feature branch (`git checkout -b feature/amazing-feature`)
3. 💾 Commit your changes (`git commit -m 'Add amazing feature'`)
4. 🚀 Push to the branch (`git push origin feature/amazing-feature`)
5. 🔔 Open a Pull Request

### **Areas for Contribution**
- Database expansion and curation
- Additional typing schemes
- Performance optimizations
- Visualization enhancements
- Documentation improvements

## 🐛 Issue Reporting

Found a bug? Have a feature request? Please let us know:

1. **Search existing issues** to avoid duplicates
2. **Create a new issue** with:
   - Clear descriptive title
   - Steps to reproduce (for bugs)
   - Expected vs actual behavior
   - Environment details (OS, Python version, etc.)

**Issue Categories:**
- 🐛 Bug Report
- 🚀 Feature Request
- 📚 Documentation
- 💡 Enhancement

## ⚠️ Limitations & Considerations

- **Species-specific:** Optimized exclusively for *E. coli*
- **Input requirement:** Requires assembled genomes (FASTA format)
- **Database dependency:** Accuracy depends on completeness of reference databases
- **Lineage coverage:** Curated database focuses on globally significant clones
- **Predictive nature:** Genotypic predictions require phenotypic confirmation for clinical decisions
- **Evolutionary dynamics:** Bacterial evolution continuously generates new variants; supplement with recent publications and local surveillance data

## 📜 License & Third-Party Components

### **EcoliTyper License**
This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

### **Third-Party Tool Licenses**
EcoliTyper integrates several third-party tools, each with their own licenses:

| Tool | License | Integration Type |
|------|---------|------------------|
| **mlst** | GPL-3.0 | Direct integration |
| **SerotypeFinder** | Apache 2.0 | Database integration |
| **ABRicate** | GPL-2.0 | Direct integration |
| **NCBI AMRFinderPlus** | Public Domain | Direct integration |
| **ezClermont** | MIT | Direct integration |
| **CH Typing databases** | Custom (research use) | Database integration |

### **Database Licenses**
- **PubMedST**: Free for academic use (https://pubmlst.org/)
- **CARD**: ODbL (Open Database License)
- **ResFinder**: Free for academic use
- **VFDB**: Free for academic use
- **PlasmidFinder**: Free for academic use

**Note:** For commercial use, verify license compatibility for all integrated components.

## 👥 Authors & Affiliations

### **Primary Authors**
- **Brown Beckley** - *Creator & Lead Developer*  
  Department of Medical Biochemistry, University of Ghana Medical School, Accra, Ghana  
  Department of Biochemistry and Biotechnology, Kwame Nkrumah University of Science and Technology, Kumasi, Ghana  
  📧 [brownbeckley94@gmail.com](mailto:brownbeckley94@gmail.com)

- **Dr. Vincent Amarh** - *Supervisor & Advisor*  
  Department of Medical Biochemistry, University of Ghana Medical School, Accra, Ghana

### **Correspondence**
**Brown Beckley**  
📧 [brownbeckley94@gmail.com](mailto:brownbeckley94@gmail.com)  
🌐 [GitHub Profile](https://github.com/bbeckley-hub)  
🔗 [EcoliTyper Repository](https://github.com/bbeckley-hub/EcoliTyper)

## 🙏 Acknowledgments

This work stands on the shoulders of the open-source bioinformatics community:

### **Tool Developers**
- Torsten Seemann (mlst, ABRicate)
- Center for Genomic Epidemiology (SerotypeFinder, CH typing databases)
- NCBI Pathogen Detection Team (AMRFinderPlus)
- EzClermont developers

### **Database Maintainers**
- PubMedST for MLST databases
- CARD, ResFinder, VFDB, PlasmidFinder database teams
- All contributors to the curated databases used in ABRicate

### **Support & Testing**
- Colleagues who provided invaluable feedback during development
- Public genome databases (NCBI RefSeq) for validation datasets
- Early adopters who tested and improved the pipeline

### **Funding**
This research received no specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

## 🔮 Future Development Roadmap

### **Short-term Goals (2025)**
- Regular Updates
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

## 📞 Support & Community

- **Documentation:** [GitHub Wiki](https://github.com/bbeckley-hub/EcoliTyper/wiki)
- **Questions:** GitHub Discussions or Issues
- **Collaborations:** Direct email contact
- **Updates:** Watch the GitHub repository for releases
- **Community:** Join the conversation on GitHub Discussions

---

<div align="center">

### **⭐ Star us on GitHub if you find EcoliTyper useful!**

*Transforming fragmented genomic surveillance into integrated public health intelligence* 🧬✨

**"From sequences to surveillance in one command"**

---
</div>

**Join the Fight Against Antimicrobial Resistance**

Antimicrobial resistance (AMR) represents one of the most significant global health threats of our time. We invite researchers, clinicians, and public health professionals to collaborate with us in:

- Expanding and validating our *E. coli* database
- Sharing regional epidemiological data
- Developing standardized typing methodologies
- Advancing AMR surveillance and intervention strategies

**Together, we can enhance global AMR monitoring and develop more effective treatment strategies.**


