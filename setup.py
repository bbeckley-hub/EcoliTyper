#!/usr/bin/env python3
"""
EcoliTyper Setup Configuration
Complete E. coli typing pipeline setup
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
"""

from setuptools import setup, find_packages
import pathlib

# Read the contents of README.md
here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="ecolityper",
    version="1.0.0",
    author="Brown Beckley",
    author_email="brownbeckley94@gmail.com",
    description="Comprehensive E. coli Typing Pipeline: MLST, Serotyping, CH Typing, Phylogrouping, AMR, and Virulence Analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bbeckley-hub/EcoliTyper",
    project_urls={
        "Bug Reports": "https://github.com/bbeckley-hub/EcoliTyper/issues",
        "Source": "https://github.com/bbeckley-hub/EcoliTyper",
        "Documentation": "https://github.com/bbeckley-hub/EcoliTyper",
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    keywords="bioinformatics, ecoli, typing, mlst, serotyping, amr, virulence, genomics",
    packages=find_packages(include=["ecoliTyper", "ecoliTyper.*"]),
    python_requires=">=3.8, <3.12",
    install_requires=[
        "pandas>=1.5.0",
        "biopython>=1.80", 
        "psutil>=5.9.0",
        "requests>=2.28.0",
        "tqdm>=4.64.0",
        "click>=8.0.0",
        "tabulate",
    ],
    entry_points={
        "console_scripts": [
            "ecolityper=ecoliTyper.ecolityper:main",
        ],
    },
    include_package_data=True,
    package_data={
        "ecoliTyper": [
            # MLST module
            "modules/mlst_module/*.py",
            "modules/mlst_module/db/**/*",
            "modules/mlst_module/bin/*",
            "modules/mlst_module/perl5/**/*",
            "modules/mlst_module/scripts/*",
            # Serotyping module
            "modules/serotypefinder_module/*.py", 
            "modules/serotypefinder_module/serotypefinder_db/**/*",
            # CH Typing module
            "modules/CHTyper_module/*.py",
            "modules/CHTyper_module/chtyper_db/**/*",
            # Other modules
            "modules/phylogrouping_module/*.py",
            "modules/Abricate_module/*.py",
            "modules/Amrfinder_module/*.py",
            "modules/Ecoli_lineage/*.py",
            # Core files
            "core/*.py",
        ],
    },
    zip_safe=False,
)
