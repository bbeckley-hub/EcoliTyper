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
    version="1.3.0",
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
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "pandas>=1.5.0",
        "biopython>=1.80", 
        "psutil>=5.9.0",
        "requests>=2.28.0",
        "tqdm>=4.64.0",
        "click>=8.0.0",
        "tabulate>=0.9.0",
        "ezclermont>=0.7.0",
        "cgecore>=1.5.6",
        # HTML parsing for summary/visualization
        "beautifulsoup4>=4.11.0",
        "lxml>=4.9.0",
        # Visualization
        "matplotlib>=3.5.0",
        "seaborn>=0.12.0",
        "scipy>=1.10.1",
    ],
    extras_require={
        'full': [
            "plotly>=5.10.0",
            "scipy>=1.9.0",
        ],
        'visualization': [
            "plotly>=5.10.0",
            "scipy>=1.9.0",
        ]
    },
    entry_points={
        "console_scripts": [
            "ecolityper=ecoliTyper.ecolityper:main",
        ],
    },
    include_package_data=True,
    package_data={
        '': ['**/*'],  # Include EVERYTHING recursively
    },
    zip_safe=False,
)