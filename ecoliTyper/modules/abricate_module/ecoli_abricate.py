#!/usr/bin/env python3
"""
EcoliTyper ABRicate Standalone Module
Comprehensive ABRicate analysis for E. coli with merged per-sample outputs.

Outputs per genome (2 files):
  - {genome}_comprehensive_abricate_report.html
  - abricate_raw.tsv                    (all databases merged)

Batch-level outputs:
  - ecoli_{db}_abricate_summary.tsv     (per database)
  - ecoli_{db}_summary_report.html      (per database)
  - ecoli_abricate_master_summary.json  (master summary)

Author: Beckley Brown <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-07-23
Send a quick mail for any issues or further explanations.
"""

import subprocess
import sys
import os
import glob
import logging
import tempfile
import random
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any
import argparse
import re
from datetime import datetime
import psutil
import math
import json
from collections import defaultdict, Counter


class EcoliAbricateExecutor:
    """ABRicate executor for E. coli with merged per-sample reporting."""

    def __init__(self, cpus: int = None, min_identity: int = 80, min_coverage: int = 80):
        self.logger = self._setup_logging()
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)
        self.min_identity = min_identity
        self.min_coverage = min_coverage

        self.required_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ecoh', 'ecoli_vf', 'bacmet2'
        ]

        self.high_risk_genes = {
            'blaCTX-M', 'blaSHV', 'blaTEM', 'BlaAmpC1_Ecoli', 'BlaAmpC2_Ecoli',
            'BlaPenicillin_Binding_Protein_Ecoli', 'BlaampH_Ecoli', 'blaTEM-105',
            'blaTEM-150', 'blaTEM-135', 'blaCTX-M-14', 'blaCTX-M-1', 'blaTEM-1',
            '(Bla)blaTEM-150', '(Bla)blaCTX-M-14', '(Bla)blaTEM-135', '(Bla)blaCTX-M-1',
            'blaKPC', 'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP',
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            '(Col)mcr-1.1', 'MCR', 'CTX',
            'qnrA', 'qnrB', 'qnrC', 'qnrD', 'qnrS', 'qnrVC', 'aac(6\')-Ib-cr',
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'rmtG', 'rmtH',
            'AGlyaadA2', 'AGlyaph3-Ia', 'AGlystrA', 'AGlystrB', 'AGlyaac3-IId',
            'aadA2', 'aph3-Ia', 'strA', 'strB', 'aac3-IId', 'aadA1-pm', 'ant3-Ia',
            'APH3-Ib', 'APH3-Ia', 'APH6-Id', 'AAC3-IId', 'ANT3-Ia',
            'PhecatA2', 'PhefloR', 'catA2', 'floR', 'catII',
            'Sulsul2', 'Sulsul3', 'Sulsul1', 'sul2', 'sul3', 'sul1',
            'TettetA', 'TettetR', 'TettetC', 'tetA', 'tetR', 'tetC',
            'TmtdfrA14', 'TmtdfrA1', 'dfrA14', 'dfrA1',
            'cfr', 'optrA', 'poxtA', 'vat', 'vga', 'mphB',
        }

        self.virulence_genes = {
            'stx1', 'stx2', 'stx1A', 'stx1B', 'stx2A', 'stx2B', 'stxA', 'stx2d1A',
            'cnf1', 'hlyA', 'hlyB', 'hlyC', 'hlyD', 'astA', 'east1',
            'eltA', 'eltB', 'estIa', 'toxB',
            'eae', 'tir', 'paa', 'fdeC',
            'afaA', 'afaB', 'afaC', 'afaD', 'afaE',
            'papA', 'papB', 'papC', 'papD', 'papE', 'papF', 'papG', 'papH', 'papI', 'papJ', 'papK', 'papX',
            'sfaA', 'sfaB', 'sfaC', 'sfaD', 'sfaE', 'sfaF', 'sfaG', 'sfaH', 'sfaS', 'sfaX', 'sfaY',
            'fimA', 'fimB', 'fimC', 'fimD', 'fimE', 'fimF', 'fimG', 'fimH', 'fimI',
            'faeC', 'faeD', 'faeE', 'faeF', 'faeG', 'faeH', 'faeI', 'faeJ',
            'fanC', 'fanD', 'fanE', 'fanF', 'fanG', 'fanH',
            'fyuA', 'irp1', 'irp2', 'ybtA', 'ybtE', 'ybtP', 'ybtQ', 'ybtS', 'ybtT', 'ybtU', 'ybtX',
            'iucA', 'iucB', 'iucC', 'iucD', 'iutA',
            'chuA', 'chuS', 'chuT', 'chuU', 'chuV', 'chuW', 'chuX', 'chuY',
            'iroB', 'iroC', 'iroD', 'iroE', 'iroN',
            'shuA', 'shuS', 'shuT', 'shuX', 'shuY',
            'sitA', 'sitB', 'sitC', 'sitD',
            'escC', 'escD', 'escE', 'escF', 'escG', 'escI', 'escJ', 'escL', 'escN', 'escO',
            'escP', 'escR', 'escS', 'escT', 'escU', 'escV',
            'espA', 'espB', 'espD', 'espF', 'espG', 'espH', 'espJ', 'espK', 'espZ',
            'nleA', 'nleB1', 'nleB2', 'nleC', 'nleD', 'nleE', 'nleF', 'nleG', 'nleH1', 'nleH2',
            'map', 'espFu', 'tccP', 'tccP2',
            'cesA', 'cesB', 'cesD', 'cesD2', 'cesF', 'cesL', 'cesT',
            'kpsD', 'kpsM', 'kpsT',
            'ibeA', 'ibeB', 'ibeC',
            'vat', 'senB', 'cif', 'ospG',
            'iss2', 'ompA',
            'efa1', 'stcE', 'pic', 'sepD', 'sepL', 'sepQ', 'sepZ',
            'aslA', 'cdtA', 'cdtB', 'cdtC', 'cdiA', 'cdiB',
            'entA', 'entB', 'entC', 'entD', 'entE', 'entF', 'entS',
            'fepA', 'fepB', 'fepC', 'fepD', 'fepE', 'fepG', 'fes',
            'csgA', 'csgB', 'csgC', 'csgD', 'csgE', 'csgF', 'csgG',
            'gspC', 'gspD', 'gspE', 'gspF', 'gspG', 'gspH', 'gspI', 'gspJ', 'gspK', 'gspL', 'gspM',
            'fliC', 'fliD', 'fliE', 'fliF', 'fliG', 'fliH', 'fliI', 'fliJ', 'fliK', 'fliL', 'fliM',
            'fliN', 'fliO', 'fliP', 'fliQ', 'fliR', 'fliS', 'fliT', 'fliY', 'fliZ',
            'flgA', 'flgB', 'flgC', 'flgD', 'flgE', 'flgF', 'flgG', 'flgH', 'flgI', 'flgJ', 'flgK', 'flgL', 'flgN',
            'flhA', 'flhB', 'flhC', 'flhD', 'flhE',
            'motA', 'motB', 'cheA', 'cheB', 'cheR', 'cheW', 'cheY', 'cheZ',
        }

        self.critical_resistance_genes = {
            'blaCTX-M', 'blaSHV', 'blaTEM', 'blaKPC', 'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP',
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            '(Col)mcr-1.1', 'MCR',
        }

        self.critical_virulence_genes = {
            'stx1', 'stx2', 'stx1A', 'stx1B', 'stx2A', 'stx2B', 'stxA', 'stx2d1A',
            'cnf1', 'hlyA', 'hlyB', 'hlyC', 'hlyD', 'eae',
        }

        self.metadata = {
            "tool_name": "EcoliTyper ABRicate",
            "version": "1.4.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }

        self.science_quotes = [
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.",
             "author": "Albert Einstein"},
            {"text": "Nothing in life is to be feared, it is only to be understood.",
             "author": "Marie Curie"},
            {"text": "The microscope opens a new world to the investigator.",
             "author": "Robert Koch"},
            {"text": "In science, the credit goes to the man who convinces the world, not to the man to whom the idea first occurs.",
             "author": "Francis Darwin"},
            {"text": "The good thing about science is that it's true whether or not you believe it.",
             "author": "Neil deGrasse Tyson"},
            {"text": "Science knows no country, because knowledge belongs to humanity.",
             "author": "Louis Pasteur"},
            {"text": "Somewhere, something incredible is waiting to be known.",
             "author": "Carl Sagan"},
            {"text": "Research is what I'm doing when I don't know what I'm doing.",
             "author": "Wernher von Braun"},
            {"text": "The science of today is the technology of tomorrow.",
             "author": "Edward Teller"},
            {"text": "EcoliTyper transforms raw E. coli sequences into actionable genomic surveillance.",
             "author": "Brown Beckley"},
        ]

    def _setup_logging(self):
        """Configure and return a logger instance."""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)

    def _get_available_ram(self) -> float:
        """Return available RAM in gigabytes."""
        try:
            return psutil.virtual_memory().available / (1024 ** 3)
        except Exception as e:
            self.logger.warning("Could not detect RAM: %s", e)
            return 8

    def _calculate_optimal_cpus(self, user_cpus: int = None) -> int:
        """Return the CPU core count to use for this run."""
        if user_cpus is not None:
            self._log_resource_info(user_cpus)
            return user_cpus
        try:
            total_physical_cores = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            if total_physical_cores <= 4:
                optimal_cpus = total_physical_cores
            elif total_physical_cores <= 8:
                optimal_cpus = total_physical_cores - 1
            elif total_physical_cores <= 16:
                optimal_cpus = max(8, total_physical_cores - 1)
            elif total_physical_cores <= 32:
                optimal_cpus = max(16, total_physical_cores - 2)
            else:
                optimal_cpus = min(32, int(total_physical_cores * 0.95))
            optimal_cpus = max(1, min(optimal_cpus, total_physical_cores))
            self._log_resource_info(optimal_cpus, total_physical_cores)
            return optimal_cpus
        except Exception as e:
            self.logger.warning("Could not detect CPU cores, using maximum available: %s", e)
            return os.cpu_count() or 4

    def _log_resource_info(self, cpus: int, total_cores: int = None):
        """Log resource allocation information."""
        self.logger.info("Available RAM: %.1f GB", self.available_ram)
        if total_cores:
            self.logger.info("System CPU cores: %d", total_cores)
            utilization = (cpus / total_cores) * 100
            self.logger.info("Using CPU cores: %d (%.1f%% of available cores)", cpus, utilization)
        else:
            self.logger.info("Using user-specified CPU cores: %d", cpus)
        if cpus == 1:
            self.logger.info("Performance: single-core")
        elif cpus <= 4:
            self.logger.info("Performance: multi-core (small system)")
        elif cpus <= 8:
            self.logger.info("Performance: high-speed mode")
        else:
            self.logger.info("Performance: maximum-speed mode")

    def check_abricate_installed(self) -> bool:
        """Check that ABRicate is installed and at least version 1.2.0."""
        try:
            result = subprocess.run(
                ['abricate', '--version'],
                capture_output=True, text=True, check=True
            )
            version_line = result.stdout.strip()
            self.logger.info("ABRicate version: %s", version_line)
            version_match = re.search(r'(\d+\.\d+\.\d+)', version_line)
            if version_match:
                version_str = version_match.group(1)
                if version_str >= "1.2.0":
                    self.logger.info("ABRicate version meets requirement (>=1.2.0)")
                    return True
                self.logger.error("ABRicate version too old: %s. Required >=1.2.0", version_str)
                return False
            self.logger.info("ABRicate installed (version check skipped)")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error(
                "ABRicate not found. Install with: conda install -c bioconda abricate"
            )
            return False

    def setup_abricate_databases(self):
        """Ensure each required database is available, setting up missing ones."""
        self.logger.info("Setting up ABRicate databases for E. coli analysis...")
        available_dbs = []
        missing_dbs = []
        try:
            check_result = subprocess.run(
                ['abricate', '--list'],
                capture_output=True, text=True, check=True
            )
            for db in self.required_databases:
                if db in check_result.stdout:
                    self.logger.info("Database available: %s", db)
                    available_dbs.append(db)
                else:
                    self.logger.warning("Database not available: %s", db)
                    missing_dbs.append(db)
            for db in missing_dbs:
                self.logger.info("Attempting to setup database: %s", db)
                try:
                    subprocess.run(
                        ['abricate', '--setupdb', '--db', db],
                        capture_output=True, text=True, check=True
                    )
                    self.logger.info("Database setup completed: %s", db)
                    available_dbs.append(db)
                except subprocess.CalledProcessError as e:
                    self.logger.error("Failed to setup database %s: %s", db, e.stderr)
            self.required_databases = available_dbs
            self.logger.info("Using databases: %s", ", ".join(self.required_databases))
        except subprocess.CalledProcessError as e:
            self.logger.error("Error checking ABRicate databases: %s", e.stderr)
        except Exception as e:
            self.logger.error("Unexpected error setting up databases: %s", e)

    def run_abricate_single_db(self, genome_file: str, database: str,
                               temp_dir: str) -> Dict[str, Any]:
        """Run ABRicate for one genome against one database and return the raw content."""
        genome_name = Path(genome_file).stem
        raw_file = os.path.join(temp_dir, f"abricate_{database}.txt")
        cmd = [
            'abricate',
            genome_file,
            '--db', database,
            '--minid', str(self.min_identity),
            '--mincov', str(self.min_coverage),
        ]
        self.logger.info(
            "Running ABRicate: %s --db %s (minid=%d%%, mincov=%d%%)",
            genome_name, database, self.min_identity, self.min_coverage
        )
        try:
            with open(raw_file, 'w') as outfile:
                subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE,
                               text=True, check=True)
            hits = self._parse_abricate_output(raw_file)
            for hit in hits:
                if not hit.get('database'):
                    hit['database'] = database
            with open(raw_file, 'r') as f:
                raw_content = f.read()
            return {
                'database': database,
                'genome': genome_name,
                'hits': hits,
                'hit_count': len(hits),
                'raw_content': raw_content,
                'status': 'success',
            }
        except subprocess.CalledProcessError as e:
            self.logger.error("ABRicate failed for %s on %s: %s",
                              database, genome_name, e.stderr)
            return {
                'database': database,
                'genome': genome_name,
                'hits': [],
                'hit_count': 0,
                'raw_content': '',
                'status': 'failed',
            }

    def _parse_abricate_output(self, abricate_file: str) -> List[Dict]:
        """Parse an ABRicate tabular output file into a list of hit dictionaries."""
        hits = []
        try:
            with open(abricate_file, 'r') as f:
                lines = f.readlines()
            if not lines:
                return hits
            headers = []
            data_lines = []
            for line in lines:
                if line.startswith('#FILE') and not headers:
                    headers = line.strip().replace('#', '').split('\t')
                elif line.strip() and not line.startswith('#'):
                    data_lines.append(line.strip())
            if not headers:
                self.logger.warning("No headers found in %s", abricate_file)
                return hits
            expected_columns = len(headers)
            for line_num, line in enumerate(data_lines, 1):
                parts = line.split('\t')
                if len(parts) > expected_columns:
                    combined_parts = parts[:expected_columns - 1]
                    combined_parts.append('\t'.join(parts[expected_columns - 1:]))
                    parts = combined_parts
                elif len(parts) < expected_columns:
                    parts.extend([''] * (expected_columns - len(parts)))
                if len(parts) == expected_columns:
                    hit = {}
                    for i, header in enumerate(headers):
                        hit[header] = parts[i] if i < len(parts) else ''
                    processed_hit = {
                        'file': hit.get('FILE', ''),
                        'sequence': hit.get('SEQUENCE', ''),
                        'start': hit.get('START', ''),
                        'end': hit.get('END', ''),
                        'strand': hit.get('STRAND', ''),
                        'gene': hit.get('GENE', ''),
                        'coverage': hit.get('COVERAGE', ''),
                        'coverage_map': hit.get('COVERAGE_MAP', ''),
                        'gaps': hit.get('GAPS', ''),
                        'coverage_percent': hit.get('%COVERAGE', ''),
                        'identity_percent': hit.get('%IDENTITY', ''),
                        'database': hit.get('DATABASE', ''),
                        'accession': hit.get('ACCESSION', ''),
                        'product': hit.get('PRODUCT', ''),
                        'resistance': hit.get('RESISTANCE', ''),
                    }
                    hits.append(processed_hit)
                else:
                    self.logger.warning(
                        "Line %d has %d parts, expected %d: %s",
                        line_num, len(parts), expected_columns, line[:100] + "..."
                    )
        except Exception as e:
            self.logger.error("Error parsing %s: %s", abricate_file, e)
        self.logger.info("Parsed %d hits from %s", len(hits), abricate_file)
        return hits

    def _write_merged_raw_tsv(self, results: Dict[str, Dict], output_dir: str):
        """Merge all per-database ABRicate outputs into a single abricate_raw.tsv."""
        output_file = os.path.join(output_dir, "abricate_raw.tsv")
        header_written = False
        data_lines = 0

        with open(output_file, 'w') as f:
            for db, result in results.items():
                raw = result.get('raw_content', '') or ''
                if not raw:
                    continue
                for line in raw.splitlines():
                    if not line.strip():
                        continue
                    if line.startswith('#'):
                        if not header_written:
                            f.write(line + '\n')
                            header_written = True
                    else:
                        f.write(line + '\n')
                        data_lines += 1

        self.logger.info("Merged raw TSV: %s (%d hits across %d databases)",
                         output_file, data_lines, len(results))

    def _get_ascii_art(self) -> str:
        """Return the ECOLITYPER ASCII art banner."""
        return """
███████╗ ██████╗ ██████╗ ██╗     ██╗████████╗██╗   ██╗██████╗ ███████╗██████╗ 
██╔════╝██╔════╝██╔═══██╗██║     ██║╚══██╔══╝╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
█████╗  ██║     ██║   ██║██║     ██║   ██║    ╚████╔╝ ██████╔╝█████╗  ██████╔╝
██╔══╝  ██║     ██║   ██║██║     ██║   ██║     ╚██╔╝  ██╔═══╝ ██╔══╝  ██╔══██╗
███████╗╚██████╗╚██████╔╝███████╗██║   ██║      ██║   ██║     ███████╗██║  ██║
╚══════╝ ╚═════╝ ╚═════╝ ╚══════╝╚═╝   ╚═╝      ╚═╝   ╚═╝     ╚══════╝╚═╝  ╚═╝
        """

    def analyze_ecoli_resistance(self, all_hits: List[Dict]) -> Dict[str, Any]:
        """Aggregate resistance and virulence hits into a per-genome risk profile."""
        analysis = {
            'esbl_status': 'negative',
            'carbapenemase_status': 'negative',
            'colistin_resistance': 'negative',
            'critical_resistance_genes': [],
            'high_risk_resistance_genes': [],
            'critical_virulence_genes': [],
            'high_risk_virulence_genes': [],
            'moderate_risk_genes': [],
            'resistance_classes': {},
            'total_critical_resistance': 0,
            'total_high_risk_resistance': 0,
            'total_critical_virulence': 0,
            'total_high_risk_virulence': 0,
            'total_hits': len(all_hits),
        }

        for hit in all_hits:
            gene = hit['gene']
            gene_base = gene.split('-')[0] if '-' in gene else gene

            if any(crit_gene in gene_base for crit_gene in self.critical_resistance_genes):
                if any(esbl in gene_base for esbl in ['blaCTX-M', 'blaSHV', 'blaTEM']):
                    analysis['esbl_status'] = 'positive'
                    risk_level = 'ESBL'
                elif any(carba in gene_base for carba in ['blaKPC', 'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP']):
                    analysis['carbapenemase_status'] = 'positive'
                    risk_level = 'CARBAPENEMASE'
                elif any(mcr in gene for mcr in ['mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5',
                                                  'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10', 'MCR']):
                    analysis['colistin_resistance'] = 'positive'
                    risk_level = 'COLISTIN-RES'
                else:
                    risk_level = 'CRITICAL'
                analysis['critical_resistance_genes'].append({
                    'gene': gene, 'product': hit['product'], 'database': hit['database'],
                    'coverage': hit['coverage_percent'], 'identity': hit['identity_percent'],
                    'risk_level': risk_level,
                })
            elif any(hr_gene in gene_base for hr_gene in self.high_risk_genes):
                analysis['high_risk_resistance_genes'].append({
                    'gene': gene, 'product': hit['product'], 'database': hit['database'],
                    'coverage': hit['coverage_percent'], 'identity': hit['identity_percent'],
                    'risk_level': 'HIGH-RISK',
                })

            if any(crit_vf in gene_base for crit_vf in self.critical_virulence_genes):
                analysis['critical_virulence_genes'].append({
                    'gene': gene, 'product': hit['product'], 'database': hit['database'],
                    'coverage': hit['coverage_percent'], 'identity': hit['identity_percent'],
                    'risk_level': 'CRITICAL-VIRULENCE',
                })
            elif any(vf_gene in gene_base for vf_gene in self.virulence_genes):
                is_high_risk_virulence = any(hr_vf in gene_base for hr_vf in [
                    'stx', 'cnf', 'hly', 'ast', 'east', 'elt', 'est', 'eae', 'tir', 'paa',
                    'afa', 'pap', 'sfa', 'fim', 'fae', 'fan', 'fyu', 'irp', 'ybt',
                    'iuc', 'iut', 'chu', 'iro', 'shu', 'esc', 'esp', 'nle', 'map',
                    'kps', 'ibe', 'vat', 'sen', 'cif', 'efa', 'stc', 'pic',
                ])
                if is_high_risk_virulence:
                    analysis['high_risk_virulence_genes'].append({
                        'gene': gene, 'product': hit['product'], 'database': hit['database'],
                        'coverage': hit['coverage_percent'], 'identity': hit['identity_percent'],
                        'risk_level': 'HIGH-VIRULENCE',
                    })
                else:
                    analysis['moderate_risk_genes'].append({
                        'gene': gene, 'product': hit['product'], 'database': hit['database'],
                        'coverage': hit['coverage_percent'], 'identity': hit['identity_percent'],
                        'risk_level': 'MODERATE',
                    })

            resistance_class = self._classify_resistance(hit['product'])
            if resistance_class:
                if resistance_class not in analysis['resistance_classes']:
                    analysis['resistance_classes'][resistance_class] = []
                if gene not in [g['gene'] for g in analysis['resistance_classes'][resistance_class]]:
                    analysis['resistance_classes'][resistance_class].append({
                        'gene': gene,
                        'product': hit['product'],
                    })

        analysis['total_critical_resistance'] = len(analysis['critical_resistance_genes'])
        analysis['total_high_risk_resistance'] = len(analysis['high_risk_resistance_genes'])
        analysis['total_critical_virulence'] = len(analysis['critical_virulence_genes'])
        analysis['total_high_risk_virulence'] = len(analysis['high_risk_virulence_genes'])
        return analysis

    def _classify_resistance(self, product: str) -> str:
        """Assign a resistance class label based on the product description."""
        product_lower = product.lower()
        if any(term in product_lower for term in ['beta-lactam', 'esbl', 'carbapenem', 'penicillin']):
            return 'Beta-lactam resistance'
        elif any(term in product_lower for term in ['aminoglycoside', 'aac', 'aad', 'aph', 'str']):
            return 'Aminoglycoside resistance'
        elif any(term in product_lower for term in ['fluoroquinolone', 'quinolone', 'qnr']):
            return 'Fluoroquinolone resistance'
        elif any(term in product_lower for term in ['macrolide', 'mph']):
            return 'Macrolide resistance'
        elif any(term in product_lower for term in ['tetracycline', 'tet']):
            return 'Tetracycline resistance'
        elif any(term in product_lower for term in ['sulfonamide', 'sul']):
            return 'Sulfonamide resistance'
        elif any(term in product_lower for term in ['trimethoprim', 'dfr']):
            return 'Trimethoprim resistance'
        elif any(term in product_lower for term in ['chloramphenicol', 'cat', 'cm']):
            return 'Chloramphenicol resistance'
        elif any(term in product_lower for term in ['colistin', 'polymyxin', 'mcr']):
            return 'Polymyxin resistance'
        elif any(term in product_lower for term in ['efflux', 'multidrug']):
            return 'Efflux pumps'
        return 'Other resistance'

    def _get_shared_css(self) -> str:
        """Return the shared CSS block used by both comprehensive and batch HTML reports."""
        return """
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff; padding: 20px; min-height: 100vh;
        }
        .container { max-width: 1500px; margin: 0 auto; }
        .header { text-align: center; margin-bottom: 30px; }
        .ascii-container {
            background: rgba(0, 0, 0, 0.7); padding: 20px; border-radius: 15px;
            margin-bottom: 20px; box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3);
        }
        .ascii-art {
            font-family: 'Courier New', monospace; font-size: 10px; line-height: 1.1;
            white-space: pre; color: #00ff00;
            text-shadow: 0 0 10px rgba(0, 255, 0, 0.5); overflow-x: auto;
        }
        .quote-container {
            background: rgba(255, 255, 255, 0.1); backdrop-filter: blur(10px);
            padding: 20px; border-radius: 10px; margin-bottom: 30px; text-align: center;
            min-height: 100px; display: flex; flex-direction: column;
            justify-content: center; box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(255, 255, 255, 0.2); transition: opacity 0.5s ease-in-out;
        }
        .quote-text { font-size: 18px; font-style: italic; margin-bottom: 10px; color: #ffffff; }
        .quote-author { font-size: 14px; color: #fbbf24; font-weight: bold; }
        .report-section {
            background: rgba(255, 255, 255, 0.95); color: #1f2937;
            padding: 25px; border-radius: 10px; margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }
        .report-section h2 {
            color: #1e3a8a; border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px; margin-bottom: 20px; font-size: 24px;
        }
        .metrics-grid {
            display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 20px; margin-top: 15px;
        }
        .metric-card {
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white; padding: 20px; border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }
        .metric-label { font-size: 14px; opacity: 0.9; margin-bottom: 5px; }
        .metric-value { font-size: 24px; font-weight: bold; }
        .summary-table {
            width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 14px;
        }
        .summary-table th {
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white; padding: 12px; text-align: left; font-weight: bold;
            position: sticky; top: 0;
        }
        .summary-table td { padding: 10px; border-bottom: 1px solid #e5e7eb; }
        .summary-table tr:nth-child(even) { background-color: #f8fafc; }
        .summary-table tr:hover { background-color: #e0f2fe; }
        .critical { background-color: #fee2e2; font-weight: bold; }
        .high-risk { background-color: #fef3c7; }
        .present { background-color: #d1fae5; }
        .product-cell {
            white-space: normal !important; word-wrap: break-word;
            max-width: 500px; min-width: 200px;
        }
        .table-responsive { width: 100%; overflow-x: auto; margin: 20px 0; }
        .risk-badge {
            display: inline-block; background: #dc2626; color: white;
            padding: 5px 10px; border-radius: 15px; margin: 2px; font-size: 0.9em;
            font-weight: bold;
        }
        .warning-badge {
            display: inline-block; background: #f59e0b; color: black;
            padding: 5px 10px; border-radius: 15px; margin: 2px; font-size: 0.9em;
            font-weight: bold;
        }
        .safe-badge {
            display: inline-block; background: #16a34a; color: white;
            padding: 5px 10px; border-radius: 15px; margin: 2px; font-size: 0.9em;
            font-weight: bold;
        }
        .stats-grid {
            display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px; margin-bottom: 20px;
        }
        .stat-card {
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white; padding: 15px; border-radius: 8px;
            text-align: center; box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }
        .stat-value { font-size: 24px; font-weight: bold; margin-bottom: 5px; }
        .stat-label { font-size: 12px; opacity: 0.9; }
        .percentage-badge {
            display: inline-block; background: #e3f2fd; color: #1976d2;
            padding: 2px 8px; border-radius: 12px; font-size: 0.85em;
            margin-left: 5px; border: 1px solid #bbdefb;
        }
        .frequency-cell { font-weight: 600; color: #2c3e50; }
        .footer {
            text-align: center; margin-top: 30px; padding: 20px;
            background: rgba(0, 0, 0, 0.3); border-radius: 10px; font-size: 14px;
        }
        .timestamp { color: #fbbf24; font-weight: bold; }
        .authorship {
            margin-top: 15px; padding: 15px;
            background: rgba(255, 255, 255, 0.1); border-radius: 8px; font-size: 12px;
        }
        .authorship a { color: #fbbf24; text-decoration: none; }
        .authorship a:hover { text-decoration: underline; }
        @media (max-width: 768px) {
            .ascii-art { font-size: 6px; }
            .metrics-grid { grid-template-columns: 1fr; }
            .summary-table { font-size: 12px; }
            .summary-table th, .summary-table td { padding: 6px; }
        }
        """

    def _get_quote_js(self) -> str:
        """Return the rotating-quote JavaScript block."""
        return f"""
        <script>
            const quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            const quoteContainer = document.getElementById('quoteContainer');
            const quoteText = document.getElementById('quoteText');
            const quoteAuthor = document.getElementById('quoteAuthor');
            function displayQuote() {{
                quoteContainer.style.opacity = '0';
                setTimeout(() => {{
                    const q = quotes[currentQuote];
                    quoteText.textContent = '"' + q.text + '"';
                    quoteAuthor.textContent = '\\u2014 ' + q.author;
                    quoteContainer.style.opacity = '1';
                    currentQuote = (currentQuote + 1) % quotes.length;
                }}, 500);
            }}
            document.addEventListener('DOMContentLoaded', function() {{
                const q = quotes[currentQuote];
                quoteText.textContent = '"' + q.text + '"';
                quoteAuthor.textContent = '\\u2014 ' + q.author;
                currentQuote = (currentQuote + 1) % quotes.length;
                setInterval(displayQuote, 10000);
            }});
        </script>
        """

    def _get_footer_html(self, title: str) -> str:
        """Return the shared footer block for HTML reports."""
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        return f"""
        <div class="footer">
            <p><strong>{title}</strong></p>
            <p class="timestamp">Generated: {current_time}</p>
            <div class="authorship">
                <p><strong>Technical Support &amp; Inquiries:</strong></p>
                <p>Author: Brown Beckley | GitHub: <a href="https://github.com/bbeckley-hub" target="_blank">bbeckley-hub</a></p>
                <p>Email: <a href="mailto:brownbeckley94@gmail.com">brownbeckley94@gmail.com</a></p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
                <p style="margin-top:10px;font-size:11px;opacity:0.8;">
                    Powered by ABRicate (Seemann) | Thresholds: minID={self.min_identity}%, minCov={self.min_coverage}%
                </p>
            </div>
        </div>
        """

    def _render_hit_table(self, title: str, genes: List[Dict], row_class: str,
                          badge_class: str, show_risk_level: bool) -> str:
        """Render an HTML table for a category of gene hits."""
        if not genes:
            return ""
        risk_col = "<th>Risk Level</th>" if show_risk_level else ""
        html = f"""
        <div class="report-section">
            <h2>{title}</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Product</th>
                            <th>Database</th>
                            <th>Coverage</th>
                            <th>Identity</th>
                            {risk_col}
                        </tr>
                    </thead>
                    <tbody>
        """
        for gene_info in genes:
            risk_cell = ""
            if show_risk_level:
                badge = badge_class if badge_class else "warning-badge"
                risk_cell = f'<td><span class="{badge}">{gene_info["risk_level"]}</span></td>'
            html += f"""
                        <tr class="{row_class}">
                            <td><strong>{gene_info['gene']}</strong></td>
                            <td class="product-cell">{gene_info['product']}</td>
                            <td>{gene_info['database']}</td>
                            <td>{gene_info['coverage']}%</td>
                            <td>{gene_info['identity']}%</td>
                            {risk_cell}
                        </tr>
            """
        html += """
                    </tbody>
                </table>
            </div>
        </div>
        """
        return html

    def create_comprehensive_html_report(self, genome_name: str, results: Dict,
                                         output_dir: str):
        """Create the single per-genome comprehensive HTML report."""
        all_hits = []
        for db_result in results.values():
            all_hits.extend(db_result['hits'])

        analysis = self.analyze_ecoli_resistance(all_hits)
        random_quote = random.choice(self.science_quotes)

        alerts_html = ""
        critical_alerts = []
        if analysis['esbl_status'] == 'positive':
            critical_alerts.append("ESBL DETECTED")
        if analysis['carbapenemase_status'] == 'positive':
            critical_alerts.append("CARBAPENEMASE DETECTED")
        if analysis['colistin_resistance'] == 'positive':
            critical_alerts.append("COLISTIN RESISTANCE DETECTED")
        if critical_alerts:
            badges = "".join(f'<span class="risk-badge">{a}</span>' for a in critical_alerts)
            alerts_html = f"""
        <div class="report-section" style="border-left: 4px solid #dc2626;">
            <h2 style="color: #dc2626;">CRITICAL RESISTANCE ALERTS</h2>
            <div style="margin: 10px 0;">{badges}</div>
        </div>
        """

        virulence_alerts_html = ""
        if analysis['critical_virulence_genes']:
            vf_genes = sorted(set(g['gene'] for g in analysis['critical_virulence_genes']))
            vf_badges = "".join(f'<span class="warning-badge">{g}</span>' for g in vf_genes)
            virulence_alerts_html = f"""
        <div class="report-section" style="border-left: 4px solid #f59e0b;">
            <h2 style="color: #f59e0b;">HIGH-VIRULENCE ALERT</h2>
            <p>Critical virulence determinants detected: {", ".join(vf_genes)}</p>
            <div style="margin: 10px 0;">{vf_badges}</div>
        </div>
        """

        resistance_classes_html = ""
        if analysis['resistance_classes']:
            blocks = ""
            for class_name, genes in analysis['resistance_classes'].items():
                gene_list = ", ".join(sorted(g['gene'] for g in genes))
                blocks += f"""
                <div style="margin: 10px 0; padding: 10px; background: #f8f9fa; border-radius: 8px;">
                    <strong style="color: #667eea;">{class_name}</strong> ({len(genes)} genes)
                    <br><span style="color: #666; font-size: 0.9em;">{gene_list}</span>
                </div>
                """
            resistance_classes_html = f"""
        <div class="report-section">
            <h2>Resistance Classes Detected</h2>
            {blocks}
        </div>
        """

        tables_html = ""
        tables_html += self._render_hit_table(
            "CRITICAL Resistance Genes", analysis['critical_resistance_genes'],
            "critical", "risk-badge", True)
        tables_html += self._render_hit_table(
            "High-Risk Resistance Genes", analysis['high_risk_resistance_genes'],
            "high-risk", "warning-badge", True)
        tables_html += self._render_hit_table(
            "CRITICAL Virulence Genes", analysis['critical_virulence_genes'],
            "critical", "risk-badge", True)
        tables_html += self._render_hit_table(
            "High-Risk Virulence Genes", analysis['high_risk_virulence_genes'],
            "high-risk", "warning-badge", True)
        tables_html += self._render_hit_table(
            "Moderate Risk Genes", analysis['moderate_risk_genes'],
            "present", "safe-badge", True)

        db_summary_rows = ""
        for db, result in results.items():
            status_icon = "OK" if result['status'] == 'success' else "FAIL"
            db_summary_rows += f"""
                        <tr class="present">
                            <td>{db}</td>
                            <td>{result['hit_count']}</td>
                            <td>{status_icon} {result['status']}</td>
                        </tr>
            """

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliTyper ABRicate - Comprehensive Report - {genome_name}</title>
    <style>{self._get_shared_css()}</style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self._get_ascii_art()}</div>
            </div>
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">&mdash; {random_quote['author']}</div>
            </div>
        </div>

        <div class="report-section">
            <h2>E. coli AMR Summary</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Total Genes</div>
                    <div class="metric-value">{analysis['total_hits']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Critical Resistance</div>
                    <div class="metric-value">{analysis['total_critical_resistance']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Critical Virulence</div>
                    <div class="metric-value">{analysis['total_critical_virulence']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Analysis Date</div>
                    <div class="metric-value">{datetime.now().strftime('%Y-%m-%d')}</div>
                </div>
            </div>
            <p style="margin-top:15px;"><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>Thresholds:</strong> minID={self.min_identity}%, minCov={self.min_coverage}%</p>
        </div>

        {alerts_html}
        {virulence_alerts_html}
        {resistance_classes_html}
        {tables_html}

        <div class="report-section">
            <h2>Database Results Summary</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Database</th>
                            <th>Hits</th>
                            <th>Status</th>
                        </tr>
                    </thead>
                    <tbody>
                        {db_summary_rows}
                    </tbody>
                </table>
            </div>
        </div>

        {self._get_footer_html("EcoliTyper - Comprehensive ABRicate Report")}
    </div>
    {self._get_quote_js()}
</body>
</html>"""

        html_file = os.path.join(output_dir, f"{genome_name}_comprehensive_abricate_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("Comprehensive HTML report generated: %s", html_file)

    def create_database_summaries(self, all_results: Dict[str, Any], output_base: str):
        """Write per-database TSV and HTML batch summaries."""
        self.logger.info("Creating E. coli database summary files and HTML reports...")

        db_results = {}
        for genome_name, genome_result in all_results.items():
            for db, db_result in genome_result['results'].items():
                if db not in db_results:
                    db_results[db] = []
                for hit in db_result['hits']:
                    hit_with_genome = hit.copy()
                    hit_with_genome['genome'] = genome_name
                    db_results[db].append(hit_with_genome)

        for db, hits in db_results.items():
            if hits:
                summary_file = os.path.join(output_base, f"ecoli_{db}_abricate_summary.tsv")
                headers = list(hits[0].keys())
                with open(summary_file, 'w') as f:
                    f.write('\t'.join(headers) + '\n')
                    for hit in hits:
                        row = [str(hit.get(header, '')) for header in headers]
                        f.write('\t'.join(row) + '\n')
                self.logger.info("Created %s summary: %s (%d hits)",
                                 db, summary_file, len(hits))
                self._create_database_summary_html(db, hits, output_base)
            else:
                self.logger.info("No hits for database %s, skipping summary", db)

    def _create_database_summary_html(self, database: str, hits: List[Dict], output_base: str):
        """Write the batch-level per-database HTML summary across all genomes."""
        random_quote = random.choice(self.science_quotes)
        unique_genomes = list(set(hit['genome'] for hit in hits))
        total_genomes = len(unique_genomes)

        genes_per_genome = {}
        for hit in hits:
            genome = hit['genome']
            if genome not in genes_per_genome:
                genes_per_genome[genome] = set()
            genes_per_genome[genome].add(hit['gene'])

        genes_by_genome_rows = ""
        for genome in sorted(unique_genomes):
            genes = genes_per_genome.get(genome, set())
            gene_list = ", ".join(sorted(genes))
            genes_by_genome_rows += f"""
                        <tr class="present">
                            <td><strong>{genome}</strong></td>
                            <td>{len(genes)}</td>
                            <td>{gene_list}</td>
                        </tr>
            """

        gene_frequency = {}
        for hit in hits:
            gene = hit['gene']
            if gene not in gene_frequency:
                gene_frequency[gene] = set()
            gene_frequency[gene].add(hit['genome'])

        frequency_rows = ""
        for gene, genomes in sorted(gene_frequency.items(),
                                    key=lambda x: len(x[1]), reverse=True):
            genome_list = ", ".join(sorted(genomes))
            count = len(genomes)
            percentage = (count / total_genomes * 100) if total_genomes > 0 else 0
            frequency_rows += f"""
                        <tr>
                            <td><strong>{gene}</strong></td>
                            <td class="frequency-cell">
                                {count}
                                <span class="percentage-badge">({percentage:.1f}%)</span>
                            </td>
                            <td>{genome_list}</td>
                        </tr>
            """

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliTyper ABRicate - Batch Summary - {database.upper()}</title>
    <style>{self._get_shared_css()}</style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self._get_ascii_art()}</div>
            </div>
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">&mdash; {random_quote['author']}</div>
            </div>
        </div>

        <div class="report-section">
            <h2>Batch ABRicate Summary - {database.upper()} Database</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{len(hits)}</div>
                    <div class="stat-label">TOTAL HITS</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{total_genomes}</div>
                    <div class="stat-label">GENOMES ANALYZED</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{len(set(hit['gene'] for hit in hits))}</div>
                    <div class="stat-label">UNIQUE GENES</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{datetime.now().strftime('%Y-%m-%d')}</div>
                    <div class="stat-label">ANALYSIS DATE</div>
                </div>
            </div>
        </div>

        <div class="report-section">
            <h2>Genes by Genome</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Genome</th>
                            <th>Gene Count</th>
                            <th>Genes Detected</th>
                        </tr>
                    </thead>
                    <tbody>
                        {genes_by_genome_rows}
                    </tbody>
                </table>
            </div>
        </div>

        <div class="report-section">
            <h2>Gene Frequency</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Frequency</th>
                            <th>Genomes</th>
                        </tr>
                    </thead>
                    <tbody>
                        {frequency_rows}
                    </tbody>
                </table>
            </div>
        </div>

        {self._get_footer_html(f"EcoliTyper - {database.upper()} Batch Summary")}
    </div>
    {self._get_quote_js()}
</body>
</html>"""

        html_file = os.path.join(output_base, f"ecoli_{database}_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("Database summary HTML report: %s", html_file)

    def create_master_json_summary(self, all_results: Dict[str, Any], output_base: str):
        """Write the single master JSON summarising every database across all genomes."""
        self.logger.info("Creating master JSON summary for all databases and genomes...")

        master_report = {
            "metadata": {
                **self.metadata,
                "analysis_type": "E. coli Comprehensive Antimicrobial Resistance & Virulence Analysis",
                "databases_analyzed": self.required_databases,
                "total_genomes": len(all_results),
                "analysis_parameters": {
                    "minimum_identity": self.min_identity,
                    "minimum_coverage": self.min_coverage,
                },
            },
            "overall_summary": {
                "total_genomes": len(all_results),
                "total_hits": 0,
                "critical_resistance_genes": 0,
                "critical_virulence_genes": 0,
                "esbl_positive_genomes": 0,
                "carbapenemase_positive_genomes": 0,
                "colistin_resistant_genomes": 0,
            },
            "database_summaries": {},
            "genome_summaries": {},
            "risk_assessment": {
                "critical_resistance_genes_found": [],
                "critical_virulence_genes_found": [],
                "high_risk_genes_found": [],
                "resistance_classes_found": defaultdict(list),
            },
            "detailed_results": {},
        }

        for genome_name, genome_result in all_results.items():
            all_genome_hits = []
            for db_result in genome_result['results'].values():
                all_genome_hits.extend(db_result['hits'])

            analysis = self.analyze_ecoli_resistance(all_genome_hits)

            master_report["overall_summary"]["total_hits"] += analysis['total_hits']
            master_report["overall_summary"]["critical_resistance_genes"] += analysis['total_critical_resistance']
            master_report["overall_summary"]["critical_virulence_genes"] += analysis['total_critical_virulence']
            if analysis['esbl_status'] == 'positive':
                master_report["overall_summary"]["esbl_positive_genomes"] += 1
            if analysis['carbapenemase_status'] == 'positive':
                master_report["overall_summary"]["carbapenemase_positive_genomes"] += 1
            if analysis['colistin_resistance'] == 'positive':
                master_report["overall_summary"]["colistin_resistant_genomes"] += 1

            master_report["genome_summaries"][genome_name] = {
                "total_hits": analysis['total_hits'],
                "critical_resistance": analysis['total_critical_resistance'],
                "critical_virulence": analysis['total_critical_virulence'],
                "esbl_status": analysis['esbl_status'],
                "carbapenemase_status": analysis['carbapenemase_status'],
                "colistin_resistance": analysis['colistin_resistance'],
                "database_hits": {db: result['hit_count']
                                  for db, result in genome_result['results'].items()},
            }

            for gene_info in analysis['critical_resistance_genes']:
                if gene_info['gene'] not in master_report["risk_assessment"]["critical_resistance_genes_found"]:
                    master_report["risk_assessment"]["critical_resistance_genes_found"].append(gene_info['gene'])
            for gene_info in analysis['critical_virulence_genes']:
                if gene_info['gene'] not in master_report["risk_assessment"]["critical_virulence_genes_found"]:
                    master_report["risk_assessment"]["critical_virulence_genes_found"].append(gene_info['gene'])
            for gene_info in analysis['high_risk_resistance_genes'] + analysis['high_risk_virulence_genes']:
                if gene_info['gene'] not in master_report["risk_assessment"]["high_risk_genes_found"]:
                    master_report["risk_assessment"]["high_risk_genes_found"].append(gene_info['gene'])
            for class_name, genes in analysis['resistance_classes'].items():
                for gene_info in genes:
                    if gene_info['gene'] not in [g['gene']
                                                 for g in master_report["risk_assessment"]["resistance_classes_found"][class_name]]:
                        master_report["risk_assessment"]["resistance_classes_found"][class_name].append(gene_info)

            master_report["detailed_results"][genome_name] = {
                "total_hits": genome_result['total_hits'],
                "database_results": genome_result['results'],
            }

        for db in self.required_databases:
            db_hits = []
            for genome_name, genome_result in all_results.items():
                if db in genome_result['results']:
                    for hit in genome_result['results'][db]['hits']:
                        hit_with_genome = hit.copy()
                        hit_with_genome['genome'] = genome_name
                        db_hits.append(hit_with_genome)
            if db_hits:
                unique_genes = len(set(h['gene'] for h in db_hits))
                unique_genomes = len(set(h['genome'] for h in db_hits))
                master_report["database_summaries"][db] = {
                    "total_hits": len(db_hits),
                    "unique_genes": unique_genes,
                    "unique_genomes": unique_genomes,
                    "hits_per_genome": round(len(db_hits) / unique_genomes, 2) if unique_genomes else 0,
                    "most_frequent_genes": dict(Counter([h['gene'] for h in db_hits]).most_common(10)),
                }
            else:
                master_report["database_summaries"][db] = {
                    "total_hits": 0,
                    "unique_genes": 0,
                    "unique_genomes": 0,
                    "hits_per_genome": 0,
                    "most_frequent_genes": {},
                }

        master_report["risk_assessment"]["critical_resistance_genes_found"].sort()
        master_report["risk_assessment"]["critical_virulence_genes_found"].sort()
        master_report["risk_assessment"]["high_risk_genes_found"].sort()
        master_report["risk_assessment"]["resistance_classes_found"] = dict(
            master_report["risk_assessment"]["resistance_classes_found"]
        )

        master_report["performance_metrics"] = {
            "cpu_cores_used": self.cpus,
            "available_ram_gb": round(self.available_ram, 2),
            "processing_mode": "MAXIMUM SPEED" if self.cpus > 8 else "STANDARD",
            "thresholds": {
                "minimum_identity": self.min_identity,
                "minimum_coverage": self.min_coverage,
            },
        }

        master_file = os.path.join(output_base, "ecoli_abricate_master_summary.json")
        with open(master_file, 'w') as f:
            json.dump(master_report, f, indent=4, default=str)
        self.logger.info("Created master JSON summary: %s", master_file)
        return master_file

    def process_single_genome(self, genome_file: str,
                              output_base: str = "ecoli_abricate_results") -> Dict[str, Any]:
        """Process one genome. Emits exactly 2 files: comprehensive HTML + merged raw TSV."""
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        self.logger.info("=== PROCESSING E. COLI GENOME: %s ===", genome_name)
        os.makedirs(results_dir, exist_ok=True)

        results = {}

        with tempfile.TemporaryDirectory() as temp_dir:
            for db in self.required_databases:
                result = self.run_abricate_single_db(genome_file, db, temp_dir)
                results[db] = result
                status_icon = "OK" if result['status'] == 'success' else "FAIL"
                self.logger.info("%s %s: %d hits", status_icon, db, result['hit_count'])

            self._write_merged_raw_tsv(results, results_dir)

        for db_result in results.values():
            db_result.pop('raw_content', None)

        self.create_comprehensive_html_report(genome_name, results, results_dir)

        return {
            'genome': genome_name,
            'results': results,
            'total_hits': sum(r['hit_count'] for r in results.values()),
        }

    def process_multiple_genomes(self, genome_pattern: str,
                                 output_base: str = "ecoli_abricate_results") -> Dict[str, Any]:
        """Process multiple genomes in parallel and write all reports."""
        print("\n" + "=" * 80)
        print(self._get_ascii_art())
        print("=" * 80)
        print("EcoliTyper ABRicate")
        print(f"   Thresholds: minID={self.min_identity}%, minCov={self.min_coverage}%")
        print("=" * 80)

        if not self.check_abricate_installed():
            raise RuntimeError("ABRicate not properly installed")

        self.setup_abricate_databases()

        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta",
                          f"{genome_pattern}.fa", f"{genome_pattern}.fna",
                          f"{genome_pattern}.faa"]
        genome_files = []
        for pattern in fasta_patterns:
            genome_files.extend(glob.glob(pattern))
        genome_files = list(set(genome_files))

        if not genome_files:
            raise FileNotFoundError(
                f"No FASTA files found matching pattern: {genome_pattern}"
            )

        self.logger.info("Found %d E. coli genomes: %s",
                         len(genome_files), [Path(f).name for f in genome_files])

        os.makedirs(output_base, exist_ok=True)
        all_results = {}

        if len(genome_files) > 1 and self.cpus > 1:
            self.logger.info("Using parallel processing with %d CPU cores", self.cpus)
            with ThreadPoolExecutor(max_workers=self.cpus) as executor:
                future_to_genome = {
                    executor.submit(self.process_single_genome, genome, output_base): genome
                    for genome in genome_files
                }
                for future in as_completed(future_to_genome):
                    genome = future_to_genome[future]
                    try:
                        result = future.result()
                        all_results[Path(genome).stem] = result
                        self.logger.info("Completed: %s (%d total hits)",
                                         result['genome'], result['total_hits'])
                    except Exception as e:
                        self.logger.error("Failed: %s - %s", genome, e)
        else:
            for genome in genome_files:
                try:
                    result = self.process_single_genome(genome, output_base)
                    all_results[Path(genome).stem] = result
                    self.logger.info("Completed: %s (%d total hits)",
                                     result['genome'], result['total_hits'])
                except Exception as e:
                    self.logger.error("Failed: %s - %s", genome, e)

        self.create_database_summaries(all_results, output_base)
        master_json_file = self.create_master_json_summary(all_results, output_base)

        self.logger.info("=== E. COLI ANALYSIS COMPLETE ===")
        self.logger.info("Processed %d genomes", len(all_results))
        self.logger.info("Results saved to: %s", output_base)
        self.logger.info("Master JSON summary: %s", master_json_file)

        return all_results


def main():
    """Command-line entry point for E. coli ABRicate analysis."""
    parser = argparse.ArgumentParser(
        description='EcoliTyper ABRicate Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python ecoli_abricate.py "*.fna"
  python ecoli_abricate.py "ECOLI_*.fasta"
  python ecoli_abricate.py "*.fa" --output my_ecoli_results
  python ecoli_abricate.py "*.fna" --cpus 4
  python ecoli_abricate.py "*.fna" --minid 90 --mincov 85

Resource management:
  CPU cores are auto-detected and scaled based on total physical cores.
  Use --cpus to force a specific count.

Supported FASTA extensions: .fasta, .fa, .fna, .faa

OUTPUTS:
  Per genome (2 files each):
    - {genome}_comprehensive_abricate_report.html
    - abricate_raw.tsv                    (all databases merged)

  Batch-level:
    - ecoli_{db}_abricate_summary.tsv     (per database)
    - ecoli_{db}_summary_report.html      (per database)
    - ecoli_abricate_master_summary.json  (master summary)
        """
    )
    parser.add_argument('pattern',
                        help='File pattern for E. coli genomes (e.g., "*.fasta")')
    parser.add_argument('--cpus', '-c', type=int, default=None,
                        help='Number of CPU cores to use (default: auto-detect)')
    parser.add_argument('--output', '-o', default='ecoli_abricate_results',
                        help='Output directory (default: ecoli_abricate_results)')
    parser.add_argument('--minid', type=int, default=80,
                        help='Minimum identity threshold for ABRicate (default: 80)')
    parser.add_argument('--mincov', type=int, default=80,
                        help='Minimum coverage threshold for ABRicate (default: 80)')
    args = parser.parse_args()

    if not 0 <= args.minid <= 100:
        parser.error("--minid must be between 0 and 100")
    if not 0 <= args.mincov <= 100:
        parser.error("--mincov must be between 0 and 100")

    executor = EcoliAbricateExecutor(
        cpus=args.cpus,
        min_identity=args.minid,
        min_coverage=args.mincov,
    )

    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)

        executor.logger.info("\n" + "=" * 50)
        executor.logger.info("EcoliTyper ABRicate FINAL SUMMARY")
        executor.logger.info("=" * 50)

        total_critical_resistance = 0
        total_critical_virulence = 0
        total_esbl = 0
        total_carba = 0
        total_colistin = 0

        for genome_name, result in results.items():
            all_genome_hits = []
            for db_result in result['results'].values():
                all_genome_hits.extend(db_result['hits'])
            analysis = executor.analyze_ecoli_resistance(all_genome_hits)
            executor.logger.info(
                "%s: %d total hits, %d critical resistance, %d critical virulence, "
                "ESBL=%s, Carba=%s, Colistin=%s",
                genome_name, result['total_hits'],
                analysis['total_critical_resistance'],
                analysis['total_critical_virulence'],
                analysis['esbl_status'],
                analysis['carbapenemase_status'],
                analysis['colistin_resistance'],
            )
            total_critical_resistance += analysis['total_critical_resistance']
            total_critical_virulence += analysis['total_critical_virulence']
            if analysis['esbl_status'] == 'positive':
                total_esbl += 1
            if analysis['carbapenemase_status'] == 'positive':
                total_carba += 1
            if analysis['colistin_resistance'] == 'positive':
                total_colistin += 1

        executor.logger.info("\n" + "=" * 50)
        executor.logger.info("DATABASE USAGE SUMMARY")
        executor.logger.info("=" * 50)
        executor.logger.info("Used databases: %s",
                             ", ".join(executor.required_databases))
        executor.logger.info("Thresholds: minID=%d%%, minCov=%d%%",
                             executor.min_identity, executor.min_coverage)

        executor.logger.info("\n" + "=" * 50)
        executor.logger.info("OUTPUT FILES SUMMARY")
        executor.logger.info("=" * 50)
        executor.logger.info("Per genome (2 files each):")
        executor.logger.info("  - {genome}_comprehensive_abricate_report.html")
        executor.logger.info("  - abricate_raw.tsv  (all databases merged)")
        executor.logger.info("Batch-level:")
        executor.logger.info("  - ecoli_{db}_abricate_summary.tsv  (per database)")
        executor.logger.info("  - ecoli_{db}_summary_report.html   (per database)")
        executor.logger.info("  - ecoli_abricate_master_summary.json  (master)")
        executor.logger.info("Location: %s", args.output)

        executor.logger.info("\n" + "=" * 50)
        executor.logger.info("PERFORMANCE SUMMARY")
        executor.logger.info("=" * 50)
        executor.logger.info("CPU cores utilized: %d cores", executor.cpus)
        executor.logger.info("Available RAM: %.1f GB", executor.available_ram)
        executor.logger.info("Total E. coli genomes processed: %d", len(results))
        executor.logger.info("Total critical resistance genes found: %d", total_critical_resistance)
        executor.logger.info("Total critical virulence genes found: %d", total_critical_virulence)
        executor.logger.info("ESBL-positive genomes: %d", total_esbl)
        executor.logger.info("Carbapenemase-positive genomes: %d", total_carba)
        executor.logger.info("Colistin-resistant genomes: %d", total_colistin)

        quote = random.choice(executor.science_quotes)
        executor.logger.info("\n%s - %s", quote['text'], quote['author'])

    except Exception as e:
        executor.logger.error("E. coli analysis failed: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()