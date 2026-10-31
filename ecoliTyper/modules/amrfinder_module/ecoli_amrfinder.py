#!/usr/bin/env python3
"""
EcoliTyper AMRfinderPlus - E. coli AMR Analysis with Dynamic Database
Comprehensive AMR analysis for E. coli with HTML reporting.

Per-genome outputs (2 files each):
  - {genome}_amrfinder_report.html   (AMR genes + point mutations combined)
  - {genome}_amrfinder_raw.tsv       (merged: Source = AMR / MUTATION)

Batch-level outputs:
  - ecoli_amrfinder_summary.tsv
  - ecoli_amrfinder_statistics_summary.tsv
  - ecoli_amrfinder_summary_report.html
  - ecoli_amrfinder_master_summary.json
  - mutation_summary.tsv
  - mutation_summary.html
  - mutation_master_summary.json

Author: Beckley Brown <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2025 / Updated 2026-07-22
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
from typing import List, Dict, Any, Optional
import argparse
import re
from datetime import datetime
import psutil
import math
import json
from collections import defaultdict


class EcoliAMRfinderPlus:
    """AMRfinderPlus executor for E. coli with merged per-sample reporting."""

    def __init__(self, cpus: int = None):
        self.logger = self._setup_logging()
        self.module_dir = os.path.dirname(os.path.abspath(__file__))
        self.bundled_amrfinder = os.path.join(self.module_dir, "bin", "amrfinder")
        self.bundled_update = os.path.join(self.module_dir, "bin", "amrfinder_update")
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)
        self.bundled_database = self._get_latest_database()
        if self.bundled_database is None:
            self.logger.warning("No AMRfinderPlus database found. Run with --update-db to download.")
        db_version = self._get_database_version() if self.bundled_database else "Unknown"

        self.metadata = {
            "tool_name": "EcoliTyper AMRfinderPlus",
            "version": "1.4.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "amrfinder_version": "4.2.7",
            "database_version": db_version,
        }

        self.high_risk_genes = {
            'blaCTX-M-14', 'blaCTX-M-1', 'blaTEM-1', 'blaEC', 'blaCTX-M', 'blaTEM', 'blaSHV',
            'mcr-1.1', 'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', '(Col)mcr-1.1', 'MCR',
            'blaKPC', 'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP',
            'qnrA', 'qnrB', 'qnrC', 'qnrD', 'qnrS', 'qnrVC',
            'aac(3)-IId', 'aac(6\')-Ib-cr', 'aadA1', 'aadA2',
            'aph(3\'\')-Ib', 'aph(3\')-Ia', 'aph(6)-Id',
            'tet(A)', 'tet(C)', 'sul1', 'sul2', 'sul3',
            'dfrA1', 'dfrA14', 'catA2', 'cmlA', 'floR',
            'ermA', 'ermB', 'ermC', 'mphA',
            'acrF', 'emrD', 'emrE', 'mdtM',
            'armA', 'rmtB', 'cfr', 'optrA', 'poxtA', 'CTX',
        }

        self.critical_risk_genes = {
            'mcr-1.1', 'mcr-1', 'blaCTX-M-14', 'blaCTX-M-1', 'blaKPC',
            'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP', 'cfr',
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
            {"text": "The good thing about science is that it's true whether or not you believe in it.",
             "author": "Neil deGrasse Tyson"},
            {"text": "Science knows no country, because knowledge belongs to humanity.",
             "author": "Louis Pasteur"},
            {"text": "EcoliTyper transforms raw E. coli sequences into actionable genomic surveillance.",
             "author": "Brown Beckley"},
            {"text": "Every E. coli genome carries a story; the tools we build determine whether that story is heard.",
             "author": "Brown Beckley"},
        ]

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
                optimal_cpus = max(16, total_physical_cores - 1)
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
        self.logger.info("Strategy: process multiple samples concurrently with optimal core allocation")

    def _get_latest_database(self) -> Optional[str]:
        """Find the latest dated database folder under data/amrfinder_db/."""
        db_root = os.path.join(self.module_dir, "data", "amrfinder_db")
        if not os.path.exists(db_root):
            self.logger.warning("Database root directory not found: %s", db_root)
            return None
        candidates = []
        for item in os.listdir(db_root):
            full_path = os.path.join(db_root, item)
            if os.path.isdir(full_path) and item.startswith('20'):
                candidates.append(item)
        if not candidates:
            self.logger.warning("No database folder starting with '20' found.")
            return None
        latest = sorted(candidates)[-1]
        latest_path = os.path.join(db_root, latest)
        self.logger.info("Using latest database: %s", latest_path)
        return latest_path

    def _get_database_version(self) -> str:
        """Read version.txt from the database folder, or fall back to folder name."""
        if not self.bundled_database:
            return "Unknown"
        version_file = os.path.join(self.bundled_database, "version.txt")
        if os.path.exists(version_file):
            with open(version_file, 'r') as f:
                return f.read().strip()
        return os.path.basename(self.bundled_database)

    def update_database(self) -> bool:
        """Download the latest AMRfinderPlus database using the bundled updater."""
        if not os.path.exists(self.bundled_update):
            self.logger.error("amrfinder_update not found at %s", self.bundled_update)
            return False
        if not os.access(self.bundled_update, os.X_OK):
            self.logger.warning("amrfinder_update not executable, fixing permissions...")
            os.chmod(self.bundled_update, 0o755)
        db_dir = os.path.join(self.module_dir, "data", "amrfinder_db")
        os.makedirs(db_dir, exist_ok=True)
        self.logger.info("Updating AMRfinderPlus database...")
        try:
            cmd = [self.bundled_update, "--database", db_dir]
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.logger.info("Database update completed successfully.")
            self.bundled_database = self._get_latest_database()
            if self.bundled_database:
                self.metadata['database_version'] = self._get_database_version()
                self.logger.info("New database version: %s", self.metadata['database_version'])
                return True
            self.logger.error("Database update succeeded but no database folder found.")
            return False
        except subprocess.CalledProcessError as e:
            self.logger.error("Database update failed: %s", e)
            self.logger.error("STDERR: %s", e.stderr)
            return False

    def check_amrfinder_installed(self) -> bool:
        """Check that the bundled AMRfinderPlus binary and database are available."""
        try:
            if not os.path.exists(self.bundled_amrfinder):
                self.logger.error("Bundled AMRfinderPlus not found at: %s", self.bundled_amrfinder)
                return False
            if not os.access(self.bundled_amrfinder, os.X_OK):
                self.logger.warning("Bundled AMRfinderPlus not executable, fixing permissions...")
                os.chmod(self.bundled_amrfinder, 0o755)
            result = subprocess.run(
                [self.bundled_amrfinder, '--version'],
                capture_output=True, text=True, check=True
            )
            self.logger.info("Bundled AMRfinderPlus version: %s", result.stdout.strip())
            if self.bundled_database and os.path.exists(self.bundled_database):
                self.logger.info("Bundled database found: %s", self.bundled_database)
                db_version_file = os.path.join(self.bundled_database, "version.txt")
                if os.path.exists(db_version_file):
                    with open(db_version_file, 'r') as f:
                        self.logger.info("Database version: %s", f.read().strip())
                return True
            self.logger.warning("Bundled database not found at expected location.")
            self.logger.info("Run with --update-db to download the latest database.")
            return False
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.logger.error("Bundled AMRfinderPlus check failed: %s", e)
            return False

    def run_amrfinder_single_genome(self, genome_file: str, output_dir: str,
                                    min_identity: float = None,
                                    min_coverage: float = None,
                                    report_mutations: bool = True) -> Dict[str, Any]:
        """Run AMRfinderPlus on one genome. Emits exactly 2 files per sample:
             - {genome}_amrfinder_raw.tsv  (AMR + MUTATION merged)
             - {genome}_amrfinder_report.html  (comprehensive AMR + mutations)
        """
        genome_name = Path(genome_file).stem
        merged_raw_file = os.path.join(output_dir, f"{genome_name}_amrfinder_raw.tsv")

        if not os.path.exists(self.bundled_amrfinder):
            self.logger.error("Bundled AMRfinderPlus not found at: %s", self.bundled_amrfinder)
            return {
                'genome': genome_name, 'output_file': merged_raw_file,
                'hits': [], 'hit_count': 0,
                'mutations': [], 'mutation_count': 0,
                'status': 'failed', 'error': 'AMRfinder binary not found',
            }

        with tempfile.TemporaryDirectory() as temp_dir:
            amr_output_file = os.path.join(temp_dir, f"{genome_name}_amrfinder.txt")
            mut_file = os.path.join(temp_dir, f"{genome_name}_mutations.tsv") if report_mutations else None

            cmd = [
                self.bundled_amrfinder,
                '-n', genome_file,
                '-O', 'Escherichia',
                '--output', amr_output_file,
            ]

            if self.bundled_database and os.path.exists(self.bundled_database):
                cmd.extend(['--database', self.bundled_database])
                self.logger.info("Using bundled database: %s", self.bundled_database)
            else:
                self.logger.warning("Using default AMRfinderPlus database location")

            if min_identity is not None:
                cmd.extend(['--ident_min', str(min_identity)])
                self.logger.info("Using minimum identity: %s", min_identity)
            if min_coverage is not None:
                cmd.extend(['--coverage_min', str(min_coverage)])
                self.logger.info("Using minimum coverage: %s", min_coverage)

            if report_mutations and mut_file:
                cmd.extend(['--mutation_all', mut_file])
                self.logger.info("Will report point mutations to temp file (merged into %s)", merged_raw_file)

            self.logger.info("Running AMRfinderPlus on %s (using %d CPU cores)", genome_name, self.cpus)
            self.logger.debug("Command: %s", ' '.join(cmd))

            try:
                subprocess.run(cmd, capture_output=True, text=True, check=True)

                hits = self._parse_amrfinder_output(amr_output_file)
                mutations = (self._parse_mutations_file(mut_file)
                             if mut_file and os.path.exists(mut_file) else [])

                self._write_merged_raw_tsv(amr_output_file, mut_file, merged_raw_file)
                self._create_amrfinder_html_report(genome_name, hits, mutations, output_dir)

                return {
                    'genome': genome_name,
                    'output_file': merged_raw_file,
                    'hits': hits,
                    'hit_count': len(hits),
                    'mutations': mutations,
                    'mutation_count': len(mutations),
                    'status': 'success',
                }

            except subprocess.CalledProcessError as e:
                self.logger.error("AMRfinderPlus failed for %s: %s", genome_name, e.stderr)
                return {
                    'genome': genome_name,
                    'output_file': merged_raw_file,
                    'hits': [], 'hit_count': 0,
                    'mutations': [], 'mutation_count': 0,
                    'status': 'failed',
                }

    def _write_merged_raw_tsv(self, amr_file: str, mut_file: Optional[str], output_file: str):
        """Merge AMRFinderPlus main output and --mutation_all output into one TSV.

        A leading 'Source' column distinguishes AMR gene hits from MUTATION rows.
        Lines appearing identically in both files are written only once (as AMR).
        """
        amr_lines_seen = set()

        with open(output_file, 'w') as out:
            if amr_file and os.path.exists(amr_file):
                with open(amr_file, 'r') as f:
                    amr_lines = f.readlines()
                if amr_lines:
                    header = amr_lines[0].rstrip('\n')
                    out.write('Source\t' + header + '\n')
                    for line in amr_lines[1:]:
                        if not line.strip():
                            continue
                        amr_lines_seen.add(line.rstrip('\n'))
                        out.write('AMR\t' + line)

            if mut_file and os.path.exists(mut_file):
                with open(mut_file, 'r') as f:
                    mut_lines = f.readlines()
                if mut_lines:
                    for line in mut_lines[1:]:
                        if not line.strip():
                            continue
                        if line.rstrip('\n') in amr_lines_seen:
                            continue
                        out.write('MUTATION\t' + line)

    def _parse_amrfinder_output(self, amrfinder_file: str) -> List[Dict]:
        """Parse an AMRFinderPlus tabular output file into structured dicts."""
        hits = []
        try:
            with open(amrfinder_file, 'r') as f:
                lines = f.readlines()
            if not lines or len(lines) < 2:
                return hits
            headers = lines[0].strip().split('\t')
            for line_num, line in enumerate(lines[1:], 2):
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < len(headers):
                    self.logger.warning(
                        "Line %d has %d parts, expected %d",
                        line_num, len(parts), len(headers)
                    )
                    continue
                hit = {}
                for i, header in enumerate(headers):
                    hit[header] = parts[i] if i < len(parts) else ''
                processed_hit = {
                    'protein_id': hit.get('Protein id', ''),
                    'contig_id': hit.get('Contig id', ''),
                    'start': hit.get('Start', ''),
                    'stop': hit.get('Stop', ''),
                    'strand': hit.get('Strand', ''),
                    'gene_symbol': hit.get('Element symbol', ''),
                    'sequence_name': hit.get('Element name', ''),
                    'scope': hit.get('Scope', ''),
                    'element_type': hit.get('Type', ''),
                    'element_subtype': hit.get('Subtype', ''),
                    'class': hit.get('Class', ''),
                    'subclass': hit.get('Subclass', ''),
                    'method': hit.get('Method', ''),
                    'target_length': hit.get('Target length', ''),
                    'ref_length': hit.get('Reference sequence length', ''),
                    'coverage': hit.get('% Coverage of reference', '').replace('%', ''),
                    'identity': hit.get('% Identity to reference', '').replace('%', ''),
                    'alignment_length': hit.get('Alignment length', ''),
                    'accession': hit.get('Closest reference accession', ''),
                    'closest_name': hit.get('Closest reference name', ''),
                    'hmm_id': hit.get('HMM accession', ''),
                    'hmm_description': hit.get('HMM description', ''),
                }
                hits.append(processed_hit)
        except Exception as e:
            self.logger.error("Error parsing %s: %s", amrfinder_file, e)
        self.logger.info("Parsed %d hits from %s", len(hits), amrfinder_file)
        return hits

    def _parse_mutations_file(self, mut_file: str) -> List[Dict]:
        """Parse the mutations TSV file (same format as the main AMR output)."""
        return self._parse_amrfinder_output(mut_file)

    def _analyze_ecoli_amr_results(self, hits: List[Dict]) -> Dict[str, Any]:
        """Aggregate AMR hits into a risk and mechanism profile."""
        analysis = {
            'total_genes': len(hits),
            'resistance_classes': {},
            'total_classes': 0,
            'high_risk_genes': 0,
            'critical_risk_genes': 0,
            'high_risk_list': [],
            'critical_risk_list': [],
            'resistance_mechanisms': {
                'esbl': [],
                'carbapenemase': [],
                'colistin_resistance': [],
                'fluoroquinolone_resistance': [],
                'aminoglycoside_resistance': [],
                'efflux_pumps': [],
                'other_amr': [],
            },
        }
        for hit in hits:
            gene_symbol = hit.get('gene_symbol', '')
            resistance_class = hit.get('class', '')
            self._categorize_resistance_mechanism(gene_symbol, resistance_class, analysis)
            if gene_symbol in self.critical_risk_genes:
                analysis['critical_risk_genes'] += 1
                if gene_symbol not in analysis['critical_risk_list']:
                    analysis['critical_risk_list'].append(gene_symbol)
            if gene_symbol in self.high_risk_genes:
                analysis['high_risk_genes'] += 1
                if gene_symbol not in analysis['high_risk_list']:
                    analysis['high_risk_list'].append(gene_symbol)
            if resistance_class:
                if resistance_class not in analysis['resistance_classes']:
                    analysis['resistance_classes'][resistance_class] = []
                if gene_symbol not in analysis['resistance_classes'][resistance_class]:
                    analysis['resistance_classes'][resistance_class].append(gene_symbol)
        analysis['total_classes'] = len(analysis['resistance_classes'])
        return analysis

    def _categorize_resistance_mechanism(self, gene_symbol: str, resistance_class: str,
                                         analysis: Dict[str, Any]):
        """Assign a gene to one of the resistance mechanism buckets."""
        esbl_genes = {'blaCTX-M', 'blaTEM', 'blaSHV', 'blaCTX-M-14', 'blaCTX-M-1', 'blaTEM-1'}
        carbapenemase_genes = {'blaKPC', 'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP'}
        colistin_genes = {'mcr-1.1', 'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', '(Col)mcr-1.1'}
        fluoroquinolone_genes = {'qnrA', 'qnrB', 'qnrC', 'qnrD', 'qnrS', 'qnrVC'}
        aminoglycoside_genes = {
            'aac(3)-IId', 'aac(6\')-Ib-cr', 'aadA1', 'aadA2',
            'aph(3\'\')-Ib', 'aph(3\')-Ia', 'aph(6)-Id',
        }
        efflux_pump_genes = {'acrF', 'emrD', 'emrE', 'mdtM'}

        if gene_symbol in esbl_genes:
            analysis['resistance_mechanisms']['esbl'].append(gene_symbol)
        elif gene_symbol in carbapenemase_genes:
            analysis['resistance_mechanisms']['carbapenemase'].append(gene_symbol)
        elif gene_symbol in colistin_genes:
            analysis['resistance_mechanisms']['colistin_resistance'].append(gene_symbol)
        elif gene_symbol in fluoroquinolone_genes:
            analysis['resistance_mechanisms']['fluoroquinolone_resistance'].append(gene_symbol)
        elif gene_symbol in aminoglycoside_genes:
            analysis['resistance_mechanisms']['aminoglycoside_resistance'].append(gene_symbol)
        elif gene_symbol in efflux_pump_genes:
            analysis['resistance_mechanisms']['efflux_pumps'].append(gene_symbol)
        else:
            analysis['resistance_mechanisms']['other_amr'].append(gene_symbol)

    def _get_shared_css(self) -> str:
        """Return the shared CSS block used by all HTML reports."""
        return """
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff; padding: 20px; min-height: 100vh;
        }
        .container { max-width: 1400px; margin: 0 auto; }
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
        .metric-card.critical-card {
            background: linear-gradient(135deg, #dc2626 0%, #b91c1c 100%);
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
        .stat-card.critical-card {
            background: linear-gradient(135deg, #dc2626 0%, #b91c1c 100%);
        }
        .stat-value { font-size: 24px; font-weight: bold; margin-bottom: 5px; }
        .stat-label { font-size: 12px; opacity: 0.9; }
        .percentage-badge {
            display: inline-block; background: #e3f2fd; color: #1976d2;
            padding: 2px 8px; border-radius: 12px; font-size: 0.85em;
            margin-left: 5px; border: 1px solid #bbdefb;
        }
        .frequency-cell { font-weight: 600; color: #2c3e50; }
        .mech-item {
            margin: 10px 0; padding: 12px; border-radius: 6px;
            color: #1f2937; font-size: 14px;
        }
        .mech-esbl { background: #fff3cd; }
        .mech-carba { background: #fee2e2; }
        .mech-colistin { background: #fee2e2; }
        .mech-fluoro { background: #dbeafe; }
        .mech-amino { background: #dbeafe; }
        .mech-efflux { background: #e5e7eb; }
        .mech-other { background: #f3f4f6; }
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
                    Powered by AMRFinderPlus {self.metadata['amrfinder_version']}
                    (database {self.metadata['database_version']}).
                </p>
            </div>
        </div>
        """

    def _build_critical_alert_html(self, analysis: Dict[str, Any]) -> str:
        """Return the critical-risk alert card, or an empty string if none."""
        if analysis['critical_risk_genes'] == 0:
            return ""
        badges = "".join(f'<span class="risk-badge">{g}</span>'
                         for g in analysis['critical_risk_list'])
        return f"""
        <div class="report-section" style="border-left: 4px solid #dc2626;">
            <h2 style="color: #dc2626;">CRITICAL RESISTANCE ALERTS</h2>
            <p><strong>{analysis['critical_risk_genes']} critical risk genes found:</strong>
            last-resort antibiotic resistance detected.</p>
            <div style="margin: 10px 0;">{badges}</div>
        </div>
        """

    def _build_high_risk_alert_html(self, analysis: Dict[str, Any]) -> str:
        """Return the high-risk alert card (only when no critical genes present)."""
        if analysis['high_risk_genes'] == 0 or analysis['critical_risk_genes'] > 0:
            return ""
        badges = "".join(f'<span class="warning-badge">{g}</span>'
                         for g in analysis['high_risk_list'])
        return f"""
        <div class="report-section" style="border-left: 4px solid #f59e0b;">
            <h2 style="color: #f59e0b;">HIGH-RISK AMR GENES DETECTED</h2>
            <p><strong>{analysis['high_risk_genes']} high-risk genes found:</strong></p>
            <div style="margin: 10px 0;">{badges}</div>
        </div>
        """

    def _build_mechanism_html(self, analysis: Dict[str, Any]) -> str:
        """Return the resistance mechanism breakdown card, or an empty string."""
        mech = analysis['resistance_mechanisms']
        blocks = []
        if mech['esbl']:
            blocks.append(f'<div class="mech-item mech-esbl"><strong>ESBL Genes:</strong> '
                          f'{", ".join(mech["esbl"])}</div>')
        if mech['carbapenemase']:
            blocks.append(f'<div class="mech-item mech-carba"><strong>Carbapenemase Genes (CRITICAL):</strong> '
                          f'{", ".join(mech["carbapenemase"])}</div>')
        if mech['colistin_resistance']:
            blocks.append(f'<div class="mech-item mech-colistin"><strong>Colistin Resistance (CRITICAL):</strong> '
                          f'{", ".join(mech["colistin_resistance"])}</div>')
        if mech['fluoroquinolone_resistance']:
            blocks.append(f'<div class="mech-item mech-fluoro"><strong>Fluoroquinolone Resistance:</strong> '
                          f'{", ".join(mech["fluoroquinolone_resistance"])}</div>')
        if mech['aminoglycoside_resistance']:
            blocks.append(f'<div class="mech-item mech-amino"><strong>Aminoglycoside Resistance:</strong> '
                          f'{", ".join(mech["aminoglycoside_resistance"])}</div>')
        if mech['efflux_pumps']:
            blocks.append(f'<div class="mech-item mech-efflux"><strong>Efflux Pumps:</strong> '
                          f'{", ".join(mech["efflux_pumps"])}</div>')
        if not blocks:
            return ""
        return f"""
        <div class="report-section">
            <h2>Resistance Mechanism Breakdown</h2>
            {''.join(blocks)}
        </div>
        """

    def _build_resistance_classes_html(self, analysis: Dict[str, Any]) -> str:
        """Return the resistance classes card, or an empty string."""
        if not analysis['resistance_classes']:
            return ""
        rows = ""
        for class_name, genes in analysis['resistance_classes'].items():
            rows += f"""
                        <tr class="present">
                            <td><strong>{class_name}</strong></td>
                            <td>{len(genes)}</td>
                            <td>{", ".join(genes)}</td>
                        </tr>
            """
        return f"""
        <div class="report-section">
            <h2>Resistance Classes Detected</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr><th>Resistance Class</th><th>Gene Count</th><th>Genes</th></tr>
                    </thead>
                    <tbody>{rows}</tbody>
                </table>
            </div>
        </div>
        """

    def _build_amr_table_html(self, analysis: Dict[str, Any], hits: List[Dict]) -> str:
        """Return the detailed AMR genes table card, or a no-hits card."""
        if not hits:
            return """
        <div class="report-section">
            <h2>No AMR Genes Detected</h2>
            <p>No antimicrobial resistance genes found in this E. coli genome.</p>
        </div>
        """
        rows = ""
        for hit in hits:
            row_class = "present"
            gene_symbol = hit.get('gene_symbol', '')
            if gene_symbol in analysis['critical_risk_list']:
                row_class = "critical"
            elif gene_symbol in analysis['high_risk_list']:
                row_class = "high-risk"
            rows += f"""
                        <tr class="{row_class}">
                            <td><strong>{gene_symbol}</strong></td>
                            <td class="product-cell">{hit.get('sequence_name', '')}</td>
                            <td>{hit.get('class', '')}</td>
                            <td>{hit.get('subclass', '')}</td>
                            <td>{hit.get('coverage', '')}%</td>
                            <td>{hit.get('identity', '')}%</td>
                            <td>{hit.get('scope', '')}</td>
                        </tr>
            """
        return f"""
        <div class="report-section">
            <h2>Detailed AMR Genes Detected</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Gene Symbol</th><th>Sequence Name</th><th>Class</th>
                            <th>Subclass</th><th>Coverage</th><th>Identity</th><th>Scope</th>
                        </tr>
                    </thead>
                    <tbody>{rows}</tbody>
                </table>
            </div>
        </div>
        """

    def _build_mutations_table_html(self, mutations: List[Dict]) -> str:
        """Return the point mutations table card, or an empty string."""
        if not mutations:
            return ""
        rows = ""
        for m in mutations:
            rows += f"""
                        <tr class="present">
                            <td><strong>{m.get('gene_symbol', '')}</strong></td>
                            <td class="product-cell">{m.get('sequence_name', '')}</td>
                            <td>{m.get('class', '')}</td>
                            <td>{m.get('subclass', '')}</td>
                            <td>{m.get('contig_id', '')}</td>
                            <td>{m.get('start', '')}</td>
                            <td>{m.get('stop', '')}</td>
                            <td>{m.get('strand', '')}</td>
                            <td>{m.get('coverage', '')}</td>
                            <td>{m.get('identity', '')}</td>
                            <td>{m.get('accession', '')}</td>
                        </tr>
            """
        return f"""
        <div class="report-section">
            <h2>Point Mutations Detected</h2>
            <p>All point mutations detected by AMRFinderPlus (including synonymous variants).</p>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Gene Symbol</th><th>Mutation</th><th>Class</th><th>Subclass</th>
                            <th>Contig</th><th>Start</th><th>Stop</th><th>Strand</th>
                            <th>Coverage (%)</th><th>Identity (%)</th><th>Accession</th>
                        </tr>
                    </thead>
                    <tbody>{rows}</tbody>
                </table>
            </div>
        </div>
        """

    def _create_amrfinder_html_report(self, genome_name: str, hits: List[Dict],
                                      mutations: List[Dict], output_dir: str):
        """Write the combined per-sample AMRFinderPlus HTML report."""
        analysis = self._analyze_ecoli_amr_results(hits)
        random_quote = random.choice(self.science_quotes)

        critical_card = self._build_critical_alert_html(analysis)
        high_risk_card = self._build_high_risk_alert_html(analysis)
        mech_card = self._build_mechanism_html(analysis)
        classes_card = self._build_resistance_classes_html(analysis)
        amr_table = self._build_amr_table_html(analysis, hits)
        mutations_table = self._build_mutations_table_html(mutations)

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliTyper AMRfinderPlus - Report - {genome_name}</title>
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
                    <div class="metric-label">Total AMR Genes</div>
                    <div class="metric-value">{analysis['total_genes']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">High Risk Genes</div>
                    <div class="metric-value">{analysis['high_risk_genes']}</div>
                </div>
                <div class="metric-card critical-card">
                    <div class="metric-label">Critical Risk</div>
                    <div class="metric-value">{analysis['critical_risk_genes']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Point Mutations</div>
                    <div class="metric-value">{len(mutations)}</div>
                </div>
            </div>
            <p style="margin-top:15px;"><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>AMRfinderPlus:</strong> {self.metadata['amrfinder_version']}</p>
            <p><strong>Database:</strong> {self.metadata['database_version']}</p>
            <p><strong>Analysis Date:</strong> {datetime.now().strftime('%Y-%m-%d')}</p>
        </div>

        {critical_card}
        {high_risk_card}
        {mech_card}
        {classes_card}
        {amr_table}
        {mutations_table}

        {self._get_footer_html("EcoliTyper - AMRFinderPlus Per-Sample Report")}
    </div>
    {self._get_quote_js()}
</body>
</html>"""

        html_file = os.path.join(output_dir, f"{genome_name}_amrfinder_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("AMRfinderPlus HTML report generated: %s", html_file)

    def create_amr_summary(self, all_results: Dict[str, Any], output_base: str):
        """Write the batch AMR TSV, statistics TSV, master JSON, and summary HTML."""
        self.logger.info("Creating E. coli AMR summary files and HTML reports...")

        summary_file = os.path.join(output_base, "ecoli_amrfinder_summary.tsv")
        with open(summary_file, 'w') as f:
            f.write("Genome\tProtein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\t"
                    "Scope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t"
                    "% Coverage of reference\t% Identity to reference\tAlignment length\t"
                    "Closest reference accession\tClosest reference name\tHMM accession\tHMM description\n")
            for genome_name, result in all_results.items():
                for hit in result['hits']:
                    row = [
                        genome_name,
                        hit.get('protein_id', ''),
                        hit.get('contig_id', ''),
                        hit.get('start', ''),
                        hit.get('stop', ''),
                        hit.get('strand', ''),
                        hit.get('gene_symbol', ''),
                        hit.get('sequence_name', ''),
                        hit.get('scope', ''),
                        hit.get('element_type', ''),
                        hit.get('element_subtype', ''),
                        hit.get('class', ''),
                        hit.get('subclass', ''),
                        hit.get('method', ''),
                        hit.get('target_length', ''),
                        hit.get('ref_length', ''),
                        hit.get('coverage', ''),
                        hit.get('identity', ''),
                        hit.get('alignment_length', ''),
                        hit.get('accession', ''),
                        hit.get('closest_name', ''),
                        hit.get('hmm_id', ''),
                        hit.get('hmm_description', ''),
                    ]
                    f.write('\t'.join(str(x) for x in row) + '\n')
        self.logger.info("AMR summary file created: %s", summary_file)

        stats_file = os.path.join(output_base, "ecoli_amrfinder_statistics_summary.tsv")
        with open(stats_file, 'w') as f:
            f.write("Genome\tTotal_AMR_Genes\tHigh_Risk_Genes\tCritical_Risk_Genes\t"
                    "Resistance_Classes\tGene_List\n")
            for genome_name, result in all_results.items():
                genes = list(set(hit.get('gene_symbol', '')
                                 for hit in result['hits'] if hit.get('gene_symbol')))
                high_risk_count = sum(1 for g in genes if g in self.high_risk_genes)
                critical_risk_count = sum(1 for g in genes if g in self.critical_risk_genes)
                classes = list(set(hit.get('class', '')
                                   for hit in result['hits'] if hit.get('class')))
                f.write(f"{genome_name}\t{result['hit_count']}\t{high_risk_count}\t"
                        f"{critical_risk_count}\t{','.join(classes)}\t{','.join(genes)}\n")
        self.logger.info("AMR statistics summary created: %s", stats_file)

        self.create_json_summaries(all_results, output_base)
        self._create_summary_html_report(all_results, output_base)
        self.create_mutation_summary(all_results, output_base)

    def create_mutation_summary(self, all_results: Dict[str, Any], output_base: str):
        """Write mutation batch TSV, HTML, and JSON summaries."""
        self.logger.info("Creating mutation batch summaries...")
        all_mutations = []
        genome_mutation_counts = {}
        for genome_name, result in all_results.items():
            muts = result.get('mutations', []) or []
            genome_mutation_counts[genome_name] = len(muts)
            for m in muts:
                m_copy = m.copy()
                m_copy['genome'] = genome_name
                all_mutations.append(m_copy)

        if not all_mutations:
            self.logger.info("No mutations found in any genome; skipping mutation summaries.")
            return

        tsv_file = os.path.join(output_base, "mutation_summary.tsv")
        with open(tsv_file, 'w') as f:
            fields = ['genome', 'gene_symbol', 'sequence_name', 'class', 'subclass',
                      'contig_id', 'start', 'stop', 'strand', 'coverage', 'identity', 'accession']
            f.write('\t'.join(fields) + '\n')
            for m in all_mutations:
                f.write('\t'.join(str(m.get(k, '')) for k in fields) + '\n')
        self.logger.info("Mutation TSV: %s", tsv_file)

        self._create_mutation_summary_html(all_mutations, genome_mutation_counts, output_base)
        self._create_mutation_json_summaries(all_mutations, genome_mutation_counts, output_base)

    def _create_mutation_summary_html(self, all_mutations: List[Dict],
                                      genome_counts: Dict[str, int], output_base: str):
        """Write the batch mutation HTML summary."""
        gene_freq = {}
        for m in all_mutations:
            gene = m.get('gene_symbol', 'unknown')
            mutation = m.get('sequence_name', '')
            key = f"{gene}_{mutation}" if mutation else gene
            if key not in gene_freq:
                gene_freq[key] = {'count': 0, 'genomes': set(), 'gene': gene,
                                  'mutation': mutation,
                                  'class': m.get('class', ''),
                                  'subclass': m.get('subclass', '')}
            gene_freq[key]['count'] += 1
            gene_freq[key]['genomes'].add(m.get('genome', ''))
        for k in gene_freq:
            gene_freq[k]['genomes'] = ', '.join(sorted(gene_freq[k]['genomes']))
        sorted_freq = sorted(gene_freq.values(), key=lambda x: x['count'], reverse=True)

        random_quote = random.choice(self.science_quotes)
        total_genomes_with_mutations = sum(1 for c in genome_counts.values() if c > 0)

        rows = ""
        for item in sorted_freq:
            rows += f"""
                        <tr>
                            <td><strong>{item['gene']}</strong></td>
                            <td class="product-cell">{item['mutation']}</td>
                            <td>{item['count']}</td>
                            <td>{item['genomes']}</td>
                            <td>{item['class']}</td>
                            <td>{item['subclass']}</td>
                        </tr>
            """

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliTyper - Mutation Batch Summary</title>
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
            <h2>Mutation Summary Across All Genomes</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{len(all_mutations)}</div>
                    <div class="stat-label">TOTAL MUTATION EVENTS</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{total_genomes_with_mutations}</div>
                    <div class="stat-label">GENOMES WITH MUTATIONS</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{len(sorted_freq)}</div>
                    <div class="stat-label">UNIQUE MUTATIONS</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{datetime.now().strftime('%Y-%m-%d')}</div>
                    <div class="stat-label">ANALYSIS DATE</div>
                </div>
            </div>
        </div>

        <div class="report-section">
            <h2>Mutation Frequency by Gene / Mutation</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Gene</th><th>Mutation</th><th>Count</th>
                            <th>Genomes</th><th>Class</th><th>Subclass</th>
                        </tr>
                    </thead>
                    <tbody>{rows}</tbody>
                </table>
            </div>
        </div>

        {self._get_footer_html("EcoliTyper - AMRFinderPlus Mutation Batch Summary")}
    </div>
    {self._get_quote_js()}
</body>
</html>"""

        out_file = os.path.join(output_base, "mutation_summary.html")
        with open(out_file, 'w') as f:
            f.write(html)
        self.logger.info("Mutation HTML summary: %s", out_file)

    def _create_mutation_json_summaries(self, all_mutations: List[Dict],
                                        genome_counts: Dict[str, int], output_base: str):
        """Write the batch mutation JSON summary."""
        genome_summary = {g: {'total_mutations': c} for g, c in genome_counts.items()}
        gene_mutation_map = defaultdict(lambda: {'count': 0, 'genomes': set(), 'details': []})
        for m in all_mutations:
            gene = m.get('gene_symbol', 'unknown')
            mut_name = m.get('sequence_name', '')
            key = f"{gene}_{mut_name}"
            gene_mutation_map[key]['count'] += 1
            gene_mutation_map[key]['genomes'].add(m.get('genome', ''))
            gene_mutation_map[key]['details'].append({
                'genome': m.get('genome'), 'gene': gene, 'mutation': mut_name,
                'class': m.get('class'), 'subclass': m.get('subclass'),
                'contig': m.get('contig_id'), 'start': m.get('start'), 'stop': m.get('stop'),
            })
        for v in gene_mutation_map.values():
            v['genomes'] = list(v['genomes'])
        master_json = {
            'metadata': {
                'tool': 'EcoliTyper AMRfinderPlus Mutation Module',
                'version': self.metadata['version'],
                'analysis_date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                'total_genomes_analyzed': len(genome_counts),
                'total_mutations_detected': len(all_mutations),
            },
            'genome_summary': genome_summary,
            'mutation_frequency': {k: {'count': v['count'], 'genomes': v['genomes']}
                                   for k, v in gene_mutation_map.items()},
            'all_mutations': all_mutations,
        }
        json_file = os.path.join(output_base, "mutation_master_summary.json")
        with open(json_file, 'w') as f:
            json.dump(master_json, f, indent=2, default=str)
        self.logger.info("Mutation master JSON: %s", json_file)

    def create_json_summaries(self, all_results: Dict[str, Any], output_base: str):
        """Write the single master JSON summarising every genome's AMR profile."""
        self.logger.info("Creating master JSON summary...")

        master_summary = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'amrfinder_version': self.metadata['amrfinder_version'],
                'database_version': self.metadata['database_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_genomes': len(all_results),
            },
            'genome_summaries': {},
            'cross_genome_patterns': {},
        }

        all_hits_by_gene = defaultdict(lambda: {'count': 0, 'genomes': set()})
        genomes_with_critical = 0
        genomes_with_high_risk = 0

        for genome_name, result in all_results.items():
            hits = result['hits']
            genes = [hit.get('gene_symbol', '') for hit in hits if hit.get('gene_symbol', '')]
            unique_genes = set(genes)
            critical_genes = [g for g in unique_genes if g in self.critical_risk_genes]
            high_risk_genes = [g for g in unique_genes
                               if g in self.high_risk_genes and g not in self.critical_risk_genes]
            if critical_genes:
                genomes_with_critical += 1
            if high_risk_genes:
                genomes_with_high_risk += 1

            master_summary['genome_summaries'][genome_name] = {
                'total_hits': result['hit_count'],
                'mutation_count': result.get('mutation_count', 0),
                'unique_genes': len(unique_genes),
                'critical_genes': critical_genes,
                'high_risk_genes': high_risk_genes,
                'genes': list(unique_genes),
                'status': result['status'],
            }

            for gene in unique_genes:
                all_hits_by_gene[gene]['count'] += 1
                all_hits_by_gene[gene]['genomes'].add(genome_name)

        cross_genome_data = {}
        for gene, data in all_hits_by_gene.items():
            cross_genome_data[gene] = {
                'frequency': data['count'],
                'genomes': list(data['genomes']),
                'risk_level': ('CRITICAL' if gene in self.critical_risk_genes
                               else 'HIGH' if gene in self.high_risk_genes
                               else 'STANDARD'),
            }

        master_summary['cross_genome_patterns'] = {
            'total_unique_genes': len(all_hits_by_gene),
            'genomes_with_critical': genomes_with_critical,
            'genomes_with_high_risk': genomes_with_high_risk,
            'gene_frequency': cross_genome_data,
        }

        master_json_file = os.path.join(output_base, "ecoli_amrfinder_master_summary.json")
        with open(master_json_file, 'w') as f:
            json.dump(master_summary, f, indent=2, default=str)
        self.logger.info("Master JSON summary created: %s", master_json_file)

    def _create_summary_html_report(self, all_results: Dict[str, Any], output_base: str):
        """Write the batch HTML summary across all genomes."""
        total_genomes = len(all_results)
        total_hits = sum(r['hit_count'] for r in all_results.values())
        total_mutations = sum(r.get('mutation_count', 0) for r in all_results.values())

        critical_genes_found = set()
        high_risk_genes_found = set()
        genomes_with_critical = 0
        genomes_with_high_risk = 0

        genes_per_genome = {}
        gene_frequency = {}

        for genome_name, result in all_results.items():
            genome_genes = set()
            for hit in result['hits']:
                gene = hit.get('gene_symbol', '')
                if gene:
                    genome_genes.add(gene)
                    if gene not in gene_frequency:
                        gene_frequency[gene] = set()
                    gene_frequency[gene].add(genome_name)
            genes_per_genome[genome_name] = genome_genes

            has_critical = any(g in genome_genes for g in self.critical_risk_genes)
            has_high_risk = any(g in genome_genes for g in self.high_risk_genes)
            if has_critical:
                genomes_with_critical += 1
                critical_genes_found.update(genome_genes & self.critical_risk_genes)
            if has_high_risk:
                genomes_with_high_risk += 1
                high_risk_genes_found.update(genome_genes & self.high_risk_genes)

        random_quote = random.choice(self.science_quotes)

        critical_card = ""
        if critical_genes_found:
            badges = "".join(f'<span class="risk-badge">{g}</span>'
                             for g in sorted(critical_genes_found))
            critical_card = f"""
        <div class="report-section" style="border-left: 4px solid #dc2626;">
            <h2 style="color: #dc2626;">CRITICAL RISK AMR GENES ACROSS ALL GENOMES</h2>
            <p><strong>{len(critical_genes_found)} unique critical genes found in
            {genomes_with_critical} genomes:</strong> immediate attention required.</p>
            <div style="margin: 10px 0;">{badges}</div>
        </div>
        """

        high_risk_card = ""
        if high_risk_genes_found and not critical_genes_found:
            badges = "".join(f'<span class="warning-badge">{g}</span>'
                             for g in sorted(high_risk_genes_found))
            high_risk_card = f"""
        <div class="report-section" style="border-left: 4px solid #f59e0b;">
            <h2 style="color: #f59e0b;">HIGH-RISK AMR GENES ACROSS ALL GENOMES</h2>
            <p><strong>{len(high_risk_genes_found)} unique high-risk genes in
            {genomes_with_high_risk} genomes:</strong></p>
            <div style="margin: 10px 0;">{badges}</div>
        </div>
        """

        genes_by_genome_rows = ""
        for genome in sorted(genes_per_genome.keys()):
            genes = sorted(genes_per_genome.get(genome, set()))
            row_class = "present"
            if any(g in genes for g in self.critical_risk_genes):
                row_class = "critical"
            elif any(g in genes for g in self.high_risk_genes):
                row_class = "high-risk"
            gene_list = ", ".join(genes) if genes else "None"
            genes_by_genome_rows += f"""
                        <tr class="{row_class}">
                            <td><strong>{genome}</strong></td>
                            <td>{len(genes)}</td>
                            <td>{gene_list}</td>
                        </tr>
            """

        frequency_rows = ""
        for gene, genomes in sorted(gene_frequency.items(),
                                    key=lambda x: len(x[1]), reverse=True):
            frequency = len(genomes)
            percentage = (frequency / total_genomes * 100) if total_genomes > 0 else 0
            if gene in self.critical_risk_genes:
                row_class = "critical"
            elif gene in self.high_risk_genes:
                row_class = "high-risk"
            else:
                row_class = "present"
            genome_list = ", ".join(sorted(genomes))
            frequency_rows += f"""
                        <tr class="{row_class}">
                            <td><strong>{gene}</strong></td>
                            <td class="frequency-cell">{frequency}
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
    <title>EcoliTyper AMRfinderPlus - Batch Summary</title>
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
            <h2>Batch AMRfinderPlus Summary</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{total_hits}</div>
                    <div class="stat-label">TOTAL AMR HITS</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{total_genomes}</div>
                    <div class="stat-label">GENOMES ANALYZED</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{len(gene_frequency)}</div>
                    <div class="stat-label">UNIQUE GENES</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{total_mutations}</div>
                    <div class="stat-label">TOTAL MUTATIONS</div>
                </div>
            </div>
            <p style="margin-top:15px;"><strong>AMRfinderPlus:</strong> {self.metadata['amrfinder_version']}</p>
            <p><strong>Database:</strong> {self.metadata['database_version']}</p>
            <p><strong>Analysis Date:</strong> {self.metadata['analysis_date']}</p>
        </div>

        {critical_card}
        {high_risk_card}

        <div class="report-section">
            <h2>Genes by Genome</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr><th>Genome</th><th>Gene Count</th><th>Genes Detected</th></tr>
                    </thead>
                    <tbody>{genes_by_genome_rows}</tbody>
                </table>
            </div>
        </div>

        <div class="report-section">
            <h2>Gene Frequency</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr><th>Gene</th><th>Frequency</th><th>Genomes</th></tr>
                    </thead>
                    <tbody>{frequency_rows}</tbody>
                </table>
            </div>
        </div>

        {self._get_footer_html("EcoliTyper - AMRFinderPlus Batch Summary")}
    </div>
    {self._get_quote_js()}
</body>
</html>"""

        html_file = os.path.join(output_base, "ecoli_amrfinder_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("AMRfinderPlus summary HTML report: %s", html_file)

    def process_single_genome(self, genome_file: str,
                              output_base: str = "ecoli_amrfinder_results",
                              min_identity: float = None,
                              min_coverage: float = None,
                              report_mutations: bool = True) -> Dict[str, Any]:
        """Process one genome with AMRfinderPlus."""
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        self.logger.info("Processing genome: %s", genome_name)
        os.makedirs(results_dir, exist_ok=True)

        if not self.check_amrfinder_installed():
            self.logger.error("Bundled AMRfinderPlus not available!")
            return {
                'genome': genome_name, 'hits': [], 'hit_count': 0,
                'mutations': [], 'mutation_count': 0, 'status': 'failed',
                'error': 'Bundled AMRfinderPlus not available',
            }

        result = self.run_amrfinder_single_genome(
            genome_file, results_dir,
            min_identity=min_identity,
            min_coverage=min_coverage,
            report_mutations=report_mutations,
        )
        self.logger.info("%s: %d AMR hits, %d mutations",
                         genome_name, result['hit_count'], result.get('mutation_count', 0))
        return result

    def process_multiple_genomes(self, genome_pattern: str,
                                 output_base: str = "ecoli_amrfinder_results",
                                 min_identity: float = None,
                                 min_coverage: float = None,
                                 report_mutations: bool = True) -> Dict[str, Any]:
        """Process multiple genomes in parallel and write all reports."""
        print("\n" + "=" * 80)
        print(self._get_ascii_art())
        print("=" * 80)
        print("EcoliTyper AMRfinderPlus")
        print("=" * 80)

        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta",
                          f"{genome_pattern}.fa", f"{genome_pattern}.fna",
                          f"{genome_pattern}.faa"]
        genome_files = []
        for pattern in fasta_patterns:
            genome_files.extend(glob.glob(pattern))
        genome_files = list(set(genome_files))

        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching pattern: {genome_pattern}")

        self.logger.info("Found %d E. coli genomes: %s",
                         len(genome_files), [Path(f).name for f in genome_files])
        os.makedirs(output_base, exist_ok=True)

        all_results = {}
        max_concurrent = max(1, min(self.cpus, len(genome_files),
                                    int(self.available_ram / 1.5)))
        self.logger.info("Using %d concurrent genome processing jobs", max_concurrent)
        self.logger.info("Bundled AMRfinderPlus: %s", self.bundled_amrfinder)
        self.logger.info("Bundled database: %s", self.bundled_database)

        with ThreadPoolExecutor(max_workers=max_concurrent) as executor:
            future_to_genome = {
                executor.submit(self.process_single_genome, genome, output_base,
                                min_identity, min_coverage, report_mutations): genome
                for genome in genome_files
            }
            for future in as_completed(future_to_genome):
                genome = future_to_genome[future]
                try:
                    result = future.result()
                    all_results[result['genome']] = result
                    self.logger.info("Completed: %s (%d AMR hits, %d mutations)",
                                     result['genome'], result['hit_count'],
                                     result.get('mutation_count', 0))
                except Exception as e:
                    self.logger.error("Failed: %s - %s", genome, e)
                    all_results[Path(genome).stem] = {
                        'genome': Path(genome).stem, 'hits': [], 'hit_count': 0,
                        'mutations': [], 'mutation_count': 0, 'status': 'failed',
                    }

        self.create_amr_summary(all_results, output_base)

        self.logger.info("=== E. COLI AMR ANALYSIS COMPLETE ===")
        self.logger.info("Processed %d genomes", len(all_results))
        self.logger.info("Results saved to: %s", output_base)
        self.logger.info("Bundled AMRfinderPlus used: %s", self.bundled_amrfinder)
        self.logger.info("Bundled database used: %s", self.bundled_database)
        return all_results


def main():
    """Command-line entry point for E. coli AMR analysis."""
    parser = argparse.ArgumentParser(
        description='EcoliTyper AMRfinderPlus Analysis - E. coli AMR',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python ecoli_amrfinder.py "*.fna"
  python ecoli_amrfinder.py "*.fna" --min-identity 0.95 --min-coverage 0.9
  python ecoli_amrfinder.py "*.fna" --skip-mutations
  python ecoli_amrfinder.py "*.fa" --cpus 16
  python ecoli_amrfinder.py --update-db
  python ecoli_amrfinder.py --db-version

Supported FASTA extensions: .fasta, .fa, .fna, .faa

PER-GENOME OUTPUTS (2 files each):
  - {genome}_amrfinder_report.html   (AMR genes + point mutations combined)
  - {genome}_amrfinder_raw.tsv       (Source column: AMR / MUTATION)
        """
    )
    parser.add_argument('pattern', nargs='?',
                        help='File pattern for E. coli genomes (e.g., "*.fasta")')
    parser.add_argument('--cpus', '-c', type=int, default=None,
                        help='Number of CPU cores to use (default: auto-detect)')
    parser.add_argument('--output', '-o', default='ecoli_amrfinder_results',
                        help='Output directory (default: ecoli_amrfinder_results)')
    parser.add_argument('--min-identity', type=float, default=None,
                        help='Minimum identity (0..1) for hits')
    parser.add_argument('--min-coverage', type=float, default=None,
                        help='Minimum coverage of reference (0..1)')
    parser.add_argument('--skip-mutations', action='store_true',
                        help='Skip point mutation reporting (reported by default)')
    parser.add_argument('--update-db', action='store_true',
                        help='Update AMRfinderPlus database and exit')
    parser.add_argument('--db-version', action='store_true',
                        help='Show current database version and exit')
    args = parser.parse_args()

    if args.update_db or args.db_version:
        executor = EcoliAMRfinderPlus(cpus=args.cpus)
        if args.update_db:
            print("Updating AMRfinderPlus database...")
            success = executor.update_database()
            print("Database updated successfully." if success else "Database update failed.")
            sys.exit(0)
        if args.db_version:
            print(f"Database version: {executor.metadata['database_version']}")
            print(f"Database path: {executor.bundled_database or 'Not found'}")
            sys.exit(0)

    if not args.pattern:
        parser.error("Please provide a file pattern for genomes (or use --update-db / --db-version)")

    executor = EcoliAMRfinderPlus(cpus=args.cpus)

    try:
        results = executor.process_multiple_genomes(
            args.pattern, args.output,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            report_mutations=not args.skip_mutations,
        )

        executor.logger.info("\n" + "=" * 50)
        executor.logger.info("EcoliTyper AMRfinderPlus FINAL SUMMARY")
        executor.logger.info("=" * 50)

        total_hits = 0
        total_mutations = 0
        high_risk_count = 0
        critical_risk_count = 0
        for genome_name, result in results.items():
            total_hits += result['hit_count']
            total_mutations += result.get('mutation_count', 0)
            genes = [h.get('gene_symbol') for h in result['hits'] if h.get('gene_symbol')]
            high_risk_count += sum(1 for g in genes if g in executor.high_risk_genes)
            critical_risk_count += sum(1 for g in genes if g in executor.critical_risk_genes)
            executor.logger.info("%s: %d AMR hits, %d mutations",
                                 genome_name, result['hit_count'],
                                 result.get('mutation_count', 0))

        executor.logger.info("\nE. COLI SUMMARY STATISTICS:")
        executor.logger.info("   Total genomes processed: %d", len(results))
        executor.logger.info("   Total AMR hits: %d", total_hits)
        executor.logger.info("   Total point mutations: %d", total_mutations)
        executor.logger.info("   High-risk genes detected: %d", high_risk_count)
        executor.logger.info("   Critical risk genes detected: %d", critical_risk_count)
        executor.logger.info("   Average AMR hits per genome: %.1f",
                             total_hits / len(results) if results else 0)

        executor.logger.info("\nPER-GENOME OUTPUTS (2 files each):")
        executor.logger.info("   - {genome}_amrfinder_report.html   (AMR + mutations combined)")
        executor.logger.info("   - {genome}_amrfinder_raw.tsv       (Source column: AMR / MUTATION)")

        executor.logger.info("\nBATCH-LEVEL OUTPUTS:")
        executor.logger.info("   Comprehensive AMR data: %s/ecoli_amrfinder_summary.tsv", args.output)
        executor.logger.info("   Statistics summary: %s/ecoli_amrfinder_statistics_summary.tsv", args.output)
        executor.logger.info("   Master JSON summary: %s/ecoli_amrfinder_master_summary.json", args.output)
        executor.logger.info("   Summary HTML report: %s/ecoli_amrfinder_summary_report.html", args.output)
        executor.logger.info("   Mutation TSV: %s/mutation_summary.tsv", args.output)
        executor.logger.info("   Mutation HTML summary: %s/mutation_summary.html", args.output)
        executor.logger.info("   Mutation master JSON: %s/mutation_master_summary.json", args.output)

        executor.logger.info("\nPERFORMANCE SUMMARY:")
        executor.logger.info("   CPU cores utilized: %d cores", executor.cpus)
        executor.logger.info("   Available RAM: %.1f GB", executor.available_ram)
        executor.logger.info("   Bundled AMRfinderPlus: %s", executor.metadata['amrfinder_version'])
        executor.logger.info("   Bundled database: %s", executor.metadata['database_version'])

        if critical_risk_count > 0:
            executor.logger.info("\nCRITICAL RISK ALERT: last-resort antibiotic resistance genes detected.")
            executor.logger.info("   Immediate clinical attention and infection control measures required.")

        quote = random.choice(executor.science_quotes)
        executor.logger.info("\n%s - %s", quote['text'], quote['author'])

    except Exception as e:
        executor.logger.error("E. coli AMR analysis failed: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()