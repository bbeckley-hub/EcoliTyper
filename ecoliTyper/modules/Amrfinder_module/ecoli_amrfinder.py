#!/usr/bin/env python3
"""
EcoliTyper AMRfinderPlus - E. coli AMR Analysis with Dynamic Database
Comprehensive AMR analysis for E. coli with beautiful HTML reporting - MAXIMUM SPEED VERSION
Author: Beckley Brown <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
Date: 2025 / Updated 2026-06-20
Send a quick mail for any issues or further explanations.
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any, Optional
import argparse
import re
from datetime import datetime
import psutil
import math
import json
import random
from collections import defaultdict

class EcoliAMRfinderPlus:
    """AMRfinderPlus executor for E. coli with DYNAMIC database detection and update capability"""
    
    def __init__(self, cpus: int = None):
        # Setup logging FIRST
        self.logger = self._setup_logging()
        
        # Get module directory and set bundled paths
        self.module_dir = os.path.dirname(os.path.abspath(__file__))
        self.bundled_amrfinder = os.path.join(self.module_dir, "bin", "amrfinder")
        self.bundled_update = os.path.join(self.module_dir, "bin", "amrfinder_update")
        
        # Initialize available_ram before calculating cpus
        self.available_ram = self._get_available_ram()
        
        # Then calculate resources - MAXIMUM SPEED MODE
        self.cpus = self._calculate_optimal_cpus(cpus)
        
        # DYNAMIC DATABASE: find the latest dated folder (starts with 20)
        self.bundled_database = self._get_latest_database()
        
        # If no database found, log warning but do not raise (analysis will be skipped later)
        if self.bundled_database is None:
            self.logger.warning("No AMRfinderPlus database found. Please run with --update-db to download.")
        
        # Read database version or set to Unknown
        db_version = self._get_database_version() if self.bundled_database else "Unknown"
        
        self.metadata = {
            "tool_name": "EcoliTyper AMRfinderPlus",
            "version": "1.2.1",   
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "amrfinder_version": "4.2.7",
            "database_version": db_version
        }
        
        # Comprehensive high-risk and critical gene sets
        self.high_risk_genes = {
            # Critical Beta-lactamases
            'blaCTX-M-14', 'blaCTX-M-1', 'blaTEM-1', 'blaEC', 'blaCTX-M', 'blaTEM', 'blaSHV',
            
            # Colistin resistance (CRITICAL - last resort antibiotic)
            'mcr-1.1', 'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', '(Col)mcr-1.1', 'MCR',
            
            # Carbapenemases (HIGHEST RISK)
            'blaKPC', 'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP',
            
            # Fluoroquinolone resistance
            'qnrA', 'qnrB', 'qnrC', 'qnrD', 'qnrS', 'qnrVC',
            
            # Aminoglycoside resistance
            'aac(3)-IId', 'aac(6\')-Ib-cr', 'aadA1', 'aadA2', 
            'aph(3\'\')-Ib', 'aph(3\')-Ia', 'aph(6)-Id',
            
            # Tetracycline resistance
            'tet(A)', 'tet(C)',
            
            # Sulfonamide resistance
            'sul1', 'sul2', 'sul3',
            
            # Trimethoprim resistance
            'dfrA1', 'dfrA14',
            
            # Chloramphenicol resistance
            'catA2', 'cmlA', 'floR',
            
            # Macrolide resistance
            'ermA', 'ermB', 'ermC', 'mphA',
            
            # Multi-drug efflux pumps
            'acrF', 'emrD', 'emrE', 'mdtM',
            
            # Other high-risk markers
            'armA', 'rmtB', 'cfr', 'optrA', 'poxtA', 'CTX'
        }

        # CRITICAL RISK genes - highest priority
        self.critical_risk_genes = {
            'mcr-1.1', 'mcr-1', 'blaCTX-M-14', 'blaCTX-M-1', 'blaKPC', 
            'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP', 'cfr'
        }
        
        self.science_quotes = [
            "The important thing is not to stop questioning. Curiosity has its own reason for existence. - Albert Einstein",
            "Nothing in life is to be feared, it is only to be understood. - Marie Curie", 
            "The microscope opens a new world to the investigator. - Robert Koch",
            "In science, the credit goes to the man who convinces the world, not to the man to whom the idea first occurs. - Francis Darwin",
            "The good thing about science is that it's true whether or not you believe in it. - Neil deGrasse Tyson",
            "Science knows no country, because knowledge belongs to humanity. - Louis Pasteur"
        ]
    
    def _get_ascii_art(self) -> str:
        """Return the ECOLITYPER ASCII art banner"""
        return """
███████╗ ██████╗ ██████╗ ██╗     ██╗████████╗██╗   ██╗██████╗ ███████╗██████╗ 
██╔════╝██╔════╝██╔═══██╗██║     ██║╚══██╔══╝╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
█████╗  ██║     ██║   ██║██║     ██║   ██║    ╚████╔╝ ██████╔╝█████╗  ██████╔╝
██╔══╝  ██║     ██║   ██║██║     ██║   ██║     ╚██╔╝  ██╔═══╝ ██╔══╝  ██╔══██╗
███████╗╚██████╗╚██████╔╝███████╗██║   ██║      ██║   ██║     ███████╗██║  ██║
╚══════╝ ╚═════╝ ╚═════╝ ╚══════╝╚═╝   ╚═╝      ╚═╝   ╚═╝     ╚══════╝╚═╝  ╚═╝
        """
    
    def _setup_logging(self):
        """Setup logging - must be called first in __init__"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)
    
    def _get_available_ram(self) -> int:
        """Get available RAM in GB"""
        try:
            ram_gb = psutil.virtual_memory().available / (1024 ** 3)
            return ram_gb
        except Exception as e:
            self.logger.warning(f"Could not detect RAM: {e}")
            return 8  # Assume 8GB as fallback
    
    def _calculate_optimal_cpus(self, user_cpus: int = None) -> int:
        """Calculate optimal number of CPU cores for MAXIMUM SPEED"""
        if user_cpus is not None:
            self._log_resource_info(user_cpus)
            return user_cpus
            
        try:
            # Get total PHYSICAL CPU cores (not logical threads)
            total_physical_cores = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            
            # MAXIMUM SPEED RULES - AGGRESSIVE CPU USAGE 
            if total_physical_cores <= 4:
                optimal_cpus = total_physical_cores  # Use ALL cores on small systems
            elif total_physical_cores <= 8:
                optimal_cpus = total_physical_cores - 1  # Use 7/8, 6/7, etc.
            elif total_physical_cores <= 16:
                optimal_cpus = max(8, total_physical_cores - 1)  # Use 15/16, 14/15, etc.
            elif total_physical_cores <= 32:
                optimal_cpus = max(16, total_physical_cores - 1)  # Use 31/32, 30/31, etc.
            else:
                optimal_cpus = min(32, int(total_physical_cores * 0.95))  # Use 95% on huge systems
            
            # Ensure at least 1 CPU and not more than available cores
            optimal_cpus = max(1, min(optimal_cpus, total_physical_cores))
            
            self._log_resource_info(optimal_cpus, total_physical_cores)
            return optimal_cpus
            
        except Exception as e:
            # Fallback to using all available cores for maximum speed
            self.logger.warning(f"Could not detect CPU cores, using maximum available: {e}")
            return os.cpu_count() or 4
    
    def _log_resource_info(self, cpus: int, total_cores: int = None):
        """Log resource allocation information - KEEPING EcoliTyper STYLING"""
        self.logger.info(f"Available RAM: {self.available_ram:.1f} GB")
        
        if total_cores:
            self.logger.info(f"System CPU cores: {total_cores}")
            utilization = (cpus / total_cores) * 100
            self.logger.info(f"Using CPU cores: {cpus} ({utilization:.1f}% of available cores)")
        else:
            self.logger.info(f"Using user-specified CPU cores: {cpus}")
        
        # Performance recommendations - MAXIMUM SPEED FOCUS 
        if cpus == 1:
            self.logger.info("💡 Performance: Single-core (max speed for 1-core systems)")
        elif cpus <= 4:
            self.logger.info("💡 Performance: Multi-core (max speed for small systems)")
        elif cpus <= 8:
            self.logger.info("💡 Performance: High-speed multi-core mode")
        elif cpus <= 16:
            self.logger.info("💡 Performance: Ultra-speed multi-core mode 🚀")
        elif cpus <= 32:
            self.logger.info("💡 Performance: MAXIMUM SPEED MULTI-CORE MODE 🚀🔥")
        else:
            self.logger.info("💡 Performance: EXTREME SPEED MULTI-CORE MODE 🚀🔥💨")
        
        # Strategy note - UPDATED for concurrent processing
        self.logger.info("📝 STRATEGY: Processing MULTIPLE samples concurrently with optimal core allocation for maximum throughput")
    
    def _get_latest_database(self) -> Optional[str]:
        """Find the latest dated database folder in data/amrfinder_db/ (starts with 20)"""
        db_root = os.path.join(self.module_dir, "data", "amrfinder_db")
        if not os.path.exists(db_root):
            self.logger.warning(f"Database root directory not found: {db_root}")
            return None
        # Find all subdirectories starting with '20'
        candidates = []
        for item in os.listdir(db_root):
            full_path = os.path.join(db_root, item)
            if os.path.isdir(full_path) and item.startswith('20'):
                candidates.append(item)
        if not candidates:
            self.logger.warning("No database folder starting with '20' found.")
            return None
        # Sort lexicographically (YYYY-MM-DD works) and take the latest
        latest = sorted(candidates)[-1]
        latest_path = os.path.join(db_root, latest)
        self.logger.info(f"Using latest database: {latest_path}")
        return latest_path
    
    def _get_database_version(self) -> str:
        """Read version.txt from the database folder or fallback to folder name"""
        if not self.bundled_database:
            return "Unknown"
        version_file = os.path.join(self.bundled_database, "version.txt")
        if os.path.exists(version_file):
            with open(version_file, 'r') as f:
                return f.read().strip()
        # Fallback to folder name
        return os.path.basename(self.bundled_database)
    
    def update_database(self) -> bool:
        """Download the latest AMRfinderPlus database using bundled amrfinder_update"""
        if not os.path.exists(self.bundled_update):
            self.logger.error(f"amrfinder_update not found at {self.bundled_update}")
            return False
        if not os.access(self.bundled_update, os.X_OK):
            self.logger.warning("amrfinder_update not executable, fixing permissions...")
            os.chmod(self.bundled_update, 0o755)
        db_dir = os.path.join(self.module_dir, "data", "amrfinder_db")
        os.makedirs(db_dir, exist_ok=True)
        self.logger.info("Updating AMRfinderPlus database...")
        try:
            cmd = [self.bundled_update, "--database", db_dir]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.logger.info("Database update completed successfully.")
            # Re‑detect latest database
            self.bundled_database = self._get_latest_database()
            if self.bundled_database:
                self.metadata['database_version'] = self._get_database_version()
                self.logger.info(f"New database version: {self.metadata['database_version']}")
                return True
            else:
                self.logger.error("Database update succeeded but no database folder found.")
                return False
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Database update failed: {e}")
            self.logger.error(f"STDERR: {e.stderr}")
            return False

    def check_amrfinder_installed(self) -> bool:
        """Check if bundled AMRfinderPlus is available and database exists"""
        try:
            if not os.path.exists(self.bundled_amrfinder):
                self.logger.error(f"Bundled AMRfinderPlus not found at: {self.bundled_amrfinder}")
                return False
            
            if not os.access(self.bundled_amrfinder, os.X_OK):
                self.logger.warning(f"Bundled AMRfinderPlus not executable, fixing permissions...")
                os.chmod(self.bundled_amrfinder, 0o755)
            
            # Test the bundled version
            result = subprocess.run(
                [self.bundled_amrfinder, '--version'], 
                capture_output=True, 
                text=True, 
                check=True
            )
            version_line = result.stdout.strip()
            self.logger.info(f"Bundled AMRfinderPlus version: {version_line}")
            
            # Check database
            if self.bundled_database and os.path.exists(self.bundled_database):
                self.logger.info(f"✅ Bundled database found: {self.bundled_database}")
                db_version_file = os.path.join(self.bundled_database, "version.txt")
                if os.path.exists(db_version_file):
                    with open(db_version_file, 'r') as f:
                        db_version = f.read().strip()
                        self.logger.info(f"✅ Database version: {db_version}")
                else:
                    self.logger.info(f"✅ Database folder: {os.path.basename(self.bundled_database)}")
                return True
            else:
                self.logger.warning(f"⚠️ Bundled database not found at expected location.")
                self.logger.info("Please run with --update-db to download the latest database.")
                return False
            
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.logger.error(f"Bundled AMRfinderPlus check failed: {e}")
            return False

    def run_amrfinder_single_genome(self, genome_file: str, output_dir: str,
                                    min_identity: float = None, min_coverage: float = None,
                                    report_mutations: bool = True) -> Dict[str, Any]:
        """Run AMRfinderPlus on a single E. coli genome - with optional mutation reporting and thresholds"""
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"{genome_name}_amrfinder.txt")
        
        # Check if bundled binary exists
        if not os.path.exists(self.bundled_amrfinder):
            self.logger.error(f"Bundled AMRfinderPlus not found at: {self.bundled_amrfinder}")
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'mutations_file': None,
                'status': 'failed',
                'error': 'AMRfinder binary not found'
            }
        
        # AMRfinderPlus uses THREADS - allocate ALL available cores for maximum speed
        run_threads = self.cpus
        
        # Build command with BUNDLED resources
        cmd = [
            self.bundled_amrfinder,
            '-n', genome_file,  # Nucleotide mode
            '-O', 'Escherichia',  # Organism (E. coli)
            '--output', output_file,
            '--plus'
        ]
        
        # Add dynamic database if available
        if self.bundled_database and os.path.exists(self.bundled_database):
            cmd.extend(['--database', self.bundled_database])
            self.logger.info(f"Using bundled database: {self.bundled_database}")
        else:
            self.logger.warning("Using default AMRfinderPlus database location")
        
        # Add min identity and min coverage if provided
        if min_identity is not None:
            cmd.extend(['--ident_min', str(min_identity)])
            self.logger.info(f"Using minimum identity: {min_identity}")
        if min_coverage is not None:
            cmd.extend(['--coverage_min', str(min_coverage)])
            self.logger.info(f"Using minimum coverage: {min_coverage}")
        
        # Mutation output file
        mut_file = None
        if report_mutations:
            mut_file = os.path.join(output_dir, f"{genome_name}_mutations.tsv")
            cmd.extend(['--mutation_all', mut_file])
            self.logger.info(f"Will report point mutations to {mut_file}")
        
        self.logger.info("🚀 MAXIMUM SPEED: Running AMRfinderPlus on %s (using ALL %d CORES)", genome_name, run_threads)
        self.logger.debug(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse results for reporting
            hits = self._parse_amrfinder_output(output_file)
            
            # Create individual HTML report
            self._create_amrfinder_html_report(genome_name, hits, output_dir)
            
            # Create individual JSON report
            self._create_amrfinder_json_report(genome_name, hits, output_dir)
            
            # If mutations were requested and file exists, create mutation HTML report
            if mut_file and os.path.exists(mut_file):
                self._create_mutation_html_report(genome_name, mut_file, output_dir)
            
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': hits,
                'hit_count': len(hits),
                'mutations_file': mut_file,
                'status': 'success'
            }
            
        except subprocess.CalledProcessError as e:
            self.logger.error("AMRfinderPlus failed for %s: %s", genome_name, e.stderr)
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'mutations_file': None,
                'status': 'failed'
            }
    
    def _parse_amrfinder_output(self, amrfinder_file: str) -> List[Dict]:
        """Parse AMRfinderPlus 4.2.4 output file into structured data"""
        hits = []
        try:
            with open(amrfinder_file, 'r') as f:
                lines = f.readlines()
                
            if not lines or len(lines) < 2:
                return hits
                
            # Parse header - AMRfinderPlus 4.2.4 uses new headers
            headers = lines[0].strip().split('\t')
            
            # Parse data lines
            for line_num, line in enumerate(lines[1:], 2):
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split('\t')
                if len(parts) >= len(headers):
                    # Create hit with original headers
                    hit = {}
                    for i, header in enumerate(headers):
                        if i < len(parts):
                            hit[header] = parts[i]
                        else:
                            hit[header] = ''
                    
                    # Map to consistent field names - BOTH OLD AND NEW HEADERS
                    processed_hit = {
                        'Protein id': hit.get('Protein id', ''),
                        'Contig id': hit.get('Contig id', ''),
                        'Start': hit.get('Start', ''),
                        'Stop': hit.get('Stop', ''),
                        'Strand': hit.get('Strand', ''),
                        'Element symbol': hit.get('Element symbol', ''),
                        'Element name': hit.get('Element name', ''),
                        'Scope': hit.get('Scope', ''),
                        'Type': hit.get('Type', ''),
                        'Subtype': hit.get('Subtype', ''),
                        'Class': hit.get('Class', ''),
                        'Subclass': hit.get('Subclass', ''),
                        'Method': hit.get('Method', ''),
                        'Target length': hit.get('Target length', ''),
                        'Reference sequence length': hit.get('Reference sequence length', ''),
                        '% Coverage of reference': hit.get('% Coverage of reference', ''),
                        '% Identity to reference': hit.get('% Identity to reference', ''),
                        'Alignment length': hit.get('Alignment length', ''),
                        'Closest reference accession': hit.get('Closest reference accession', ''),
                        'Closest reference name': hit.get('Closest reference name', ''),
                        'HMM accession': hit.get('HMM accession', ''),
                        'HMM description': hit.get('HMM description', ''),
                        
                        # Old header compatibility
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
                        
                        '_original_headers': headers,
                        '_original_values': parts
                    }
                    hits.append(processed_hit)
                else:
                    self.logger.warning("Line %d has %d parts, expected %d: %s", 
                                      line_num, len(parts), len(headers), line[:100] + "...")
                    
        except Exception as e:
            self.logger.error("Error parsing %s: %s", amrfinder_file, e)
            
        self.logger.info("Parsed %d AMR hits from %s", len(hits), amrfinder_file)
        return hits
    
    def _parse_mutations_file(self, mut_file: str) -> List[Dict]:
        """Parse mutations TSV file (same format as AMRfinder output)"""
        return self._parse_amrfinder_output(mut_file)
    
    def _create_mutation_html_report(self, genome_name: str, mutations_file: str, output_dir: str):
        """Create a beautiful HTML report for point mutations"""
        mutations = self._parse_mutations_file(mutations_file)
        if not mutations:
            self.logger.info(f"No mutations found for {genome_name}, skipping mutation HTML.")
            return
        
        random_quote = random.choice(self.science_quotes)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        html = f"""<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>EcoliTyper - Mutation Report: {genome_name}</title>
        <style>
            * {{ margin: 0; padding: 0; box-sizing: border-box; }}
            body {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                color: #ffffff;
                padding: 20px;
                min-height: 100vh;
            }}
            .container {{ max-width: 1400px; margin: 0 auto; }}
            .header {{ text-align: center; margin-bottom: 30px; }}
            .ascii-container {{
                background: rgba(0, 0, 0, 0.7);
                padding: 20px;
                border-radius: 15px;
                margin-bottom: 20px;
                box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
                border: 2px solid rgba(0, 255, 0, 0.3);
            }}
            .ascii-art {{
                font-family: 'Courier New', monospace;
                font-size: 10px;
                line-height: 1.1;
                white-space: pre;
                color: #00ff00;
                text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
                overflow-x: auto;
            }}
            .quote-container {{
                background: rgba(255, 255, 255, 0.1);
                backdrop-filter: blur(10px);
                padding: 20px;
                border-radius: 10px;
                margin-bottom: 30px;
                text-align: center;
                box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
                border: 1px solid rgba(255, 255, 255, 0.2);
            }}
            .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; }}
            .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
            .report-section {{
                background: rgba(255, 255, 255, 0.95);
                color: #1f2937;
                padding: 25px;
                border-radius: 10px;
                margin-bottom: 20px;
                box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            }}
            .report-section h2 {{
                color: #1e3a8a;
                border-bottom: 3px solid #3b82f6;
                padding-bottom: 10px;
                margin-bottom: 20px;
                font-size: 24px;
            }}
            .table-responsive {{ width: 100%; overflow-x: auto; margin: 20px 0; }}
            .mutation-table {{
                width: 100%;
                border-collapse: collapse;
                font-size: 13px;
                min-width: 1000px;
            }}
            .mutation-table th {{
                background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
                color: white;
                padding: 12px;
                text-align: left;
            }}
            .mutation-table td {{
                padding: 10px;
                border-bottom: 1px solid #e5e7eb;
            }}
            .mutation-table tr:nth-child(even) {{ background-color: #f8fafc; }}
            .footer {{
                text-align: center;
                margin-top: 30px;
                padding: 20px;
                background: rgba(0, 0, 0, 0.3);
                border-radius: 10px;
                font-size: 14px;
            }}
            .timestamp {{ color: #fbbf24; font-weight: bold; }}
            .authorship {{ margin-top: 15px; padding: 15px; background: rgba(255, 255, 255, 0.1); border-radius: 8px; font-size: 12px; }}
        </style>
    </head>
    <body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self._get_ascii_art()}</div>
            </div>
            <div class="quote-container">
                <div class="quote-text">"{random_quote}"</div>
            </div>
        </div>

        <div class="report-section">
            <h2>🧬 Point Mutation Report: {genome_name}</h2>
            <p>All point mutations detected by AMRfinderPlus (including synonymous variants).</p>
            <div class="table-responsive">
                <table class="mutation-table">
                    <thead>
                        <tr><th>Gene Symbol</th><th>Mutation</th><th>Class</th><th>Subclass</th>
                        <th>Contig</th><th>Start</th><th>Stop</th><th>Strand</th>
                        <th>Coverage (%)</th><th>Identity (%)</th><th>Accession</th></tr>
                    </thead>
                    <tbody>
    """
        for m in mutations:
            html += f"""
                        <tr>
                            <td>{m.get('gene_symbol', '')}</td>
                            <td>{m.get('element_name', '')}</td>
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
        html += f"""
                    </tbody>
                </table>
            </div>
        </div>

        <div class="footer">
            <p><strong>ECOLITYPER</strong> - Mutation Analysis Module</p>
            <p class="timestamp">Generated: {current_time}</p>
            <div class="authorship">
                <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
        </div>
    </div>
    </body>
    </html>"""
        out_file = os.path.join(output_dir, f"{genome_name}_mutations.html")
        with open(out_file, 'w') as f:
            f.write(html)
        self.logger.info(f"✓ Mutation HTML report: {out_file}")
    
    def _create_amrfinder_html_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        """Create comprehensive HTML report for AMRfinderPlus results with beautiful styling and ASCII art"""
        
        # Analyze AMR results for E. coli
        analysis = self._analyze_ecoli_amr_results(hits)
        
        # JavaScript for rotating quotes
        quotes_js = f"""
        <script>
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            
            // Rotate every 10 seconds
            setInterval(rotateQuote, 10000);
            
            // Initial display
            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
            }});
        </script>
        """
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>EcoliTyper AMRfinderPlus Analysis Report</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        body {{ 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
            margin: 0; 
            padding: 0; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}
        .container {{ 
            max-width: 1200px; 
            margin: 0 auto; 
            padding: 20px; 
        }}
        .header {{ 
            background: rgba(255, 255, 255, 0.95); 
            padding: 30px; 
            border-radius: 15px; 
            margin-bottom: 30px; 
            box-shadow: 0 8px 32px rgba(0,0,0,0.1);
            backdrop-filter: blur(10px);
        }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            text-align: center;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 12px;
            line-height: 1.2;
            white-space: pre;
            color: #00ff00;
            text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
            overflow-x: auto;
        }}
        .card {{ 
            background: rgba(255, 255, 255, 0.95); 
            padding: 25px; 
            margin: 20px 0; 
            border-radius: 12px; 
            box-shadow: 0 4px 20px rgba(0,0,0,0.1);
            backdrop-filter: blur(10px);
        }}
        .gene-table, .class-table {{ 
            width: 100%; 
            border-collapse: collapse; 
            margin: 20px 0; 
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .gene-table th, .gene-table td, .class-table th, .class-table td {{ 
            padding: 15px; 
            text-align: left; 
            border-bottom: 1px solid #e0e0e0; 
        }}
        .gene-table th, .class-table th {{ 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            font-weight: 600;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .success {{ color: #28a745; font-weight: 600; }}
        .warning {{ color: #ffc107; font-weight: 600; }}
        .error {{ color: #dc3545; font-weight: 600; }}
        .summary-stats {{ 
            display: flex; 
            justify-content: space-around; 
            margin: 20px 0; 
            flex-wrap: wrap;
        }}
        .stat-card {{ 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px; 
            border-radius: 12px; 
            text-align: center; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .critical-stat-card {{
            background: linear-gradient(135deg, #dc3545 0%, #c82333 100%);
            color: white;
            padding: 20px; 
            border-radius: 12px; 
            text-align: center; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            color: white;
            padding: 20px;
            border-radius: 12px;
            margin: 20px 0;
            text-align: center;
            font-style: italic;
            border-left: 4px solid #fff;
        }}
        .footer {{
            background: rgba(0, 0, 0, 0.8);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin-top: 40px;
        }}
        .footer a {{
            color: #667eea;
            text-decoration: none;
        }}
        .footer a:hover {{
            text-decoration: underline;
        }}
        .resistance-badge {{
            display: inline-block;
            background: #667eea;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .high-risk {{ background: #dc3545; }}
        .critical-risk {{ background: #8b0000; font-weight: bold; }}
        .medium-risk {{ background: #ffc107; color: black; }}
        .low-risk {{ background: #28a745; }}
        .present {{ background-color: #d4edda; }}
        .critical-row {{ background-color: #f8d7da; font-weight: bold; border-left: 4px solid #dc3545; }}
        .high-risk-row {{ background-color: #fff3cd; border-left: 4px solid #ffc107; }}
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">
{self._get_ascii_art()}
                </div>
            </div>
            <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 EcoliTyper AMRfinderPlus Analysis Report</h1>
            <p style="color: #666; font-size: 1.2em;">Comprehensive E. coli Antimicrobial Resistance Analysis</p>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
"""
        
        # CRITICAL RISK ALERT - Show first if critical genes detected
        if analysis['critical_risk_genes'] > 0:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #dc3545; background: #f8d7da;">
            <h2 style="color: #dc3545;">🚨 CRITICAL RISK AMR GENES DETECTED</h2>
            <p><strong>{analysis['critical_risk_genes']} CRITICAL RISK antimicrobial resistance genes found:</strong></p>
            <div style="margin: 10px 0;">
                <p style="color: #721c24; font-weight: bold;">
                    ⚠️ These genes confer resistance to last-resort antibiotics and represent 
                    a serious public health concern requiring immediate attention.
                </p>
"""
            for gene in analysis['critical_risk_list']:
                html_content += f'<span class="resistance-badge critical-risk" style="font-size: 1.1em;">🚨 {gene}</span>'
            html_content += """
            </div>
        </div>
"""
        
        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">📊 E. coli AMR Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total AMR Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['total_genes']}</p>
                </div>
                <div class="stat-card">
                    <h3>High Risk Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['high_risk_genes']}</p>
                </div>
                <div class="critical-stat-card">
                    <h3>Critical Risk</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['critical_risk_genes']}</p>
                </div>
            </div>
            <p><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>AMRfinderPlus Version:</strong> {self.metadata['amrfinder_version']}</p>
            <p><strong>Database Version:</strong> {self.metadata['database_version']}</p>
        </div>
"""
        
        # High-risk genes warning (non-critical)
        if analysis['high_risk_genes'] > 0 and analysis['critical_risk_genes'] == 0:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #ffc107;">
            <h2 style="color: #856404;">⚠️ High-Risk AMR Genes Detected</h2>
            <p><strong>{analysis['high_risk_genes']} high-risk antimicrobial resistance genes found:</strong></p>
            <div style="margin: 10px 0;">
"""
            for gene in analysis['high_risk_list']:
                html_content += f'<span class="resistance-badge high-risk">{gene}</span>'
            html_content += """
            </div>
        </div>
"""
        
        # Resistance Mechanism Breakdown
        if any(analysis['resistance_mechanisms'].values()):
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🔬 Resistance Mechanism Breakdown</h2>
"""
            
            mechanisms = analysis['resistance_mechanisms']
            if mechanisms['esbl']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #fff3cd; border-radius: 5px;">
                <strong>ESBL Genes:</strong> {', '.join(mechanisms['esbl'])}
            </div>
"""
            if mechanisms['carbapenemase']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #f8d7da; border-radius: 5px;">
                <strong>Carbapenemase Genes (CRITICAL):</strong> {', '.join(mechanisms['carbapenemase'])}
            </div>
"""
            if mechanisms['colistin_resistance']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #f8d7da; border-radius: 5px;">
                <strong>Colistin Resistance (CRITICAL):</strong> {', '.join(mechanisms['colistin_resistance'])}
            </div>
"""
            if mechanisms['fluoroquinolone_resistance']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #d1ecf1; border-radius: 5px;">
                <strong>Fluoroquinolone Resistance:</strong> {', '.join(mechanisms['fluoroquinolone_resistance'])}
            </div>
"""
            if mechanisms['aminoglycoside_resistance']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #d1ecf1; border-radius: 5px;">
                <strong>Aminoglycoside Resistance:</strong> {', '.join(mechanisms['aminoglycoside_resistance'])}
            </div>
"""
            if mechanisms['efflux_pumps']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #e2e3e5; border-radius: 5px;">
                <strong>Efflux Pumps:</strong> {', '.join(mechanisms['efflux_pumps'])}
            </div>
"""
            
            html_content += """
        </div>
"""
        
        # Resistance classes summary
        if analysis['resistance_classes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🧪 Resistance Classes Detected</h2>
            <table class="class-table">
                <thead>
                    <tr>
                        <th>Resistance Class</th>
                        <th>Gene Count</th>
                        <th>Genes</th>
                    </tr>
                </thead>
                <tbody>
"""
            
            for class_name, genes in analysis['resistance_classes'].items():
                gene_list = ", ".join(genes)
                html_content += f"""
                    <tr>
                        <td><strong>{class_name}</strong></td>
                        <td>{len(genes)}</td>
                        <td>{gene_list}</td>
                    </tr>
"""
            
            html_content += """
                </tbody>
            </table>
        </div>
"""
        
        # Detailed AMR genes table
        if hits:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🔬 Detailed AMR Genes Detected</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Gene Symbol</th>
                        <th>Sequence Name</th>
                        <th>Class</th>
                        <th>Subclass</th>
                        <th>Coverage</th>
                        <th>Identity</th>
                        <th>Scope</th>
                    </tr>
                </thead>
                <tbody>
"""
            
            for hit in hits:
                # Determine row class based on risk level
                row_class = "present"
                gene_symbol = hit.get('gene_symbol', '')
                if gene_symbol in analysis['critical_risk_list']:
                    row_class = "critical-row"
                elif gene_symbol in analysis['high_risk_list']:
                    row_class = "high-risk-row"
                
                html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{gene_symbol}</strong></td>
                        <td title="{hit.get('sequence_name', '')}">{hit.get('sequence_name', '')[:1000]}{'...' if len(hit.get('sequence_name', '')) > 1000 else ''}</td>
                        <td>{hit.get('class', '')}</td>
                        <td>{hit.get('subclass', '')}</td>
                        <td>{hit.get('coverage', '')}%</td>
                        <td>{hit.get('identity', '')}%</td>
                        <td>{hit.get('scope', '')}</td>
                    </tr>
"""
            
            html_content += """
                </tbody>
            </table>
        </div>
"""
        else:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">✅ No AMR Genes Detected</h2>
            <p>No antimicrobial resistance genes found in this E. coli genome.</p>
        </div>
"""
        
        # Footer
        html_content += f"""
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #667eea; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
            </p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write HTML report
        html_file = os.path.join(output_dir, f"{genome_name}_amrfinder_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info("E. coli AMRfinderPlus HTML report generated: %s", html_file)
    
    def _create_amrfinder_json_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        """Create JSON report for AMRfinderPlus results"""
        analysis = self._analyze_ecoli_amr_results(hits)
        
        json_data = {
            'metadata': {
                'genome': genome_name,
                'analysis_date': self.metadata['analysis_date'],
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'amrfinder_version': self.metadata['amrfinder_version'],
                'database_version': self.metadata['database_version']
            },
            'summary': {
                'total_genes': analysis['total_genes'],
                'high_risk_genes': analysis['high_risk_genes'],
                'critical_risk_genes': analysis['critical_risk_genes'],
                'high_risk_list': analysis['high_risk_list'],
                'critical_risk_list': analysis['critical_risk_list'],
                'resistance_classes': analysis['resistance_classes'],
                'resistance_mechanisms': analysis['resistance_mechanisms']
            },
            'hits': hits
        }
        
        json_file = os.path.join(output_dir, f"{genome_name}_amrfinder_report.json")
        with open(json_file, 'w') as f:
            json.dump(json_data, f, indent=2)
        
        self.logger.info("E. coli AMRfinderPlus JSON report generated: %s", json_file)
    
    def _analyze_ecoli_amr_results(self, hits: List[Dict]) -> Dict[str, Any]:
        """Analyze AMR results specifically for E. coli with enhanced risk assessment"""
        
        analysis = {
            'total_genes': len(hits),
            'resistance_classes': {},
            'total_classes': 0,
            'high_risk_genes': 0,
            'critical_risk_genes': 0,
            'high_risk_list': [],
            'critical_risk_list': [],
            'resistance_mechanisms': {
                'esbl': [],          # Extended Spectrum Beta-lactamases
                'carbapenemase': [], # Carbapenemases
                'colistin_resistance': [], # Colistin resistance
                'fluoroquinolone_resistance': [], # Fluoroquinolone resistance
                'aminoglycoside_resistance': [], # Aminoglycoside resistance
                'efflux_pumps': [],  # Multi-drug efflux pumps
                'other_amr': []      # Other resistance mechanisms
            }
        }
        
        for hit in hits:
            gene_symbol = hit.get('gene_symbol', '')
            resistance_class = hit.get('class', '')
            
            # Categorize resistance mechanism
            self._categorize_resistance_mechanism(gene_symbol, resistance_class, analysis)
            
            # Check for critical risk genes
            if gene_symbol in self.critical_risk_genes:
                analysis['critical_risk_genes'] += 1
                if gene_symbol not in analysis['critical_risk_list']:
                    analysis['critical_risk_list'].append(gene_symbol)
            
            # Check for high-risk genes (includes critical ones)
            if gene_symbol in self.high_risk_genes:
                analysis['high_risk_genes'] += 1
                if gene_symbol not in analysis['high_risk_list']:
                    analysis['high_risk_list'].append(gene_symbol)
            
            # Group by resistance class
            if resistance_class:
                if resistance_class not in analysis['resistance_classes']:
                    analysis['resistance_classes'][resistance_class] = []
                if gene_symbol not in analysis['resistance_classes'][resistance_class]:
                    analysis['resistance_classes'][resistance_class].append(gene_symbol)
        
        analysis['total_classes'] = len(analysis['resistance_classes'])
        return analysis

    def _categorize_resistance_mechanism(self, gene_symbol: str, resistance_class: str, analysis: Dict[str, Any]):
        """Categorize genes by resistance mechanism"""
        
        # ESBL genes
        esbl_genes = {'blaCTX-M', 'blaTEM', 'blaSHV', 'blaCTX-M-14', 'blaCTX-M-1', 'blaTEM-1'}
        
        # Carbapenemase genes
        carbapenemase_genes = {'blaKPC', 'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP'}
        
        # Colistin resistance genes
        colistin_genes = {'mcr-1.1', 'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', '(Col)mcr-1.1'}
        
        # Fluoroquinolone resistance genes
        fluoroquinolone_genes = {'qnrA', 'qnrB', 'qnrC', 'qnrD', 'qnrS', 'qnrVC'}
        
        # Aminoglycoside resistance genes
        aminoglycoside_genes = {
            'aac(3)-IId', 'aac(6\')-Ib-cr', 'aadA1', 'aadA2', 
            'aph(3\'\')-Ib', 'aph(3\')-Ia', 'aph(6)-Id'
        }
        
        # Efflux pumps
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
    
    def create_amr_summary(self, all_results: Dict[str, Any], output_base: str):
        """Create comprehensive AMR summary files and HTML reports for all E. coli samples"""
        self.logger.info("Creating E. coli AMR summary files and HTML reports...")
        
        # Create TSV summary files
        summary_file = os.path.join(output_base, "ecoli_amrfinder_summary.tsv")
        
        with open(summary_file, 'w') as f:
            # Write header with NEW headers
            f.write("Genome\tProtein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description\n")
            
            # Write data for all genomes
            for genome_name, result in all_results.items():
                for hit in result['hits']:
                    row = [
                        genome_name,
                        hit.get('Protein id', ''),
                        hit.get('Contig id', ''),
                        hit.get('Start', ''),
                        hit.get('Stop', ''),
                        hit.get('Strand', ''),
                        hit.get('Element symbol', ''),
                        hit.get('Element name', ''),
                        hit.get('Scope', ''),
                        hit.get('Type', ''),
                        hit.get('Subtype', ''),
                        hit.get('Class', ''),
                        hit.get('Subclass', ''),
                        hit.get('Method', ''),
                        hit.get('Target length', ''),
                        hit.get('Reference sequence length', ''),
                        hit.get('% Coverage of reference', ''),
                        hit.get('% Identity to reference', ''),
                        hit.get('Alignment length', ''),
                        hit.get('Closest reference accession', ''),
                        hit.get('Closest reference name', ''),
                        hit.get('HMM accession', ''),
                        hit.get('HMM description', '')
                    ]
                    f.write('\t'.join(str(x) for x in row) + '\n')
        
        self.logger.info("✓ E. coli AMR summary file created: %s", summary_file)
        
        # Create statistics summary
        stats_file = os.path.join(output_base, "ecoli_amrfinder_statistics_summary.tsv")
        with open(stats_file, 'w') as f:
            f.write("Genome\tTotal_AMR_Genes\tHigh_Risk_Genes\tCritical_Risk_Genes\tResistance_Classes\tGene_List\n")
            
            for genome_name, result in all_results.items():
                # Get unique genes
                genes = list(set(hit.get('gene_symbol', '') for hit in result['hits'] if hit.get('gene_symbol')))
                gene_list = ",".join(genes)
                
                # Count high-risk and critical genes
                high_risk_count = sum(1 for gene in genes if gene in self.high_risk_genes)
                critical_risk_count = sum(1 for gene in genes if gene in self.critical_risk_genes)
                
                # Get resistance classes
                classes = list(set(hit.get('class', '') for hit in result['hits'] if hit.get('class')))
                class_list = ",".join(classes)
                
                f.write(f"{genome_name}\t{result['hit_count']}\t{high_risk_count}\t{critical_risk_count}\t{class_list}\t{gene_list}\n")
        
        self.logger.info("✓ E. coli AMR statistics summary created: %s", stats_file)
        
        # Create JSON summaries
        self.create_json_summaries(all_results, output_base)
        
        # Create comprehensive HTML summary report
        self._create_summary_html_report(all_results, output_base)
        
        # Create mutation batch summary if any mutations exist
        self.create_mutation_summary(all_results, output_base)
    
    def create_mutation_summary(self, all_results: Dict[str, Any], output_base: str):
        """Create mutation batch summary across all genomes (TSV, HTML, JSON)"""
        self.logger.info("Creating mutation batch summaries...")
        all_mutations = []
        genome_mutation_counts = {}
        
        for genome_name, result in all_results.items():
            if 'mutations_file' in result and result['mutations_file'] and os.path.exists(result['mutations_file']):
                muts = self._parse_mutations_file(result['mutations_file'])
                if muts:
                    genome_mutation_counts[genome_name] = len(muts)
                    for m in muts:
                        m_copy = m.copy()
                        m_copy['genome'] = genome_name
                        all_mutations.append(m_copy)
                else:
                    genome_mutation_counts[genome_name] = 0
            else:
                genome_mutation_counts[genome_name] = 0
        
        if not all_mutations:
            self.logger.info("No mutations found in any genome; skipping mutation summaries.")
            return
        
        # TSV summary
        tsv_file = os.path.join(output_base, "mutation_summary.tsv")
        with open(tsv_file, 'w') as f:
            fieldnames = ['genome', 'gene_symbol', 'element_name', 'class', 'subclass',
                          'contig_id', 'start', 'stop', 'strand', 'coverage', 'identity', 'accession']
            f.write('\t'.join(fieldnames) + '\n')
            for m in all_mutations:
                row = [m.get('genome', ''),
                       m.get('gene_symbol', ''),
                       m.get('element_name', ''),
                       m.get('class', ''),
                       m.get('subclass', ''),
                       m.get('contig_id', ''),
                       m.get('start', ''),
                       m.get('stop', ''),
                       m.get('strand', ''),
                       m.get('coverage', ''),
                       m.get('identity', ''),
                       m.get('accession', '')]
                f.write('\t'.join(str(x) for x in row) + '\n')
        self.logger.info(f"✓ Mutation TSV: {tsv_file}")
        
        # HTML summary
        self._create_mutation_summary_html(all_mutations, genome_mutation_counts, output_base)
        # JSON summary
        self._create_mutation_json_summaries(all_mutations, genome_mutation_counts, output_base)
    
    def _create_mutation_summary_html(self, all_mutations: List[Dict], genome_counts: Dict[str, int], output_base: str):
        """Create HTML summary for mutations across all genomes"""
        # Group by gene and mutation
        gene_freq = {}
        for m in all_mutations:
            gene = m.get('gene_symbol', 'unknown')
            mutation = m.get('element_name', '')
            key = f"{gene}_{mutation}" if mutation else gene
            if key not in gene_freq:
                gene_freq[key] = {'count': 0, 'genomes': set(), 'gene': gene, 'mutation': mutation,
                                'class': m.get('class',''), 'subclass': m.get('subclass','')}
            gene_freq[key]['count'] += 1
            gene_freq[key]['genomes'].add(m.get('genome',''))
        for k in gene_freq:
            gene_freq[k]['genomes'] = ', '.join(sorted(gene_freq[k]['genomes']))
        sorted_freq = sorted(gene_freq.values(), key=lambda x: x['count'], reverse=True)
        
        random_quote = random.choice(self.science_quotes)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        html = f"""<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>EcoliTyper - Mutation Batch Summary</title>
        <style>
            * {{ margin: 0; padding: 0; box-sizing: border-box; }}
            body {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                color: #ffffff;
                padding: 20px;
                min-height: 100vh;
            }}
            .container {{ max-width: 1400px; margin: 0 auto; }}
            .header {{ text-align: center; margin-bottom: 30px; }}
            .ascii-container {{
                background: rgba(0, 0, 0, 0.7);
                padding: 20px;
                border-radius: 15px;
                margin-bottom: 20px;
                border: 2px solid rgba(0, 255, 0, 0.3);
            }}
            .ascii-art {{
                font-family: 'Courier New', monospace;
                font-size: 10px;
                line-height: 1.1;
                white-space: pre;
                color: #00ff00;
                overflow-x: auto;
            }}
            .quote-container {{
                background: rgba(255, 255, 255, 0.1);
                backdrop-filter: blur(10px);
                padding: 20px;
                border-radius: 10px;
                margin-bottom: 30px;
                text-align: center;
            }}
            .report-section {{
                background: rgba(255, 255, 255, 0.95);
                color: #1f2937;
                padding: 25px;
                border-radius: 10px;
                margin-bottom: 20px;
            }}
            .report-section h2 {{
                color: #1e3a8a;
                border-bottom: 3px solid #3b82f6;
                padding-bottom: 10px;
                margin-bottom: 20px;
            }}
            .summary-table {{
                width: 100%;
                border-collapse: collapse;
                margin-top: 20px;
                font-size: 14px;
            }}
            .summary-table th {{
                background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
                color: white;
                padding: 12px;
                text-align: left;
            }}
            .summary-table td {{
                padding: 12px;
                border-bottom: 1px solid #e5e7eb;
            }}
            .summary-table tr:nth-child(even) {{ background-color: #f8fafc; }}
            .table-responsive {{ overflow-x: auto; margin: 20px 0; }}
            .footer {{
                text-align: center;
                margin-top: 30px;
                padding: 20px;
                background: rgba(0, 0, 0, 0.3);
                border-radius: 10px;
            }}
            .timestamp {{ color: #fbbf24; }}
            .authorship {{ margin-top: 15px; padding: 15px; background: rgba(255, 255, 255, 0.1); border-radius: 8px; font-size: 12px; }}
        </style>
    </head>
    <body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self._get_ascii_art()}</div>
            </div>
            <div class="quote-container">
                <div class="quote-text">"{random_quote}"</div>
            </div>
        </div>

        <div class="report-section">
            <h2>🧬 Mutation Summary Across All Genomes</h2>
            <p>Total genomes with mutations: {len([c for c in genome_counts.values() if c > 0])} / {len(genome_counts)}<br>
            Total mutation events: {len(all_mutations)}</p>
        </div>

        <div class="report-section">
            <h2>📊 Mutation Frequency by Gene/Mutation</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead><tr><th>Gene</th><th>Mutation</th><th>Count</th><th>Genomes</th><th>Class</th><th>Subclass</th></tr></thead>
                    <tbody>
    """
        for item in sorted_freq:
            html += f"""
                        <tr>
                            <td><strong>{item['gene']}</strong></td>
                            <td>{item['mutation']}</td>
                            <td>{item['count']}</td>
                            <td class="sequence-cell">{item['genomes']}</td>
                            <td>{item['class']}</td>
                            <td>{item['subclass']}</td>
                        </tr>
    """
        html += f"""
                    </tbody>
                </table>
            </div>
        </div>

        <div class="footer">
            <p><strong>ECOLITYPER</strong> - Mutation Batch Summary Module</p>
            <p class="timestamp">Generated: {current_time}</p>
            <div class="authorship">
                <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
        </div>
    </div>
    </body>
    </html>"""
        out_file = os.path.join(output_base, "mutation_summary.html")
        with open(out_file, 'w') as f:
            f.write(html)
        self.logger.info(f"✓ Mutation HTML summary: {out_file}")
    
    def _create_mutation_json_summaries(self, all_mutations: List[Dict], genome_counts: Dict[str, int], output_base: str):
        """Create JSON summaries for mutations"""
        genome_summary = {genome: {'total_mutations': count} for genome, count in genome_counts.items()}
        gene_mutation_map = defaultdict(lambda: {'count': 0, 'genomes': set(), 'details': []})
        for m in all_mutations:
            gene = m.get('gene_symbol', 'unknown')
            mut_name = m.get('element_name', '')
            key = f"{gene}_{mut_name}"
            gene_mutation_map[key]['count'] += 1
            gene_mutation_map[key]['genomes'].add(m.get('genome',''))
            gene_mutation_map[key]['details'].append({
                'genome': m.get('genome'),
                'gene': gene,
                'mutation': mut_name,
                'class': m.get('class'),
                'subclass': m.get('subclass'),
                'contig': m.get('contig_id'),
                'start': m.get('start'),
                'stop': m.get('stop')
            })
        for v in gene_mutation_map.values():
            v['genomes'] = list(v['genomes'])
        master_json = {
            'metadata': {
                'tool': 'EcoliTyper AMRfinderPlus Mutation Module',
                'version': self.metadata['version'],
                'analysis_date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                'total_genomes_analyzed': len(genome_counts),
                'total_mutations_detected': len(all_mutations)
            },
            'genome_summary': genome_summary,
            'mutation_frequency': {k: {'count': v['count'], 'genomes': v['genomes']} for k, v in gene_mutation_map.items()},
            'all_mutations': all_mutations
        }
        json_file = os.path.join(output_base, "mutation_master_summary.json")
        with open(json_file, 'w') as f:
            json.dump(master_json, f, indent=2)
        self.logger.info(f"✓ Mutation master JSON: {json_file}")
    
    def create_json_summaries(self, all_results: Dict[str, Any], output_base: str):
        """Create JSON summary files (master and per‑genome)"""
        self.logger.info("Creating JSON summaries...")
        
        # Create master JSON summary
        master_summary = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'amrfinder_version': self.metadata['amrfinder_version'],
                'database_version': self.metadata['database_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_genomes': len(all_results)
            },
            'genome_summaries': {},
            'cross_genome_patterns': {}
        }
        
        # Collect all data for cross-genome analysis
        all_hits_by_gene = defaultdict(lambda: {'count': 0, 'genomes': set()})
        genomes_with_critical = 0
        genomes_with_high_risk = 0
        
        for genome_name, result in all_results.items():
            # Create genome-specific summary
            hits = result['hits']
            genes = [hit.get('gene_symbol', '') for hit in hits if hit.get('gene_symbol', '')]
            unique_genes = set(genes)
            
            critical_genes = [g for g in unique_genes if g in self.critical_risk_genes]
            high_risk_genes = [g for g in unique_genes if g in self.high_risk_genes and g not in self.critical_risk_genes]
            
            if critical_genes:
                genomes_with_critical += 1
            if high_risk_genes:
                genomes_with_high_risk += 1
            
            # Add to genome summaries
            master_summary['genome_summaries'][genome_name] = {
                'total_hits': result['hit_count'],
                'unique_genes': len(unique_genes),
                'critical_genes': critical_genes,
                'high_risk_genes': high_risk_genes,
                'genes': list(unique_genes),
                'status': result['status']
            }
            
            # Update gene frequency
            for gene in unique_genes:
                all_hits_by_gene[gene]['count'] += 1
                all_hits_by_gene[gene]['genomes'].add(genome_name)
        
        # Prepare cross-genome patterns
        cross_genome_data = {}
        for gene, data in all_hits_by_gene.items():
            cross_genome_data[gene] = {
                'frequency': data['count'],
                'genomes': list(data['genomes']),
                'risk_level': 'CRITICAL' if gene in self.critical_risk_genes else 'HIGH' if gene in self.high_risk_genes else 'STANDARD'
            }
        
        master_summary['cross_genome_patterns'] = {
            'total_unique_genes': len(all_hits_by_gene),
            'genomes_with_critical': genomes_with_critical,
            'genomes_with_high_risk': genomes_with_high_risk,
            'gene_frequency': cross_genome_data
        }
        
        # Write master JSON
        master_json_file = os.path.join(output_base, "ecoli_amrfinder_master_summary.json")
        with open(master_json_file, 'w') as f:
            json.dump(master_summary, f, indent=2)
        
        self.logger.info("✓ Master JSON summary created: %s", master_json_file)
        
        # Create individual genome JSON files in their directories
        for genome_name, result in all_results.items():
            genome_dir = os.path.join(output_base, genome_name)
            if os.path.exists(genome_dir):
                json_file = os.path.join(genome_dir, f"{genome_name}_amrfinder_summary.json")
                with open(json_file, 'w') as f:
                    json.dump({
                        'metadata': {
                            'genome': genome_name,
                            'analysis_date': self.metadata['analysis_date']
                        },
                        'summary': {
                            'total_hits': result['hit_count'],
                            'genes': list(set(hit.get('gene_symbol', '') for hit in result['hits'] if hit.get('gene_symbol', '')))
                        },
                        'hits': result['hits'][:10000]  # Limit to first 10000 hits to keep file manageable
                    }, f, indent=2)
    
    def _create_summary_html_report(self, all_results: Dict[str, Any], output_base: str):
        """Create comprehensive HTML summary report with pattern discovery - UPDATED STYLING and ASCII art"""
        
        # Collect all data for pattern analysis
        all_hits = []
        for genome_name, result in all_results.items():
            for hit in result['hits']:
                hit_with_genome = hit.copy()
                hit_with_genome['genome'] = genome_name
                all_hits.append(hit_with_genome)
        
        # Calculate statistics
        total_genomes = len(all_results)
        total_hits = len(all_hits)
        
        # Track critical and high-risk genes across all genomes
        critical_genes_found = set()
        high_risk_genes_found = set()
        genomes_with_critical = 0
        genomes_with_high_risk = 0
        
        # Calculate genes per genome and gene frequency
        genes_per_genome = {}
        gene_frequency = {}
        
        for genome_name, result in all_results.items():
            genome_genes = set()
            for hit in result['hits']:
                gene = hit.get('gene_symbol', '')
                if gene:
                    genome_genes.add(gene)
                    
                    # Track gene frequency
                    if gene not in gene_frequency:
                        gene_frequency[gene] = set()
                    gene_frequency[gene].add(genome_name)
            
            genes_per_genome[genome_name] = genome_genes
            
            # Check for critical and high-risk genes
            has_critical = any(gene in genome_genes for gene in self.critical_risk_genes)
            has_high_risk = any(gene in genome_genes for gene in self.high_risk_genes)
            
            if has_critical:
                genomes_with_critical += 1
                critical_genes_found.update(genome_genes.intersection(self.critical_risk_genes))
            
            if has_high_risk:
                genomes_with_high_risk += 1
                high_risk_genes_found.update(genome_genes.intersection(self.high_risk_genes))
        
        # JavaScript for rotating quotes
        quotes_js = f"""
        <script>
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            
            // Rotate every 10 seconds
            setInterval(rotateQuote, 10000);
            
            // Initial display
            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
            }});
        </script>
        """
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>EcoliTyper AMRfinderPlus - Summary Report</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        body {{ 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
            margin: 0; 
            padding: 0; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}
        .container {{ 
            max-width: 1400px; 
            margin: 0 auto; 
            padding: 20px; 
        }}
        .header {{ 
            background: rgba(255, 255, 255, 0.95); 
            padding: 30px; 
            border-radius: 15px; 
            margin-bottom: 30px; 
            box-shadow: 0 8px 32px rgba(0,0,0,0.1);
            backdrop-filter: blur(10px);
        }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            text-align: center;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 12px;
            line-height: 1.2;
            white-space: pre;
            color: #00ff00;
            text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
            overflow-x: auto;
        }}
        .card {{ 
            background: rgba(255, 255, 255, 0.95); 
            padding: 25px; 
            margin: 20px 0; 
            border-radius: 12px; 
            box-shadow: 0 4px 20px rgba(0,0,0,0.1);
            backdrop-filter: blur(10px);
        }}
        .gene-table {{ 
            width: 100%; 
            border-collapse: collapse; 
            margin: 20px 0; 
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .gene-table th, .gene-table td {{ 
            padding: 15px; 
            text-align: left; 
            border-bottom: 1px solid #e0e0e0; 
        }}
        .gene-table th {{ 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            font-weight: 600;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .summary-stats {{ 
            display: flex; 
            justify-content: space-around; 
            margin: 20px 0; 
            flex-wrap: wrap;
        }}
        .stat-card {{ 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px; 
            border-radius: 12px; 
            text-align: center; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .critical-stat-card {{
            background: linear-gradient(135deg, #dc3545 0%, #c82333 100%);
            color: white;
            padding: 20px; 
            border-radius: 12px; 
            text-align: center; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            color: white;
            padding: 20px;
            border-radius: 12px;
            margin: 20px 0;
            text-align: center;
            font-style: italic;
            border-left: 4px solid #fff;
        }}
        .footer {{
            background: rgba(0, 0, 0, 0.8);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin-top: 40px;
        }}
        .footer a {{
            color: #667eea;
            text-decoration: none;
        }}
        .footer a:hover {{
            text-decoration: underline;
        }}
        .risk-badge {{
            display: inline-block;
            background: #dc3545;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .warning-badge {{
            display: inline-block;
            background: #ffc107;
            color: black;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .safe-badge {{
            display: inline-block;
            background: #28a745;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .present {{ background-color: #d4edda; }}
        .critical {{ background-color: #f8d7da; font-weight: bold; }}
        .high-risk {{ background-color: #fff3cd; }}
        .gene-list-container {{
            font-size: 0.9em;
            line-height: 1.4;
            word-wrap: break-word;
        }}
        .genome-list-simple {{
            font-size: 0.9em;
            line-height: 1.4;
            white-space: normal;
        }}
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">
{self._get_ascii_art()}
                </div>
            </div>
            <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 EcoliTyper AMRfinderPlus - Summary Report</h1>
            <p style="color: #666; font-size: 1.2em;">Comprehensive E. coli Antimicrobial Resistance Analysis Across All Genomes</p>
            <p style="color: #666; font-size: 1.1em;">AMRfinderPlus {self.metadata['amrfinder_version']} | Database: {self.metadata['database_version']}</p>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
"""
        
        # CRITICAL RISK ALERT - Show first if critical genes detected
        if critical_genes_found:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #dc3545; background: #f8d7da;">
            <h2 style="color: #dc3545;">🚨 CRITICAL RISK AMR GENES ACROSS ALL GENOMES</h2>
            <p><strong>{len(critical_genes_found)} unique critical risk genes found in {genomes_with_critical} genomes:</strong></p>
            <div style="margin: 10px 0;">
                <p style="color: #721c24; font-weight: bold;">
                    ⚠️ IMMEDIATE ATTENTION REQUIRED: These genes confer resistance to last-resort antibiotics
                </p>
"""
            for gene in sorted(critical_genes_found):
                html_content += f'<span class="risk-badge">🚨 {gene}</span>'
            html_content += """
            </div>
        </div>
"""
        
        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">📊 Overall Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Genomes</h3>
                    <p style="font-size: 2em; margin: 0;">{total_genomes}</p>
                </div>
                <div class="stat-card">
                    <h3>Total AMR Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{total_hits}</p>
                </div>
                <div class="critical-stat-card">
                    <h3>High-Risk Genomes</h3>
                    <p style="font-size: 2em; margin: 0;">{genomes_with_high_risk}</p>
                </div>
            </div>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>AMRfinderPlus:</strong> {self.metadata['amrfinder_version']}</p>
            <p><strong>Database:</strong> {self.metadata['database_version']}</p>
        </div>
"""
        
        # High-risk genes summary (non-critical)
        if high_risk_genes_found and not critical_genes_found:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #ffc107;">
            <h2 style="color: #856404;">⚠️ High-Risk AMR Genes Detected</h2>
            <p><strong>{len(high_risk_genes_found)} unique high-risk genes found across {genomes_with_high_risk} genomes:</strong></p>
            <div style="margin: 10px 0;">
"""
            for gene in sorted(high_risk_genes_found):
                html_content += f'<span class="warning-badge">{gene}</span>'
            html_content += """
            </div>
        </div>
"""
        
        # Enhanced Genes by Genome table (Showing ALL genes)
        html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🔍 Genes by Genome</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Genome</th>
                        <th>Gene Count</th>
                        <th>Genes Detected</th>
                    </tr>
                </thead>
                <tbody>
"""
        
        for genome in sorted(genes_per_genome.keys()):
            genes = sorted(genes_per_genome.get(genome, set()))
            
            # Create simple comma-separated gene list
            gene_list = ", ".join(genes)
            
            # Determine row class based on risk level
            row_class = "present"
            if any(gene in genes for gene in self.critical_risk_genes):
                row_class = "critical"
            elif any(gene in genes for gene in self.high_risk_genes):
                row_class = "high-risk"
            
            html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{genome}</strong></td>
                        <td>{len(genes)}</td>
                        <td>{gene_list}</td>
                    </tr>
"""
        
        html_content += """
                </tbody>
            </table>
        </div>
        
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">📈 Gene Frequency</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Frequency</th>
                        <th>Prevalence</th>
                        <th>Risk Level</th>
                        <th>Genomes</th>
                    </tr>
                </thead>
                <tbody>
"""
        
        # Calculate gene frequency with SIMPLIFIED STYLING
        for gene, genomes in sorted(gene_frequency.items(), key=lambda x: len(x[1]), reverse=True):
            frequency = len(genomes)
            genome_list = ", ".join(sorted(genomes))
            frequency_percent = (frequency / total_genomes) * 100 if total_genomes > 0 else 0
            
            # Determine risk level - KEEP YOUR LOGIC
            if gene in self.critical_risk_genes:
                risk_level = '<span class="risk-badge">CRITICAL</span>'
                row_class = "critical"
            elif gene in self.high_risk_genes:
                risk_level = '<span class="warning-badge">HIGH</span>'
                row_class = "high-risk"
            else:
                risk_level = '<span class="safe-badge">Standard</span>'
                row_class = "present"
            
            # Determine prevalence badge
            if frequency_percent >= 75:
                prevalence_badge = '<span class="risk-badge">Very High</span>'
            elif frequency_percent >= 50:
                prevalence_badge = '<span class="warning-badge">High</span>'
            elif frequency_percent >= 25:
                prevalence_badge = '<span class="warning-badge">Medium</span>'
            elif frequency_percent >= 10:
                prevalence_badge = '<span class="safe-badge">Low</span>'
            else:
                prevalence_badge = '<span class="safe-badge">Rare</span>'
            
            html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{gene}</strong></td>
                        <td>{frequency} ({frequency_percent:.1f}%)</td>
                        <td>{prevalence_badge}</td>
                        <td>{risk_level}</td>
                        <td class="genome-list-simple">{genome_list}</td>
                    </tr>
"""
        
        html_content += """
                </tbody>
            </table>
        </div>
        
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">📁 Generated Files</h2>
            <ul style="color: #666; font-size: 1.1em;">
                <li><strong>ecoli_amrfinder_summary.tsv</strong> - Complete AMR data for all genomes</li>
                <li><strong>ecoli_amrfinder_statistics_summary.tsv</strong> - Statistical summary</li>
                <li><strong>ecoli_amrfinder_master_summary.json</strong> - Master JSON summary</li>
                <li><strong>mutation_summary.tsv</strong> - All point mutations across genomes</li>
                <li><strong>mutation_summary.html</strong> - Mutation summary report</li>
                <li><strong>mutation_master_summary.json</strong> - Mutation JSON summary</li>
                <li><strong>Individual genome HTML reports</strong> - Detailed analysis per genome</li>
                <li><strong>Individual genome JSON reports</strong> - JSON data per genome</li>
                <li><strong>This summary report</strong> - Cross-genome analysis with pattern discovery</li>
            </ul>
        </div>
        
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #667eea; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
            </p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write summary HTML report
        html_file = os.path.join(output_base, "ecoli_amrfinder_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info("✓ E. coli AMRfinderPlus summary HTML report created: %s", html_file)   
    
    def process_single_genome(self, genome_file: str, output_base: str = "ecoli_amrfinder_results",
                              min_identity: float = None, min_coverage: float = None,
                              report_mutations: bool = True) -> Dict[str, Any]:
        """Process a single E. coli genome with AMRfinderPlus"""
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        
        self.logger.info("=== PROCESSING E. COLI GENOME: %s ===", genome_name)
        
        # Create output directory
        os.makedirs(results_dir, exist_ok=True)
        
        # Check bundled AMRfinderPlus before running
        if not self.check_amrfinder_installed():
            self.logger.error("Bundled AMRfinderPlus not available!")
            return {
                'genome': genome_name,
                'hits': [],
                'hit_count': 0,
                'mutations_file': None,
                'status': 'failed',
                'error': 'Bundled AMRfinderPlus not available'
            }
        
        # Run AMRfinderPlus
        result = self.run_amrfinder_single_genome(genome_file, results_dir,
                                                  min_identity=min_identity,
                                                  min_coverage=min_coverage,
                                                  report_mutations=report_mutations)
        
        status_icon = "✓" if result['status'] == 'success' else "✗"
        self.logger.info("%s %s: %d AMR hits", status_icon, genome_name, result['hit_count'])
        
        return result
    
    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "ecoli_amrfinder_results",
                                 min_identity: float = None, min_coverage: float = None,
                                 report_mutations: bool = True) -> Dict[str, Any]:
        """Process multiple E. coli genomes using wildcard pattern - MAXIMUM SPEED"""
        
        # Print ASCII art banner at start (console)
        print("\n" + "="*80)
        print(self._get_ascii_art())
        print("="*80)
        print("🧬 EcoliTyper AMRfinderPlus - MAXIMUM SPEED MODE")
        print("="*80)
        
        # Find genome files (support all FASTA extensions)
        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta", f"{genome_pattern}.fa", 
                         f"{genome_pattern}.fna", f"{genome_pattern}.faa"]
        
        genome_files = []
        for pattern in fasta_patterns:
            genome_files.extend(glob.glob(pattern))
        
        # Remove duplicates
        genome_files = list(set(genome_files))
        
        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching pattern: {genome_pattern}")
        
        self.logger.info("Found %d E. coli genomes: %s", len(genome_files), [Path(f).name for f in genome_files])
        
        # Create output directory
        os.makedirs(output_base, exist_ok=True)
        
        # Process genomes with threading - MAXIMUM SPEED CONFIGURATION 
        all_results = {}
        
        # Calculate optimal concurrent genomes - BE AGGRESSIVE FOR SPEED
        # Use all available CPU cores for concurrent processing
        max_concurrent = max(1, min(self.cpus, len(genome_files), int(self.available_ram / 1.5)))  # 1.5GB per genome
        
        self.logger.info("🚀 MAXIMUM SPEED: Using %d concurrent genome processing jobs", max_concurrent)
        self.logger.info("   Using BUNDLED AMRfinderPlus: %s", self.bundled_amrfinder)
        self.logger.info("   Using BUNDLED database: %s", self.bundled_database)
        
        with ThreadPoolExecutor(max_workers=max_concurrent) as executor:
            # Submit all tasks
            future_to_genome = {
                executor.submit(self.process_single_genome, genome, output_base,
                                min_identity, min_coverage, report_mutations): genome 
                for genome in genome_files
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_genome):
                genome = future_to_genome[future]
                try:
                    result = future.result()
                    all_results[result['genome']] = result
                    self.logger.info("✓ COMPLETED: %s (%d AMR hits)", result['genome'], result['hit_count'])
                except Exception as e:
                    self.logger.error("✗ FAILED: %s - %s", genome, e)
                    all_results[Path(genome).stem] = {
                        'genome': Path(genome).stem,
                        'hits': [],
                        'hit_count': 0,
                        'mutations_file': None,
                        'status': 'failed'
                    }
        
        # Create AMR summary files and HTML reports after processing all genomes
        self.create_amr_summary(all_results, output_base)
        
        self.logger.info("=== E. COLI AMR ANALYSIS COMPLETE ===")
        self.logger.info("Processed %d genomes", len(all_results))
        self.logger.info("Results saved to: %s", output_base)
        self.logger.info("Bundled AMRfinderPlus used: %s", self.bundled_amrfinder)
        self.logger.info("Bundled database used: %s", self.bundled_database)
        
        return all_results


def main():
    """Command line interface for E. coli AMR analysis"""
    parser = argparse.ArgumentParser(
        description='EcoliTyper AMRfinderPlus Analysis - E. coli Antimicrobial Resistance - MAXIMUM SPEED VERSION',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on all E. coli FASTA files (auto-detect optimal CPU cores - MAXIMUM SPEED)
  python ecoli_amrfinder.py "*.fna"
  
  # Run with custom identity and coverage thresholds
  python ecoli_amrfinder.py "*.fna" --min-identity 0.95 --min-coverage 0.9
  
  # Skip mutation reporting (mutations are reported by default)
  python ecoli_amrfinder.py "*.fna" --skip-mutations
  
  # Force specific number of CPU cores
  python ecoli_amrfinder.py "*.fa" --cpus 16

  # Update database only
  python ecoli_amrfinder.py --update-db

  # Show current database version
  python ecoli_amrfinder.py --db-version

MAXIMUM SPEED RESOURCE MANAGEMENT:
  • 1-4 cores: Uses ALL CPU cores (100% utilization)
  • 5-8 cores: Uses (cores-1) for optimal performance  
  • 9-16 cores: Uses (cores-2) for high performance
  • 17-32 cores: Uses (cores-4) for maximum throughput
  • 32+ cores: Uses 95% of cores (capped at 32)

Supported FASTA extensions: .fasta, .fa, .fna, .faa
        """
    )
    
    # Make pattern optional (nargs='?') so that --update-db and --db-version work without it
    parser.add_argument('pattern', nargs='?', help='File pattern for E. coli genomes (e.g., "*.fasta", "genomes/*.fna")')
    parser.add_argument('--cpus', '-c', type=int, default=None, 
                       help='Number of CPU cores to use (default: auto-detect optimal for MAXIMUM SPEED)')
    parser.add_argument('--output', '-o', default='ecoli_amrfinder_results', 
                       help='Output directory (default: ecoli_amrfinder_results)')
    parser.add_argument('--min-identity', type=float, default=None,
                       help='Minimum identity (0..1) for hits. Default: AMRfinder auto threshold')
    parser.add_argument('--min-coverage', type=float, default=None,
                       help='Minimum coverage of reference (0..1). Default: 0.5')
    parser.add_argument('--skip-mutations', action='store_true',
                       help='Skip point mutation reporting (mutations are reported by default)')
    parser.add_argument('--update-db', action='store_true', 
                       help='Update AMRfinderPlus database to latest version and exit')
    parser.add_argument('--db-version', action='store_true', 
                       help='Show current database version and exit')
    
    args = parser.parse_args()
    
    # Handle database operations without requiring pattern
    if args.update_db or args.db_version:
        executor = EcoliAMRfinderPlus(cpus=args.cpus)
        if args.update_db:
            print("Updating AMRfinderPlus database...")
            success = executor.update_database()
            if success:
                print("Database updated successfully.")
            else:
                print("Database update failed.")
            sys.exit(0)
        if args.db_version:
            print(f"Database version: {executor.metadata['database_version']}")
            print(f"Database path: {executor.bundled_database or 'Not found'}")
            sys.exit(0)
    
    # For analysis, pattern is required
    if not args.pattern:
        parser.error("Please provide a file pattern for genomes (or use --update-db / --db-version)")
    
    executor = EcoliAMRfinderPlus(cpus=args.cpus)
    
    try:
        results = executor.process_multiple_genomes(args.pattern, args.output,
                                                    min_identity=args.min_identity,
                                                    min_coverage=args.min_coverage,
                                                    report_mutations=not args.skip_mutations)
        
        # Print summary
        executor.logger.info("\n" + "="*50)
        executor.logger.info("🧬 EcoliTyper AMRfinderPlus FINAL SUMMARY")
        executor.logger.info("="*50)
        
        total_hits = 0
        high_risk_count = 0
        critical_risk_count = 0
        
        for genome_name, result in results.items():
            total_hits += result['hit_count']
            
            # Count high-risk and critical genes
            genes = [hit.get('gene_symbol') for hit in result['hits'] if hit.get('gene_symbol')]
            high_risk_count += sum(1 for gene in genes if gene in executor.high_risk_genes)
            critical_risk_count += sum(1 for gene in genes if gene in executor.critical_risk_genes)
            
            executor.logger.info("✓ %s: %d AMR hits", genome_name, result['hit_count'])
        
        executor.logger.info("\n📊 E. COLI SUMMARY STATISTICS:")
        executor.logger.info("   Total genomes processed: %d", len(results))
        executor.logger.info("   Total AMR hits: %d", total_hits)
        executor.logger.info("   High-risk genes detected: %d", high_risk_count)
        executor.logger.info("   CRITICAL RISK genes detected: %d", critical_risk_count)
        executor.logger.info("   Average AMR hits per genome: %.1f", total_hits / len(results) if results else 0)
        
        # Show summary file locations
        executor.logger.info("\n📁 SUMMARY FILES CREATED:")
        executor.logger.info("   Comprehensive AMR data: %s/ecoli_amrfinder_summary.tsv", args.output)
        executor.logger.info("   Statistics summary: %s/ecoli_amrfinder_statistics_summary.tsv", args.output)
        executor.logger.info("   Master JSON summary: %s/ecoli_amrfinder_master_summary.json", args.output)
        executor.logger.info("   Summary HTML report: %s/ecoli_amrfinder_summary_report.html", args.output)
        executor.logger.info("   Mutation TSV: %s/mutation_summary.tsv", args.output)
        executor.logger.info("   Mutation HTML summary: %s/mutation_summary.html", args.output)
        executor.logger.info("   Mutation master JSON: %s/mutation_master_summary.json", args.output)
        executor.logger.info("   Individual genome reports in: %s/*/", args.output)
        
        # Performance summary
        executor.logger.info("\n⚡ MAXIMUM SPEED PERFORMANCE SUMMARY:")
        executor.logger.info("   CPU cores utilized: %d cores", executor.cpus)
        executor.logger.info("   Available RAM: %.1f GB", executor.available_ram)
        executor.logger.info("   Processing mode: MAXIMUM SPEED CONCURRENT MODE 🚀")
        executor.logger.info("   Strategy: Process multiple genomes concurrently with optimal core allocation")
        executor.logger.info("   Bundled AMRfinderPlus: %s", executor.metadata['amrfinder_version'])
        executor.logger.info("   Bundled database: %s", executor.metadata['database_version'])
        
        # Critical risk warning if detected
        if critical_risk_count > 0:
            executor.logger.info("\n🚨 CRITICAL RISK ALERT: Last-resort antibiotic resistance genes detected!")
            executor.logger.info("   Immediate clinical attention and infection control measures required.")
        
        import random
        executor.logger.info("\n💡 %s", random.choice(executor.science_quotes))
        
    except Exception as e:
        executor.logger.error("E. coli AMR analysis failed: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()