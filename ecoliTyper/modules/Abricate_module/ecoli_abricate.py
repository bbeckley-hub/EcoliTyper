#!/usr/bin/env python3
"""
EcoliTyper ABRicate Standalone Module
Comprehensive ABRicate analysis for E. coli with HTML reporting - MAXIMUM SPEED VERSION
Author: Beckley Brown <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
Date: 2025
Send a quick mail for any issues or further explanations.
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any
import argparse
import re
from datetime import datetime
import psutil
import math
import json
import pandas as pd

class EcoliAbricateExecutor:
    """ABRicate executor for E. coli with comprehensive HTML reporting - MAXIMUM SPEED"""
    
    def __init__(self, cpus: int = None):
        # Setup logging FIRST
        self.logger = self._setup_logging()
        
        # Initialize available_ram before calculating cpus
        self.available_ram = self._get_available_ram()
        
        # Then calculate resources - MAXIMUM SPEED MODE
        self.cpus = self._calculate_optimal_cpus(cpus)
        
        # E. coli specific databases
        self.required_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot', 
            'plasmidfinder', 'megares', 'ecoh', 'ecoli_vf'
        ]
        
        # COMPREHENSIVE E. coli high-risk resistance genes
        self.high_risk_genes = {
            # 🔴 CRITICAL ESBL genes 
            'blaCTX-M', 'blaSHV', 'blaTEM', 'BlaAmpC1_Ecoli', 'BlaAmpC2_Ecoli', 
            'BlaPenicillin_Binding_Protein_Ecoli', 'BlaampH_Ecoli', 'blaTEM-105',
            'blaTEM-150', 'blaTEM-135', 'blaCTX-M-14', 'blaCTX-M-1', 'blaTEM-1',
            '(Bla)blaTEM-150', '(Bla)blaCTX-M-14', '(Bla)blaTEM-135', '(Bla)blaCTX-M-1',
            
            # 🔴 CRITICAL Carbapenemases
            'blaKPC', 'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP',
            
            # 🔴 CRITICAL Colistin resistance
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            '(Col)mcr-1.1', 'MCR', 'CTX',
            
            # 🔴 CRITICAL Fluoroquinolone resistance
            'qnrA', 'qnrB', 'qnrC', 'qnrD', 'qnrS', 'qnrVC', 'aac(6\')-Ib-cr',
            
            # 🟡 HIGH RISK Aminoglycoside resistance
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'rmtG', 'rmtH',
            'AGlyaadA2', 'AGlyaph3-Ia', 'AGlystrA', 'AGlystrB', 'AGlyaac3-IId',
            'aadA2', 'aph3-Ia', 'strA', 'strB', 'aac3-IId', 'aadA1-pm', 'ant3-Ia',
            'APH3-Ib', 'APH3-Ia', 'APH6-Id', 'AAC3-IId', 'ANT3-Ia',
            
            # 🟡 HIGH RISK Other resistance classes
            'PhecatA2', 'PhefloR', 'catA2', 'floR', 'catII',
            'Sulsul2', 'Sulsul3', 'Sulsul1', 'sul2', 'sul3', 'sul1',
            'TettetA', 'TettetR', 'TettetC', 'tetA', 'tetR', 'tetC',
            'TmtdfrA14', 'TmtdfrA1', 'dfrA14', 'dfrA1',
            
            # 🟡 HIGH RISK Other high-risk
            'cfr', 'optrA', 'poxtA', 'vat', 'vga', 'mphB'
        }
        
        # COMPREHENSIVE E. coli virulence factors
        self.virulence_genes = {
            # 🔴 CRITICAL RISK - Toxins & Lethal Factors
            'stx1', 'stx2', 'stx1A', 'stx1B', 'stx2A', 'stx2B', 'stxA', 'stx2d1A',
            'cnf1', 'hlyA', 'hlyB', 'hlyC', 'hlyD', 'astA', 'east1', 
            'eltA', 'eltB', 'estIa', 'toxB',
            
            # 🟡 HIGH RISK - Adhesins & Colonization Factors
            'eae', 'tir', 'paa', 'fdeC', 
            'afaA', 'afaB', 'afaC', 'afaD', 'afaE',
            'papA', 'papB', 'papC', 'papD', 'papE', 'papF', 'papG', 'papH', 'papI', 'papJ', 'papK', 'papX',
            'sfaA', 'sfaB', 'sfaC', 'sfaD', 'sfaE', 'sfaF', 'sfaG', 'sfaH', 'sfaS', 'sfaX', 'sfaY',
            'fimA', 'fimB', 'fimC', 'fimD', 'fimE', 'fimF', 'fimG', 'fimH', 'fimI',
            'faeC', 'faeD', 'faeE', 'faeF', 'faeG', 'faeH', 'faeI', 'faeJ',
            'fanC', 'fanD', 'fanE', 'fanF', 'fanG', 'fanH',
            
            # 🟡 HIGH RISK - Iron Acquisition Systems
            'fyuA', 'irp1', 'irp2', 'ybtA', 'ybtE', 'ybtP', 'ybtQ', 'ybtS', 'ybtT', 'ybtU', 'ybtX',
            'iucA', 'iucB', 'iucC', 'iucD', 'iutA',
            'chuA', 'chuS', 'chuT', 'chuU', 'chuV', 'chuW', 'chuX', 'chuY',
            'iroB', 'iroC', 'iroD', 'iroE', 'iroN',
            'shuA', 'shuS', 'shuT', 'shuX', 'shuY',
            'sitA', 'sitB', 'sitC', 'sitD',
            
            # 🟡 HIGH RISK - Type III Secretion System (T3SS) & Effectors
            'escC', 'escD', 'escE', 'escF', 'escG', 'escI', 'escJ', 'escL', 'escN', 'escO', 'escP', 'escR', 'escS', 'escT', 'escU', 'escV',
            'espA', 'espB', 'espD', 'espF', 'espG', 'espH', 'espJ', 'espK', 'espZ',
            'nleA', 'nleB1', 'nleB2', 'nleC', 'nleD', 'nleE', 'nleF', 'nleG', 'nleH1', 'nleH2',
            'map', 'tir', 'espFu', 'tccP', 'tccP2',
            'cesA', 'cesB', 'cesD', 'cesD2', 'cesF', 'cesL', 'cesT',
            
            # 🟡 HIGH RISK - Capsule & Immune Evasion
            'kpsD', 'kpsM', 'kpsT',
            'ibeA', 'ibeB', 'ibeC',
            'vat', 'senB', 'cif', 'ospG',
            'iss2', 'ompA',
            
            # 🟡 HIGH RISK - Other Virulence Mechanisms
            'efa1', 'stcE', 'pic', 'sepD', 'sepL', 'sepQ', 'sepZ',
            'aslA', 'cdtA', 'cdtB', 'cdtC', 'cdiA', 'cdiB',
            
            # 🔵 MODERATE RISK - Enterobactin & Basic Systems  
            'entA', 'entB', 'entC', 'entD', 'entE', 'entF', 'entS',
            'fepA', 'fepB', 'fepC', 'fepD', 'fepE', 'fepG', 'fes',
            'csgA', 'csgB', 'csgC', 'csgD', 'csgE', 'csgF', 'csgG',
            'gspC', 'gspD', 'gspE', 'gspF', 'gspG', 'gspH', 'gspI', 'gspJ', 'gspK', 'gspL', 'gspM',
            'fliC', 'fliD', 'fliE', 'fliF', 'fliG', 'fliH', 'fliI', 'fliJ', 'fliK', 'fliL', 'fliM', 
            'fliN', 'fliO', 'fliP', 'fliQ', 'fliR', 'fliS', 'fliT', 'fliY', 'fliZ',
            'flgA', 'flgB', 'flgC', 'flgD', 'flgE', 'flgF', 'flgG', 'flgH', 'flgI', 'flgJ', 'flgK', 'flgL', 'flgN',
            'flhA', 'flhB', 'flhC', 'flhD', 'flhE',
            'motA', 'motB', 'cheA', 'cheB', 'cheR', 'cheW', 'cheY', 'cheZ'
        }

        # Enhanced risk categories for better reporting
        self.critical_resistance_genes = {
            'blaCTX-M', 'blaSHV', 'blaTEM', 'blaKPC', 'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP',
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            '(Col)mcr-1.1', 'MCR'
        }
        
        self.critical_virulence_genes = {
            'stx1', 'stx2', 'stx1A', 'stx1B', 'stx2A', 'stx2B', 'stxA', 'stx2d1A',
            'cnf1', 'hlyA', 'hlyB', 'hlyC', 'hlyD', 'eae'
        }
        
        self.metadata = {
            "tool_name": "EcoliTyper ABRicate",
            "version": "1.0.0", 
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        
        self.science_quotes = [
            "“The important thing is not to stop questioning. Curiosity has its own reason for existence.” - Albert Einstein",
            "“Nothing in life is to be feared, it is only to be understood.” - Marie Curie", 
            "“The microscope opens a new world to the investigator.” - Robert Koch",
            "“In science, the credit goes to the man who convinces the world, not to the man to whom the idea first occurs.” - Francis Darwin",
            "“The good thing about science is that it's true whether or not you believe in it.” - Neil deGrasse Tyson",
            "“Science knows no country, because knowledge belongs to humanity.” - Louis Pasteur"
        ]
    
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
                optimal_cpus = max(8, total_physical_cores - 2)  # Use 14/16, 13/15, etc.
            elif total_physical_cores <= 32:
                optimal_cpus = max(16, total_physical_cores - 4)  # Use 28/32, 27/31, etc.
            else:
                optimal_cpus = min(32, int(total_physical_cores * 0.85))  # Use 85% on huge systems
            
            # Ensure at least 1 CPU and not more than available cores
            optimal_cpus = max(1, min(optimal_cpus, total_physical_cores))
            
            self._log_resource_info(optimal_cpus, total_physical_cores)
            return optimal_cpus
            
        except Exception as e:
            # Fallback to using all available cores for maximum speed
            self.logger.warning(f"Could not detect CPU cores, using maximum available: {e}")
            return os.cpu_count() or 4
    
    def _log_resource_info(self, cpus: int, total_cores: int = None):
        """Log resource allocation information"""
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
            self.logger.info("💡 Performance: High-speed mode")
        else:
            self.logger.info("💡 Performance: MAXIMUM SPEED MODE 🚀")

    def check_abricate_installed(self) -> bool:
        """Check if ABRicate is installed and meets version requirements"""
        try:
            result = subprocess.run(['abricate', '--version'], 
                                  capture_output=True, text=True, check=True)
            version_line = result.stdout.strip()
            self.logger.info("ABRicate version: %s", version_line)
            
            version_match = re.search(r'(\d+\.\d+\.\d+)', version_line)
            if version_match:
                version_str = version_match.group(1)
                if version_str >= "1.0.1":
                    self.logger.info("✓ ABRicate version meets requirement (>=1.0.1)")
                    return True
                else:
                    self.logger.error("ABRicate version too old: %s. Required >=1.0.1", version_str)
                    return False
            self.logger.info("✓ ABRicate installed (version check skipped)")
            return True
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error("ABRicate not found. Please install with: conda install -c bioconda abricate")
            return False
    
    def setup_abricate_databases(self):
        """Setup ABRicate databases if they don't exist"""
        self.logger.info("Setting up ABRicate databases for E. coli analysis...")
        
        available_dbs = []
        missing_dbs = []
        
        try:
            # Check which databases exist
            check_result = subprocess.run(['abricate', '--list'], 
                                        capture_output=True, text=True, check=True)
            
            for db in self.required_databases:
                if db in check_result.stdout:
                    self.logger.info("✓ Database available: %s", db)
                    available_dbs.append(db)
                else:
                    self.logger.warning("Database not available: %s", db)
                    missing_dbs.append(db)
            
            # Setup missing databases
            for db in missing_dbs:
                self.logger.info("Attempting to setup database: %s", db)
                try:
                    result = subprocess.run(
                        ['abricate', '--setupdb', '--db', db],
                        capture_output=True, text=True, check=True
                    )
                    self.logger.info("✓ Database setup completed: %s", db)
                    available_dbs.append(db)
                except subprocess.CalledProcessError as e:
                    self.logger.error("Failed to setup database %s: %s", db, e.stderr)
            
            # Update required databases to only include available ones
            self.required_databases = available_dbs
            self.logger.info("Using databases: %s", ", ".join(self.required_databases))
            
        except subprocess.CalledProcessError as e:
            self.logger.error("Error checking ABRicate databases: %s", e.stderr)
        except Exception as e:
            self.logger.error("Unexpected error setting up databases: %s", e)
    
    def run_abricate_single_db(self, genome_file: str, database: str, output_dir: str) -> Dict[str, Any]:
        """Run ABRicate on a single genome with specific database"""
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"abricate_{database}.txt")
        
        cmd = [
            'abricate', 
            genome_file, 
            '--db', database,
            '--minid', '80',
            '--mincov', '80'
        ]
        
        self.logger.info("Running ABRicate: %s --db %s", genome_name, database)
        
        try:
            with open(output_file, 'w') as outfile:
                result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True, check=True)
            
            # Parse results for reporting
            hits = self._parse_abricate_output(output_file)
            
            # Create individual database HTML report
            self._create_database_html_report(genome_name, database, hits, output_dir)
            
            return {
                'database': database,
                'genome': genome_name,
                'output_file': output_file,
                'hits': hits,
                'hit_count': len(hits),
                'status': 'success'
            }
            
        except subprocess.CalledProcessError as e:
            self.logger.error("ABRicate failed for %s on %s: %s", database, genome_name, e.stderr)
            return {
                'database': database,
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'status': 'failed'
            }
    
    def _parse_abricate_output(self, abricate_file: str) -> List[Dict]:
        """Parse ABRicate output file - ROBUST VERSION that handles tabs in fields"""
        hits = []
        try:
            with open(abricate_file, 'r') as f:
                lines = f.readlines()
                
            if not lines:
                return hits
                
            # Find header line
            headers = []
            data_lines = []
            
            for line in lines:
                if line.startswith('#FILE') and not headers:
                    # This is the header line
                    headers = line.strip().replace('#', '').split('\t')
                elif line.strip() and not line.startswith('#'):
                    data_lines.append(line.strip())
            
            if not headers:
                self.logger.warning("No headers found in %s", abricate_file)
                return hits
            
            # Expected column count based on standard ABRicate output
            expected_columns = len(headers)
                
            # Parse data lines with robust tab handling
            for line_num, line in enumerate(data_lines, 1):
                # Split by tab but be careful about fields that contain tabs
                parts = line.split('\t')
                
                # Handle cases where there are more parts than headers due to tabs in fields
                if len(parts) > expected_columns:
                    # Combine extra fields into the last column (usually PRODUCT)
                    combined_parts = parts[:expected_columns-1]  # Take all but the last expected column
                    combined_parts.append('\t'.join(parts[expected_columns-1:]))  # Combine the rest into PRODUCT
                    parts = combined_parts
                elif len(parts) < expected_columns:
                    # Pad with empty strings if fewer columns
                    parts.extend([''] * (expected_columns - len(parts)))
                
                if len(parts) == expected_columns:
                    hit = {}
                    for i, header in enumerate(headers):
                        hit[header] = parts[i] if i < len(parts) else ''
                    
                    # Map to consistent field names
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
                        'resistance': hit.get('RESISTANCE', '')
                    }
                    hits.append(processed_hit)
                else:
                    self.logger.warning("Line %d has %d parts, expected %d: %s", 
                                      line_num, len(parts), expected_columns, line[:100] + "...")
                    
        except Exception as e:
            self.logger.error("Error parsing %s: %s", abricate_file, e)
            
        self.logger.info("Parsed %d hits from %s", len(hits), abricate_file)
        return hits
    
    def _create_database_html_report(self, genome_name: str, database: str, hits: List[Dict], output_dir: str):
        """Create individual HTML report for each database with beautiful styling"""
        
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
    <title>EcoliTyper ABRicate - {database.upper()} Database</title>
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
        .present {{ background-color: #d4edda; }}
        .high-risk {{ background-color: #fff3cd; }}
        .critical {{ background-color: #f8d7da; font-weight: bold; }}
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 EcoliTyper ABRicate - {database.upper()} Database</h1>
            <p style="color: #666; font-size: 1.2em;">Genome: {genome_name} | Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
        
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">📊 Database Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{len(hits)}</p>
                </div>
                <div class="stat-card">
                    <h3>Database</h3>
                    <p style="font-size: 1.5em; margin: 0;">{database.upper()}</p>
                </div>
            </div>
        </div>
"""
        
        if hits:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🔍 Genes Detected</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Product</th>
                        <th>Coverage</th>
                        <th>Identity</th>
                        <th>Accession</th>
                    </tr>
                </thead>
                <tbody>
"""
            
            for hit in hits:
                # Determine row class based on gene risk
                row_class = "present"
                gene_base = hit['gene'].split('-')[0]  # Get base gene name (handle blaTEM-1)
                
                if any(hr_gene in gene_base for hr_gene in self.high_risk_genes):
                    row_class = "critical"
                elif any(vf_gene in gene_base for vf_gene in self.virulence_genes):
                    row_class = "high-risk"
                
                # Truncate very long product descriptions for display
                product_display = hit['product']
                if len(product_display) > 150:
                    product_display = product_display[:147] + "..."
                
                html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{hit['gene']}</strong></td>
                        <td title="{hit['product']}">{product_display}</td>
                        <td>{hit['coverage_percent']}%</td>
                        <td>{hit['identity_percent']}%</td>
                        <td>{hit['accession']}</td>
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
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">✅ No Genes Detected</h2>
            <p>No significant hits found in the {database.upper()} database.</p>
        </div>
"""
        
        html_content += f"""
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #667eea; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
                Analysis performed using EcoliTyper ABRicate v1.0.1
            </p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write individual database HTML report
        html_file = os.path.join(output_dir, f"abricate_{database}_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info("Individual database report: %s", html_file)
    
    def analyze_ecoli_resistance(self, all_hits: List[Dict]) -> Dict[str, Any]:
        """Enhanced E. coli resistance analysis with comprehensive risk assessment"""
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
            'total_hits': len(all_hits)
        }
        
        for hit in all_hits:
            gene = hit['gene']
            gene_base = gene.split('-')[0] if '-' in gene else gene
            
            # Check for CRITICAL resistance patterns
            if any(crit_gene in gene_base for crit_gene in self.critical_resistance_genes):
                if any(esbl in gene_base for esbl in ['blaCTX-M', 'blaSHV', 'blaTEM']):
                    analysis['esbl_status'] = 'positive'
                    risk_level = 'ESBL'
                elif any(carba in gene_base for carba in ['blaKPC', 'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP']):
                    analysis['carbapenemase_status'] = 'positive'
                    risk_level = 'CARBAPENEMASE'
                elif any(mcr in gene for mcr in ['mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10', 'MCR']):
                    analysis['colistin_resistance'] = 'positive'
                    risk_level = 'COLISTIN-RES'
                else:
                    risk_level = 'CRITICAL'
                
                analysis['critical_resistance_genes'].append({
                    'gene': gene,
                    'product': hit['product'],
                    'database': hit['database'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent'],
                    'risk_level': risk_level
                })
            
            # Check for HIGH RISK resistance genes
            elif any(hr_gene in gene_base for hr_gene in self.high_risk_genes):
                analysis['high_risk_resistance_genes'].append({
                    'gene': gene,
                    'product': hit['product'],
                    'database': hit['database'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent'],
                    'risk_level': 'HIGH-RISK'
                })
            
            # Check for CRITICAL virulence genes
            if any(crit_vf in gene_base for crit_vf in self.critical_virulence_genes):
                analysis['critical_virulence_genes'].append({
                    'gene': gene,
                    'product': hit['product'],
                    'database': hit['database'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent'],
                    'risk_level': 'CRITICAL-VIRULENCE'
                })
            
            # Check for HIGH RISK virulence genes
            elif any(vf_gene in gene_base for vf_gene in self.virulence_genes):
                # Check if it's high risk virulence (not moderate)
                is_high_risk_virulence = any(hr_vf in gene_base for hr_vf in [
                    'stx', 'cnf', 'hly', 'ast', 'east', 'elt', 'est', 'eae', 'tir', 'paa',
                    'afa', 'pap', 'sfa', 'fim', 'fae', 'fan', 'fyu', 'irp', 'ybt',
                    'iuc', 'iut', 'chu', 'iro', 'shu', 'esc', 'esp', 'nle', 'map',
                    'kps', 'ibe', 'vat', 'sen', 'cif', 'efa', 'stc', 'pic'
                ])
                
                if is_high_risk_virulence:
                    analysis['high_risk_virulence_genes'].append({
                        'gene': gene,
                        'product': hit['product'],
                        'database': hit['database'],
                        'coverage': hit['coverage_percent'],
                        'identity': hit['identity_percent'],
                        'risk_level': 'HIGH-VIRULENCE'
                    })
                else:
                    analysis['moderate_risk_genes'].append({
                        'gene': gene,
                        'product': hit['product'],
                        'database': hit['database'],
                        'coverage': hit['coverage_percent'],
                        'identity': hit['identity_percent'],
                        'risk_level': 'MODERATE'
                    })
            
            # Track resistance classes
            resistance_class = self._classify_resistance(hit['product'])
            if resistance_class:
                if resistance_class not in analysis['resistance_classes']:
                    analysis['resistance_classes'][resistance_class] = []
                if gene not in [g['gene'] for g in analysis['resistance_classes'][resistance_class]]:
                    analysis['resistance_classes'][resistance_class].append({
                        'gene': gene,
                        'product': hit['product']
                    })
        
        # Calculate totals
        analysis['total_critical_resistance'] = len(analysis['critical_resistance_genes'])
        analysis['total_high_risk_resistance'] = len(analysis['high_risk_resistance_genes'])
        analysis['total_critical_virulence'] = len(analysis['critical_virulence_genes'])
        analysis['total_high_risk_virulence'] = len(analysis['high_risk_virulence_genes'])
        
        return analysis

    def _classify_resistance(self, product: str) -> str:
        """Enhanced resistance classification"""
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
        else:
            return 'Other resistance'
    
    def create_comprehensive_html_report(self, genome_name: str, results: Dict, output_dir: str):
        """Create comprehensive HTML report for E. coli with beautiful styling"""
        
        # Collect all hits
        all_hits = []
        for db_result in results.values():
            all_hits.extend(db_result['hits'])
        
        # Analyze E. coli resistance
        analysis = self.analyze_ecoli_resistance(all_hits)
        
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
    <title>EcoliTyper ABRicate Analysis Report</title>
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
        .success {{ color: #28a745; font-weight: 600; }}
        .warning {{ color: #ffc107; font-weight: 600; }}
        .error {{ color: #dc3545; font-weight: 600; }}
        .critical {{ background-color: #f8d7da; font-weight: bold; }}
        .high-risk {{ background-color: #fff3cd; }}
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
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 EcoliTyper ABRicate Analysis Report</h1>
            <p style="color: #666; font-size: 1.2em;">Comprehensive E. coli Antimicrobial Resistance & Virulence Analysis</p>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
        
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">📊 E. coli AMR Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['total_hits']}</p>
                </div>
                <div class="stat-card">
                    <h3>Critical Resistance</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['total_critical_resistance']}</p>
                </div>
                <div class="stat-card">
                    <h3>Critical Virulence</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['total_critical_virulence']}</p>
                </div>
            </div>
            <p><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
        </div>
"""
        
        # Critical resistance alerts
        critical_alerts = []
        if analysis['esbl_status'] == 'positive':
            critical_alerts.append("🟡 ESBL DETECTED")
        if analysis['carbapenemase_status'] == 'positive':
            critical_alerts.append("🔴 CARBAPENEMASE DETECTED")
        if analysis['colistin_resistance'] == 'positive':
            critical_alerts.append("🔴 COLISTIN RESISTANCE DETECTED")
        
        if critical_alerts:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #dc3545;">
            <h2 style="color: #dc3545;">⚠️ CRITICAL RESISTANCE ALERTS</h2>
            <div style="margin: 10px 0;">
"""
            for alert in critical_alerts:
                html_content += f'<span class="risk-badge">{alert}</span>'
            html_content += """
            </div>
        </div>
"""
        
        # Resistance classes summary
        if analysis['resistance_classes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🧪 Resistance Classes Detected</h2>
            <div style="margin: 20px 0;">
"""
            
            for class_name, genes in analysis['resistance_classes'].items():
                gene_list = ", ".join([g['gene'] for g in genes])
                html_content += f"""
                <div style="margin: 10px 0; padding: 10px; background: #f8f9fa; border-radius: 8px;">
                    <strong style="color: #667eea;">{class_name}</strong> ({len(genes)} genes)
                    <br><span style="color: #666; font-size: 0.9em;">{gene_list}</span>
                </div>
"""
            
            html_content += "</div></div>"
        
        # Critical resistance genes table
        if analysis['critical_resistance_genes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🔴 CRITICAL Resistance Genes</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Product</th>
                        <th>Database</th>
                        <th>Coverage</th>
                        <th>Identity</th>
                        <th>Risk Level</th>
                    </tr>
                </thead>
                <tbody>
"""
            
            for gene_info in analysis['critical_resistance_genes']:
                product_display = gene_info['product']
                if len(product_display) > 100:
                    product_display = gene_info['product'][:97] + "..."
                
                html_content += f"""
                    <tr class="critical">
                        <td><strong>{gene_info['gene']}</strong></td>
                        <td title="{gene_info['product']}">{product_display}</td>
                        <td>{gene_info['database']}</td>
                        <td>{gene_info['coverage']}%</td>
                        <td>{gene_info['identity']}%</td>
                        <td><span class="risk-badge">{gene_info['risk_level']}</span></td>
                    </tr>
"""
            
            html_content += """
                </tbody>
            </table>
        </div>
"""
        
        # High-risk resistance genes table
        if analysis['high_risk_resistance_genes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🟡 High-Risk Resistance Genes</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Product</th>
                        <th>Database</th>
                        <th>Coverage</th>
                        <th>Identity</th>
                        <th>Risk Level</th>
                    </tr>
                </thead>
                <tbody>
"""
            
            for gene_info in analysis['high_risk_resistance_genes']:
                product_display = gene_info['product']
                if len(product_display) > 100:
                    product_display = gene_info['product'][:97] + "..."
                
                html_content += f"""
                    <tr class="high-risk">
                        <td><strong>{gene_info['gene']}</strong></td>
                        <td title="{gene_info['product']}">{product_display}</td>
                        <td>{gene_info['database']}</td>
                        <td>{gene_info['coverage']}%</td>
                        <td>{gene_info['identity']}%</td>
                        <td><span class="warning-badge">{gene_info['risk_level']}</span></td>
                    </tr>
"""
            
            html_content += """
                </tbody>
            </table>
        </div>
"""
        
        # Critical virulence genes table
        if analysis['critical_virulence_genes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🔴 CRITICAL Virulence Factors</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Product</th>
                        <th>Database</th>
                        <th>Coverage</th>
                        <th>Identity</th>
                        <th>Risk Level</th>
                    </tr>
                </thead>
                <tbody>
"""
            
            for gene_info in analysis['critical_virulence_genes']:
                product_display = gene_info['product']
                if len(product_display) > 100:
                    product_display = gene_info['product'][:97] + "..."
                
                html_content += f"""
                    <tr class="critical">
                        <td><strong>{gene_info['gene']}</strong></td>
                        <td title="{gene_info['product']}">{product_display}</td>
                        <td>{gene_info['database']}</td>
                        <td>{gene_info['coverage']}%</td>
                        <td>{gene_info['identity']}%</td>
                        <td><span class="risk-badge">{gene_info['risk_level']}</span></td>
                    </tr>
"""
            
            html_content += """
                </tbody>
            </table>
        </div>
"""
        
        # High-risk virulence genes table
        if analysis['high_risk_virulence_genes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🟡 High-Risk Virulence Factors</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Product</th>
                        <th>Database</th>
                        <th>Coverage</th>
                        <th>Identity</th>
                        <th>Risk Level</th>
                    </tr>
                </thead>
                <tbody>
"""
            
            for gene_info in analysis['high_risk_virulence_genes']:
                product_display = gene_info['product']
                if len(product_display) > 100:
                    product_display = gene_info['product'][:97] + "..."
                
                html_content += f"""
                    <tr class="high-risk">
                        <td><strong>{gene_info['gene']}</strong></td>
                        <td title="{gene_info['product']}">{product_display}</td>
                        <td>{gene_info['database']}</td>
                        <td>{gene_info['coverage']}%</td>
                        <td>{gene_info['identity']}%</td>
                        <td><span class="warning-badge">{gene_info['risk_level']}</span></td>
                    </tr>
"""
            
            html_content += """
                </tbody>
            </table>
        </div>
"""
        
        # Database summary
        html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🗃️ Database Results Summary</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Database</th>
                        <th>Hits</th>
                        <th>Status</th>
                    </tr>
                </thead>
                <tbody>
"""
        
        for db, result in results.items():
            status_icon = "✅" if result['status'] == 'success' else "❌"
            html_content += f"""
                    <tr>
                        <td>{db}</td>
                        <td>{result['hit_count']}</td>
                        <td>{status_icon} {result['status']}</td>
                    </tr>
"""
        
        html_content += """
                </tbody>
            </table>
        </div>
        
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #667eea; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
                Analysis performed using EcoliTyper ABRicate v1.0.1
            </p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write comprehensive HTML report
        html_file = os.path.join(output_dir, f"{genome_name}_comprehensive_abricate_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info("Comprehensive E. coli HTML report generated: %s", html_file)
    
    def create_database_summaries(self, all_results: Dict[str, Any], output_base: str):
        """Create ABRicate summary files and HTML reports for each database across all genomes"""
        self.logger.info("Creating E. coli database summary files and HTML reports...")
        
        # Group results by database
        db_results = {}
        for genome_name, genome_result in all_results.items():
            for db, db_result in genome_result['results'].items():
                if db not in db_results:
                    db_results[db] = []
                
                # Add hits with genome name
                for hit in db_result['hits']:
                    hit_with_genome = hit.copy()
                    hit_with_genome['genome'] = genome_name
                    db_results[db].append(hit_with_genome)
        
        # Create summary file and HTML report for each database
        for db, hits in db_results.items():
            if hits:
                # Create TSV summary
                summary_file = os.path.join(output_base, f"ecoli_{db}_abricate_summary.tsv")
                
                # Get headers from first hit
                headers = list(hits[0].keys())
                
                with open(summary_file, 'w') as f:
                    # Write header
                    f.write('\t'.join(headers) + '\n')
                    
                    # Write data
                    for hit in hits:
                        row = [str(hit.get(header, '')) for header in headers]
                        f.write('\t'.join(row) + '\n')
                
                self.logger.info("✓ Created %s summary: %s (%d hits)", db, summary_file, len(hits))
                
                # Create HTML summary report for this database
                self._create_database_summary_html(db, hits, output_base)
            else:
                self.logger.info("No hits for database %s, skipping summary", db)
    
    def _create_database_summary_html(self, database: str, hits: List[Dict], output_base: str):
        """Create HTML summary report for a specific database across all genomes"""
        
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
        
        # Count unique genomes
        unique_genomes = list(set(hit['genome'] for hit in hits))
        
        # Count genes per genome
        genes_per_genome = {}
        for hit in hits:
            genome = hit['genome']
            if genome not in genes_per_genome:
                genes_per_genome[genome] = set()
            genes_per_genome[genome].add(hit['gene'])
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>EcoliTyper ABRicate - {database.upper()} Database Summary</title>
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
        .present {{ background-color: #d4edda; }}
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 EcoliTyper ABRicate - {database.upper()} Database Summary</h1>
            <p style="color: #666; font-size: 1.2em;">Cross-genome analysis of {database.upper()} database results</p>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
        
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">📊 Database Overview</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Hits</h3>
                    <p style="font-size: 2em; margin: 0;">{len(hits)}</p>
                </div>
                <div class="stat-card">
                    <h3>Genomes</h3>
                    <p style="font-size: 2em; margin: 0;">{len(unique_genomes)}</p>
                </div>
                <div class="stat-card">
                    <h3>Unique Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{len(set(hit['gene'] for hit in hits))}</p>
                </div>
            </div>
            <p><strong>Database:</strong> {database.upper()}</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
        </div>
        
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
        
        for genome in sorted(unique_genomes):
            genes = genes_per_genome.get(genome, set())
            gene_list = ", ".join(sorted(genes))
            html_content += f"""
                    <tr class="present">
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
                        <th>Genomes</th>
                    </tr>
                </thead>
                <tbody>
"""
        
        # Calculate gene frequency
        gene_frequency = {}
        for hit in hits:
            gene = hit['gene']
            if gene not in gene_frequency:
                gene_frequency[gene] = set()
            gene_frequency[gene].add(hit['genome'])
        
        for gene, genomes in sorted(gene_frequency.items(), key=lambda x: len(x[1]), reverse=True):
            genome_list = ", ".join(sorted(genomes))
            html_content += f"""
                    <tr>
                        <td><strong>{gene}</strong></td>
                        <td>{len(genomes)}</td>
                        <td>{genome_list}</td>
                    </tr>
"""
        
        html_content += """
                </tbody>
            </table>
        </div>
        
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #667eea; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
                Analysis performed using EcoliTyper ABRicate v1.0.1
            </p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write database summary HTML report
        html_file = os.path.join(output_base, f"ecoli_{database}_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info("Database summary HTML report: %s", html_file)
    
    def process_single_genome(self, genome_file: str, output_base: str = "ecoli_abricate_results") -> Dict[str, Any]:
        """Process a single E. coli genome with all databases and HTML reporting"""
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        
        self.logger.info("=== PROCESSING E. COLI GENOME: %s ===", genome_name)
        
        # Create output directory
        os.makedirs(results_dir, exist_ok=True)
        
        databases = self.required_databases
        
        # Run ABRicate on all databases
        results = {}
        for db in databases:
            result = self.run_abricate_single_db(genome_file, db, results_dir)
            results[db] = result
            status_icon = "✓" if result['status'] == 'success' else "✗"
            self.logger.info("%s %s: %d hits", status_icon, db, result['hit_count'])
        
        # Create comprehensive HTML report
        self.create_comprehensive_html_report(genome_name, results, results_dir)
        
        return {
            'genome': genome_name,
            'results': results,
            'total_hits': sum(r['hit_count'] for r in results.values())
        }
    
    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "ecoli_abricate_results") -> Dict[str, Any]:
        """Process multiple E. coli genomes using wildcard pattern - MAXIMUM SPEED"""
        
        # Check ABRicate installation
        if not self.check_abricate_installed():
            raise RuntimeError("ABRicate not properly installed")
        
        # Setup databases
        self.setup_abricate_databases()
        
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
        
        # Process genomes with parallel execution - MAXIMUM SPEED CONFIGURATION
        all_results = {}
        
        if len(genome_files) > 1 and self.cpus > 1:
            # Use ThreadPoolExecutor for parallel processing of multiple genomes
            self.logger.info("Using parallel processing with %d CPU cores (MAXIMUM SPEED)", self.cpus)
            
            with ThreadPoolExecutor(max_workers=self.cpus) as executor:
                # Submit all genomes for processing
                future_to_genome = {
                    executor.submit(self.process_single_genome, genome, output_base): genome 
                    for genome in genome_files
                }
                
                # Collect results as they complete
                for future in as_completed(future_to_genome):
                    genome = future_to_genome[future]
                    try:
                        result = future.result()
                        all_results[Path(genome).stem] = result
                        self.logger.info("✓ Completed: %s (%d total hits)", result['genome'], result['total_hits'])
                    except Exception as e:
                        self.logger.error("✗ Failed: %s - %s", genome, e)
        else:
            # Process genomes sequentially
            for genome in genome_files:
                try:
                    result = self.process_single_genome(genome, output_base)
                    all_results[Path(genome).stem] = result
                    self.logger.info("✓ Completed: %s (%d total hits)", result['genome'], result['total_hits'])
                except Exception as e:
                    self.logger.error("✗ Failed: %s - %s", genome, e)
        
        # Create database summary files and HTML reports after processing all genomes
        self.create_database_summaries(all_results, output_base)
        
        self.logger.info("=== E. COLI ANALYSIS COMPLETE ===")
        self.logger.info("Processed %d genomes", len(all_results))
        self.logger.info("Results saved to: %s", output_base)
        
        return all_results


def main():
    """Command line interface for E. coli ABRicate analysis"""
    parser = argparse.ArgumentParser(
        description='EcoliTyper ABRicate Analysis - MAXIMUM SPEED VERSION',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on all E. coli FASTA files (auto-detect optimal CPU cores - MAXIMUM SPEED)
  python ecoli_abricate.py "*.fna"
  
  # Run on specific pattern with auto CPU detection
  python ecoli_abricate.py "ECOLI_*.fasta"
  
  # Run with custom output directory and auto CPUs
  python ecoli_abricate.py "*.fa" --output my_ecoli_results

  # Force specific number of CPU cores
  python ecoli_abricate.py "*.fna" --cpus 4

MAXIMUM SPEED RESOURCE MANAGEMENT:
  • 1-4 cores: Uses ALL CPU cores (100% utilization)
  • 5-8 cores: Uses (cores-1) for optimal performance  
  • 9-16 cores: Uses (cores-2) for high performance
  • 17-32 cores: Uses (cores-4) for maximum throughput
  • 32+ cores: Uses 85% of cores (capped at 32)

Supported FASTA extensions: .fasta, .fa, .fna, .faa
        """
    )
    
    parser.add_argument('pattern', help='File pattern for E. coli genomes (e.g., "*.fasta", "genomes/*.fna")')
    parser.add_argument('--cpus', '-c', type=int, default=None, 
                       help='Number of CPU cores to use (default: auto-detect optimal for MAXIMUM SPEED)')
    parser.add_argument('--output', '-o', default='ecoli_abricate_results', 
                       help='Output directory (default: ecoli_abricate_results)')
    
    args = parser.parse_args()
    
    executor = EcoliAbricateExecutor(cpus=args.cpus)
    
    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)
        
        # Print summary
        executor.logger.info("\n" + "="*50)
        executor.logger.info("🧬 EcoliTyper ABRicate FINAL SUMMARY")
        executor.logger.info("="*50)
        
        total_critical_resistance = 0
        total_critical_virulence = 0
        
        for genome_name, result in results.items():
            # Collect all hits for this genome
            all_genome_hits = []
            for db_result in result['results'].values():
                all_genome_hits.extend(db_result['hits'])
            
            # Analyze for this genome
            analysis = executor.analyze_ecoli_resistance(all_genome_hits)
            
            executor.logger.info("✓ %s: %d total hits, %d critical resistance, %d critical virulence", 
                               genome_name, result['total_hits'], analysis['total_critical_resistance'], analysis['total_critical_virulence'])
            
            total_critical_resistance += analysis['total_critical_resistance']
            total_critical_virulence += analysis['total_critical_virulence']
        
        # Database usage summary
        executor.logger.info("\n" + "="*50)
        executor.logger.info("🗃️  DATABASE USAGE SUMMARY")
        executor.logger.info("="*50)
        executor.logger.info("Used databases: %s", ", ".join(executor.required_databases))
        
        # Performance summary
        executor.logger.info("\n" + "="*50)
        executor.logger.info("⚡ MAXIMUM SPEED PERFORMANCE SUMMARY")
        executor.logger.info("="*50)
        executor.logger.info("CPU cores utilized: %d cores", executor.cpus)
        executor.logger.info("Available RAM: %.1f GB", executor.available_ram)
        executor.logger.info("Total E. coli genomes processed: %d", len(results))
        executor.logger.info("Total critical resistance genes found: %d", total_critical_resistance)
        executor.logger.info("Total critical virulence genes found: %d", total_critical_virulence)
        executor.logger.info("Processing mode: MAXIMUM SPEED 🚀")
        
        import random
        executor.logger.info("\n💡 %s", random.choice(executor.science_quotes))
        
    except Exception as e:
        executor.logger.error("E. coli analysis failed: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
