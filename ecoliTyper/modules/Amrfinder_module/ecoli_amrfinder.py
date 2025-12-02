#!/usr/bin/env python3
"""
EcoliTyper AMRfinderPlus - E. coli AMR Analysis
Comprehensive AMR analysis for E. coli with beautiful HTML reporting - MAXIMUM SPEED VERSION
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

class EcoliAMRfinderPlus:
    """AMRfinderPlus executor for E. coli with comprehensive HTML reporting - MAXIMUM SPEED"""
    
    def __init__(self, cpus: int = None):
        # Setup logging FIRST
        self.logger = self._setup_logging()
        
        # Initialize available_ram before calculating cpus
        self.available_ram = self._get_available_ram()
        
        # Then calculate resources - MAXIMUM SPEED MODE
        self.cpus = self._calculate_optimal_cpus(cpus)
        
        self.metadata = {
            "tool_name": "EcoliTyper AMRfinderPlus",
            "version": "1.0.0", 
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
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
        """Log resource allocation information - KEEPING EcoliTyper STYLING"""
        self.logger.info(f"Available RAM: {self.available_ram:.1f} GB")
        
        if total_cores:
            self.logger.info(f"System CPU cores: {total_cores}")
            utilization = (cpus / total_cores) * 100
            self.logger.info(f"Using CPU cores: {cpus} ({utilization:.1f}% of available cores)")
        else:
            self.logger.info(f"Using user-specified CPU cores: {cpus}")
        
        # Performance recommendations - MAXIMUM SPEED FOCUS (EcoliTyper style)
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

    def run_amrfinder_single_genome(self, genome_file: str, output_dir: str) -> Dict[str, Any]:
        """Run AMRfinderPlus on a single E. coli genome - MAXIMUM SPEED WITH ALL CORES"""
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"{genome_name}_amrfinder.txt")
        
        # AMRfinderPlus uses THREADS - allocate ALL available cores for maximum speed
        # Since we're processing one sample at a time, we can use all system resources
        run_threads = self.cpus
        
        cmd = [
            'amrfinder',
            '--nucleotide', genome_file,
            '--output', output_file,
            '--threads', str(run_threads),  # Use ALL cores for this single sample
            '--plus'
        ]
        
        self.logger.info("🚀 MAXIMUM SPEED: Running AMRfinderPlus on %s (using ALL %d CORES)", genome_name, run_threads)
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse results for reporting
            hits = self._parse_amrfinder_output(output_file)
            
            # Create individual HTML report
            self._create_amrfinder_html_report(genome_name, hits, output_dir)
            
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': hits,
                'hit_count': len(hits),
                'status': 'success'
            }
            
        except subprocess.CalledProcessError as e:
            self.logger.error("AMRfinderPlus failed for %s: %s", genome_name, e.stderr)
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'status': 'failed'
            }
    
    def _parse_amrfinder_output(self, amrfinder_file: str) -> List[Dict]:
        """Parse AMRfinderPlus output file into structured data"""
        hits = []
        try:
            with open(amrfinder_file, 'r') as f:
                lines = f.readlines()
                
            if not lines or len(lines) < 2:
                return hits
                
            # Parse header
            headers = lines[0].strip().split('\t')
            
            # Parse data lines
            for line_num, line in enumerate(lines[1:], 2):
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split('\t')
                if len(parts) >= len(headers):
                    hit = {}
                    for i, header in enumerate(headers):
                        if i < len(parts):
                            hit[header] = parts[i]
                        else:
                            hit[header] = ''
                    
                    # Map to consistent field names
                    processed_hit = {
                        'protein_id': hit.get('Protein identifier', ''),
                        'contig_id': hit.get('Contig id', ''),
                        'start': hit.get('Start', ''),
                        'stop': hit.get('Stop', ''),
                        'strand': hit.get('Strand', ''),
                        'gene_symbol': hit.get('Gene symbol', ''),
                        'sequence_name': hit.get('Sequence name', ''),
                        'scope': hit.get('Scope', ''),
                        'element_type': hit.get('Element type', ''),
                        'element_subtype': hit.get('Element subtype', ''),
                        'class': hit.get('Class', ''),
                        'subclass': hit.get('Subclass', ''),
                        'method': hit.get('Method', ''),
                        'target_length': hit.get('Target length', ''),
                        'ref_length': hit.get('Reference sequence length', ''),
                        'coverage': hit.get('% Coverage of reference sequence', ''),
                        'identity': hit.get('% Identity to reference sequence', ''),
                        'alignment_length': hit.get('Alignment length', ''),
                        'accession': hit.get('Accession of closest sequence', ''),
                        'closest_name': hit.get('Name of closest sequence', ''),
                        'hmm_id': hit.get('HMM id', ''),
                        'hmm_description': hit.get('HMM description', '')
                    }
                    hits.append(processed_hit)
                else:
                    self.logger.warning("Line %d has %d parts, expected %d: %s", 
                                      line_num, len(parts), len(headers), line[:100] + "...")
                    
        except Exception as e:
            self.logger.error("Error parsing %s: %s", amrfinder_file, e)
            
        self.logger.info("Parsed %d AMR hits from %s", len(hits), amrfinder_file)
        return hits
    
    def _create_amrfinder_html_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        """Create comprehensive HTML report for AMRfinderPlus results with beautiful styling"""
        
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
                if hit['gene_symbol'] in analysis['critical_risk_list']:
                    row_class = "critical-row"
                elif hit['gene_symbol'] in analysis['high_risk_list']:
                    row_class = "high-risk-row"
                
                html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{hit['gene_symbol']}</strong></td>
                        <td title="{hit['sequence_name']}">{hit['sequence_name'][:80]}{'...' if len(hit['sequence_name']) > 80 else ''}</td>
                        <td>{hit['class']}</td>
                        <td>{hit['subclass']}</td>
                        <td>{hit['coverage']}%</td>
                        <td>{hit['identity']}%</td>
                        <td>{hit['scope']}</td>
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
                Analysis performed using EcoliTyper AMRfinderPlus v3.12.8
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
            gene_symbol = hit['gene_symbol']
            resistance_class = hit['class']
            
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
            # Write header
            f.write("Genome\tGene_Symbol\tSequence_Name\tClass\tSubclass\tCoverage\tIdentity\tScope\tElement_Type\tAccession\tContig\tStart\tStop\n")
            
            # Write data for all genomes
            for genome_name, result in all_results.items():
                for hit in result['hits']:
                    row = [
                        genome_name,
                        hit.get('gene_symbol', ''),
                        hit.get('sequence_name', ''),
                        hit.get('class', ''),
                        hit.get('subclass', ''),
                        hit.get('coverage', ''),
                        hit.get('identity', ''),
                        hit.get('scope', ''),
                        hit.get('element_type', ''),
                        hit.get('accession', ''),
                        hit.get('contig_id', ''),
                        hit.get('start', ''),
                        hit.get('stop', '')
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
        
        # Create comprehensive HTML summary report for ecoli_amrfinder_summary.tsv
        self._create_summary_html_report(all_results, output_base)
    
    def _create_summary_html_report(self, all_results: Dict[str, Any], output_base: str):
        """Create comprehensive HTML summary report with pattern discovery"""
        
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
            padding: 12px; 
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
        .resistance-badge {{
            display: inline-block;
            background: #dc3545;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .critical-resistance-badge {{
            display: inline-block;
            background: #8b0000;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
            font-weight: bold;
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
        .success-badge {{
            display: inline-block;
            background: #28a745;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        /* IMPROVED GENE FREQUENCY COLOR SCHEME */
        .frequency-high {{ background-color: #f8d7da; font-weight: bold; border-left: 4px solid #dc3545; }}
        .frequency-medium-high {{ background-color: #ffeaa7; border-left: 4px solid #fdcb6e; }}
        .frequency-medium {{ background-color: #fff3cd; border-left: 4px solid #ffc107; }}
        .frequency-low-medium {{ background-color: #d1ecf1; border-left: 4px solid #17a2b8; }}
        .frequency-low {{ background-color: #d4edda; border-left: 4px solid #28a745; }}
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 EcoliTyper AMRfinderPlus - Summary Report</h1>
            <p style="color: #666; font-size: 1.2em;">Comprehensive E. coli Antimicrobial Resistance Analysis Across All Genomes</p>
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
                html_content += f'<span class="critical-resistance-badge">🚨 {gene}</span>'
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
                    <h3>Critical Risk Genomes</h3>
                    <p style="font-size: 2em; margin: 0;">{genomes_with_critical}</p>
                </div>
            </div>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
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
                html_content += f'<span class="resistance-badge">{gene}</span>'
            html_content += """
            </div>
        </div>
"""
        
        # Genes by Genome table (Pattern Discovery)
        html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🔍 Genes by Genome</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Genome</th>
                        <th>Gene Count</th>
                        <th>Critical Genes</th>
                        <th>High Risk Genes</th>
                    </tr>
                </thead>
                <tbody>
"""
        
        for genome in sorted(genes_per_genome.keys()):
            genes = genes_per_genome.get(genome, set())
            critical_genes = [g for g in genes if g in self.critical_risk_genes]
            high_risk_genes = [g for g in genes if g in self.high_risk_genes and g not in self.critical_risk_genes]
            
            critical_display = ", ".join(critical_genes) if critical_genes else "None"
            high_risk_display = ", ".join(high_risk_genes) if high_risk_genes else "None"
            
            # Highlight rows with critical genes
            row_class = "critical-row" if critical_genes else "high-risk-row" if high_risk_genes else ""
            
            html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{genome}</strong></td>
                        <td>{len(genes)}</td>
                        <td>{critical_display}</td>
                        <td>{high_risk_display}</td>
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
        
        # Calculate gene frequency with IMPROVED color highlighting
        for gene, genomes in sorted(gene_frequency.items(), key=lambda x: len(x[1]), reverse=True):
            frequency = len(genomes)
            genome_list = ", ".join(sorted(genomes))
            frequency_percent = (frequency / total_genomes) * 100
            
            # Determine risk level
            if gene in self.critical_risk_genes:
                risk_level = '<span class="critical-resistance-badge">CRITICAL</span>'
            elif gene in self.high_risk_genes:
                risk_level = '<span class="resistance-badge">HIGH</span>'
            else:
                risk_level = '<span class="success-badge">Standard</span>'
            
            # IMPROVED COLOR SCHEME: Better visual distinction
            if frequency_percent >= 75:
                frequency_class = "frequency-high"
                prevalence_badge = '<span class="resistance-badge">Very High</span>'
            elif frequency_percent >= 50:
                frequency_class = "frequency-medium-high"
                prevalence_badge = '<span class="warning-badge">High</span>'
            elif frequency_percent >= 25:
                frequency_class = "frequency-medium"
                prevalence_badge = '<span class="warning-badge">Medium</span>'
            elif frequency_percent >= 10:
                frequency_class = "frequency-low-medium"
                prevalence_badge = '<span class="success-badge">Low</span>'
            else:
                frequency_class = "frequency-low"
                prevalence_badge = '<span class="success-badge">Rare</span>'
            
            html_content += f"""
                    <tr class="{frequency_class}">
                        <td><strong>{gene}</strong></td>
                        <td>{frequency} ({frequency_percent:.1f}%)</td>
                        <td>{prevalence_badge}</td>
                        <td>{risk_level}</td>
                        <td>{genome_list}</td>
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
                <li><strong>Individual genome HTML reports</strong> - Detailed analysis per genome</li>
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
                Analysis performed using EcoliTyper AMRfinderPlus v3.12.8
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
    
    def process_single_genome(self, genome_file: str, output_base: str = "ecoli_amrfinder_results") -> Dict[str, Any]:
        """Process a single E. coli genome with AMRfinderPlus"""
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        
        self.logger.info("=== PROCESSING E. COLI GENOME: %s ===", genome_name)
        
        # Create output directory
        os.makedirs(results_dir, exist_ok=True)
        
        # Run AMRfinderPlus
        result = self.run_amrfinder_single_genome(genome_file, results_dir)
        
        status_icon = "✓" if result['status'] == 'success' else "✗"
        self.logger.info("%s %s: %d AMR hits", status_icon, genome_name, result['hit_count'])
        
        return result
    
    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "ecoli_amrfinder_results") -> Dict[str, Any]:
        """Process multiple E. coli genomes using wildcard pattern - MAXIMUM SPEED"""
        
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
        max_concurrent = max(1, min(self.cpus, len(genome_files), int(self.available_ram / 2.5)))  # 2.5GB per genome
        
        self.logger.info("🚀 MAXIMUM SPEED: Using %d concurrent genome processing jobs", max_concurrent)
        self.logger.info("   Each AMRfinderPlus instance uses %d threads internally", self.cpus)
        self.logger.info("   This provides maximum throughput for multiple genome analysis")
        
        with ThreadPoolExecutor(max_workers=max_concurrent) as executor:
            # Submit all tasks
            future_to_genome = {
                executor.submit(self.process_single_genome, genome, output_base): genome 
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
                        'status': 'failed'
                    }
        
        # Create AMR summary files and HTML reports after processing all genomes
        self.create_amr_summary(all_results, output_base)
        
        self.logger.info("=== E. COLI AMR ANALYSIS COMPLETE ===")
        self.logger.info("Processed %d genomes", len(all_results))
        self.logger.info("Results saved to: %s", output_base)
        
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
  
  # Run on specific pattern with auto CPU detection
  python ecoli_amrfinder.py "ECOLI_*.fasta"
  
  # Force specific number of CPU cores
  python ecoli_amrfinder.py "*.fa" --cpus 16

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
    parser.add_argument('--output', '-o', default='ecoli_amrfinder_results', 
                       help='Output directory (default: ecoli_amrfinder_results)')
    
    args = parser.parse_args()
    
    executor = EcoliAMRfinderPlus(cpus=args.cpus)
    
    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)
        
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
        executor.logger.info("   Summary HTML report: %s/ecoli_amrfinder_summary_report.html", args.output)
        
        # Performance summary
        executor.logger.info("\n⚡ MAXIMUM SPEED PERFORMANCE SUMMARY:")
        executor.logger.info("   CPU cores utilized: %d cores", executor.cpus)
        executor.logger.info("   Available RAM: %.1f GB", executor.available_ram)
        executor.logger.info("   Processing mode: MAXIMUM SPEED CONCURRENT MODE 🚀")
        executor.logger.info("   Strategy: Process multiple genomes concurrently with optimal core allocation")
        
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
