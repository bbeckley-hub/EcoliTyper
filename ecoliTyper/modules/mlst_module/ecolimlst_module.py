#!/usr/bin/env python3
"""
EcoliTyper Module for MLST Analysis of E. coli
Comprehensive MLST analysis with beautiful HTML reporting
Author: Beckley Brown <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
Date: 2025
Send a quick mail for any issues or further explanations.
"""

import os
import sys
import json
import glob
import argparse
import subprocess
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
from datetime import datetime

class EcoliTyper:
    def __init__(self, database_dir: Path, script_dir: Path):
        self.database_dir = database_dir
        self.script_dir = script_dir
        self.mlst_bin = script_dir / "mlst"
        
        self.metadata = {
            "tool_name": "EcoliTyper MLST Analysis",
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
        
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files using glob patterns"""
        if os.path.isfile(input_path):
            return [Path(input_path)]
        
        fasta_patterns = [
            input_path,
            f"{input_path}/*.fna", f"{input_path}/*.fasta",
            f"{input_path}/*.fa", f"{input_path}/*.fn",
            f"{input_path}/*.fna.gz", f"{input_path}/*.fasta.gz",
            f"{input_path}/*.fa.gz", f"{input_path}/*.gb",
            f"{input_path}/*.gbk", f"{input_path}/*.gbff"
        ]
        
        fasta_files = []
        for pattern in fasta_patterns:
            matched_files = glob.glob(pattern)
            for file_path in matched_files:
                path = Path(file_path)
                if path.is_file():
                    fasta_files.append(path)
        
        return sorted(list(set(fasta_files)))

    def run_mlst_single(self, input_file: Path, output_dir: Path, scheme: str = "ecoli_achtman_4") -> Dict:
        """Run MLST analysis for a single file"""
        print(f"🧬 Processing: {input_file.name}")
        
        # Create sample-specific output directory
        sample_output_dir = output_dir / input_file.stem
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save raw MLST output first
        raw_output_file = sample_output_dir / "mlst_raw_output.txt"
        
        # Run MLST command
        mlst_cmd = [
            "perl", str(self.mlst_bin),
            str(input_file),
            "--scheme", scheme,
            "--csv",
            "--nopath"
        ]
        
        try:
            # Run and capture output
            result = subprocess.run(mlst_cmd, capture_output=True, text=True, check=True)
            
            # Save raw output
            with open(raw_output_file, 'w') as f:
                f.write("STDOUT:\n")
                f.write(result.stdout)
                f.write("\nSTDERR:\n")
                f.write(result.stderr)
            
            # Parse the CSV output (it's comma-separated!)
            mlst_results = self.parse_mlst_csv(result.stdout, input_file.name)
            
            # Generate only 3 output files
            self.generate_output_files(mlst_results, sample_output_dir)
            
            print(f"✅ Completed: {input_file.name} -> ST{mlst_results.get('st', 'ND')}")
            return mlst_results
            
        except subprocess.CalledProcessError as e:
            print(f"❌ MLST failed for {input_file.name}")
            error_result = self.get_fallback_results(input_file.name)
            self.generate_output_files(error_result, sample_output_dir)
            return error_result

    def parse_mlst_csv(self, stdout: str, sample_name: str) -> Dict:
        """Parse MLST CSV output - it's comma-separated!"""
        lines = stdout.strip().split('\n')
        if not lines:
            return self.get_empty_results(sample_name)
        
        # Find the result line (usually the last line with data)
        result_line = None
        for line in reversed(lines):
            if line.strip() and ',' in line and not line.startswith('['):
                result_line = line.strip()
                break
        
        if not result_line:
            return self.get_empty_results(sample_name)
        
        # Split by COMMA, not tab!
        parts = result_line.split(',')
        
        if len(parts) < 3:
            return self.get_empty_results(sample_name)
        
        # Extract components - format: filename,scheme,ST,allele1,allele2,...
        filename = parts[0]
        scheme = parts[1]
        st = parts[2]
        
        # Extract alleles from remaining parts
        alleles = {}
        allele_parts = []
        
        for i in range(3, len(parts)):
            allele_str = parts[i]
            if '(' in allele_str and ')' in allele_str:
                # Format: arcC(1)
                gene = allele_str.split('(')[0]
                allele = allele_str.split('(')[1].rstrip(')')
                alleles[gene] = allele
                allele_parts.append(f"{gene}({allele})")
        
        allele_profile = '-'.join(allele_parts) if allele_parts else ""
        
        return {
            "sample": sample_name,
            "st": st,
            "scheme": scheme,
            "alleles": alleles,
            "allele_profile": allele_profile,
            "confidence": "HIGH" if st and st != '-' and st != 'ND' else "LOW"
        }

    def get_empty_results(self, sample_name: str) -> Dict:
        """Return empty results structure"""
        return {
            "sample": sample_name,
            "st": "ND",
            "scheme": "ecoli_achtman_4",
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW"
        }

    def get_fallback_results(self, sample_name: str) -> Dict:
        """Fallback when MLST fails"""
        return {
            "sample": sample_name,
            "st": "UNKNOWN",
            "scheme": "ecoli_achtman_4",
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW",
            "error": "MLST analysis failed"
        }

    def generate_output_files(self, mlst_results: Dict, output_dir: Path):
        """Generate only 3 output files: HTML, TXT, and TSV"""
        # 1. Beautiful HTML Report
        self.generate_html_report(mlst_results, output_dir)
        
        # 2. Detailed Text Report
        self.generate_text_report(mlst_results, output_dir)
        
        # 3. Simple TSV Report
        self.generate_tsv_report(mlst_results, output_dir)

    def generate_text_report(self, mlst_results: Dict, output_dir: Path):
        """Generate detailed text report"""
        report = f"""MLST Analysis Report
===================

Sample: {mlst_results['sample']}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

MLST TYPING RESULTS:
-------------------
Sequence Type (ST): {mlst_results['st']}
Scheme: {mlst_results['scheme']}
Confidence: {mlst_results['confidence']}

Allele Profile:
{mlst_results['allele_profile']}

Detailed Alleles:
"""
        for gene, allele in mlst_results['alleles'].items():
            report += f"- {gene}: {allele}\n"
        
        with open(output_dir / "mlst_report.txt", 'w') as f:
            f.write(report)

    def generate_tsv_report(self, mlst_results: Dict, output_dir: Path):
        """Generate simple TSV report"""
        tsv_content = f"Sample\tST\tScheme\tAllele_Profile\tConfidence\n"
        tsv_content += f"{mlst_results['sample']}\t{mlst_results['st']}\t{mlst_results['scheme']}\t{mlst_results['allele_profile']}\t{mlst_results['confidence']}\n"
        
        with open(output_dir / "mlst_report.tsv", 'w') as f:
            f.write(tsv_content)

    def generate_html_report(self, mlst_results: Dict, output_dir: Path):
        """Generate HTML report with beautiful purple styling"""
        
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
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliTyper - MLST Analysis Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        
        .container {{
            max-width: 1200px;
            margin: 0 auto;
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
            color: #333;
            padding: 25px;
            border-radius: 12px;
            margin-bottom: 20px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.1);
            backdrop-filter: blur(10px);
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
        
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-top: 15px;
        }}
        
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
        }}
        
        .allele-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(150px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }}
        
        .allele-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
            text-align: center;
            font-weight: bold;
        }}
        
        .confidence-high {{
            color: #28a745;
            font-weight: bold;
        }}
        
        .confidence-low {{
            color: #dc3545;
            font-weight: bold;
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
        
        h1, h2, h3 {{
            color: #333;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        
        .st-badge {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 10px 20px;
            border-radius: 25px;
            font-size: 1.5em;
            font-weight: bold;
            display: inline-block;
            margin: 10px 0;
        }}
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 EcoliTyper MLST Analysis Report</h1>
            <p style="color: #666; font-size: 1.2em;">Comprehensive MLST typing for Escherichia coli</p>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
        
        <div class="card">
            <h2>📊 Sample Information</h2>
            <div class="metrics-grid">
                <div class="stat-card">
                    <div style="font-size: 14px; opacity: 0.9;">Sample Name</div>
                    <div style="font-size: 1.5em; font-weight: bold;">{mlst_results['sample']}</div>
                </div>
                <div class="stat-card">
                    <div style="font-size: 14px; opacity: 0.9;">Analysis Date</div>
                    <div style="font-size: 1.2em; font-weight: bold;">{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
                </div>
                <div class="stat-card">
                    <div style="font-size: 14px; opacity: 0.9;">MLST Scheme</div>
                    <div style="font-size: 1.2em; font-weight: bold;">{mlst_results['scheme'].title()}</div>
                </div>
            </div>
        </div>
        
        <div class="card">
            <h2>🔬 MLST Typing Results</h2>
            <div class="metrics-grid">
                <div class="stat-card">
                    <div style="font-size: 14px; opacity: 0.9;">Sequence Type</div>
                    <div class="st-badge">ST{mlst_results['st']}</div>
                </div>
                <div class="stat-card">
                    <div style="font-size: 14px; opacity: 0.9;">Confidence</div>
                    <div style="font-size: 1.5em; font-weight: bold;" class="confidence-{mlst_results['confidence'].lower()}">{mlst_results['confidence']}</div>
                </div>
            </div>
            
            <h3>🧬 Allele Profile</h3>
            <div style="background: #f8f9fa; padding: 15px; border-radius: 8px; margin: 15px 0; border-left: 4px solid #667eea;">
                <code style="font-size: 16px; color: #667eea; font-weight: bold;">{mlst_results['allele_profile']}</code>
            </div>
            
            <h3>🔍 Individual Alleles</h3>
            <div class="allele-grid">
'''
        
        # Add allele cards
        for gene, allele in mlst_results['alleles'].items():
            html_content += f'''                <div class="allele-card">
                    <div style="font-size: 12px; opacity: 0.9;">{gene}</div>
                    <div style="font-size: 18px;">{allele}</div>
                </div>
'''
        
        html_content += f'''            </div>
        </div>
        
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #667eea; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
                Analysis performed using EcoliTyper MLST(ecoli_achtman_4) Scheme
            </p>
        </div>
    </div>
</body>
</html>'''
        
        with open(output_dir / "mlst_report.html", 'w', encoding='utf-8') as f:
            f.write(html_content)

    def create_mlst_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create comprehensive MLST summary files for all samples"""
        print("📊 Creating MLST summary files...")
        
        # Create TSV summary
        self.create_mlst_tsv_summary(all_results, output_dir)
        
        # Create HTML summary
        self.create_mlst_html_summary(all_results, output_dir)
        
        print("✅ MLST summary files created successfully!")

    def create_mlst_tsv_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create TSV summary file with all samples"""
        summary_file = output_dir / "mlst_summary.tsv"
        
        with open(summary_file, 'w') as f:
            # Write header
            f.write("Sample\tST\tAllele_Profile\t")
            
            # Get all unique gene names from all samples
            all_genes = set()
            for result in all_results.values():
                all_genes.update(result['alleles'].keys())
            
            # Write gene headers
            for gene in sorted(all_genes):
                f.write(f"{gene}\t")
            f.write("\n")
            
            # Write data for each sample
            for sample_name, result in all_results.items():
                f.write(f"{sample_name}\t{result['st']}\t{result['allele_profile']}\t")
                
                # Write allele values for each gene
                for gene in sorted(all_genes):
                    allele = result['alleles'].get(gene, '')
                    f.write(f"{allele}\t")
                f.write("\n")
        
        print(f"📄 TSV summary created: {summary_file}")

    def create_mlst_html_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create HTML summary with beautiful styling"""
        summary_file = output_dir / "mlst_summary.html"
        
        # Get all unique gene names
        all_genes = set()
        for result in all_results.values():
            all_genes.update(result['alleles'].keys())
        sorted_genes = sorted(all_genes)
        
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
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliTyper - MLST Summary Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 0;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}
        
        .container {{
            max-width: 1800px;
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
            color: #333;
            padding: 25px;
            border-radius: 12px;
            margin-bottom: 20px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.1);
            backdrop-filter: blur(10px);
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
        
        .summary-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        
        .summary-table th {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: 600;
        }}
        
        .summary-table td {{
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
        }}
        
        .summary-table tr:nth-child(even) {{
            background-color: #f8fafc;
        }}
        
        .summary-table tr:hover {{
            background-color: #f1f5f9;
        }}
        
        .st-cell {{
            font-weight: bold;
            color: #667eea;
        }}
        
        .stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 20px;
        }}
        
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
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
        
        h1, h2 {{
            color: #333;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 EcoliTyper - MLST Summary Report</h1>
            <p style="color: #666; font-size: 1.2em;">Comprehensive MLST analysis summary for all samples</p>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
        
        <div class="card">
            <h2>📊 Summary Statistics</h2>
            <div class="stats">
                <div class="stat-card">
                    <div style="font-size: 2em; font-weight: bold;">{len(all_results)}</div>
                    <div>Samples Processed</div>
                </div>
                <div class="stat-card">
                    <div style="font-size: 2em; font-weight: bold;">{len(set(result['st'] for result in all_results.values()))}</div>
                    <div>Unique STs</div>
                </div>
                <div class="stat-card">
                    <div style="font-size: 2em; font-weight: bold;">{len(sorted_genes)}</div>
                    <div>MLST Genes</div>
                </div>
            </div>
        </div>
        
        <div class="card">
            <h2>🔍 MLST Results</h2>
            <div style="overflow-x: auto;">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>ST</th>
                            <th>Allele Profile</th>
'''
        
        # Add gene headers
        for gene in sorted_genes:
            html_content += f'                            <th>{gene}</th>\n'
        
        html_content += '''                        </tr>
                    </thead>
                    <tbody>
'''
        
        # Add data rows
        for sample_name, result in all_results.items():
            html_content += f'''                        <tr>
                            <td><strong>{sample_name}</strong></td>
                            <td class="st-cell">ST{result['st']}</td>
                            <td>{result['allele_profile']}</td>
'''
            
            # Add allele values for each gene
            for gene in sorted_genes:
                allele = result['alleles'].get(gene, '')
                html_content += f'                            <td>{allele}</td>\n'
            
            html_content += '                        </tr>\n'
        
        html_content += f'''                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #667eea; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
                Analysis performed using EcoliTyper MLST(ecoli_achtman_4) Scheme
            </p>
        </div>
    </div>
</body>
</html>'''
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"📊 HTML summary created: {summary_file}")

    def run_mlst_batch(self, input_path: str, output_dir: Path, scheme: str = "ecoli_achtman_4") -> Dict[str, Dict]:
        """Run MLST analysis for multiple files"""
        print("🔍 Searching for FASTA files...")
        fasta_files = self.find_fasta_files(input_path)
        
        if not fasta_files:
            print("❌ No FASTA files found!")
            return {}
        
        print(f"📁 Found {len(fasta_files)} FASTA files")
        
        results = {}
        for fasta_file in fasta_files:
            result = self.run_mlst_single(fasta_file, output_dir, scheme)
            results[fasta_file.name] = result
        
        # Create summary files after processing all samples
        self.create_mlst_summary(results, output_dir)
        
        return results

def main():
    parser = argparse.ArgumentParser(
        description='EcoliTyper - MLST Analyzer for E. coli',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file analysis
  python ecolimlst_module.py -i genome.fasta -o results -db /path/to/db -sc /path/to/script
  
  # Batch analysis
  python ecolimlst_module.py -i "*.fasta" -o results -db /path/to/db -sc /path/to/script --batch

Author: Brown Beckley <brownbeckley94@gmail.com>
GitHub: https://github.com/bbeckley-hub
Affiliation: University of Ghana Medical School
        """
    )
    
    parser.add_argument('-i', '--input', required=True, 
                       help='Input FASTA file or directory (supports wildcards)')
    parser.add_argument('-o', '--output-dir', required=True, 
                       help='Output directory')
    parser.add_argument('-db', '--database-dir', required=True,
                       help='Database directory')
    parser.add_argument('-sc', '--script-dir', required=True,
                       help='Script directory (contains mlst binary)')
    parser.add_argument('-s', '--scheme', default='ecoli_achtman_4',
                       help='MLST scheme (default: ecoli_achtman_4)')
    parser.add_argument('--batch', action='store_true',
                       help='Process multiple files')
    
    args = parser.parse_args()
    
    analyzer = EcoliTyper(
        database_dir=Path(args.database_dir),
        script_dir=Path(args.script_dir)
    )
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("🧬 EcoliTyper MLST Analysis")
    print("=" * 50)
    print(f"Input: {args.input}")
    print(f"Output: {output_dir}")
    print(f"Scheme: {args.scheme}")
    print("=" * 50)
    
    if args.batch:
        results = analyzer.run_mlst_batch(args.input, output_dir, args.scheme)
        print(f"✅ Batch MLST completed! Processed {len(results)} samples")
    else:
        input_file = Path(args.input)
        if input_file.exists():
            result = analyzer.run_mlst_single(input_file, output_dir, args.scheme)
            print(f"✅ MLST completed for {input_file.name}: ST{result.get('st', 'ND')}")
        else:
            print(f"❌ Input file not found: {args.input}")

if __name__ == "__main__":
    main()
