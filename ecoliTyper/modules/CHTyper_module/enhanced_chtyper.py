#!/usr/bin/env python3
"""
EcoliTyper - Enhanced CHTyper Wrapper
Comprehensive E. coli typing with FumC and FimH with auto BLAST detection
Email: <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
Date: 2025
Send a quick mail for any issues or further explanations.
"""

import os
import sys
import json
import argparse
import subprocess
import shutil
import glob
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any
import pandas as pd
import multiprocessing as mp

class EnhancedCHTyper:
    def __init__(self, db_path: str = "chtyper_db", threads: int = 4):
        self.db_path = Path(db_path)
        self.threads = threads
        self.results = []
        self.metadata = {
            "tool_name": "EcoliTyper CHTyper",
            "version": "1.0.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        
        # Auto-detect BLAST path
        self.blast_path = self._detect_blast_path()
        if not self.blast_path:
            raise RuntimeError("BLAST not found. Please install BLAST or ensure it's in your PATH")
        
        self.science_quotes = [
            "“The important thing is not to stop questioning. Curiosity has its own reason for existence.” - Albert Einstein",
            "“Nothing in life is to be feared, it is only to be understood.” - Marie Curie",
            "“The microscope opens a new world to the investigator.” - Robert Koch",
            "“In science, the credit goes to the man who convinces the world, not to the man to whom the idea first occurs.” - Francis Darwin",
            "“The good thing about science is that it's true whether or not you believe in it.” - Neil deGrasse Tyson",
            "“Science knows no country, because knowledge belongs to humanity.” - Louis Pasteur"
        ]
        
        self.ascii_art = """
 ██████╗██╗  ██╗████████╗██╗   ██╗██████╗ ███████╗██████╗ 
██╔════╝██║  ██║╚══██╔══╝╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
██║     ███████║   ██║    ╚████╔╝ ██████╔╝█████╗  ██████╔╝
██║     ██╔══██║   ██║     ╚██╔╝  ██╔═══╝ ██╔══╝  ██╔══██╗
╚██████╗██║  ██║   ██║      ██║   ██║     ███████╗██║  ██║
 ╚═════╝╚═╝  ╚═╝   ╚═╝      ╚═╝   ╚═╝     ╚══════╝╚═╝  ╚═╝
        """
    
    def _detect_blast_path(self) -> str:
        """Auto-detect BLAST installation"""
        try:
            result = subprocess.run(["which", "blastn"], capture_output=True, text=True)
            if result.returncode == 0:
                return result.stdout.strip()
        except:
            pass
        return ""
    
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files using glob patterns or direct paths"""
        fasta_files = []
        
        if '*' in input_path or '?' in input_path or '[' in input_path:
            matches = glob.glob(input_path)
            for match in matches:
                path = Path(match)
                if path.is_file() and path.suffix.lower() in ['.fasta', '.fna', '.fa', '.fsa']:
                    fasta_files.append(path)
        else:
            input_path = Path(input_path)
            if input_path.is_file():
                if input_path.suffix.lower() in ['.fasta', '.fna', '.fa', '.fsa']:
                    fasta_files = [input_path]
            elif input_path.is_dir():
                for ext in ['*.fasta', '*.fna', '*.fa', '.fsa']:
                    fasta_files.extend(input_path.glob(ext))
                    fasta_files.extend(input_path.glob(ext.upper()))
        
        if not fasta_files:
            raise ValueError(f"No FASTA files found matching: {input_path}")
            
        print(f"Found {len(fasta_files)} FASTA file(s) for analysis")
        return fasta_files
    
    def run_chtyper_analysis(self, fasta_file: Path, output_base: Path) -> Dict[str, Any]:
        """Run CHTyper on a single FASTA file"""
        try:
            sample_name = fasta_file.stem
            sample_output_dir = output_base / sample_name
            
            # Create sample-specific output directory
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            
            print(f"🔬 Analyzing {sample_name}...")
            
            # Copy the input file to the sample directory
            sample_fasta = sample_output_dir / fasta_file.name
            shutil.copy2(fasta_file, sample_fasta)
            
            # Use CHTyper-fixed.py (the working version)
            chyper_script = "CHTyper-fixed.py"
            if not Path(chyper_script).exists():
                return self._create_error_result(sample_name, str(fasta_file), "CHTyper script not found")
            
            # Build CHTyper command with auto-detected BLAST path
            cmd = [
                "python3", chyper_script,
                "-i", str(sample_fasta),
                "-o", str(sample_output_dir),
                "-p", str(self.db_path),
                "-b", self.blast_path,
                "-l", "0.6",  # Minimum coverage
                "-t", "0.9"   # Threshold
            ]
            
            # Run CHTyper
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                error_msg = f"Command failed with return code {result.returncode}"
                if result.stderr:
                    error_msg += f": {result.stderr}"
                return self._create_error_result(sample_name, str(fasta_file), error_msg)
            
            # Parse results from output files
            return self._parse_sample_results(sample_name, str(fasta_file), sample_output_dir)
            
        except Exception as e:
            return self._create_error_result(sample_name, str(fasta_file), str(e))
    
    def _parse_sample_results(self, sample_name: str, fasta_path: str, output_dir: Path) -> Dict[str, Any]:
        """Parse results from CHTyper output files - FIXED PARSING"""
        try:
            # Parse tabular results
            tab_file = output_dir / "results_tab.txt"
            if tab_file.exists():
                with open(tab_file, 'r') as f:
                    content = f.read()
                
                # Parse FumC and FimH types CORRECTLY
                fumc_type = "Unknown"
                fimh_type = "Unknown"
                fumc_identity = "0.00"
                fimh_identity = "0.00"
                fumc_coverage = "0/0"
                fimh_coverage = "0/0"
                
                lines = content.strip().split('\n')
                for line in lines:
                    if line.startswith('fumC') and '\t' in line:
                        parts = line.split('\t')
                        if len(parts) >= 6:
                            fumc_type = parts[0]  # fumC11 (the actual type)
                            fumc_identity = parts[1]  # 100.00
                            fumc_coverage = parts[2]  # 469/469
                    elif line.startswith('fimH') and '\t' in line:
                        parts = line.split('\t')
                        if len(parts) >= 6:
                            fimh_type = parts[0]  # fimH27 (the actual type)
                            fimh_identity = parts[1]  # 100.00
                            fimh_coverage = parts[2]  # 489/489
                
                # Parse detailed results for additional info
                detailed_results = ""
                txt_file = output_dir / "results.txt"
                if txt_file.exists():
                    with open(txt_file, 'r') as f:
                        detailed_results = f.read()
                
                return {
                    "sample_id": sample_name,
                    "file_path": fasta_path,
                    "fumc_type": fumc_type,
                    "fimh_type": fimh_type,
                    "fumc_identity": fumc_identity,
                    "fimh_identity": fimh_identity,
                    "fumc_coverage": fumc_coverage,
                    "fimh_coverage": fimh_coverage,
                    "status": "Completed",
                    "output_directory": str(output_dir),
                    "warnings": [],
                    "detailed_results": detailed_results
                }
            else:
                return self._create_error_result(sample_name, fasta_path, "No results_tab.txt file found")
                
        except Exception as e:
            return self._create_error_result(sample_name, fasta_path, f"Error parsing results: {str(e)}")
    
    def _create_error_result(self, sample_name: str, fasta_path: str, error_msg: str) -> Dict[str, Any]:
        """Create error result structure"""
        return {
            "sample_id": sample_name,
            "file_path": fasta_path,
            "fumc_type": "Unknown",
            "fimh_type": "Unknown",
            "fumc_identity": "0.00",
            "fimh_identity": "0.00",
            "fumc_coverage": "0/0",
            "fimh_coverage": "0/0",
            "status": f"Error: {error_msg}",
            "output_directory": "",
            "warnings": [error_msg],
            "detailed_results": ""
        }
    
    def process_batch(self, input_path: str, main_output_dir: Path) -> List[Dict[str, Any]]:
        """Process all FASTA files in batch"""
        fasta_files = self.find_fasta_files(input_path)
        
        print(f"🔄 Processing {len(fasta_files)} samples using {self.threads} threads...")
        print(f"🔍 Using BLAST at: {self.blast_path}")
        
        # Prepare arguments for parallel processing
        args = [(fasta_file, main_output_dir) for fasta_file in fasta_files]
        
        # Use multiprocessing for parallel execution
        with mp.Pool(processes=self.threads) as pool:
            results = pool.starmap(self.run_chtyper_analysis, args)
        
        self.results = results
        return results
    
    def generate_html_report(self, output_dir: Path) -> str:
        """Generate comprehensive HTML report with rotating science quotes"""
        # JavaScript for rotating quotes
        quotes_js = """
        <script>
            let quotes = %s;
            let currentQuote = 0;
            
            function rotateQuote() {
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) %% quotes.length;
            }
            
            // Rotate every 10 seconds
            setInterval(rotateQuote, 10000);
            
            // Initial display
            document.addEventListener('DOMContentLoaded', function() {
                rotateQuote();
            });
        </script>
        """ % json.dumps(self.science_quotes)
        
        html_content = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>EcoliTyper CHTyper Analysis Report</title>
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
                table {{ 
                    width: 100%; 
                    border-collapse: collapse; 
                    margin: 20px 0; 
                    background: white;
                    border-radius: 8px;
                    overflow: hidden;
                    box-shadow: 0 2px 10px rgba(0,0,0,0.1);
                }}
                th, td {{ 
                    padding: 15px; 
                    text-align: left; 
                    border-bottom: 1px solid #e0e0e0; 
                }}
                th {{ 
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
                .type-badge {{
                    display: inline-block;
                    background: #667eea;
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
                    <h1 style="color: #333; margin: 0; font-size: 2.5em;">🔬 EcoliTyper CHTyper Analysis Report</h1>
                    <p style="color: #666; font-size: 1.2em;">Comprehensive FumC and FimH Typing Results</p>
                </div>
                
                <div class="quote-container">
                    <div id="science-quote" style="font-size: 1.1em;"></div>
                </div>
                
                <div class="card">
                    <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">📊 Analysis Summary</h2>
                    <div class="summary-stats">
                        <div class="stat-card">
                            <h3>Total Samples</h3>
                            <p style="font-size: 2em; margin: 0;">{len(self.results)}</p>
                        </div>
                        <div class="stat-card">
                            <h3>Successful</h3>
                            <p style="font-size: 2em; margin: 0;" class="success">{len([r for r in self.results if r['status'] == 'Completed'])}</p>
                        </div>
                        <div class="stat-card">
                            <h3>Failed</h3>
                            <p style="font-size: 2em; margin: 0;" class="error">{len([r for r in self.results if r['status'] != 'Completed'])}</p>
                        </div>
                    </div>
                    <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
                    <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
                    <p><strong>Threads Used:</strong> {self.threads}</p>
                    <p><strong>BLAST Path:</strong> {self.blast_path}</p>
                </div>
                
                <div class="card">
                    <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🧬 Type Distribution</h2>
                    <div style="margin: 20px 0;">
        """
        
        # Count types for distribution
        fumc_types = {}
        fimh_types = {}
        for result in self.results:
            if result['status'] == 'Completed':
                if result['fumc_type'] != 'Unknown':
                    fumc_types[result['fumc_type']] = fumc_types.get(result['fumc_type'], 0) + 1
                if result['fimh_type'] != 'Unknown':
                    fimh_types[result['fimh_type']] = fimh_types.get(result['fimh_type'], 0) + 1
        
        # Add FumC type badges
        if fumc_types:
            html_content += "<h4>FumC Types:</h4>"
            for fumc_type, count in sorted(fumc_types.items()):
                html_content += f'<span class="type-badge">{fumc_type} ({count})</span>'
        
        # Add FimH type badges
        if fimh_types:
            html_content += "<h4 style='margin-top: 15px;'>FimH Types:</h4>"
            for fimh_type, count in sorted(fimh_types.items()):
                html_content += f'<span class="type-badge">{fimh_type} ({count})</span>'
        
        html_content += """
                    </div>
                </div>
                
                <div class="card">
                    <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🧪 CHTyper Results</h2>
                    <table>
                        <thead>
                            <tr>
                                <th>Sample ID</th>
                                <th>FumC Type</th>
                                <th>FumC Identity</th>
                                <th>FumC Coverage</th>
                                <th>FimH Type</th>
                                <th>FimH Identity</th>
                                <th>FimH Coverage</th>
                                <th>Status</th>
                            </tr>
                        </thead>
                        <tbody>
        """
        
        for result in self.results:
            status_class = "success" if result["status"] == "Completed" else "error"
            
            html_content += f"""
                            <tr>
                                <td><strong>{result['sample_id']}</strong></td>
                                <td><strong style="color: #667eea;">{result['fumc_type']}</strong></td>
                                <td>{result['fumc_identity']}%</td>
                                <td>{result['fumc_coverage']}</td>
                                <td><strong style="color: #667eea;">{result['fimh_type']}</strong></td>
                                <td>{result['fimh_identity']}%</td>
                                <td>{result['fimh_coverage']}</td>
                                <td class="{status_class}">{result['status']}</td>
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
                        Analysis performed using EcoliTyper CGE CHTyper v1.0.0
                    </p>
                </div>
            </div>
        </body>
        </html>
        """
        
        html_file = output_dir / "chtyper_results.html"
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        return str(html_file)
    
    def generate_tsv_report(self, output_dir: Path) -> str:
        """Generate TSV report"""
        data = []
        for result in self.results:
            data.append({
                'Sample_ID': result['sample_id'],
                'FumC_Type': result['fumc_type'],
                'FumC_Identity': result['fumc_identity'],
                'FumC_Coverage': result['fumc_coverage'],
                'FimH_Type': result['fimh_type'],
                'FimH_Identity': result['fimh_identity'],
                'FimH_Coverage': result['fimh_coverage'],
                'Status': result['status'],
                'File_Path': result['file_path']
            })
        
        df = pd.DataFrame(data)
        tsv_file = output_dir / "chtyper_results.tsv"
        df.to_csv(tsv_file, sep='\t', index=False)
        return str(tsv_file)

def main():
    parser = argparse.ArgumentParser(
        description="EcoliTyper Enhanced CHTyper - Auto BLAST detection with batch processing"
    )
    
    parser.add_argument('-i', '--input', required=True, 
                       help='Input FASTA file, directory, or glob pattern')
    parser.add_argument('-db', '--database', default='chtyper_db',
                       help='Path to CHTyper database')
    parser.add_argument('-o', '--output', required=True,
                       help='Main output directory')
    parser.add_argument('-t', '--threads', type=int, default=4,
                       help='Number of threads to use')
    
    args = parser.parse_args()
    
    # Create main output directory
    main_output_dir = Path(args.output) / "chtyper_results"
    main_output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Initialize enhanced CHTyper
        finder = EnhancedCHTyper(args.database, args.threads)
        
        print(finder.ascii_art)
        print("🔬 EcoliTyper Enhanced CHTyper")
        print("=" * 50)
        print(f"Input: {args.input}")
        print(f"Output: {main_output_dir}")
        print(f"Threads: {args.threads}")
        print(f"Auto-detected BLAST: {finder.blast_path}")
        print("=" * 50)
        
        # Process all samples
        results = finder.process_batch(args.input, main_output_dir)
        
        # Generate reports
        print("\n📊 Generating reports...")
        html_file = finder.generate_html_report(main_output_dir)
        tsv_file = finder.generate_tsv_report(main_output_dir)
        
        # Summary
        print("\n✅ Analysis Complete!")
        print(f"📊 Samples processed: {len(results)}")
        print(f"📁 Results directory: {main_output_dir}")
        print(f"📄 HTML Report: {html_file}")
        print(f"📊 TSV Report: {tsv_file}")
        
        successful = len([r for r in results if r['status'] == 'Completed'])
        print(f"🎯 Success rate: {successful}/{len(results)} ({successful/len(results)*100:.1f}%)")
        
        # Print unique types found
        fumc_types = list(set([r['fumc_type'] for r in results if r['status'] == 'Completed' and r['fumc_type'] != 'Unknown']))
        fimh_types = list(set([r['fimh_type'] for r in results if r['status'] == 'Completed' and r['fimh_type'] != 'Unknown']))
        print(f"🧬 FumC types found: {', '.join(fumc_types) if fumc_types else 'None'}")
        print(f"🧬 FimH types found: {', '.join(fimh_types) if fimh_types else 'None'}")
        
        import random
        print(f"\n💡 {random.choice(finder.science_quotes)}")
        
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
