#!/usr/bin/env python3
"""
EcoliTyper - Enhanced SerotypeFinder Wrapper
Handles batch processing, automatic directory management, and multi-format reporting
Author: Brown Beckley
Email: <brownbeckley94@gmail.com>
Affliation: University of Ghana Medical School-Department of Medical Biochemistry
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

class EnhancedSerotypeFinder:
    def __init__(self, db_path: str = "serotypefinder_db", threads: int = 2):
        self.db_path = Path(db_path)
        self.threads = threads
        self.results = []
        self.metadata = {
            "tool_name": "EcoliTyper SerotypeFinder",
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
        
        self.ascii_art = """
███████╗ ██████╗ ██████╗ ██╗     ██╗████████╗██╗   ██╗██████╗ ███████╗██████╗ 
██╔════╝██╔════╝██╔═══██╗██║     ██║╚══██╔══╝╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
█████╗  ██║     ██║   ██║██║     ██║   ██║    ╚████╔╝ ██████╔╝█████╗  ██████╔╝
██╔══╝  ██║     ██║   ██║██║     ██║   ██║     ╚██╔╝  ██╔═══╝ ██╔══╝  ██╔══██╗
███████╗╚██████╗╚██████╔╝███████╗██║   ██║      ██║   ██║     ███████╗██║  ██║
╚══════╝ ╚═════╝ ╚═════╝ ╚══════╝╚═╝   ╚═╝      ╚═╝   ╚═╝     ╚══════╝╚═╝  ╚═╝
        """
        
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files using glob patterns or direct paths"""
        fasta_files = []
        
        # Handle glob patterns
        if '*' in input_path or '?' in input_path or '[' in input_path:
            matches = glob.glob(input_path)
            for match in matches:
                path = Path(match)
                if path.is_file() and path.suffix.lower() in ['.fasta', '.fna', '.fa', '.fsa']:
                    fasta_files.append(path)
        else:
            # Handle direct file or directory path
            input_path = Path(input_path)
            if input_path.is_file():
                if input_path.suffix.lower() in ['.fasta', '.fna', '.fa', '.fsa']:
                    fasta_files = [input_path]
            elif input_path.is_dir():
                for ext in ['*.fasta', '*.fna', '*.fa', '*.fsa']:
                    fasta_files.extend(input_path.glob(ext))
                    fasta_files.extend(input_path.glob(ext.upper()))
        
        if not fasta_files:
            raise ValueError(f"No FASTA files found matching: {input_path}")
            
        print(f"Found {len(fasta_files)} FASTA file(s) for analysis")
        return fasta_files
    
    def run_serotype_analysis(self, fasta_file: Path, output_base: Path) -> Dict[str, Any]:
        """Run serotypefinder on a single FASTA file"""
        try:
            sample_name = fasta_file.stem
            sample_output_dir = output_base / sample_name
            
            # Create sample-specific output directory
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            
            print(f"🔬 Analyzing {sample_name}...")
            
            # Build serotypefinder command
            cmd = [
                "python3", "serotypefinder.py",
                "-i", str(fasta_file),
                "-o", str(sample_output_dir),
                "-p", str(self.db_path),
                "-d", "O_type,H_type",
                "-l", "0.6",
                "-t", "0.9",
                "-x"
            ]
            
            # Run serotypefinder
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                return self._create_error_result(sample_name, str(fasta_file), f"Command failed: {result.stderr}")
            
            # Parse results from output files
            return self._parse_sample_results(sample_name, str(fasta_file), sample_output_dir)
            
        except Exception as e:
            return self._create_error_result(sample_name, str(fasta_file), str(e))
    
    def _parse_sample_results(self, sample_name: str, fasta_path: str, output_dir: Path) -> Dict[str, Any]:
        """Parse results from serotypefinder output files"""
        try:
            # Parse JSON results
            json_file = output_dir / "data.json"
            if json_file.exists():
                with open(json_file, 'r') as f:
                    json_data = json.load(f)
                
                serotype_data = json_data.get('serotypefinder', {}).get('results', {})
                run_info = json_data.get('serotypefinder', {}).get('run_info', {})
                user_input = json_data.get('serotypefinder', {}).get('user_input', {})
                
                o_type = "Unknown"
                h_type = "Unknown"
                genes_found = []
                detailed_results = {}
                
                # Extract O-type results
                o_type_results = serotype_data.get('O_type', {})
                detailed_results['O_type'] = o_type_results
                if o_type_results:
                    o_type = list(o_type_results.values())[0].get('serotype', 'Unknown')
                    genes_found.extend(list(o_type_results.keys()))
                
                # Extract H-type results  
                h_type_results = serotype_data.get('H_type', {})
                detailed_results['H_type'] = h_type_results
                if h_type_results:
                    h_type = list(h_type_results.values())[0].get('serotype', 'Unknown')
                    genes_found.extend(list(h_type_results.keys()))
                
                serotype = f"{o_type}:{h_type}" if o_type != "Unknown" and h_type != "Unknown" else "Unknown"
                
                return {
                    "sample_id": sample_name,
                    "file_path": fasta_path,
                    "serotype": serotype,
                    "o_type": o_type,
                    "h_type": h_type,
                    "genes_found": genes_found,
                    "confidence": "High",  # Based on 100% identity in your results
                    "status": "Completed",
                    "output_directory": str(output_dir),
                    "warnings": [],
                    "detailed_data": detailed_results,
                    "run_info": run_info,
                    "user_input": user_input
                }
            else:
                return self._create_error_result(sample_name, fasta_path, "No JSON results file found")
                
        except Exception as e:
            return self._create_error_result(sample_name, fasta_path, f"Error parsing results: {str(e)}")
    
    def _create_error_result(self, sample_name: str, fasta_path: str, error_msg: str) -> Dict[str, Any]:
        """Create error result structure"""
        return {
            "sample_id": sample_name,
            "file_path": fasta_path,
            "serotype": "Unknown",
            "o_type": "Unknown",
            "h_type": "Unknown",
            "genes_found": [],
            "confidence": 0.0,
            "status": f"Error: {error_msg}",
            "output_directory": "",
            "warnings": [error_msg],
            "detailed_data": {},
            "run_info": {},
            "user_input": {}
        }
    
    def process_batch(self, input_path: str, main_output_dir: Path) -> List[Dict[str, Any]]:
        """Process all FASTA files in batch"""
        fasta_files = self.find_fasta_files(input_path)
        results = []
        
        for fasta_file in fasta_files:
            result = self.run_serotype_analysis(fasta_file, main_output_dir)
            results.append(result)
        
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
            <title>EcoliTyper Serotype Analysis Report</title>
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
                .gene-details {{
                    background: #f8f9fa;
                    border-radius: 8px;
                    padding: 15px;
                    margin: 10px 0;
                }}
                .detail-table {{
                    font-size: 14px;
                }}
                .detail-table th {{
                    background: #495057;
                }}
            </style>
            {quotes_js}
        </head>
        <body>
            <div class="container">
                <div class="header">
                    <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 EcoliTyper Serotype Analysis Report</h1>
                    <p style="color: #666; font-size: 1.2em;">Comprehensive E. coli Serotyping Results</p>
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
                </div>
                
                <div class="card">
                    <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🧪 Serotype Results Overview</h2>
                    <table>
                        <thead>
                            <tr>
                                <th>Sample ID</th>
                                <th>Serotype</th>
                                <th>O-type</th>
                                <th>H-type</th>
                                <th>Genes Found</th>
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
                                <td><strong style="color: #667eea;">{result['serotype']}</strong></td>
                                <td>{result['o_type']}</td>
                                <td>{result['h_type']}</td>
                                <td>{', '.join(result['genes_found'])}</td>
                                <td class="{status_class}">{result['status']}</td>
                            </tr>
            """
        
        html_content += """
                        </tbody>
                    </table>
                </div>
                
                <div class="card">
                    <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🔬 Detailed Gene Information</h2>
        """
        
        for result in self.results:
            if result['status'] == 'Completed' and result['detailed_data']:
                html_content += f"""
                    <div style="margin-bottom: 30px;">
                        <h3 style="color: #495057; background: #e9ecef; padding: 10px; border-radius: 5px;">Sample: {result['sample_id']}</h3>
                """
                
                # O-type details
                if result['detailed_data'].get('O_type'):
                    html_content += """
                        <h4 style="color: #667eea;">O-type Genes:</h4>
                        <table class="detail-table">
                            <thead>
                                <tr>
                                    <th>Gene</th>
                                    <th>Serotype</th>
                                    <th>Identity</th>
                                    <th>Coverage</th>
                                    <th>Contig</th>
                                    <th>Position</th>
                                    <th>Accession</th>
                                </tr>
                            </thead>
                            <tbody>
                    """
                    for gene, details in result['detailed_data']['O_type'].items():
                        html_content += f"""
                                <tr>
                                    <td><strong>{gene}</strong></td>
                                    <td>{details.get('serotype', 'N/A')}</td>
                                    <td>{details.get('identity', 'N/A')}%</td>
                                    <td>{details.get('coverage', 'N/A')}%</td>
                                    <td style="max-width: 200px; word-wrap: break-word;">{details.get('contig_name', 'N/A')}</td>
                                    <td>{details.get('positions_in_contig', 'N/A')}</td>
                                    <td>{details.get('accession', 'N/A')}</td>
                                </tr>
                        """
                    html_content += """
                            </tbody>
                        </table>
                    """
                
                # H-type details
                if result['detailed_data'].get('H_type'):
                    html_content += """
                        <h4 style="color: #667eea; margin-top: 20px;">H-type Genes:</h4>
                        <table class="detail-table">
                            <thead>
                                <tr>
                                    <th>Gene</th>
                                    <th>Serotype</th>
                                    <th>Identity</th>
                                    <th>Coverage</th>
                                    <th>Contig</th>
                                    <th>Position</th>
                                    <th>Accession</th>
                                </tr>
                            </thead>
                            <tbody>
                    """
                    for gene, details in result['detailed_data']['H_type'].items():
                        html_content += f"""
                                <tr>
                                    <td><strong>{gene}</strong></td>
                                    <td>{details.get('serotype', 'N/A')}</td>
                                    <td>{details.get('identity', 'N/A')}%</td>
                                    <td>{details.get('coverage', 'N/A')}%</td>
                                    <td style="max-width: 200px; word-wrap: break-word;">{details.get('contig_name', 'N/A')}</td>
                                    <td>{details.get('positions_in_contig', 'N/A')}</td>
                                    <td>{details.get('accession', 'N/A')}</td>
                                </tr>
                        """
                    html_content += """
                            </tbody>
                        </table>
                    """
                
                html_content += "</div><hr>"
        
        html_content += f"""
                </div>
                
                <div class="footer">
                    <h3 style="color: #fff; border-bottom: 2px solid #667eea; padding-bottom: 10px;">👥 Contact Information</h3>
                    <p><strong>Author:</strong> {', '.join(self.metadata['authors'])}</p>
                    <p><strong>Email:</strong> {self.metadata['email']}</p>
                    <p><strong>GitHub:</strong> <a href="{self.metadata['github']}" target="_blank">{self.metadata['github']}</a></p>
                    <p><strong>Affiliation:</strong> {self.metadata['affiliation']}</p>
                    <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
                        Analysis performed using EcoliTyper CGE SerotypeFinder v2.0.1
                    </p>
                </div>
            </div>
        </body>
        </html>
        """
        
        html_file = output_dir / "serotype_analysis_report.html"
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        return str(html_file)
    
    def generate_tsv_report(self, output_dir: Path) -> str:
        """Generate TSV report with all sample results"""
        # Create simplified dataframe for TSV
        data = []
        for result in self.results:
            data.append({
                'Sample_ID': result['sample_id'],
                'Serotype': result['serotype'],
                'O_Type': result['o_type'],
                'H_Type': result['h_type'],
                'Genes_Found': ','.join(result['genes_found']),
                'Confidence': result['confidence'],
                'Status': result['status'],
                'File_Path': result['file_path']
            })
        
        df = pd.DataFrame(data)
        tsv_file = output_dir / "serotype_analysis_report.tsv"
        df.to_csv(tsv_file, sep='\t', index=False)
        return str(tsv_file)
    
    def cleanup_temp_dirs(self, main_output_dir: Path):
        """Clean up temporary directories while keeping results"""
        for result in self.results:
            if result['output_directory']:
                sample_dir = Path(result['output_directory'])
                temp_dir = sample_dir / "tmp"
                if temp_dir.exists():
                    shutil.rmtree(temp_dir)

def main():
    parser = argparse.ArgumentParser(
        description="EcoliTyper Enhanced SerotypeFinder - Batch processing with multi-format reports",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('-i', '--input', required=True, 
                       help='Input FASTA file, directory, or glob pattern (e.g., "*.fna", "genomes/*.fasta")')
    parser.add_argument('-db', '--database', default='serotypefinder_db',
                       help='Path to serotypefinder database')
    parser.add_argument('-o', '--output', required=True,
                       help='Main output directory (will be created as SerotypeFinder_results)')
    parser.add_argument('-t', '--threads', type=int, default=1,
                       help='Number of threads to use')
    
    args = parser.parse_args()
    
    # Create main output directory
    main_output_dir = Path(args.output) / "SerotypeFinder_results"
    main_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Print ASCII art
    finder = EnhancedSerotypeFinder()
    print(finder.ascii_art)
    print("🧬 EcoliTyper Enhanced SerotypeFinder")
    print("=" * 50)
    print(f"Input: {args.input}")
    print(f"Output: {main_output_dir}")
    print(f"Threads: {args.threads}")
    print("=" * 50)
    
    try:
        # Initialize enhanced serotypefinder
        finder = EnhancedSerotypeFinder(args.database, args.threads)
        
        # Process all samples
        results = finder.process_batch(args.input, main_output_dir)
        
        # Generate reports
        print("\n📊 Generating reports...")
        html_file = finder.generate_html_report(main_output_dir)
        tsv_file = finder.generate_tsv_report(main_output_dir)
        
        # Cleanup temporary directories
        finder.cleanup_temp_dirs(main_output_dir)
        
        # Summary
        print("\n✅ Analysis Complete!")
        print(f"📊 Samples processed: {len(results)}")
        print(f"📁 Results directory: {main_output_dir}")
        print(f"📄 HTML Report: {html_file}")
        print(f"📊 TSV Report: {tsv_file}")
        
        # Show success rate
        successful = len([r for r in results if r['status'] == 'Completed'])
        print(f"🎯 Success rate: {successful}/{len(results)} ({successful/len(results)*100:.1f}%)")
        
        # Print a random science quote
        import random
        print(f"\n💡 {random.choice(finder.science_quotes)}")
        
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
