#!/usr/bin/env python3
"""
EcoliTyper - Enhanced ezClermont Wrapper
E. coli phylogrouping using ezClermont with beautiful reporting
Author: Brown Beckley
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
import re
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any
import pandas as pd
import multiprocessing as mp

class EnhancedEzClermont:
    def __init__(self, threads: int = 4):
        self.threads = threads
        self.results = []
        self.metadata = {
            "tool_name": "EcoliTyper Phylogrouping",
            "version": "1.0.0", 
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        
        # Check if ezClermont is available
        self.ezclermont_path = self._detect_ezclermont()
        if not self.ezclermont_path:
            raise RuntimeError("ezClermont not found. Please install ezClermont or ensure it's in your PATH")
        
        self.science_quotes = [
            "“The important thing is not to stop questioning. Curiosity has its own reason for existence.” - Albert Einstein",
            "“Nothing in life is to be feared, it is only to be understood.” - Marie Curie", 
            "“The microscope opens a new world to the investigator.” - Robert Koch",
            "“In science, the credit goes to the man who convinces the world, not to the man to whom the idea first occurs.” - Francis Darwin",
            "“The good thing about science is that it's true whether or not you believe in it.” - Neil deGrasse Tyson",
            "“Science knows no country, because knowledge belongs to humanity.” - Louis Pasteur"
        ]
        
        self.ascii_art = """
██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ 
██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗
██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║  
██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║ 
██║     ██║  ██║   ██║   ███████╗╚██████╔╝
╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝  
        """
    
    def _detect_ezclermont(self) -> str:
        """Auto-detect ezClermont installation"""
        try:
            result = subprocess.run(["which", "ezclermont"], capture_output=True, text=True)
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
    
    def run_ezclermont_analysis(self, fasta_file: Path, output_base: Path) -> Dict[str, Any]:
        """Run ezClermont on a single FASTA file - IMPROVED STATUS HANDLING"""
        try:
            sample_name = fasta_file.stem
            sample_output_dir = output_base / sample_name
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            
            print(f"🔬 Analyzing {sample_name}...")
            
            # Copy input file to sample directory
            sample_fasta = sample_output_dir / fasta_file.name
            shutil.copy2(fasta_file, sample_fasta)
            
            # Build command using the COPIED file in the sample directory
            cmd = [
                self.ezclermont_path, 
                str(sample_fasta.absolute()),
                "-e", 
                sample_name
            ]
            
            # Run with timeout in the sample output directory
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                timeout=300, 
                cwd=str(sample_output_dir.absolute())
            )
            
            # Save raw output for debugging
            output_file = sample_output_dir / "ezclermont_raw_output.txt"
            with open(output_file, 'w') as f:
                f.write("STDOUT:\n" + result.stdout)
                if result.stderr:
                    f.write("\nSTDERR:\n" + result.stderr)
                f.write(f"\nReturn code: {result.returncode}")
            
            # ALWAYS try to parse results regardless of return code
            parsed_result = self._parse_sample_results(sample_name, str(fasta_file), sample_output_dir, result.stdout)
            
            # IMPROVED STATUS HANDLING:
            if result.returncode != 0:
                # Check if we successfully parsed a valid Clermont type despite the exit code
                if parsed_result["clermont_type"] != "Unknown":
                    # Successfully got result despite non-zero exit code - treat as completed
                    parsed_result["status"] = "Completed"
                    parsed_result["warnings"].append(f"ezClermont returned exit code {result.returncode} but analysis completed successfully")
                else:
                    # Couldn't parse result and non-zero exit code - treat as warning
                    parsed_result["status"] = f"Error: ezClermont failed with exit code {result.returncode}"
                    parsed_result["warnings"].append(f"ezClermont exited with code {result.returncode}")
                
                # Save stderr for debugging
                if result.stderr:
                    stderr_file = sample_output_dir / "stderr.log"
                    with open(stderr_file, 'w') as f:
                        f.write(result.stderr)
            else:
                parsed_result["status"] = "Completed"
            
            return parsed_result
                
        except subprocess.TimeoutExpired:
            return self._create_error_result(sample_name, str(fasta_file), "Analysis timed out after 5 minutes")
        except Exception as e:
            return self._create_error_result(sample_name, str(fasta_file), str(e))
    
    def _parse_sample_results(self, sample_name: str, fasta_path: str, output_dir: Path, output_text: str) -> Dict[str, Any]:
        """Parse results from ezClermont output - UPDATED BASED ON ACTUAL OUTPUT"""
        try:
            # Initialize with defaults
            result_data = {
                "sample_id": sample_name,
                "file_path": fasta_path,
                "clermont_type": "Unknown",
                "tspe4": "Unknown",
                "arpa": "Unknown", 
                "chu": "Unknown",
                "yjaa": "Unknown",
                "trpba_control": "Unknown",
                "vira": "Unknown",
                "status": "Completed",
                "output_directory": str(output_dir),
                "warnings": [],
                "raw_output": output_text
            }
            
            lines = output_text.strip().split('\n')
            
            # Parse line by line
            for line in lines:
                line = line.strip()
                
                # Parse control genes
                if line.startswith('trpBA_control:'):
                    result_data["trpba_control"] = self._extract_gene_result(line)
                elif line.startswith('virA:'):
                    result_data["vira"] = self._extract_gene_result(line)
                
                # Parse Quadriplex PCR results
                elif line.startswith('TspE4:'):
                    result_data["tspe4"] = self._extract_gene_result(line)
                elif line.startswith('arpA:'):
                    result_data["arpa"] = self._extract_gene_result(line)
                elif line.startswith('chu:'):
                    result_data["chu"] = self._extract_gene_result(line)
                elif line.startswith('yjaA:'):
                    result_data["yjaa"] = self._extract_gene_result(line)
                
                # Parse Clermont type from formatted output
                elif line.startswith('Clermont type:'):
                    result_data["clermont_type"] = line.split(':')[1].strip()
                
                # Parse final tab-separated result line: "test_sample	B1"
                elif '\t' in line and len(line.split('\t')) == 2:
                    parts = line.split('\t')
                    if parts[1] in ['A', 'B1', 'B2', 'C', 'D', 'E', 'F', 'G']:
                        result_data["clermont_type"] = parts[1]
            
            # Special handling for the tuple format at the end: "('B1', 'TspE4: +\\narpA: +\\nchu: -\\nyjaA: -\\n')"
            for line in lines:
                if line.startswith("('") and line.endswith("')"):
                    # Extract from format: ('B1', 'TspE4: +\\narpA: +\\nchu: -\\nyjaA: -\\n')
                    tuple_match = re.search(r"\('([^']+)',\s*'([^']+)'\)", line)
                    if tuple_match:
                        clermont_type = tuple_match.group(1)
                        markers_text = tuple_match.group(2).replace('\\n', '\n')
                        
                        # Use the tuple result as final authority
                        result_data["clermont_type"] = clermont_type
                        
                        # Parse markers from tuple if we haven't found them already
                        if result_data["tspe4"] == "Unknown":
                            tspe4_match = re.search(r'TspE4:\s*([+-])', markers_text)
                            if tspe4_match:
                                result_data["tspe4"] = tspe4_match.group(1)
                        if result_data["arpa"] == "Unknown":
                            arpa_match = re.search(r'arpA:\s*([+-])', markers_text)
                            if arpa_match:
                                result_data["arpa"] = arpa_match.group(1)
                        if result_data["chu"] == "Unknown":
                            chu_match = re.search(r'chu:\s*([+-])', markers_text)
                            if chu_match:
                                result_data["chu"] = chu_match.group(1)
                        if result_data["yjaa"] == "Unknown":
                            yjaa_match = re.search(r'yjaA:\s*([+-])', markers_text)
                            if yjaa_match:
                                result_data["yjaa"] = yjaa_match.group(1)
            
            # Final validation and warnings
            self._validate_results(result_data)
            
            return result_data
                
        except Exception as e:
            return self._create_error_result(sample_name, fasta_path, f"Error parsing results: {str(e)}")
    
    def _extract_gene_result(self, line: str) -> str:
        """Extract + or - result from a gene line"""
        if ':' in line:
            parts = line.split(':')
            if len(parts) > 1:
                value = parts[1].strip()
                if value == '+':
                    return '+'
                elif value == '-':
                    return '-'
        return "Unknown"
    
    def _validate_results(self, result_data: Dict[str, Any]):
        """Validate parsed results and add warnings if needed"""
        warnings = []
        
        if result_data["clermont_type"] == "Unknown":
            warnings.append("Could not determine Clermont phylogroup")
        
        unknown_markers = []
        for marker in ['tspe4', 'arpa', 'chu', 'yjaa', 'trpba_control', 'vira']:
            if result_data[marker] == "Unknown":
                unknown_markers.append(marker)
        
        if unknown_markers:
            warnings.append(f"Could not parse markers: {', '.join(unknown_markers)}")
        
        result_data["warnings"] = warnings
    
    def _create_error_result(self, sample_name: str, fasta_path: str, error_msg: str) -> Dict[str, Any]:
        """Create error result structure"""
        return {
            "sample_id": sample_name,
            "file_path": fasta_path,
            "clermont_type": "Unknown",
            "tspe4": "Unknown",
            "arpa": "Unknown",
            "chu": "Unknown", 
            "yjaa": "Unknown",
            "trpba_control": "Unknown",
            "vira": "Unknown",
            "status": f"Error: {error_msg}",
            "output_directory": "",
            "warnings": [error_msg],
            "raw_output": ""
        }
    
    def process_batch(self, input_path: str, main_output_dir: Path) -> List[Dict[str, Any]]:
        """Process all FASTA files in batch"""
        fasta_files = self.find_fasta_files(input_path)
        
        print(f"🔄 Processing {len(fasta_files)} samples using {self.threads} threads...")
        print(f"🔍 Using ezClermont at: {self.ezclermont_path}")
        
        # Prepare arguments for parallel processing
        args = [(fasta_file, main_output_dir) for fasta_file in fasta_files]
        
        # Use multiprocessing for parallel execution
        with mp.Pool(processes=self.threads) as pool:
            results = pool.starmap(self.run_ezclermont_analysis, args)
        
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
            <title>EcoliTyper Phylogrouping Analysis Report</title>
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
                .info-box {{
                    background: #e7f3ff;
                    border-left: 4px solid #2196F3;
                    padding: 15px;
                    margin: 15px 0;
                    border-radius: 4px;
                }}
            </style>
            {quotes_js}
        </head>
        <body>
            <div class="container">
                <div class="header">
                    <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 EcoliTyper Phylogrouping Analysis Report</h1>
                    <p style="color: #666; font-size: 1.2em;">Comprehensive E. coli Clermont Phylotyping Results</p>
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
                            <p style="font-size: 2em; margin: 0;" class="error">{len([r for r in self.results if r['status'].startswith('Error')])}</p>
                        </div>
                    </div>
                    <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
                    <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
                    <p><strong>Threads Used:</strong> {self.threads}</p>
                    <p><strong>ezClermont Path:</strong> {self.ezclermont_path}</p>
                </div>
                
                <div class="card">
                    <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🧬 Clermont Phylogrouping Method</h2>
                    <div class="info-box">
                        <p><strong>Note:</strong> The Clermont phylogrouping is based on the original Clermont algorithm that uses specific gene markers including <em>chuA</em>, <em>yjaA</em>, <em>TspE4.C2</em>, and <em>arpA</em> to classify E. coli into phylogroups A, B1, B2, C, D, E, F, and G.</p>
                        <p>The algorithm determines phylogroups based on the presence/absence patterns of these key genetic markers following the established Clermont typing scheme.</p>
                    </div>
                </div>
                
                <div class="card">
                    <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🧪 Clermont Type Distribution</h2>
                    <div style="margin: 20px 0;">
        """
        
        # Count Clermont types for distribution
        clermont_types = {}
        for result in self.results:
            if result['status'] == 'Completed' and result['clermont_type'] != 'Unknown':
                clermont_types[result['clermont_type']] = clermont_types.get(result['clermont_type'], 0) + 1
        
        # Add Clermont type badges
        for clermont_type, count in sorted(clermont_types.items()):
            html_content += f'<span class="type-badge">{clermont_type} ({count})</span>'
        
        html_content += """
                    </div>
                </div>
                
                <div class="card">
                    <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🔬 Phylogrouping Results</h2>
                    <table>
                        <thead>
                            <tr>
                                <th>Sample ID</th>
                                <th>Clermont Type</th>
                                <th>Status</th>
                            </tr>
                        </thead>
                        <tbody>
        """
        
        for result in self.results:
            if result["status"] == "Completed":
                status_class = "success"
            elif result["status"].startswith("Error"):
                status_class = "error"
            else:
                status_class = "warning"
            
            html_content += f"""
                            <tr>
                                <td><strong>{result['sample_id']}</strong></td>
                                <td><strong style="color: #667eea;">{result['clermont_type']}</strong></td>
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
                        Analysis performed using EcoliTyper ezClermont V0.7.0
                    </p>
                </div>
            </div>
        </body>
        </html>
        """
        
        html_file = output_dir / "phylogrouping_results.html"
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        return str(html_file)
    
    def generate_tsv_report(self, output_dir: Path) -> str:
        """Generate TSV report - SIMPLIFIED WITHOUT GENE COLUMNS"""
        data = []
        for result in self.results:
            data.append({
                'Sample_ID': result['sample_id'],
                'Clermont_Type': result['clermont_type'],
                'Status': result['status'],
                'File_Path': result['file_path']
            })
        
        df = pd.DataFrame(data)
        tsv_file = output_dir / "phylogrouping_results.tsv"
        df.to_csv(tsv_file, sep='\t', index=False)
        return str(tsv_file)

def main():
    parser = argparse.ArgumentParser(
        description="EcoliTyper Enhanced ezClermont - Phylogrouping with batch processing"
    )
    
    parser.add_argument('-i', '--input', required=True, 
                       help='Input FASTA file, directory, or glob pattern')
    parser.add_argument('-o', '--output', required=True,
                       help='Main output directory')
    parser.add_argument('-t', '--threads', type=int, default=4,
                       help='Number of threads to use')
    
    args = parser.parse_args()
    
    # Create main output directory
    main_output_dir = Path(args.output) / "phylogrouping_results"
    main_output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Initialize enhanced ezClermont
        finder = EnhancedEzClermont(args.threads)
        
        print(finder.ascii_art)
        print("🧬 EcoliTyper Enhanced ezClermont")
        print("=" * 50)
        print(f"Input: {args.input}")
        print(f"Output: {main_output_dir}")
        print(f"Threads: {args.threads}")
        print(f"Auto-detected ezClermont: {finder.ezclermont_path}")
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
        
        # Print unique Clermont types found
        clermont_types = list(set([r['clermont_type'] for r in results if r['status'] == 'Completed' and r['clermont_type'] != 'Unknown']))
        print(f"🧬 Clermont types found: {', '.join(clermont_types) if clermont_types else 'None'}")
        
        import random
        print(f"\n💡 {random.choice(finder.science_quotes)}")
        
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()