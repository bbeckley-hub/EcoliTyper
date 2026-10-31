#!/usr/bin/env python3
"""
EcoliTyper - Enhanced SerotypeFinder Wrapper
Handles batch processing, automatic directory management, and multi-format reporting.

Per-sample output (1 file each):
    <sample>/<sample>_serotype_report.html

Batch-level outputs (at the top of -o):
    serotype_analysis_report.html
    serotype_analysis_report.tsv
    serotype_analysis_report.json

All serotypefinder.py intermediates (data.json, results files, tmp/) live
inside a TemporaryDirectory that self-destructs on exit — no intermediates
are left on disk.

Author: Brown Beckley
Email: <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
Date: 2025-12-16 / 2026-07-22
Send a quick mail for any issues or further explanations.
"""

import os
import sys
import json
import argparse
import subprocess
import glob
import tempfile
import random
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any
import pandas as pd


class EnhancedSerotypeFinder:
    """Enhanced SerotypeFinder wrapper with per-sample HTML reporting."""

    def __init__(self, db_path: str = "serotypefinder_db", threads: int = 2):
        self.db_path = Path(db_path)
        self.threads = threads
        self.results = []
        self.metadata = {
            "tool_name": "EcoliTyper SerotypeFinder",
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
            {"text": "The good thing about science is that it's true whether or not you believe in it.",
             "author": "Neil deGrasse Tyson"},
            {"text": "Science knows no country, because knowledge belongs to humanity.",
             "author": "Louis Pasteur"},
            {"text": "Somewhere, something incredible is waiting to be known.",
             "author": "Carl Sagan"},
            {"text": "Research is what I'm doing when I don't know what I'm doing.",
             "author": "Wernher von Braun"},
            {"text": "EcoliTyper transforms raw E. coli sequences into actionable genomic surveillance.",
             "author": "Brown Beckley"},
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
        """Find all FASTA files using glob patterns or direct paths."""
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
                for ext in ['*.fasta', '*.fna', '*.fa', '*.fsa']:
                    fasta_files.extend(input_path.glob(ext))
                    fasta_files.extend(input_path.glob(ext.upper()))

        if not fasta_files:
            raise ValueError(f"No FASTA files found matching: {input_path}")

        print(f"Found {len(fasta_files)} FASTA file(s) for analysis")
        return fasta_files

    def run_serotype_analysis(self, fasta_file: Path, output_base: Path) -> Dict[str, Any]:
        """Run serotypefinder on one genome. Emits exactly 1 file per sample:
             - <sample>/<sample>_serotype_report.html
        All serotypefinder.py intermediates live inside a TemporaryDirectory
        that self-destructs on exit.
        """
        sample_name = fasta_file.stem
        try:
            print(f"Analyzing {sample_name}...")

            with tempfile.TemporaryDirectory(prefix=f"serotype_{sample_name}_") as temp_dir:
                temp_path = Path(temp_dir)

                script_dir = os.path.dirname(__file__)
                serotypefinder_script = os.path.join(script_dir, "serotypefinder.py")

                cmd = [
                    sys.executable,
                    serotypefinder_script,
                    "-i", str(fasta_file),
                    "-o", str(temp_path),
                    "-p", str(self.db_path),
                    "-d", "O_type,H_type",
                    "-l", "0.6",
                    "-t", "0.9",
                    "-x",
                ]

                result = subprocess.run(cmd, capture_output=True, text=True)

                if result.returncode != 0:
                    parsed = self._create_error_result(
                        sample_name, str(fasta_file),
                        f"Command failed: {result.stderr}",
                    )
                else:
                    parsed = self._parse_sample_results_from_temp(
                        sample_name, str(fasta_file), temp_path
                    )

            sample_output_dir = output_base / sample_name
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            self._generate_per_sample_html(parsed, sample_output_dir)
            parsed["output_directory"] = str(sample_output_dir)

            return parsed

        except Exception as e:
            error = self._create_error_result(sample_name, str(fasta_file), str(e))
            try:
                sample_output_dir = output_base / sample_name
                sample_output_dir.mkdir(parents=True, exist_ok=True)
                self._generate_per_sample_html(error, sample_output_dir)
                error["output_directory"] = str(sample_output_dir)
            except Exception:
                pass
            return error

    def _parse_sample_results_from_temp(self, sample_name: str, fasta_path: str,
                                        temp_dir: Path) -> Dict[str, Any]:
        """Parse serotypefinder output files from the temporary working directory."""
        try:
            json_file = temp_dir / "data.json"
            if not json_file.exists():
                return self._create_error_result(
                    sample_name, fasta_path, "No JSON results file found"
                )

            with open(json_file, 'r') as f:
                json_data = json.load(f)

            serotype_data = json_data.get('serotypefinder', {}).get('results', {})
            if not isinstance(serotype_data, dict):
                serotype_data = {}

            run_info = json_data.get('serotypefinder', {}).get('run_info', {})
            user_input = json_data.get('serotypefinder', {}).get('user_input', {})

            o_type = "Unknown"
            h_type = "Unknown"
            genes_found = []
            detailed_results = {'O_type': {}, 'H_type': {}}

            # --- Parse O-type ---
            o_type_results = serotype_data.get('O_type', {})
            if isinstance(o_type_results, dict):
                detailed_results['O_type'] = o_type_results
                if o_type_results:
                    first_result = next(iter(o_type_results.values())) if o_type_results else {}
                    if isinstance(first_result, dict):
                        o_type = first_result.get('serotype', 'Unknown')
                        genes_found.extend(list(o_type_results.keys()))
            elif isinstance(o_type_results, list):
                for item in o_type_results:
                    if isinstance(item, dict):
                        detailed_results['O_type'][item.get('gene', 'unknown')] = item
                        if item.get('serotype') and o_type == "Unknown":
                            o_type = item.get('serotype')
                        genes_found.append(item.get('gene', ''))
                if o_type == "Unknown" and o_type_results:
                    o_type = o_type_results[0].get('serotype', 'Unknown')
            else:
                detailed_results['O_type'] = {'raw': str(o_type_results) if o_type_results else 'None'}

            # --- Parse H-type ---
            h_type_results = serotype_data.get('H_type', {})
            if isinstance(h_type_results, dict):
                detailed_results['H_type'] = h_type_results
                if h_type_results:
                    first_result = next(iter(h_type_results.values())) if h_type_results else {}
                    if isinstance(first_result, dict):
                        h_type = first_result.get('serotype', 'Unknown')
                        genes_found.extend(list(h_type_results.keys()))
            elif isinstance(h_type_results, list):
                for item in h_type_results:
                    if isinstance(item, dict):
                        detailed_results['H_type'][item.get('gene', 'unknown')] = item
                        if item.get('serotype') and h_type == "Unknown":
                            h_type = item.get('serotype')
                        genes_found.append(item.get('gene', ''))
                if h_type == "Unknown" and h_type_results:
                    h_type = h_type_results[0].get('serotype', 'Unknown')
            else:
                detailed_results['H_type'] = {'raw': str(h_type_results) if h_type_results else 'None'}

            if o_type == "Unknown" and h_type == "Unknown":
                serotype = "Unknown"
            else:
                serotype = f"{o_type}:{h_type}"

            return {
                "sample_id": sample_name,
                "file_path": fasta_path,
                "serotype": serotype,
                "o_type": o_type,
                "h_type": h_type,
                "genes_found": genes_found,
                "confidence": "High" if genes_found else "Low",
                "status": "Completed",
                "output_directory": "",
                "warnings": [],
                "detailed_data": detailed_results,
                "run_info": run_info,
                "user_input": user_input,
            }

        except Exception as e:
            return self._create_error_result(
                sample_name, fasta_path, f"Error parsing results: {str(e)}"
            )

    def _create_error_result(self, sample_name: str, fasta_path: str,
                             error_msg: str) -> Dict[str, Any]:
        """Create error result structure."""
        return {
            "sample_id": sample_name,
            "file_path": fasta_path,
            "serotype": "Unknown",
            "o_type": "Unknown",
            "h_type": "Unknown",
            "genes_found": [],
            "confidence": "Low",
            "status": f"Error: {error_msg}",
            "output_directory": "",
            "warnings": [error_msg],
            "detailed_data": {},
            "run_info": {},
            "user_input": {},
        }

    def process_batch(self, input_path: str, main_output_dir: Path) -> List[Dict[str, Any]]:
        """Process all FASTA files in batch."""
        fasta_files = self.find_fasta_files(input_path)
        results = []

        for fasta_file in fasta_files:
            result = self.run_serotype_analysis(fasta_file, main_output_dir)
            results.append(result)

        self.results = results
        return results

    # ------------------------------------------------------------------
    # Shared HTML helpers
    # ------------------------------------------------------------------

    def _get_shared_css(self) -> str:
        """Return shared CSS used by per-sample and batch HTML reports."""
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
            display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px; margin-top: 15px;
        }
        .metric-card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white; padding: 20px; border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15); text-align: center;
        }
        .metric-label { font-size: 13px; opacity: 0.9; margin-bottom: 5px; text-transform: uppercase; letter-spacing: 0.5px; }
        .metric-value { font-size: 22px; font-weight: bold; }
        .primary-result {
            color: white; padding: 30px; border-radius: 12px;
            margin: 20px 0; box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            text-align: center;
        }
        .primary-label {
            font-size: 13px; opacity: 0.9; letter-spacing: 3px;
            margin-bottom: 12px; text-transform: uppercase;
        }
        .primary-value { font-size: 56px; font-weight: bold; margin-bottom: 15px; }
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
        .summary-table a { color: #2563eb; text-decoration: none; }
        .summary-table a:hover { text-decoration: underline; }
        .product-cell {
            white-space: normal !important; word-wrap: break-word;
            max-width: 400px; min-width: 180px;
        }
        .stats-grid {
            display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px; margin-bottom: 20px;
        }
        .stat-card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white; padding: 15px; border-radius: 8px;
            text-align: center; box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }
        .stat-value { font-size: 24px; font-weight: bold; margin-bottom: 5px; }
        .stat-label { font-size: 12px; opacity: 0.9; }
        .marker-badge {
            display: inline-block; padding: 4px 14px; border-radius: 12px;
            font-weight: bold; font-size: 14px; color: white; margin: 2px;
        }
        .marker-present { background: #16a34a; }
        .marker-absent { background: #dc2626; }
        .marker-unknown { background: #6b7280; }
        .marker-info { background: #2563eb; }
        .warning-box {
            background: #fef3c7; border-left: 4px solid #f59e0b;
            padding: 15px; margin: 15px 0; border-radius: 6px; color: #92400e;
        }
        .info-box {
            background: #eff6ff; border-left: 4px solid #3b82f6;
            padding: 15px; margin: 15px 0; border-radius: 6px; color: #1e3a8a;
            font-size: 0.95em; line-height: 1.6;
        }
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
        .table-responsive { width: 100%; overflow-x: auto; margin: 20px 0; }
        @media (max-width: 768px) {
            .ascii-art { font-size: 6px; }
            .metrics-grid { grid-template-columns: 1fr; }
            .primary-value { font-size: 32px; }
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
            <p><strong>EcoliTyper</strong> - {title}</p>
            <p class="timestamp">Generated: {current_time}</p>
            <div class="authorship">
                <p><strong>Technical Support &amp; Inquiries:</strong></p>
                <p>Author: Brown Beckley | GitHub: <a href="https://github.com/bbeckley-hub" target="_blank">bbeckley-hub</a></p>
                <p>Email: <a href="mailto:brownbeckley94@gmail.com">brownbeckley94@gmail.com</a></p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
                <p style="margin-top:10px;font-size:11px;opacity:0.8;">
                    Powered by SerotypeFinder (Joensen et al.).<br>
                    Joensen KG, et al. Rapid and easy in silico serotyping of Escherichia coli using whole genome sequencing (WGS) data.
                    J Clin Microbiol. 2015;53(8):2410-2426. doi:10.1128/JCM.00008-15
                </p>
            </div>
        </div>
        """

    def _render_gene_table_html(self, gene_dict: Dict, gene_label: str) -> str:
        """Render a gene detail table (for O-type or H-type)."""
        if not gene_dict:
            return f'<p style="color:#6b7280;">No {gene_label} gene details available.</p>'

        rows = ''
        for gene, details in gene_dict.items():
            if isinstance(details, dict):
                rows += f"""
                    <tr>
                        <td><strong>{gene}</strong></td>
                        <td>{details.get('serotype', 'N/A')}</td>
                        <td>{details.get('identity', 'N/A')}%</td>
                        <td>{details.get('coverage', 'N/A')}%</td>
                        <td class="product-cell">{details.get('contig_name', 'N/A')}</td>
                        <td>{details.get('positions_in_contig', 'N/A')}</td>
                        <td>{details.get('accession', 'N/A')}</td>
                    </tr>
                """
            else:
                rows += f'<tr><td colspan="7"><em>{details}</em></td></tr>'

        return f"""
            <div class="table-responsive">
                <table class="summary-table">
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
                    <tbody>{rows}</tbody>
                </table>
            </div>
        """

    # ------------------------------------------------------------------
    # Per-sample HTML
    # ------------------------------------------------------------------

    def _generate_per_sample_html(self, result: Dict[str, Any], output_dir: Path):
        """Write the single per-sample serotype HTML report."""
        sample_name = result['sample_id']
        random_quote = random.choice(self.science_quotes)
        is_success = result['status'] == 'Completed' and result['serotype'] != 'Unknown'

        if is_success:
            primary_value = result['serotype']
            primary_label = 'Serotype (O : H)'
            primary_bg = 'linear-gradient(135deg, #10b981 0%, #059669 100%)'
        else:
            primary_value = 'Not Determined'
            primary_label = 'Status'
            primary_bg = 'linear-gradient(135deg, #6b7280 0%, #4b5563 100%)'

        gene_count = len(result['genes_found'])
        confidence = result['confidence']

        failure_card = ''
        if not is_success:
            failure_card = f"""
        <div class="report-section" style="border-left: 4px solid #dc2626;">
            <h2 style="color: #dc2626;">Analysis Error</h2>
            <p>{result['status']}</p>
        </div>
            """

        warnings_html = ''
        if result.get('warnings'):
            items = ''.join(f'<li>{w}</li>' for w in result['warnings'])
            warnings_html = f'<div class="warning-box"><strong>Warnings:</strong><ul>{items}</ul></div>'

        detailed = result.get('detailed_data', {})
        o_type_html = self._render_gene_table_html(detailed.get('O_type', {}), 'O-type')
        h_type_html = self._render_gene_table_html(detailed.get('H_type', {}), 'H-type')

        gene_badges = ''
        if result['genes_found']:
            gene_badges = ' '.join(
                f'<span class="marker-badge marker-info">{g}</span>'
                for g in result['genes_found']
            )
        else:
            gene_badges = '<span style="color:#6b7280;">None detected</span>'

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliTyper - {sample_name} Serotype Report</title>
    <style>{self._get_shared_css()}</style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">&mdash; {random_quote['author']}</div>
            </div>
        </div>

        <div class="report-section">
            <h2>Sample Information</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Sample Name</div>
                    <div class="metric-value" style="font-size:16px;">{sample_name}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Analysis Date</div>
                    <div class="metric-value" style="font-size:14px;">{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Status</div>
                    <div class="metric-value" style="font-size:16px;">{result['status']}</div>
                </div>
            </div>
        </div>

        {failure_card}

        <div class="report-section">
            <h2>Primary Serotype Result</h2>
            <div class="primary-result" style="background: {primary_bg};">
                <div class="primary-label">{primary_label}</div>
                <div class="primary-value">{primary_value}</div>
            </div>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">O-type</div>
                    <div class="metric-value">{result['o_type']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">H-type</div>
                    <div class="metric-value">{result['h_type']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Genes Found</div>
                    <div class="metric-value">{gene_count}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Confidence</div>
                    <div class="metric-value" style="font-size:18px;">{confidence}</div>
                </div>
            </div>
        </div>

        <div class="report-section">
            <h2>Detected Serotype Genes</h2>
            <div style="margin: 10px 0;">{gene_badges}</div>
            <div class="info-box">
                <strong>Method:</strong> SerotypeFinder (Joensen et al.) predicts O and H antigens
                by matching genome sequences against curated reference gene databases at &ge;90%
                identity and &ge;60% coverage.
            </div>
            {warnings_html}
        </div>

        <div class="report-section">
            <h2>O-type Gene Details</h2>
            {o_type_html}
        </div>

        <div class="report-section">
            <h2>H-type Gene Details</h2>
            {h_type_html}
        </div>

        {self._get_footer_html("Serotype Per-Sample Report")}
    </div>
    {self._get_quote_js()}
</body>
</html>"""

        html_file = output_dir / f"{sample_name}_serotype_report.html"
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

    # ------------------------------------------------------------------
    # Batch summaries
    # ------------------------------------------------------------------

    def generate_json_report(self, output_dir: Path) -> str:
        """Generate structured JSON report with comprehensive analysis data."""
        json_report = {
            "metadata": self.metadata,
            "analysis_summary": {
                "total_samples": len(self.results),
                "successful_analyses": len([r for r in self.results if r['status'] == 'Completed']),
                "failed_analyses": len([r for r in self.results if r['status'] != 'Completed']),
                "unique_serotypes_found": len(set(
                    r['serotype'] for r in self.results if r['serotype'] != "Unknown"
                )),
                "unique_o_types": len(set(
                    r['o_type'] for r in self.results if r['o_type'] != "Unknown"
                )),
                "unique_h_types": len(set(
                    r['h_type'] for r in self.results if r['h_type'] != "Unknown"
                )),
            },
            "serotype_distribution": self._calculate_serotype_distribution(),
            "detailed_results": self.results,
            "run_parameters": {
                "database_path": str(self.db_path),
                "threads": self.threads,
                "analysis_timestamp": datetime.now().isoformat(),
            },
        }

        # Add per_sample_report path per sample for convenience
        for r in json_report["detailed_results"]:
            r["per_sample_report"] = (
                f"{r['sample_id']}/{r['sample_id']}_serotype_report.html"
            )

        json_file = output_dir / "serotype_analysis_report.json"
        with open(json_file, 'w') as f:
            json.dump(json_report, f, indent=4, default=str)

        return str(json_file)

    def _calculate_serotype_distribution(self) -> Dict[str, Any]:
        """Calculate distribution statistics for serotypes."""
        distribution = {
            "serotype_counts": {},
            "o_type_counts": {},
            "h_type_counts": {},
            "gene_frequency": {},
        }

        for result in self.results:
            if result['status'] == 'Completed':
                serotype = result['serotype']
                o_type = result['o_type']
                h_type = result['h_type']

                distribution["serotype_counts"][serotype] = (
                    distribution["serotype_counts"].get(serotype, 0) + 1
                )

                if o_type != "Unknown":
                    distribution["o_type_counts"][o_type] = (
                        distribution["o_type_counts"].get(o_type, 0) + 1
                    )

                if h_type != "Unknown":
                    distribution["h_type_counts"][h_type] = (
                        distribution["h_type_counts"].get(h_type, 0) + 1
                    )

                for gene in result['genes_found']:
                    distribution["gene_frequency"][gene] = (
                        distribution["gene_frequency"].get(gene, 0) + 1
                    )

        for key in ["serotype_counts", "o_type_counts", "h_type_counts", "gene_frequency"]:
            distribution[key] = dict(sorted(
                distribution[key].items(), key=lambda x: x[1], reverse=True,
            ))

        return distribution

    def generate_html_report(self, output_dir: Path) -> str:
        """Generate the batch HTML summary with links to per-sample reports."""
        random_quote = random.choice(self.science_quotes)

        successful = len([r for r in self.results if r['status'] == 'Completed'])
        failed = len([r for r in self.results if r['status'] != 'Completed'])
        unique_serotypes = len(set(
            r['serotype'] for r in self.results if r['serotype'] != "Unknown"
        ))

        serotype_counts = {}
        for r in self.results:
            if r['status'] == 'Completed' and r['serotype'] != 'Unknown':
                serotype_counts[r['serotype']] = serotype_counts.get(r['serotype'], 0) + 1

        serotype_badges_html = ''.join(
            f'<span class="marker-badge marker-info" style="margin:4px;">{s} ({c})</span>'
            for s, c in sorted(serotype_counts.items())
        )
        if not serotype_badges_html:
            serotype_badges_html = '<span style="color:#6b7280;">None detected</span>'

        rows = ''
        for result in self.results:
            status_class = 'marker-present' if result['status'] == 'Completed' else 'marker-absent'
            link = f"{result['sample_id']}/{result['sample_id']}_serotype_report.html"
            genes_display = ', '.join(result['genes_found']) if result['genes_found'] else '&mdash;'
            rows += f"""
                        <tr>
                            <td><strong><a href="{link}">{result['sample_id']}</a></strong></td>
                            <td><strong style="color: #2563eb;">{result['serotype']}</strong></td>
                            <td>{result['o_type']}</td>
                            <td>{result['h_type']}</td>
                            <td class="product-cell">{genes_display}</td>
                            <td><span class="marker-badge {status_class}">{result['status']}</span></td>
                        </tr>
            """

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliTyper - Serotype Batch Summary</title>
    <style>{self._get_shared_css()}</style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">&mdash; {random_quote['author']}</div>
            </div>
        </div>

        <div class="report-section">
            <h2>Batch Serotyping Summary</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{len(self.results)}</div>
                    <div class="stat-label">SAMPLES PROCESSED</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{successful}</div>
                    <div class="stat-label">SUCCESSFUL</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{failed}</div>
                    <div class="stat-label">FAILED</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{unique_serotypes}</div>
                    <div class="stat-label">UNIQUE SEROTYPES</div>
                </div>
            </div>
            <p style="margin-top:15px;"><strong>Analysis Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>Database:</strong> {self.db_path}</p>
        </div>

        <div class="report-section">
            <h2>Serotype Distribution</h2>
            <div style="margin: 15px 0;">{serotype_badges_html}</div>
            <div class="info-box">
                <strong>Method:</strong> SerotypeFinder (Joensen et al., 2015) predicts E. coli O and H
                antigens by matching genome sequences against curated reference gene databases.
                Click any sample below to open its detailed per-sample report.
            </div>
        </div>

        <div class="report-section">
            <h2>Per-Sample Results</h2>
            <div class="table-responsive">
                <table class="summary-table">
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
                        {rows}
                    </tbody>
                </table>
            </div>
        </div>

        {self._get_footer_html("Serotype Batch Summary")}
    </div>
    {self._get_quote_js()}
</body>
</html>"""

        html_file = output_dir / "serotype_analysis_report.html"
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        return str(html_file)

    def generate_tsv_report(self, output_dir: Path) -> str:
        """Generate TSV report with all sample results."""
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
                'File_Path': result['file_path'],
            })

        df = pd.DataFrame(data)
        tsv_file = output_dir / "serotype_analysis_report.tsv"
        df.to_csv(tsv_file, sep='\t', index=False)
        return str(tsv_file)


def main():
    parser = argparse.ArgumentParser(
        description="EcoliTyper Enhanced SerotypeFinder - Batch processing with multi-format reports",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Per-sample output (1 file each):
  <sample>/<sample>_serotype_report.html

Batch-level outputs (at the top of -o):
  serotype_analysis_report.html
  serotype_analysis_report.tsv
  serotype_analysis_report.json
        """,
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input FASTA file, directory, or glob pattern (e.g., "*.fna", "genomes/*.fasta")')
    parser.add_argument('-db', '--database', default='serotypefinder_db',
                        help='Path to serotypefinder database')
    parser.add_argument('-o', '--output', required=True,
                        help='Main output directory')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of threads to use')

    args = parser.parse_args()

    main_output_dir = Path(args.output)
    main_output_dir.mkdir(parents=True, exist_ok=True)

    try:
        finder = EnhancedSerotypeFinder(args.database, args.threads)

        print(finder.ascii_art)
        print("EcoliTyper Enhanced SerotypeFinder")
        print("=" * 50)
        print(f"Input: {args.input}")
        print(f"Output: {main_output_dir}")
        print(f"Threads: {args.threads}")
        print("=" * 50)

        results = finder.process_batch(args.input, main_output_dir)

        print("\nGenerating reports...")
        json_file = finder.generate_json_report(main_output_dir)
        html_file = finder.generate_html_report(main_output_dir)
        tsv_file = finder.generate_tsv_report(main_output_dir)

        print("\nAnalysis Complete.")
        print(f"Samples processed: {len(results)}")
        print(f"Results directory: {main_output_dir}")
        print("\nPer-sample outputs (1 file each):")
        print("  <sample>/<sample>_serotype_report.html")
        print("\nBatch-level outputs:")
        print(f"  {html_file}")
        print(f"  {tsv_file}")
        print(f"  {json_file}")

        successful = len([r for r in results if r['status'] == 'Completed'])
        if len(results) > 0:
            print(f"\nSuccess rate: {successful}/{len(results)} ({successful/len(results)*100:.1f}%)")

        serotypes = list(set(
            r['serotype'] for r in results
            if r['status'] == 'Completed' and r['serotype'] != 'Unknown'
        ))
        print(f"Serotypes found: {', '.join(serotypes) if serotypes else 'None'}")

        quote = random.choice(finder.science_quotes)
        print(f"\n{quote['text']} - {quote['author']}")

    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()