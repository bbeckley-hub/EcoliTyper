#!/usr/bin/env python3
"""
EcoliTyper - Enhanced CHTyper Wrapper
Comprehensive E. coli typing with FumC and FimH with auto BLAST detection

Per-sample output (1 file each):
  - {sample}_chtyper_report.html   (FumC + FimH typing, metrics, detailed output)

Batch-level outputs:
  - chtyper_results.html
  - chtyper_results.tsv
  - chtyper_summary.json

Email: <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
Date: 2025 / Updated 2026-07-22
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
import html
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any
import pandas as pd
import multiprocessing as mp
from collections import Counter


class EnhancedCHTyper:
    """CHTyper wrapper with auto BLAST detection and single-file per-sample reports."""

    def __init__(self, db_path: str = "chtyper_db", threads: int = 4):
        self.db_path = Path(db_path)
        self.threads = threads
        self.results = []
        self.metadata = {
            "tool_name": "EcoliTyper CHTyper",
            "version": "1.4.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "citation": (
                "Roer L, Tchesnokova V, Allesøe R, et al. Development of a Web Tool for "
                "Escherichia coli Subtyping Based on FumC and FimH (CH Typing). "
                "Front Microbiol. 2017;8:2511."
            ),
        }

        self.blast_path = self._detect_blast_path()
        if not self.blast_path:
            raise RuntimeError("BLAST not found. Please install BLAST or ensure it's in your PATH")

        self.science_quotes = [
            {"text": "The good thing about science is that it's true whether or not you believe in it.",
             "author": "Neil deGrasse Tyson"},
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.",
             "author": "Albert Einstein"},
            {"text": "Nothing in life is to be feared, it is only to be understood.",
             "author": "Marie Curie"},
            {"text": "The microscope opens a new world to the investigator.",
             "author": "Robert Koch"},
            {"text": "In science, the credit goes to the man who convinces the world, not to the man to whom the idea first occurs.",
             "author": "Francis Darwin"},
            {"text": "Science knows no country, because knowledge belongs to humanity.",
             "author": "Louis Pasteur"},
            {"text": "Somewhere, something incredible is waiting to be known.",
             "author": "Carl Sagan"},
            {"text": "Research is what I'm doing when I don't know what I'm doing.",
             "author": "Wernher von Braun"},
        ]

        self.ascii_art = """
███████╗ ██████╗ ██████╗ ██╗     ██╗████████╗██╗   ██╗██████╗ ███████╗██████╗ 
██╔════╝██╔════╝██╔═══██╗██║     ██║╚══██╔══╝╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
█████╗  ██║     ██║   ██║██║     ██║   ██║    ╚████╔╝ ██████╔╝█████╗  ██████╔╝
██╔══╝  ██║     ██║   ██║██║     ██║   ██║     ╚██╔╝  ██╔═══╝ ██╔══╝  ██╔══██╗
███████╗╚██████╗╚██████╔╝███████╗██║   ██║      ██║   ██║     ███████╗██║  ██║
╚══════╝ ╚═════╝ ╚═════╝ ╚══════╝╚═╝   ╚═╝      ╚═╝   ╚═╝     ╚══════╝╚═╝  ╚═╝
        """

    def _detect_blast_path(self) -> str:
        """Auto-detect BLAST installation."""
        try:
            result = subprocess.run(["which", "blastn"], capture_output=True, text=True)
            if result.returncode == 0:
                return result.stdout.strip()
        except Exception:
            pass
        return ""

    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files matching a path, directory, or glob pattern."""
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

    def run_chtyper_analysis(self, fasta_file: Path, output_base: Path) -> Dict[str, Any]:
        """Run CHTyper on one genome. Emits exactly 1 file per sample:
             - {sample}_chtyper_report.html
        All CHTyper intermediates live in a TemporaryDirectory that self-destructs.
        """
        sample_name = fasta_file.stem
        try:
            print(f"Analyzing {sample_name}...")

            chyper_script = "CHTyper-fixed.py"
            if not Path(chyper_script).exists():
                parsed = self._create_error_result(
                    sample_name, str(fasta_file), "CHTyper script not found"
                )
            else:
                with tempfile.TemporaryDirectory() as temp_dir:
                    temp_path = Path(temp_dir)

                    cmd = [
                        sys.executable, chyper_script,
                        "-i", str(fasta_file),
                        "-o", str(temp_path),
                        "-p", str(self.db_path),
                        "-b", self.blast_path,
                        "-l", "0.6",
                        "-t", "0.9",
                    ]

                    result = subprocess.run(cmd, capture_output=True, text=True)

                    if result.returncode != 0:
                        error_msg = f"Command failed with return code {result.returncode}"
                        if result.stderr:
                            error_msg += f": {result.stderr}"
                        parsed = self._create_error_result(
                            sample_name, str(fasta_file), error_msg
                        )
                    else:
                        parsed = self._parse_sample_results_from_temp(
                            sample_name, str(fasta_file), temp_path
                        )

            sample_output_dir = output_base / sample_name
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            self._generate_per_sample_html(parsed, sample_output_dir)
            parsed['output_directory'] = str(sample_output_dir)

            return parsed

        except Exception as e:
            try:
                error_result = self._create_error_result(
                    sample_name, str(fasta_file), str(e)
                )
                sample_output_dir = output_base / sample_name
                sample_output_dir.mkdir(parents=True, exist_ok=True)
                self._generate_per_sample_html(error_result, sample_output_dir)
                error_result['output_directory'] = str(sample_output_dir)
                return error_result
            except Exception:
                return self._create_error_result(sample_name, str(fasta_file), str(e))

    def _parse_sample_results_from_temp(self, sample_name: str, fasta_path: str,
                                        temp_dir: Path) -> Dict[str, Any]:
        """Parse CHTyper output files from the temporary working directory."""
        try:
            tab_file = temp_dir / "results_tab.txt"
            if not tab_file.exists():
                return self._create_error_result(
                    sample_name, fasta_path, "No results_tab.txt file found"
                )

            with open(tab_file, 'r') as f:
                content = f.read()

            fumc_type = "Unknown"
            fimh_type = "Unknown"
            fumc_identity = "0.00"
            fimh_identity = "0.00"
            fumc_coverage = "0/0"
            fimh_coverage = "0/0"

            for line in content.strip().split('\n'):
                if line.startswith('fumC') and '\t' in line:
                    parts = line.split('\t')
                    if len(parts) >= 6:
                        fumc_type = parts[0]
                        fumc_identity = parts[1]
                        fumc_coverage = parts[2]
                elif line.startswith('fimH') and '\t' in line:
                    parts = line.split('\t')
                    if len(parts) >= 6:
                        fimh_type = parts[0]
                        fimh_identity = parts[1]
                        fimh_coverage = parts[2]

            detailed_results = ""
            txt_file = temp_dir / "results.txt"
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
                "output_directory": "",
                "warnings": [],
                "detailed_results": detailed_results,
            }
        except Exception as e:
            return self._create_error_result(
                sample_name, fasta_path, f"Error parsing results: {str(e)}"
            )

    def _create_error_result(self, sample_name: str, fasta_path: str,
                             error_msg: str) -> Dict[str, Any]:
        """Return a placeholder result when a sample fails."""
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
            "detailed_results": "",
        }

    def _get_shared_css(self) -> str:
        """Return the shared CSS block used by per-sample and batch HTML reports."""
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
        .primary-result {
            color: white; padding: 30px; border-radius: 12px;
            margin: 20px 0; box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            text-align: center;
        }
        .primary-label {
            font-size: 13px; opacity: 0.9; letter-spacing: 3px;
            margin-bottom: 12px; text-transform: uppercase;
        }
        .primary-value { font-size: 48px; font-weight: bold; margin-bottom: 15px; }
        .metrics-grid {
            display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px; margin-top: 20px;
        }
        .metric-card {
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white; padding: 20px; border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }
        .metric-card.fumc-card {
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
        }
        .metric-card.fimh-card {
            background: linear-gradient(135deg, #10b981 0%, #059669 100%);
        }
        .metric-label { font-size: 13px; opacity: 0.9; margin-bottom: 5px; }
        .metric-value { font-size: 22px; font-weight: bold; }
        .confidence-badge {
            display: inline-block; padding: 6px 18px; border-radius: 20px;
            font-size: 14px; font-weight: bold; text-transform: uppercase;
            letter-spacing: 1px;
        }
        .confidence-high { background-color: #16a34a; color: white; }
        .confidence-medium { background-color: #f59e0b; color: white; }
        .confidence-low { background-color: #dc2626; color: white; }
        .confidence-failed { background-color: #6b7280; color: white; }
        .detailed-output {
            background: #f8fafc; padding: 20px; border-radius: 8px;
            font-family: 'Courier New', monospace; font-size: 13px;
            color: #1f2937; overflow-x: auto; white-space: pre-wrap;
            border: 1px solid #e5e7eb; max-height: 600px; overflow-y: auto;
        }
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
        .type-badge {
            display: inline-block; padding: 4px 12px; border-radius: 20px;
            font-weight: bold; font-size: 13px; background: #3b82f6; color: white;
            margin: 2px;
        }
        .success { color: #16a34a; font-weight: bold; }
        .error { color: #dc2626; font-weight: bold; }
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
            .primary-value { font-size: 28px; }
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
                    Powered by CHTyper and BLAST (Camacho et al., 2009).<br>
                    Coverage threshold: 0.6 | Identity threshold: 0.9<br>
                    {self.metadata['citation']}
                </p>
            </div>
        </div>
        """

    def _generate_per_sample_html(self, result: Dict[str, Any], output_dir: Path):
        """Write the single per-sample CHTyper HTML report."""
        sample_name = result['sample_id']
        random_quote = random.choice(self.science_quotes)
        is_success = result['status'] == 'Completed'

        fumc_type = result['fumc_type']
        fimh_type = result['fimh_type']

        try:
            fumc_id_f = float(result['fumc_identity'])
            fimh_id_f = float(result['fimh_identity'])
        except (ValueError, TypeError):
            fumc_id_f = fimh_id_f = 0.0

        if not is_success:
            confidence = "Failed"
            confidence_class = "confidence-failed"
        elif fumc_id_f >= 99.0 and fimh_id_f >= 99.0:
            confidence = "High"
            confidence_class = "confidence-high"
        elif fumc_id_f >= 90.0 and fimh_id_f >= 90.0:
            confidence = "Medium"
            confidence_class = "confidence-medium"
        else:
            confidence = "Low"
            confidence_class = "confidence-low"

        if is_success:
            primary_value = f"{fumc_type}-{fimh_type}"
            primary_label = "CH Type"
            primary_bg = "linear-gradient(135deg, #10b981 0%, #059669 100%)"
        else:
            primary_value = "Analysis Failed"
            primary_label = "Status"
            primary_bg = "linear-gradient(135deg, #dc2626 0%, #b91c1c 100%)"

        failure_card = ""
        if not is_success:
            failure_card = f"""
        <div class="report-section" style="border-left: 4px solid #dc2626;">
            <h2 style="color: #dc2626;">Analysis Error</h2>
            <p>{html.escape(result['status'])}</p>
        </div>
            """

        detailed_card = ""
        if result.get('detailed_results'):
            detailed_card = f"""
        <div class="report-section">
            <h2>Detailed CHTyper Output</h2>
            <pre class="detailed-output">{html.escape(result['detailed_results'])}</pre>
        </div>
            """

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliTyper CHTyper - Report - {sample_name}</title>
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
                    <div class="metric-value">{sample_name}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Analysis Date</div>
                    <div class="metric-value">{datetime.now().strftime('%Y-%m-%d')}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Confidence</div>
                    <div class="metric-value">
                        <span class="confidence-badge {confidence_class}">{confidence}</span>
                    </div>
                </div>
            </div>
        </div>

        {failure_card}

        <div class="report-section">
            <h2>CH Typing Result</h2>
            <div class="primary-result" style="background: {primary_bg};">
                <div class="primary-label">{primary_label}</div>
                <div class="primary-value">{primary_value}</div>
                <div class="confidence-badge {confidence_class}">{confidence} Confidence</div>
            </div>

            <div class="metrics-grid">
                <div class="metric-card fumc-card">
                    <div class="metric-label">FumC Type</div>
                    <div class="metric-value">{fumc_type}</div>
                </div>
                <div class="metric-card fumc-card">
                    <div class="metric-label">FumC Identity</div>
                    <div class="metric-value">{result['fumc_identity']}%</div>
                </div>
                <div class="metric-card fumc-card">
                    <div class="metric-label">FumC Coverage</div>
                    <div class="metric-value">{result['fumc_coverage']}</div>
                </div>
                <div class="metric-card fimh-card">
                    <div class="metric-label">FimH Type</div>
                    <div class="metric-value">{fimh_type}</div>
                </div>
                <div class="metric-card fimh-card">
                    <div class="metric-label">FimH Identity</div>
                    <div class="metric-value">{result['fimh_identity']}%</div>
                </div>
                <div class="metric-card fimh-card">
                    <div class="metric-label">FimH Coverage</div>
                    <div class="metric-value">{result['fimh_coverage']}</div>
                </div>
            </div>
        </div>

        {detailed_card}

        {self._get_footer_html("CHTyper Per-Sample Report")}
    </div>
    {self._get_quote_js()}
</body>
</html>"""

        html_file = output_dir / f"{sample_name}_chtyper_report.html"
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

    def process_batch(self, input_path: str, main_output_dir: Path) -> List[Dict[str, Any]]:
        """Process all FASTA files in batch."""
        fasta_files = self.find_fasta_files(input_path)

        print(f"Processing {len(fasta_files)} samples using {self.threads} threads...")
        print(f"Using BLAST at: {self.blast_path}")

        args = [(fasta_file, main_output_dir) for fasta_file in fasta_files]

        with mp.Pool(processes=self.threads) as pool:
            results = pool.starmap(self.run_chtyper_analysis, args)

        self.results = results
        return results

    def generate_json_report(self, output_dir: Path) -> str:
        """Write the batch CHTyper JSON summary."""
        total_samples = len(self.results)
        successful_samples = len([r for r in self.results if r['status'] == 'Completed'])

        def parse_identity(identity_str):
            try:
                return float(str(identity_str).replace('%', ''))
            except Exception:
                return 0.0

        def parse_coverage(coverage_str):
            try:
                if '/' in str(coverage_str):
                    num, denom = str(coverage_str).split('/')
                    return float(num) / float(denom) * 100
                return 0.0
            except Exception:
                return 0.0

        fumc_types = Counter()
        fimh_types = Counter()
        for result in self.results:
            if result['status'] == 'Completed':
                if result['fumc_type'] != 'Unknown':
                    fumc_types[result['fumc_type']] += 1
                if result['fimh_type'] != 'Unknown':
                    fimh_types[result['fimh_type']] += 1

        fumc_distribution = dict(sorted(fumc_types.items(), key=lambda x: x[1], reverse=True))
        fimh_distribution = dict(sorted(fimh_types.items(), key=lambda x: x[1], reverse=True))

        fumc_identities, fimh_identities = [], []
        fumc_coverages, fimh_coverages = [], []
        for result in self.results:
            if result['status'] == 'Completed' and result['fumc_type'] != 'Unknown':
                fumc_identities.append(parse_identity(result['fumc_identity']))
                fumc_coverages.append(parse_coverage(result['fumc_coverage']))
            if result['status'] == 'Completed' and result['fimh_type'] != 'Unknown':
                fimh_identities.append(parse_identity(result['fimh_identity']))
                fimh_coverages.append(parse_coverage(result['fimh_coverage']))

        avg = lambda x: sum(x) / len(x) if x else 0.0

        detailed_results = []
        for result in self.results:
            detailed_results.append({
                "sample_id": result["sample_id"],
                "fumc_type": result["fumc_type"],
                "fimh_type": result["fimh_type"],
                "fumc_identity": result["fumc_identity"],
                "fimh_identity": result["fimh_identity"],
                "fumc_coverage": result["fumc_coverage"],
                "fimh_coverage": result["fimh_coverage"],
                "status": result["status"],
                "file_path": result["file_path"],
                "per_sample_report": (
                    f"{result['sample_id']}/{result['sample_id']}_chtyper_report.html"
                ),
                "has_detailed_results": bool(result.get("detailed_results", "")),
            })

        json_report = {
            "metadata": {
                **self.metadata,
                "analysis_type": "CHTyper (FumC/FimH Typing)",
                "blast_path": self.blast_path,
                "database_path": str(self.db_path),
                "threads_used": self.threads,
                "analysis_parameters": {
                    "minimum_coverage": 0.6,
                    "identity_threshold": 0.9,
                },
            },
            "analysis_summary": {
                "total_samples": total_samples,
                "successful_analyses": successful_samples,
                "failed_analyses": total_samples - successful_samples,
                "success_rate": round(successful_samples / total_samples * 100, 2) if total_samples > 0 else 0,
                "unique_fumc_types": len(fumc_distribution),
                "unique_fimh_types": len(fimh_distribution),
                "most_common_fumc_type": list(fumc_distribution.items())[0][0] if fumc_distribution else "None",
                "most_common_fimh_type": list(fimh_distribution.items())[0][0] if fimh_distribution else "None",
                "quality_metrics": {
                    "average_fumc_identity": round(avg(fumc_identities), 2),
                    "average_fimh_identity": round(avg(fimh_identities), 2),
                    "average_fumc_coverage": round(avg(fumc_coverages), 2),
                    "average_fimh_coverage": round(avg(fimh_coverages), 2),
                },
            },
            "type_distributions": {
                "fumc_types": fumc_distribution,
                "fimh_types": fimh_distribution,
            },
            "methodology": {
                "description": "CHTyper analysis for E. coli using FumC and FimH genes",
                "genes_analyzed": ["fumC", "fimH"],
                "purpose": "FumC-FimH typing for E. coli subtyping and epidemiological studies",
                "algorithm": "BLAST-based sequence typing",
                "database": "CHTyper database containing FumC and FimH allele sequences",
                "reference": self.metadata['citation'],
            },
            "detailed_results": detailed_results,
            "performance_statistics": {
                "high_confidence_results": len([
                    r for r in self.results
                    if parse_identity(r['fumc_identity']) >= 99.0
                    and parse_identity(r['fimh_identity']) >= 99.0
                ]),
                "medium_confidence_results": len([
                    r for r in self.results
                    if 90.0 <= parse_identity(r['fumc_identity']) < 99.0
                    and 90.0 <= parse_identity(r['fimh_identity']) < 99.0
                ]),
                "low_confidence_results": len([
                    r for r in self.results
                    if parse_identity(r['fumc_identity']) < 90.0
                    or parse_identity(r['fimh_identity']) < 90.0
                ]),
            },
            "combined_typing_summary": {
                "unique_combinations": len(set([
                    (r['fumc_type'], r['fimh_type'])
                    for r in self.results
                    if r['status'] == 'Completed'
                    and r['fumc_type'] != 'Unknown'
                    and r['fimh_type'] != 'Unknown'
                ])),
                "most_common_combination": (
                    max(
                        [(r['fumc_type'], r['fimh_type']) for r in self.results
                         if r['status'] == 'Completed'
                         and r['fumc_type'] != 'Unknown'
                         and r['fimh_type'] != 'Unknown'],
                        key=lambda x: list((r['fumc_type'], r['fimh_type'])
                                           for r in self.results).count(x),
                        default=("None", "None"),
                    )
                ),
            },
        }

        json_file = output_dir / "chtyper_summary.json"
        with open(json_file, 'w') as f:
            json.dump(json_report, f, indent=4, default=str)
        return str(json_file)

    def generate_html_report(self, output_dir: Path) -> str:
        """Write the batch CHTyper HTML summary with links to per-sample reports."""
        random_quote = random.choice(self.science_quotes)
        successful = len([r for r in self.results if r['status'] == 'Completed'])
        failed = len(self.results) - successful

        fumc_types = Counter()
        fimh_types = Counter()
        for result in self.results:
            if result['status'] == 'Completed':
                if result['fumc_type'] != 'Unknown':
                    fumc_types[result['fumc_type']] += 1
                if result['fimh_type'] != 'Unknown':
                    fimh_types[result['fimh_type']] += 1

        fumc_badges = "".join(
            f'<span class="type-badge">{t} ({c})</span>'
            for t, c in sorted(fumc_types.items())
        ) or '<span style="color: #6b7280;">None detected</span>'
        fimh_badges = "".join(
            f'<span class="type-badge">{t} ({c})</span>'
            for t, c in sorted(fimh_types.items())
        ) or '<span style="color: #6b7280;">None detected</span>'

        table_rows = ""
        for result in self.results:
            status_class = "success" if result["status"] == "Completed" else "error"
            link = f"{result['sample_id']}/{result['sample_id']}_chtyper_report.html"
            table_rows += f"""
                        <tr>
                            <td><a href="{link}"><strong>{result['sample_id']}</strong></a></td>
                            <td><strong style="color: #2563eb;">{result['fumc_type']}</strong></td>
                            <td>{result['fumc_identity']}%</td>
                            <td>{result['fumc_coverage']}</td>
                            <td><strong style="color: #059669;">{result['fimh_type']}</strong></td>
                            <td>{result['fimh_identity']}%</td>
                            <td>{result['fimh_coverage']}</td>
                            <td class="{status_class}">{result['status']}</td>
                        </tr>
            """

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliTyper CHTyper - Batch Summary</title>
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
            <h2>Batch CHTyper Summary</h2>
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
                    <div class="stat-value">{len(fumc_types)}</div>
                    <div class="stat-label">UNIQUE FUMC TYPES</div>
                </div>
            </div>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>Threads Used:</strong> {self.threads}</p>
            <p><strong>BLAST Path:</strong> {self.blast_path}</p>
        </div>

        <div class="report-section">
            <h2>Type Distribution</h2>
            <p><strong>FumC Types:</strong></p>
            <div style="margin: 10px 0 20px 0;">{fumc_badges}</div>
            <p><strong>FimH Types:</strong></p>
            <div style="margin: 10px 0;">{fimh_badges}</div>
        </div>

        <div class="report-section">
            <h2>CHTyper Results</h2>
            <div class="table-responsive">
                <table class="summary-table">
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
                    <tbody>{table_rows}</tbody>
                </table>
            </div>
        </div>

        {self._get_footer_html("CHTyper Batch Summary Report")}
    </div>
    {self._get_quote_js()}
</body>
</html>"""

        html_file = output_dir / "chtyper_results.html"
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        return str(html_file)

    def generate_tsv_report(self, output_dir: Path) -> str:
        """Write the CHTyper TSV summary."""
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
                'File_Path': result['file_path'],
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

    main_output_dir = Path(args.output)
    main_output_dir.mkdir(parents=True, exist_ok=True)

    try:
        finder = EnhancedCHTyper(args.database, args.threads)

        print(finder.ascii_art)
        print("EcoliTyper Enhanced CHTyper")
        print("=" * 50)
        print(f"Input: {args.input}")
        print(f"Output: {main_output_dir}")
        print(f"Threads: {args.threads}")
        print(f"Auto-detected BLAST: {finder.blast_path}")
        print("=" * 50)

        results = finder.process_batch(args.input, main_output_dir)

        print("\nGenerating reports...")
        json_file = finder.generate_json_report(main_output_dir)
        html_file = finder.generate_html_report(main_output_dir)
        tsv_file = finder.generate_tsv_report(main_output_dir)

        print("\nAnalysis complete.")
        print(f"Samples processed: {len(results)}")
        print(f"Results directory: {main_output_dir}")
        print("\nPer-sample outputs (1 file each):")
        print("  - {sample}/{sample}_chtyper_report.html")
        print("\nBatch-level outputs:")
        print(f"  - {json_file}")
        print(f"  - {html_file}")
        print(f"  - {tsv_file}")

        successful = len([r for r in results if r['status'] == 'Completed'])
        print(f"\nSuccess rate: {successful}/{len(results)} ({successful / len(results) * 100:.1f}%)")

        fumc_types = list(set([r['fumc_type'] for r in results
                               if r['status'] == 'Completed'
                               and r['fumc_type'] != 'Unknown']))
        fimh_types = list(set([r['fimh_type'] for r in results
                               if r['status'] == 'Completed'
                               and r['fimh_type'] != 'Unknown']))
        print(f"FumC types found: {', '.join(fumc_types) if fumc_types else 'None'}")
        print(f"FimH types found: {', '.join(fimh_types) if fimh_types else 'None'}")

        quote = random.choice(finder.science_quotes)
        print(f"\n{quote['text']} - {quote['author']}")

    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()