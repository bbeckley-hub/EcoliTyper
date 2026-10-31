#!/usr/bin/env python3
"""
EcoliTyper Module for MLST Analysis of E. coli
Comprehensive MLST analysis with beautiful HTML reporting.

Per-sample output (1 file each):
    <sample>/<sample>_mlst_report.html

Batch-level outputs (at the top of -o):
    mlst_summary.html
    mlst_summary.tsv
    mlst_summary.json

Author: Beckley Brown <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
Date: 2025-12-16 / 2026-07-22
Send a quick mail for any issues or further explanations.
"""

import os
import sys
import json
import glob
import shutil
import tempfile
import argparse
import subprocess
import random
import html as html_module
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
from datetime import datetime
from collections import Counter


class EcoliTyper:
    """EcoliTyper MLST analyzer with per-sample HTML reporting."""

    def __init__(self, database_dir: Path, script_dir: Path):
        # Resolve both to absolute paths. The subprocess changes cwd to a
        # TemporaryDirectory, so any relative path (mlst binary, MLST_DATADIR)
        self.database_dir = Path(database_dir).resolve()
        self.script_dir = Path(script_dir).resolve()
        self.mlst_bin = self.script_dir / "mlst"

        # Fail fast with a clear message instead of 15 opaque STUNKNOWN results.
        if not self.mlst_bin.exists():
            raise FileNotFoundError(
                f"mlst binary not found at {self.mlst_bin}\n"
                f"Check that -sc/--script-dir points to a directory containing 'mlst'."
            )
        if not self.database_dir.exists():
            raise FileNotFoundError(
                f"MLST database directory not found at {self.database_dir}\n"
                f"Check that -db/--database-dir points to the mlst database root."
            )

        self.ascii_art = """
███████╗ ██████╗ ██████╗ ██╗     ██╗████████╗██╗   ██╗██████╗ ███████╗██████╗ 
██╔════╝██╔════╝██╔═══██╗██║     ██║╚══██╔══╝╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
█████╗  ██║     ██║   ██║██║     ██║   ██║    ╚████╔╝ ██████╔╝█████╗  ██████╔╝
██╔══╝  ██║     ██║   ██║██║     ██║   ██║     ╚██╔╝  ██╔═══╝ ██╔══╝  ██╔══██╗
███████╗╚██████╗╚██████╔╝███████╗██║   ██║      ██║   ██║     ███████╗██║  ██║
╚══════╝ ╚═════╝ ╚═════╝ ╚══════╝╚═╝   ╚═╝      ╚═╝   ╚═╝     ╚══════╝╚═╝  ╚═╝
        """

        self.metadata = {
            "tool_name": "EcoliTyper MLST Analysis",
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
        ]

    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files. Handles three input shapes:
          - a single file path
          - a directory (all matching files inside)
          - a glob pattern (e.g. \"*.fna\") — must be quoted by the shell
        """
        if os.path.isfile(input_path):
            return [Path(input_path)]

        fasta_patterns = [
            input_path,
            f"{input_path}/*.fna", f"{input_path}/*.fasta",
            f"{input_path}/*.fa", f"{input_path}/*.fn",
            f"{input_path}/*.fna.gz", f"{input_path}/*.fasta.gz",
            f"{input_path}/*.fa.gz", f"{input_path}/*.gb",
            f"{input_path}/*.gbk", f"{input_path}/*.gbff",
        ]

        fasta_files = []
        for pattern in fasta_patterns:
            matched_files = glob.glob(pattern)
            for file_path in matched_files:
                path = Path(file_path)
                if path.is_file():
                    fasta_files.append(path)

        return sorted(list(set(fasta_files)))

    def run_mlst_single(self, input_file: Path, output_dir: Path,
                        scheme: str = "ecoli_achtman_4") -> Dict:
        """Run MLST on one genome. Emits exactly 1 file per sample:
             - <sample>/<sample>_mlst_report.html
        mlst --csv writes to stdout, so no file intermediates; still uses a
        TemporaryDirectory as cwd for safety against stray BLAST logs.
        """
        sample_name = input_file.stem
        print(f"Processing: {sample_name}")

        mlst_cmd = [
            "perl", str(self.mlst_bin),
            str(input_file.absolute()),
            "--scheme", scheme,
            "--csv",
            "--nopath",
        ]

        env = os.environ.copy()
        env["MLST_DATADIR"] = str(self.database_dir)

        raw_stdout = ""
        raw_stderr = ""

        try:
            with tempfile.TemporaryDirectory(prefix=f"mlst_{sample_name}_") as temp_dir:
                result = subprocess.run(
                    mlst_cmd,
                    capture_output=True,
                    text=True,
                    check=False,
                    env=env,
                    cwd=temp_dir,
                )
                raw_stdout = result.stdout or ""
                raw_stderr = result.stderr or ""
                returncode = result.returncode

            raw_combined = raw_stdout
            if raw_stderr:
                raw_combined += "\n--- STDERR ---\n" + raw_stderr

            if returncode == 0:
                mlst_results = self.parse_mlst_csv(raw_stdout, input_file.name)
                mlst_results["status"] = "Completed"
                mlst_results["raw_output"] = raw_combined
            else:
                mlst_results = self.get_fallback_results(input_file.name)
                mlst_results["status"] = f"Error: mlst exited with code {returncode}"
                mlst_results["raw_output"] = raw_combined
                if raw_stderr:
                    mlst_results["warnings"] = [raw_stderr.strip()[:300]]

            sample_output_dir = output_dir / sample_name
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            self.generate_output_files(mlst_results, sample_output_dir)

            st_display = mlst_results.get("st", "ND")
            if st_display in ("UNKNOWN", "ND", "-", ""):
                print(f"Completed: {input_file.name} -> ST not assigned")
            else:
                print(f"Completed: {input_file.name} -> ST{st_display}")
            return mlst_results

        except Exception as e:
            error_result = self.get_fallback_results(input_file.name)
            error_result["status"] = f"Error: {str(e)}"
            error_result["raw_output"] = (
                f"stdout:\n{raw_stdout}\n\nstderr:\n{raw_stderr}"
            )
            sample_output_dir = output_dir / sample_name
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            try:
                self.generate_output_files(error_result, sample_output_dir)
            except Exception:
                pass
            return error_result

    def parse_mlst_csv(self, stdout: str, sample_name: str) -> Dict:
        """Parse MLST CSV output (comma-separated)."""
        lines = stdout.strip().split('\n')
        if not lines:
            return self.get_empty_results(sample_name)

        result_line = None
        for line in reversed(lines):
            if line.strip() and ',' in line and not line.startswith('['):
                result_line = line.strip()
                break

        if not result_line:
            return self.get_empty_results(sample_name)

        parts = result_line.split(',')
        if len(parts) < 3:
            return self.get_empty_results(sample_name)

        filename = parts[0]
        scheme = parts[1]
        st = parts[2]

        alleles = {}
        allele_parts = []

        for i in range(3, len(parts)):
            allele_str = parts[i]
            if '(' in allele_str and ')' in allele_str:
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
            "confidence": "HIGH" if st and st != '-' and st != 'ND' else "LOW",
        }

    def get_empty_results(self, sample_name: str) -> Dict:
        """Return empty results structure."""
        return {
            "sample": sample_name,
            "st": "ND",
            "scheme": "ecoli_achtman_4",
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW",
        }

    def get_fallback_results(self, sample_name: str) -> Dict:
        """Fallback when MLST fails."""
        return {
            "sample": sample_name,
            "st": "UNKNOWN",
            "scheme": "ecoli_achtman_4",
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW",
            "status": "Error: MLST analysis failed",
            "raw_output": "",
            "warnings": ["MLST analysis failed"],
        }

    def generate_output_files(self, mlst_results: Dict, output_dir: Path):
        """Generate per-sample output: HTML report only (TSV is batch-level)."""
        self.generate_html_report(mlst_results, output_dir)

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
        .report-section h3 {
            color: #1e40af; margin-top: 20px; margin-bottom: 10px; font-size: 18px;
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
        .allele-grid {
            display: grid; grid-template-columns: repeat(auto-fill, minmax(150px, 1fr));
            gap: 15px; margin-top: 15px;
        }
        .allele-card {
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white; padding: 15px; border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1); text-align: center;
        }
        .allele-card .allele-gene { font-size: 12px; opacity: 0.9; margin-bottom: 5px; font-weight: bold; }
        .allele-card .allele-num { font-size: 20px; font-weight: bold; }
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
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white; padding: 15px; border-radius: 8px;
            text-align: center; box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }
        .stat-value { font-size: 24px; font-weight: bold; margin-bottom: 5px; }
        .stat-label { font-size: 12px; opacity: 0.9; }
        .st-badge {
            display: inline-block; background: #2563eb; color: white;
            padding: 5px 14px; border-radius: 15px; margin: 3px;
            font-size: 13px; font-weight: bold;
        }
        .confidence-badge {
            display: inline-block; padding: 4px 14px; border-radius: 12px;
            font-weight: bold; font-size: 14px; color: white;
        }
        .confidence-high { background: #16a34a; }
        .confidence-low { background: #dc2626; }
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
                    Powered by Seemann's mlst and PubMLST.<br>
                </p>
            </div>
        </div>
        """

    # ------------------------------------------------------------------
    # Per-sample HTML
    # ------------------------------------------------------------------

    def generate_html_report(self, mlst_results: Dict, output_dir: Path):
        """Write the single per-sample MLST HTML report."""
        sample_name = mlst_results["sample"]
        if sample_name.endswith((".fna", ".fasta", ".fa", ".fn", ".gz")):
            sample_name = Path(sample_name).stem
            if sample_name.endswith((".fna", ".fasta", ".fa", ".fn")):
                sample_name = Path(sample_name).stem

        random_quote = random.choice(self.science_quotes)

        st = mlst_results.get("st", "ND")
        status = mlst_results.get("status", "Completed")
        confidence = mlst_results.get("confidence", "LOW")

        is_success = (
            status == "Completed"
            and st not in ("ND", "UNKNOWN", "-", "")
        )

        if is_success:
            primary_label = "Sequence Type"
            primary_value = f"ST{st}"
            primary_bg = "linear-gradient(135deg, #10b981 0%, #059669 100%)"
        else:
            primary_label = "Status"
            primary_value = "Not Assigned"
            primary_bg = "linear-gradient(135deg, #6b7280 0%, #374151 100%)"

        confidence_class = (
            "confidence-high" if confidence == "HIGH" else "confidence-low"
        )

        allele_cards = ""
        for gene, allele in mlst_results.get("alleles", {}).items():
            allele_cards += f"""
                <div class="allele-card">
                    <div class="allele-gene">{gene}</div>
                    <div class="allele-num">{allele}</div>
                </div>
            """
        if not allele_cards:
            allele_cards = '<p style="color:#6b7280;">No alleles assigned.</p>'

        allele_profile = mlst_results.get("allele_profile", "") or "&mdash;"

        warnings_html = ""
        if mlst_results.get("warnings"):
            items = "".join(f"<li>{html_module.escape(w)}</li>" for w in mlst_results["warnings"])
            warnings_html = f'<div class="warning-box"><strong>Warnings:</strong><ul>{items}</ul></div>'

        failure_card = ""
        if not is_success:
            if status == "Completed":
                failure_card = """
        <div class="report-section" style="border-left: 4px solid #f59e0b;">
            <h2 style="color: #f59e0b;">Sequence Type Not Determined</h2>
            <p>mlst ran successfully but no sequence type could be assigned. This usually
            means one or more housekeeping genes failed to match the reference database.
            See raw output below for details.</p>
        </div>
                """
            else:
                failure_card = f"""
        <div class="report-section" style="border-left: 4px solid #dc2626;">
            <h2 style="color: #dc2626;">Analysis Error</h2>
            <p>{html_module.escape(status)}</p>
        </div>
                """

        raw_output_html = ""
        if mlst_results.get("raw_output"):
            escaped_raw = html_module.escape(mlst_results["raw_output"])
            raw_output_html = f"""
        <div class="report-section">
            <h2>Troubleshooting: Raw mlst Output</h2>
            <details>
                <summary style="cursor:pointer; color:#2563eb; font-weight:bold; font-size:14px;">
                    Show combined stdout + stderr (click to expand)
                </summary>
                <pre style="background:#f3f4f6; padding:15px; border-radius:6px;
                            max-height:400px; overflow-y:auto; font-size:12px;
                            color:#1f2937; white-space:pre-wrap;
                            margin-top:10px; border:1px solid #d1d5db;">{escaped_raw}</pre>
            </details>
        </div>
        """

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliTyper - {sample_name} MLST Report</title>
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
                    <div class="metric-label">MLST Scheme</div>
                    <div class="metric-value" style="font-size:16px;">{mlst_results.get('scheme', 'ecoli_achtman_4')}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Status</div>
                    <div class="metric-value" style="font-size:14px;">{html_module.escape(status)}</div>
                </div>
            </div>
        </div>

        {failure_card}

        <div class="report-section">
            <h2>Primary MLST Result</h2>
            <div class="primary-result" style="background: {primary_bg};">
                <div class="primary-label">{primary_label}</div>
                <div class="primary-value">{primary_value}</div>
                <span class="confidence-badge {confidence_class}">{confidence} Confidence</span>
            </div>
        </div>

        <div class="report-section">
            <h2>Allele Profile</h2>
            <p style="color:#6b7280;margin-bottom:10px;">
                Full allele profile in standard MLST notation (gene&rarr;allele number).
            </p>
            <div style="background:#f8f9fa; padding:15px; border-radius:8px;
                        margin:15px 0; border-left:4px solid #667eea;">
                <code style="font-size:14px; color:#1e3a8a; font-weight:bold; word-break:break-all;">
                    {allele_profile}
                </code>
            </div>
            <h3>Individual Alleles</h3>
            <div class="allele-grid">
                {allele_cards}
            </div>
        </div>

        <div class="report-section">
            <h2>Interpretation</h2>
            <div class="info-box">
                <strong>Method:</strong> Multi-Locus Sequence Typing (MLST) assigns a Sequence
                Type (ST) number by comparing the allelic profile of seven E. coli housekeeping
                genes (<em>adk</em>, <em>fumC</em>, <em>gyrB</em>, <em>icd</em>, <em>mdh</em>,
                <em>purA</em>, <em>recA</em>) against the Achtman scheme on PubMLST.
                STs sharing 4 or more identical alleles belong to the same clonal complex.
            </div>
            {warnings_html}
        </div>

        {raw_output_html}

        {self._get_footer_html("MLST Per-Sample Report")}
    </div>
    {self._get_quote_js()}
</body>
</html>"""

        html_file = output_dir / f"{sample_name}_mlst_report.html"
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

    # ------------------------------------------------------------------
    # Batch summaries
    # ------------------------------------------------------------------

    def create_mlst_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create comprehensive MLST summary files for all samples."""
        print("Creating MLST summary files...")
        self.create_mlst_tsv_summary(all_results, output_dir)
        self.create_mlst_html_summary(all_results, output_dir)
        self.create_mlst_json_summary(all_results, output_dir)
        print("MLST summary files created successfully.")

    def create_mlst_json_summary(self, all_results: Dict[str, Dict], output_dir: Path) -> str:
        """Generate structured JSON summary report."""
        total_samples = len(all_results)
        successful_samples = len(
            [r for r in all_results.values() if r["st"] not in ["ND", "UNKNOWN", "-"]]
        )

        st_distribution = {}
        for result in all_results.values():
            st = result["st"]
            if st not in ["ND", "UNKNOWN", "-"]:
                st_distribution[st] = st_distribution.get(st, 0) + 1

        st_distribution = dict(
            sorted(st_distribution.items(), key=lambda x: x[1], reverse=True)
        )
        unique_sts = list(st_distribution.keys())

        allele_frequencies = {}
        for result in all_results.values():
            for gene, allele in result["alleles"].items():
                if gene not in allele_frequencies:
                    allele_frequencies[gene] = {}
                allele_frequencies[gene][allele] = (
                    allele_frequencies[gene].get(allele, 0) + 1
                )

        for gene in allele_frequencies:
            allele_frequencies[gene] = dict(
                sorted(allele_frequencies[gene].items(), key=lambda x: x[1], reverse=True)
            )

        all_genes = set()
        for result in all_results.values():
            all_genes.update(result["alleles"].keys())

        detailed_results = []
        for sample_name, result in all_results.items():
            stem = Path(sample_name).stem
            detailed_results.append({
                "sample_id": sample_name,
                "sequence_type": f"ST{result['st']}",
                "st": result["st"],
                "allele_profile": result["allele_profile"],
                "confidence": result["confidence"],
                "alleles": result["alleles"],
                "file_name": sample_name,
                "status": result.get("status", "Completed"),
                "per_sample_report": f"{stem}/{stem}_mlst_report.html",
            })

        json_report = {
            "metadata": {
                **self.metadata,
                "analysis_type": "MLST",
                "mlst_scheme": "ecoli_achtman_4",
                "database_directory": str(self.database_dir),
                "total_files_analyzed": total_samples,
            },
            "analysis_summary": {
                "total_samples": total_samples,
                "successful_typing": successful_samples,
                "success_rate": round(successful_samples / total_samples * 100, 2)
                if total_samples > 0 else 0,
                "unique_sequence_types": len(unique_sts),
                "most_common_st": list(st_distribution.items())[0][0]
                if st_distribution else "None",
                "st_distribution": st_distribution,
            },
            "mlst_scheme_info": {
                "name": "ecoli_achtman_4",
                "description": "Escherichia coli MLST scheme based on 7 housekeeping genes",
                "reference": "Wirth T, Falush D, Lan R, et al. Sex and virulence in Escherichia coli: an evolutionary perspective. Mol Microbiol. 2006;60(5):1136-51.",
                "genes": sorted(all_genes),
            },
            "allele_statistics": {
                "gene_count": len(all_genes),
                "allele_frequencies": allele_frequencies,
                "most_common_alleles": {
                    gene: list(freq.items())[0][0] if freq else "None"
                    for gene, freq in allele_frequencies.items()
                },
            },
            "detailed_results": detailed_results,
            "quality_metrics": {
                "high_confidence_samples": len(
                    [r for r in all_results.values() if r["confidence"] == "HIGH"]
                ),
                "low_confidence_samples": len(
                    [r for r in all_results.values() if r["confidence"] == "LOW"]
                ),
                "failed_samples": total_samples - successful_samples,
            },
        }

        json_file = output_dir / "mlst_summary.json"
        with open(json_file, 'w') as f:
            json.dump(json_report, f, indent=4, default=str)

        print(f"JSON summary created: {json_file}")
        return str(json_file)

    def create_mlst_tsv_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create TSV summary file with all samples."""
        summary_file = output_dir / "mlst_summary.tsv"

        with open(summary_file, 'w') as f:
            f.write("Sample\tST\tAllele_Profile\tConfidence\t")

            all_genes = set()
            for result in all_results.values():
                all_genes.update(result["alleles"].keys())

            for gene in sorted(all_genes):
                f.write(f"{gene}\t")
            f.write("\n")

            for sample_name, result in all_results.items():
                f.write(
                    f"{sample_name}\t{result['st']}\t{result['allele_profile']}\t"
                    f"{result['confidence']}\t"
                )
                for gene in sorted(all_genes):
                    allele = result["alleles"].get(gene, "")
                    f.write(f"{allele}\t")
                f.write("\n")

        print(f"TSV summary created: {summary_file}")

    def create_mlst_html_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create batch HTML summary with links to per-sample reports."""
        summary_file = output_dir / "mlst_summary.html"

        all_genes = set()
        for result in all_results.values():
            all_genes.update(result["alleles"].keys())
        sorted_genes = sorted(all_genes)

        st_counts = Counter()
        for result in all_results.values():
            st_counts[result["st"]] += 1

        random_quote = random.choice(self.science_quotes)

        successful = sum(
            1 for r in all_results.values()
            if r["st"] not in ["ND", "UNKNOWN", "-"]
        )
        unique_sts = sum(1 for st in st_counts.keys() if st not in ["ND", "UNKNOWN", "-"])
        high_conf = sum(1 for r in all_results.values() if r["confidence"] == "HIGH")

        st_badges = "".join(
            f'<span class="st-badge">ST{st} ({count})</span>'
            for st, count in sorted(st_counts.items(), key=lambda x: x[1], reverse=True)
            if st not in ["ND", "UNKNOWN", "-"]
        )
        if not st_badges:
            st_badges = '<span style="color:#6b7280;">None detected</span>'

        header_cells = "".join(f"<th>{g}</th>" for g in sorted_genes)
        table_rows = ""
        for sample_name, result in all_results.items():
            stem = Path(sample_name).stem
            link = f"{stem}/{stem}_mlst_report.html"
            conf_class = "confidence-high" if result["confidence"] == "HIGH" else "confidence-low"
            st_display = (
                f"ST{result['st']}" if result["st"] not in ["ND", "UNKNOWN", "-"]
                else "&mdash;"
            )
            allele_cells = "".join(
                f"<td>{result['alleles'].get(g, '')}</td>" for g in sorted_genes
            )
            table_rows += f"""
                        <tr>
                            <td><strong><a href="{link}">{sample_name}</a></strong></td>
                            <td><strong style="color:#2563eb;">{st_display}</strong></td>
                            <td style="font-family:monospace;font-size:12px;">{result['allele_profile'] or '&mdash;'}</td>
                            <td><span class="confidence-badge {conf_class}">{result['confidence']}</span></td>
                            {allele_cells}
                        </tr>
            """

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliTyper - MLST Batch Summary</title>
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
            <h2>Batch MLST Summary</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{len(all_results)}</div>
                    <div class="stat-label">SAMPLES PROCESSED</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{successful}</div>
                    <div class="stat-label">TYPED SUCCESSFULLY</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{unique_sts}</div>
                    <div class="stat-label">UNIQUE SEQUENCE TYPES</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{high_conf}</div>
                    <div class="stat-label">HIGH CONFIDENCE</div>
                </div>
            </div>
            <p style="margin-top:15px;"><strong>Analysis Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>Scheme:</strong> ecoli_achtman_4</p>
        </div>

        <div class="report-section">
            <h2>Sequence Type Distribution</h2>
            <div style="margin: 15px 0;">{st_badges}</div>
            <div class="info-box">
                <strong>Method:</strong> MLST assigns sequence types (ST) by matching seven
                housekeeping genes against the Achtman scheme on PubMLST. Click any sample
                below to open its full per-sample report.
            </div>
        </div>

        <div class="report-section">
            <h2>Per-Sample Results</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>ST</th>
                            <th>Allele Profile</th>
                            <th>Confidence</th>
                            {header_cells}
                        </tr>
                    </thead>
                    <tbody>
                        {table_rows}
                    </tbody>
                </table>
            </div>
        </div>

        {self._get_footer_html("MLST Batch Summary")}
    </div>
    {self._get_quote_js()}
</body>
</html>"""

        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        print(f"HTML summary created: {summary_file}")

    def run_mlst_batch(self, input_path: str, output_dir: Path,
                       scheme: str = "ecoli_achtman_4") -> Dict[str, Dict]:
        """Run MLST analysis for one or more files. Handles single-file,
        directory, and glob input uniformly via find_fasta_files()."""
        print(self.ascii_art)
        print("EcoliTyper MLST Analysis")
        print("=" * 50)
        print(f"Input: {input_path}")
        print(f"Output: {output_dir}")
        print(f"Scheme: {scheme}")
        print(f"mlst binary: {self.mlst_bin}")
        print(f"Database:    {self.database_dir}")
        print("=" * 50)

        print("Searching for FASTA files...")
        fasta_files = self.find_fasta_files(input_path)

        if not fasta_files:
            print(f"No FASTA files found for input: {input_path}")
            return {}

        print(f"Found {len(fasta_files)} FASTA files")

        results = {}
        for fasta_file in fasta_files:
            result = self.run_mlst_single(fasta_file, output_dir, scheme)
            results[fasta_file.name] = result

        self.create_mlst_summary(results, output_dir)
        return results


def update_mlst_database(db_dir: Path, script_dir: Path,
                          scheme: str = "ecoli_achtman_4") -> bool:
    """
    Download and install the PubMLST E. coli MLST scheme using mlstdb.

    Runs mlstdb inside a temporary directory, locates the E. coli scheme under
    a set of known PubMLST folder names (ecoli_achtman_4, ecoli, etc.), copies
    it into db_dir/pubmlst/<scheme> using the fixed scheme name as the on-disk
    folder name, prunes any other schemes, and rebuilds the BLAST database.
    """
    if not shutil.which('mlstdb'):
        print("Error: 'mlstdb' not found. Install with: pip install mlstdb")
        print("       Then configure PubMLST access: mlstdb connect --db pubmlst --api-key")
        return False

    db_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix="ecolityper_mlst_update_"))
    print(f"Temporary download directory: {temp_dir}")

    try:
        print("Running mlstdb update (full download, may take 10-20 minutes)...")
        print("Progress will be shown below:\n")

        result = subprocess.run(
            ['mlstdb', 'update'],
            input='y\n',
            text=True,
            cwd=temp_dir,
            env=os.environ.copy(),
        )

        if result.returncode != 0:
            print(f"\nmlstdb update failed with exit code {result.returncode}")
            return False

        print("\nDownload completed.")

        temp_pubmlst = temp_dir / "pubmlst"
        if not temp_pubmlst.exists():
            print("pubmlst directory not found in downloaded data.")
            return False

        scheme_src = None
        if (temp_pubmlst / scheme).exists():
            scheme_src = temp_pubmlst / scheme
        else:
            for candidate in ("ecoli_achtman_4", "ecoli", "ecoli_achtman"):
                if (temp_pubmlst / candidate).exists():
                    scheme_src = temp_pubmlst / candidate
                    print(f"Source folder '{candidate}' will be installed as '{scheme}'.")
                    break
            if scheme_src is None:
                available = ", ".join(
                    sorted(p.name for p in temp_pubmlst.iterdir() if p.is_dir())
                ) or "(none)"
                print(f"E. coli scheme not found in downloaded data. Available: {available}")
                return False

        target_pubmlst = db_dir / "pubmlst"
        target_pubmlst.mkdir(parents=True, exist_ok=True)
        target_scheme = target_pubmlst / scheme
        if target_scheme.exists():
            shutil.rmtree(target_scheme)
        shutil.copytree(scheme_src, target_scheme)
        print(f"Installed scheme as {target_scheme}")

        blast_src = temp_dir / "blast"
        if blast_src.exists():
            target_blast = db_dir / "blast"
            if target_blast.exists():
                shutil.rmtree(target_blast)
            shutil.copytree(blast_src, target_blast)
            print(f"Copied BLAST database to {target_blast}")

        removed = 0
        for item in target_pubmlst.iterdir():
            if item.is_dir() and item.name != scheme:
                shutil.rmtree(item)
                removed += 1
                print(f"  Removed scheme: {item.name}")
        if removed:
            print(f"Removed {removed} other scheme(s).")
        else:
            print("No extra schemes to remove.")

        rebuild_script = script_dir.parent / "scripts" / "mlst-make_blast_db"
        if rebuild_script.exists():
            print("Rebuilding BLAST database...")
            rebuild_env = os.environ.copy()
            rebuild_env["MLST_DATADIR"] = str(db_dir)
            try:
                subprocess.run(
                    ['perl', str(rebuild_script)],
                    check=True,
                    cwd=str(script_dir.parent),
                    env=rebuild_env,
                )
                print("BLAST database rebuilt successfully.")
            except subprocess.CalledProcessError as e:
                print(f"BLAST rebuild failed ({e}); mlstdb output may suffice.")
        else:
            print(f"mlst-make_blast_db not found at {rebuild_script}")
            print("    The mlstdb download already produced a blast/ directory,")
            print("    so analysis can proceed; rebuild manually if you see errors.")

        print(f"MLST database updated and pruned to '{scheme}' in {db_dir}")
        return True

    except subprocess.CalledProcessError as e:
        print(f"Update failed: {e}")
        return False
    except Exception as e:
        print(f"Unexpected error: {e}")
        return False
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)
        print(f"Removed temporary directory: {temp_dir}")


def main():
    parser = argparse.ArgumentParser(
        description='EcoliTyper - MLST Analyzer for E. coli',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file analysis
  python ecolimlst_module.py -i genome.fasta -o results -db /path/to/db -sc /path/to/script

  # Batch analysis (quote globs — the shell must NOT expand them)
  python ecolimlst_module.py -i "*.fna" -o results -db /path/to/db -sc /path/to/script

  # Equivalent — pass a directory
  python ecolimlst_module.py -i . -o results -db /path/to/db -sc /path/to/script

  # Pull a fresh PubMLST E. coli MLST database (full download)
  python ecolimlst_module.py -db /path/to/db -sc /path/to/script --pull-db

  # Update the PubMLST E. coli MLST database (refresh in place)
  python ecolimlst_module.py -db /path/to/db -sc /path/to/script --update-db

Per-sample output (1 file each):
  <sample>/<sample>_mlst_report.html

Batch-level outputs (at the top of -o):
  mlst_summary.html
  mlst_summary.tsv
  mlst_summary.json

Author: Brown Beckley <brownbeckley94@gmail.com>
GitHub: https://github.com/bbeckley-hub
Affiliation: University of Ghana Medical School
        """,
    )

    parser.add_argument('-i', '--input',
                        help='Input FASTA file, directory, or glob pattern. '
                             'Quote globs (e.g. -i "*.fna") or pass the directory (e.g. -i .)')
    parser.add_argument('-o', '--output-dir',
                        help='Output directory')
    parser.add_argument('-db', '--database-dir', required=True,
                        help='Database directory')
    parser.add_argument('-sc', '--script-dir', required=True,
                        help='Script directory (contains mlst binary)')
    parser.add_argument('-s', '--scheme', default='ecoli_achtman_4',
                        help='MLST scheme (default: ecoli_achtman_4)')
    parser.add_argument('--batch', action='store_true',
                        help='(Deprecated — batch mode is auto-detected from input glob/directory)')
    parser.add_argument('--pull-db', action='store_true',
                        help='Download a fresh PubMLST E. coli MLST database and exit')
    parser.add_argument('--update-db', action='store_true',
                        help='Update the PubMLST E. coli MLST database and exit')

    args = parser.parse_args()

    db_dir = Path(args.database_dir)
    script_dir = Path(args.script_dir)

    if args.pull_db or args.update_db:
        success = update_mlst_database(db_dir, script_dir, args.scheme)
        sys.exit(0 if success else 1)

    if not args.input or not args.output_dir:
        parser.error("When not using --pull-db or --update-db, -i and -o are required.")

    try:
        analyzer = EcoliTyper(database_dir=db_dir, script_dir=script_dir)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)


    results = analyzer.run_mlst_batch(args.input, output_dir, args.scheme)

    if not results:
        print(f"No FASTA files found for input: {args.input}")
        print('Hint: quote globs (e.g. -i "*.fna") or pass the directory (e.g. -i .)')
        sys.exit(1)

    print(f"MLST completed. Processed {len(results)} sample(s)")


if __name__ == "__main__":
    main()