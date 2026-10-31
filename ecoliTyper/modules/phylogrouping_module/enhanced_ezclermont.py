#!/usr/bin/env python3
"""
EcoliTyper - Enhanced ezClermont Wrapper
E. coli phylogrouping using ezClermont with beautiful reporting.

Per-sample output (1 file each):
    <sample>/<sample>_phylogroup_report.html

Batch-level outputs (at the top of -o):
    phylogrouping_results.html
    phylogrouping_results.tsv
    phylogrouping_results.json

All intermediate files (raw ezClermont output, stderr logs) live inside a
TemporaryDirectory that self-destructs on exit — no intermediates left on disk.

Author: Brown Beckley
Email: <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
Date: 2025-12-16 / 2026-07-22
Send a quick mail for any issues or further explanations.
"""

import os
import sys
import json
import html as html_module
import argparse
import subprocess
import glob
import re
import tempfile
import random
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any
import pandas as pd
import multiprocessing as mp


class EnhancedEzClermont:
    """Enhanced ezClermont wrapper with per-sample HTML reporting."""

    def __init__(self, threads: int = 4):
        self.threads = threads
        self.results = []
        self.metadata = {
            "tool_name": "EcoliTyper Phylogrouping",
            "version": "1.4.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }

        self.ezclermont_path = self._detect_ezclermont()
        if not self.ezclermont_path:
            raise RuntimeError(
                "ezClermont not found. Please install ezClermont or ensure it's in your PATH"
            )

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

    def _detect_ezclermont(self) -> str:
        """Auto-detect ezClermont installation."""
        try:
            result = subprocess.run(["which", "ezclermont"], capture_output=True, text=True)
            if result.returncode == 0:
                return result.stdout.strip()
        except Exception:
            pass
        return ""

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

    def run_ezclermont_analysis(self, fasta_file: Path, output_base: Path) -> Dict[str, Any]:
        """Run ezClermont on one genome. Emits exactly 1 file per sample:
             - <sample>/<sample>_phylogroup_report.html
        All intermediates live inside a TemporaryDirectory that self-destructs.
        """
        sample_name = fasta_file.stem
        try:
            print(f"Analyzing {sample_name}...")

            with tempfile.TemporaryDirectory(prefix=f"ezclermont_{sample_name}_") as temp_dir:
                temp_path = Path(temp_dir)

                cmd = [
                    self.ezclermont_path,
                    str(fasta_file.absolute()),
                ]

                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=300,
                    cwd=str(temp_path),
                )

                parsed_result = self._parse_sample_results(
                    sample_name,
                    str(fasta_file),
                    str(output_base / sample_name),
                    result.stdout or "",
                    result.stderr or "",
                )

                if result.returncode != 0:
                    if parsed_result["clermont_type"] != "Unknown":
                        parsed_result["status"] = "Completed"
                        parsed_result["warnings"].append(
                            f"ezClermont returned exit code {result.returncode} "
                            f"but analysis completed successfully"
                        )
                    else:
                        parsed_result["status"] = (
                            f"Error: ezClermont failed with exit code {result.returncode}"
                        )
                        if result.stderr:
                            parsed_result["warnings"].append(
                                f"stderr: {result.stderr.strip()[:300]}"
                            )
                else:
                    parsed_result["status"] = "Completed"

            sample_output_dir = output_base / sample_name
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            self._generate_per_sample_html(parsed_result, sample_output_dir)
            parsed_result["output_directory"] = str(sample_output_dir)

            return parsed_result

        except subprocess.TimeoutExpired:
            error = self._create_error_result(
                sample_name, str(fasta_file), "Analysis timed out after 5 minutes"
            )
            sample_output_dir = output_base / sample_name
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            self._generate_per_sample_html(error, sample_output_dir)
            error["output_directory"] = str(sample_output_dir)
            return error
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

    def _parse_sample_results(self, sample_name: str, fasta_path: str,
                              output_dir: str, stdout_text: str,
                              stderr_text: str = "") -> Dict[str, Any]:
        """Parse results from ezClermont combined stdout + stderr.

        ezClermont writes its final `GCA_xxx\\tB2` result line to stdout, but
        the marker / progress information (TspE4, arpA, chu, yjaA, trpBA, virA)
        goes to stderr on this build. Both streams are therefore combined
        before parsing so no marker evidence is lost.
        """
        try:
            stdout_text = stdout_text or ""
            stderr_text = stderr_text or ""

            combined_text = stdout_text
            if stderr_text:
                combined_text += "\n" + stderr_text

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
                "output_directory": output_dir,
                "warnings": [],
                "raw_output": combined_text,
            }

            text_to_parse = combined_text.strip()

            if not text_to_parse:
                result_data["warnings"].append(
                    "ezClermont produced no output on stdout or stderr"
                )
                self._validate_results(result_data)
                return result_data

            marker_patterns = [
                ("trpba_control", r"trpba[_\s]?control\b"),
                ("vira", r"\bvir\s*a\b"),
                ("tspe4", r"\btsp\s*e4(?:\.c2)?\b"),
                ("arpa", r"\barp\s*a\b"),
                ("chu", r"\bchu(?:a)?\b"),
                ("yjaa", r"\byja\s*a\b"),
            ]

            for line in text_to_parse.split('\n'):
                line = line.strip()
                if not line:
                    continue
                lower = line.lower()

                for marker_key, prefix_pattern in marker_patterns:
                    if result_data[marker_key] == "Unknown" and \
                       re.search(prefix_pattern, line, re.IGNORECASE):
                        extracted = self._extract_gene_result(line)
                        if extracted != "Unknown":
                            result_data[marker_key] = extracted

                if lower.startswith("clermont type"):
                    parts = line.split(':', 1)
                    if len(parts) == 2:
                        value = parts[1].strip()
                        if value in ['A', 'B1', 'B2', 'C', 'D', 'E', 'F', 'G']:
                            result_data["clermont_type"] = value

                elif '\t' in line:
                    parts = line.split('\t')
                    if len(parts) == 2:
                        value = parts[1].strip()
                        if value in ['A', 'B1', 'B2', 'C', 'D', 'E', 'F', 'G']:
                            result_data["clermont_type"] = value

                if line.startswith("('") and line.endswith("')"):
                    tuple_match = re.search(r"\('([^']+)',\s*'([^']*)'\)", line)
                    if tuple_match:
                        tuple_type = tuple_match.group(1).strip()
                        if tuple_type in ['A', 'B1', 'B2', 'C', 'D', 'E', 'F', 'G']:
                            result_data["clermont_type"] = tuple_type
                        markers_text = tuple_match.group(2).replace('\\n', '\n')
                        for marker_key, regex_pat in [
                            ('tspe4', r'TspE4\s*:\s*([+-])'),
                            ('arpa', r'arpA\s*:\s*([+-])'),
                            ('chu', r'chu\s*:\s*([+-])'),
                            ('yjaa', r'yjaA\s*:\s*([+-])'),
                        ]:
                            if result_data[marker_key] == "Unknown":
                                m = re.search(regex_pat, markers_text)
                                if m:
                                    result_data[marker_key] = m.group(1)

            if result_data["clermont_type"] != "Unknown":
                known_markers = sum(
                    1 for m in ['tspe4', 'arpa', 'chu', 'yjaa']
                    if result_data[m] != "Unknown"
                )
                if known_markers == 0:
                    stale_type = result_data["clermont_type"]
                    result_data["clermont_type"] = "Unknown"
                    result_data["warnings"].append(
                        f"Assigned type '{stale_type}' had no supporting marker evidence — "
                        f"treated as Unknown. See raw output below for details."
                    )

            self._validate_results(result_data)
            return result_data

        except Exception as e:
            return self._create_error_result(
                sample_name, fasta_path, f"Error parsing results: {str(e)}"
            )

    def _extract_gene_result(self, line: str) -> str:
        """Extract + or - from a marker line.

        The trailing sign may be followed by whitespace or end-of-string —
        NOT a word character (which is why the previous `\\b` anchor failed
        at end-of-line for `+`/`-`).

        Handles all of:
          • "TspE4: +"
          • "TspE4:+"
          • "TspE4: +   (extra note)"
          • "TspE4 = +"
        """
        m = re.search(r'[:=]\s*([+-])(?:\s|$)', line)
        if m:
            return m.group(1)
        return "Unknown"

    def _validate_results(self, result_data: Dict[str, Any]):
        """Validate parsed results and add warnings if needed."""
        warnings = list(result_data.get("warnings", []))

        if result_data["clermont_type"] == "Unknown":
            warnings.append("Could not determine Clermont phylogroup")

        unknown_markers = []
        for marker in ['tspe4', 'arpa', 'chu', 'yjaa', 'trpba_control', 'vira']:
            if result_data[marker] == "Unknown":
                unknown_markers.append(marker)

        if unknown_markers:
            warnings.append(f"Could not parse markers: {', '.join(unknown_markers)}")

        result_data["warnings"] = warnings

    def _create_error_result(self, sample_name: str, fasta_path: str,
                             error_msg: str) -> Dict[str, Any]:
        """Create error result structure."""
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
            "raw_output": "",
        }

    def process_batch(self, input_path: str, main_output_dir: Path) -> List[Dict[str, Any]]:
        """Process all FASTA files in batch."""
        fasta_files = self.find_fasta_files(input_path)

        print(f"Processing {len(fasta_files)} samples using {self.threads} threads...")
        print(f"Using ezClermont at: {self.ezclermont_path}")

        args = [(fasta_file, main_output_dir) for fasta_file in fasta_files]

        with mp.Pool(processes=self.threads) as pool:
            results = pool.starmap(self.run_ezclermont_analysis, args)

        self.results = results
        return results

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
            font-weight: bold; font-size: 15px; color: white;
        }
        .marker-present { background: #16a34a; }
        .marker-absent { background: #dc2626; }
        .marker-unknown { background: #6b7280; }
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
                    Powered by ezClermont. Clermont phylotyping: Clermont O, Christenson JK, Denamur E, Gordon DM.
                    Environ Microbiol Rep. 2013;5(1):58-65.
                </p>
            </div>
        </div>
        """

    def _generate_per_sample_html(self, result: Dict[str, Any], output_dir: Path):
        """Write the single per-sample phylogroup HTML report."""
        sample_name = result['sample_id']
        random_quote = random.choice(self.science_quotes)
        is_success = result['status'] == 'Completed' and result['clermont_type'] != 'Unknown'
        clermont_type = result['clermont_type']

        type_backgrounds = {
            'A': 'linear-gradient(135deg, #16a34a 0%, #15803d 100%)',
            'B1': 'linear-gradient(135deg, #2563eb 0%, #1e40af 100%)',
            'B2': 'linear-gradient(135deg, #dc2626 0%, #991b1b 100%)',
            'C': 'linear-gradient(135deg, #f59e0b 0%, #b45309 100%)',
            'D': 'linear-gradient(135deg, #7c3aed 0%, #5b21b6 100%)',
            'E': 'linear-gradient(135deg, #ec4899 0%, #9d174d 100%)',
            'F': 'linear-gradient(135deg, #06b6d4 0%, #0e7490 100%)',
            'G': 'linear-gradient(135deg, #84cc16 0%, #4d7c0f 100%)',
            'Unknown': 'linear-gradient(135deg, #6b7280 0%, #374151 100%)',
        }
        primary_bg = type_backgrounds.get(clermont_type, type_backgrounds['Unknown'])

        if is_success:
            primary_label = 'Clermont Phylogroup'
            primary_value = clermont_type
        else:
            primary_label = 'Status'
            primary_value = 'Not Assigned'

        markers = [
            ('TspE4.C2', result['tspe4']),
            ('arpA', result['arpa']),
            ('chuA', result['chu']),
            ('yjaA', result['yjaa']),
            ('trpBA control', result['trpba_control']),
            ('virA', result['vira']),
        ]

        marker_rows = ''
        for name, value in markers:
            if value == '+':
                badge = '<span class="marker-badge marker-present">Present (+)</span>'
            elif value == '-':
                badge = '<span class="marker-badge marker-absent">Absent (&minus;)</span>'
            else:
                badge = '<span class="marker-badge marker-unknown">Unknown (?)</span>'
            marker_rows += f'<tr><td><strong>{name}</strong></td><td>{badge}</td></tr>'

        warnings_html = ''
        if result.get('warnings'):
            items = ''.join(f'<li>{w}</li>' for w in result['warnings'])
            warnings_html = f'<div class="warning-box"><strong>Warnings:</strong><ul>{items}</ul></div>'

        failure_card = ''
        if not is_success:
            if result['status'] == 'Completed':
                failure_card = '''
        <div class="report-section" style="border-left: 4px solid #f59e0b;">
            <h2 style="color: #f59e0b;">Phylogroup Not Determined</h2>
            <p>ezClermont ran successfully but its output could not be parsed into a valid
            phylogroup assignment. See the warnings and raw output below for troubleshooting.</p>
        </div>
                '''
            else:
                failure_card = f'''
        <div class="report-section" style="border-left: 4px solid #dc2626;">
            <h2 style="color: #dc2626;">Analysis Error</h2>
            <p>{result['status']}</p>
        </div>
                '''

        raw_output_html = ''
        if result.get('raw_output'):
            escaped_raw = html_module.escape(result['raw_output'])
            raw_output_html = f'''
        <div class="report-section">
            <h2>Troubleshooting: Raw ezClermont Output</h2>
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
        '''

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliTyper - {sample_name} Phylogroup Report</title>
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
            <h2>Primary Phylogroup Result</h2>
            <div class="primary-result" style="background: {primary_bg};">
                <div class="primary-label">{primary_label}</div>
                <div class="primary-value">{primary_value}</div>
            </div>
        </div>

        <div class="report-section">
            <h2>Clermont Marker Genotypes</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr><th>Marker</th><th>Result</th></tr>
                    </thead>
                    <tbody>
                        {marker_rows}
                    </tbody>
                </table>
            </div>
            <div class="info-box">
                <strong>Method:</strong> Clermont phylotyping uses presence/absence of the
                <em>chuA</em>, <em>yjaA</em>, <em>TspE4.C2</em>, and <em>arpA</em> markers
                to assign E. coli isolates to phylogroups A, B1, B2, C, D, E, F, or G.
                <em>trpBA</em> and <em>virA</em> serve as control markers for the assay.
            </div>
            {warnings_html}
        </div>

        {raw_output_html}

        {self._get_footer_html("Phylogroup Per-Sample Report")}
    </div>
    {self._get_quote_js()}
</body>
</html>"""

        html_file = output_dir / f"{sample_name}_phylogroup_report.html"
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

    def generate_json_report(self, output_dir: Path) -> str:
        """Generate structured JSON summary report."""
        successful = len([r for r in self.results if r['status'] == 'Completed'])
        failed = len([r for r in self.results if r['status'].startswith('Error')])

        clermont_distribution = {}
        for result in self.results:
            if result['status'] == 'Completed' and result['clermont_type'] != 'Unknown':
                clermont_type = result['clermont_type']
                clermont_distribution[clermont_type] = clermont_distribution.get(clermont_type, 0) + 1

        clermont_distribution = dict(sorted(
            clermont_distribution.items(),
            key=lambda x: x[1],
            reverse=True,
        ))

        marker_results = {
            "tspe4": {"+": 0, "-": 0, "Unknown": 0},
            "arpa": {"+": 0, "-": 0, "Unknown": 0},
            "chu": {"+": 0, "-": 0, "Unknown": 0},
            "yjaa": {"+": 0, "-": 0, "Unknown": 0},
        }

        for result in self.results:
            if result['status'] == 'Completed':
                for marker in marker_results.keys():
                    marker_value = result.get(marker, "Unknown")
                    marker_results[marker][marker_value] = marker_results[marker].get(marker_value, 0) + 1

        simplified_results = []
        for result in self.results:
            simplified_results.append({
                "sample_id": result["sample_id"],
                "clermont_type": result["clermont_type"],
                "status": result["status"],
                "file_path": result["file_path"],
                "per_sample_report": (
                    f"{result['sample_id']}/{result['sample_id']}_phylogroup_report.html"
                ),
                "markers": {
                    "tspe4": result["tspe4"],
                    "arpa": result["arpa"],
                    "chu": result["chu"],
                    "yjaa": result["yjaa"],
                },
            })

        json_report = {
            "metadata": {
                **self.metadata,
                "ezclermont_path": self.ezclermont_path,
                "threads_used": self.threads,
            },
            "analysis_summary": {
                "total_samples": len(self.results),
                "successful_analyses": successful,
                "failed_analyses": failed,
                "success_rate": round(successful / len(self.results) * 100, 1) if len(self.results) > 0 else 0,
                "unique_clermont_types": len(clermont_distribution),
                "clermont_type_distribution": clermont_distribution,
            },
            "marker_summary": marker_results,
            "clermont_type_descriptions": {
                "A": "Phylogroup A: Typically commensal strains",
                "B1": "Phylogroup B1: Often associated with environmental strains",
                "B2": "Phylogroup B2: Frequently contains extra-intestinal pathogenic E. coli (ExPEC)",
                "C": "Phylogroup C: Intermediate between B2 and D",
                "D": "Phylogroup D: Contains many pathogenic strains including uropathogenic E. coli",
                "E": "Phylogroup E: Recently identified group, often associated with specific pathotypes",
                "F": "Phylogroup F: Associated with certain pathogenic strains",
                "G": "Phylogroup G: Newly identified group with distinct characteristics",
                "Unknown": "Could not be classified into known Clermont phylogroups",
            },
            "detailed_results": simplified_results,
            "analysis_notes": {
                "method": "Clermont phylotyping based on chuA, yjaA, TspE4.C2, and arpA genes",
                "algorithm": "Original Clermont algorithm (2000, 2013)",
                "reference": "Clermont O, Christenson JK, Denamur E, Gordon DM. The Clermont Escherichia coli phylo-typing method revisited: improvement of specificity and detection of new phylo-groups. Environ Microbiol Rep. 2013;5(1):58-65.",
                "version": "ezClermont V0.7.0",
            },
        }

        json_file = output_dir / "phylogrouping_results.json"
        with open(json_file, 'w') as f:
            json.dump(json_report, f, indent=4, default=str)

        return str(json_file)

    def generate_html_report(self, output_dir: Path) -> str:
        """Generate comprehensive batch HTML report with links to per-sample reports."""
        random_quote = random.choice(self.science_quotes)

        successful = len([r for r in self.results if r['status'] == 'Completed'])
        failed = len([r for r in self.results if r['status'].startswith('Error')])

        clermont_types = {}
        for result in self.results:
            if result['status'] == 'Completed' and result['clermont_type'] != 'Unknown':
                clermont_types[result['clermont_type']] = clermont_types.get(result['clermont_type'], 0) + 1

        type_badges_html = ''.join(
            f'<span class="marker-badge marker-present" style="margin:4px;">{t} ({c})</span>'
            for t, c in sorted(clermont_types.items())
        )
        if not type_badges_html:
            type_badges_html = '<span style="color:#6b7280;">None detected</span>'

        rows = ''
        for result in self.results:
            status_class = 'marker-present' if result['status'] == 'Completed' else 'marker-absent'
            link = f"{result['sample_id']}/{result['sample_id']}_phylogroup_report.html"
            rows += f"""
                        <tr>
                            <td><strong><a href="{link}">{result['sample_id']}</a></strong></td>
                            <td><strong style="color: #2563eb;">{result['clermont_type']}</strong></td>
                            <td><span class="marker-badge {status_class}">{result['status']}</span></td>
                        </tr>
            """

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliTyper - Phylogrouping Batch Summary</title>
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
            <h2>Batch Phylogrouping Summary</h2>
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
                    <div class="stat-value">{len(clermont_types)}</div>
                    <div class="stat-label">UNIQUE PHYLOGROUPS</div>
                </div>
            </div>
            <p style="margin-top:15px;"><strong>Analysis Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>Threads Used:</strong> {self.threads}</p>
            <p><strong>ezClermont Path:</strong> {self.ezclermont_path}</p>
        </div>

        <div class="report-section">
            <h2>Clermont Phylogroup Distribution</h2>
            <div style="margin: 15px 0;">{type_badges_html}</div>
            <div class="info-box">
                <strong>Method:</strong> Clermont phylotyping is based on the original Clermont
                algorithm using <em>chuA</em>, <em>yjaA</em>, <em>TspE4.C2</em>, and <em>arpA</em>
                markers. Click any sample below to open its full per-sample report.
            </div>
        </div>

        <div class="report-section">
            <h2>Per-Sample Results</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr><th>Sample ID</th><th>Clermont Type</th><th>Status</th></tr>
                    </thead>
                    <tbody>
                        {rows}
                    </tbody>
                </table>
            </div>
        </div>

        {self._get_footer_html("Phylogrouping Batch Summary")}
    </div>
    {self._get_quote_js()}
</body>
</html>"""

        html_file = output_dir / "phylogrouping_results.html"
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        return str(html_file)

    def generate_tsv_report(self, output_dir: Path) -> str:
        """Generate simplified TSV report (no per-marker columns)."""
        data = []
        for result in self.results:
            data.append({
                'Sample_ID': result['sample_id'],
                'Clermont_Type': result['clermont_type'],
                'Status': result['status'],
                'File_Path': result['file_path'],
            })

        df = pd.DataFrame(data)
        tsv_file = output_dir / "phylogrouping_results.tsv"
        df.to_csv(tsv_file, sep='\t', index=False)
        return str(tsv_file)


def main():
    parser = argparse.ArgumentParser(
        description="EcoliTyper Enhanced ezClermont - Phylogrouping with batch processing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Per-sample output (1 file each):
  <sample>/<sample>_phylogroup_report.html

Batch-level outputs (at the top of -o):
  phylogrouping_results.html
  phylogrouping_results.tsv
  phylogrouping_results.json
        """,
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input FASTA file, directory, or glob pattern')
    parser.add_argument('-o', '--output', required=True,
                        help='Main output directory')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='Number of threads to use')

    args = parser.parse_args()

    main_output_dir = Path(args.output)
    main_output_dir.mkdir(parents=True, exist_ok=True)

    try:
        finder = EnhancedEzClermont(args.threads)

        print(finder.ascii_art)
        print("EcoliTyper Enhanced ezClermont")
        print("=" * 50)
        print(f"Input: {args.input}")
        print(f"Output: {main_output_dir}")
        print(f"Threads: {args.threads}")
        print(f"Auto-detected ezClermont: {finder.ezclermont_path}")
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
        print("  <sample>/<sample>_phylogroup_report.html")
        print("\nBatch-level outputs:")
        print(f"  {html_file}")
        print(f"  {tsv_file}")
        print(f"  {json_file}")

        successful = len([r for r in results if r['status'] == 'Completed'])
        if len(results) > 0:
            print(f"\nSuccess rate: {successful}/{len(results)} ({successful/len(results)*100:.1f}%)")

        clermont_types = list(set(
            r['clermont_type'] for r in results
            if r['status'] == 'Completed' and r['clermont_type'] != 'Unknown'
        ))
        print(f"Clermont types found: {', '.join(clermont_types) if clermont_types else 'None'}")

        quote = random.choice(finder.science_quotes)
        print(f"\n{quote['text']} - {quote['author']}")

    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()