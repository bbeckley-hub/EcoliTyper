#!/usr/bin/env python3
"""
EcoliTyper Main Orchestrator - Complete E. coli Typing Pipeline (Temp‑dir version)/ HPC-friendly/ Docker
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
MIT License
2026-07-23
Version: 1.3.0
"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
import signal
import threading
import tempfile
import time
import traceback
import logging
from pathlib import Path
from typing import Dict, List, Optional

class Colors:
    RESET = "\033[0m"
    BOLD = "\033[1m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    WHITE = "\033[37m"
    BRIGHT_RED = "\033[91m"
    BRIGHT_GREEN = "\033[92m"
    BRIGHT_YELLOW = "\033[93m"
    BRIGHT_BLUE = "\033[94m"
    BRIGHT_MAGENTA = "\033[95m"
    BRIGHT_CYAN = "\033[96m"
    BRIGHT_WHITE = "\033[97m"


class ColoredHelpFormatter(argparse.RawTextHelpFormatter, argparse.HelpFormatter):
    def _format_action(self, action):
        action_str = super()._format_action(action)
        if not action_str:
            return action_str
        lines = action_str.split('\n')
        colored_lines = []
        for line in lines:
            if line.strip():
                if line.lstrip().startswith('-'):
                    parts = line.split('  ', 1)
                    if len(parts) == 2:
                        options = parts[0].strip()
                        help_text = parts[1]
                        colored_line = f"  {Colors.BRIGHT_CYAN}{options}{Colors.RESET}  {help_text}"
                    else:
                        colored_line = f"  {Colors.BRIGHT_CYAN}{line.strip()}{Colors.RESET}"
                else:
                    colored_line = f"  {Colors.YELLOW}{line}{Colors.RESET}"
                colored_lines.append(colored_line)
            else:
                colored_lines.append(line)
        return '\n'.join(colored_lines)

    def start_section(self, heading):
        heading = f"{Colors.BOLD}{Colors.BRIGHT_GREEN}{heading}{Colors.RESET}"
        super().start_section(heading)


try:
    from .core.banner import EcoliTyperBanner
except (ImportError, SystemError):
    sys.path.insert(0, str(Path(__file__).parent))
    from core.banner import EcoliTyperBanner


class EcoliTyperOrchestrator:
    def __init__(self):
        self.banner = EcoliTyperBanner()
        self.base_dir = Path(__file__).parent
        self.fasta_files = []
        self.interrupted = False
        self.output_lock = threading.Lock()
        self.start_time = None
        self.total_duration = None
        self.keep_temp = False
        self.user_output_dir = None
        self.logger = None

        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)

    def setup_logging(self, output_dir: Path):
        log_file = output_dir / "ecolityper_run.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[logging.FileHandler(log_file, mode='w')]
        )
        self.logger = logging.getLogger("EcoliTyper")

    def _log_info(self, msg: str):
        if self.logger:
            self.logger.info(msg)

    def _log_error(self, msg: str):
        if self.logger:
            self.logger.error(msg)

    def _log_warning(self, msg: str):
        if self.logger:
            self.logger.warning(msg)

    def _format_duration(self, seconds: float) -> str:
        if seconds < 60:
            return f"{seconds:.1f} seconds"
        elif seconds < 3600:
            minutes = int(seconds // 60)
            secs = seconds % 60
            return f"{minutes} minutes {secs:.0f} seconds"
        else:
            hours = int(seconds // 3600)
            minutes = int((seconds % 3600) // 60)
            secs = seconds % 60
            return f"{hours} hours {minutes} minutes {secs:.0f} seconds"

    def _signal_handler(self, signum, frame):
        self.interrupted = True
        if self.start_time:
            self.total_duration = time.time() - self.start_time
            self._log_error(f"Total analysis ran for: {self._format_duration(self.total_duration)}")
        self._log_error(f"Analysis interrupted by user (signal {signum})")
        self.banner.display_error("Analysis interrupted by user – cleaning up...")
        sys.exit(1)

    def update_amr_database(self, force: bool = False) -> bool:
        amr_module_path = self.base_dir / "modules" / "amrfinder_module"
        amr_script = amr_module_path / "ecoli_amrfinder.py"
        if not amr_script.exists():
            self.banner.display_error(f"AMR script not found at: {amr_script}")
            return False
        self.banner.display_info("Updating AMRfinderPlus database...")
        cmd = [sys.executable, str(amr_script), "--force-update" if force else "--update-db"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
        if result.returncode == 0:
            self.banner.display_success("AMR database updated successfully.")
            version_cmd = [sys.executable, str(amr_script), "--db-version"]
            version_result = subprocess.run(version_cmd, capture_output=True, text=True, cwd=amr_module_path)
            if version_result.returncode == 0:
                self.banner.display_info(f"New database version: {version_result.stdout.strip()}")
            return True
        else:
            self.banner.display_error("AMR database update failed.")
            if result.stderr:
                self._log_error(result.stderr)
            return False

    def ensure_amr_database(self) -> bool:
        amr_module_path = self.base_dir / "modules" / "amrfinder_module"
        amr_script = amr_module_path / "ecoli_amrfinder.py"
        if not amr_script.exists():
            self.banner.display_error("AMR script not found, cannot check database.")
            return False
        cmd = [sys.executable, str(amr_script), "--db-version"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
        if result.returncode == 0 and "Unknown" not in result.stdout and "No database" not in result.stdout:
            self.banner.display_success(f"AMR database already present: {result.stdout.strip()}")
            return True
        else:
            self.banner.display_warning("AMR database not found or outdated. Attempting automatic update...")
            return self.update_amr_database(force=False)

    def find_fasta_files(self, input_path: str) -> List[Path]:
        self.banner.display_info(f"Searching for files with pattern: {input_path}")
        if '*' in input_path or '?' in input_path:
            matched_files = glob.glob(input_path)
            fasta_files = [Path(f) for f in matched_files if Path(f).is_file() and
                           f.lower().endswith(('.fna', '.fasta', '.fa', '.fsa')) and
                           not Path(f).name.startswith('.')]
            self.banner.display_success(f"Found {len(fasta_files)} FASTA files")
            return sorted(fasta_files)
        input_path_obj = Path(input_path)
        if input_path_obj.is_file() and input_path_obj.suffix.lower() in ['.fna', '.fasta', '.fa', '.fsa']:
            self.banner.display_success(f"Found single FASTA file: {input_path_obj.name}")
            return [input_path_obj]
        if input_path_obj.is_dir():
            patterns = [f"{input_path}/*.fna", f"{input_path}/*.fasta", f"{input_path}/*.fa", f"{input_path}/*.fsa"]
            fasta_files = []
            for pattern in patterns:
                for file_path in glob.glob(pattern):
                    path = Path(file_path)
                    if path.is_file() and not path.name.startswith('.'):
                        fasta_files.append(path)
            fasta_files = sorted(list(set(fasta_files)))
            if fasta_files:
                self.banner.display_success(f"Found {len(fasta_files)} FASTA files in directory")
            else:
                self.banner.display_warning(f"No FASTA files found in directory: {input_path}")
            return fasta_files
        self.banner.display_error(f"Input path not found: {input_path}")
        return []

    def get_file_pattern(self, fasta_files: List[Path]) -> str:
        if not fasta_files:
            return "*.fna"
        extensions = set(f.suffix.lower() for f in fasta_files)
        if len(extensions) == 1:
            ext = list(extensions)[0]
            return f"*{ext}"
        return "*"

    def run_module_in_temp(self, module_name: str, fasta_files: List[Path],
                           cmd: List[str], result_subdir: str, target_dir: str,
                           extra_result_files: List[str] = None,
                           env: Dict[str, str] = None) -> bool:
        module_orig = self.base_dir / "modules" / module_name
        if not module_orig.exists():
            self.banner.display_error(f"Module directory not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix=f"ecolityper_{module_name}_"))
        self._log_info(f"Temporary directory for {module_name}: {temp_dir}")

        try:
            shutil.copytree(module_orig, temp_dir / module_name, dirs_exist_ok=True)
            work_dir = temp_dir / module_name
            for f in fasta_files:
                shutil.copy2(f, work_dir / f.name)

            self.banner.display_info(f"Running {module_name} analysis...")
            self._log_info(f"Command: {' '.join(cmd)}")

            result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True, env=env)
            if result.stdout:
                self._log_info(result.stdout[:500])
            if result.stderr:
                self._log_warning(result.stderr[:500])

            if result.returncode != 0:
                self.banner.display_error(f"{module_name} failed with return code {result.returncode}")
                return False

            src_dir = work_dir / result_subdir
            if src_dir.exists():
                dst_dir = self.user_output_dir / target_dir
                if dst_dir.exists():
                    shutil.rmtree(dst_dir)
                shutil.copytree(src_dir, dst_dir)
                self.banner.display_success(f"Results copied to {dst_dir}")

            if extra_result_files:
                for extra in extra_result_files:
                    src_extra = work_dir / extra
                    if src_extra.exists():
                        shutil.copy2(src_extra, self.user_output_dir / extra)
                        self._log_info(f"Copied {extra} to output directory")

            for ext in ['*.html', '*.tsv', '*.json']:
                for root_file in work_dir.glob(ext):
                    if root_file.is_file():
                        shutil.copy2(root_file, self.user_output_dir / root_file.name)
                        self._log_info(f"Copied {root_file.name} to output directory")

            return True
        except Exception as e:
            self.banner.display_error(f"Exception in {module_name}: {e}")
            self._log_error(traceback.format_exc())
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self._log_info(f"Removed temporary directory: {temp_dir}")

    def run_fasta_qc_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "ecolityper_fasta_qc.py", pattern, "-o", "ecolityper_qc_results", "-c", str(threads)]
        extra = ["FASTA_QC_summary.html", "FASTA_QC_summary.tsv", "FASTA_QC_summary.json"]
        return self.run_module_in_temp("fasta_qc_module", fasta_files, cmd,
                                       "ecolityper_qc_results", "fasta_qc_results", extra)

    def run_mlst_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "ecolimlst_module.py", "-i", pattern, "-o", "ecolityper_mlst_results",
               "-db", "db", "-sc", "bin", "--batch"]
        extra = ["mlst_summary.html", "mlst_summary.tsv", "mlst_summary.json"]
        return self.run_module_in_temp("mlst_module", fasta_files, cmd,
                                       "ecolityper_mlst_results", "mlst_results", extra)

    def run_serotyping_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "enhanced_serotypefinder.py", "-i", pattern, "-o", "ecolityper_serotype_results"]
        extra = ["serotype_analysis_report.html", "serotype_analysis_report.tsv", "serotype_analysis_report.json"]
        return self.run_module_in_temp("serotypefinder_module", fasta_files, cmd,
                                       "ecolityper_serotype_results/SerotypeFinder_results",
                                       "serotyping_results", extra)

    def run_chtyper_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "enhanced_chtyper.py", "-i", pattern, "-o", "ecolityper_chtyper_results"]
        extra = ["chtyper_results.html", "chtyper_results.tsv", "chtyper_summary.json"]
        return self.run_module_in_temp("chtyper_module", fasta_files, cmd,
                                       "ecolityper_chtyper_results/chtyper_results",
                                       "chtyper_results", extra)

    def run_phylogrouping_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "enhanced_ezclermont.py", "-i", pattern, "-o", "ecolityper_phylo_results"]
        extra = ["phylogrouping_results.html", "phylogrouping_results.tsv", "phylogrouping_results.json"]
        return self.run_module_in_temp("phylogrouping_module", fasta_files, cmd,
                                       "ecolityper_phylo_results/phylogrouping_results",
                                       "phylogrouping_results", extra)

    def run_abricate_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int,
                              min_identity: int = 80, min_coverage: int = 80) -> bool:
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "ecoli_abricate.py", pattern]
        if min_identity is not None:
            cmd.extend(["--minid", str(min_identity)])
        if min_coverage is not None:
            cmd.extend(["--mincov", str(min_coverage)])
        extra = [
            "ecoli_abricate_master_summary.json",
            "ecoli_argannot_abricate_summary.tsv", "ecoli_argannot_summary.json", "ecoli_argannot_summary_report.html",
            "ecoli_bacmet2_abricate_summary.tsv", "ecoli_bacmet2_summary.json", "ecoli_bacmet2_summary_report.html",
            "ecoli_card_abricate_summary.tsv", "ecoli_card_summary.json", "ecoli_card_summary_report.html",
            "ecoli_ecoh_abricate_summary.tsv", "ecoli_ecoh_summary.json", "ecoli_ecoh_summary_report.html",
            "ecoli_ecoli_vf_abricate_summary.tsv", "ecoli_ecoli_vf_summary.json", "ecoli_ecoli_vf_summary_report.html",
            "ecoli_megares_abricate_summary.tsv", "ecoli_megares_summary.json", "ecoli_megares_summary_report.html",
            "ecoli_resfinder_abricate_summary.tsv", "ecoli_resfinder_summary.json", "ecoli_resfinder_summary_report.html",
            "ecoli_ncbi_abricate_summary.tsv", "ecoli_ncbi_summary.json", "ecoli_ncbi_summary_report.html",
            "ecoli_plasmidfinder_abricate_summary.tsv", "ecoli_plasmidfinder_summary.json", "ecoli_plasmidfinder_summary_report.html",
            "ecoli_vfdb_abricate_summary.tsv", "ecoli_vfdb_summary.json", "ecoli_vfdb_summary_report.html"
        ]
        return self.run_module_in_temp("abricate_module", fasta_files, cmd,
                                       "ecoli_abricate_results", "abricate_results", extra)

    def run_amrfinder_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int,
                               min_identity: float = None, min_coverage: float = None,
                               skip_mutations: bool = False, force_update: bool = False) -> bool:
        if force_update:
            self.update_amr_database(force=True)
        elif not self.ensure_amr_database():
            self.banner.display_error("AMR database is missing and could not be updated automatically.")
            return False

        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "ecoli_amrfinder.py", pattern]
        if min_identity is not None:
            cmd.extend(["--min-identity", str(min_identity)])
        if min_coverage is not None:
            cmd.extend(["--min-coverage", str(min_coverage)])
        if skip_mutations:
            cmd.append("--skip-mutations")
        extra = [
            "ecoli_amrfinder_summary_report.html", "ecoli_amrfinder_summary.tsv", "ecoli_amrfinder_statistics_summary.tsv",
            "ecoli_amrfinder_master_summary.json", "mutation_summary.html", "mutation_summary.tsv", "mutation_master_summary.json"
        ]
        return self.run_module_in_temp("amrfinder_module", fasta_files, cmd,
                                       "ecoli_amrfinder_results", "amrfinder_results", extra)

    def _copy_required_files_to_temp(self, src_output_dir: Path, temp_work_dir: Path):
        required_files = [
            ("fasta_qc_results", ["FASTA_QC_summary.html", "FASTA_QC_summary.tsv", "FASTA_QC_summary.json"]),
            ("mlst_results", ["mlst_summary.html", "mlst_summary.tsv", "mlst_summary.json"]),
            ("serotyping_results", ["serotype_analysis_report.html", "serotype_analysis_report.tsv", "serotype_analysis_report.json"]),
            ("chtyper_results", ["chtyper_results.html", "chtyper_results.tsv", "chtyper_summary.json"]),
            ("phylogrouping_results", ["phylogrouping_results.html", "phylogrouping_results.tsv", "phylogrouping_results.json"]),
            ("abricate_results", [
                "ecoli_argannot_abricate_summary.tsv", "ecoli_bacmet2_abricate_summary.tsv",
                "ecoli_card_abricate_summary.tsv", "ecoli_ecoh_abricate_summary.tsv",
                "ecoli_ecoli_vf_abricate_summary.tsv", "ecoli_megares_abricate_summary.tsv",
                "ecoli_ncbi_abricate_summary.tsv", "ecoli_plasmidfinder_abricate_summary.tsv",
                "ecoli_vfdb_abricate_summary.tsv", "ecoli_resfinder_abricate_summary.tsv",
                "ecoli_argannot_summary_report.html", "ecoli_bacmet2_summary_report.html",
                "ecoli_card_summary_report.html", "ecoli_ecoh_summary_report.html",
                "ecoli_ecoli_vf_summary_report.html", "ecoli_megares_summary_report.html",
                "ecoli_ncbi_summary_report.html", "ecoli_plasmidfinder_summary_report.html",
                "ecoli_vfdb_summary_report.html", "ecoli_resfinder_summary_report.html"
            ]),
            ("amrfinder_results", [
                "ecoli_amrfinder_summary_report.html", "ecoli_amrfinder_summary.tsv",
                "ecoli_amrfinder_statistics_summary.tsv", "mutation_summary.html", "mutation_summary.tsv"
            ]),
        ]
        for dirname, filenames in required_files:
            src_dir = src_output_dir / dirname
            if not src_dir.exists():
                self._log_warning(f"Source directory not found: {src_dir}")
                continue
            for fname in filenames:
                src_file = src_dir / fname
                if src_file.exists():
                    shutil.copy2(src_file, temp_work_dir / fname)
                    self._log_info(f"Copied {fname} to temporary directory")
                else:
                    self._log_warning(f"Required file not found: {src_file}")

    def run_sample_centric_analysis(self, output_dir: Path) -> bool:
        module_name = "sample_centric_module"
        module_orig = self.base_dir / "modules" / module_name
        if not module_orig.exists():
            self.banner.display_error(f"Module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix=f"ecolityper_{module_name}_"))
        try:
            shutil.copytree(module_orig, temp_dir / module_name, dirs_exist_ok=True)
            work_dir = temp_dir / module_name

            self._copy_required_files_to_temp(output_dir, work_dir)

            cmd = [sys.executable, "genius_ecoli_ultimate_sample_centric_reporter.py", "-i", "."]
            self.banner.display_info("Running Sample-Centric Hybrid Reporter...")
            result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)

            if result.returncode != 0:
                self.banner.display_error(f"Sample-centric reporter failed: {result.stderr}")
                if result.stdout:
                    self._log_info(result.stdout)
                return False

            results_dir = work_dir / "GENIUS_ECOLI_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
            if results_dir.exists():
                dst = output_dir / "GENIUS_ECOLI_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(results_dir, dst)
                self.banner.display_success(f"Sample-centric reports copied to {dst}")
            return True

        except Exception as e:
            self.banner.display_error(f"Sample-centric module error: {e}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_summary_analysis(self, output_dir: Path) -> bool:
        module_name = "gene_centric_module"
        module_orig = self.base_dir / "modules" / module_name
        if not module_orig.exists():
            self.banner.display_error(f"Module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix=f"ecolityper_{module_name}_"))
        try:
            shutil.copytree(module_orig, temp_dir / module_name, dirs_exist_ok=True)
            work_dir = temp_dir / module_name

            self._copy_required_files_to_temp(output_dir, work_dir)

            cmd = [sys.executable, "genius_ultimate_gene_centric_reporter.py", "-i", "."]
            result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)

            if result.returncode != 0:
                self.banner.display_error(f"Summary report failed: {result.stderr}")
                if result.stdout:
                    self._log_info(result.stdout)
                return False

            results_dir = work_dir / "GENIUS_ECOLI_ULTIMATE_GENE_CENTRIC_REPORTS"
            if results_dir.exists():
                dst = output_dir / "GENIUS_ECOLI_ULTIMATE_GENE_CENTRIC_REPORTS"
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(results_dir, dst)
                self.banner.display_success(f"Gene-centric reports copied to {dst}")
            return True

        except Exception as e:
            self.banner.display_error(f"Gene-centric module error: {e}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_visualization_analysis(self, output_dir: Path) -> bool:
        module_name = "visualization_module"
        module_orig = self.base_dir / "modules" / module_name
        if not module_orig.exists():
            self.banner.display_error(f"Module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix=f"ecolityper_{module_name}_"))
        try:
            shutil.copytree(module_orig, temp_dir / module_name, dirs_exist_ok=True)
            work_dir = temp_dir / module_name

            self._copy_required_files_to_temp(output_dir, work_dir)

            cmd = [sys.executable, "visualization_reporter.py"]
            result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)

            if result.returncode != 0:
                self.banner.display_warning(f"Visualization had issues: {result.stderr}")

            viz_dir = work_dir / "ECOLI_VISUALIZATIONS"
            if viz_dir.exists():
                dst = output_dir / "visualization_results"
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(viz_dir, dst / "ECOLI_VISUALIZATIONS")
                self.banner.display_success(f"Visualizations copied to {dst}")
            return True

        except Exception as e:
            self.banner.display_error(f"Visualization error: {e}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_lineage_analysis(self, output_dir: Path) -> bool:
        module_name = "ecoli_lineage"
        module_orig = self.base_dir / "modules" / module_name
        if not module_orig.exists():
            self.banner.display_error(f"Lineage module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix="ecolityper_lineage_"))
        try:
            shutil.copytree(module_orig, temp_dir / module_name, dirs_exist_ok=True)
            work_dir = temp_dir / module_name

            cmd = [sys.executable, "ecoli_html_reference.py"]
            result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)

            if result.returncode != 0:
                self.banner.display_warning(f"Lineage generation had warnings: {result.stderr}")

            html_src = work_dir / "ecoli_comprehensive_reference.html"
            if html_src.exists():
                dst_dir = output_dir / "lineage_results"
                dst_dir.mkdir(exist_ok=True)
                shutil.copy2(html_src, dst_dir / "ecoli_comprehensive_reference.html")
                self.banner.display_success("Lineage reference generated")
            return True

        except Exception as e:
            self.banner.display_error(f"Lineage error: {e}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_sequential_analyses(self, fasta_files: List[Path], output_dir: Path, threads: int,
                                skip_modules: Dict[str, bool],
                                amr_min_identity: float, amr_min_coverage: float,
                                amr_skip_mutations: bool, amr_force_update: bool,
                                abricate_min_identity: int, abricate_min_coverage: int) -> Dict[str, bool]:
        results = {}

        if not skip_modules.get('fasta_qc', False) and not self.interrupted:
            self.banner.display_module_header("FASTA QC", "Sequence Quality Control & Statistics")
            results["FASTA QC"] = self.run_fasta_qc_analysis(fasta_files, output_dir, threads)
            if results["FASTA QC"]:
                self.banner.display_success("✅ FASTA QC completed successfully!")
            else:
                self.banner.display_warning("⚠️ FASTA QC completed with issues")

        if not skip_modules.get('mlst', False) and not self.interrupted:
            self.banner.display_module_header("MLST", "Multi-Locus Sequence Typing")
            results["MLST"] = self.run_mlst_analysis(fasta_files, output_dir, threads)
            if results["MLST"]:
                self.banner.display_success("✅ MLST completed successfully!")
            else:
                self.banner.display_warning("⚠️ MLST completed with issues")

        if not skip_modules.get('serotyping', False) and not self.interrupted:
            self.banner.display_module_header("Serotyping", "O and H antigen determination")
            results["Serotyping"] = self.run_serotyping_analysis(fasta_files, output_dir, threads)
            if results["Serotyping"]:
                self.banner.display_success("✅ Serotyping completed successfully!")
            else:
                self.banner.display_warning("⚠️ Serotyping completed with issues")

        if not skip_modules.get('chtyper', False) and not self.interrupted:
            self.banner.display_module_header("CH Typing", "FumC and FimH typing")
            results["CH Typing"] = self.run_chtyper_analysis(fasta_files, output_dir, threads)
            if results["CH Typing"]:
                self.banner.display_success("✅ CH Typing completed successfully!")
            else:
                self.banner.display_warning("⚠️ CH Typing completed with issues")

        if not skip_modules.get('phylogrouping', False) and not self.interrupted:
            self.banner.display_module_header("Phylogrouping", "ezClermont phylogrouping algorithm")
            results["Phylogrouping"] = self.run_phylogrouping_analysis(fasta_files, output_dir, threads)
            if results["Phylogrouping"]:
                self.banner.display_success("✅ Phylogrouping completed successfully!")
            else:
                self.banner.display_warning("⚠️ Phylogrouping completed with issues")

        if not skip_modules.get('abricate', False) and not self.interrupted:
            self.banner.display_module_header("ABRicate", "Resistance, Virulence, and Plasmid gene screening")
            results["ABRicate"] = self.run_abricate_analysis(fasta_files, output_dir, threads,
                                                             abricate_min_identity, abricate_min_coverage)
            if results["ABRicate"]:
                self.banner.display_success("✅ ABRicate completed successfully!")
            else:
                self.banner.display_warning("⚠️ ABRicate completed with issues")

        if not skip_modules.get('amrfinder', False) and not self.interrupted:
            self.banner.display_module_header("AMRfinderPlus", "NCBI AMR gene detection with optional thresholds and mutation reporting")
            results["AMRfinderPlus"] = self.run_amrfinder_analysis(fasta_files, output_dir, threads,
                                                                   amr_min_identity, amr_min_coverage,
                                                                   amr_skip_mutations, amr_force_update)
            if results["AMRfinderPlus"]:
                self.banner.display_success("✅ AMRfinderPlus completed successfully!")
            else:
                self.banner.display_warning("⚠️ AMRfinderPlus completed with issues")

        return results

    def run_complete_analysis(self, input_path: str, output_dir: str, threads: int = 1,
                              skip_modules: Dict[str, bool] = None,
                              update_amr_db_only: bool = False,
                              force_update_amr_db: bool = False,
                              amr_min_identity: float = None,
                              amr_min_coverage: float = None,
                              amr_skip_mutations: bool = False,
                              amr_force_update: bool = False,
                              abricate_min_identity: int = 80,
                              abricate_min_coverage: int = 80,
                              keep_temp: bool = False):
        if skip_modules is None:
            skip_modules = {}
        self.keep_temp = keep_temp

        if update_amr_db_only:
            self.update_amr_database(force=False)
            return
        if force_update_amr_db:
            self.update_amr_database(force=True)
            return

        self.start_time = time.time()
        try:
            self.banner.display_startup_sequence()
            self.banner.display_banner(show_quote=True, show_author=True)

            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            self.user_output_dir = output_path
            self.setup_logging(output_path)

            self.fasta_files = self.find_fasta_files(input_path)
            if not self.fasta_files:
                self.banner.display_error("No FASTA files found! Analysis stopped.")
                return

            extensions = set(f.suffix.lower() for f in self.fasta_files)
            self.banner.display_success(f"Starting analysis of {len(self.fasta_files)} E. coli genomes")
            self.banner.display_info(f"File formats detected: {', '.join(extensions)}")

            self.banner.display_module_header("Analysis Plan", "Modules to be executed")
            plan = [
                ("FASTA QC", not skip_modules.get('fasta_qc', False)),
                ("MLST", not skip_modules.get('mlst', False)),
                ("Serotyping", not skip_modules.get('serotyping', False)),
                ("CH Typing", not skip_modules.get('chtyper', False)),
                ("Phylogrouping", not skip_modules.get('phylogrouping', False)),
                ("ABRicate", not skip_modules.get('abricate', False)),
                ("AMRfinderPlus", not skip_modules.get('amrfinder', False)),
                ("Sample-Centric Reporter", not skip_modules.get('samplecentric', False)),
                ("Gene-Centric Reporter", not skip_modules.get('summary', False)),
                ("Lineage Reference", not skip_modules.get('lineage', False)),
                ("Visualizations", not skip_modules.get('visualization', False))
            ]
            for analysis, enabled in plan:
                status = "✅ ENABLED" if enabled else "⏸️  SKIPPED"
                print(f"   {status} - {analysis}")

            analysis_results = self.run_sequential_analyses(
                self.fasta_files, output_path, threads, skip_modules,
                amr_min_identity, amr_min_coverage, amr_skip_mutations, amr_force_update,
                abricate_min_identity, abricate_min_coverage
            )

            if not skip_modules.get('samplecentric', False) and not self.interrupted:
                self.banner.display_module_header("Sample-Centric Reporter", "Interactive Isolate Boxes for AMR, Virulence, Plasmids, Bacmet & Mutations")
                self.run_sample_centric_analysis(output_path)

            if not skip_modules.get('summary', False) and not self.interrupted:
                self.banner.display_module_header("Gene-Centric Reporter", "HTML summary reports (GENIUS Reporter)")
                self.run_summary_analysis(output_path)

            if not skip_modules.get('lineage', False) and not self.interrupted:
                self.banner.display_module_header("Lineage Reference", "E. coli lineage reference database")
                self.run_lineage_analysis(output_path)

            if not skip_modules.get('visualization', False) and not self.interrupted:
                self.banner.display_module_header("Visualizations", "Charts and visualizations")
                self.run_visualization_analysis(output_path)

            self.total_duration = time.time() - self.start_time
            self.banner.display_footer(analysis_time=self._format_duration(self.total_duration),
                                       samples_processed=len(self.fasta_files))
            self.banner.display_success("🎉 Analysis completed successfully!")
        except KeyboardInterrupt:
            self.banner.display_error("Analysis interrupted by user")
        except Exception as e:
            self.banner.display_error(f"Critical error: {str(e)}")
            traceback.print_exc()
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description=f"{Colors.BRIGHT_YELLOW}EcoliTyper: Complete E. coli Typing Pipeline (Temp‑dir version){Colors.RESET}",
        formatter_class=ColoredHelpFormatter,
        epilog=f"""
{Colors.BOLD}{Colors.BRIGHT_GREEN}EXAMPLES:{Colors.RESET}

  {Colors.BRIGHT_CYAN}ecolityper -i genome.fna -o results/{Colors.RESET}
  {Colors.BRIGHT_CYAN}ecolityper -i "*.fna" -o batch_results --threads 8{Colors.RESET}
  {Colors.BRIGHT_CYAN}ecolityper -i "*.fasta" -o analysis --threads 16 --skip-lineage{Colors.RESET}
  {Colors.BRIGHT_CYAN}ecolityper -i "genome*.fa" -o results/ --threads 4{Colors.RESET}
  {Colors.BRIGHT_CYAN}ecolityper --update-amr-db{Colors.RESET}                # Update AMR database (incremental)
  {Colors.BRIGHT_CYAN}ecolityper --force-update-amr-db{Colors.RESET}          # Force full AMR database update
  {Colors.BRIGHT_CYAN}ecolityper -i "*.fna" -o results --amr-min-identity 0.95 --amr-min-coverage 0.9 --skip-amr-mutations{Colors.RESET}
  {Colors.BRIGHT_CYAN}ecolityper -i "*.fna" -o results --abricate-minid 90 --abricate-mincov 85{Colors.RESET}

{Colors.BOLD}{Colors.BRIGHT_GREEN}SUPPORTED FASTA FORMATS:{Colors.RESET} {Colors.YELLOW}.fna, .fasta, .fa, .fsa{Colors.RESET}

{Colors.BOLD}{Colors.BRIGHT_GREEN}ANALYSIS MODULES:{Colors.RESET}
  {Colors.GREEN}• FASTA QC (Sequence Quality Control & Statistics){Colors.RESET}
  {Colors.GREEN}• MLST (Multi-Locus Sequence Typing){Colors.RESET}
  {Colors.GREEN}• Serotyping (O and H antigen determination){Colors.RESET}
  {Colors.GREEN}• CH Typing (FumC and FimH typing){Colors.RESET}
  {Colors.GREEN}• Phylogrouping (ezClermont algorithm){Colors.RESET}
  {Colors.GREEN}• ABRicate (Resistance/Virulence/Plasmid screening) – with configurable identity/coverage thresholds{Colors.RESET}
  {Colors.GREEN}• AMRfinderPlus (NCBI AMR gene detection) – with optional thresholds and mutation reporting{Colors.RESET}
  {Colors.GREEN}• Sample-Centric Reporter (Interactive Isolate Boxes for AMR, Virulence, Plasmids, Bacmet & Mutations){Colors.RESET}
  {Colors.GREEN}• Gene-Centric Reporter (Integrated summary reports){Colors.RESET}
  {Colors.GREEN}• Lineage reference database{Colors.RESET}
  {Colors.GREEN}• Visualizations (Charts and visualizations){Colors.RESET}

{Colors.BOLD}{Colors.BRIGHT_GREEN}PRIOR SETUP REQUIRED:{Colors.RESET}
  {Colors.YELLOW}• abricate --setupdb{Colors.RESET}
  {Colors.YELLOW}• ecolityper --update-amr-db{Colors.RESET}

{Colors.BOLD}{Colors.BRIGHT_GREEN}OUTPUT:{Colors.RESET}
  Comprehensive results for all analyses in organized directories.
  A detailed log file `ecolityper_run.log` is written to the output directory.

{Colors.BOLD}{Colors.BRIGHT_GREEN}CONTACT:{Colors.RESET}
  Brown Beckley <brownbeckley94@gmail.com>
  University of Ghana Medical School - Department of Medical Biochemistry

  {Colors.BRIGHT_YELLOW}⭐ Drop a Star On Github If You Find This Tool Useful!!! ⭐{Colors.RESET}
"""
    )
    parser.add_argument('-i', '--input', help='Input FASTA file(s) - can use glob patterns like "*.fna"')
    parser.add_argument('-o', '--output', help='Output directory for all results')
    parser.add_argument('-t', '--threads', type=int, default=2, help='Number of threads (default: 2)')
    parser.add_argument('--keep-temp', action='store_true', help='Do not delete temporary directories (for debugging)')

    parser.add_argument('--update-amr-db', action='store_true', help='Update AMRfinderPlus database (incremental) and exit')
    parser.add_argument('--force-update-amr-db', action='store_true', help='Force complete AMR database update and exit')

    parser.add_argument('--amr-min-identity', type=float, help='Minimum identity for AMR hits (0..1)')
    parser.add_argument('--amr-min-coverage', type=float, help='Minimum coverage for AMR hits (0..1)')
    parser.add_argument('--skip-amr-mutations', action='store_true', help='Disable point mutation reporting in AMR (enabled by default)')
    parser.add_argument('--amr-force-update', action='store_true', help='Force update AMR database before analysis')

    parser.add_argument('--abricate-minid', type=int, default=80, help='Minimum identity for ABRicate hits (0-100, default: 80)')
    parser.add_argument('--abricate-mincov', type=int, default=80, help='Minimum coverage for ABRicate hits (0-100, default: 80)')

    skip_group = parser.add_argument_group(f'{Colors.BOLD}{Colors.BRIGHT_MAGENTA}Skip Options (disable specific analyses){Colors.RESET}')
    skip_group.add_argument('--skip-fasta-qc', action='store_true', help='Skip FASTA QC analysis')
    skip_group.add_argument('--skip-amrfinder', action='store_true', help='Skip AMRfinderPlus analysis')
    skip_group.add_argument('--skip-abricate', action='store_true', help='Skip ABRicate analysis')
    skip_group.add_argument('--skip-mlst', action='store_true', help='Skip MLST analysis')
    skip_group.add_argument('--skip-serotyping', action='store_true', help='Skip serotyping analysis')
    skip_group.add_argument('--skip-chtyper', action='store_true', help='Skip CH typing analysis')
    skip_group.add_argument('--skip-phylogrouping', action='store_true', help='Skip phylogrouping analysis')
    skip_group.add_argument('--skip-lineage', action='store_true', help='Skip lineage reference generation')
    skip_group.add_argument('--skip-summary', action='store_true', help='Skip gene-centric summary report generation')
    skip_group.add_argument('--skip-samplecentric', action='store_true', help='Skip sample-centric hybrid reporter generation')
    skip_group.add_argument('--skip-visualization', action='store_true', help='Skip visualization generation')

    args = parser.parse_args()

    # Validate thresholds
    if args.abricate_minid < 0 or args.abricate_minid > 100:
        parser.error("--abricate-minid must be between 0 and 100")
    if args.abricate_mincov < 0 or args.abricate_mincov > 100:
        parser.error("--abricate-mincov must be between 0 and 100")

    if args.update_amr_db or args.force_update_amr_db:
        orch = EcoliTyperOrchestrator()
        if args.force_update_amr_db:
            orch.update_amr_database(force=True)
        else:
            orch.update_amr_database(force=False)
        sys.exit(0)

    if not args.input or not args.output:
        parser.error("Both -i/--input and -o/--output are required for analysis (or use --update-amr-db / --force-update-amr-db)")

    skip_modules = {
        'fasta_qc': args.skip_fasta_qc,
        'amrfinder': args.skip_amrfinder,
        'abricate': args.skip_abricate,
        'mlst': args.skip_mlst,
        'serotyping': args.skip_serotyping,
        'chtyper': args.skip_chtyper,
        'phylogrouping': args.skip_phylogrouping,
        'lineage': args.skip_lineage,
        'summary': args.skip_summary,
        'samplecentric': args.skip_samplecentric,
        'visualization': args.skip_visualization
    }

    orchestrator = EcoliTyperOrchestrator()
    orchestrator.run_complete_analysis(
        input_path=args.input,
        output_dir=args.output,
        threads=args.threads,
        skip_modules=skip_modules,
        update_amr_db_only=False,
        force_update_amr_db=False,
        amr_min_identity=args.amr_min_identity,
        amr_min_coverage=args.amr_min_coverage,
        amr_skip_mutations=args.skip_amr_mutations,
        amr_force_update=args.amr_force_update,
        abricate_min_identity=args.abricate_minid,
        abricate_min_coverage=args.abricate_mincov,
        keep_temp=args.keep_temp
    )


if __name__ == "__main__":
    main()