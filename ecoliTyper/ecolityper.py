#!/usr/bin/env python3
"""
EcoliTyper Main Orchestrator - Complete E. coli Typing Pipeline
HPC-friendly / Docker-ready

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
MIT License
Version: 1.4.0
Date: 2026-10-24 / 2026-11-04
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
import urllib.request
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Set


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
    """ANSI-coloured argparse formatter for readable --help output."""

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
                        colored_line = (
                            f"  {Colors.BRIGHT_CYAN}{options}{Colors.RESET}  {help_text}"
                        )
                    else:
                        colored_line = (
                            f"  {Colors.BRIGHT_CYAN}{line.strip()}{Colors.RESET}"
                        )
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


# ----------------------------------------------------------------------
# Canonical output directory names
# ----------------------------------------------------------------------

DIR_FASTA_QC = "fasta_qc_results"
DIR_MLST = "mlst_results"
DIR_SEROTYPE = "serotype_results"
DIR_CHTYPER = "chtyper_results"
DIR_PHYLOGROUP = "phylogroup_results"
DIR_ABRICATE = "abricate_results"
DIR_AMRFINDER = "amrfinder_results"
DIR_MGE = "mge_results"
DIR_LINEAGE = "lineage_results"
DIR_VISUALIZATION = "visualization_results"
DIR_FINAL = "EcoliTyper_final_report"
DIR_GENE_CENTRIC = "GENIUS_ECOLI_ULTIMATE_GENE_CENTRIC_REPORTS"
DIR_SAMPLE_CENTRIC = "GENIUS_ECOLI_ULTIMATE_SAMPLE_CENTRIC_REPORTS"

ALL_MODULE_DIRS = [
    DIR_FASTA_QC, DIR_MLST, DIR_SEROTYPE, DIR_CHTYPER, DIR_PHYLOGROUP,
    DIR_ABRICATE, DIR_AMRFINDER, DIR_MGE, DIR_LINEAGE, DIR_VISUALIZATION,
]

CITATION_LINES = [
    "Beckley, B., Amarh, V. EcoliTyper: a species-optimized computational pipeline for comprehensive genotyping and surveillance of Escherichia coli. BMC Bioinformatics 27, 182 (2026).",
    "https://doi.org/10.1186/s12859-026-06529-6",
]


class EcoliTyperOrchestrator:
    """Main orchestrator for the EcoliTyper pipeline.

    Runs analysis modules in isolated temporary directories, collects
    batch summaries into a consolidated final report, and prunes
    intermediate files according to the selected report level.

    All subprocess stdout/stderr is captured in full and written to the
    run log. When a subprocess fails, its complete output is also echoed
    to the console so the failure is visible immediately.
    """

    def __init__(self):
        self.banner = EcoliTyperBanner()
        self.base_dir = Path(__file__).parent
        self.fasta_files: List[Path] = []
        self.interrupted = False
        self.output_lock = threading.Lock()
        self.start_time = None
        self.total_duration = None
        self.keep_temp = False
        self.user_output_dir: Optional[Path] = None
        self.logger = None
        self.report_level = 'summary'
        self.temp_dirs: List[str] = []

        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)

    # ------------------------------------------------------------------
    # Logging
    # ------------------------------------------------------------------

    def setup_logging(self, output_dir: Path):
        """Initialise file logging for the current run."""
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

    def _log_subprocess_result(self, name: str, result) -> bool:
        """Log full subprocess output; echo to console on failure.

        Writes the complete stdout and stderr to the run log regardless
        of the outcome. On non-zero exit code, both streams are also
        echoed to the console so the user sees the failure immediately
        without opening the log file.

        Returns True when returncode == 0, else False.
        """
        if result.stdout:
            self._log_info(f"=== {name} STDOUT (full) ===\n{result.stdout}")
        if result.stderr:
            self._log_warning(f"=== {name} STDERR (full) ===\n{result.stderr}")

        if result.returncode != 0:
            print(f"\n--- {name} FAILED (return code {result.returncode}) ---")
            if result.stdout:
                print(f"STDOUT:\n{result.stdout}")
            if result.stderr:
                print(f"STDERR:\n{result.stderr}")
            print(f"--- END {name} ---\n")
            return False
        return True

    def _format_duration(self, seconds: float) -> str:
        """Render a duration in seconds as a human-readable string."""
        if seconds < 60:
            return f"{seconds:.1f} seconds"
        elif seconds < 3600:
            return f"{int(seconds // 60)} minutes {seconds % 60:.0f} seconds"
        return (f"{int(seconds // 3600)} hours "
                f"{int((seconds % 3600) // 60)} minutes "
                f"{seconds % 60:.0f} seconds")

    def _signal_handler(self, signum, frame):
        """Remove temporary directories and exit on Ctrl-C / SIGTERM."""
        self.interrupted = True
        if self.start_time:
            self.total_duration = time.time() - self.start_time
            self._log_error(
                f"Total analysis ran for: {self._format_duration(self.total_duration)}"
            )
        self._log_error(f"Analysis interrupted by user (signal {signum})")
        self.banner.display_error("Analysis interrupted by user - cleaning up...")
        for temp_dir in self.temp_dirs:
            shutil.rmtree(temp_dir, ignore_errors=True)
        sys.exit(1)

    def _register_temp_dir(self, path: Path):
        """Track a temporary directory for SIGINT cleanup."""
        self.temp_dirs.append(str(path))

    # ------------------------------------------------------------------
    # Paths
    # ------------------------------------------------------------------

    def get_module_path(self, module_name: str) -> Path:
        """Resolve a module directory, honouring environment overrides."""
        env_override = os.environ.get(f'ECOLITYPER_{module_name.upper()}_MODULE_PATH')
        if env_override:
            return Path(env_override)
        if hasattr(sys, 'prefix'):
            share_path = (Path(sys.prefix) / "share" / "ecolityper"
                          / "modules" / module_name)
            if share_path.exists():
                return share_path
        return self.base_dir / "modules" / module_name

    def get_shared_data_root(self) -> Path:
        """Writable user-local data root for read-only module installs."""
        root = Path(os.environ.get(
            "ECOLITYPER_SHARED_DATA",
            Path.home() / ".local" / "share" / "ecolityper",
        ))
        root.mkdir(parents=True, exist_ok=True)
        return root

    def get_mlst_db_dir(self, for_update: bool = False) -> Path:
        """Return the MLST database directory.

        Resolution order:
          1. $ECOLITYPER_MLST_DB env override
          2. User-local module db/ if populated
          3. Packaged module db/ as fallback

        When for_update=True, always return a writable path.
        """
        env_db = os.environ.get("ECOLITYPER_MLST_DB")
        if env_db:
            env_path = Path(env_db)
            if for_update:
                env_path.mkdir(parents=True, exist_ok=True)
                return env_path
            if self._mlst_db_has_scheme(env_path):
                return env_path

        user_module = Path(os.environ.get(
            "ECOLITYPER_MLST_MODULE_PATH",
            Path.home() / ".local" / "share" / "ecolityper" / "mlst_module",
        ))
        user_db = user_module / "db"

        if for_update:
            user_db.mkdir(parents=True, exist_ok=True)
            return user_db

        if self._mlst_db_has_scheme(user_db):
            self.banner.display_info(f"Using user-local MLST database: {user_db}")
            return user_db

        module_db = self.get_module_path("mlst_module") / "db"
        if module_db.exists():
            return module_db

        self.banner.display_warning(
            "No populated MLST database found. Run `ecolityper --pull-mlst-db`."
        )
        user_db.mkdir(parents=True, exist_ok=True)
        return user_db

    @staticmethod
    def _mlst_db_has_scheme(db_dir: Path) -> bool:
        """Check whether the MLST db dir contains a usable E. coli scheme."""
        for scheme in ("ecoli_achtman_4", "ecoli"):
            if (db_dir / "pubmlst" / scheme).exists():
                return True
        return False

    # ------------------------------------------------------------------
    # MGE database resolution + management
    # ------------------------------------------------------------------

    def get_mge_data_dir(self, for_update: bool = False) -> Optional[Path]:
        """Return the mobileOG-db data directory.

        Resolution order:
          1. $ECOLITYPER_MGE_DATA_DIR env override
          2. User-local ~/.local/share/ecolityper/mge_module/mobileOG-db/
             (populated by `ecolityper --pull-mge-db`)
          3. Packaged <module>/mobileOG-db/ (Docker / Singularity)
        """
        env_dir = os.environ.get("ECOLITYPER_MGE_DATA_DIR")
        if env_dir:
            env_path = Path(env_dir)
            if for_update:
                env_path.mkdir(parents=True, exist_ok=True)
                return env_path
            if self._mge_db_has_data(env_path):
                return env_path

        user_dir = self.get_shared_data_root() / "mge_module" / "mobileOG-db"
        if for_update:
            user_dir.mkdir(parents=True, exist_ok=True)
            return user_dir
        if self._mge_db_has_data(user_dir):
            self._log_info(f"Using user-local mobileOG-db: {user_dir}")
            return user_dir

        module_db = self.get_module_path("mge_module") / "mobileOG-db"
        if self._mge_db_has_data(module_db):
            return module_db

        return None

    @staticmethod
    def _mge_db_has_data(db_dir: Optional[Path]) -> bool:
        """Return True if the directory contains the beatrix DIAMOND index."""
        if db_dir is None or not db_dir.exists():
            return False
        candidates = [
            db_dir / "beatrix-1-6_v1_all" / "mobileOG-db-beatrix-1.6.dmnd",
            db_dir / "mobileOG-db-beatrix-1.6.dmnd",
            db_dir / "beatrix-1-6_v1_all.dmnd",
        ]
        return any(c.exists() for c in candidates)

    def pull_mge_database(self) -> bool:
        """Download the mobileOG-db payload from the GitHub release.

        Docker and Singularity users do not need this — the data is baked
        into the image and get_mge_data_dir() finds it there automatically.
        Conda users run `ecolityper --pull-mge-db` once after install.
        """
        module_path = self.get_module_path("mge_module")

        if module_path.exists() and os.access(module_path, os.W_OK):
            target_root = module_path
            self.banner.display_info(f"Using writable module directory: {target_root}")
        else:
            target_root = self.get_shared_data_root() / "mge_module"
            target_root.mkdir(parents=True, exist_ok=True)
            self.banner.display_info(
                f"Module dir is read-only; using user-local path: {target_root}"
            )

        repo_url = os.environ.get(
            "ECOLITYPER_MGE_PULL_URL",
            "https://github.com/bbeckley-hub/EcoliTyper/releases/download/"
            "v1.4.0/ecolityper-mge-data-v1.4.0.zip",
        )
        temp_zip = Path(tempfile.mktemp(suffix=".zip"))
        extract_dir = None

        try:
            self.banner.display_info(f"Downloading mobileOG-db from {repo_url}")
            self.banner.display_info(
                "Payload is approximately 300 MB, this may take a few minutes..."
            )
            urllib.request.urlretrieve(repo_url, temp_zip)
            self.banner.display_success("Download complete")

            extract_dir = Path(tempfile.mkdtemp(prefix="ecolityper_mge_extract_"))
            with zipfile.ZipFile(temp_zip, "r") as zip_ref:
                zip_ref.extractall(extract_dir)
            self._log_info(f"Extracted to {extract_dir}")

            beatrix_dir = None
            for root, dirs, _ in os.walk(extract_dir):
                if "beatrix-1-6_v1_all" in dirs:
                    beatrix_dir = Path(root) / "beatrix-1-6_v1_all"
                    break

            if beatrix_dir is None:
                self.banner.display_error(
                    "Could not find 'beatrix-1-6_v1_all' directory "
                    "in downloaded zip"
                )
                return False

            dmnd = beatrix_dir / "mobileOG-db-beatrix-1.6.dmnd"
            if not dmnd.exists():
                self.banner.display_error(
                    f"DIAMOND index missing from archive: {dmnd.name}"
                )
                return False

            target_beatrix = target_root / "mobileOG-db" / "beatrix-1-6_v1_all"
            if target_beatrix.exists():
                shutil.rmtree(target_beatrix)
            target_beatrix.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(beatrix_dir, target_beatrix)

            dmnd_target = target_beatrix / "mobileOG-db-beatrix-1.6.dmnd"
            if not dmnd_target.exists():
                self.banner.display_error(
                    f"Install verification failed: {dmnd_target} not present"
                )
                return False

            size_mb = dmnd_target.stat().st_size / (1024 * 1024)
            self.banner.display_success(
                f"mobileOG-db installed to {target_root / 'mobileOG-db'} "
                f"(DIAMOND index: {size_mb:.0f} MB)"
            )

            os.environ["ECOLITYPER_MGE_DATA_DIR"] = str(target_root / "mobileOG-db")
            os.environ["ECOLITYPER_MGE_MODULE_PATH"] = str(target_root)
            return True

        except Exception as e:
            self.banner.display_error(
                f"Failed to download/extract mobileOG-db: {e}"
            )
            return False
        finally:
            if temp_zip.exists():
                temp_zip.unlink()
            if extract_dir and extract_dir.exists():
                shutil.rmtree(extract_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # AMR database management
    # ------------------------------------------------------------------

    def update_amr_database(self, force: bool = False) -> bool:
        """Update the AMRFinderPlus database (incremental or full)."""
        amr_module_path = self.base_dir / "modules" / "amrfinder_module"
        amr_script = amr_module_path / "ecoli_amrfinder.py"
        if not amr_script.exists():
            self.banner.display_error(f"AMR script not found at: {amr_script}")
            return False
        self.banner.display_info("Updating AMRfinderPlus database...")
        cmd = [sys.executable, str(amr_script),
               "--force-update" if force else "--update-db"]
        result = subprocess.run(cmd, capture_output=True, text=True,
                                cwd=amr_module_path)
        if not self._log_subprocess_result("AMR database update", result):
            self.banner.display_error("AMR database update failed.")
            return False
        self.banner.display_success("AMR database updated successfully.")
        version_cmd = [sys.executable, str(amr_script), "--db-version"]
        version_result = subprocess.run(version_cmd, capture_output=True,
                                        text=True, cwd=amr_module_path)
        if version_result.returncode == 0:
            self.banner.display_info(
                f"New database version: {version_result.stdout.strip()}"
            )
        return True

    def ensure_amr_database(self) -> bool:
        """Verify the AMRFinderPlus database is present, updating if needed."""
        amr_module_path = self.base_dir / "modules" / "amrfinder_module"
        amr_script = amr_module_path / "ecoli_amrfinder.py"
        if not amr_script.exists():
            self.banner.display_error("AMR script not found, cannot check database.")
            return False
        cmd = [sys.executable, str(amr_script), "--db-version"]
        result = subprocess.run(cmd, capture_output=True, text=True,
                                cwd=amr_module_path)
        if (result.returncode == 0
                and "Unknown" not in result.stdout
                and "No database" not in result.stdout):
            self.banner.display_success(
                f"AMR database already present: {result.stdout.strip()}"
            )
            return True
        self.banner.display_warning(
            "AMR database not found or outdated. Attempting automatic update..."
        )
        return self.update_amr_database(force=False)

    # ------------------------------------------------------------------
    # MLST database management
    # ------------------------------------------------------------------

    def pull_mlst_database(self) -> bool:
        """Download the prebuilt E. coli MLST database from GitHub."""
        module_path = self.get_module_path("mlst_module")
        system_path = module_path

        if os.access(module_path, os.W_OK):
            target_path = module_path
            self.banner.display_info(f"Using writable module directory: {target_path}")
        else:
            target_path = Path(os.environ.get(
                "ECOLITYPER_MLST_MODULE_PATH",
                Path.home() / ".local" / "share" / "ecolityper" / "mlst_module",
            ))
            target_path.mkdir(parents=True, exist_ok=True)
            self.banner.display_info(
                f"System MLST module is read-only; using user-local path: {target_path}"
            )

        repo_url = "https://github.com/bbeckley-hub/mlst/archive/refs/heads/master.zip"
        temp_zip = Path(tempfile.mktemp(suffix=".zip"))
        extract_dir = None

        try:
            self.banner.display_info("Downloading MLST data from GitHub...")
            urllib.request.urlretrieve(repo_url, temp_zip)

            extract_dir = Path(tempfile.mkdtemp(prefix="ecolityper_mlst_extract_"))
            with zipfile.ZipFile(temp_zip, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)

            src_root = extract_dir / "mlst-master"
            if not src_root.exists():
                self.banner.display_error("Extracted repository root not found.")
                return False

            for dirname in ['bin', 'db', 'perl5', 'scripts']:
                src = src_root / dirname
                if src.exists():
                    dst = target_path / dirname
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                    self._log_info(f"Copied {dirname} to {dst}")
                else:
                    self.banner.display_warning(
                        f"'{dirname}' not found in repository, skipping."
                    )

            if target_path != system_path:
                for script in ['ecolimlst_module.py', 'mlst_database.sh']:
                    src_script = system_path / script
                    if src_script.exists():
                        shutil.copy2(src_script, target_path / script)
                        self._log_info(f"Copied {script} to {target_path / script}")

            target_db = target_path / "db"
            if not self._mlst_db_has_scheme(target_db):
                self.banner.display_error(
                    f"E. coli MLST scheme missing after pull. Expected "
                    f"pubmlst/ecoli_achtman_4 or pubmlst/ecoli under {target_db}"
                )
                return False

            self.banner.display_success(f"MLST module pulled to {target_path}")
            os.environ['ECOLITYPER_MLST_MODULE_PATH'] = str(target_path)
            os.environ['ECOLITYPER_MLST_DB'] = str(target_db)
            return True

        except Exception as e:
            self.banner.display_error(f"Failed to download/extract MLST data: {e}")
            return False
        finally:
            if temp_zip.exists():
                temp_zip.unlink()
            if extract_dir and extract_dir.exists():
                shutil.rmtree(extract_dir, ignore_errors=True)

    def update_mlst_database(self) -> bool:
        """Refresh the E. coli MLST scheme from PubMLST (API key required)."""
        module_path = self.get_module_path("mlst_module")
        if os.access(module_path, os.W_OK):
            target_module = module_path
        else:
            target_module = Path(os.environ.get(
                "ECOLITYPER_MLST_MODULE_PATH",
                Path.home() / ".local" / "share" / "ecolityper" / "mlst_module",
            ))

        db_dir = target_module / "db"
        script_dir = target_module / "bin"
        scheme_map = db_dir / "scheme_species_map.tab"

        if not scheme_map.exists():
            self.banner.display_error(
                f"MLST module tree is not populated at {target_module}"
            )
            self.banner.display_info(
                "Run `ecolityper --pull-mlst-db` first to install the module, "
                "then run `ecolityper --update-mlst-db` to refresh the scheme."
            )
            return False

        if not script_dir.exists():
            self.banner.display_error(f"MLST bin/ not found at {script_dir}")
            return False

        os.environ['ECOLITYPER_MLST_MODULE_PATH'] = str(target_module)
        os.environ['ECOLITYPER_MLST_DB'] = str(db_dir)

        try:
            from ecoliTyper.modules.mlst_module.ecolimlst_module import (
                update_mlst_database as _update,
            )
        except ImportError:
            sys.path.insert(0, str(self.base_dir / "modules" / "mlst_module"))
            from ecolimlst_module import update_mlst_database as _update

        self.banner.display_info(f"Using MLST database directory: {db_dir}")
        return _update(db_dir, script_dir, scheme="ecoli_achtman_4")

    # ------------------------------------------------------------------
    # File discovery
    # ------------------------------------------------------------------

    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Locate FASTA files from a glob, single file, or directory."""
        self.banner.display_info(f"Searching for files with pattern: {input_path}")
        extensions = ('.fna', '.fasta', '.fa', '.fsa')

        if '*' in input_path or '?' in input_path:
            matched_files = glob.glob(input_path)
            fasta_files = [Path(f) for f in matched_files
                           if Path(f).is_file()
                           and f.lower().endswith(extensions)
                           and not Path(f).name.startswith('.')]
            self.banner.display_success(f"Found {len(fasta_files)} FASTA files")
            return sorted(fasta_files)

        p = Path(input_path)
        if p.is_file() and p.suffix.lower() in extensions:
            self.banner.display_success(f"Found single FASTA file: {p.name}")
            return [p]

        if p.is_dir():
            fasta_files = []
            for ext in extensions:
                fasta_files.extend(p.glob(f"*{ext}"))
                fasta_files.extend(p.glob(f"*{ext.upper()}"))
            fasta_files = sorted(set(
                f for f in fasta_files
                if f.is_file() and not f.name.startswith('.')
            ))
            if fasta_files:
                self.banner.display_success(
                    f"Found {len(fasta_files)} FASTA files in directory"
                )
            else:
                self.banner.display_warning(
                    f"No FASTA files found in directory: {input_path}"
                )
            return fasta_files

        self.banner.display_error(f"Input path not found: {input_path}")
        return []

    def get_file_pattern(self, fasta_files: List[Path]) -> str:
        """Return a glob pattern that matches the extensions in use."""
        if not fasta_files:
            return "*.fna"
        extensions = set(f.suffix.lower() for f in fasta_files)
        if len(extensions) == 1:
            return f"*{list(extensions)[0]}"
        return "*"

    # ------------------------------------------------------------------
    # Species gate
    # ------------------------------------------------------------------

    def _load_species_gate(self, output_dir: Path) -> Optional[Set[str]]:
        """Read the species-gate pass list produced by FASTA QC."""
        gate_file = output_dir / DIR_FASTA_QC / "FASTA_QC_species_passed.txt"
        if not gate_file.exists():
            return None
        with open(gate_file, 'r', encoding='utf-8') as f:
            names = {line.strip() for line in f if line.strip()}
        self._log_info(
            f"Species gate loaded from {gate_file}: {len(names)} sample(s) passed"
        )
        return names

    def _filter_by_species_gate(self, fasta_files: List[Path],
                                gate: Optional[Set[str]]) -> List[Path]:
        """Keep only FASTA files that passed the species gate."""
        if gate is None:
            self.banner.display_warning(
                "Species gate not available (FASTA QC failed or was skipped). "
                "Running ALL samples through downstream modules."
            )
            return fasta_files

        passed = [f for f in fasta_files if f.name in gate]
        excluded = [f for f in fasta_files if f.name not in gate]

        if excluded:
            self.banner.display_info(
                f"Species gate: {len(passed)} of {len(fasta_files)} samples passed"
            )
            for f in excluded:
                self.banner.display_warning(
                    f"  Excluded (not E. coli/Shigella): {f.name}"
                )
        else:
            self.banner.display_success(
                f"Species gate: all {len(passed)} sample(s) passed species check"
            )
        return passed

    # ------------------------------------------------------------------
    # Copy helpers
    # ------------------------------------------------------------------

    def _copy_module_results(self, src: Path, dst: Path,
                             skip_dirs: Optional[Set[str]] = None) -> None:
        """Copy module results to user output, respecting report_level.

        Batch-level files (top-level files in src) are always copied.
        Per-sample folders (subdirectories) are only copied when
        report_level == 'full'. This avoids copy-then-delete churn on
        large batches.
        """
        if not src.exists():
            return
        if dst.exists():
            shutil.rmtree(dst)
        dst.mkdir(parents=True, exist_ok=True)

        skip_dirs = skip_dirs or set()
        for item in src.iterdir():
            if item.name in skip_dirs:
                continue
            if item.is_file():
                shutil.copy2(item, dst / item.name)
            elif item.is_dir() and self.report_level == 'full':
                shutil.copytree(item, dst / item.name)

    # ------------------------------------------------------------------
    # Generic module runner
    # ------------------------------------------------------------------

    def run_module_in_temp(self, module_name: str, fasta_files: List[Path],
                           cmd: List[str], result_subdir: str, target_dir: str,
                           extra_result_files: Optional[List[str]] = None,
                           env: Optional[Dict[str, str]] = None) -> bool:
        """Run a module inside an isolated temp directory and copy results back."""
        module_orig = self.base_dir / "modules" / module_name
        if not module_orig.exists():
            self.banner.display_error(f"Module directory not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix=f"ecolityper_{module_name}_"))
        self._register_temp_dir(temp_dir)
        self._log_info(f"Temporary directory for {module_name}: {temp_dir}")

        try:
            shutil.copytree(module_orig, temp_dir / module_name, dirs_exist_ok=True)
            work_dir = temp_dir / module_name
            for f in fasta_files:
                shutil.copy2(f, work_dir / f.name)

            self.banner.display_info(f"Running {module_name} analysis...")
            self._log_info(f"Command: {' '.join(cmd)}")

            result = subprocess.run(cmd, cwd=work_dir, capture_output=True,
                                    text=True, env=env)

            if not self._log_subprocess_result(module_name, result):
                self.banner.display_error(
                    f"{module_name} failed with return code {result.returncode}"
                )
                return False

            src_dir = work_dir / result_subdir
            if src_dir.exists():
                dst_dir = self.user_output_dir / target_dir
                self._copy_module_results(src_dir, dst_dir)
                self.banner.display_success(f"Results copied to {dst_dir}")

            if extra_result_files:
                for extra in extra_result_files:
                    src_extra = work_dir / extra
                    if src_extra.exists():
                        shutil.copy2(src_extra, self.user_output_dir / extra)
                        self._log_info(f"Copied {extra} to output directory")

            return True
        except Exception as e:
            self.banner.display_error(f"Exception in {module_name}: {e}")
            self._log_error(traceback.format_exc())
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self._log_info(f"Removed temporary directory: {temp_dir}")

    # ------------------------------------------------------------------
    # Individual modules
    # ------------------------------------------------------------------

    def run_fasta_qc_analysis(self, fasta_files: List[Path], output_dir: Path,
                              threads: int) -> bool:
        """Run the FASTA QC module (assembly statistics + species check)."""
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "ecolityper_fasta_qc.py", pattern,
               "-o", "ecolityper_qc_results", "-c", str(threads)]
        return self.run_module_in_temp(
            "fasta_qc_module", fasta_files, cmd,
            "ecolityper_qc_results", DIR_FASTA_QC,
        )

    def run_mlst_analysis(self, fasta_files: List[Path], output_dir: Path,
                          threads: int) -> bool:
        """Run the MLST module against the resolved E. coli scheme."""
        pattern = self.get_file_pattern(fasta_files)
        env = os.environ.copy()

        module_db = self.get_module_path("mlst_module") / "db"
        if module_db.exists():
            db_arg = "db"
        else:
            resolved_db = self.get_mlst_db_dir()
            if not resolved_db.exists():
                self.banner.display_error(
                    "No MLST database found. Run `ecolityper --pull-mlst-db` first."
                )
                return False
            db_arg = str(resolved_db)
            env['ECOLITYPER_MLST_DB'] = str(resolved_db)

        cmd = [sys.executable, "ecolimlst_module.py",
               "-i", pattern, "-o", "mlst_results",
               "-db", db_arg, "-sc", "bin"]
        return self.run_module_in_temp(
            "mlst_module", fasta_files, cmd,
            "mlst_results", DIR_MLST,
            env=env,
        )

    def run_serotyping_analysis(self, fasta_files: List[Path], output_dir: Path,
                                threads: int) -> bool:
        """Run the serotyping module (O and H antigen prediction)."""
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "enhanced_serotypefinder.py",
               "-i", pattern, "-o", "serotype_results", "-t", str(threads)]
        return self.run_module_in_temp(
            "serotypefinder_module", fasta_files, cmd,
            "serotype_results", DIR_SEROTYPE,
        )

    def run_chtyper_analysis(self, fasta_files: List[Path], output_dir: Path,
                              threads: int) -> bool:
        """Run the CH typing module (FumC and FimH alleles)."""
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "enhanced_chtyper.py",
               "-i", pattern, "-o", "chtyper_results", "-t", str(threads)]
        return self.run_module_in_temp(
            "chtyper_module", fasta_files, cmd,
            "chtyper_results", DIR_CHTYPER,
        )

    def run_phylogrouping_analysis(self, fasta_files: List[Path], output_dir: Path,
                                    threads: int) -> bool:
        """Run the ezClermont phylogrouping module."""
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "enhanced_ezclermont.py",
               "-i", pattern, "-o", "phylogroup_results", "-t", str(threads)]
        return self.run_module_in_temp(
            "phylogrouping_module", fasta_files, cmd,
            "phylogroup_results", DIR_PHYLOGROUP,
        )

    def run_abricate_analysis(self, fasta_files: List[Path], output_dir: Path,
                              threads: int, min_identity: int = 80,
                              min_coverage: int = 80) -> bool:
        """Run ABRicate across all configured databases."""
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "ecoli_abricate.py", pattern,
               "--minid", str(min_identity), "--mincov", str(min_coverage)]
        return self.run_module_in_temp(
            "abricate_module", fasta_files, cmd,
            "ecoli_abricate_results", DIR_ABRICATE,
        )

    def run_amrfinder_analysis(self, fasta_files: List[Path], output_dir: Path,
                               threads: int, min_identity: Optional[float] = None,
                               min_coverage: Optional[float] = None,
                               skip_mutations: bool = False,
                               force_update: bool = False) -> bool:
        """Run the AMRFinderPlus module with optional thresholds."""
        if force_update:
            self.update_amr_database(force=True)
        elif not self.ensure_amr_database():
            self.banner.display_error(
                "AMR database is missing and could not be updated automatically."
            )
            return False

        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "ecoli_amrfinder.py", pattern]
        if min_identity is not None:
            cmd.extend(["--min-identity", str(min_identity)])
        if min_coverage is not None:
            cmd.extend(["--min-coverage", str(min_coverage)])
        if skip_mutations:
            cmd.append("--skip-mutations")
        return self.run_module_in_temp(
            "amrfinder_module", fasta_files, cmd,
            "ecoli_amrfinder_results", DIR_AMRFINDER,
        )

    def run_mge_analysis(self, fasta_files: List[Path], output_dir: Path,
                         threads: int) -> bool:
        """Run the mobileOG-based MGE profiler.

        Resolves the mobileOG-db payload via get_mge_data_dir():
          - Docker/Singularity: found inside the packaged mge_module/
          - Conda + `--pull-mge-db`: found at ~/.local/share/ecolityper/
        Either way, the beatrix directory is symlinked (or copied as a
        fallback) into the temp working directory so MGEFinder.py sees the
        layout it expects:
            mobileOG-db/beatrix-1-6_v1_all/mobileOG-db-beatrix-1.6.dmnd
        """
        module_orig = self.get_module_path("mge_module")
        if not module_orig.exists():
            self.banner.display_error(f"mge_module not found: {module_orig}")
            return False

        mge_data = self.get_mge_data_dir()
        if mge_data is None:
            self.banner.display_error(
                "mobileOG-db not found. Run `ecolityper --pull-mge-db` once, "
                "or use the Docker/Singularity image (which has it bundled)."
            )
            return False
        self._log_info(f"Resolved mobileOG-db: {mge_data}")

        temp_dir = Path(tempfile.mkdtemp(prefix="ecolityper_mge_"))
        self._register_temp_dir(temp_dir)
        self._log_info(f"Temporary directory for mge: {temp_dir}")

        try:
            # Copy the module tree, excluding any existing mobileOG-db —
            # the resolved copy below is authoritative.
            shutil.copytree(
                module_orig, temp_dir, dirs_exist_ok=True,
                ignore=shutil.ignore_patterns(
                    "mge_results", "annotations", "*.fna", "*.fasta",
                    "mobileOG-db", "__pycache__", "*.pyc",
                ),
            )

            # Link or copy the resolved mobileOG-db into the temp tree
            temp_mge = temp_dir / "mobileOG-db"
            temp_mge.mkdir(parents=True, exist_ok=True)

            beatrix_src = mge_data / "beatrix-1-6_v1_all"
            if not beatrix_src.exists() and mge_data.name == "beatrix-1-6_v1_all":
                beatrix_src = mge_data
            if not beatrix_src.exists():
                self.banner.display_error(
                    f"beatrix-1-6_v1_all not found under {mge_data}"
                )
                return False

            beatrix_dst = temp_mge / "beatrix-1-6_v1_all"
            if beatrix_dst.exists() or beatrix_dst.is_symlink():
                if beatrix_dst.is_symlink() or beatrix_dst.is_file():
                    beatrix_dst.unlink()
                else:
                    shutil.rmtree(beatrix_dst)

            try:
                os.symlink(beatrix_src.resolve(), beatrix_dst)
                self._log_info(f"Symlinked {beatrix_src} -> {beatrix_dst}")
            except OSError as exc:
                self._log_warning(
                    f"Symlink failed ({exc}); falling back to copy"
                )
                shutil.copytree(beatrix_src, beatrix_dst)
                self._log_info(f"Copied {beatrix_src} -> {beatrix_dst}")

            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            self.banner.display_info(
                f"Running MGE analysis on {len(fasta_files)} genome(s)..."
            )
            script = temp_dir / "MGEFinder.py"
            if not script.exists():
                self.banner.display_error("MGEFinder.py not found in mge_module")
                return False

            mge_inputs = [f.name for f in fasta_files]
            cmd = [sys.executable, str(script), "-i"] + mge_inputs + [
                "-d", "mge_results",
                "-db_dir", "mobileOG-db/beatrix-1-6_v1_all",
            ]
            # Log the extension pattern rather than the full filename list.
            # MGEFinder still receives `cmd` with every filename.
            display_pattern = self.get_file_pattern(fasta_files)
            self._log_info(
                f"Command: MGEFinder.py -i {display_pattern} "
                f"-d mge_results -db_dir mobileOG-db/beatrix-1-6_v1_all"
            )

            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True,
                                    text=True)

            if not self._log_subprocess_result("MGE", result):
                self.banner.display_error(
                    f"MGE analysis failed with return code {result.returncode}"
                )
                return False

            src = temp_dir / "mge_results"
            if not src.exists():
                self.banner.display_error(
                    f"MGE results directory not found: {src}"
                )
                return False

            dst = self.user_output_dir / DIR_MGE
            self._copy_module_results(src, dst, skip_dirs={"annotations"})

            expected_files = [
                "ecolityper_mge_summary.tsv",
                "ecolityper_mge_all_hits.tsv",
                "ecolityper_mge_summary.json",
                "ecolityper_mge_summary.html",
            ]
            missing = [f for f in expected_files if not (dst / f).exists()]
            if missing:
                self.banner.display_error(
                    f"MGE ran but these summary files are missing: "
                    f"{', '.join(missing)}"
                )
                return False

            self.banner.display_success(
                f"MGE results copied to {dst} (annotations skipped)"
            )
            return True
        except Exception as e:
            self.banner.display_error(f"MGE exception: {e}")
            self._log_error(traceback.format_exc())
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Reporter inputs
    # ------------------------------------------------------------------

    def _copy_required_files_to_temp(self, src_output_dir: Path,
                                     temp_work_dir: Path) -> None:
        """Copy every batch file the two reporters require into a work dir."""
        required_files = [
            (DIR_FASTA_QC, ["FASTA_QC_summary.html", "FASTA_QC_summary.tsv",
                            "FASTA_QC_summary.json"]),
            (DIR_MLST, ["mlst_summary.html", "mlst_summary.tsv",
                        "mlst_summary.json"]),
            (DIR_SEROTYPE, ["serotype_analysis_report.html",
                            "serotype_analysis_report.tsv",
                            "serotype_analysis_report.json"]),
            (DIR_CHTYPER, ["chtyper_results.html", "chtyper_results.tsv",
                           "chtyper_summary.json"]),
            (DIR_PHYLOGROUP, ["phylogrouping_results.html",
                              "phylogrouping_results.tsv",
                              "phylogrouping_results.json"]),
            (DIR_ABRICATE, [
                "ecoli_argannot_abricate_summary.tsv",
                "ecoli_bacmet2_abricate_summary.tsv",
                "ecoli_card_abricate_summary.tsv",
                "ecoli_ecoh_abricate_summary.tsv",
                "ecoli_ecoli_vf_abricate_summary.tsv",
                "ecoli_megares_abricate_summary.tsv",
                "ecoli_ncbi_abricate_summary.tsv",
                "ecoli_plasmidfinder_abricate_summary.tsv",
                "ecoli_vfdb_abricate_summary.tsv",
                "ecoli_resfinder_abricate_summary.tsv",
                "ecoli_argannot_summary_report.html",
                "ecoli_bacmet2_summary_report.html",
                "ecoli_card_summary_report.html",
                "ecoli_ecoh_summary_report.html",
                "ecoli_ecoli_vf_summary_report.html",
                "ecoli_megares_summary_report.html",
                "ecoli_ncbi_summary_report.html",
                "ecoli_plasmidfinder_summary_report.html",
                "ecoli_vfdb_summary_report.html",
                "ecoli_resfinder_summary_report.html",
                "ecoli_abricate_master_summary.json",
            ]),
            (DIR_AMRFINDER, [
                "ecoli_amrfinder_summary_report.html",
                "ecoli_amrfinder_summary.tsv",
                "ecoli_amrfinder_statistics_summary.tsv",
                "ecoli_amrfinder_master_summary.json",
                "mutation_summary.html",
                "mutation_summary.tsv",
                "mutation_master_summary.json",
            ]),
            (DIR_MGE, [
                "ecolityper_mge_summary.html",
                "ecolityper_mge_all_hits.tsv",
                "ecolityper_mge_summary.tsv",
                "ecolityper_mge_summary.json",
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
                    self._log_info(f"Copied {fname} to reporter work dir")
                else:
                    self._log_warning(f"Required file not found: {src_file}")

    # ------------------------------------------------------------------
    # Comprehensive + gene-centric reports
    # ------------------------------------------------------------------

    def run_comprehensive_and_gene_centric(self, output_dir: Path) -> bool:
        """Run comprehensive_report.py then the gene-centric reporter.

        Both scripts execute in the same temporary directory so the
        gene-centric reporter finds the comprehensive TSV as its master
        typing input.
        """
        module_name = "gene_centric_module"
        module_orig = self.base_dir / "modules" / module_name
        if not module_orig.exists():
            self.banner.display_error(f"Module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix=f"ecolityper_{module_name}_"))
        self._register_temp_dir(temp_dir)
        try:
            shutil.copytree(module_orig, temp_dir / module_name, dirs_exist_ok=True)
            work_dir = temp_dir / module_name

            self._copy_required_files_to_temp(output_dir, work_dir)

            # --- Step 1: comprehensive report ---
            comp_script = work_dir / "comprehensive_report.py"
            if not comp_script.exists():
                self.banner.display_error("comprehensive_report.py not found")
                return False

            self.banner.display_info(
                "Step 1/2: Running comprehensive report "
                "(generates ecolityper_comprehensive_report.*)"
            )
            comp_result = subprocess.run(
                [sys.executable, str(comp_script)],
                cwd=work_dir, capture_output=True, text=True,
            )
            if not self._log_subprocess_result("Comprehensive report", comp_result):
                self.banner.display_error(
                    "Comprehensive report failed."
                )
                return False

            copied_comp = 0
            for ext in [".html", ".tsv", ".json"]:
                src = work_dir / f"ecolityper_comprehensive_report{ext}"
                if src.exists():
                    shutil.copy2(src, output_dir / src.name)
                    copied_comp += 1
            if copied_comp == 0:
                self.banner.display_error(
                    "Comprehensive report ran but produced no output files."
                )
                return False
            self.banner.display_success(
                f"Comprehensive report complete ({copied_comp}/3 files)"
            )

            # --- Step 2: gene-centric reporter ---
            gene_script = work_dir / "genius_ultimate_gene_centric_reporter.py"
            if not gene_script.exists():
                self.banner.display_error(
                    "genius_ultimate_gene_centric_reporter.py not found"
                )
                return False

            self.banner.display_info(
                "Step 2/2: Running gene-centric reporter "
                "(consumes comprehensive report)"
            )
            gene_result = subprocess.run(
                [sys.executable, str(gene_script), "-i", "."],
                cwd=work_dir, capture_output=True, text=True,
            )
            if not self._log_subprocess_result("Gene-centric reporter", gene_result):
                self.banner.display_error(
                    "Gene-centric reporter failed."
                )
                return False

            results_dir = work_dir / DIR_GENE_CENTRIC
            if not results_dir.exists():
                self.banner.display_error(
                    f"{DIR_GENE_CENTRIC} not produced by gene-centric reporter"
                )
                return False

            dst = output_dir / DIR_GENE_CENTRIC
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(results_dir, dst)
            self.banner.display_success(f"Gene-centric reports copied to {dst}")
            return True

        except Exception as e:
            self.banner.display_error(
                f"Comprehensive + gene-centric exception: {e}"
            )
            self._log_error(traceback.format_exc())
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Sample-centric reporter
    # ------------------------------------------------------------------

    def run_sample_centric_analysis(self, output_dir: Path) -> bool:
        """Run the sample-centric reporter with robust output discovery.

        The reporter is known to occasionally exit cleanly without producing
        output, write to a directory name different from the expected one,
        or bail silently when an input file is missing. This method:

          1. Logs full stdout/stderr to the run log (via the shared helper)
          2. Snapshots work_dir before/after to identify new files
          3. Finds any HTML file the reporter created, regardless of folder
          4. Preserves the temp dir on failure (for inspection)
        """
        module_name = "sample_centric_module"
        module_orig = self.base_dir / "modules" / module_name
        if not module_orig.exists():
            self.banner.display_error(f"Module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix=f"ecolityper_{module_name}_"))
        self._register_temp_dir(temp_dir)
        try:
            shutil.copytree(module_orig, temp_dir / module_name, dirs_exist_ok=True)
            work_dir = temp_dir / module_name

            self._copy_required_files_to_temp(output_dir, work_dir)
            for ext in [".html", ".tsv", ".json"]:
                src = output_dir / f"ecolityper_comprehensive_report{ext}"
                if src.exists():
                    shutil.copy2(src, work_dir / src.name)
                else:
                    self.banner.display_warning(
                        f"ecolityper_comprehensive_report{ext} not found - "
                        f"sample-centric reporter may fail"
                    )

            script = work_dir / "genius_ecoli_ultimate_sample_centric_reporter.py"
            if not script.exists():
                self.banner.display_error(
                    "genius_ecoli_ultimate_sample_centric_reporter.py not found"
                )
                return False

            def _snapshot(root: Path):
                return {p for p in root.rglob("*") if p.is_file()}

            before = _snapshot(work_dir)
            self._log_info(
                f"Input files in sample-centric work dir: "
                f"{sorted(p.name for p in before if p.parent == work_dir)}"
            )

            self.banner.display_info("Running Sample-Centric Reporter...")
            result = subprocess.run(
                [sys.executable, str(script), "-i", "."],
                cwd=work_dir, capture_output=True, text=True,
                input="",
            )

            self._log_info(f"Sample-centric rc={result.returncode}")
            if not self._log_subprocess_result("Sample-centric", result):
                self.banner.display_error(
                    f"Sample-centric reporter exited with rc={result.returncode}"
                )
                if not self.keep_temp:
                    self._log_info(
                        f"Preserving temp dir for debugging: {temp_dir}"
                    )
                    self.keep_temp = True
                return False

            after = _snapshot(work_dir)
            new_files = after - before

            self._log_info(
                f"Sample-centric produced {len(new_files)} new file(s) in work dir"
            )
            for f in sorted(new_files):
                self._log_info(f"  NEW: {f.relative_to(work_dir)}")

            candidates = []

            expected = work_dir / DIR_SAMPLE_CENTRIC
            if expected.exists() and expected.is_dir():
                candidates.append(expected)

            if not candidates:
                dirs_with_html = set()
                for f in new_files:
                    if f.suffix == ".html":
                        dirs_with_html.add(f.parent)
                for d in dirs_with_html:
                    if d != work_dir and d not in candidates:
                        candidates.append(d)
                        self._log_info(
                            f"Reporter wrote to unexpected dir: "
                            f"{d.relative_to(work_dir)}"
                        )

            if not candidates:
                for subdir in work_dir.iterdir():
                    if subdir.is_dir() and subdir.name not in {"__pycache__"}:
                        htmls = list(subdir.glob("*.html"))
                        if htmls:
                            candidates.append(subdir)
                            self._log_info(
                                f"Fallback: found HTMLs in "
                                f"{subdir.relative_to(work_dir)}"
                            )
                            break

            if not candidates:
                self.banner.display_error(
                    "Sample-centric reporter produced no recognisable output. "
                    "Work dir contents logged to ecolityper_run.log for inspection."
                )
                self._log_error("=== FULL WORK DIR CONTENTS AFTER RUN ===")
                for p in sorted(work_dir.rglob("*")):
                    if p.is_file():
                        self._log_error(
                            f"  {p.relative_to(work_dir)} "
                            f"({p.stat().st_size} bytes)"
                        )
                self._log_error("=== END WORK DIR CONTENTS ===")

                self.banner.display_warning(
                    f"Temp dir preserved for debugging: {temp_dir}\n"
                    f"   Inspect with: ls -R {temp_dir}\n"
                    f"   Then check: grep 'Sample-centric' "
                    f"{self.user_output_dir / 'ecolityper_run.log'}"
                )
                self.keep_temp = True
                return False

            results_dir = candidates[0]
            dst = output_dir / DIR_SAMPLE_CENTRIC
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(results_dir, dst)

            html_count = len(list(dst.glob("*.html")))
            self.banner.display_success(
                f"Sample-centric reports copied to {dst} "
                f"({html_count} HTML files)"
            )
            return True

        except Exception as e:
            self.banner.display_error(f"Sample-centric reporter exception: {e}")
            self._log_error(traceback.format_exc())
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Lineage + visualization
    # ------------------------------------------------------------------

    def run_visualization_analysis(self, output_dir: Path) -> bool:
        """Run the visualization module (dashboard + static plots)."""
        module_name = "visualization_module"
        module_orig = self.base_dir / "modules" / module_name
        if not module_orig.exists():
            self.banner.display_error(f"Module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix=f"ecolityper_{module_name}_"))
        self._register_temp_dir(temp_dir)
        try:
            shutil.copytree(module_orig, temp_dir / module_name, dirs_exist_ok=True)
            work_dir = temp_dir / module_name

            self._copy_required_files_to_temp(output_dir, work_dir)

            for ext in [".html", ".tsv", ".json"]:
                src = output_dir / f"ecolityper_comprehensive_report{ext}"
                if src.exists():
                    shutil.copy2(src, work_dir / src.name)

            script = work_dir / "visualization_reporter.py"
            if not script.exists():
                self.banner.display_error("visualization_reporter.py not found")
                return False

            result = subprocess.run(
                [sys.executable, str(script)],
                cwd=work_dir, capture_output=True, text=True,
            )
            if not self._log_subprocess_result("Visualization", result):
                self.banner.display_warning(
                    f"Visualization had issues (rc={result.returncode})"
                )

            viz_dir = work_dir / "ECOLI_VISUALIZATIONS"
            if viz_dir.exists():
                dst = output_dir / DIR_VISUALIZATION
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(viz_dir, dst)
                self.banner.display_success(f"Visualizations copied to {dst}")
                return True
            return False

        except Exception as e:
            self.banner.display_error(f"Visualization error: {e}")
            self._log_error(traceback.format_exc())
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_lineage_analysis(self, output_dir: Path) -> bool:
        """Generate the E. coli lineage reference HTML."""
        module_name = "ecoli_lineage"
        module_orig = self.base_dir / "modules" / module_name
        if not module_orig.exists():
            self.banner.display_error(f"Module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix="ecolityper_lineage_"))
        self._register_temp_dir(temp_dir)
        try:
            shutil.copytree(module_orig, temp_dir / module_name, dirs_exist_ok=True)
            work_dir = temp_dir / module_name

            script = work_dir / "ecoli_html_reference.py"
            if not script.exists():
                self.banner.display_error("ecoli_html_reference.py not found")
                return False

            result = subprocess.run(
                [sys.executable, str(script)],
                cwd=work_dir, capture_output=True, text=True,
            )
            if not self._log_subprocess_result("Lineage", result):
                self.banner.display_warning(
                    f"Lineage generation had warnings (rc={result.returncode})"
                )

            html_src = work_dir / "ecoli_comprehensive_reference.html"
            if html_src.exists():
                dst_dir = output_dir / DIR_LINEAGE
                dst_dir.mkdir(exist_ok=True)
                shutil.copy2(html_src,
                             dst_dir / "ecoli_comprehensive_reference.html")
                self.banner.display_success("Lineage reference generated")
                return True
            return False

        except Exception as e:
            self.banner.display_error(f"Lineage error: {e}")
            self._log_error(traceback.format_exc())
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Final report + de-duplication
    # ------------------------------------------------------------------

    def copy_summary_results_to_final_directory(self, output_dir: Path) -> None:
        """Copy the top-level deliverables into EcoliTyper_final_report/.

        MGE standalone files are not copied — they live in mge_results/
        and are embedded inside the gene-centric / sample-centric reports.
        """
        try:
            self.banner.display_info(
                "Copying summary results to final directory..."
            )
            final_dir = output_dir / DIR_FINAL
            final_dir.mkdir(parents=True, exist_ok=True)

            for ext in [".html", ".json", ".tsv"]:
                src = output_dir / f"ecolityper_comprehensive_report{ext}"
                if src.exists():
                    shutil.copy2(src, final_dir / src.name)
                    self.banner.display_info(f"  * Copied: {src.name}")

            gc_dir = output_dir / DIR_GENE_CENTRIC
            if gc_dir.exists() and gc_dir.is_dir():
                target = final_dir / DIR_GENE_CENTRIC
                if target.exists():
                    shutil.rmtree(target)
                shutil.copytree(gc_dir, target)
                self.banner.display_info(
                    f"  * Copied: {DIR_GENE_CENTRIC}/ directory"
                )

            sc_dir = output_dir / DIR_SAMPLE_CENTRIC
            if sc_dir.exists() and sc_dir.is_dir():
                target = final_dir / DIR_SAMPLE_CENTRIC
                if target.exists():
                    shutil.rmtree(target)
                shutil.copytree(sc_dir, target)
                self.banner.display_info(
                    f"  * Copied: {DIR_SAMPLE_CENTRIC}/ directory"
                )

            self.banner.display_success(f"All results copied to: {final_dir}")
        except Exception as e:
            self.banner.display_error(f"Error copying summary results: {e}")

    def _cleanup_top_level_duplicates(self, output_path: Path) -> None:
        """Remove top-level copies that were duplicated into the final report."""
        removed_files = 0
        removed_dirs = 0

        for ext in [".html", ".json", ".tsv"]:
            f = output_path / f"ecolityper_comprehensive_report{ext}"
            if f.exists():
                f.unlink()
                removed_files += 1

        for dirname in [DIR_GENE_CENTRIC, DIR_SAMPLE_CENTRIC]:
            d = output_path / dirname
            if d.exists():
                shutil.rmtree(d, ignore_errors=True)
                removed_dirs += 1

        if removed_files or removed_dirs:
            self._log_info(
                f"Cleaned up {removed_files} duplicate file(s) and "
                f"{removed_dirs} duplicate dir(s) from top level"
            )

    def _finalize_output_by_report_level(self, output_path: Path) -> None:
        """Prune intermediate files based on report level.

        final   : Only EcoliTyper_final_report/ is kept.
        summary : EcoliTyper_final_report/ + batch files in module dirs.
        full    : No pruning.
        """
        if self.report_level == 'full':
            self.banner.display_info(
                "Report level: full (all intermediate files kept)"
            )
            return

        final_dir = output_path / DIR_FINAL
        final_exists = final_dir.exists() and any(final_dir.iterdir())
        if not final_exists:
            self.banner.display_warning(
                "Final report directory is missing or empty - skipping prune "
                "to avoid losing all output."
            )
            return

        if self.report_level == 'final':
            removed = 0
            for module_dir_name in ALL_MODULE_DIRS:
                module_dir = output_path / module_dir_name
                if not module_dir.exists():
                    continue
                shutil.rmtree(module_dir, ignore_errors=True)
                removed += 1
            self.banner.display_success(
                f"Report level: final - removed {removed} module directories"
            )
            return

        if self.report_level == 'summary':
            removed = 0
            for module_dir_name in ALL_MODULE_DIRS:
                module_dir = output_path / module_dir_name
                if not module_dir.exists():
                    continue
                for item in list(module_dir.iterdir()):
                    if item.is_dir():
                        shutil.rmtree(item, ignore_errors=True)
                        removed += 1
            if removed:
                self.banner.display_success(
                    f"Report level: summary - removed {removed} per-sample folders"
                )
            else:
                self.banner.display_success(
                    "Report level: summary - only batch summaries retained"
                )

    # ------------------------------------------------------------------
    # Sequential analyses
    # ------------------------------------------------------------------

    def run_sequential_analyses(self, fasta_files: List[Path], output_dir: Path,
                                threads: int, skip_modules: Dict[str, bool],
                                amr_min_identity: Optional[float],
                                amr_min_coverage: Optional[float],
                                amr_skip_mutations: bool,
                                amr_force_update: bool,
                                abricate_min_identity: int,
                                abricate_min_coverage: int) -> Dict[str, bool]:
        """Run every enabled analysis module in sequence.

        FASTA QC runs first so its species gate can filter which samples
        flow into the downstream modules.
        """
        results: Dict[str, bool] = {}
        downstream_files = fasta_files

        if not skip_modules.get('fasta_qc', False) and not self.interrupted:
            self.banner.display_module_header(
                "FASTA QC", "Sequence Quality Control & Species Confirmation"
            )
            results["FASTA QC"] = self.run_fasta_qc_analysis(
                fasta_files, output_dir, threads
            )
            if results["FASTA QC"]:
                self.banner.display_success("FASTA QC completed")
                gate = self._load_species_gate(output_dir)
                downstream_files = self._filter_by_species_gate(
                    fasta_files, gate
                )
                if not downstream_files:
                    self.banner.display_error(
                        "No samples passed the species gate - downstream "
                        "modules will be skipped."
                    )
                    for name in ["MLST", "Serotyping", "CH Typing",
                                 "Phylogrouping", "ABRicate",
                                 "AMRFinderPlus", "MGE"]:
                        results[name] = False
                    return results
            else:
                self.banner.display_warning(
                    "FASTA QC failed - species gate unavailable; "
                    "running all samples through downstream modules"
                )

        if not skip_modules.get('mlst', False) and not self.interrupted:
            self.banner.display_module_header(
                "MLST", "Multi-Locus Sequence Typing"
            )
            results["MLST"] = self.run_mlst_analysis(
                downstream_files, output_dir, threads
            )

        if not skip_modules.get('serotyping', False) and not self.interrupted:
            self.banner.display_module_header(
                "Serotyping", "O and H antigen determination"
            )
            results["Serotyping"] = self.run_serotyping_analysis(
                downstream_files, output_dir, threads
            )

        if not skip_modules.get('chtyper', False) and not self.interrupted:
            self.banner.display_module_header(
                "CH Typing", "FumC and FimH typing"
            )
            results["CH Typing"] = self.run_chtyper_analysis(
                downstream_files, output_dir, threads
            )

        if not skip_modules.get('phylogrouping', False) and not self.interrupted:
            self.banner.display_module_header(
                "Phylogrouping", "ezClermont phylogrouping algorithm"
            )
            results["Phylogrouping"] = self.run_phylogrouping_analysis(
                downstream_files, output_dir, threads
            )

        if not skip_modules.get('abricate', False) and not self.interrupted:
            self.banner.display_module_header(
                "ABRicate",
                "Resistance, Virulence, and Plasmid gene screening"
            )
            results["ABRicate"] = self.run_abricate_analysis(
                downstream_files, output_dir, threads,
                abricate_min_identity, abricate_min_coverage,
            )

        if not skip_modules.get('amrfinder', False) and not self.interrupted:
            self.banner.display_module_header(
                "AMRfinderPlus",
                "NCBI AMR gene detection with optional thresholds and "
                "mutation reporting",
            )
            results["AMRFinderPlus"] = self.run_amrfinder_analysis(
                downstream_files, output_dir, threads,
                amr_min_identity, amr_min_coverage,
                amr_skip_mutations, amr_force_update,
            )

        if not skip_modules.get('mge', False) and not self.interrupted:
            self.banner.display_module_header(
                "MGE Profiling",
                "Mobile Genetic Element Detection (mobileOG-db)"
            )
            results["MGE"] = self.run_mge_analysis(
                downstream_files, output_dir, threads
            )

        return results

    # ------------------------------------------------------------------
    # Full pipeline
    # ------------------------------------------------------------------

    def run_complete_analysis(self, input_path: str, output_dir: str,
                              threads: int = 1,
                              skip_modules: Optional[Dict[str, bool]] = None,
                              update_amr_db_only: bool = False,
                              force_update_amr_db: bool = False,
                              amr_min_identity: Optional[float] = None,
                              amr_min_coverage: Optional[float] = None,
                              amr_skip_mutations: bool = False,
                              amr_force_update: bool = False,
                              abricate_min_identity: int = 80,
                              abricate_min_coverage: int = 80,
                              keep_temp: bool = False,
                              report_level: str = 'summary'):
        """Execute the full pipeline end-to-end.

        Handles database-only modes first, then runs every enabled analysis,
        both reporters, the lineage reference, and the visualization module
        before consolidating and pruning the output directory.
        """
        if skip_modules is None:
            skip_modules = {}
        self.keep_temp = keep_temp
        self.report_level = report_level

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
            self.banner.display_success(
                f"Starting analysis of {len(self.fasta_files)} E. coli genomes"
            )
            self.banner.display_info(
                f"File formats detected: {', '.join(extensions)}"
            )

            # Warn early if MGE is enabled but the DB is missing
            if (not skip_modules.get('mge', False)
                    and self.get_mge_data_dir() is None):
                self.banner.display_warning(
                    "MGE profiling is enabled but mobileOG-db is not installed. "
                    "Run `ecolityper --pull-mge-db` first, or add --skip-mge."
                )

            self.banner.display_module_header(
                "Analysis Plan", "Modules to be executed"
            )
            plan = [
                ("FASTA QC", not skip_modules.get('fasta_qc', False)),
                ("MLST", not skip_modules.get('mlst', False)),
                ("Serotyping", not skip_modules.get('serotyping', False)),
                ("CH Typing", not skip_modules.get('chtyper', False)),
                ("Phylogrouping", not skip_modules.get('phylogrouping', False)),
                ("ABRicate", not skip_modules.get('abricate', False)),
                ("AMRfinderPlus", not skip_modules.get('amrfinder', False)),
                ("MGE Profiling", not skip_modules.get('mge', False)),
                ("Comprehensive + Gene-Centric Reports",
                 not skip_modules.get('summary', False)),
                ("Sample-Centric Reporter",
                 not skip_modules.get('samplecentric', False)),
                ("Lineage Reference", not skip_modules.get('lineage', False)),
                ("Visualizations", not skip_modules.get('visualization', False)),
            ]
            for analysis, enabled in plan:
                status = "✅ ENABLED" if enabled else "⏸️  SKIPPED"
                print(f"   {status} - {analysis}")
            print(f"\n   Report level: {self.report_level}")
            sys.stdout.flush()

            analysis_results = self.run_sequential_analyses(
                self.fasta_files, output_path, threads, skip_modules,
                amr_min_identity, amr_min_coverage,
                amr_skip_mutations, amr_force_update,
                abricate_min_identity, abricate_min_coverage,
            )

            if not skip_modules.get('summary', False) and not self.interrupted:
                self.banner.display_module_header(
                    "Comprehensive + Gene-Centric Reports",
                    "Typing integration and gene-centric interactive report",
                )
                analysis_results["Comprehensive + Gene-Centric Reports"] = \
                    self.run_comprehensive_and_gene_centric(output_path)

            if not skip_modules.get('samplecentric', False) and not self.interrupted:
                self.banner.display_module_header(
                    "Sample-Centric Reporter",
                    "Interactive isolate boxes (AMR, Virulence, Plasmids, "
                    "BACMET, Mutations, MGE)",
                )
                analysis_results["Sample-Centric Reporter"] = \
                    self.run_sample_centric_analysis(output_path)

            if not skip_modules.get('lineage', False) and not self.interrupted:
                self.banner.display_module_header(
                    "Lineage Reference", "E. coli lineage reference database"
                )
                analysis_results["Lineage Reference"] = \
                    self.run_lineage_analysis(output_path)

            if not skip_modules.get('visualization', False) and not self.interrupted:
                self.banner.display_module_header(
                    "Visualizations", "Charts and interactive dashboards"
                )
                analysis_results["Visualizations"] = \
                    self.run_visualization_analysis(output_path)

            self.copy_summary_results_to_final_directory(output_path)
            self._cleanup_top_level_duplicates(output_path)
            self._finalize_output_by_report_level(output_path)

            self.total_duration = time.time() - self.start_time
            self.banner.display_footer(
                analysis_time=self._format_duration(self.total_duration),
                samples_processed=len(self.fasta_files),
            )

            successful = sum(analysis_results.values())
            total = len(analysis_results)
            if successful == total:
                self.banner.display_success(
                    f"All {total} analyses completed successfully!"
                )
            else:
                self.banner.display_warning(
                    f"{successful}/{total} analyses completed successfully"
                )

            print("\n📁 FINAL OUTPUT STRUCTURE:")
            for item in sorted(output_path.glob("*")):
                if item.is_dir():
                    if item.name == DIR_FINAL:
                        print(f"   📂 {item.name}/")
                        for sub_item in sorted(item.glob("*")):
                            if sub_item.is_file():
                                print(f"      📄 {sub_item.name}")
                            elif sub_item.is_dir():
                                sub_count = len(list(sub_item.rglob("*")))
                                print(
                                    f"      📂 {sub_item.name}/ "
                                    f"({sub_count} items)"
                                )
                    else:
                        count = len(list(item.rglob("*")))
                        print(f"   📂 {item.name}/ ({count} items)")
                else:
                    print(f"   📄 {item.name}")

            print("\n📚 Please cite our EcoliTyper paper:")
            for line in CITATION_LINES:
                print(f"   {line}")

        except KeyboardInterrupt:
            self.banner.display_error("Analysis interrupted by user")
        except Exception as e:
            self.banner.display_error(f"Critical error: {str(e)}")
            traceback.print_exc()
            sys.exit(1)


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            f"{Colors.BRIGHT_YELLOW}EcoliTyper: Complete E. coli Typing Pipeline"
            f"{Colors.RESET}"
        ),
        formatter_class=ColoredHelpFormatter,
        epilog=f"""
{Colors.BOLD}{Colors.BRIGHT_GREEN}EXAMPLES:{Colors.RESET}

  {Colors.BRIGHT_CYAN}ecolityper -i genome.fna -o results/{Colors.RESET}
  {Colors.BRIGHT_CYAN}ecolityper -i "*.fna" -o batch_results --threads 8{Colors.RESET}
  {Colors.BRIGHT_CYAN}ecolityper -i "*.fasta" -o analysis --threads 16 --skip-visualization{Colors.RESET}
  {Colors.BRIGHT_CYAN}ecolityper -i "*.fna" -o results --report-level final{Colors.RESET}   # minimal disk usage

  {Colors.BRIGHT_CYAN}ecolityper --update-amr-db{Colors.RESET}                # Update AMR database (incremental)
  {Colors.BRIGHT_CYAN}ecolityper --force-update-amr-db{Colors.RESET}          # Force full AMR database update

  {Colors.BRIGHT_CYAN}ecolityper -i "*.fna" -o results --abricate-minid 90 --abricate-mincov 85{Colors.RESET}
  {Colors.BRIGHT_CYAN}ecolityper -i "*.fna" -o results --amr-min-identity 0.95 --amr-min-coverage 0.9 --skip-amr-mutations{Colors.RESET}

{Colors.BOLD}{Colors.BRIGHT_GREEN}MLST DATABASE MANAGEMENT:{Colors.RESET}
  {Colors.BRIGHT_CYAN}ecolityper --pull-mlst-db{Colors.RESET}
     Download a pre-built E. coli MLST database from GitHub.
     Recommended for first-time users (no API key needed).

  {Colors.BRIGHT_CYAN}ecolityper --update-mlst-db{Colors.RESET}
     Update the E. coli MLST database directly from PubMLST.
     Requires one-time API key setup (see instructions below).
     Download may take 10-20 minutes depending on internet speed.

{Colors.BOLD}{Colors.BRIGHT_GREEN}MGE DATABASE MANAGEMENT:{Colors.RESET}
  {Colors.BRIGHT_CYAN}ecolityper --pull-mge-db{Colors.RESET}
     Download the mobileOG-db payload (~300 MB) from GitHub.
     Required once after conda install if you want MGE profiling.
     Docker and Singularity users do NOT need this —
     the data is already baked into the image.

     Location: ~/.local/share/ecolityper/mge_module/mobileOG-db/
     Env override: ECOLIITYPER_MGE_DATA_DIR
     Override URL: ECOLIITYPER_MGE_PULL_URL

{Colors.BOLD}{Colors.BRIGHT_GREEN}PubMLST API Key Setup (for --update-mlst-db):{Colors.RESET}
  1. Register/login at https://pubmlst.org/site-accounts
  2. Generate an API key from your profile (under "API keys").
  3. Run: mlstdb connect --db pubmlst --api-key
  4. Paste your API key when prompted.
  5. After setup, run ecolityper --update-mlst-db whenever you want to refresh the database.

{Colors.BOLD}{Colors.BRIGHT_GREEN}HPC / Docker Environment Variables:{Colors.RESET}
  All overrides are optional. Useful when the module directory is read-only
  (Docker images, system-wide conda, shared HPC installs) or when databases
  live on a shared filesystem (Lustre, NFS, BeeGFS).

  ECOLIITYPER_SHARED_DATA         Writable root for user-local data
                                   (default: ~/.local/share/ecolityper/)

  ECOLIITYPER_MLST_DB             MLST database directory
  ECOLIITYPER_MLST_MODULE_PATH    MLST module root (bin/, db/, perl5/, scripts/)
  ECOLIITYPER_MLST_PULL_URL       Override default GitHub pull URL

  ECOLIITYPER_MGE_DATA_DIR        mobileOG-db data directory
  ECOLIITYPER_MGE_MODULE_PATH     MGE module root
  ECOLIITYPER_MGE_PULL_URL        Override default GitHub pull URL

  Example (HPC with shared DBs on Lustre):
    export ECOLIITYPER_MLST_DB=/lustre/scratch/dbs/mlst/pubmlst
    export ECOLIITYPER_MGE_DATA_DIR=/lustre/scratch/dbs/mge/mobileOG-db
    ecolityper -i "*.fna" -o results --report-level final

{Colors.BOLD}{Colors.BRIGHT_GREEN}REPORT LEVELS:{Colors.RESET}
  {Colors.YELLOW}final{Colors.RESET}   - Only {DIR_FINAL}/ is kept. All module dirs are deleted.
            Best for archival and very large batches.
  {Colors.YELLOW}summary{Colors.RESET} - (default) {DIR_FINAL}/ + batch TSVs/HTMLs/JSONs in each
            module directory. Per-sample folders are removed.
  {Colors.YELLOW}full{Colors.RESET}    - Everything is kept, including all per-sample HTMLs.
            Best for small batches or when debugging.

{Colors.BOLD}{Colors.BRIGHT_GREEN}SUPPORTED FASTA FORMATS:{Colors.RESET} .fna, .fasta, .fa, .fsa

{Colors.BOLD}{Colors.BRIGHT_GREEN}ANALYSIS MODULES:{Colors.RESET}
  * FASTA QC (Sequence Quality Control + Species Confirmation)
  * MLST (Multi-Locus Sequence Typing)
  * Serotyping (O and H antigen determination)
  * CH Typing (FumC and FimH typing)
  * Phylogrouping (ezClermont algorithm)
  * ABRicate (Resistance/Virulence/Plasmid screening)
  * AMRfinderPlus (NCBI AMR gene detection + mutations)
  * MGE Profiling (mobileOG-db)
  * Comprehensive + Gene-Centric Reports
  * Sample-Centric Reporter
  * Lineage reference
  * Visualizations

{Colors.BOLD}{Colors.BRIGHT_GREEN}PRIOR SETUP REQUIRED:{Colors.RESET}
  {Colors.YELLOW}* abricate --setupdb{Colors.RESET}
  {Colors.YELLOW}* ecolityper --update-amr-db{Colors.RESET}
  {Colors.YELLOW}* ecolityper --pull-mge-db{Colors.RESET}   (conda users only)

{Colors.BOLD}{Colors.BRIGHT_GREEN}CONTACT:{Colors.RESET}
  Brown Beckley <brownbeckley94@gmail.com>
  University of Ghana Medical School - Department of Medical Biochemistry

  {Colors.BRIGHT_YELLOW}Drop a Star on GitHub if you find this tool useful!{Colors.RESET}
"""
    )
    parser.add_argument('-i', '--input',
                        help='Input FASTA file(s) - can use glob patterns like "*.fna"')
    parser.add_argument('-o', '--output',
                        help='Output directory for all results')
    parser.add_argument('-t', '--threads', type=int, default=2,
                        help='Number of threads (default: 2)')
    parser.add_argument('--keep-temp', action='store_true',
                        help='Do not delete temporary directories (for debugging)')

    parser.add_argument('--update-amr-db', action='store_true',
                        help='Update AMRfinderPlus database (incremental) and exit')
    parser.add_argument('--force-update-amr-db', action='store_true',
                        help='Force complete AMR database update and exit')

    parser.add_argument('--pull-mlst-db', action='store_true',
                        help='Pull E. coli MLST database from GitHub '
                             '(recommended for first-time setup)')
    parser.add_argument('--update-mlst-db', action='store_true',
                        help='Update E. coli MLST database from PubMLST '
                             '(API key required)')

    parser.add_argument('--pull-mge-db', action='store_true',
                        help='Pull the mobileOG-db payload (MGE profiling data) '
                             'from GitHub (~300 MB, no API key)')

    parser.add_argument('--amr-min-identity', type=float,
                        help='Minimum identity for AMR hits (0..1)')
    parser.add_argument('--amr-min-coverage', type=float,
                        help='Minimum coverage for AMR hits (0..1)')
    parser.add_argument('--skip-amr-mutations', action='store_true',
                        help='Disable point mutation reporting in AMR '
                             '(enabled by default)')
    parser.add_argument('--amr-force-update', action='store_true',
                        help='Force update AMR database before analysis')

    parser.add_argument('--abricate-minid', type=int, default=80,
                        help='Minimum identity for ABRicate hits (default: 80)')
    parser.add_argument('--abricate-mincov', type=int, default=80,
                        help='Minimum coverage for ABRicate hits (default: 80)')

    parser.add_argument('--report-level',
                        choices=['final', 'summary', 'full'],
                        default='summary',
                        help='Controls which files are kept after the run. '
                             'final = only EcoliTyper_final_report/ ; '
                             'summary = (default) final report + batch files, '
                             'per-sample folders removed ; '
                             'full = everything including per-sample folders.')

    skip_group = parser.add_argument_group(
        f'{Colors.BOLD}{Colors.BRIGHT_MAGENTA}'
        f'Skip Options (disable specific analyses)'
        f'{Colors.RESET}'
    )
    skip_group.add_argument('--skip-fasta-qc', action='store_true',
                            help='Skip FASTA QC analysis')
    skip_group.add_argument('--skip-amrfinder', action='store_true',
                            help='Skip AMRfinderPlus analysis')
    skip_group.add_argument('--skip-abricate', action='store_true',
                            help='Skip ABRicate analysis')
    skip_group.add_argument('--skip-mlst', action='store_true',
                            help='Skip MLST analysis')
    skip_group.add_argument('--skip-serotyping', action='store_true',
                            help='Skip serotyping analysis')
    skip_group.add_argument('--skip-chtyper', action='store_true',
                            help='Skip CH typing analysis')
    skip_group.add_argument('--skip-phylogrouping', action='store_true',
                            help='Skip phylogrouping analysis')
    skip_group.add_argument('--skip-mge', action='store_true',
                            help='Skip MGE profiling (mobileOG)')
    skip_group.add_argument('--skip-lineage', action='store_true',
                            help='Skip lineage reference generation')
    skip_group.add_argument('--skip-summary', action='store_true',
                            help='Skip comprehensive + gene-centric report generation')
    skip_group.add_argument('--skip-samplecentric', action='store_true',
                            help='Skip sample-centric reporter generation')
    skip_group.add_argument('--skip-visualization', action='store_true',
                            help='Skip visualization generation')

    args = parser.parse_args()

    if args.abricate_minid < 0 or args.abricate_minid > 100:
        parser.error("--abricate-minid must be between 0 and 100")
    if args.abricate_mincov < 0 or args.abricate_mincov > 100:
        parser.error("--abricate-mincov must be between 0 and 100")

    # Database-management-only paths exit early
    if args.update_amr_db or args.force_update_amr_db:
        orch = EcoliTyperOrchestrator()
        if args.force_update_amr_db:
            orch.update_amr_database(force=True)
        else:
            orch.update_amr_database(force=False)
        sys.exit(0)

    if args.pull_mlst_db or args.update_mlst_db:
        orch = EcoliTyperOrchestrator()
        temp_log_dir = Path(tempfile.mkdtemp(prefix="ecolityper_mlst_db_"))
        orch.setup_logging(temp_log_dir)
        if args.pull_mlst_db:
            success = orch.pull_mlst_database()
        else:
            success = orch.update_mlst_database()
        if not orch.keep_temp:
            shutil.rmtree(temp_log_dir, ignore_errors=True)
        sys.exit(0 if success else 1)

    if args.pull_mge_db:
        orch = EcoliTyperOrchestrator()
        temp_log_dir = Path(tempfile.mkdtemp(prefix="ecolityper_mge_db_"))
        orch.setup_logging(temp_log_dir)
        success = orch.pull_mge_database()
        if not orch.keep_temp:
            shutil.rmtree(temp_log_dir, ignore_errors=True)
        sys.exit(0 if success else 1)

    if not args.input or not args.output:
        parser.error(
            "Both -i/--input and -o/--output are required for analysis "
            "(or use --update-amr-db / --force-update-amr-db / "
            "--pull-mlst-db / --update-mlst-db / --pull-mge-db)"
        )

    if (args.report_level == 'final'
            and args.skip_summary
            and args.skip_samplecentric):
        parser.error(
            "--report-level final requires at least one reporter to run "
            "(--skip-summary and --skip-samplecentric both set would leave "
            "the final report empty and pruning would delete all output)."
        )

    skip_modules = {
        'fasta_qc': args.skip_fasta_qc,
        'amrfinder': args.skip_amrfinder,
        'abricate': args.skip_abricate,
        'mlst': args.skip_mlst,
        'serotyping': args.skip_serotyping,
        'chtyper': args.skip_chtyper,
        'phylogrouping': args.skip_phylogrouping,
        'mge': args.skip_mge,
        'lineage': args.skip_lineage,
        'summary': args.skip_summary,
        'samplecentric': args.skip_samplecentric,
        'visualization': args.skip_visualization,
    }

    orchestrator = EcoliTyperOrchestrator()
    orchestrator.run_complete_analysis(
        input_path=args.input,
        output_dir=args.output,
        threads=args.threads,
        skip_modules=skip_modules,
        amr_min_identity=args.amr_min_identity,
        amr_min_coverage=args.amr_min_coverage,
        amr_skip_mutations=args.skip_amr_mutations,
        amr_force_update=args.amr_force_update,
        abricate_min_identity=args.abricate_minid,
        abricate_min_coverage=args.abricate_mincov,
        keep_temp=args.keep_temp,
        report_level=args.report_level,
    )


if __name__ == "__main__":
    main()