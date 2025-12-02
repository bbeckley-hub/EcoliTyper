#!/usr/bin/env python3
"""
EcoliTyper Main Orchestrator - Complete E. coli Typing Pipeline
Comprehensive E. coli analysis: MLST, Serotyping, CH Typing, Phylogrouping, Abricate, AMRfinderPlus
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
Date: 2025
Send a quick mail for any issues or further explanations.
"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
import signal
import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from typing import Dict, List, Set

# Import banner
try:
    from .core.banner import EcoliTyperBanner
except (ImportError, SystemError):
    # Fallback import
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from core.banner import EcoliTyperBanner

class EcoliTyperOrchestrator:
    """EcoliTyper orchestrator with comprehensive cleanup and interrupt handling"""
    
    def __init__(self):
        self.banner = EcoliTyperBanner()
        self.base_dir = Path(__file__).parent
        self.fasta_files = []
        self.interrupted = False
        self.output_lock = threading.Lock()  # Prevent output mixing
        
        # Setup interrupt handler for automatic cleanup
        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)
    
    def _signal_handler(self, signum, frame):
        """Handle interrupt signals with automatic cleanup"""
        self.interrupted = True
        self.banner.display_error(f"Analysis interrupted by user (signal {signum})")
        self.banner.display_info("Starting automatic cleanup...")
        self._emergency_cleanup()
        sys.exit(1)
    
    def _emergency_cleanup(self):
        """Emergency cleanup when interrupted"""
        try:
            # Clean all module directories
            modules = [
                "mlst_module", "serotypefinder_module", 
                "CHTyper_module", "phylogrouping_module",
                "Abricate_module", "Amrfinder_module"
            ]
            
            for module in modules:
                module_path = self.base_dir / "modules" / module
                if module_path.exists():
                    self.cleanup_module_directory(module_path, self.fasta_files)
            
            self.banner.display_success("Emergency cleanup completed!")
        except Exception as e:
            self.banner.display_error(f"Emergency cleanup failed: {str(e)}")
    
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files using glob patterns - SAME LOGIC AS YOUR MODULES"""
        self.banner.display_info(f"Searching for files with pattern: {input_path}")
        
        # Handle quoted wildcards properly
        if '*' in input_path or '?' in input_path:
            matched_files = glob.glob(input_path)
            fasta_files = [Path(f) for f in matched_files if Path(f).is_file() and 
                          f.lower().endswith(('.fna', '.fasta', '.fa', '.fsa')) and
                          not Path(f).name.startswith('.')]
            self.banner.display_success(f"Found {len(fasta_files)} FASTA files")
            return sorted(fasta_files)
        
        # Handle direct file path
        input_path_obj = Path(input_path)
        if input_path_obj.is_file() and input_path_obj.suffix.lower() in ['.fna', '.fasta', '.fa', '.fsa']:
            self.banner.display_success(f"Found single FASTA file: {input_path_obj.name}")
            return [input_path_obj]
        
        # Handle directory
        if input_path_obj.is_dir():
            patterns = [
                f"{input_path}/*.fna", f"{input_path}/*.fasta",
                f"{input_path}/*.fa", f"{input_path}/*.fsa"
            ]
            fasta_files = []
            for pattern in patterns:
                matched_files = glob.glob(pattern)
                for file_path in matched_files:
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
        """Get the correct file pattern based on actual file extensions"""
        if not fasta_files:
            return "*.fna"
        
        # Get all unique extensions from the input files
        extensions = set(f.suffix.lower() for f in fasta_files)
        
        # If all files have the same extension, use that
        if len(extensions) == 1:
            ext = list(extensions)[0]
            return f"*{ext}"
        
        # If mixed extensions, use a pattern that matches all FASTA files
        return "*"

    def cleanup_module_directory(self, module_path: Path, fasta_files: List[Path]):
        """COMPREHENSIVE cleanup of module directory after analysis"""
        try:
            with self.output_lock:
                self.banner.display_info(f"Cleaning up {module_path.name}...")
            
            # 1. Remove copied input files
            for fasta_file in fasta_files:
                temp_file = module_path / fasta_file.name
                if temp_file.exists():
                    temp_file.unlink()
            
            # 2. Remove common output directories
            output_dirs = [
                "mlst_results", "results", "SerotypeFinder_results",
                "chtyper_results", "phylogrouping_results",
                "ecoli_abricate_results", "ecoli_amrfinder_results"
            ]
            for output_dir in output_dirs:
                dir_path = module_path / output_dir
                if dir_path.exists():
                    shutil.rmtree(dir_path)
            
            # 3. Remove any other common temporary files
            temp_patterns = ["*.txt", "*.log", "*.tmp", "temp_*", "*.html", "*.tsv"]
            for pattern in temp_patterns:
                for temp_file in module_path.glob(pattern):
                    if temp_file.is_file():
                        temp_file.unlink()
            
            with self.output_lock:
                self.banner.display_success(f"✅ {module_path.name} cleaned up successfully")
            
        except Exception as e:
            with self.output_lock:
                self.banner.display_warning(f"⚠️  Partial cleanup issue in {module_path.name}: {str(e)}")

    def run_mlst_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run MLST analysis - PRODUCTION VERSION"""
        mlst_module_path = self.base_dir / "modules" / "mlst_module"
        
        try:
            with self.output_lock:
                self.banner.start_analysis_timer("mlst")
                self.banner.display_module_header("MLST Analysis", "Multi-Locus Sequence Typing for E. coli")
            
            mlst_script = mlst_module_path / "ecolimlst_module.py"
            
            if not mlst_script.exists():
                with self.output_lock:
                    self.banner.display_error(f"MLST script not found at: {mlst_script}")
                return False
            
            # Clean any existing results first
            results_dir = mlst_module_path / "results"
            if results_dir.exists():
                shutil.rmtree(results_dir)
            
            # Copy files to MLST module directory
            with self.output_lock:
                self.banner.display_info(f"Copied {len(fasta_files)} files to MLST module")
            
            for fasta_file in fasta_files:
                target_file = mlst_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            # Build the correct command
            if len(fasta_files) == 1:
                # Single file - use direct filename
                fasta_file = fasta_files[0]
                cmd = [
                    sys.executable, str(mlst_script),
                    "-i", fasta_file.name,
                    "-o", "results",
                    "-db", "db",
                    "-sc", "bin",
                    "--batch"
                ]
                with self.output_lock:
                    self.banner.display_info(f"Running MLST analysis on: {fasta_file.name}")
            else:
                # Multiple files - use pattern
                file_pattern = self.get_file_pattern(fasta_files)
                cmd = [
                    sys.executable, str(mlst_script),
                    "-i", file_pattern,
                    "-o", "results", 
                    "-db", "db",
                    "-sc", "bin",
                    "--batch"
                ]
                with self.output_lock:
                    self.banner.display_info(f"Running MLST analysis with pattern: {file_pattern}")
        
            # Run the analysis
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=mlst_module_path)
            
            # Check if results were actually generated
            mlst_success = False
            if results_dir.exists():
                # Check for summary files to verify success
                summary_files = list(results_dir.glob("*summary*"))
                sample_dirs = [d for d in results_dir.iterdir() if d.is_dir()]
                
                if summary_files or sample_dirs:
                    mlst_success = True
                    with self.output_lock:
                        self.banner.stop_analysis_timer("mlst")
                        self.banner.display_success("MLST analysis completed!")
                    
                    # Copy results to output directory
                    mlst_target = output_dir / "mlst_results"
                    if mlst_target.exists():
                        shutil.rmtree(mlst_target)
                    shutil.copytree(results_dir, mlst_target)
                    
                    with self.output_lock:
                        self.banner.display_success(f"MLST results copied to: {mlst_target}")
                    
                    # Verify we have actual ST results, not just UNKNOWN
                    summary_tsv = mlst_target / "mlst_summary.tsv"
                    if summary_tsv.exists():
                        with open(summary_tsv, 'r') as f:
                            content = f.read()
                            if "UNKNOWN" in content or "ND" in content:
                                with self.output_lock:
                                    self.banner.display_warning("MLST analysis completed but some samples have unknown ST")
                else:
                    with self.output_lock:
                        self.banner.display_warning("MLST analysis ran but produced no result files")
            else:
                with self.output_lock:
                    self.banner.display_warning("MLST analysis produced no output directory")
            
            return mlst_success
                
        except Exception as e:
            with self.output_lock:
                self.banner.display_error(f"MLST analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup
            self.cleanup_module_directory(mlst_module_path, fasta_files)

    def run_serotyping_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run serotyping analysis - WITH CLEANUP"""
        sero_module_path = self.base_dir / "modules" / "serotypefinder_module"
        
        try:
            with self.output_lock:
                self.banner.start_analysis_timer("serotyping")
                self.banner.display_module_header("Serotyping Analysis", "O and H antigen determination")
            
            sero_script = sero_module_path / "enhanced_serotypefinder.py"
            
            if not sero_script.exists():
                with self.output_lock:
                    self.banner.display_error(f"Serotyping script not found at: {sero_script}")
                return False
            
            # Copy files to serotyping module directory
            for fasta_file in fasta_files:
                target_file = sero_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            with self.output_lock:
                self.banner.display_info(f"Copied {len(fasta_files)} files to serotyping module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            
            # Build command - use direct command list without shell=True
            cmd = [
                sys.executable, str(sero_script),
                "-i", file_pattern,
                "-o", "Serotype"
            ]
            
            with self.output_lock:
                self.banner.display_info(f"Running serotyping analysis with pattern: {file_pattern}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=sero_module_path)
            
            if result.returncode == 0:
                with self.output_lock:
                    self.banner.stop_analysis_timer("serotyping")
                    self.banner.display_success("Serotyping analysis completed!")
                
                # Copy results to output directory
                sero_source = sero_module_path / "Serotype" / "SerotypeFinder_results"
                sero_target = output_dir / "serotyping_results"
                
                if sero_source.exists():
                    if sero_target.exists():
                        shutil.rmtree(sero_target)
                    shutil.copytree(sero_source, sero_target)
                    
                    with self.output_lock:
                        self.banner.display_success(f"Serotyping results copied to: {sero_target}")
                
                return True
            else:
                with self.output_lock:
                    self.banner.display_warning("Serotyping analysis had warnings")
                    if result.stderr:
                        self.banner.display_info(f"Serotyping stderr: {result.stderr[:200]}...")
                return True
                
        except Exception as e:
            with self.output_lock:
                self.banner.display_error(f"Serotyping analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.cleanup_module_directory(sero_module_path, fasta_files)

    def run_chtyper_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run CH typing analysis - WITH CLEANUP"""
        chtyper_module_path = self.base_dir / "modules" / "CHTyper_module"
        
        try:
            with self.output_lock:
                self.banner.start_analysis_timer("ch_typing")
                self.banner.display_module_header("CH Typing Analysis", "CH (FumC and FimH) typing analysis")
            
            chtyper_script = chtyper_module_path / "enhanced_chtyper.py"
            
            if not chtyper_script.exists():
                with self.output_lock:
                    self.banner.display_error(f"CHTyper script not found at: {chtyper_script}")
                return False
            
            # Copy files to CHTyper module directory
            for fasta_file in fasta_files:
                target_file = chtyper_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            with self.output_lock:
                self.banner.display_info(f"Copied {len(fasta_files)} files to CHTyper module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            
            # Build command - use direct command list
            cmd = [
                sys.executable, str(chtyper_script),
                "-i", file_pattern,
                "-o", "CH_results"
            ]
            
            with self.output_lock:
                self.banner.display_info(f"Running CH typing analysis with pattern: {file_pattern}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=chtyper_module_path)
            
            if result.returncode == 0:
                with self.output_lock:
                    self.banner.stop_analysis_timer("ch_typing")
                    self.banner.display_success("CH typing analysis completed!")
                
                # Copy results to output directory
                chtyper_source = chtyper_module_path / "CH_results" / "chtyper_results"
                chtyper_target = output_dir / "chtyper_results"
                
                if chtyper_source.exists():
                    if chtyper_target.exists():
                        shutil.rmtree(chtyper_target)
                    shutil.copytree(chtyper_source, chtyper_target)
                    
                    with self.output_lock:
                        self.banner.display_success(f"CH typing results copied to: {chtyper_target}")
                
                return True
            else:
                with self.output_lock:
                    self.banner.display_warning("CH typing analysis had warnings")
                    if result.stderr:
                        self.banner.display_info(f"CH typing stderr: {result.stderr[:200]}...")
                return True
                
        except Exception as e:
            with self.output_lock:
                self.banner.display_error(f"CH typing analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.cleanup_module_directory(chtyper_module_path, fasta_files)

    def run_phylogrouping_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run phylogrouping analysis - WITH CLEANUP"""
        phylo_module_path = self.base_dir / "modules" / "phylogrouping_module"
        
        try:
            with self.output_lock:
                self.banner.start_analysis_timer("phylogrouping")
                self.banner.display_module_header("Phylogrouping Analysis", "zClermont phylogrouping algorithm")
            
            phylo_script = phylo_module_path / "enhanced_ezclermont.py"
            
            if not phylo_script.exists():
                with self.output_lock:
                    self.banner.display_error(f"Phylogrouping script not found at: {phylo_script}")
                return False
            
            # Copy files to phylogrouping module directory
            for fasta_file in fasta_files:
                target_file = phylo_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            with self.output_lock:
                self.banner.display_info(f"Copied {len(fasta_files)} files to phylogrouping module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            
            # Build command - use direct command list
            cmd = [
                sys.executable, str(phylo_script),
                "-i", file_pattern,
                "-o", "Phylo"
            ]
            
            with self.output_lock:
                self.banner.display_info(f"Running phylogrouping analysis with pattern: {file_pattern}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=phylo_module_path)
            
            if result.returncode == 0:
                with self.output_lock:
                    self.banner.stop_analysis_timer("phylogrouping")
                    self.banner.display_success("Phylogrouping analysis completed!")
                
                # Copy results to output directory
                phylo_source = phylo_module_path / "Phylo" / "phylogrouping_results"
                phylo_target = output_dir / "phylogrouping_results"
                
                if phylo_source.exists():
                    if phylo_target.exists():
                        shutil.rmtree(phylo_target)
                    shutil.copytree(phylo_source, phylo_target)
                    
                    with self.output_lock:
                        self.banner.display_success(f"Phylogrouping results copied to: {phylo_target}")
                
                return True
            else:
                with self.output_lock:
                    self.banner.display_warning("Phylogrouping analysis had warnings")
                    if result.stderr:
                        self.banner.display_info(f"Phylogrouping stderr: {result.stderr[:200]}...")
                return True
                
        except Exception as e:
            with self.output_lock:
                self.banner.display_error(f"Phylogrouping analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.cleanup_module_directory(phylo_module_path, fasta_files)

    def run_abricate_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run Abricate analysis - WITH CLEANUP"""
        abricate_module_path = self.base_dir / "modules" / "Abricate_module"
        
        try:
            with self.output_lock:
                self.banner.start_analysis_timer("abricate")
                self.banner.display_module_header("ABRicate Analysis", "Resistance, Virulence, and Plasmid gene screening")
            
            abricate_script = abricate_module_path / "ecoli_abricate.py"
            
            if not abricate_script.exists():
                with self.output_lock:
                    self.banner.display_error(f"ABRicate script not found at: {abricate_script}")
                return False
            
            # Copy files to Abricate module directory
            for fasta_file in fasta_files:
                target_file = abricate_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            with self.output_lock:
                self.banner.display_info(f"Copied {len(fasta_files)} files to ABRicate module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            
            # Build command - use direct command list
            cmd = [
                sys.executable, str(abricate_script),
                file_pattern
            ]
            
            with self.output_lock:
                self.banner.display_info(f"Running ABRicate analysis with pattern: {file_pattern}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=abricate_module_path)
            
            if result.returncode == 0:
                with self.output_lock:
                    self.banner.stop_analysis_timer("abricate")
                    self.banner.display_success("ABRicate analysis completed!")
                
                # Copy results to output directory
                abricate_source = abricate_module_path / "ecoli_abricate_results"
                abricate_target = output_dir / "abricate_results"
                
                if abricate_source.exists():
                    if abricate_target.exists():
                        shutil.rmtree(abricate_target)
                    shutil.copytree(abricate_source, abricate_target)
                    
                    with self.output_lock:
                        self.banner.display_success(f"ABRicate results copied to: {abricate_target}")
                
                return True
            else:
                with self.output_lock:
                    self.banner.display_warning("ABRicate analysis had warnings")
                    if result.stderr:
                        self.banner.display_info(f"ABRicate stderr: {result.stderr[:200]}...")
                return True
                
        except Exception as e:
            with self.output_lock:
                self.banner.display_error(f"ABRicate analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.cleanup_module_directory(abricate_module_path, fasta_files)

    def run_amrfinder_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run AMRfinderPlus analysis - WITH CLEANUP (ALWAYS LAST)"""
        amr_module_path = self.base_dir / "modules" / "Amrfinder_module"
        
        try:
            with self.output_lock:
                self.banner.start_analysis_timer("amrfinder")
                self.banner.display_module_header("AMRfinderPlus Analysis", "NCBI AMR gene detection")
            
            amr_script = amr_module_path / "ecoli_amrfinder.py"
            
            if not amr_script.exists():
                with self.output_lock:
                    self.banner.display_error(f"AMRfinderPlus script not found at: {amr_script}")
                return False
            
            # Copy files to AMR module directory
            for fasta_file in fasta_files:
                target_file = amr_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            with self.output_lock:
                self.banner.display_info(f"Copied {len(fasta_files)} files to AMRfinderPlus module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            
            # Build command - use direct command list
            cmd = [
                sys.executable, str(amr_script),
                file_pattern
            ]
            
            with self.output_lock:
                self.banner.display_info(f"Running AMRfinderPlus analysis with pattern: {file_pattern}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
            
            if result.returncode == 0:
                with self.output_lock:
                    self.banner.stop_analysis_timer("amrfinder")
                    self.banner.display_success("AMRfinderPlus analysis completed!")
                
                # Copy results to output directory
                amr_source = amr_module_path / "ecoli_amrfinder_results"
                amr_target = output_dir / "amrfinder_results"
                
                if amr_source.exists():
                    if amr_target.exists():
                        shutil.rmtree(amr_target)
                    shutil.copytree(amr_source, amr_target)
                    
                    with self.output_lock:
                        self.banner.display_success(f"AMRfinderPlus results copied to: {amr_target}")
                
                return True
            else:
                with self.output_lock:
                    self.banner.display_warning("AMRfinderPlus analysis had warnings")
                    if result.stderr:
                        self.banner.display_info(f"AMRfinderPlus stderr: {result.stderr[:200]}...")
                return True
                
        except Exception as e:
            with self.output_lock:
                self.banner.display_error(f"AMRfinderPlus analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.cleanup_module_directory(amr_module_path, fasta_files)

    def run_lineage_analysis(self, output_dir: Path) -> bool:
        """Run lineage database generation"""
        try:
            with self.output_lock:
                self.banner.start_analysis_timer("lineage_db")
                self.banner.display_module_header("Lineage Database", "E. coli lineage reference generation")
            
            lineage_module_path = self.base_dir / "modules" / "Ecoli_lineage"
            lineage_script = lineage_module_path / "ecoli_html_reference.py"
            
            if not lineage_script.exists():
                with self.output_lock:
                    self.banner.display_error(f"Lineage script not found at: {lineage_script}")
                return False
            
            # Build command
            cmd = [sys.executable, str(lineage_script)]
            
            with self.output_lock:
                self.banner.display_info("Generating E. coli lineage reference database...")
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=lineage_module_path)
            
            if result.returncode == 0:
                with self.output_lock:
                    self.banner.stop_analysis_timer("lineage_db")
                    self.banner.display_success("E. coli lineage reference database generated!")
                
                # Create lineage output directory
                lineage_output = output_dir / "lineage_results"
                lineage_output.mkdir(parents=True, exist_ok=True)
                
                # Copy lineage HTML to output directory
                lineage_html = lineage_module_path / "ecoli_comprehensive_reference.html"
                if lineage_html.exists():
                    target_html = lineage_output / "ecoli_comprehensive_reference.html"
                    shutil.copy2(lineage_html, target_html)
                    
                    with self.output_lock:
                        self.banner.display_success(f"E. coli lineage reference copied to: {target_html}")
                return True
            else:
                with self.output_lock:
                    self.banner.display_warning("Lineage database generation had warnings")
                    if result.stderr:
                        self.banner.display_info(f"Lineage stderr: {result.stderr[:200]}...")
                return True
                
        except Exception as e:
            with self.output_lock:
                self.banner.display_error(f"Lineage database generation failed: {str(e)}")
            return False

    def run_parallel_analyses(self, fasta_files: List[Path], output_dir: Path, threads: int, 
                            skip_modules: Dict[str, bool]) -> Dict[str, bool]:
        """Run analyses in parallel with synchronized output"""
        # Regular analyses (run in parallel)
        analysis_functions = [
            (self.run_mlst_analysis, "MLST", not skip_modules.get('mlst', False)),
            (self.run_serotyping_analysis, "Serotyping", not skip_modules.get('serotyping', False)),
            (self.run_chtyper_analysis, "CH Typing", not skip_modules.get('chtyper', False)),
            (self.run_phylogrouping_analysis, "Phylogrouping", not skip_modules.get('phylogrouping', False)),
            (self.run_abricate_analysis, "ABRicate", not skip_modules.get('abricate', False))
        ]
        
        # Filter out skipped analyses
        active_analyses = [(func, name) for func, name, enabled in analysis_functions if enabled]
        
        if not active_analyses:
            self.banner.display_warning("All parallel analyses were skipped!")
            return {}
        
        with self.output_lock:
            self.banner.display_info(f"Running {len(active_analyses)} analyses in parallel")
        
        results = {}
        
        # Run analyses in parallel
        with ThreadPoolExecutor(max_workers=min(len(active_analyses), max(1, threads // 2))) as executor:
            future_to_analysis = {
                executor.submit(func, fasta_files, output_dir, max(1, threads // len(active_analyses))): name 
                for func, name in active_analyses
            }
            
            for future in as_completed(future_to_analysis):
                if self.interrupted:
                    break
                    
                analysis_name = future_to_analysis[future]
                
                try:
                    success = future.result()
                    results[analysis_name] = success
                    
                    with self.output_lock:
                        if success:
                            self.banner.display_success(f"✅ {analysis_name} completed")
                        else:
                            self.banner.display_error(f"❌ {analysis_name} failed")
                        
                except Exception as e:
                    with self.output_lock:
                        self.banner.display_error(f"❌ {analysis_name} failed with exception: {str(e)}")
                    results[analysis_name] = False
        
        return results

    def run_complete_analysis(self, input_path: str, output_dir: str, threads: int = 1, 
                            skip_modules: Dict[str, bool] = None):
        """Run complete EcoliTyper analysis pipeline"""
        if skip_modules is None:
            skip_modules = {}
        
        start_time = datetime.now()
        
        try:
            # Display beautiful startup sequence
            self.banner.display_startup_sequence()
            self.banner.display_banner(show_quote=True, show_author=True)
            
            # Create output directory
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            
            # Find input files
            self.fasta_files = self.find_fasta_files(input_path)
            
            if not self.fasta_files:
                self.banner.display_error("No FASTA files found! Analysis stopped.")
                return
            
            # Show file formats detected
            extensions = set(f.suffix.lower() for f in self.fasta_files)
            self.banner.display_success(f"Starting analysis of {len(self.fasta_files)} E. coli genomes")
            self.banner.display_info(f"File formats detected: {', '.join(extensions)}")
            
            # Create output structure
            subdirs = [
                "mlst_results", "serotyping_results", "chtyper_results",
                "phylogrouping_results", "abricate_results", "amrfinder_results", "lineage_results"
            ]
            for subdir in subdirs:
                (output_path / subdir).mkdir(exist_ok=True)
            
            # Display analysis plan
            self.banner.display_module_header("Analysis Plan", "Modules to be executed")
            analyses_to_run = [
                ("MLST", not skip_modules.get('mlst', False)),
                ("Serotyping", not skip_modules.get('serotyping', False)),
                ("CH Typing", not skip_modules.get('chtyper', False)),
                ("Phylogrouping", not skip_modules.get('phylogrouping', False)),
                ("ABRicate", not skip_modules.get('abricate', False)),
                ("AMRfinderPlus", not skip_modules.get('amrfinder', False)),
                ("Lineage Reference", not skip_modules.get('lineage', False))
            ]
            
            for analysis, enabled in analyses_to_run:
                status = "✅ ENABLED" if enabled else "⏸️  SKIPPED"
                print(f"   {status} - {analysis}")
            
            # Run main analyses in parallel (except AMRfinder)
            analysis_results = self.run_parallel_analyses(self.fasta_files, output_path, threads, skip_modules)
            
            # Run AMRfinderPlus analysis LAST (always sequential due to resource usage)
            if not skip_modules.get('amrfinder', False) and not self.interrupted:
                amr_success = self.run_amrfinder_analysis(self.fasta_files, output_path, threads)
                analysis_results["AMRfinderPlus"] = amr_success
            
            # Run lineage analysis if not skipped
            if not skip_modules.get('lineage', False) and not self.interrupted:
                lineage_success = self.run_lineage_analysis(output_path)
                analysis_results["Lineage Reference"] = lineage_success
            
            # Display beautiful completion footer
            successful_count = sum(analysis_results.values())
            total_count = len(analysis_results)
            
            self.banner.display_footer(samples_processed=len(self.fasta_files))
            
            # ADDED: Display citation request and random footer message
            self.banner.display_citation_request()
            self.banner.display_random_footer()
            
            # Final status
            if successful_count == total_count:
                self.banner.display_success(f"🎉 All {total_count} analyses completed successfully!")
                self.banner.display_success("🧹 All module directories have been cleaned up")
            else:
                self.banner.display_warning(
                    f"⚠️  {successful_count}/{total_count} analyses completed successfully."
                )
                
        except KeyboardInterrupt:
            self.banner.display_error("Analysis interrupted by user")
            self._emergency_cleanup()
        except Exception as e:
            self.banner.display_error(f"Critical error in analysis pipeline: {str(e)}")
            import traceback
            traceback.print_exc()

def main():
    """Main entry point for EcoliTyper"""
    parser = argparse.ArgumentParser(
        description="EcoliTyper: Complete E. coli Typing Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  ecolityper -i genome.fna -o results/
  ecolityper -i "*.fna" -o batch_results --threads 8
  ecolityper -i "*.fasta" -o analysis --threads 16 --skip-lineage
  ecolityper -i "genome*.fa" -o results/ --threads 4

Supported FASTA formats: .fna, .fasta, .fa, .fsa

Analysis Modules:
  • MLST (Multi-Locus Sequence Typing)
  • Serotyping (O and H antigen determination)
  • CH Typing (FumC and FimH typing)  
  • Phylogrouping (zClermont algorithm)
  • ABRicate (Resistance/Virulence/Plasmid screening)
  • AMRfinderPlus (NCBI AMR gene detection) 
  • Lineage reference database


  Please run the following commands prior to analysis:
  • amrfinder -u (For database update)
  • abricate --setupdb (For latest database)

Output: Comprehensive results for all analyses in organized directories
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input FASTA file(s) - can use glob patterns like "*.fna" or "*.fasta"')
    parser.add_argument('-o', '--output', required=True,
                       help='Output directory for all results')
    parser.add_argument('-t', '--threads', type=int, default=2,
                       help='Number of threads (default: 2)')
    
    # Skip options
    parser.add_argument('--skip-amrfinder', action='store_true', 
                       help='Skip AMRfinderPlus analysis')
    parser.add_argument('--skip-abricate', action='store_true',
                       help='Skip ABRicate analysis')
    parser.add_argument('--skip-mlst', action='store_true',
                       help='Skip MLST analysis')
    parser.add_argument('--skip-serotyping', action='store_true',
                       help='Skip serotyping analysis')
    parser.add_argument('--skip-chtyper', action='store_true',
                       help='Skip CH typing analysis')
    parser.add_argument('--skip-phylogrouping', action='store_true',
                       help='Skip phylogrouping analysis')
    parser.add_argument('--skip-lineage', action='store_true',
                       help='Skip lineage reference generation')
    
    args = parser.parse_args()
    
    # Create skip modules dictionary
    skip_modules = {
        'amrfinder': args.skip_amrfinder,
        'abricate': args.skip_abricate,
        'mlst': args.skip_mlst,
        'serotyping': args.skip_serotyping,
        'chtyper': args.skip_chtyper,
        'phylogrouping': args.skip_phylogrouping,
        'lineage': args.skip_lineage
    }
    
    # Create and run EcoliTyper
    ecolityper = EcoliTyperOrchestrator()
    
    try:
        ecolityper.run_complete_analysis(
            input_path=args.input,
            output_dir=args.output,
            threads=args.threads,
            skip_modules=skip_modules
        )
    except KeyboardInterrupt:
        print("\n❌ Analysis interrupted by user - automatic cleanup completed")
    except Exception as e:
        print(f"\n💥 Critical error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()