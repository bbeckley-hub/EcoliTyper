#!/usr/bin/env python3
"""
EcoliTyper Banner Module - FIXED VERSION
Beautiful ASCII art and scientific quotes for terminal display - StaphScope Style
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
Date: 2025/2026-05-15
Send a quick mail for any issues or further explanations.
"""

import random
from datetime import datetime
import sys
import time
import textwrap
import os
import shutil

class EcoliTyperBanner:
    """EcoliTyper Banner Display with Scientific Quotes"""
    
    def __init__(self):
        self.banner_art = self._get_banner_art()
        self.quotes = self._get_scientific_quotes()
        self.version = "v1.2.0"
        self.author_info = self._get_author_info()
        self.terminal_width = self._get_terminal_width()
        self.analysis_times = {}  # For timing functions
    
    def _get_terminal_width(self):
        """Get terminal width, default to 88 if cannot determine"""
        try:
            return min(100, shutil.get_terminal_size().columns - 2)
        except:
            return 88
    
    def _get_banner_art(self):
        """Return the main EcoliTyper ASCII art"""
        return r"""
    ███████╗ ██████╗ ██████╗ ██╗     ██╗████████╗██║   ██╗██████╗ ███████╗██████╗ 
    ██╔════╝██╔════╝██╔═══██╗██║     ██║╚══██╔══╝╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
    ██████╗ ██║     ██║   ██║██║     ██║   ██║    ╚████╔╝ ██████╔╝█████╗  ██████╔╝
    ██╔══╝  ██║     ██║   ██║██║     ██║   ██║     ╚██╔╝  ██╔═══╝ ██╔══╝  ██╔══██╗
    ███████╗╚██████╗╚██████╔╝███████╗██║   ██║      ██║   ██║     ███████╗██║  ██║
    ╚══════╝ ╚═════╝ ╚═════╝ ╚══════╝╚═╝   ╚═╝      ╚═╝   ╚═╝     ╚══════╝╚═╝  ╚═╝
    
        🧫 Comprehensive Escherichia coli Typing & Surveillance Platform
    """
    
    def _get_scientific_quotes(self):
        """Return collection of scientific quotes about microbiology and E. coli"""
        return [
            {
                "quote": "The important thing in science is not so much to obtain new facts as to discover new ways of thinking about them.",
                "author": "William Lawrence Bragg"
            },
            {
                "quote": "In the fields of observation chance favors only the prepared mind.",
                "author": "Louis Pasteur"
            },
            {
                "quote": "E. coli is not just a laboratory workhorse; it's a window into the fundamental processes of life.",
                "author": "Molecular Biology Principle"
            },
            {
                "quote": "The genome is the book of life, and we are learning to read it with increasing clarity.",
                "author": "Francis Collins"
            },
            {
                "quote": "Understanding bacterial pathogenesis requires knowing not just the pathogen, but its lineage and evolution.",
                "author": "Infectious Disease Research"
            },
            {
                "quote": "Every E. coli strain tells a story of adaptation, resistance, and survival.",
                "author": "Microbial Genomics"
            },
            {
                "quote": "The art of epidemiology lies in connecting genetic markers to public health outcomes.",
                "author": "Public Health Research"
            },
            {
                "quote": "DNA sequencing has revolutionized our understanding of bacterial diversity and evolution.",
                "author": "Genomics Research"
            },
            {
                "quote": "In microbiology, the smallest details often reveal the biggest truths.",
                "author": "Laboratory Wisdom"
            },
            {
                "quote": "The fight against antimicrobial resistance begins with understanding resistance mechanisms.",
                "author": "AMR Research"
            }
        ]
    
    def _get_author_info(self):
        """Return author information"""
        return {
            "name": "Brown Beckley",
            "github": "bbeckley-hub",
            "email": "brownbeckley94@gmail.com",
            "affiliation": "University of Ghana Medical School - Department of Medical Biochemistry",
            "license": "MIT"
        }
    
    def _get_colors(self):
        """Define color codes for terminal output"""
        class Colors:
            RED = '\033[91m'
            GREEN = '\033[92m'
            YELLOW = '\033[93m'
            BLUE = '\033[94m'
            MAGENTA = '\033[95m'
            CYAN = '\033[96m'
            WHITE = '\033[97m'
            BOLD = '\033[1m'
            UNDERLINE = '\033[4m'
            DIM = '\033[2m'
            END = '\033[0m'
            # Extended colors
            LIGHT_BLUE = '\033[38;5;117m'
            LIGHT_GREEN = '\033[38;5;120m'
            LIGHT_YELLOW = '\033[38;5;228m'
            ORANGE = '\033[38;5;214m'
            PURPLE = '\033[38;5;141m'
            PINK = '\033[38;5;213m'
        return Colors
    
    def display_banner(self, show_quote=True, show_author=True):
        """Display the main EcoliTyper banner"""
        C = self._get_colors()
        width = self.terminal_width
        
        # Main banner
        print(f"\n{C.CYAN}{C.BOLD}{self.banner_art}{C.END}")
        
        # Top decoration
        print(f"{C.CYAN}{'▓' * width}{C.END}")
        print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
        
        # Version and date
        version_text = f"🔬 Version: {self.version}  |  {datetime.now().strftime('%Y-%m-%d')}  |  EcoliTyper Platform 🧬"
        print(f"{C.YELLOW}{C.BOLD}{version_text.center(width)}{C.END}")
        
        print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
        print(f"{C.CYAN}{'▓' * width}{C.END}\n")
        
        if show_quote:
            quote = random.choice(self.quotes)
            
            # Quote section header
            print(f"{C.CYAN}{'▓' * width}{C.END}")
            print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
            print(f"{C.CYAN}{C.BOLD}{'💡 SCIENTIFIC INSPIRATION'.center(width)}{C.END}")
            print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
            print(f"{C.CYAN}{'▓' * width}{C.END}\n")
            
            # Quote text with wrapping
            wrapper = textwrap.TextWrapper(width=width-4, initial_indent='  ', subsequent_indent='  ')
            quote_lines = wrapper.wrap(f'"{quote["quote"]}"')
            
            for line in quote_lines:
                print(f"{C.WHITE}{line}{C.END}")
            
            print()
            author_text = f"— {quote['author']}"
            print(f"{C.YELLOW}{C.BOLD}{author_text.rjust(width-2)}{C.END}")
            
            print(f"\n{C.CYAN}{'▓' * width}{C.END}")
            print(f"{C.LIGHT_BLUE}{'─' * width}{C.END}\n")
        
        if show_author:
            # Author section header
            print(f"{C.CYAN}{'▓' * width}{C.END}")
            print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
            print(f"{C.CYAN}{C.BOLD}{'👨‍💻 DEVELOPER & CONTACT INFORMATION'.center(width)}{C.END}")
            print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
            print(f"{C.CYAN}{'▓' * width}{C.END}\n")
            
            # Author details in columns
            print(f"{C.LIGHT_BLUE}  👤 Name:{C.END}        {C.WHITE}{self.author_info['name']}{C.END}")
            print(f"{C.LIGHT_BLUE}  🐙 GitHub:{C.END}      {C.WHITE}{self.author_info['github']}{C.END}")
            print(f"{C.LIGHT_BLUE}  📧 Email:{C.END}       {C.WHITE}{self.author_info['email']}{C.END}")
            print(f"{C.LIGHT_BLUE}  🏛️  Affiliation:{C.END} {C.WHITE}{self.author_info['affiliation']}{C.END}")
            print(f"{C.LIGHT_BLUE}  📜 License:{C.END}     {C.WHITE}{self.author_info['license']}{C.END}")
            
            print(f"\n{C.CYAN}{'▓' * width}{C.END}")
            print(f"{C.LIGHT_BLUE}{'─' * width}{C.END}\n")
    
    def display_startup_sequence(self):
        """Display animated startup sequence with progress - FIXED PROGRESS BARS"""
        C = self._get_colors()
        width = self.terminal_width
        
        # Header
        print(f"\n{C.CYAN}{'▓' * width}{C.END}")
        print(f"{C.CYAN}{'█' * width}{C.END}")
        print(f"{C.CYAN}{C.BOLD}{'🚀 INITIALIZING ECOLITYPER ANALYSIS PLATFORM 🚀'.center(width)}{C.END}")
        print(f"{C.CYAN}{'█' * width}{C.END}")
        print(f"{C.CYAN}{'▓' * width}{C.END}")
        print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}\n")
        time.sleep(0.3)
        
        steps = [
            ("🗄️", "Loading E. coli genomic databases", C.BLUE),
            ("🧬", "Initializing MLST analysis engine", C.CYAN),
            ("🔍", "Configuring serotyping algorithms", C.GREEN),
            ("🧬", "Setting up CH typing analysis", C.YELLOW),
            ("🌳", "Preparing zClermont phylogrouping", C.ORANGE),
            ("🛡️", "Enabling Abricate analysis (Resistance/Virulence/Plasmids)", C.MAGENTA),
            ("💊", "Configuring AMRfinderPlus (NCBI)", C.LIGHT_BLUE),
            ("⚡", "Optimizing multi-threading capabilities", C.LIGHT_GREEN),
        ]
        
        for i, (icon, step, color) in enumerate(steps, 1):
            print(f"{color}[{i}/{len(steps)}] {icon}  {step}...{C.END}")
            
            # FIXED ANIMATION
            progress_stages = [
                '░░░░░░░░░░░░░░░░░░░░',
                '▓░░░░░░░░░░░░░░░░░░░', 
                '▓▓▓░░░░░░░░░░░░░░░░░',
                '▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓'
            ]
            
            # Show all 4 progress stages with proper animation
            for j, stage in enumerate(progress_stages):
                if j < len(progress_stages) - 1:
                    # Show intermediate stage with carriage return
                    print(f"{color}  {stage}{C.END}", end='\r', flush=True)
                    time.sleep(0.08)
                else:
                    # Final stage - print with checkmark and newline
                    print(f"{color}  {stage} {C.GREEN}✓{C.END}")
            
            time.sleep(0.1)
        
        # Footer
        print(f"\n{C.LIGHT_BLUE}{'═' * width}{C.END}")
        print(f"{C.GREEN}{'▓' * width}{C.END}")
        print(f"{C.GREEN}{C.BOLD}{'✅ ECOLITYPER READY FOR ANALYSIS! ✅'.center(width)}{C.END}")
        print(f"{C.GREEN}{'▓' * width}{C.END}")
        print(f"{C.LIGHT_GREEN}{'═' * width}{C.END}\n")
        time.sleep(0.3)
    
    # TIMING METHODS - From the orchestrator
    def start_analysis_timer(self, analysis_name):
        """Start timer for a specific analysis"""
        self.analysis_times[analysis_name] = {
            'start': datetime.now(),
            'end': None,
            'duration': None
        }
    
    def stop_analysis_timer(self, analysis_name):
        """Stop timer for a specific analysis and calculate duration"""
        if analysis_name in self.analysis_times:
            end_time = datetime.now()
            self.analysis_times[analysis_name]['end'] = end_time
            duration = end_time - self.analysis_times[analysis_name]['start']
            self.analysis_times[analysis_name]['duration'] = duration
            return duration
        return None
    
    def get_analysis_time(self, analysis_name):
        """Get formatted analysis time"""
        if analysis_name in self.analysis_times and self.analysis_times[analysis_name]['duration']:
            duration = self.analysis_times[analysis_name]['duration']
            total_seconds = duration.total_seconds()
            
            if total_seconds < 60:
                return f"{total_seconds:.1f} seconds"
            elif total_seconds < 3600:
                minutes = total_seconds // 60
                seconds = total_seconds % 60
                return f"{int(minutes)}m {int(seconds)}s"
            else:
                hours = total_seconds // 3600
                minutes = (total_seconds % 3600) // 60
                return f"{int(hours)}h {int(minutes)}m"
        return "Not completed"
    
    def display_footer(self, analysis_time=None, samples_processed=0):
        """Display analysis completion footer"""
        C = self._get_colors()
        width = self.terminal_width
        
        print()
        # Header
        print(f"{C.CYAN}{'▓' * width}{C.END}")
        print(f"{C.LIGHT_BLUE}{'━' * width}{C.END}")
        print(f"{C.CYAN}{C.BOLD}{'🎉 ANALYSIS COMPLETE 🎉'.center(width)}{C.END}")
        print(f"{C.LIGHT_BLUE}{'━' * width}{C.END}")
        print(f"{C.CYAN}{'▓' * width}{C.END}\n")
        
        # Statistics
        if analysis_time or samples_processed > 0:
            print(f"{C.CYAN}{'▓' * width}{C.END}")
            print(f"{C.LIGHT_BLUE}{'─' * width}{C.END}")
            print(f"{C.YELLOW}{C.BOLD}{'📊 STATISTICS'.center(width)}{C.END}")
            print(f"{C.LIGHT_BLUE}{'─' * width}{C.END}")
            print(f"{C.CYAN}{'▓' * width}{C.END}\n")
            
            if analysis_time:
                print(f"{C.WHITE}  ⏱️  Analysis Duration: {C.GREEN}{C.BOLD}{analysis_time}{C.END}")
            if samples_processed > 0:
                print(f"{C.WHITE}  🧫 E. coli Genomes Processed:  {C.GREEN}{C.BOLD}{samples_processed}{C.END}")
            
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"{C.WHITE}  📅 Completion Time:    {C.CYAN}{C.BOLD}{current_time}{C.END}")
            
            print(f"\n{C.CYAN}{'▓' * width}{C.END}")
            print(f"{C.LIGHT_BLUE}{'─' * width}{C.END}\n")
        
        # Support section
        print(f"{C.CYAN}{'▓' * width}{C.END}")
        print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
        print(f"{C.CYAN}{C.BOLD}{'📞 SUPPORT & INQUIRIES'.center(width)}{C.END}")
        print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
        print(f"{C.CYAN}{'▓' * width}{C.END}\n")
        
        print(f"{C.WHITE}  👤 Contact:  {C.LIGHT_BLUE}{self.author_info['name']}{C.END}")
        print(f"{C.WHITE}  🐙 GitHub:   {C.LIGHT_BLUE}{self.author_info['github']}{C.END}")
        print(f"{C.WHITE}  📧 Email:    {C.LIGHT_BLUE}{self.author_info['email']}{C.END}")
        
        print(f"\n{C.CYAN}{'▓' * width}{C.END}")
        print(f"{C.LIGHT_BLUE}{'─' * width}{C.END}\n")
    
    def display_module_header(self, module_name, description=""):
        """Display header for specific analysis modules"""
        C = self._get_colors()
        width = self.terminal_width
        
        print(f"\n{C.CYAN}{'▓' * width}{C.END}")
        print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
        
        # Module title
        title = f"🧬 {module_name.upper()}"
        print(f"{C.YELLOW}{C.BOLD}{title.center(width)}{C.END}")
        
        if description:
            print(f"{C.LIGHT_BLUE}{'─' * width}{C.END}")
            # Wrap long descriptions
            wrapper = textwrap.TextWrapper(width=width-4, initial_indent='  ', subsequent_indent='  ')
            desc_lines = wrapper.wrap(description)
            for line in desc_lines:
                print(f"{C.WHITE}{line}{C.END}")
        
        print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
        print(f"{C.CYAN}{'▓' * width}{C.END}\n")
        
        # Flush output to ensure header displays immediately
        sys.stdout.flush()
    
    def display_section_divider(self, title="", color=None):
        """Display a colorful section divider"""
        C = self._get_colors()
        section_color = color or C.CYAN
        width = self.terminal_width
        
        print(f"\n{section_color}{'═' * width}{C.END}")
        if title:
            print(f"{section_color}{C.BOLD}{title.center(width)}{C.END}")
            print(f"{section_color}{'═' * width}{C.END}\n")
    
    def display_warning(self, message):
        """Display warning message"""
        C = self._get_colors()
        width = self.terminal_width
        
        # Wrap long messages
        wrapper = textwrap.TextWrapper(width=width-4, initial_indent='  ', subsequent_indent='  ')
        message_lines = wrapper.wrap(message)
        
        print(f"{C.YELLOW}{'─' * width}{C.END}")
        if len(message_lines) == 1:
            print(f"{C.YELLOW}{C.BOLD}⚠️  WARNING:{C.END} {C.YELLOW}{message}{C.END}")
        else:
            print(f"{C.YELLOW}{C.BOLD}⚠️  WARNING:{C.END}")
            for line in message_lines:
                print(f"{C.YELLOW}  {line}{C.END}")
        print(f"{C.YELLOW}{'─' * width}{C.END}")
    
    def display_error(self, message):
        """Display error message"""
        C = self._get_colors()
        width = self.terminal_width
        
        # Wrap long messages
        wrapper = textwrap.TextWrapper(width=width-4, initial_indent='  ', subsequent_indent='  ')
        message_lines = wrapper.wrap(message)
        
        print(f"{C.RED}{'─' * width}{C.END}")
        if len(message_lines) == 1:
            print(f"{C.RED}{C.BOLD}❌ ERROR:{C.END} {C.RED}{message}{C.END}")
        else:
            print(f"{C.RED}{C.BOLD}❌ ERROR:{C.END}")
            for line in message_lines:
                print(f"{C.RED}  {line}{C.END}")
        print(f"{C.RED}{'─' * width}{C.END}")
    
    def display_success(self, message):
        """Display success message"""
        C = self._get_colors()
        width = self.terminal_width
        
        # Wrap long messages
        wrapper = textwrap.TextWrapper(width=width-4, initial_indent='  ', subsequent_indent='  ')
        message_lines = wrapper.wrap(message)
        
        print(f"{C.GREEN}{'─' * width}{C.END}")
        if len(message_lines) == 1:
            print(f"{C.GREEN}{C.BOLD}✅ SUCCESS:{C.END} {C.GREEN}{message}{C.END}")
        else:
            print(f"{C.GREEN}{C.BOLD}✅ SUCCESS:{C.END}")
            for line in message_lines:
                print(f"{C.GREEN}  {line}{C.END}")
        print(f"{C.GREEN}{'─' * width}{C.END}")
    
    def display_info(self, message):
        """Display info message"""
        C = self._get_colors()
        width = self.terminal_width
        
        # Wrap long messages
        wrapper = textwrap.TextWrapper(width=width-4, initial_indent='  ', subsequent_indent='  ')
        message_lines = wrapper.wrap(message)
        
        print(f"{C.CYAN}{'─' * width}{C.END}")
        if len(message_lines) == 1:
            print(f"{C.CYAN}{C.BOLD}💡 INFO:{C.END} {C.CYAN}{message}{C.END}")
        else:
            print(f"{C.CYAN}{C.BOLD}💡 INFO:{C.END}")
            for line in message_lines:
                print(f"{C.CYAN}  {line}{C.END}")
        print(f"{C.CYAN}{'─' * width}{C.END}")
    
    def display_progress_bar(self, iteration, total, prefix='', suffix='', length=50, fill='█'):
        """Display enhanced progress bar"""
        C = self._get_colors()
        
        # Don't display if we're at 0% or 100% and suffix is empty
        if iteration == 0 and not suffix:
            return
            
        percent = ("{0:.1f}").format(100 * (iteration / float(total)))
        filled_length = int(length * iteration // total)
        bar = fill * filled_length + '░' * (length - filled_length)
        
        # Color gradient based on progress
        if iteration < total * 0.33:
            color = C.RED
        elif iteration < total * 0.66:
            color = C.YELLOW
        else:
            color = C.GREEN
        
        # For intermediate updates, use \r to overwrite
        print(f'\r{C.CYAN}{prefix}{C.END} {color}[{bar}]{C.END} {C.BOLD}{percent}%{C.END} {C.DIM}{suffix}{C.END}', end='\r')
        
        # If this is the final iteration, print newline
        if iteration == total:
            print()
            sys.stdout.flush()
    
    # Additional methods for compatibility
    def display_citation_request(self):
        """Display citation request"""
        C = self._get_colors()
        width = self.terminal_width
        
        print(f"\n{C.MAGENTA}{'═' * width}{C.END}")
        print(f"{C.MAGENTA}{C.BOLD}{'📚 CITATION REQUEST'.center(width)}{C.END}")
        print(f"{C.MAGENTA}{'═' * width}{C.END}")
        
        messages = [
            "If you use EcoliTyper in your research, please cite:",
            "",
            "EcoliTyper: A species-optimized computational pipeline",
            "for comprehensive genotyping and surveillance of Escherichia coli.",
            "",
            "GitHub: https://github.com/bbeckley-hub/EcoliTyper"
        ]
        
        for message in messages:
            if message:
                print(f"{C.WHITE}{message.center(width)}{C.END}")
            else:
                print()
        
        print(f"{C.MAGENTA}{'═' * width}{C.END}\n")
    
    def display_random_footer(self):
        """Display random footer message about E. coli genomics"""
        C = self._get_colors()
        
        footer_messages = [
            "🔬 Advancing E. coli genomics research to combat antimicrobial resistance worldwide.",
            "🌍 Contributing to global AMR surveillance through comprehensive E. coli typing.",
            "💡 Harnessing genomic data for better understanding of E. coli epidemiology.",
            "🦠 Bridging genomics and clinical practice in infectious disease management.",
            "🧪 Pioneering open-source tools for accessible bacterial genomics research.",
        ]
        
        message = random.choice(footer_messages)
        width = self.terminal_width
        
        print(f"\n{C.CYAN}{'═' * width}{C.END}")
        print(f"{C.CYAN}{C.BOLD}{'✨ ' + message + ' ✨'.center(width)}{C.END}")
        print(f"{C.CYAN}{'═' * width}{C.END}\n")

def main():
    """Test the fixed banner display"""
    banner = EcoliTyperBanner()
    
    # Display full banner with startup sequence
    banner.display_startup_sequence()
    banner.display_banner(show_quote=True, show_author=True)
    
    # Test module headers
    banner.display_module_header("MLST Analysis", "Multi-Locus Sequence Typing for E. coli")
    banner.display_info("Copied 1 files to MLST module")
    banner.display_info("Running MLST analysis with pattern: '*.fna'")
    banner.display_success("MLST analysis completed!")
    
    print()
    banner.display_module_header("Serotyping Analysis", "O and H antigen determination")
    banner.display_info("Copied 1 files to serotyping module")
    banner.display_success("Serotyping analysis completed!")
    
    # Test progress bar
    print("\n🔄 Progress demonstration:")
    for i in range(101):
        banner.display_progress_bar(i, 100, prefix='Analyzing E. coli:', suffix='Complete', length=50)
        time.sleep(0.02)
    
    # Test timing methods
    banner.start_analysis_timer("Test Analysis")
    time.sleep(1)
    banner.stop_analysis_timer("Test Analysis")
    print(f"\nAnalysis time: {banner.get_analysis_time('Test Analysis')}")
    
    # Display footer
    banner.display_footer(analysis_time="2 minutes, 15 seconds", samples_processed=8)
    
    # Test citation and footer
    banner.display_citation_request()
    banner.display_random_footer()

if __name__ == "__main__":
    main()