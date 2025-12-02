#!/usr/bin/env python3
"""
EcoliTyper Banner Module
Beautiful ASCII art and scientific quotes for terminal display
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
Date: 2025
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
    """EcoliTyper Banner Display with Scientific Quotes and Time Tracking"""
    
    def __init__(self):
        self.banner_art = self._get_banner_art()
        self.quotes = self._get_scientific_quotes()
        self.footer_messages = self._get_footer_messages()
        self.version = "v1.0.0"
        self.author_info = self._get_author_info()
        self.analysis_times = {}  # Track analysis times
    
    def _get_optimal_width(self, content_type="general"):
        """Get optimal width based on content type"""
        try:
            terminal_width = shutil.get_terminal_size().columns
            
            if content_type == "quote":
                return min(max(90, terminal_width - 20), 110)
            elif content_type == "author":
                return min(max(85, terminal_width - 25), 100)
            elif content_type == "footer":
                return min(max(95, terminal_width - 15), 120)
            elif content_type == "citation":
                return min(max(80, terminal_width - 20), 100)
            else:
                return 80
        except:
            return {
                "quote": 100,
                "author": 90, 
                "footer": 100,
                "citation": 85,
                "general": 80
            }[content_type]
    
    def _get_banner_art(self):
        """Return the main EcoliTyper ASCII art"""
        return r"""
    ███████╗ ██████╗ ██████╗ ██╗     ██╗████████╗██╗   ██╗██████╗ ███████╗██████╗ 
    ██╔════╝██╔════╝██╔═══██╗██║     ██║╚══██╔══╝╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
    ██████╗ ██║     ██║   ██║██║     ██║   ██║    ╚████╔╝ ██████╔╝█████╗  ██████╔╝
    ██╔══╝  ██║     ██║   ██║██║     ██║   ██║     ╚██╔╝  ██╔═══╝ ██╔══╝  ██╔══██╗
    ███████╗╚██████╗╚██████╔╝███████╗██║   ██║      ██║   ██║     ███████╗██║  ██║
    ╚══════╝ ╚═════╝ ╚═════╝ ╚══════╝╚═╝   ╚═╝      ╚═╝   ╚═╝     ╚══════╝╚═╝  ╚═╝
    
        🧫 Comprehensive E. coli Typing: MLST • Serotyping • CH Typing • Phylogrouping • Resistance Gene Profiling 
                                       • Virulence & Plasmid Profiling • Lineage Db • AMR Screening
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

    def _get_footer_messages(self):
        """Return collection of footer messages about E. coli genomics"""
        return [
            "🔬 Advancing E. coli genomics research to combat antimicrobial resistance worldwide.",
            "🌍 Contributing to global AMR surveillance through comprehensive E. coli typing.",
            "💡 Harnessing genomic data for better understanding of E. coli epidemiology.",
            "🦠 Bridging genomics and clinical practice in infectious disease management.",
            "🧪 Pioneering open-source tools for accessible bacterial genomics research.",
            "📊 Transforming raw sequences into actionable public health intelligence.",
            "🔍 Uncovering the hidden stories within E. coli genomes for better outcomes.",
            "💊 Fighting antimicrobial resistance one genome at a time through precise typing.",
            "🧬 Decoding E. coli diversity to improve global health security.",
            "⚕️ From sequence to insight: Advancing clinical microbiology through genomics."
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
            END = '\033[0m'
        return Colors
    
    def _create_box(self, width, char='═', corner_left='╔', corner_right='╗'):
        """Create a horizontal box line"""
        return f"{corner_left}{char * (width - 2)}{corner_right}"
    
    def _create_text_line(self, text, width, color=None, align='left', padding=2):
        """Create a text line within the box"""
        C = self._get_colors()
        text_color = color or C.WHITE
        
        # Calculate available width (subtract borders and padding)
        available_width = width - (padding * 2) - 2  # -2 for borders
        
        # Handle text wrapping
        if len(text) > available_width:
            text = textwrap.fill(text, width=available_width)
        
        lines = text.split('\n') if '\n' in text else [text]
        
        formatted_lines = []
        for line in lines:
            if align == 'center':
                formatted_text = line.center(available_width)
            elif align == 'right':
                formatted_text = line.rjust(available_width)
            else:  # left align
                formatted_text = line.ljust(available_width)
            
            formatted_lines.append(f"║{' ' * padding}{text_color}{formatted_text}{C.END}{' ' * padding}║")
        
        return formatted_lines
    
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
    
    def display_banner(self, show_quote=True, show_author=True):
        """Display the main EcoliTyper banner"""
        C = self._get_colors()
        
        print(f"{C.CYAN}{C.BOLD}{self.banner_art}{C.END}")
        print(f"{C.YELLOW}{C.BOLD}         Version: {self.version}{C.END}")
        print()
        
        if show_quote:
            quote = random.choice(self.quotes)
            quote_width = self._get_optimal_width("quote")
            
            # Quote box
            print(f"{C.GREEN}{self._create_box(quote_width, '═', '╔', '╗')}{C.END}")
            
            # Title
            title_lines = self._create_text_line("💡 SCIENTIFIC QUOTE", quote_width, C.CYAN + C.BOLD, 'left')
            for line in title_lines:
                print(line)
            
            print(f"{C.GREEN}{self._create_box(quote_width, '─', '╠', '╣')}{C.END}")
            
            # Quote text
            quote_lines = self._create_text_line(quote['quote'], quote_width, C.WHITE, 'left')
            for line in quote_lines:
                print(line)
            
            # Author
            author_lines = self._create_text_line(f"— {quote['author']}", quote_width, C.YELLOW, 'right')
            for line in author_lines:
                print(line)
            
            print(f"{C.GREEN}{self._create_box(quote_width, '═', '╚', '╝')}{C.END}")
            print()
        
        if show_author:
            author_width = self._get_optimal_width("author")
            
            # Author info box
            print(f"{C.GREEN}{self._create_box(author_width, '═', '╔', '╗')}{C.END}")
            
            # Title
            title_lines = self._create_text_line("👨‍💻 AUTHOR INFORMATION", author_width, C.MAGENTA + C.BOLD, 'left')
            for line in title_lines:
                print(line)
            
            print(f"{C.GREEN}{self._create_box(author_width, '─', '╠', '╣')}{C.END}")
            
            # Author details
            details = [
                f"Name: {self.author_info['name']}",
                f"GitHub: {self.author_info['github']}",
                f"Email: {self.author_info['email']}",
                f"Affiliation: {self.author_info['affiliation']}",
                f"License: {self.author_info['license']}"
            ]
            
            for detail in details:
                detail_lines = self._create_text_line(detail, author_width, C.WHITE, 'left')
                for line in detail_lines:
                    print(line)
            
            print(f"{C.GREEN}{self._create_box(author_width, '═', '╚', '╝')}{C.END}")
        print()
    
    def display_startup_sequence(self):
        """Display animated startup sequence"""
        C = self._get_colors()
        
        print(f"{C.CYAN}{C.BOLD}🚀 Initializing EcoliTyper Analysis Platform...{C.END}")
        time.sleep(0.5)
        
        steps = [
            "Loading E. coli genomic databases...",
            "Initializing MLST analysis engine...",
            "Configuring serotyping algorithms...",
            "Setting up CH typing analysis...",
            "Preparing zClermont phylogrouping...",
            "Enabling Abricate analysis (Resistance/Virulence/Plasmids)...",
            "Configuring AMRfinderPlus (NCBI)...",
            "Building lineage database reference...",
            "Optimizing multi-threading capabilities...",
        ]
        
        for i, step in enumerate(steps, 1):
            print(f"{C.YELLOW}[{i}/{len(steps)}] {step}{C.END}")
            time.sleep(0.3)
        
        print(f"{C.GREEN}{C.BOLD}✅ EcoliTyper ready for comprehensive analysis!{C.END}")
        print()
    
    def display_citation_request(self):
        """Display citation request box"""
        C = self._get_colors()
        
        citation_width = self._get_optimal_width("citation")
        
        print()
        print(f"{C.MAGENTA}{self._create_box(citation_width, '═', '╔', '╗')}{C.END}")
        
        # Title
        title_lines = self._create_text_line("📚 CITATION REQUEST", citation_width, C.MAGENTA + C.BOLD, 'center')
        for line in title_lines:
            print(line)
        
        print(f"{C.MAGENTA}{self._create_box(citation_width, '─', '╠', '╣')}{C.END}")
        
        # Citation message
        messages = [
            "If you use EcoliTyper in your research, please cite:",
            "",
            "EcoliTyper: A species-optimized computational pipeline",
            "for comprehensive genotyping and surveillance of Escherichia coli.",
            "",
            "Citation details will be updated upon manuscript",
            "acceptance. Please check GitHub for updates.",
            "",
            "GitHub: https://github.com/bbeckley-hub/EcoliTyper"
        ]
        
        for message in messages:
            message_lines = self._create_text_line(message, citation_width, C.WHITE, 'center')
            for line in message_lines:
                print(line)
        
        print(f"{C.MAGENTA}{self._create_box(citation_width, '═', '╚', '╝')}{C.END}")
        print()
    
    def display_random_footer(self):
        """Display random footer message about E. coli genomics"""
        C = self._get_colors()
        
        message = random.choice(self.footer_messages)
        
        # Simple centered message without box
        print(f"{C.CYAN}{C.BOLD}✨ {message}{C.END}")
        print()
    
    def display_footer(self, samples_processed=0):
        """Display analysis completion footer"""
        C = self._get_colors()
        
        footer_width = self._get_optimal_width("footer")
        
        print()
        print(f"{C.MAGENTA}{self._create_box(footer_width, '═', '╔', '╗')}{C.END}")
        
        # Title
        title_lines = self._create_text_line("ANALYSIS COMPLETE - ECOLITYPER", footer_width, C.MAGENTA + C.BOLD, 'center')
        for line in title_lines:
            print(line)
        
        print(f"{C.MAGENTA}{self._create_box(footer_width, '═', '╠', '╣')}{C.END}")
        
        # Analysis details - Show total time instead of individual timings
        total_duration = sum(
            (t['duration'].total_seconds() if t['duration'] else 0) 
            for t in self.analysis_times.values()
        )
        
        if total_duration > 0:
            if total_duration < 60:
                total_str = f"{total_duration:.1f} seconds"
            elif total_duration < 3600:
                minutes = total_duration // 60
                seconds = total_duration % 60
                total_str = f"{int(minutes)} minutes {int(seconds)} seconds"
            else:
                hours = total_duration // 3600
                minutes = (total_duration % 3600) // 60
                total_str = f"{int(hours)} hours {int(minutes)} minutes"
            
            time_lines = self._create_text_line(f"⏱️  Total Analysis Time: {total_str}", footer_width, C.CYAN, 'left')
            for line in time_lines:
                print(line)
        
        if samples_processed > 0:
            samples_lines = self._create_text_line(f"🧫 E. coli Genomes Processed: {samples_processed}", footer_width, C.GREEN, 'left')
            for line in samples_lines:
                print(line)
        
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        time_lines = self._create_text_line(f"📅 Completion Time: {current_time}", footer_width, C.YELLOW, 'left')
        for line in time_lines:
            print(line)
        
        print(f"{C.MAGENTA}{self._create_box(footer_width, '═', '╠', '╣')}{C.END}")
        
        # Support information
        support_lines = self._create_text_line("TECHNICAL SUPPORT & INQUIRIES:", footer_width, C.WHITE + C.BOLD, 'left')
        for line in support_lines:
            print(line)
        
        author_lines = self._create_text_line(f"Author: {self.author_info['name']}", footer_width, C.CYAN, 'left')
        for line in author_lines:
            print(line)
        
        github_lines = self._create_text_line(f"GitHub: {self.author_info['github']}", footer_width, C.BLUE, 'left')
        for line in github_lines:
            print(line)
        
        email_lines = self._create_text_line(f"Email: {self.author_info['email']}", footer_width, C.GREEN, 'left')
        for line in email_lines:
            print(line)
        
        affiliation_lines = self._create_text_line(f"Affiliation: {self.author_info['affiliation']}", footer_width, C.WHITE, 'left')
        for line in affiliation_lines:
            print(line)
        
        print(f"{C.MAGENTA}{self._create_box(footer_width, '═', '╚', '╝')}{C.END}")
        print()
    
    def display_module_header(self, module_name, description=""):
        """Display header for specific analysis modules"""
        C = self._get_colors()
        
        module_width = self._get_optimal_width("general")
        
        print()
        print(f"{C.BLUE}{self._create_box(module_width, '═', '╔', '╗')}{C.END}")
        
        # Module name
        module_lines = self._create_text_line(f"🧬 {module_name.upper()}", module_width, C.BLUE + C.BOLD, 'left')
        for line in module_lines:
            print(line)
        
        if description:
            print(f"{C.BLUE}{self._create_box(module_width, '─', '╠', '╣')}{C.END}")
            # Description
            desc_lines = self._create_text_line(description, module_width, C.WHITE, 'left')
            for line in desc_lines:
                print(line)
        
        print(f"{C.BLUE}{self._create_box(module_width, '═', '╚', '╝')}{C.END}")
        print()
    
    def display_warning(self, message):
        """Display warning message"""
        C = self._get_colors()
        print(f"{C.YELLOW}⚠️  WARNING: {message}{C.END}")
    
    def display_error(self, message):
        """Display error message"""
        C = self._get_colors()
        print(f"{C.RED}❌ ERROR: {message}{C.END}")
    
    def display_success(self, message):
        """Display success message"""
        C = self._get_colors()
        print(f"{C.GREEN}✅ SUCCESS: {message}{C.END}")
    
    def display_info(self, message):
        """Display info message"""
        C = self._get_colors()
        print(f"{C.CYAN}💡 INFO: {message}{C.END}")
    
    def display_progress_bar(self, iteration, total, prefix='', suffix='', length=50, fill='█'):
        """Display progress bar"""
        C = self._get_colors()
        percent = ("{0:.1f}").format(100 * (iteration / float(total)))
        filled_length = int(length * iteration // total)
        bar = fill * filled_length + '-' * (length - filled_length)
        print(f'\r{C.CYAN}{prefix} |{bar}| {percent}% {suffix}{C.END}', end='\r')
        if iteration == total:
            print()

def main():
    """Test the banner display with all modules"""
    banner = EcoliTyperBanner()
    
    # Display full banner with startup sequence
    banner.display_startup_sequence()
    banner.display_banner(show_quote=True, show_author=True)
    
    # Simulate ALL analysis modules with timing
    modules = [
        ("MLST Analysis", "mlst", "Multi-Locus Sequence Typing for E. coli"),
        ("Serotyping", "serotyping", "O and H antigen determination"),
        ("CH Typing", "ch_typing", "CH typing analysis"),
        ("Phylogrouping", "phylogrouping", "zClermont phylogrouping algorithm"),
        ("Abricate Analysis", "abricate", "Resistance, Virulence, and Plasmid gene screening"),
        ("AMRfinderPlus", "amrfinder", "NCBI AMR gene detection"),
        ("Lineage Database", "lineage_db", "Lineage database querying and reference")
    ]
    
    for module_name, module_key, description in modules:
        banner.start_analysis_timer(module_key)
        banner.display_module_header(module_name, description)
        # Simulate processing time
        time.sleep(0.5 + random.random())
        banner.stop_analysis_timer(module_key)
        banner.display_success(f"{module_name} completed!")
    
    # Test progress bar
    print("Overall progress demonstration:")
    for i in range(101):
        banner.display_progress_bar(i, 100, prefix='EcoliTyper Progress:', suffix='Complete', length=40)
        time.sleep(0.02)
    print()
    
    # Display footer with total time only
    banner.display_footer(samples_processed=8)
    
    # NEW: Display citation request and random footer
    banner.display_citation_request()
    banner.display_random_footer()

if __name__ == "__main__":
    main()