#!/usr/bin/env python3
"""
GENIUS COMPREHENSIVE E.COLI REPORTER - ULTIMATE EDITION
Advanced HTML Parser with Gene-Centric Cross-Genome Analysis
Author: Beckley Brown <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 1.0.0 - Ultimate Gene-Centric Edition-Department of Medical Biochemistry
Date: 2025-12-16
"""

import os
import sys
import json
import re
import glob
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any, Optional
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

# HTML parsing
from bs4 import BeautifulSoup

# Visualization
try:
    import plotly.graph_objects as go
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

class UltimateHTMLParser:
    """Ultimate HTML parser for all EcoliTyper reports"""
    
    def __init__(self):
        self.abricate_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ecoh', 'ecoli_vf'
        ]
    
    def normalize_sample_id(self, sample_id: str) -> str:
        """Normalize sample ID"""
        sample = str(sample_id)
        extensions = ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv']
        for ext in extensions:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        
        return sample.strip()
    
    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
        """Parse HTML table"""
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            
            if not tables or table_index >= len(tables):
                return pd.DataFrame()
            
            table = tables[table_index]
            rows = table.find_all('tr')
            
            headers = []
            for th in rows[0].find_all(['th', 'td']):
                headers.append(th.get_text().strip())
            
            data = []
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    if len(row_data) == len(headers):
                        data.append(row_data)
            
            if not data:
                return pd.DataFrame()
            
            return pd.DataFrame(data, columns=headers)
            
        except Exception as e:
            print(f"  ⚠️ Table parsing error: {e}")
            return pd.DataFrame()
    
    def parse_mlst_report(self, file_path: Path) -> Dict[str, Dict]:
        """Parse MLST HTML report"""
        print(f"  🧬 Parsing MLST: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            
            mlst_table = None
            for table in tables:
                if table.find(string=re.compile(r'Sample|ST|Allele', re.I)):
                    mlst_table = table
                    break
            
            if not mlst_table:
                mlst_table = soup.find('table')
            
            if not mlst_table:
                return {}
            
            data = []
            rows = mlst_table.find_all('tr')
            
            if len(rows) < 2:
                return {}
            
            headers = []
            header_cells = rows[0].find_all(['th', 'td'])
            for cell in header_cells:
                headers.append(cell.get_text().strip())
            
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    if len(row_data) >= 2:
                        data.append(row_data)
            
            if not data:
                return {}
            
            df = pd.DataFrame(data)
            if len(df.columns) > len(headers):
                df = df.iloc[:, :len(headers)]
            df.columns = headers[:len(df.columns)]
            
            df.columns = [col.strip() for col in df.columns]
            
            if 'Sample' in df.columns:
                df['normalized_sample'] = df['Sample'].apply(self.normalize_sample_id)
            elif 'sample' in df.columns.lower():
                sample_col = [col for col in df.columns if 'sample' in col.lower()][0]
                df['normalized_sample'] = df[sample_col].apply(self.normalize_sample_id)
            else:
                df['normalized_sample'] = df.iloc[:, 0].apply(self.normalize_sample_id)
            
            results = {}
            for _, row in df.iterrows():
                sample = row['normalized_sample']
                
                st = 'ND'
                if 'ST' in df.columns:
                    st = str(row['ST']) if pd.notna(row['ST']) else 'ND'
                elif 'st' in [col.lower() for col in df.columns]:
                    st_col = [col for col in df.columns if col.lower() == 'st'][0]
                    st = str(row[st_col]) if pd.notna(row[st_col]) else 'ND'
                
                allele_profile = 'ND'
                if 'Allele Profile' in df.columns:
                    allele_profile = str(row['Allele Profile']) if pd.notna(row['Allele Profile']) else 'ND'
                
                results[sample] = {
                    'ST': st,
                    'Allele_Profile': allele_profile
                }
            
            print(f"    ✓ Found {len(results)} samples")
            return results
            
        except Exception as e:
            print(f"    ❌ Error parsing MLST: {e}")
            return {}
    
    def parse_serotype_report(self, file_path: Path) -> Dict[str, Dict]:
        """Parse Serotype HTML report"""
        print(f"  🧬 Parsing Serotype: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            df = self.parse_html_table(html_content, 0)
            if df.empty:
                return {}
            
            df.columns = [col.strip() for col in df.columns]
            
            sample_col = None
            for col in df.columns:
                if 'sample' in col.lower() or 'id' in col.lower():
                    sample_col = col
                    break
            
            if not sample_col:
                sample_col = df.columns[0]
            
            df['normalized_sample'] = df[sample_col].apply(self.normalize_sample_id)
            
            results = {}
            for _, row in df.iterrows():
                sample = row['normalized_sample']
                
                serotype = 'ND'
                o_type = 'ND'
                h_type = 'ND'
                
                if 'Serotype' in df.columns:
                    serotype = str(row['Serotype']) if pd.notna(row['Serotype']) else 'ND'
                
                if 'O-type' in df.columns:
                    o_type = str(row['O-type']) if pd.notna(row['O-type']) else 'ND'
                elif 'O_Type' in df.columns:
                    o_type = str(row['O_Type']) if pd.notna(row['O_Type']) else 'ND'
                
                if 'H-type' in df.columns:
                    h_type = str(row['H-type']) if pd.notna(row['H-type']) else 'ND'
                elif 'H_Type' in df.columns:
                    h_type = str(row['H_Type']) if pd.notna(row['H_Type']) else 'ND'
                
                results[sample] = {
                    'Serotype': serotype,
                    'O_Type': o_type,
                    'H_Type': h_type
                }
            
            print(f"    ✓ Found {len(results)} samples")
            return results
            
        except Exception as e:
            print(f"    ❌ Error parsing Serotype: {e}")
            return {}
    
    def parse_chtyper_report(self, file_path: Path) -> Dict[str, Dict]:
        """Parse CHTyper HTML report"""
        print(f"  🧬 Parsing CHTyper: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            df = self.parse_html_table(html_content, 0)
            if df.empty:
                return {}
            
            df.columns = [col.strip() for col in df.columns]
            
            sample_col = None
            for col in df.columns:
                if 'sample' in col.lower() or 'id' in col.lower():
                    sample_col = col
                    break
            
            if not sample_col:
                sample_col = df.columns[0]
            
            df['normalized_sample'] = df[sample_col].apply(self.normalize_sample_id)
            
            results = {}
            for _, row in df.iterrows():
                sample = row['normalized_sample']
                
                fumc_type = 'ND'
                fimh_type = 'ND'
                
                if 'FumC Type' in df.columns:
                    fumc_type = str(row['FumC Type']) if pd.notna(row['FumC Type']) else 'ND'
                if 'FimH Type' in df.columns:
                    fimh_type = str(row['FimH Type']) if pd.notna(row['FimH Type']) else 'ND'
                
                ch_type = f"{fumc_type}:{fimh_type}"
                
                results[sample] = {
                    'FumC_Type': fumc_type,
                    'FimH_Type': fimh_type,
                    'CH_Type': ch_type
                }
            
            print(f"    ✓ Found {len(results)} samples")
            return results
            
        except Exception as e:
            print(f"    ❌ Error parsing CHTyper: {e}")
            return {}
    
    def parse_phylogrouping_report(self, file_path: Path) -> Dict[str, Dict]:
        """Parse Phylogrouping HTML report"""
        print(f"  🧬 Parsing Phylogrouping: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            df = self.parse_html_table(html_content, 0)
            if df.empty:
                return {}
            
            df.columns = [col.strip() for col in df.columns]
            
            sample_col = None
            for col in df.columns:
                if 'sample' in col.lower() or 'id' in col.lower():
                    sample_col = col
                    break
            
            if not sample_col:
                sample_col = df.columns[0]
            
            df['normalized_sample'] = df[sample_col].apply(self.normalize_sample_id)
            
            results = {}
            for _, row in df.iterrows():
                sample = row['normalized_sample']
                
                clermont_type = 'ND'
                if 'Clermont Type' in df.columns:
                    clermont_type = str(row['Clermont Type']) if pd.notna(row['Clermont Type']) else 'ND'
                elif 'Type' in df.columns:
                    clermont_type = str(row['Type']) if pd.notna(row['Type']) else 'ND'
                
                results[sample] = {
                    'Clermont_Type': clermont_type
                }
            
            print(f"    ✓ Found {len(results)} samples")
            return results
            
        except Exception as e:
            print(f"    ❌ Error parsing Phylogrouping: {e}")
            return {}
    
    def parse_amrfinder_report(self, file_path: Path) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        """Parse AMRfinder HTML report"""
        print(f"  🧬 Parsing AMRfinder: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            
            if len(tables) < 2:
                return {}, {}
            
            # Parse second table: Gene Frequency (contains genomes)
            gene_frequencies = {}
            df2 = self.parse_html_table(str(tables[1]), 0)
            if not df2.empty and 'Gene' in df2.columns:
                for _, row in df2.iterrows():
                    gene = str(row['Gene']).strip()
                    frequency = str(row.get('Frequency', '0')).strip()
                    
                    genomes = []
                    if 'Genomes' in df2.columns and pd.notna(row.get('Genomes')):
                        genomes_str = str(row['Genomes'])
                        genomes = [self.normalize_sample_id(g.strip()) 
                                  for g in genomes_str.split(',') if g.strip()]
                    
                    count = 0
                    match = re.search(r'(\d+)', frequency)
                    if match:
                        count = int(match.group(1))
                    
                    gene_frequencies[gene] = {
                        'frequency': frequency,
                        'count': count,
                        'genomes': genomes,
                        'database': 'amrfinder'
                    }
            
            # Parse first table: Genes by Genome (for reverse mapping)
            genes_by_genome = {}
            df1 = self.parse_html_table(str(tables[0]), 0)
            if not df1.empty:
                for _, row in df1.iterrows():
                    if 'Genome' not in df1.columns:
                        continue
                    sample = self.normalize_sample_id(row['Genome'])
                    genes = []
                    if pd.notna(row.get('Genes Detected')):
                        gene_str = str(row['Genes Detected'])
                        genes = [g.strip() for g in gene_str.split(',') if g.strip()]
                    genes_by_genome[sample] = genes
            
            print(f"    ✓ Found {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return genes_by_genome, gene_frequencies
            
        except Exception as e:
            print(f"    ❌ Error parsing AMRfinder: {e}")
            return {}, {}
    
    def parse_abricate_database_report(self, file_path: Path) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        """Parse ABRicate database HTML report"""
        print(f"  🧬 Parsing ABRicate: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            
            if len(tables) < 2:
                return {}, {}
            
            # Determine database name
            db_name = 'unknown'
            title_tag = soup.find('title')
            if title_tag:
                title_text = title_tag.get_text()
                match = re.search(r'(\w+)\s+Database', title_text, re.I)
                if match:
                    db_name = match.group(1).lower()
                else:
                    db_name = file_path.stem.lower().replace('ecoli_', '').replace('_summary_report', '')
            
            # Parse second table: Gene Frequency (contains genomes)
            gene_frequencies = {}
            df2 = self.parse_html_table(str(tables[1]), 0)
            if not df2.empty and 'Gene' in df2.columns:
                for _, row in df2.iterrows():
                    gene = str(row['Gene']).strip()
                    frequency = str(row.get('Frequency', '0')).strip()
                    
                    genomes = []
                    if 'Genomes' in df2.columns and pd.notna(row.get('Genomes')):
                        genomes_str = str(row['Genomes'])
                        genomes = [self.normalize_sample_id(g.strip()) 
                                  for g in genomes_str.split(',') if g.strip()]
                    
                    count = 0
                    match = re.search(r'(\d+)', frequency)
                    if match:
                        count = int(match.group(1))
                    
                    gene_frequencies[gene] = {
                        'frequency': frequency,
                        'count': count,
                        'genomes': genomes,
                        'database': db_name
                    }
            
            # Parse first table: Genes by Genome (for reverse mapping)
            genes_by_genome = {}
            df1 = self.parse_html_table(str(tables[0]), 0)
            if not df1.empty:
                sample_col = None
                for col in df1.columns:
                    if 'genome' in col.lower() or 'sample' in col.lower():
                        sample_col = col
                        break
                
                if sample_col:
                    for _, row in df1.iterrows():
                        sample = self.normalize_sample_id(row[sample_col])
                        genes = []
                        
                        genes_col = None
                        for col in df1.columns:
                            if 'genes' in col.lower() or 'detected' in col.lower():
                                genes_col = col
                                break
                        
                        if genes_col and pd.notna(row.get(genes_col)):
                            gene_str = str(row[genes_col])
                            genes = [g.strip() for g in gene_str.split(',') if g.strip()]
                        
                        genes_by_genome[sample] = genes
            
            print(f"    ✓ {db_name.upper()}: {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return genes_by_genome, gene_frequencies
            
        except Exception as e:
            print(f"    ❌ Error parsing ABRicate report: {e}")
            return {}, {}


class UltimateDataAnalyzer:
    """Analyzes data for ultimate gene-centric reporting"""
    
    def __init__(self):
        self.critical_amr_genes = {
            'blaCTX-M', 'blaSHV', 'blaTEM', 'blaKPC', 'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP',
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            'qnrA', 'qnrB', 'qnrC', 'qnrD', 'qnrS', 'aac(6\')-Ib-cr'
        }
        
        self.critical_virulence_genes = {
            'stx1', 'stx2', 'stx1A', 'stx1B', 'stx2A', 'stx2B',
            'cnf1', 'hlyA', 'hlyB', 'hlyC', 'hlyD', 'eae'
        }
    
    def create_gene_centric_tables(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        """Create gene-centric tables showing genes with their genomes"""
        gene_centric = {
            'amr_databases': {},
            'virulence_databases': {},
            'combined_gene_frequencies': []
        }
        
        # Process AMRfinder
        if 'amrfinder' in integrated_data.get('gene_frequencies', {}):
            amr_data = integrated_data['gene_frequencies']['amrfinder']
            gene_list = []
            
            for gene, data in amr_data.items():
                gene_list.append({
                    'gene': gene,
                    'database': 'AMRfinder',
                    'frequency': data.get('frequency', '0'),
                    'count': data.get('count', 0),
                    'genomes': data.get('genomes', [])
                })
            
            gene_centric['amr_databases']['amrfinder'] = sorted(gene_list, key=lambda x: x['count'], reverse=True)
        
        # Process ABRicate databases
        if 'abricate' in integrated_data.get('gene_frequencies', {}):
            abricate_data = integrated_data['gene_frequencies']['abricate']
            
            for db_name, db_genes in abricate_data.items():
                gene_list = []
                
                for gene, data in db_genes.items():
                    gene_list.append({
                        'gene': gene,
                        'database': db_name.upper(),
                        'frequency': data.get('frequency', '0'),
                        'count': data.get('count', 0),
                        'genomes': data.get('genomes', [])
                    })
                
                # Sort by count and store
                if gene_list:
                    gene_list.sort(key=lambda x: x['count'], reverse=True)
                    
                    if db_name in ['vfdb', 'ecoli_vf']:
                        gene_centric['virulence_databases'][db_name] = gene_list
                    else:
                        gene_centric['amr_databases'][db_name] = gene_list
        
        # Create combined gene frequencies for pattern discovery
        all_genes = []
        
        for db_type in ['amr_databases', 'virulence_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    all_genes.append(gene_data)
        
        # Sort combined list by count
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        gene_centric['combined_gene_frequencies'] = all_genes
        
        return gene_centric
    
    def create_cross_genome_patterns(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        """Create cross-genome patterns"""
        patterns = {
            'st_distribution': Counter(),
            'serotype_distribution': Counter(),
            'phylogroup_distribution': Counter(),
            'ch_type_distribution': Counter(),
            'st_serotype_combinations': defaultdict(list),
            'gene_cooccurrence': defaultdict(Counter),
            'high_risk_combinations': []
        }
        
        samples_data = integrated_data.get('samples', {})
        gene_centric = integrated_data.get('gene_centric', {})
        
        # Collect all genes per sample for co-occurrence
        sample_genes = defaultdict(list)
        for db_type in ['amr_databases', 'virulence_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    for genome in gene_data['genomes']:
                        if gene_data['gene'] not in sample_genes[genome]:
                            sample_genes[genome].append(gene_data['gene'])
        
        # Analyze each sample
        for sample, data in samples_data.items():
            st = data.get('mlst', {}).get('ST', 'ND')
            serotype = data.get('serotype', {}).get('Serotype', 'ND')
            phylogroup = data.get('phylogrouping', {}).get('Clermont_Type', 'ND')
            ch_type = data.get('chtyper', {}).get('CH_Type', 'ND')
            
            # Basic distributions
            if st != 'ND':
                patterns['st_distribution'][st] += 1
            if serotype != 'ND':
                patterns['serotype_distribution'][serotype] += 1
            if phylogroup != 'ND':
                patterns['phylogroup_distribution'][phylogroup] += 1
            if ch_type != 'ND':
                patterns['ch_type_distribution'][ch_type] += 1
            
            # ST-Serotype combinations
            if st != 'ND' and serotype != 'ND':
                patterns['st_serotype_combinations'][f"ST{st} - {serotype}"].append(sample)
            
            # Gene co-occurrence
            genes = sample_genes.get(sample, [])
            for i, gene1 in enumerate(genes):
                for gene2 in genes[i+1:]:
                    patterns['gene_cooccurrence'][gene1][gene2] += 1
            
            # High-risk combinations
            amr_genes = data.get('amr_genes', [])
            virulence_genes = data.get('virulence_genes', [])
            
            critical_amr = [g for g in amr_genes if any(crit in str(g).lower() for crit in self.critical_amr_genes)]
            critical_vir = [g for g in virulence_genes if any(crit in str(g).lower() for crit in self.critical_virulence_genes)]
            
            if critical_amr and critical_vir:
                patterns['high_risk_combinations'].append({
                    'sample': sample,
                    'st': st,
                    'serotype': serotype,
                    'critical_amr_genes': critical_amr,
                    'critical_virulence_genes': critical_vir
                })
        
        return patterns


class UltimateHTMLGenerator:
    """Generates ultimate HTML reports with gene-centric approach"""
    
    def __init__(self, data_analyzer: UltimateDataAnalyzer):
        self.data_analyzer = data_analyzer
        self.tab_colors = {
            'summary': '#4CAF50',
            'samples': '#2196F3',
            'mlst': '#FF9800',
            'serotype': '#9C27B0',
            'chtype': '#009688',
            'phylogroup': '#795548',
            'amr': '#F44336',
            'virulence': '#E91E63',
            'patterns': '#FF5722',
            'export': '#3F51B5'
        }
    
    def generate_main_report(self, integrated_data: Dict[str, Any], output_dir: Path) -> str:
        """Generate the ultimate HTML report"""
        print("\n🎨 Generating ULTIMATE HTML report...")
        
        # Extract data
        samples_data = integrated_data.get('samples', {})
        patterns = integrated_data.get('patterns', {})
        gene_centric = integrated_data.get('gene_centric', {})
        metadata = integrated_data.get('metadata', {})
        
        # Create HTML
        html = self._create_ultimate_html(
            metadata=metadata,
            samples_data=samples_data,
            patterns=patterns,
            gene_centric=gene_centric
        )
        
        # Save HTML file
        output_file = output_dir / "genius_ultimate_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)
    
    def _create_ultimate_html(self, **kwargs) -> str:
        """Create ultimate HTML with all sections"""
        
        # CSS Styles - UPDATED TO REMOVE FIXED WIDTHS AND ALLOW FULL EXPANSION
        css = """
        <style>
        :root {
            --summary-color: #4CAF50;
            --samples-color: #2196F3;
            --mlst-color: #FF9800;
            --serotype-color: #9C27B0;
            --chtype-color: #009688;
            --phylogroup-color: #795548;
            --amr-color: #F44336;
            --virulence-color: #E91E63;
            --patterns-color: #FF5722;
            --export-color: #3F51B5;
        }
        
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background: #f5f5f5;
            min-width: 1200px; /* Ensure minimum width for tables */
        }
        
        .container {
            max-width: none; /* REMOVED max-width constraint */
            margin: 0 auto;
            padding: 20px;
            width: 100%;
            overflow-x: auto; /* Allow horizontal scrolling if needed */
        }
        
        .main-header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
            margin-bottom: 30px;
            text-align: center;
        }
        
        .main-header h1 {
            font-size: 2.8em;
            margin-bottom: 10px;
            color: white;
        }
        
        .metadata-bar {
            background: rgba(255,255,255,0.1);
            padding: 15px;
            border-radius: 10px;
            margin: 20px 0;
            display: flex;
            justify-content: space-around;
            flex-wrap: wrap;
            gap: 15px;
            backdrop-filter: blur(10px);
        }
        
        .metadata-item {
            display: flex;
            align-items: center;
            gap: 8px;
            font-size: 0.95em;
        }
        
        .dashboard-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }
        
        .dashboard-card {
            background: white;
            padding: 25px;
            border-radius: 12px;
            box-shadow: 0 5px 20px rgba(0,0,0,0.1);
            text-align: center;
            transition: all 0.3s ease;
            cursor: pointer;
            border-left: 5px solid;
            position: relative;
            overflow: hidden;
        }
        
        .dashboard-card::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 3px;
            background: linear-gradient(90deg, transparent, rgba(255,255,255,0.8), transparent);
        }
        
        .dashboard-card:hover {
            transform: translateY(-10px);
            box-shadow: 0 15px 30px rgba(0,0,0,0.2);
        }
        
        .card-summary { border-left-color: var(--summary-color); }
        .card-samples { border-left-color: var(--samples-color); }
        .card-mlst { border-left-color: var(--mlst-color); }
        .card-serotype { border-left-color: var(--serotype-color); }
        .card-chtype { border-left-color: var(--chtype-color); }
        .card-phylogroup { border-left-color: var(--phylogroup-color); }
        .card-amr { border-left-color: var(--amr-color); }
        .card-virulence { border-left-color: var(--virulence-color); }
        .card-patterns { border-left-color: var(--patterns-color); }
        .card-export { border-left-color: var(--export-color); }
        
        .card-number {
            font-size: 3em;
            font-weight: bold;
            margin: 15px 0;
            background: linear-gradient(90deg, #667eea, #764ba2);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        
        .tab-navigation {
            display: flex;
            gap: 5px;
            margin-bottom: 20px;
            flex-wrap: wrap;
            background: white;
            padding: 15px;
            border-radius: 12px;
            box-shadow: 0 5px 20px rgba(0,0,0,0.1);
            position: sticky;
            top: 10px;
            z-index: 100;
        }
        
        .tab-button {
            padding: 12px 25px;
            background: #f5f5f5;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            font-weight: 600;
            color: #666;
            transition: all 0.3s ease;
            display: flex;
            align-items: center;
            gap: 8px;
            position: relative;
            overflow: hidden;
        }
        
        .tab-button::after {
            content: '';
            position: absolute;
            bottom: 0;
            left: 50%;
            right: 50%;
            height: 3px;
            background: currentColor;
            transition: all 0.3s ease;
        }
        
        .tab-button:hover::after {
            left: 10%;
            right: 10%;
        }
        
        .tab-button.active {
            color: white;
        }
        
        .tab-button.active::after {
            left: 10%;
            right: 10%;
        }
        
        .tab-button.summary.active { background: var(--summary-color); }
        .tab-button.samples.active { background: var(--samples-color); }
        .tab-button.mlst.active { background: var(--mlst-color); }
        .tab-button.serotype.active { background: var(--serotype-color); }
        .tab-button.chtype.active { background: var(--chtype-color); }
        .tab-button.phylogroup.active { background: var(--phylogroup-color); }
        .tab-button.amr.active { background: var(--amr-color); }
        .tab-button.virulence.active { background: var(--virulence-color); }
        .tab-button.patterns.active { background: var(--patterns-color); }
        .tab-button.export.active { background: var(--export-color); }
        
        .tab-content {
            display: none;
            background: white;
            padding: 30px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
            margin-bottom: 30px;
            animation: fadeIn 0.5s ease;
            width: 100%;
            overflow-x: auto; /* Allow horizontal scrolling for tables */
        }
        
        .tab-content.active {
            display: block;
        }
        
        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(20px); }
            to { opacity: 1; transform: translateY(0); }
        }
        
        .section-header {
            color: #2c3e50;
            margin-bottom: 25px;
            padding-bottom: 15px;
            border-bottom: 3px solid;
            font-size: 1.8em;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        
        .summary-header { border-color: var(--summary-color); }
        .samples-header { border-color: var(--samples-color); }
        .mlst-header { border-color: var(--mlst-color); }
        .serotype-header { border-color: var(--serotype-color); }
        .chtype-header { border-color: var(--chtype-color); }
        .phylogroup-header { border-color: var(--phylogroup-color); }
        .amr-header { border-color: var(--amr-color); }
        .virulence-header { border-color: var(--virulence-color); }
        .patterns-header { border-color: var(--patterns-color); }
        .export-header { border-color: var(--export-color); }
        
        .data-table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 0.95em;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            border-radius: 8px;
            overflow: hidden;
            table-layout: auto; /* Allow columns to expand based on content */
        }
        
        .data-table th {
            background: #2c3e50;
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
            position: sticky;
            top: 0;
            white-space: nowrap;
        }
        
        .data-table td {
            padding: 12px;
            border-bottom: 1px solid #e0e0e0;
            vertical-align: top;
            word-wrap: break-word;
            word-break: break-word;
            white-space: normal; /* Allow text wrapping */
        }
        
        .data-table tr:hover {
            background: #f8f9fa;
        }
        
        .scrollable-table {
            max-height: none; /* REMOVED height restriction */
            overflow-y: auto;
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            margin: 20px 0;
            width: 100%;
        }
        
        .search-box {
            width: 100%;
            padding: 12px;
            margin-bottom: 20px;
            border: 2px solid #e0e0e0;
            border-radius: 8px;
            font-size: 1em;
            transition: all 0.3s ease;
        }
        
        .search-box:focus {
            outline: none;
            border-color: #667eea;
            box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.1);
        }
        
        .badge {
            display: inline-block;
            padding: 5px 15px;
            border-radius: 20px;
            font-size: 0.85em;
            font-weight: 600;
            margin: 2px;
        }
        
        .badge-low { background: #4CAF50; color: white; }
        .badge-medium { background: #FF9800; color: black; }
        .badge-high { background: #F44336; color: white; }
        .badge-critical { background: #9C27B0; color: white; }
        
        .alert-box {
            padding: 20px;
            border-radius: 10px;
            margin: 20px 0;
            display: flex;
            align-items: center;
            gap: 20px;
            border-left: 5px solid;
        }
        
        .alert-success { background: #d4edda; color: #155724; border-left-color: #28a745; }
        .alert-warning { background: #fff3cd; color: #856404; border-left-color: #ffc107; }
        .alert-danger { background: #f8d7da; color: #721c24; border-left-color: #dc3545; }
        .alert-info { background: #d1ecf1; color: #0c5460; border-left-color: #17a2b8; }
        
        .action-buttons {
            display: flex;
            gap: 10px;
            margin: 20px 0;
            flex-wrap: wrap;
        }
        
        .action-btn {
            padding: 10px 20px;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            font-weight: 600;
            display: flex;
            align-items: center;
            gap: 8px;
            transition: all 0.3s ease;
        }
        
        .action-btn:hover {
            transform: translateY(-2px);
            box-shadow: 0 5px 15px rgba(0,0,0,0.2);
        }
        
        .btn-primary { background: #667eea; color: white; }
        .btn-success { background: #28a745; color: white; }
        .btn-danger { background: #dc3545; color: white; }
        .btn-warning { background: #ffc107; color: black; }
        
        .database-section {
            margin: 30px 0;
            padding: 25px;
            border-radius: 12px;
            background: #f8f9fa;
            box-shadow: 0 3px 15px rgba(0,0,0,0.08);
        }
        
        .database-header {
            font-size: 1.4em;
            color: #2c3e50;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid #667eea;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        
        .print-section-btn {
            background: #667eea;
            color: white;
            border: none;
            border-radius: 5px;
            padding: 8px 15px;
            cursor: pointer;
            display: flex;
            align-items: center;
            gap: 5px;
            font-size: 0.9em;
        }
        
        .print-section-btn:hover {
            background: #764ba2;
        }
        
        .genome-list {
            display: flex;
            flex-wrap: wrap;
            gap: 5px;
            margin-top: 5px;
            max-width: none; /* REMOVED width restriction */
        }
        
        .genome-tag {
            background: #e3f2fd;
            color: #1976d2;
            padding: 3px 10px;
            border-radius: 12px;
            font-size: 0.85em;
            border: 1px solid #bbdefb;
            white-space: nowrap;
            margin: 2px;
        }
        
        .footer {
            text-align: center;
            padding: 30px;
            color: white;
            margin-top: 40px;
            border-radius: 15px;
            background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
        }
        
        /* REMOVED fixed width constraints from tables */
        .data-table td:nth-child(5), /* Genomes column */
        .data-table td:nth-child(4), /* Samples column */
        .data-table td:nth-child(3) { /* Common STs column */
            min-width: 300px; /* Minimum width but can expand */
            max-width: none; /* No maximum width */
        }
        
        @media print {
            body * {
                visibility: hidden;
            }
            .tab-content.active,
            .tab-content.active * {
                visibility: visible;
            }
            .tab-content.active {
                position: absolute;
                left: 0;
                top: 0;
                width: 100%;
                padding: 20px;
                box-shadow: none;
                border-radius: 0;
            }
            .print-section-btn,
            .tab-navigation,
            .dashboard-grid,
            .search-box,
            .action-buttons {
                display: none !important;
            }
            
            /* Ensure tables print properly */
            .data-table {
                page-break-inside: auto;
            }
            .data-table tr {
                page-break-inside: avoid;
                page-break-after: auto;
            }
            .data-table td, .data-table th {
                page-break-inside: avoid;
            }
        }
        
        @media (max-width: 768px) {
            .container {
                padding: 10px;
            }
            
            .main-header {
                padding: 20px;
            }
            
            .main-header h1 {
                font-size: 2em;
            }
            
            .tab-button {
                padding: 10px 15px;
                font-size: 0.9em;
            }
            
            .dashboard-grid {
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            }
            
            .data-table {
                font-size: 0.85em;
            }
            
            body {
                min-width: auto;
                overflow-x: auto;
            }
        }
        </style>
        """
        
        # JavaScript
        js = """
        <script>
        // Tab switching
        function switchTab(tabName) {
            // Hide all tabs
            document.querySelectorAll('.tab-content').forEach(tab => {
                tab.classList.remove('active');
            });
            
            // Remove active class from all buttons
            document.querySelectorAll('.tab-button').forEach(button => {
                button.classList.remove('active');
            });
            
            // Show selected tab
            document.getElementById(tabName + '-tab').classList.add('active');
            
            // Activate selected button
            event.currentTarget.classList.add('active');
            
            // Update URL hash
            window.location.hash = tabName;
        }
        
        // Search functionality
        function searchTable(tableId, searchId) {
            const input = document.getElementById(searchId);
            const filter = input.value.toUpperCase();
            const table = document.getElementById(tableId);
            const rows = table.getElementsByTagName('tr');
            
            for (let i = 1; i < rows.length; i++) {
                const cells = rows[i].getElementsByTagName('td');
                let found = false;
                
                for (let j = 0; j < cells.length; j++) {
                    const cell = cells[j];
                    if (cell) {
                        const txtValue = cell.textContent || cell.innerText;
                        if (txtValue.toUpperCase().indexOf(filter) > -1) {
                            found = true;
                            break;
                        }
                    }
                }
                
                rows[i].style.display = found ? '' : 'none';
            }
        }
        
        // Print current section
        function printSection(sectionId) {
            const content = document.getElementById(sectionId);
            const printWindow = window.open('', '_blank');
            printWindow.document.write('<html><head><title>Print Section</title>');
            printWindow.document.write('<style>' + document.querySelector('style').textContent + '</style>');
            printWindow.document.write('</head><body>');
            printWindow.document.write(content.innerHTML);
            printWindow.document.write('</body></html>');
            printWindow.document.close();
            printWindow.print();
        }
        
        // Export table to CSV
        function exportTableToCSV(tableId, filename) {
            const table = document.getElementById(tableId);
            const rows = table.querySelectorAll('tr');
            const csv = [];
            
            for (let i = 0; i < rows.length; i++) {
                const row = [], cols = rows[i].querySelectorAll('td, th');
                
                for (let j = 0; j < cols.length; j++) {
                    row.push('"' + (cols[j].innerText || '').replace(/"/g, '""') + '"');
                }
                
                csv.push(row.join(','));
            }
            
            const csvFile = new Blob([csv.join('\\n')], {type: 'text/csv'});
            const downloadLink = document.createElement('a');
            downloadLink.download = filename;
            downloadLink.href = window.URL.createObjectURL(csvFile);
            downloadLink.style.display = 'none';
            document.body.appendChild(downloadLink);
            downloadLink.click();
            document.body.removeChild(downloadLink);
        }
        
        // Initialize from URL hash
        document.addEventListener('DOMContentLoaded', function() {
            const hash = window.location.hash.substring(1);
            if (hash) {
                const tabButton = document.querySelector(`.tab-button.${hash}`);
                if (tabButton) {
                    tabButton.click();
                }
            } else {
                // Show first tab
                document.querySelector('.tab-button').click();
            }
        });
        </script>
        """
        
        # Build HTML
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GENIUS E.coli Ultimate Report</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    {css}
    {js}
</head>
<body>
    <div class="container">
        <!-- Main Header -->
        <div class="main-header">
            <h1><i class="fas fa-dna"></i> GENIUS E.coli Ultimate Analysis Report</h1>
            <p>Gene-Centric Cross-Genome Analysis with Complete Genome Lists</p>
            
            <div class="metadata-bar">
                <div class="metadata-item">
                    <i class="fas fa-calendar"></i>
                    <span>Generated: {kwargs['metadata'].get('analysis_date', 'Unknown')}</span>
                </div>
                <div class="metadata-item">
                    <i class="fas fa-database"></i>
                    <span>Samples: {len(kwargs['samples_data'])}</span>
                </div>
                <div class="metadata-item">
                    <i class="fas fa-user-md"></i>
                    <span>Tool: GENIUS Ultimate v1.0.0</span>
                </div>
                <div class="metadata-item">
                    <i class="fas fa-university"></i>
                    <span>University of Ghana Medical School</span>
                </div>
            </div>
        </div>
        
        <!-- Dashboard -->
        <div class="dashboard-grid">
            <div class="dashboard-card card-summary" onclick="switchTab('summary')">
                <div class="card-number">{len(kwargs['samples_data'])}</div>
                <div class="card-label">Total Samples</div>
                <i class="fas fa-vial fa-2x" style="color: var(--summary-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-mlst" onclick="switchTab('mlst')">
                <div class="card-number">{len(kwargs['patterns'].get('st_distribution', {}))}</div>
                <div class="card-label">Unique STs</div>
                <i class="fas fa-code-branch fa-2x" style="color: var(--mlst-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-serotype" onclick="switchTab('serotype')">
                <div class="card-number">{len(kwargs['patterns'].get('serotype_distribution', {}))}</div>
                <div class="card-label">Serotypes</div>
                <i class="fas fa-tag fa-2x" style="color: var(--serotype-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-amr" onclick="switchTab('amr')">
                <div class="card-number">{len(kwargs['gene_centric'].get('amr_databases', {}).get('amrfinder', [])) if kwargs['gene_centric'].get('amr_databases', {}).get('amrfinder') else 0}</div>
                <div class="card-label">AMR Genes</div>
                <i class="fas fa-biohazard fa-2x" style="color: var(--amr-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-virulence" onclick="switchTab('virulence')">
                <div class="card-number">{sum(len(genes) for genes in kwargs['gene_centric'].get('virulence_databases', {}).values())}</div>
                <div class="card-label">Virulence Genes</div>
                <i class="fas fa-virus fa-2x" style="color: var(--virulence-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-patterns" onclick="switchTab('patterns')">
                <div class="card-number">{len(kwargs['patterns'].get('high_risk_combinations', []))}</div>
                <div class="card-label">High-Risk Combos</div>
                <i class="fas fa-project-diagram fa-2x" style="color: var(--patterns-color); margin-top: 10px;"></i>
            </div>
        </div>
        
        <!-- Tab Navigation -->
        <div class="tab-navigation">
            <button class="tab-button summary active" onclick="switchTab('summary')">
                <i class="fas fa-chart-pie"></i> Summary
            </button>
            <button class="tab-button samples" onclick="switchTab('samples')">
                <i class="fas fa-list-alt"></i> Sample Overview
            </button>
            <button class="tab-button mlst" onclick="switchTab('mlst')">
                <i class="fas fa-code-branch"></i> MLST Analysis
            </button>
            <button class="tab-button serotype" onclick="switchTab('serotype')">
                <i class="fas fa-tag"></i> Serotype Analysis
            </button>
            <button class="tab-button chtype" onclick="switchTab('chtype')">
                <i class="fas fa-project-diagram"></i> CH Type
            </button>
            <button class="tab-button phylogroup" onclick="switchTab('phylogroup')">
                <i class="fas fa-sitemap"></i> Phylogrouping
            </button>
            <button class="tab-button amr" onclick="switchTab('amr')">
                <i class="fas fa-biohazard"></i> AMR Genes
            </button>
            <button class="tab-button virulence" onclick="switchTab('virulence')">
                <i class="fas fa-virus"></i> Virulence Genes
            </button>
            <button class="tab-button patterns" onclick="switchTab('patterns')">
                <i class="fas fa-project-diagram"></i> Pattern Discovery
            </button>
            <button class="tab-button export" onclick="switchTab('export')">
                <i class="fas fa-download"></i> Export
            </button>
        </div>
        
        <!-- Summary Tab -->
        <div id="summary-tab" class="tab-content active">
            <h2 class="section-header summary-header">
                <i class="fas fa-chart-pie"></i> Executive Summary
                <button class="print-section-btn" onclick="printSection('summary-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_summary_section(kwargs)}
        </div>
        
        <!-- Sample Overview Tab -->
        <div id="samples-tab" class="tab-content">
            <h2 class="section-header samples-header">
                <i class="fas fa-list-alt"></i> Complete Sample Overview
                <button class="print-section-btn" onclick="printSection('samples-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_sample_overview_section(kwargs)}
        </div>
        
        <!-- MLST Analysis Tab -->
        <div id="mlst-tab" class="tab-content">
            <h2 class="section-header mlst-header">
                <i class="fas fa-code-branch"></i> MLST Analysis
                <button class="print-section-btn" onclick="printSection('mlst-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_mlst_section(kwargs)}
        </div>
        
        <!-- Serotype Analysis Tab -->
        <div id="serotype-tab" class="tab-content">
            <h2 class="section-header serotype-header">
                <i class="fas fa-tag"></i> Serotype Analysis
                <button class="print-section-btn" onclick="printSection('serotype-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_serotype_section(kwargs)}
        </div>
        
        <!-- CH Type Tab -->
        <div id="chtype-tab" class="tab-content">
            <h2 class="section-header chtype-header">
                <i class="fas fa-project-diagram"></i> CH Type Analysis (FumC/FimH)
                <button class="print-section-btn" onclick="printSection('chtype-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_chtype_section(kwargs)}
        </div>
        
        <!-- Phylogrouping Tab -->
        <div id="phylogroup-tab" class="tab-content">
            <h2 class="section-header phylogroup-header">
                <i class="fas fa-sitemap"></i> Phylogrouping Analysis
                <button class="print-section-btn" onclick="printSection('phylogroup-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_phylogroup_section(kwargs)}
        </div>
        
        <!-- AMR Genes Tab -->
        <div id="amr-tab" class="tab-content">
            <h2 class="section-header amr-header">
                <i class="fas fa-biohazard"></i> Antimicrobial Resistance Genes
                <button class="print-section-btn" onclick="printSection('amr-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_amr_section(kwargs)}
        </div>
        
        <!-- Virulence Genes Tab -->
        <div id="virulence-tab" class="tab-content">
            <h2 class="section-header virulence-header">
                <i class="fas fa-virus"></i> Virulence Genes
                <button class="print-section-btn" onclick="printSection('virulence-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_virulence_section(kwargs)}
        </div>
        
        <!-- Pattern Discovery Tab -->
        <div id="patterns-tab" class="tab-content">
            <h2 class="section-header patterns-header">
                <i class="fas fa-project-diagram"></i> Cross-Genome Pattern Discovery
                <button class="print-section-btn" onclick="printSection('patterns-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_pattern_discovery_section(kwargs)}
        </div>
        
        <!-- Export Tab -->
        <div id="export-tab" class="tab-content">
            <h2 class="section-header export-header">
                <i class="fas fa-download"></i> Export Data
                <button class="print-section-btn" onclick="printSection('export-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_export_section(kwargs)}
        </div>
        
        <!-- Footer -->
        <div class="footer">
            <h3>GENIUS E.coli Ultimate Reporter v1.0.0</h3>
            <p>University of Ghana Medical School | Brown Beckley <brownbeckley94@gmail.com></p>
            <p>Generated on {kwargs['metadata'].get('analysis_date', 'Unknown')}</p>
        </div>
    </div>
</body>
</html>
        """
        
        return html
    
    def _generate_summary_section(self, kwargs: Dict) -> str:
        """Generate summary section"""
        samples_data = kwargs['samples_data']
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        
        total_samples = len(samples_data)
        total_amr_genes = sum(len(genes) for genes in gene_centric.get('amr_databases', {}).values())
        total_virulence_genes = sum(len(genes) for genes in gene_centric.get('virulence_databases', {}).values())
        critical_findings = len(patterns.get('high_risk_combinations', []))
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Analysis Overview</h3>
                <p>This ultimate gene-centric report analyzes <strong>{total_samples}</strong> E. coli genomes. 
                Instead of listing genes per genome, we show each gene with all genomes that contain it - 
                making it easy to track gene distribution across samples.</p>
            </div>
        </div>
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="switchTab('amr')">
                <i class="fas fa-biohazard"></i> View AMR Genes
            </button>
            <button class="action-btn btn-success" onclick="switchTab('virulence')">
                <i class="fas fa-virus"></i> View Virulence Genes
            </button>
            <button class="action-btn btn-danger" onclick="switchTab('patterns')">
                <i class="fas fa-exclamation-triangle"></i> Check High-Risk Combos
            </button>
        </div>
        
        <h3><i class="fas fa-chart-bar"></i> Key Statistics</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Metric</th>
                        <th>Count</th>
                        <th>Details</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td>Total Samples Analyzed</td>
                        <td><strong>{total_samples}</strong></td>
                        <td>Complete genomic analysis with all databases</td>
                    </tr>
                    <tr>
                        <td>Unique Sequence Types</td>
                        <td><strong>{len(patterns.get('st_distribution', {}))}</strong></td>
                        <td>MLST typing results</td>
                    </tr>
                    <tr>
                        <td>Unique Serotypes</td>
                        <td><strong>{len(patterns.get('serotype_distribution', {}))}</strong></td>
                        <td>O and H antigen typing</td>
                    </tr>
                    <tr>
                        <td>AMR Genes Identified</td>
                        <td><strong>{total_amr_genes}</strong></td>
                        <td>Across all AMR databases (AMRfinder, CARD, ResFinder, etc.)</td>
                    </tr>
                    <tr>
                        <td>Virulence Genes Identified</td>
                        <td><strong>{total_virulence_genes}</strong></td>
                        <td>From virulence databases (VFDB, Ecoli_VF)</td>
                    </tr>
                    <tr>
                        <td>High-Risk Combinations</td>
                        <td><span class="badge {'badge-critical' if critical_findings > 0 else 'badge-low'}">{critical_findings}</span></td>
                        <td>Samples with both critical AMR and virulence genes</td>
                    </tr>
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-lightbulb"></i> Report Features</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">
            <div class="database-section">
                <h4><i class="fas fa-gene"></i> Gene-Centric Approach</h4>
                <p>Each gene is shown with ALL genomes that contain it. No more searching through sample lists!</p>
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-print"></i> Section-Specific Printing</h4>
                <p>Print button on each section header prints only that section. Report any issues by quick mail!</p>
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-search"></i> Comprehensive Search</h4>
                <p>Search any table by gene name, sample name, ST, serotype, or any field.</p>
            </div>
        </div>
        """
        
        return html
    
    def _generate_sample_overview_section(self, kwargs: Dict) -> str:
        """Generate complete sample overview"""
        samples_data = kwargs['samples_data']
        
        html = f"""
        <input type="text" class="search-box" id="search-samples" 
               onkeyup="searchTable('samples-table', 'search-samples')" 
               placeholder="🔍 Search samples by any field...">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table', 'sample_overview.csv')">
                <i class="fas fa-download"></i> Export to CSV
            </button>
            <button class="action-btn btn-success" onclick="document.getElementById('search-samples').value=''; searchTable('samples-table', 'search-samples')">
                <i class="fas fa-sync"></i> Clear Search
            </button>
        </div>
        
        <div class="scrollable-table">
            <table id="samples-table" class="data-table">
                <thead>
                    <tr>
                        <th>Sample</th>
                        <th>ST</th>
                        <th>Serotype</th>
                        <th>Phylogroup</th>
                        <th>CH Type</th>
                        <th>Virulence Gene Count</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for sample, data in samples_data.items():
            st = data.get('mlst', {}).get('ST', 'ND')
            serotype = data.get('serotype', {}).get('Serotype', 'ND')
            phylogroup = data.get('phylogrouping', {}).get('Clermont_Type', 'ND')
            ch_type = data.get('chtyper', {}).get('CH_Type', 'ND')
            amr_count = len(data.get('amr_genes', []))
            virulence_count = len(data.get('virulence_genes', []))
            
            html += f"""
                    <tr>
                        <td><strong>{sample}</strong></td>
                        <td>{st}</td>
                        <td>{serotype}</td>
                        <td>{phylogroup}</td>
                        <td>{ch_type}</td>
                        <td>{virulence_count}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_mlst_section(self, kwargs: Dict) -> str:
        """Generate MLST analysis section"""
        patterns = kwargs['patterns']
        st_dist = patterns.get('st_distribution', Counter())
        st_serotype_combos = patterns.get('st_serotype_combinations', {})
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>MLST Analysis</h3>
                <p><strong>{len(st_dist)} unique sequence types</strong> identified. Each ST is shown with its associated serotypes and sample counts.</p>
            </div>
        </div>
        
        <h3><i class="fas fa-chart-bar"></i> Sequence Type Distribution</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Sequence Type</th>
                        <th>Count</th>
                        <th>Percentage</th>
                        <th>Associated Serotypes</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        total = sum(st_dist.values())
        for st, count in st_dist.most_common():
            percentage = (count / total) * 100
            
            # Get associated serotypes for this ST
            associated_serotypes = []
            for combo, samples in st_serotype_combos.items():
                if f"ST{st} - " in combo:
                    serotype = combo.split(" - ")[1]
                    associated_serotypes.append(serotype)
            
            if associated_serotypes:
                unique_serotypes = list(set(associated_serotypes))
                serotype_list = ', '.join(unique_serotypes) if unique_serotypes else 'ND'
            else:
                serotype_list = 'ND'            
            
            
            html += f"""
                    <tr>
                        <td><strong>{st if st.startswith('ST') else f'ST{st}'}</strong></td>
                        <td>{count}</td>
                        <td>{percentage:.1f}%</td>
                        <td>{serotype_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-project-diagram"></i> ST-Serotype Combinations</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>ST-Serotype Combination</th>
                        <th>Samples</th>
                        <th>Count</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for combo, samples in sorted(st_serotype_combos.items(), key=lambda x: len(x[1]), reverse=True):
            sample_list = ', '.join(samples) if samples else 'None'
            
            html += f"""
                    <tr>
                        <td><strong>{combo}</strong></td>
                        <td>{sample_list}</td>
                        <td>{len(samples)}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_serotype_section(self, kwargs: Dict) -> str:
        """Generate serotype analysis section"""
        patterns = kwargs['patterns']
        serotype_dist = patterns.get('serotype_distribution', Counter())
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Serotype Analysis</h3>
                <p><strong>{len(serotype_dist)} unique serotypes</strong> identified across all samples.</p>
            </div>
        </div>
        
        <h3><i class="fas fa-chart-bar"></i> Serotype Distribution</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Serotype</th>
                        <th>Count</th>
                        <th>Percentage</th>
                        <th>Common STs</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        total = sum(serotype_dist.values())
        for serotype, count in serotype_dist.most_common():
            if serotype == 'ND':
                continue
                
            percentage = (count / total) * 100
            
            # Find STs with this serotype
            sts = []
            st_serotype_combos = patterns.get('st_serotype_combinations', {})
            for combo, samples in st_serotype_combos.items():
                if serotype in combo:
                    st = combo.split(" - ")[0].replace("ST", "")
                    if st not in sts:
                        sts.append(st)
            
            st_list = ', '.join([f"ST{s}" for s in sts]) if sts else 'Various'
            
            html += f"""
                    <tr>
                        <td><strong>{serotype}</strong></td>
                        <td>{count}</td>
                        <td>{percentage:.1f}%</td>
                        <td>{st_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_chtype_section(self, kwargs: Dict) -> str:
        """Generate CH type analysis section"""
        patterns = kwargs['patterns']
        ch_dist = patterns.get('ch_type_distribution', Counter())
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>CH Type (FumC/FimH) Analysis</h3>
                <p><strong>{len(ch_dist)} unique CH types</strong> identified through FumC and FimH allele typing.</p>
            </div>
        </div>
        
        <h3><i class="fas fa-chart-bar"></i> CH Type Distribution</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>CH Type (FumC:FimH)</th>
                        <th>Count</th>
                        <th>Percentage</th>
                        <th>Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        total = sum(ch_dist.values())
        for ch_type, count in ch_dist.most_common():
            if ch_type == 'ND' or 'ND:' in ch_type or ':ND' in ch_type:
                continue
                
            percentage = (count / total) * 100
            
            # Find samples with this CH type
            samples_with_chtype = []
            for sample, data in kwargs['samples_data'].items():
                if data.get('chtyper', {}).get('CH_Type') == ch_type:
                    samples_with_chtype.append(sample)
            
            sample_list = ', '.join(samples_with_chtype) if samples_with_chtype else 'None'
            
            html += f"""
                    <tr>
                        <td><strong>{ch_type}</strong></td>
                        <td>{count}</td>
                        <td>{percentage:.1f}%</td>
                        <td>{sample_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_phylogroup_section(self, kwargs: Dict) -> str:
        """Generate phylogrouping analysis section"""
        patterns = kwargs['patterns']
        phylogroup_dist = patterns.get('phylogroup_distribution', Counter())
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Phylogrouping Analysis</h3>
                <p><strong>{len(phylogroup_dist)} unique phylogroups</strong> identified using the Clermont scheme.</p>
            </div>
        </div>
        
        <h3><i class="fas fa-chart-bar"></i> Phylogroup Distribution</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Phylogroup</th>
                        <th>Count</th>
                        <th>Percentage</th>
                        <th>Common STs</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        total = sum(phylogroup_dist.values())
        for phylogroup, count in phylogroup_dist.most_common():
            if phylogroup == 'ND':
                continue
                
            percentage = (count / total) * 100
            
            # Find STs in this phylogroup
            sts_in_phylogroup = []
            for sample, data in kwargs['samples_data'].items():
                if data.get('phylogrouping', {}).get('Clermont_Type') == phylogroup:
                    st = data.get('mlst', {}).get('ST', 'ND')
                    if st != 'ND' and st not in sts_in_phylogroup:
                        sts_in_phylogroup.append(st)
            
            st_list = ', '.join([f"ST{s}" for s in sts_in_phylogroup]) if sts_in_phylogroup else 'Various'
            
            html += f"""
                    <tr>
                        <td><strong>{phylogroup}</strong></td>
                        <td>{count}</td>
                        <td>{percentage:.1f}%</td>
                        <td>{st_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_amr_section(self, kwargs: Dict) -> str:
        """Generate AMR genes section with gene-centric approach"""
        gene_centric = kwargs['gene_centric']
        amr_databases = gene_centric.get('amr_databases', {})
        
        html = """
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>AMR Gene Analysis</h3>
                <p>Each AMR gene is shown with ALL genomes that contain it. Search for specific genes or browse by database.</p>
            </div>
        </div>
        
        <input type="text" class="search-box" id="search-amr" 
               onkeyup="searchTable('amr-table', 'search-amr')" 
               placeholder="🔍 Search AMR genes by name or database...">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('amr-table', 'amr_genes.csv')">
                <i class="fas fa-download"></i> Export All AMR Genes
            </button>
        </div>
        
        <h3><i class="fas fa-shield-virus"></i> All AMR Genes Across Databases</h3>
        <div class="scrollable-table">
            <table id="amr-table" class="data-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Database</th>
                        <th>Frequency</th>
                        <th>Count</th>
                        <th>Genomes</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Combine all AMR genes from all databases
        all_amr_genes = []
        for db_name, genes in amr_databases.items():
            for gene_data in genes:
                all_amr_genes.append(gene_data)
        
        # Sort by count
        all_amr_genes.sort(key=lambda x: x['count'], reverse=True)
        
        for gene_data in all_amr_genes:
            gene = gene_data['gene']
            database = gene_data['database']
            frequency = gene_data['frequency']
            count = gene_data['count']
            genomes = gene_data.get('genomes', [])
            
            # Create genome tags - NO TRUNCATION
            genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in genomes])
            
            html += f"""
                    <tr>
                        <td><strong>{gene}</strong></td>
                        <td>{database}</td>
                        <td>{frequency}</td>
                        <td>{count}</td>
                        <td><div class="genome-list">{genome_tags}</div></td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-database"></i> AMR Databases Summary</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">
        """
        
        for db_name, genes in amr_databases.items():
            db_display = db_name.upper() if db_name != 'amrfinder' else 'AMRfinder'
            html += f"""
            <div class="database-section">
                <h4>{db_display}</h4>
                <p><strong>{len(genes)} unique genes</strong></p>
                <p>Top genes: {', '.join([g['gene'] for g in genes[:3]])}...</p>
            </div>
            """
        
        html += """
        </div>
        """
        
        return html
    
    def _generate_virulence_section(self, kwargs: Dict) -> str:
        """Generate virulence genes section with gene-centric approach"""
        gene_centric = kwargs['gene_centric']
        virulence_databases = gene_centric.get('virulence_databases', {})
        
        html = """
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Virulence Gene Analysis</h3>
                <p>Each virulence gene is shown with ALL genomes that contain it. Critical virulence genes are highlighted.</p>
            </div>
        </div>
        
        <input type="text" class="search-box" id="search-virulence" 
               onkeyup="searchTable('virulence-table', 'search-virulence')" 
               placeholder="🔍 Search virulence genes by name or database...">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('virulence-table', 'virulence_genes.csv')">
                <i class="fas fa-download"></i> Export All Virulence Genes
            </button>
        </div>
        
        <h3><i class="fas fa-virus"></i> All Virulence Genes Across Databases</h3>
        <div class="scrollable-table">
            <table id="virulence-table" class="data-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Database</th>
                        <th>Frequency</th>
                        <th>Count</th>
                        <th>Genomes</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Combine all virulence genes from all databases
        all_virulence_genes = []
        for db_name, genes in virulence_databases.items():
            for gene_data in genes:
                all_virulence_genes.append(gene_data)
        
        # Sort by count
        all_virulence_genes.sort(key=lambda x: x['count'], reverse=True)
        
        for gene_data in all_virulence_genes:
            gene = gene_data['gene']
            database = gene_data['database']
            frequency = gene_data['frequency']
            count = gene_data['count']
            genomes = gene_data.get('genomes', [])
            
            # Check if critical virulence gene
            is_critical = any(crit_gene in gene.lower() for crit_gene in self.data_analyzer.critical_virulence_genes)
            gene_display = f"<strong>{gene}</strong>" + (" ⚠️" if is_critical else "")
            
            # Create genome tags - NO TRUNCATION
            genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in genomes])
            
            html += f"""
                    <tr>
                        <td>{gene_display}</td>
                        <td>{database}</td>
                        <td>{frequency}</td>
                        <td>{count}</td>
                        <td><div class="genome-list">{genome_tags}</div></td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-database"></i> Virulence Databases Summary</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">
        """
        
        for db_name, genes in virulence_databases.items():
            db_display = db_name.upper()
            html += f"""
            <div class="database-section">
                <h4>{db_display}</h4>
                <p><strong>{len(genes)} unique virulence genes</strong></p>
                <p>Top genes: {', '.join([g['gene'] for g in genes[:3]])}...</p>
            </div>
            """
        
        html += """
        </div>
        """
        
        return html
    
    def _generate_pattern_discovery_section(self, kwargs: Dict) -> str:
        """Generate pattern discovery section"""
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        
        # Get critical genes from gene-centric data
        critical_amr_genes = []
        critical_virulence_genes = []
        
        for db_type in ['amr_databases', 'virulence_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    gene = gene_data['gene'].lower()
                    if any(crit in gene for crit in self.data_analyzer.critical_amr_genes):
                        critical_amr_genes.append(gene_data)
                    if any(crit in gene for crit in self.data_analyzer.critical_virulence_genes):
                        critical_virulence_genes.append(gene_data)
        
        html = """
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Cross-Genome Pattern Discovery</h3>
                <p>Discover associations between genes and identify high-risk combinations across all samples.</p>
            </div>
        </div>
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('high-risk-table', 'high_risk_combinations.csv')">
                <i class="fas fa-download"></i> Export High-Risk Combos
            </button>
        </div>
        """
        
        # High-risk combinations
        high_risk_combinations = patterns.get('high_risk_combinations', [])
        if high_risk_combinations:
            html += f"""
            <h3><i class="fas fa-exclamation-triangle"></i> High-Risk Combinations ({len(high_risk_combinations)})</h3>
            <div class="alert-box alert-danger">
                <i class="fas fa-radiation fa-2x"></i>
                <div>
                    <h3>⚠️ Critical Alert</h3>
                    <p><strong>{len(high_risk_combinations)} samples</strong> contain dangerous combinations of critical AMR and virulence genes.</p>
                </div>
            </div>
            
            <div class="scrollable-table">
                <table id="high-risk-table" class="data-table">
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>ST</th>
                            <th>Serotype</th>
                            <th>Critical AMR Genes</th>
                            <th>Critical Virulence Genes</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for combo in high_risk_combinations:
                # FIXED: Use comma separation for critical AMR genes
                amr_genes = ', '.join(combo['critical_amr_genes'])
                vir_genes = ', '.join(combo['critical_virulence_genes'])
                
                html += f"""
                        <tr>
                            <td><strong>{combo['sample']}</strong></td>
                            <td>{combo['st']}</td>
                            <td>{combo['serotype']}</td>
                            <td>{amr_genes}</td>
                            <td>{vir_genes}</td>
                        </tr>
                """
            
            html += """
                    </tbody>
                </table>
            </div>
            """
        
        # Critical genes summary - SHOW ALL, NO LIMIT
        html += f"""
        <h3 style="margin-top: 30px;"><i class="fas fa-skull-crossbones"></i> Critical Genes Summary</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; margin: 20px 0;">
            <div class="database-section">
                <h4>Critical AMR Genes</h4>
                <p><strong>{len(critical_amr_genes)} critical AMR genes</strong> identified</p>
                <div style="margin-top: 10px;">
        """
        
        # Show ALL critical AMR genes, no limit - COMMA SEPARATED
        if critical_amr_genes:
            # Create comma-separated list of all critical AMR genes
            amr_gene_names = [gene_data['gene'] for gene_data in critical_amr_genes]
            amr_gene_list = ', '.join(sorted(set(amr_gene_names)))
            html += f"""
                    <div style="padding: 10px; background: #fff3cd; border-radius: 5px;">
                        <strong>All Critical AMR Genes:</strong><br>
                        {amr_gene_list}
                    </div>
            """
        else:
            html += """
                    <div style="padding: 10px; background: #d4edda; border-radius: 5px;">
                        No critical AMR genes detected.
                    </div>
            """
        
        html += """
                </div>
            </div>
            
            <div class="database-section">
                <h4>Critical Virulence Genes</h4>
                <p><strong>critical virulence genes</strong> identified</p>
                <div style="margin-top: 10px;">
        """
        
        # Show ALL critical virulence genes, no limit - COMMA SEPARATED
        if critical_virulence_genes:
            # Create comma-separated list of all critical virulence genes
            vir_gene_names = [gene_data['gene'] for gene_data in critical_virulence_genes]
            vir_gene_list = ', '.join(sorted(set(vir_gene_names)))
            html += f"""
                    <div style="padding: 10px; background: #fff3cd; border-radius: 5px;">
                        <strong>All Critical Virulence Genes:</strong><br>
                        {vir_gene_list}
                    </div>
            """
        else:
            html += """
                    <div style="padding: 10px; background: #d4edda; border-radius: 5px;">
                        No critical virulence genes detected.
                    </div>
            """
        
        html += """
                </div>
            </div>
        </div>
        """
        
        return html   
    
    def _generate_export_section(self, kwargs: Dict) -> str:
        """Generate export section"""
        return """
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Export Data and Reports</h3>
                <p>Download comprehensive data in various formats for further analysis and reporting.</p>
            </div>
        </div>
        
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 30px 0;">
            <div class="dashboard-card card-export" onclick="exportTableToCSV('samples-table', 'sample_overview.csv')">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-file-csv"></i></div>
                <div class="card-label">Sample Overview CSV</div>
                <p style="font-size: 0.9em; margin-top: 10px;">All samples with ST, serotype, phylogroup, CH type</p>
            </div>
            
            <div class="dashboard-card card-export" onclick="exportTableToCSV('amr-table', 'amr_genes.csv')">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-biohazard"></i></div>
                <div class="card-label">AMR Genes CSV</div>
                <p style="font-size: 0.9em; margin-top: 10px;">All AMR genes with genomes and frequencies</p>
            </div>
            
            <div class="dashboard-card card-export" onclick="exportTableToCSV('virulence-table', 'virulence_genes.csv')">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-virus"></i></div>
                <div class="card-label">Virulence Genes CSV</div>
                <p style="font-size: 0.9em; margin-top: 10px;">All virulence genes with genomes and frequencies</p>
            </div>
            
            <div class="dashboard-card card-export" onclick="location.href='genius_ultimate_report.json'">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-file-code"></i></div>
                <div class="card-label">Complete JSON Data</div>
                <p style="font-size: 0.9em; margin-top: 10px;">All analysis data in structured JSON format</p>
            </div>
        </div>
        
        <h3><i class="fas fa-download"></i> Available Export Files</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>File</th>
                        <th>Description</th>
                        <th>Format</th>
                        <th>Contents</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td><strong>genius_ultimate_report.html</strong></td>
                        <td>This interactive HTML report</td>
                        <td>HTML</td>
                        <td>Complete analysis with all sections</td>
                    </tr>
                    <tr>
                        <td><strong>genius_ultimate_report.json</strong></td>
                        <td>Complete structured data</td>
                        <td>JSON</td>
                        <td>All analysis data for programmatic use</td>
                    </tr>
                    <tr>
                        <td><strong>sample_overview.csv</strong></td>
                        <td>Sample overview data</td>
                        <td>CSV</td>
                        <td>All samples with typing results</td>
                    </tr>
                    <tr>
                        <td><strong>amr_genes.csv</strong></td>
                        <td>AMR gene analysis</td>
                        <td>CSV</td>
                        <td>All AMR genes with genomes and frequencies</td>
                    </tr>
                    <tr>
                        <td><strong>virulence_genes.csv</strong></td>
                        <td>Virulence gene analysis</td>
                        <td>CSV</td>
                        <td>All virulence genes with genomes and frequencies</td>
                    </tr>
                    <tr>
                        <td><strong>pattern_discovery.csv</strong></td>
                        <td>Pattern discovery results</td>
                        <td>CSV</td>
                        <td>Cross-genome patterns and associations</td>
                    </tr>
                </tbody>
            </table>
        </div>
        """


class GeniusUltimateReporter:
    """MASTER CLASS: Generates ultimate gene-centric reports"""
    
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "GENIUS_ULTIMATE_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize components
        self.parser = UltimateHTMLParser()
        self.analyzer = UltimateDataAnalyzer()
        self.html_generator = UltimateHTMLGenerator(self.analyzer)
        
        # Metadata
        self.metadata = {
            "tool_name": "GENIUS E.coli Ultimate Reporter",
            "version": "1.0.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir)
        }
    
    def find_html_files(self) -> Dict[str, List[Path]]:
        """Find all HTML report files"""
        print("🔍 Searching for HTML reports...")
        
        html_files = {
            'mlst': [],
            'serotype': [],
            'chtyper': [],
            'phylogrouping': [],
            'amrfinder': [],
            'abricate': []
        }
        
        for html_file in self.input_dir.glob("**/*.html"):
            filename = html_file.name.lower()
            
            if 'mlst' in filename:
                html_files['mlst'].append(html_file)
            elif 'serotype' in filename:
                html_files['serotype'].append(html_file)
            elif 'chtyper' in filename or 'fumc' in filename or 'fimh' in filename:
                html_files['chtyper'].append(html_file)
            elif 'phylogroup' in filename or 'clermont' in filename:
                html_files['phylogrouping'].append(html_file)
            elif 'amr' in filename or 'resistance' in filename:
                html_files['amrfinder'].append(html_file)
            elif 'abricate' in filename or any(db in filename for db in self.parser.abricate_databases):
                html_files['abricate'].append(html_file)
        
        # Print findings
        for file_type, files in html_files.items():
            if files:
                print(f"  📁 {file_type.upper()}: {len(files)} files found")
        
        return html_files
    
    def integrate_all_data(self, html_files: Dict[str, List[Path]]) -> Dict[str, Any]:
        """Integrate data from all HTML reports"""
        print("\n🔗 Integrating data from all reports...")
        
        integrated_data = {
            'metadata': self.metadata,
            'samples': {},
            'patterns': {},
            'gene_centric': {}
        }
        
        # Parse all reports
        mlst_data = {}
        if html_files['mlst']:
            mlst_data = self.parser.parse_mlst_report(html_files['mlst'][0])
        
        serotype_data = {}
        if html_files['serotype']:
            serotype_data = self.parser.parse_serotype_report(html_files['serotype'][0])
        
        chtyper_data = {}
        if html_files['chtyper']:
            chtyper_data = self.parser.parse_chtyper_report(html_files['chtyper'][0])
        
        phylogrouping_data = {}
        if html_files['phylogrouping']:
            phylogrouping_data = self.parser.parse_phylogrouping_report(html_files['phylogrouping'][0])
        
        amr_by_sample, amr_gene_freq = {}, {}
        if html_files['amrfinder']:
            amr_by_sample, amr_gene_freq = self.parser.parse_amrfinder_report(html_files['amrfinder'][0])
        
        # Parse ABRicate databases
        abricate_data = {}
        abricate_gene_freq = {}
        
        for abricate_file in html_files['abricate']:
            genes_by_sample, gene_freq = self.parser.parse_abricate_database_report(abricate_file)
            
            db_name = 'unknown'
            for db in self.parser.abricate_databases:
                if db in abricate_file.name.lower():
                    db_name = db
                    break
            
            abricate_data[db_name] = genes_by_sample
            abricate_gene_freq[db_name] = gene_freq
        
        # Combine all samples
        all_samples = set()
        all_samples.update(mlst_data.keys())
        all_samples.update(serotype_data.keys())
        all_samples.update(chtyper_data.keys())
        all_samples.update(phylogrouping_data.keys())
        all_samples.update(amr_by_sample.keys())
        
        for db_data in abricate_data.values():
            all_samples.update(db_data.keys())
        
        all_samples = sorted(list(all_samples))
        
        if not all_samples:
            print("❌ No samples found in any report!")
            return {}
        
        print(f"\n📊 Found {len(all_samples)} unique samples")
        
        # Integrate data for each sample
        for sample in all_samples:
            # Get virulence genes from virulence databases
            virulence_genes = []
            for db_name in ['vfdb', 'ecoli_vf']:
                if db_name in abricate_data:
                    virulence_genes.extend(abricate_data[db_name].get(sample, []))
            
            sample_data = {
                'mlst': mlst_data.get(sample, {'ST': 'ND', 'Allele_Profile': 'ND'}),
                'serotype': serotype_data.get(sample, {'Serotype': 'ND', 'O_Type': 'ND', 'H_Type': 'ND'}),
                'chtyper': chtyper_data.get(sample, {'FumC_Type': 'ND', 'FimH_Type': 'ND', 'CH_Type': 'ND'}),
                'phylogrouping': phylogrouping_data.get(sample, {'Clermont_Type': 'ND'}),
                'amr_genes': amr_by_sample.get(sample, []),
                'virulence_genes': list(set(virulence_genes))  # Remove duplicates
            }
            
            integrated_data['samples'][sample] = sample_data
        
        # Store gene frequencies
        integrated_data['gene_frequencies'] = {
            'amrfinder': amr_gene_freq,
            'abricate': abricate_gene_freq
        }
        
        # Process gene-centric data and patterns
        print("\n🧠 Processing gene-centric analysis...")
        integrated_data['gene_centric'] = self.analyzer.create_gene_centric_tables(integrated_data)
        integrated_data['patterns'] = self.analyzer.create_cross_genome_patterns(integrated_data)
        
        return integrated_data
    
    def generate_json_report(self, integrated_data: Dict[str, Any]) -> Path:
        """Generate comprehensive JSON report"""
        print("\n📝 Generating JSON report...")
        
        output_file = self.output_dir / "genius_ultimate_report.json"
        
        # Create serializable copy
        def make_serializable(obj):
            if obj is None:
                return None
            elif isinstance(obj, (str, int, float, bool)):
                return obj
            elif isinstance(obj, (list, tuple)):
                return [make_serializable(item) for item in obj]
            elif isinstance(obj, dict):
                return {str(k): make_serializable(v) for k, v in obj.items()}
            elif isinstance(obj, set):
                return [make_serializable(item) for item in obj]
            elif isinstance(obj, (Counter, defaultdict)):
                return {str(k): make_serializable(v) for k, v in obj.items()}
            elif isinstance(obj, Path):
                return str(obj)
            elif hasattr(obj, 'isoformat'):
                return obj.isoformat()
            else:
                try:
                    return str(obj)
                except:
                    return None
        
        # Create serializable data
        serializable_data = make_serializable(integrated_data)
        
        # Write to file
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(serializable_data, f, indent=2, ensure_ascii=False)
        
        print(f"    ✅ JSON report saved: {output_file}")
        return output_file
    
    def generate_csv_reports(self, integrated_data: Dict[str, Any]):
        """Generate multiple CSV reports"""
        print("\n📊 Generating CSV reports...")
        
        # 1. Sample summary
        samples_data = []
        for sample, data in integrated_data['samples'].items():
            row = {
                'Sample': sample,
                'ST': data['mlst']['ST'],
                'Allele_Profile': data['mlst']['Allele_Profile'],
                'Serotype': data['serotype']['Serotype'],
                'O_Type': data['serotype']['O_Type'],
                'H_Type': data['serotype']['H_Type'],
                'FumC_Type': data['chtyper']['FumC_Type'],
                'FimH_Type': data['chtyper']['FimH_Type'],
                'CH_Type': data['chtyper']['CH_Type'],
                'Clermont_Type': data['phylogrouping']['Clermont_Type'],
                'AMR_Gene_Count': len(data['amr_genes']),
                'Virulence_Gene_Count': len(data['virulence_genes'])
            }
            samples_data.append(row)
        
        df_samples = pd.DataFrame(samples_data)
        samples_file = self.output_dir / "sample_overview.csv"
        df_samples.to_csv(samples_file, index=False)
        
        # 2. AMR genes (gene-centric)
        amr_data = []
        gene_centric = integrated_data.get('gene_centric', {})
        
        for db_name, genes in gene_centric.get('amr_databases', {}).items():
            for gene_info in genes:
                amr_data.append({
                    'Gene': gene_info['gene'],
                    'Database': gene_info['database'],
                    'Frequency': gene_info['frequency'],
                    'Count': gene_info['count'],
                    'Genomes': ';'.join(gene_info.get('genomes', []))
                })
        
        if amr_data:
            df_amr = pd.DataFrame(amr_data)
            amr_file = self.output_dir / "amr_genes.csv"
            df_amr.to_csv(amr_file, index=False)
        
        # 3. Virulence genes (gene-centric)
        virulence_data = []
        for db_name, genes in gene_centric.get('virulence_databases', {}).items():
            for gene_info in genes:
                virulence_data.append({
                    'Gene': gene_info['gene'],
                    'Database': gene_info['database'],
                    'Frequency': gene_info['frequency'],
                    'Count': gene_info['count'],
                    'Genomes': ';'.join(gene_info.get('genomes', []))
                })
        
        if virulence_data:
            df_virulence = pd.DataFrame(virulence_data)
            virulence_file = self.output_dir / "virulence_genes.csv"
            df_virulence.to_csv(virulence_file, index=False)
        
        # 4. Pattern discovery
        pattern_data = []
        patterns = integrated_data['patterns']
        
        # ST distribution
        for st, count in patterns.get('st_distribution', Counter()).items():
            pattern_data.append({
                'Pattern_Type': 'ST_Distribution',
                'ST': st,
                'Count': count
            })
        
        # ST-Serotype combinations
        for combo, samples in patterns.get('st_serotype_combinations', {}).items():
            pattern_data.append({
                'Pattern_Type': 'ST_Serotype_Combination',
                'Combination': combo,
                'Samples': ';'.join(samples),
                'Count': len(samples)
            })
        
        # High-risk combinations
        for combo in patterns.get('high_risk_combinations', []):
            pattern_data.append({
                'Pattern_Type': 'High_Risk_Combination',
                'Sample': combo['sample'],
                'ST': combo['st'],
                'Serotype': combo['serotype'],
                'Critical_AMR_Genes': ';'.join(combo['critical_amr_genes']),
                'Critical_Virulence_Genes': ';'.join(combo['critical_virulence_genes'])
            })
        
        if pattern_data:
            df_patterns = pd.DataFrame(pattern_data)
            patterns_file = self.output_dir / "pattern_discovery.csv"
            df_patterns.to_csv(patterns_file, index=False)
        
        print(f"    ✅ CSV reports generated: sample_overview.csv, amr_genes.csv, virulence_genes.csv, pattern_discovery.csv")
    
    def run(self):
        """Run the complete analysis"""
        print("=" * 80)
        print("🧠 GENIUS E.COLI ULTIMATE REPORTER v1.0.0")
        print("=" * 80)
        print(f"📁 Input directory: {self.input_dir}")
        
        # Find HTML files
        html_files = self.find_html_files()
        
        if not any(html_files.values()):
            print("❌ No HTML report files found!")
            return False
        
        # Integrate data
        integrated_data = self.integrate_all_data(html_files)
        if not integrated_data:
            return False
        
        # Generate reports
        print("\n" + "=" * 80)
        print("📊 GENERATING ULTIMATE REPORTS")
        print("=" * 80)
        
        # Generate JSON
        json_file = self.generate_json_report(integrated_data)
        
        # Generate CSV
        self.generate_csv_reports(integrated_data)
        
        # Generate HTML
        html_file = self.html_generator.generate_main_report(integrated_data, self.output_dir)
        
        # Print summary
        total_samples = len(integrated_data['samples'])
        patterns = integrated_data['patterns']
        high_risk = len(patterns.get('high_risk_combinations', []))
        gene_centric = integrated_data['gene_centric']
        
        total_amr_genes = sum(len(genes) for genes in gene_centric.get('amr_databases', {}).values())
        total_virulence_genes = sum(len(genes) for genes in gene_centric.get('virulence_databases', {}).values())
        
        print("\n" + "=" * 80)
        print("✅ ULTIMATE ANALYSIS COMPLETE!")
        print("=" * 80)
        print(f"📁 Output directory: {self.output_dir}")
        print(f"📄 Files generated:")
        print(f"   • genius_ultimate_report.html (Interactive report)")
        print(f"   • genius_ultimate_report.json (Complete data)")
        print(f"   • sample_overview.csv (Sample data)")
        print(f"   • amr_genes.csv (Gene-centric AMR data)")
        print(f"   • virulence_genes.csv (Gene-centric virulence data)")
        print(f"   • pattern_discovery.csv (Pattern analysis)")
        
        print(f"\n🔬 KEY FEATURES:")
        print(f"   • Gene-centric approach: Genes shown with all genomes")
        print(f"   • NO TRUNCATION: Complete genome lists for each gene")
        print(f"   • Section-specific printing: Print button on each section")
        print(f"   • ST-Serotype combinations: Clear association display")
        print(f"   • Removed all fixed-width constraints: Tables expand as needed")
        
        print(f"\n📈 ANALYSIS SUMMARY:")
        print(f"   • {total_samples} total samples analyzed")
        print(f"   • {total_amr_genes} AMR genes across all databases")
        print(f"   • {total_virulence_genes} virulence genes")
        print(f"   • {high_risk} high-risk AMR+virulence combinations")
        
        print("\n🎯 Next steps:")
        print("   1. Open genius_ultimate_report.html in your browser")
        print("   2. Use AMR and Virulence tabs to see genes with ALL their genomes")
        print("   3. Use print buttons on each section header to print specific sections")
        print("   4. Export data using the Export tab or individual CSV buttons")
        
        print("\n" + "=" * 80)
        return True


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='GENIUS E.coli Ultimate Reporter - Gene-Centric Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python genius_reporter.py -i /path/to/html/reports
  
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
        """
    )
    
    parser.add_argument('-i', '--input-dir', required=True,
                       help='Directory containing HTML report files')
    parser.add_argument('-o', '--output-dir',
                       help='Custom output directory')
    
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    
    if not input_dir.exists():
        print(f"❌ Input directory not found: {input_dir}")
        sys.exit(1)
    
    # Create and run reporter
    reporter = GeniusUltimateReporter(input_dir)
    
    if args.output_dir:
        reporter.output_dir = Path(args.output_dir)
        reporter.output_dir.mkdir(parents=True, exist_ok=True)
    
    success = reporter.run()
    
    if not success:
        print("❌ Report generation failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()