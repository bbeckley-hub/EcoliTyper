#!/usr/bin/env python3
"""
GENIUS E.COLI VISUALIZATION MODULE - ULTIMATE HTML EDITION
Parses EcoliTyper HTML reports and generates publication-ready visualizations
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 1.3.0 
Date: 2026-07-22
"""

import os
import sys
import re
import glob
import json
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any, Optional, Union
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

# HTML parsing
from bs4 import BeautifulSoup

# Visualization
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
from scipy import stats

# Set publication-quality style
plt.style.use('seaborn-v0_8-whitegrid')
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['savefig.pad_inches'] = 0.1
mpl.rcParams['figure.max_open_warning'] = 50

# Color palettes
COLOR_PALETTES = {
    'mlst': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
             '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
             '#9edae5', '#dbdb8d', '#c7c7c7', '#f7b6d2', '#c49c94'],
    'serotype': ['#4c78a8', '#f58518', '#e45756', '#72b7b2', '#54a24b',
                 '#eeca3b', '#b279a2', '#ff9da6', '#9d755d', '#bab0ac',
                 '#6a9fb5', '#d4b9da', '#00cc96', '#ab63fa', '#ffa15a'],
    'phylogroup': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                   '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
    'chtype': ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3',
               '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd',
               '#ccebc5', '#ffed6f', '#e78ac3', '#a6d854', '#ffd92f'],
    'database': ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00',
                 '#ffff33', '#a65628', '#f781bf', '#999999', '#66c2a5']
}


class UltimateHTMLParser:
    """Parses ALL EcoliTyper HTML reports with BeautifulSoup"""
    
    def __init__(self):
        self.normalized_samples = {}
        
    def normalize_sample_name(self, sample_name: str) -> str:
        """Normalize sample name by removing extensions and paths"""
        sample = str(sample_name)
        
        # Remove common extensions
        extensions = ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', 
                     '.fna.fna', '.fasta.fasta', '.fa.fa']
        for ext in extensions:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        
        # Remove path if present
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        
        # Remove any remaining extensions
        for ext in ['.fna', '.fasta', '.fa']:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        
        return sample.strip()
    
    def parse_mlst_html(self, file_path: Path) -> pd.DataFrame:
        """Parse MLST HTML report - gets ALL STs"""
        print(f"  📄 Parsing MLST HTML: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            
            # Find the main results table
            tables = soup.find_all('table')
            mlst_table = None
            
            for table in tables:
                if table.find('th') and any('ST' in th.get_text() for th in table.find_all('th')):
                    mlst_table = table
                    break
            
            if not mlst_table:
                # Try to find by text content
                for table in tables:
                    if table.get_text().strip() and 'ST' in table.get_text():
                        mlst_table = table
                        break
            
            if not mlst_table:
                raise ValueError("Could not find MLST table in HTML")
            
            # Parse table headers
            headers = []
            header_row = mlst_table.find('tr')
            if header_row:
                for th in header_row.find_all(['th', 'td']):
                    headers.append(th.get_text().strip())
            
            # Parse table rows
            rows = []
            for tr in mlst_table.find_all('tr')[1:]:  # Skip header
                cols = tr.find_all(['td', 'th'])
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    if len(row_data) >= 2:  # At least sample and ST
                        rows.append(row_data)
            
            # Create DataFrame
            if len(headers) > 0 and len(rows) > 0:
                # Ensure headers match columns
                if len(headers) > len(rows[0]):
                    headers = headers[:len(rows[0])]
                elif len(headers) < len(rows[0]):
                    headers = headers + [f'Col{i}' for i in range(len(headers), len(rows[0]))]
                
                df = pd.DataFrame(rows, columns=headers)
                
                # Identify sample column
                sample_col = None
                for col in df.columns:
                    if 'sample' in col.lower() or 'genome' in col.lower() or 'id' in col.lower():
                        sample_col = col
                        break
                
                if not sample_col and len(df.columns) > 0:
                    sample_col = df.columns[0]
                
                # Identify ST column
                st_col = None
                for col in df.columns:
                    if col.lower() == 'st' or 'sequence type' in col.lower():
                        st_col = col
                        break
                
                if not st_col:
                    # Look for column containing 'ST'
                    for col in df.columns:
                        if any(row and 'ST' in str(row) for row in df[col]):
                            st_col = col
                            break
                
                if not st_col and len(df.columns) > 1:
                    st_col = df.columns[1]
                
                # Normalize sample names
                if sample_col:
                    df['Sample'] = df[sample_col].apply(self.normalize_sample_name)
                
                # Extract ST values
                if st_col:
                    df['ST'] = df[st_col].astype(str).str.replace('ST', '').str.strip()
                
                # Keep only essential columns
                keep_cols = ['Sample', 'ST']
                for col in ['Allele Profile', 'adk', 'fumC', 'gyrB', 'icd', 'mdh', 'purA', 'recA']:
                    if col in df.columns:
                        keep_cols.append(col)
                
                df = df[[col for col in keep_cols if col in df.columns]]
                
                print(f"    ✓ Found {len(df)} samples, {df['ST'].nunique()} unique STs")
                return df
            
            return pd.DataFrame()
            
        except Exception as e:
            print(f"    ✗ Error parsing MLST HTML: {e}")
            import traceback
            traceback.print_exc()
            return pd.DataFrame()
    
    def parse_serotype_html(self, file_path: Path) -> pd.DataFrame:
        """Parse Serotype HTML report - gets ALL serotypes"""
        print(f"  📄 Parsing Serotype HTML: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            
            # Find the main results table
            tables = soup.find_all('table')
            serotype_table = None
            
            for table in tables:
                if table.find('th') and any('serotype' in th.get_text().lower() for th in table.find_all('th')):
                    serotype_table = table
                    break
            
            if not serotype_table:
                # Try to find by text content
                for table in tables:
                    if table.get_text().strip() and ('O-type' in table.get_text() or 'H-type' in table.get_text()):
                        serotype_table = table
                        break
            
            if not serotype_table:
                raise ValueError("Could not find serotype table in HTML")
            
            # Parse table headers
            headers = []
            header_row = serotype_table.find('tr')
            if header_row:
                for th in header_row.find_all(['th', 'td']):
                    headers.append(th.get_text().strip())
            
            # Parse table rows
            rows = []
            for tr in serotype_table.find_all('tr')[1:]:  # Skip header
                cols = tr.find_all(['td', 'th'])
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    if len(row_data) >= 2:
                        rows.append(row_data)
            
            # Create DataFrame
            if len(headers) > 0 and len(rows) > 0:
                # Ensure headers match columns
                if len(headers) > len(rows[0]):
                    headers = headers[:len(rows[0])]
                elif len(headers) < len(rows[0]):
                    headers = headers + [f'Col{i}' for i in range(len(headers), len(rows[0]))]
                
                df = pd.DataFrame(rows, columns=headers)
                
                # Identify sample column
                sample_col = None
                for col in df.columns:
                    if 'sample' in col.lower() or 'id' in col.lower():
                        sample_col = col
                        break
                
                if not sample_col and len(df.columns) > 0:
                    sample_col = df.columns[0]
                
                # Identify serotype columns
                serotype_col = None
                o_type_col = None
                h_type_col = None
                
                for col in df.columns:
                    col_lower = col.lower()
                    if 'serotype' in col_lower:
                        serotype_col = col
                    elif 'o-type' in col_lower or 'o_type' in col_lower:
                        o_type_col = col
                    elif 'h-type' in col_lower or 'h_type' in col_lower:
                        h_type_col = col
                
                # Normalize sample names
                if sample_col:
                    df['Sample'] = df[sample_col].apply(self.normalize_sample_name)
                
                # Extract serotype values
                if serotype_col:
                    df['Serotype'] = df[serotype_col].astype(str).str.strip()
                else:
                    # Create serotype from O and H types if available
                    if o_type_col and h_type_col:
                        df['Serotype'] = df[o_type_col].astype(str) + ':' + df[h_type_col].astype(str)
                    elif o_type_col:
                        df['Serotype'] = df[o_type_col].astype(str) + ':H-ND'
                    elif h_type_col:
                        df['Serotype'] = 'O-ND:' + df[h_type_col].astype(str)
                    else:
                        df['Serotype'] = 'O-ND:H-ND'
                
                # Extract O and H types separately
                if o_type_col:
                    df['O_Type'] = df[o_type_col].astype(str).str.strip()
                else:
                    # Extract from serotype
                    df['O_Type'] = df['Serotype'].apply(lambda x: x.split(':')[0] if ':' in str(x) else 'ND')
                
                if h_type_col:
                    df['H_Type'] = df[h_type_col].astype(str).str.strip()
                else:
                    # Extract from serotype
                    df['H_Type'] = df['Serotype'].apply(lambda x: x.split(':')[1] if ':' in str(x) else 'ND')
                
                # Keep only essential columns
                keep_cols = ['Sample', 'Serotype', 'O_Type', 'H_Type']
                df = df[[col for col in keep_cols if col in df.columns]]
                
                print(f"    ✓ Found {len(df)} samples, {df['Serotype'].nunique()} unique serotypes")
                return df
            
            return pd.DataFrame()
            
        except Exception as e:
            print(f"    ✗ Error parsing Serotype HTML: {e}")
            import traceback
            traceback.print_exc()
            return pd.DataFrame()
    
    def parse_phylogroup_html(self, file_path: Path) -> pd.DataFrame:
        """Parse Phylogroup HTML report - gets ALL phylogroups"""
        print(f"  📄 Parsing Phylogroup HTML: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            
            # Try to find distribution data first (pie chart data)
            distribution_data = {}
            dist_sections = soup.find_all(string=re.compile(r'Clermont Type Distribution|Phylogroup Distribution', re.I))
            
            # Look for distribution in text or tables
            for text in soup.stripped_strings:
                if '(' in text and ')' in text:
                    # Pattern like "B1 (21)" or "A (3)"
                    match = re.search(r'([A-Za-z0-9]+)\s*\((\d+)\)', text)
                    if match:
                        phylogroup = match.group(1)
                        count = int(match.group(2))
                        distribution_data[phylogroup] = count
            
            # Find the main results table
            tables = soup.find_all('table')
            phylogroup_table = None
            
            for table in tables:
                if table.find('th') and any('clermont' in th.get_text().lower() for th in table.find_all('th')):
                    phylogroup_table = table
                    break
            
            if not phylogroup_table:
                # Try to find by text content
                for table in tables:
                    if table.get_text().strip() and any(pg in table.get_text() for pg in ['A', 'B1', 'B2', 'C', 'D', 'E', 'F', 'G']):
                        phylogroup_table = table
                        break
            
            # Parse from table if available
            if phylogroup_table:
                # Parse table headers
                headers = []
                header_row = phylogroup_table.find('tr')
                if header_row:
                    for th in header_row.find_all(['th', 'td']):
                        headers.append(th.get_text().strip())
                
                # Parse table rows
                rows = []
                for tr in phylogroup_table.find_all('tr')[1:]:  # Skip header
                    cols = tr.find_all(['td', 'th'])
                    if cols:
                        row_data = [col.get_text().strip() for col in cols]
                        if len(row_data) >= 2:
                            rows.append(row_data)
                
                # Create DataFrame
                if len(headers) > 0 and len(rows) > 0:
                    # Ensure headers match columns
                    if len(headers) > len(rows[0]):
                        headers = headers[:len(rows[0])]
                    elif len(headers) < len(rows[0]):
                        headers = headers + [f'Col{i}' for i in range(len(headers), len(rows[0]))]
                    
                    df = pd.DataFrame(rows, columns=headers)
                    
                    # Identify sample column
                    sample_col = None
                    for col in df.columns:
                        if 'sample' in col.lower() or 'id' in col.lower():
                            sample_col = col
                            break
                    
                    if not sample_col and len(df.columns) > 0:
                        sample_col = df.columns[0]
                    
                    # Identify phylogroup column
                    pg_col = None
                    for col in df.columns:
                        col_lower = col.lower()
                        if 'clermont' in col_lower or 'phylogroup' in col_lower or 'type' in col_lower:
                            pg_col = col
                            break
                    
                    if not pg_col and len(df.columns) > 1:
                        pg_col = df.columns[1]
                    
                    # Normalize sample names
                    if sample_col:
                        df['Sample'] = df[sample_col].apply(self.normalize_sample_name)
                    
                    # Extract phylogroup values
                    if pg_col:
                        df['Phylogroup'] = df[pg_col].astype(str).str.strip()
                    
                    # Keep only essential columns
                    df = df[['Sample', 'Phylogroup']] if 'Sample' in df.columns and 'Phylogroup' in df.columns else pd.DataFrame()
            
            else:
                # Create DataFrame from distribution data if no table found
                df = pd.DataFrame({
                    'Sample': [f'Sample_{i}' for i in range(sum(distribution_data.values()))],
                    'Phylogroup': []
                })
                
                for pg, count in distribution_data.items():
                    df = pd.concat([df, pd.DataFrame({
                        'Sample': [f'Sample_{len(df) + i}' for i in range(count)],
                        'Phylogroup': [pg] * count
                    })], ignore_index=True)
            
            if not df.empty:
                print(f"    ✓ Found {len(df)} samples, {df['Phylogroup'].nunique()} unique phylogroups")
                return df
            else:
                print("    ⚠️ No phylogroup data found in HTML")
                return pd.DataFrame()
            
        except Exception as e:
            print(f"    ✗ Error parsing Phylogroup HTML: {e}")
            import traceback
            traceback.print_exc()
            return pd.DataFrame()
    
    def parse_chtyper_html(self, file_path: Path) -> pd.DataFrame:
        """Parse CHTyper HTML report - gets ALL CH types"""
        print(f"  📄 Parsing CHTyper HTML: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            
            # Try to find distribution data first
            fumc_dist = {}
            fimh_dist = {}
            
            # Look for FumC and FimH distributions
            for text in soup.stripped_strings:
                if 'FumC Types:' in text or 'FimH Types:' in text:
                    lines = text.split('\n')
                    for line in lines:
                        if '(' in line and ')' in line:
                            match = re.search(r'([a-zA-Z0-9]+)\s*\((\d+)\)', line)
                            if match:
                                gene_type = match.group(1)
                                count = int(match.group(2))
                                if 'FumC' in text:
                                    fumc_dist[gene_type] = count
                                elif 'FimH' in text:
                                    fimh_dist[gene_type] = count
            
            # Find the main results table
            tables = soup.find_all('table')
            chtyper_table = None
            
            for table in tables:
                if table.find('th') and any('fumc' in th.get_text().lower() or 'fimh' in th.get_text().lower() for th in table.find_all('th')):
                    chtyper_table = table
                    break
            
            if not chtyper_table:
                # Try to find by text content
                for table in tables:
                    if table.get_text().strip() and ('FumC' in table.get_text() and 'FimH' in table.get_text()):
                        chtyper_table = table
                        break
            
            # Parse from table if available
            if chtyper_table:
                # Parse table headers
                headers = []
                header_row = chtyper_table.find('tr')
                if header_row:
                    for th in header_row.find_all(['th', 'td']):
                        headers.append(th.get_text().strip())
                
                # Parse table rows
                rows = []
                for tr in chtyper_table.find_all('tr')[1:]:  # Skip header
                    cols = tr.find_all(['td', 'th'])
                    if cols:
                        row_data = [col.get_text().strip() for col in cols]
                        if len(row_data) >= 2:
                            rows.append(row_data)
                
                # Create DataFrame
                if len(headers) > 0 and len(rows) > 0:
                    # Ensure headers match columns
                    if len(headers) > len(rows[0]):
                        headers = headers[:len(rows[0])]
                    elif len(headers) < len(rows[0]):
                        headers = headers + [f'Col{i}' for i in range(len(headers), len(rows[0]))]
                    
                    df = pd.DataFrame(rows, columns=headers)
                    
                    # Identify sample column
                    sample_col = None
                    for col in df.columns:
                        if 'sample' in col.lower() or 'id' in col.lower():
                            sample_col = col
                            break
                    
                    if not sample_col and len(df.columns) > 0:
                        sample_col = df.columns[0]
                    
                    # Identify FumC and FimH columns
                    fumc_col = None
                    fimh_col = None
                    
                    for col in df.columns:
                        col_lower = col.lower()
                        if 'fumc' in col_lower:
                            fumc_col = col
                        elif 'fimh' in col_lower:
                            fimh_col = col
                    
                    # Normalize sample names
                    if sample_col:
                        df['Sample'] = df[sample_col].apply(self.normalize_sample_name)
                    
                    # Extract FumC and FimH values
                    if fumc_col:
                        df['FumC_Type'] = df[fumc_col].astype(str).str.strip()
                    
                    if fimh_col:
                        df['FimH_Type'] = df[fimh_col].astype(str).str.strip()
                    
                    # Create CH Type (FumC:FimH)
                    if 'FumC_Type' in df.columns and 'FimH_Type' in df.columns:
                        df['CH_Type'] = df['FumC_Type'] + ':' + df['FimH_Type']
                    elif 'FumC_Type' in df.columns:
                        df['CH_Type'] = df['FumC_Type'] + ':ND'
                    elif 'FimH_Type' in df.columns:
                        df['CH_Type'] = 'ND:' + df['FimH_Type']
                    else:
                        df['CH_Type'] = 'ND:ND'
                    
                    # Keep only essential columns
                    keep_cols = ['Sample', 'CH_Type', 'FumC_Type', 'FimH_Type']
                    df = df[[col for col in keep_cols if col in df.columns]]
            
            else:
                # Create DataFrame from distribution data if no table found
                df = pd.DataFrame({
                    'Sample': [f'Sample_{i}' for i in range(sum(fumc_dist.values()))],
                    'FumC_Type': [],
                    'FimH_Type': [],
                    'CH_Type': []
                })
                
                # This is simplified - in reality we'd need sample-level data
                print("    ⚠️ Only distribution data found, not sample-level data")
                return pd.DataFrame()
            
            if not df.empty:
                print(f"    ✓ Found {len(df)} samples, {df['CH_Type'].nunique()} unique CH types")
                return df
            else:
                print("    ⚠️ No CHTyper data found in HTML")
                return pd.DataFrame()
            
        except Exception as e:
            print(f"    ✗ Error parsing CHTyper HTML: {e}")
            import traceback
            traceback.print_exc()
            return pd.DataFrame()
    
    def parse_abricate_html(self, file_path: Path) -> Dict[str, pd.DataFrame]:
        """Parse ABRicate/AMRfinder HTML with dual tables"""
        print(f"  📄 Parsing Database HTML: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            
            # Determine database name from file name
            db_name = file_path.stem.lower()
            
            # Clean up the database name
            db_name = db_name.replace('ecoli_', '').replace('summary_report', '').replace('_', '').strip()
            
            # Handle specific known databases
            if 'amrfinder' in db_name:
                db_name = 'amrfinder'
            elif 'vfdb' in db_name:
                db_name = 'vfdb'
            elif 'card' in db_name:
                db_name = 'card'
            elif 'resfinder' in db_name:
                db_name = 'resfinder'
            elif 'plasmidfinder' in db_name:
                db_name = 'plasmidfinder'
            elif 'argannot' in db_name:
                db_name = 'argannot'
            elif 'megares' in db_name:
                db_name = 'megares'
            elif 'ecoh' in db_name:
                db_name = 'ecoh'
            elif 'ecolivf' in db_name or 'vf' in db_name:
                db_name = 'ecoli_vf'
            elif 'ncbi' in db_name:
                db_name = 'ncbi'
            
            # Find all tables
            tables = soup.find_all('table')
            if len(tables) < 2:
                print(f"    ⚠️ Only {len(tables)} table(s) found, expected at least 2")
                return {'database': db_name, 'genes_by_genome': pd.DataFrame(), 'gene_frequency': pd.DataFrame()}
            
            # Parse first table: Genes by Genome
            df1 = self._parse_html_table(tables[0])
            if not df1.empty:
                # Normalize column names
                df1.columns = [col.strip() for col in df1.columns]
                
                # Identify genome column
                genome_col = None
                for col in df1.columns:
                    if 'genome' in col.lower() or 'sample' in col.lower():
                        genome_col = col
                        break
                
                if genome_col:
                    df1['Sample'] = df1[genome_col].apply(self.normalize_sample_name)
                
                # Identify gene count column
                count_col = None
                for col in df1.columns:
                    if 'count' in col.lower() or 'number' in col.lower():
                        count_col = col
                        break
                
                if count_col:
                    df1['Gene_Count'] = pd.to_numeric(df1[count_col], errors='coerce').fillna(0)
                else:
                    # If no count column, try to calculate from other columns
                    for col in df1.columns:
                        if 'gene' in col.lower() and 'detected' in col.lower():
                            # Try to parse gene lists
                            df1['Gene_Count'] = df1[col].apply(lambda x: len(str(x).split(',')) if pd.notna(x) else 0)
                            break
                
                # Keep essential columns
                keep_cols = ['Sample', 'Gene_Count']
                if 'Genes Detected' in df1.columns:
                    df1['Genes'] = df1['Genes Detected']
                    keep_cols.append('Genes')
                
                df1 = df1[[col for col in keep_cols if col in df1.columns]]
            
            # Parse second table: Gene Frequency
            df2 = self._parse_html_table(tables[1])
            if not df2.empty:
                # Normalize column names
                df2.columns = [col.strip() for col in df2.columns]
                
                # Identify gene column
                gene_col = None
                for col in df2.columns:
                    if 'gene' in col.lower():
                        gene_col = col
                        break
                
                if gene_col:
                    df2['Gene'] = df2[gene_col].astype(str).str.strip()
                
                # Identify frequency column
                freq_col = None
                for col in df2.columns:
                    if 'frequency' in col.lower() or 'count' in col.lower():
                        freq_col = col
                        break
                
                if freq_col:
                    # Extract count from frequency string like "30 (100.0%)"
                    df2['Count'] = df2[freq_col].apply(self._extract_count_from_frequency)
                else:
                    # If no frequency column, try to extract from any column
                    for col in df2.columns:
                        if any(str(x).isdigit() for x in df2[col].head()):
                            try:
                                df2['Count'] = pd.to_numeric(df2[col], errors='coerce').fillna(0)
                                break
                            except:
                                pass
                
                # If still no Count column, create one with 1s
                if 'Count' not in df2.columns:
                    df2['Count'] = 1
                
                # Keep essential columns
                keep_cols = ['Gene', 'Count']
                if 'Frequency' in df2.columns:
                    df2['Frequency_Text'] = df2['Frequency']
                    keep_cols.append('Frequency_Text')
                
                df2 = df2[[col for col in keep_cols if col in df2.columns]]
            
            result = {
                'database': db_name,
                'genes_by_genome': df1,
                'gene_frequency': df2
            }
            
            print(f"    ✓ {db_name.upper()}: {len(df1)} genome entries, {len(df2)} gene entries")
            return result
            
        except Exception as e:
            print(f"    ✗ Error parsing database HTML: {e}")
            import traceback
            traceback.print_exc()
            return {'database': 'unknown', 'genes_by_genome': pd.DataFrame(), 'gene_frequency': pd.DataFrame()}
    
    def _parse_html_table(self, table) -> pd.DataFrame:
        """Parse a single HTML table"""
        try:
            # Parse headers
            headers = []
            header_row = table.find('tr')
            if header_row:
                for th in header_row.find_all(['th', 'td']):
                    headers.append(th.get_text().strip())
            
            # Parse rows
            rows = []
            for tr in table.find_all('tr')[1:]:  # Skip header
                cols = tr.find_all(['td', 'th'])
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    rows.append(row_data)
            
            # Create DataFrame
            if len(headers) > 0 and len(rows) > 0:
                # Ensure headers match columns
                if len(headers) > len(rows[0]):
                    headers = headers[:len(rows[0])]
                elif len(headers) < len(rows[0]):
                    headers = headers + [f'Col{i}' for i in range(len(headers), len(rows[0]))]
                
                return pd.DataFrame(rows, columns=headers)
            
            return pd.DataFrame()
            
        except Exception as e:
            print(f"    ⚠️ Error parsing table: {e}")
            return pd.DataFrame()
    
    def _extract_count_from_frequency(self, freq_str: str) -> int:
        """Extract count from frequency string like '30 (100.0%)'"""
        try:
            if pd.isna(freq_str):
                return 0
            
            freq_str = str(freq_str)
            # Look for numbers in the string (first number before space or parenthesis)
            match = re.search(r'(\d+)', freq_str)
            if match:
                return int(match.group(1))
            
            return 0
        except:
            return 0


class UltimateVisualizer:
    """Generates ALL requested visualizations from HTML data"""
    
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        
        # Create subdirectories
        self.subdirs = {
            'png': self.output_dir / 'PNG',
            'pdf': self.output_dir / 'PDF',
            'svg': self.output_dir / 'SVG',
            'data': self.output_dir / 'DATA'
        }
        
        for subdir in self.subdirs.values():
            subdir.mkdir(parents=True, exist_ok=True)
        
        # Set color cycle for consistent coloring
        plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Set3.colors)
    
    def _save_figure(self, fig, name: str, formats: List[str] = ['png', 'pdf', 'svg']):
        """Save figure in multiple formats"""
        for fmt in formats:
            if fmt == 'png':
                fig.savefig(self.subdirs['png'] / f"{name}.png", dpi=300, bbox_inches='tight')
            elif fmt == 'pdf':
                fig.savefig(self.subdirs['pdf'] / f"{name}.pdf", bbox_inches='tight')
            elif fmt == 'svg':
                fig.savefig(self.subdirs['svg'] / f"{name}.svg", bbox_inches='tight')
        
        plt.close(fig)
    
    def _get_colors(self, n_colors: int, palette: str = 'categorical') -> List[str]:
        """Get colors from palette"""
        if palette == 'mlst':
            return COLOR_PALETTES['mlst'][:n_colors]
        elif palette == 'serotype':
            return COLOR_PALETTES['serotype'][:n_colors]
        elif palette == 'phylogroup':
            return COLOR_PALETTES['phylogroup'][:n_colors]
        elif palette == 'chtype':
            return COLOR_PALETTES['chtype'][:n_colors]
        elif palette == 'database':
            return COLOR_PALETTES['database'][:n_colors]
        else:
            # Use tab20c colormap
            return plt.cm.tab20c(np.linspace(0, 1, n_colors))
    
    def create_complete_distribution(self, data: pd.DataFrame, category_column: str, 
                                     title: str, output_prefix: str):
        """
        Create BOTH pie and bar charts for complete distribution
        Shows ALL categories, no truncation
        """
        print(f"  📊 Creating distribution plots for {title}...")
        
        if data.empty or category_column not in data.columns:
            print(f"    ⚠️ No data for {title}")
            return
        
        # Calculate distribution
        distribution = data[category_column].value_counts().reset_index()
        distribution.columns = [category_column, 'Count']
        distribution['Percentage'] = (distribution['Count'] / distribution['Count'].sum() * 100).round(2)
        
        # Sort by count
        distribution = distribution.sort_values('Count', ascending=False)
        
        # Save data
        distribution.to_csv(self.subdirs['data'] / f"{output_prefix}_distribution.csv", index=False)
        
        # Get colors
        n_categories = len(distribution)
        colors = self._get_colors(n_categories, output_prefix.split('_')[0])
        
        # ========== PIE CHART ==========
        fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        
        # Create pie chart
        wedges, texts, autotexts = ax1.pie(
            distribution['Count'],
            labels=distribution[category_column],
            colors=colors,
            autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100.*sum(distribution["Count"]))})',
            startangle=90,
            textprops={'fontsize': 9}
        )
        
        ax1.set_title(f'{title} - Pie Chart\nTotal: {distribution["Count"].sum()} samples', 
                     fontsize=14, fontweight='bold', pad=20)
        
        # Improve label readability
        for text in texts:
            text.set_fontsize(8)
            text.set_fontweight('normal')
        for autotext in autotexts:
            autotext.set_fontsize(8)
            autotext.set_fontweight('bold')
            autotext.set_color('white')
        
        # ========== BAR CHART ==========
        # Horizontal bar chart
        y_pos = np.arange(len(distribution))
        bars = ax2.barh(y_pos, distribution['Count'], color=colors, edgecolor='black', linewidth=0.5)
        ax2.set_yticks(y_pos)
        
        # Create y-tick labels with counts
        y_labels = [f"{row[category_column]}\n(n={row['Count']}, {row['Percentage']}%)" 
                   for _, row in distribution.iterrows()]
        ax2.set_yticklabels(y_labels, fontsize=9)
        ax2.invert_yaxis()  # Highest count on top
        
        # Add count labels on bars
        for bar, count, pct in zip(bars, distribution['Count'], distribution['Percentage']):
            width = bar.get_width()
            ax2.text(width + max(distribution['Count']) * 0.01, bar.get_y() + bar.get_height()/2,
                    f'{count} ({pct}%)', va='center', fontsize=9, fontweight='bold')
        
        ax2.set_xlabel('Number of Samples', fontsize=12)
        ax2.set_title(f'{title} - Bar Chart\nTotal: {distribution["Count"].sum()} samples',
                     fontsize=14, fontweight='bold', pad=20)
        
        # Add grid
        ax2.grid(True, alpha=0.3, linestyle='--')
        
        plt.suptitle(title, fontsize=16, fontweight='bold', y=0.98)
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig1, f"{output_prefix}_complete")
        print(f"    ✓ Created {output_prefix}_complete.[png/pdf/svg]")
    
    def create_stacked_combinations(self, mlst_data: pd.DataFrame, 
                                    serotype_data: pd.DataFrame = None,
                                    phylogroup_data: pd.DataFrame = None):
        """
        Create 3 stacked bar plots:
        1. MLST × Serotype
        2. MLST × Phylogroup  
        3. MLST × Serotype × Phylogroup
        """
        print("\n📊 Creating stacked combination plots...")
        
        if mlst_data.empty:
            print("  ⚠️ No MLST data for stacked plots")
            return
        
        # Merge all available data
        merged_data = mlst_data.copy()
        
        if not serotype_data.empty:
            merged_data = pd.merge(merged_data, serotype_data[['Sample', 'Serotype']], 
                                  on='Sample', how='left')
        
        if not phylogroup_data.empty:
            merged_data = pd.merge(merged_data, phylogroup_data[['Sample', 'Phylogroup']], 
                                  on='Sample', how='left')
        
        # Save merged data
        merged_data.to_csv(self.subdirs['data'] / "merged_typing_data.csv", index=False)
        
        # Plot 1: MLST × Serotype
        if 'Serotype' in merged_data.columns:
            self._create_stacked_plot(
                merged_data, 
                primary_col='ST',
                secondary_col='Serotype',
                title='MLST × Serotype Relationship',
                filename='stacked_mlst_serotype'
            )
        
        # Plot 2: MLST × Phylogroup
        if 'Phylogroup' in merged_data.columns:
            self._create_stacked_plot(
                merged_data,
                primary_col='ST',
                secondary_col='Phylogroup',
                title='MLST × Phylogroup Relationship',
                filename='stacked_mlst_phylogroup'
            )
        
        # Plot 3: MLST × Serotype × Phylogroup (triple)
        if 'Serotype' in merged_data.columns and 'Phylogroup' in merged_data.columns:
            self._create_triple_stacked_plot(merged_data)
    
    def _create_stacked_plot(self, data: pd.DataFrame, primary_col: str, 
                            secondary_col: str, title: str, filename: str):
        """Create a single stacked bar plot"""
        print(f"  Creating {filename}...")
        
        # Create cross-tabulation
        cross_tab = pd.crosstab(data[primary_col], data[secondary_col])
        
        # Remove 'nan' if exists
        if np.nan in cross_tab.columns:
            cross_tab = cross_tab.drop(columns=np.nan)
        
        # Sort by total count
        cross_tab = cross_tab.loc[cross_tab.sum(axis=1).sort_values(ascending=False).index]
        
        # Sort columns by total
        cross_tab = cross_tab[cross_tab.sum().sort_values(ascending=False).index]
        
        # Plot
        fig, ax = plt.subplots(figsize=(16, 10))
        
        # Get colors
        n_colors = len(cross_tab.columns)
        colors = self._get_colors(n_colors, 'categorical')
        
        # Create stacked bar
        bottom = np.zeros(len(cross_tab))
        for i, col in enumerate(cross_tab.columns):
            ax.bar(cross_tab.index, cross_tab[col], bottom=bottom, 
                  label=str(col), color=colors[i], edgecolor='black', linewidth=0.5)
            bottom += cross_tab[col]
        
        ax.set_xlabel(primary_col.upper(), fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
        ax.set_title(f'{title}\nTotal Samples: {len(data)}', 
                    fontsize=14, fontweight='bold', pad=20)
        
        # Rotate x-tick labels for readability
        plt.xticks(rotation=45, ha='right', fontsize=10)
        plt.yticks(fontsize=10)
        
        # Add count labels on top of each bar
        for i, total in enumerate(cross_tab.sum(axis=1)):
            ax.text(i, total + max(cross_tab.sum(axis=1)) * 0.01, f'{int(total)}', 
                   ha='center', va='bottom', fontsize=9, fontweight='bold')
        
        # Add legend
        ax.legend(title=secondary_col, bbox_to_anchor=(1.05, 1), loc='upper left', 
                 fontsize=9, title_fontsize=10)
        
        # Add grid
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, filename)
        print(f"    ✓ Created {filename}.[png/pdf/svg]")
    
    def _create_triple_stacked_plot(self, data: pd.DataFrame):
        """Create triple stacked plot (MLST × Serotype × Phylogroup)"""
        print("  Creating stacked_mlst_serotype_phylogroup...")
        
        # Create grouped data
        data['ST_Serotype'] = data['ST'] + ' × ' + data['Serotype']
        grouped = data.groupby(['ST_Serotype', 'Phylogroup']).size().reset_index(name='count')
        
        # Pivot for plotting
        pivot_data = grouped.pivot(index='ST_Serotype', columns='Phylogroup', 
                                  values='count').fillna(0)
        
        # Sort by total count
        pivot_data = pivot_data.loc[pivot_data.sum(axis=1).sort_values(ascending=False).index]
        
        # Get phylogroups sorted by frequency
        phylogroups = pivot_data.sum().sort_values(ascending=False).index.tolist()
        
        # Plot
        fig, ax = plt.subplots(figsize=(20, 12))
        
        # Get colors for phylogroups
        colors = self._get_colors(len(phylogroups), 'phylogroup')
        
        # Create stacked bars
        bottom = np.zeros(len(pivot_data))
        for i, phylogroup in enumerate(phylogroups):
            if phylogroup in pivot_data.columns:
                ax.bar(range(len(pivot_data)), pivot_data[phylogroup], bottom=bottom,
                      label=phylogroup, color=colors[i], edgecolor='black', linewidth=0.5)
                bottom += pivot_data[phylogroup]
        
        # Customize plot
        ax.set_xlabel('MLST × Serotype Combination', fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
        ax.set_title('MLST × Serotype × Phylogroup Relationship\nTriple Typing Analysis', 
                    fontsize=16, fontweight='bold', pad=20)
        
        # Set x-tick labels
        tick_labels = [idx.replace(' × ', '\n×\n') for idx in pivot_data.index]
        ax.set_xticks(range(len(pivot_data)))
        ax.set_xticklabels(tick_labels, rotation=45, ha='right', fontsize=8)
        
        # Add total counts on top of bars
        for i, total in enumerate(pivot_data.sum(axis=1)):
            ax.text(i, total + max(pivot_data.sum(axis=1)) * 0.01, f'{int(total)}', 
                   ha='center', va='bottom', fontsize=8, fontweight='bold')
        
        # Add legend
        ax.legend(title='Phylogroup', bbox_to_anchor=(1.05, 1), loc='upper left', 
                 fontsize=9, title_fontsize=10)
        
        # Add grid
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, "stacked_mlst_serotype_phylogroup")
        print("    ✓ Created stacked_mlst_serotype_phylogroup.[png/pdf/svg]")
    
    def create_database_comparison(self, databases_data: Dict[str, Dict[str, pd.DataFrame]]):
        """
        Create box/violin plots comparing gene hits across databases
        """
        print("\n📊 Creating database comparison plots...")
        
        if not databases_data:
            print("  ⚠️ No database data for comparison")
            return
        
        # Prepare data for box plot
        box_data = []
        database_names = []
        
        for db_name, db_data in databases_data.items():
            if 'genes_by_genome' in db_data and not db_data['genes_by_genome'].empty:
                df = db_data['genes_by_genome']
                if 'Gene_Count' in df.columns:
                    for count in df['Gene_Count']:
                        box_data.append(float(count))
                        database_names.append(db_name)
        
        if not box_data:
            print("  ⚠️ No gene count data found")
            return
        
        # Create DataFrame
        comparison_df = pd.DataFrame({
            'Database': database_names,
            'Gene_Count': box_data
        })
        
        # Save comparison data
        comparison_df.to_csv(self.subdirs['data'] / "database_comparison.csv", index=False)
        
        # Calculate statistics
        stats_df = comparison_df.groupby('Database')['Gene_Count'].agg([
            'count', 'mean', 'median', 'std', 'min', 'max'
        ]).round(2)
        stats_df.to_csv(self.subdirs['data'] / "database_statistics.csv")
        
        # ========== BOX PLOT ==========
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        
        # Get colors
        n_databases = len(stats_df)
        colors = self._get_colors(n_databases, 'database')
        
        # Box plot
        box_data_by_db = [comparison_df[comparison_df['Database'] == db]['Gene_Count'].values 
                         for db in stats_df.index]
        
        bp = ax1.boxplot(box_data_by_db, patch_artist=True, showfliers=True, 
                        labels=stats_df.index, medianprops={'color': 'black', 'linewidth': 2})
        
        # Color boxes
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax1.set_xticklabels(stats_df.index, rotation=45, ha='right', fontsize=10)
        ax1.set_ylabel('Number of Genes Detected', fontsize=12, fontweight='bold')
        ax1.set_title('Gene Detection by Database - Box Plot', 
                     fontsize=14, fontweight='bold', pad=20)
        
        # Add grid
        ax1.grid(True, alpha=0.3, linestyle='--')
        
        # ========== VIOLIN PLOT ==========
        # Violin plot with swarm
        sns.violinplot(data=comparison_df, x='Database', y='Gene_Count', 
                      ax=ax2, palette=colors, inner='quartile', cut=0)
        
        # Add swarm plot for individual points
        sns.swarmplot(data=comparison_df, x='Database', y='Gene_Count', 
                     ax=ax2, color='black', alpha=0.5, size=3)
        
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right', fontsize=10)
        ax2.set_ylabel('Number of Genes Detected', fontsize=12, fontweight='bold')
        ax2.set_title('Gene Detection by Database - Violin Plot', 
                     fontsize=14, fontweight='bold', pad=20)
        
        # Add grid
        ax2.grid(True, alpha=0.3, linestyle='--')
        
        plt.suptitle('Database Comparison - Gene Detection Statistics', 
                    fontsize=16, fontweight='bold', y=0.98)
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, "database_boxplots")
        print("    ✓ Created database_boxplots.[png/pdf/svg]")
        
        # Create additional statistics table visualization
        self._create_database_statistics_table(stats_df)
    
    def _create_database_statistics_table(self, stats_df: pd.DataFrame):
        """Create a visual table of database statistics"""
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.axis('tight')
        ax.axis('off')
        
        # Create table
        table_data = []
        for db in stats_df.index:
            row = [
                db,
                f"{stats_df.loc[db, 'count']}",
                f"{stats_df.loc[db, 'mean']:.2f}",
                f"{stats_df.loc[db, 'median']:.2f}",
                f"{stats_df.loc[db, 'std']:.2f}",
                f"{stats_df.loc[db, 'min']:.0f}",
                f"{stats_df.loc[db, 'max']:.0f}"
            ]
            table_data.append(row)
        
        # Create table - FIXED: proper column count
        n_columns = len(['Database', 'N', 'Mean', 'Median', 'Std', 'Min', 'Max'])
        table = ax.table(cellText=table_data,
                        colLabels=['Database', 'N', 'Mean', 'Median', 'Std', 'Min', 'Max'],
                        cellLoc='center',
                        loc='center',
                        colWidths=[0.15] * n_columns)
        
        # Style table
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 2)
        
        # Color header - FIXED: use proper column count
        for i in range(n_columns):
            table[(0, i)].set_facecolor('#4C78A8')
            table[(0, i)].set_text_props(weight='bold', color='white')
        
        # Alternate row colors - FIXED: use proper column count
        for i in range(1, len(table_data) + 1):
            if i % 2 == 0:
                for j in range(n_columns):
                    table[(i, j)].set_facecolor('#f0f0f0')
        
        ax.set_title('Database Statistics Summary', fontsize=14, fontweight='bold', pad=20)
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, "database_statistics_table")
        print("    ✓ Created database_statistics_table.[png/pdf/svg]")
    
    def create_gene_distributions(self, databases_data: Dict[str, Dict[str, pd.DataFrame]]):
        """
        Create statistical distribution plots from gene frequency tables
        """
        print("\n📊 Creating gene frequency distribution plots...")
        
        if not databases_data:
            print("  ⚠️ No database data for gene distributions")
            return
        
        # Extract frequency data
        freq_data = {}
        for db_name, db_data in databases_data.items():
            if 'gene_frequency' in db_data and not db_data['gene_frequency'].empty:
                df = db_data['gene_frequency']
                if 'Count' in df.columns:
                    # Convert to numeric and drop NaN
                    counts = pd.to_numeric(df['Count'], errors='coerce').dropna().values
                    if len(counts) > 0:
                        freq_data[db_name] = counts
        
        if not freq_data:
            print("  ⚠️ No frequency data found")
            return
        
        # Create subplots
        n_databases = len(freq_data)
        n_cols = min(3, n_databases)
        n_rows = (n_databases + n_cols - 1) // n_cols
        
        if n_databases == 0:
            print("  ⚠️ No valid frequency data to plot")
            return
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
        if n_databases == 1:
            axes = np.array([axes])
        axes = axes.flatten()
        
        # Get colors
        colors = self._get_colors(n_databases, 'database')
        
        for i, (db_name, counts) in enumerate(freq_data.items()):
            if i >= len(axes):
                break
                
            ax = axes[i]
            
            # Plot histogram with KDE
            n, bins, patches = ax.hist(counts, bins=30, alpha=0.7, density=True,
                                      color=colors[i], edgecolor='black', linewidth=0.5)
            
            try:
                # Plot KDE if we have enough data points
                if len(counts) > 1:
                    kde = stats.gaussian_kde(counts)
                    x_range = np.linspace(min(counts), max(counts), 1000)
                    ax.plot(x_range, kde(x_range), 'r-', linewidth=2, label='KDE')
            except Exception as e:
                print(f"    ⚠️ Could not compute KDE for {db_name}: {e}")
            
            # Add statistics lines
            mean_val = np.mean(counts)
            median_val = np.median(counts)
            
            ax.axvline(mean_val, color='blue', linestyle='--', linewidth=2, 
                      label=f'Mean: {mean_val:.2f}')
            ax.axvline(median_val, color='green', linestyle='--', linewidth=2,
                      label=f'Median: {median_val:.2f}')
            
            # Fill area under curve if KDE was computed
            if len(counts) > 1:
                try:
                    ax.fill_between(x_range, kde(x_range), alpha=0.3, color=colors[i])
                except:
                    pass
            
            ax.set_xlabel('Gene Frequency (Count)', fontsize=11)
            ax.set_ylabel('Density', fontsize=11)
            ax.set_title(f'{db_name.upper()}\nGene Frequency Distribution', 
                        fontsize=12, fontweight='bold', pad=15)
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3, linestyle='--')
            
            # Add text box with statistics
            stats_text = f'n = {len(counts)}\nMean = {mean_val:.2f}\nMedian = {median_val:.2f}\nStd = {np.std(counts):.2f}'
            ax.text(0.95, 0.95, stats_text, transform=ax.transAxes, fontsize=9,
                   verticalalignment='top', horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Hide unused subplots
        for j in range(i + 1, len(axes)):
            axes[j].axis('off')
        
        plt.suptitle('Gene Frequency Distributions Across Databases', 
                    fontsize=16, fontweight='bold', y=0.98)
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, "gene_distributions")
        print("    ✓ Created gene_distributions.[png/pdf/svg]")
        
        # Create cumulative distribution plot
        self._create_cumulative_distribution(freq_data)
    
    def _create_cumulative_distribution(self, freq_data: Dict[str, np.ndarray]):
        """Create cumulative distribution function plot"""
        if not freq_data:
            return
            
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Get colors
        colors = self._get_colors(len(freq_data), 'database')
        
        for (db_name, counts), color in zip(freq_data.items(), colors):
            # Sort data for CDF
            sorted_counts = np.sort(counts)
            cdf = np.arange(1, len(sorted_counts) + 1) / len(sorted_counts)
            
            ax.plot(sorted_counts, cdf, label=db_name.upper(), color=color, linewidth=2)
        
        ax.set_xlabel('Gene Frequency (Count)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Cumulative Probability', fontsize=12, fontweight='bold')
        ax.set_title('Cumulative Distribution Function (CDF) of Gene Frequencies',
                    fontsize=14, fontweight='bold', pad=20)
        
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, linestyle='--')
        
        # Add reference lines
        for percentile in [25, 50, 75, 90]:
            ax.axhline(percentile/100, color='gray', linestyle=':', alpha=0.5, linewidth=1)
            ax.text(ax.get_xlim()[1] * 0.95, percentile/100, f'{percentile}%', 
                   ha='right', va='bottom', fontsize=9, color='gray')
        
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, "gene_cumulative_distribution")
        print("    ✓ Created gene_cumulative_distribution.[png/pdf/svg]")


class VisualizationReporter:
    """Main class to orchestrate the complete visualization pipeline"""
    
    def __init__(self, input_dir: Path = Path.cwd(), output_dir: Path = None):
        self.input_dir = Path(input_dir)
        
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            self.output_dir = self.input_dir / "ECOLI_VISUALIZATIONS"
        
        self.parser = UltimateHTMLParser()
        self.visualizer = UltimateVisualizer(self.output_dir)
        
        # Data storage
        self.mlst_data = pd.DataFrame()
        self.serotype_data = pd.DataFrame()
        self.phylogroup_data = pd.DataFrame()
        self.chtype_data = pd.DataFrame()
        self.databases_data = {}
    
    def run_pipeline(self):
        """Run the complete visualization pipeline"""
        print("=" * 70)
        print("🧬 GENIUS E.COLI VISUALIZATION MODULE - ULTIMATE HTML EDITION")
        print("=" * 70)
        print(f"Input directory: {self.input_dir}")
        print(f"Output directory: {self.output_dir}")
        print("-" * 70)
        
        start_time = datetime.now()
        
        # Step 1: Parse all HTML files
        self._parse_all_files()
        
        # Step 2: Generate Category 1 plots (Complete Distributions)
        print("\n" + "=" * 70)
        print("🎯 CATEGORY 1: Complete Distribution Plots")
        print("=" * 70)
        self._generate_category1_plots()
        
        # Step 3: Generate Category 2 plots (Stacked Combinations)
        print("\n" + "=" * 70)
        print("🎯 CATEGORY 2: Stacked Combination Plots")
        print("=" * 70)
        self._generate_category2_plots()
        
        # Step 4: Generate Category 3 plots (Database Statistics)
        print("\n" + "=" * 70)
        print("🎯 CATEGORY 3: Database Statistical Plots")
        print("=" * 70)
        self._generate_category3_plots()
        
        # Step 5: Generate summary report
        self._generate_summary_report(start_time)
        
        print("\n" + "=" * 70)
        print("✅ VISUALIZATION PIPELINE COMPLETED SUCCESSFULLY!")
        print("=" * 70)
    
    def _parse_all_files(self):
        """Parse all HTML files in the input directory"""
        print("\n📁 Parsing HTML files...")
        
        # Parse MLST
        mlst_files = list(self.input_dir.glob("*mlst*.html"))
        if mlst_files:
            self.mlst_data = self.parser.parse_mlst_html(mlst_files[0])
        
        # Parse Serotype
        serotype_files = list(self.input_dir.glob("*serotype*.html"))
        if serotype_files:
            self.serotype_data = self.parser.parse_serotype_html(serotype_files[0])
        
        # Parse Phylogroup
        phylogroup_files = list(self.input_dir.glob("*phylogroup*.html"))
        if phylogroup_files:
            self.phylogroup_data = self.parser.parse_phylogroup_html(phylogroup_files[0])
        
        # Parse CHTyper
        chtyper_files = list(self.input_dir.glob("*chtyper*.html"))
        if chtyper_files:
            self.chtype_data = self.parser.parse_chtyper_html(chtyper_files[0])
        
        # Parse all database files
        database_patterns = [
            "*amrfinder*.html",
            "*vfdb*.html", 
            "*card*.html",
            "*resfinder*.html",
            "*plasmidfinder*.html",
            "*argannot*.html",
            "*megares*.html",
            "*ecoh*.html",
            "*ecoli_vf*.html",
            "*ncbi*.html"
        ]
        
        # Track parsed files to avoid duplicates
        parsed_files = set()
        
        for pattern in database_patterns:
            for db_file in self.input_dir.glob(pattern):
                # Skip files already parsed as MLST, serotype, etc.
                if (db_file in mlst_files or db_file in serotype_files or 
                    db_file in phylogroup_files or db_file in chtyper_files):
                    continue
                    
                # Skip already parsed files (case-insensitive)
                db_file_lower = str(db_file).lower()
                if any(db_file_lower in pf for pf in parsed_files):
                    continue
                    
                db_data = self.parser.parse_abricate_html(db_file)
                db_name = db_data.get('database', 'unknown')
                
                # Only add if we have valid data
                if (not db_data['genes_by_genome'].empty or 
                    not db_data['gene_frequency'].empty):
                    self.databases_data[db_name] = db_data
                    parsed_files.add(db_file_lower)
    
    def _generate_category1_plots(self):
        """Generate Category 1: Complete Distribution Plots"""
        # MLST Distribution
        if not self.mlst_data.empty and 'ST' in self.mlst_data.columns:
            self.visualizer.create_complete_distribution(
                self.mlst_data, 'ST',
                "MLST Sequence Types - Complete Distribution",
                "mlst_distribution"
            )
        
        # Serotype Distribution
        if not self.serotype_data.empty and 'Serotype' in self.serotype_data.columns:
            self.visualizer.create_complete_distribution(
                self.serotype_data, 'Serotype',
                "Serotype Distribution - Complete Analysis",
                "serotype_distribution"
            )
        
        # Phylogroup Distribution
        if not self.phylogroup_data.empty and 'Phylogroup' in self.phylogroup_data.columns:
            self.visualizer.create_complete_distribution(
                self.phylogroup_data, 'Phylogroup',
                "Clermont Phylogroup Distribution",
                "phylogroup_distribution"
            )
        
        # CH Type Distribution
        if not self.chtype_data.empty and 'CH_Type' in self.chtype_data.columns:
            self.visualizer.create_complete_distribution(
                self.chtype_data, 'CH_Type',
                "CH Type (FumC:FimH) Distribution",
                "chtype_distribution"
            )
    
    def _generate_category2_plots(self):
        """Generate Category 2: Stacked Combination Plots"""
        self.visualizer.create_stacked_combinations(
            self.mlst_data,
            self.serotype_data,
            self.phylogroup_data
        )
    
    def _generate_category3_plots(self):
        """Generate Category 3: Database Statistical Plots"""
        if self.databases_data:
            # Create database comparison plots
            self.visualizer.create_database_comparison(self.databases_data)
            
            # Create gene frequency distributions
            self.visualizer.create_gene_distributions(self.databases_data)
        else:
            print("  ⚠️ No database data available for statistical plots")
    
    def _generate_summary_report(self, start_time):
        """Generate a comprehensive summary report"""
        end_time = datetime.now()
        duration = end_time - start_time
        
        report_path = self.output_dir / "visualization_report.txt"
        
        with open(report_path, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("🧬 GENIUS E.COLI VISUALIZATION MODULE - SUMMARY REPORT\n")
            f.write("=" * 70 + "\n\n")
            
            f.write("📋 EXECUTION DETAILS\n")
            f.write("-" * 40 + "\n")
            f.write(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Duration: {duration}\n")
            f.write(f"Input directory: {self.input_dir}\n")
            f.write(f"Output directory: {self.output_dir}\n\n")
            
            f.write("📊 DATA SUMMARY\n")
            f.write("-" * 40 + "\n")
            
            # MLST summary
            if not self.mlst_data.empty:
                unique_sts = self.mlst_data['ST'].nunique() if 'ST' in self.mlst_data.columns else 0
                f.write(f"MLST Data: {len(self.mlst_data)} samples, {unique_sts} unique STs\n")
            
            # Serotype summary
            if not self.serotype_data.empty:
                unique_seros = self.serotype_data['Serotype'].nunique() if 'Serotype' in self.serotype_data.columns else 0
                f.write(f"Serotype Data: {len(self.serotype_data)} samples, {unique_seros} unique serotypes\n")
            
            # Phylogroup summary
            if not self.phylogroup_data.empty:
                unique_phylos = self.phylogroup_data['Phylogroup'].nunique() if 'Phylogroup' in self.phylogroup_data.columns else 0
                f.write(f"Phylogroup Data: {len(self.phylogroup_data)} samples, {unique_phylos} unique phylogroups\n")
            
            # CH Type summary
            if not self.chtype_data.empty:
                unique_chtypes = self.chtype_data['CH_Type'].nunique() if 'CH_Type' in self.chtype_data.columns else 0
                f.write(f"CH Type Data: {len(self.chtype_data)} samples, {unique_chtypes} unique CH types\n")
            
            # Database summary
            f.write(f"Databases Analyzed: {len(self.databases_data)}\n")
            for db_name, db_data in self.databases_data.items():
                gbg_len = len(db_data.get('genes_by_genome', pd.DataFrame()))
                gf_len = len(db_data.get('gene_frequency', pd.DataFrame()))
                f.write(f"  • {db_name.upper()}: {gbg_len} genomes, {gf_len} genes\n")
            
            f.write("\n📈 GENERATED VISUALIZATIONS\n")
            f.write("-" * 40 + "\n")
            
            # List generated files
            for fmt in ['PNG', 'PDF', 'SVG', 'DATA']:
                fmt_dir = self.output_dir / fmt
                if fmt_dir.exists():
                    files = list(fmt_dir.glob("*"))
                    f.write(f"\n{fmt} Files ({len(files)}):\n")
                    for file in sorted(files):
                        f.write(f"  • {file.name}\n")
            
            f.write("\n" + "=" * 70 + "\n")
            f.write("✅ REPORT COMPLETE\n")
            f.write("=" * 70 + "\n")
        
        print(f"\n📋 Summary report saved: {report_path}")
        
        # Print summary to console
        print("\n📋 FINAL SUMMARY:")
        print("-" * 40)
        print(f"Total time: {duration}")
        print(f"Total plots generated: {len(list(self.output_dir.glob('**/*.png')))}")
        print(f"Output directory: {self.output_dir}")


def main():
    """Main function to run the visualization pipeline"""
    parser = argparse.ArgumentParser(
        description="GENIUS E.Coli Visualization Module - Parse HTML reports and generate publication-quality visualizations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python visualization_reporter.py
  python visualization_reporter.py --input /path/to/html/files --output /path/to/results
  python visualization_reporter.py --force
        
This script will:
  1. Parse all HTML files in the input directory
  2. Generate complete distribution plots for MLST, Serotype, Phylogroup, and CH Types
  3. Create stacked combination plots
  4. Generate database comparison and statistical plots
  5. Save results in PNG, PDF, SVG, and CSV formats
        """
    )
    
    parser.add_argument("--input", type=str, default=".",
                       help="Directory containing HTML files (default: current directory)")
    parser.add_argument("--output", type=str, default=None,
                       help="Output directory (default: ECOLI_VISUALIZATIONS in input directory)")
    parser.add_argument("--force", action="store_true",
                       help="Overwrite existing output directory")
    
    args = parser.parse_args()
    
    # Create reporter
    reporter = VisualizationReporter(
        input_dir=Path(args.input),
        output_dir=Path(args.output) if args.output else None
    )
    
    # Run pipeline
    try:
        reporter.run_pipeline()
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()