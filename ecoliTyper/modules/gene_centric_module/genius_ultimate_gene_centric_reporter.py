#!/usr/bin/env python3
"""
GENIUS E. COLI ULTIMATE REPORTER v4.0.0
========================================
Interactive, gene-centric HTML report generator for EcoliTyper outputs.

Combines a master typing TSV with module HTML reports (AMR, virulence,
plasmids, BACMET, mutations, MGE, QC) into a single multi-tab report with
live filtering, grouping, comparison, cluster detection, and export.

Sections
--------
Summary · Sample Overview · FASTA QC · MLST · Serotype · CH Type ·
Phylogroup · Pathotype · AMR · Virulence · Plasmids & Serotyping ·
BACMET · Mutations · MGE · Compare & Cluster · Patterns ·
Database Metrics · Citation · Funding · Call to Action · Export

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School, Department of Medical Biochemistry
License: MIT
"""

import sys
import json
import re
import io
import argparse
import pandas as pd
from pathlib import Path
from typing import Dict, List, Any
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

from bs4 import BeautifulSoup


# =============================================================================
# PARSER
# =============================================================================
class EcoliParser:
    """Parses EcoliTyper module outputs into structured dictionaries."""

    def __init__(self):
        self.abricate_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ecoh', 'ecoli_vf', 'bacmet2',
        ]

    def normalize_sample_id(self, sample_id: str) -> str:
        """Strip extension and directory prefix."""
        s = str(sample_id).strip()
        for ext in ('.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv'):
            if s.endswith(ext):
                s = s[:-len(ext)]
        if '/' in s or '\\' in s:
            s = Path(s).name
        return s.strip()

    def clean_st(self, st_str: str) -> str:
        if not st_str or st_str in ('ND', 'nan', 'None'):
            return 'ND'
        st_str = re.sub(r'^(ST)+', 'ST', str(st_str).strip(), flags=re.I)
        if not st_str.upper().startswith('ST'):
            st_str = 'ST' + st_str
        return st_str

    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if not tables or table_index >= len(tables):
                return pd.DataFrame()
            rows = tables[table_index].find_all('tr')
            if not rows:
                return pd.DataFrame()
            headers = [th.get_text().strip() for th in rows[0].find_all(['th', 'td'])]
            data = []
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                if cols:
                    row_data = [c.get_text().strip() for c in cols]
                    if len(row_data) == len(headers):
                        data.append(row_data)
            return pd.DataFrame(data, columns=headers) if data else pd.DataFrame()
        except Exception as e:
            print(f"    Table parsing error: {e}")
            return pd.DataFrame()

    def load_master_tsv(self, path: Path) -> Dict[str, Dict]:
        print(f"  Loading master typing TSV: {path.name}")
        if not path.exists():
            print(f"    File not found: {path}")
            return {}
        try:
            df = pd.read_csv(path, sep='\t', dtype=str).fillna('ND')
            results = {}
            for _, row in df.iterrows():
                sample = self.normalize_sample_id(row.get('Sample_ID', ''))
                if not sample:
                    continue
                results[sample] = {
                    'ST': self.clean_st(row.get('MLST_ST', 'ND')),
                    'Serotype': str(row.get('Serotype', 'ND')).strip(),
                    'O_Type': str(row.get('O_Type', 'ND')).strip(),
                    'H_Type': str(row.get('H_Type', 'ND')).strip(),
                    'Serotype_Confidence': str(row.get('Serotype_Confidence', 'ND')).strip(),
                    'Phylogroup': str(row.get('Phylogroup', 'ND')).strip(),
                    'FumC_Type': str(row.get('FumC_Type', 'ND')).strip(),
                    'FimH_Type': str(row.get('FimH_Type', 'ND')).strip(),
                    'CH_Type': str(row.get('CH_Type', 'ND')).strip(),
                    'FumC_Identity': str(row.get('FumC_Identity', 'ND')).strip(),
                    'FimH_Identity': str(row.get('FimH_Identity', 'ND')).strip(),
                    'FumC_Coverage': str(row.get('FumC_Coverage', 'ND')).strip(),
                    'FimH_Coverage': str(row.get('FimH_Coverage', 'ND')).strip(),
                    'Status': str(row.get('Status', 'ND')).strip(),
                }
            print(f"    Loaded {len(results)} samples")
            return results
        except Exception as e:
            print(f"    Error reading master TSV: {e}")
            return {}

    def parse_qc_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  Parsing FASTA QC: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty:
                return {}
            sample_col = next((c for c in df.columns
                               if 'filename' in c.lower() or 'sample' in c.lower()),
                              df.columns[0])
            results = {}
            for _, row in df.iterrows():
                sample = self.normalize_sample_id(row[sample_col])
                if not sample:
                    continue
                qc = {}
                for col in df.columns:
                    if col == sample_col:
                        continue
                    val = row[col]
                    if pd.isna(val) or val in ('', 'ND'):
                        qc[col] = 'ND'
                    else:
                        cleaned = str(val).replace('%', '').replace(',', '').strip()
                        try:
                            qc[col] = float(cleaned)
                        except Exception:
                            qc[col] = str(val)
                results[sample] = qc
            print(f"    Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    Error parsing QC: {e}")
            return {}

    def parse_amrfinder_report(self, file_path: Path, total_samples: int = 0):
        print(f"  Parsing AMRfinder: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}, {}
            gene_frequencies = {}
            df2 = self.parse_html_table(str(tables[1]), 0)
            if not df2.empty and 'Gene' in df2.columns:
                for _, row in df2.iterrows():
                    gene = str(row['Gene']).strip()
                    freq = str(row.get('Frequency', '0')).strip()
                    genomes = []
                    if 'Genomes' in df2.columns and pd.notna(row.get('Genomes')):
                        genomes = [self.normalize_sample_id(g.strip())
                                   for g in str(row['Genomes']).split(',') if g.strip()]
                    m = re.search(r'(\d+)', freq)
                    count = int(m.group(1)) if m else 0
                    pct = (count / total_samples * 100) if total_samples else 0
                    gene_frequencies[gene] = {
                        'count': count,
                        'percentage': round(pct, 2),
                        'frequency_display': f"{count} ({pct:.1f}%)",
                        'genomes': genomes,
                        'database': 'amrfinder',
                    }
            genes_by_genome = {}
            df1 = self.parse_html_table(str(tables[0]), 0)
            if not df1.empty and 'Genome' in df1.columns:
                for _, row in df1.iterrows():
                    sample = self.normalize_sample_id(row['Genome'])
                    genes = []
                    if pd.notna(row.get('Genes Detected')):
                        genes = [g.strip() for g in str(row['Genes Detected']).split(',') if g.strip()]
                    genes_by_genome[sample] = genes
            print(f"    {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return genes_by_genome, gene_frequencies
        except Exception as e:
            print(f"    Error parsing AMRfinder: {e}")
            return {}, {}

    def parse_abricate_database_report(self, file_path: Path, total_samples: int = 0):
        print(f"  Parsing ABRicate: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return 'unknown', {}, {}
            db_name = 'unknown'
            title_tag = soup.find('title')
            if title_tag:
                m = re.search(r'(\w+)\s+Database', title_tag.get_text(), re.I)
                if m:
                    db_name = m.group(1).lower()
            if db_name == 'unknown':
                for db in self.abricate_databases:
                    if db in file_path.name.lower():
                        db_name = db
                        break
            gene_frequencies = {}
            df2 = self.parse_html_table(str(tables[1]), 0)
            if not df2.empty and 'Gene' in df2.columns:
                for _, row in df2.iterrows():
                    gene = str(row['Gene']).strip()
                    freq = str(row.get('Frequency', '0')).strip()
                    genomes = []
                    if 'Genomes' in df2.columns and pd.notna(row.get('Genomes')):
                        genomes = [self.normalize_sample_id(g.strip())
                                   for g in str(row['Genomes']).split(',') if g.strip()]
                    m = re.search(r'(\d+)', freq)
                    count = int(m.group(1)) if m else 0
                    pct = (count / total_samples * 100) if total_samples else 0
                    gene_frequencies[gene] = {
                        'count': count,
                        'percentage': round(pct, 2),
                        'frequency_display': f"{count} ({pct:.1f}%)",
                        'genomes': genomes,
                        'database': db_name,
                    }
            genes_by_genome = {}
            df1 = self.parse_html_table(str(tables[0]), 0)
            if not df1.empty:
                sample_col = next((c for c in df1.columns
                                   if 'genome' in c.lower() or 'sample' in c.lower()), None)
                genes_col = next((c for c in df1.columns
                                  if 'genes' in c.lower() or 'detected' in c.lower()), None)
                if sample_col:
                    for _, row in df1.iterrows():
                        sample = self.normalize_sample_id(row[sample_col])
                        genes = []
                        if genes_col and pd.notna(row.get(genes_col)):
                            genes = [g.strip() for g in str(row[genes_col]).split(',') if g.strip()]
                        genes_by_genome[sample] = genes
            print(f"    {db_name.upper()}: {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return db_name, genes_by_genome, gene_frequencies
        except Exception as e:
            print(f"    Error parsing ABRicate report: {e}")
            return 'unknown', {}, {}

    def parse_mutation_summary_html(self, file_path: Path) -> Dict[str, Any]:
        print(f"  Parsing mutation summary: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            table = None
            for t in soup.find_all('table'):
                if t.find(string=re.compile(r'Gene', re.I)) and \
                   t.find(string=re.compile(r'Mutation|Element', re.I)):
                    table = t
                    break
            if not table:
                return {}
            header_row = (table.find('thead').find('tr') if table.find('thead')
                          else table.find('tr'))
            if not header_row:
                return {}
            headers = [c.get_text().strip() for c in header_row.find_all(['th', 'td']) if c.get_text().strip()]
            col_idx = {}
            for i, h in enumerate(headers):
                hl = h.lower()
                if 'gene' in hl: col_idx['gene'] = i
                elif 'mutation' in hl or 'element' in hl or 'variant' in hl: col_idx['mutation'] = i
                elif 'class' in hl: col_idx['class'] = i
                elif 'subclass' in hl: col_idx['subclass'] = i
                elif 'genome' in hl: col_idx['genomes'] = i
            if 'gene' not in col_idx or 'genomes' not in col_idx:
                return {}
            tbody = table.find('tbody')
            rows = tbody.find_all('tr') if tbody else table.find_all('tr')[1:]
            mutations_list = []
            genome_counts = defaultdict(int)
            for row in rows:
                cells = row.find_all('td')
                if len(cells) <= max(col_idx.values()):
                    continue
                gene = cells[col_idx['gene']].get_text().strip()
                mutation = cells[col_idx.get('mutation', 0)].get_text().strip() if 'mutation' in col_idx else ''
                class_name = cells[col_idx['class']].get_text().strip() if 'class' in col_idx else ''
                subclass = cells[col_idx['subclass']].get_text().strip() if 'subclass' in col_idx else ''
                genomes_str = cells[col_idx['genomes']].get_text().strip()
                genomes = [self.normalize_sample_id(g.strip()) for g in genomes_str.split(',') if g.strip()]
                if not genomes:
                    continue
                for g in genomes:
                    genome_counts[g] += 1
                mutations_list.append({
                    'gene': gene, 'mutation': mutation,
                    'class': class_name, 'subclass': subclass,
                    'count': len(genomes), 'genomes': genomes,
                })
            mutations_list.sort(key=lambda x: x['count'], reverse=True)
            print(f"    {len(mutations_list)} mutations across {len(genome_counts)} genomes")
            return {'mutations': mutations_list,
                    'genome_mutation_counts': dict(genome_counts)}
        except Exception as e:
            print(f"    Error parsing mutation summary: {e}")
            return {}

    def parse_mge_summary_html(self, file_path: Path) -> Dict[str, Any]:
        """Parse the mobileOG-db MGE summary HTML."""
        print(f"  Parsing MGE summary: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')

            stats = {}
            labels = ('Total Samples', 'With mobileOG hits', 'Without mobileOG hits',
                      'Total mobileOG hits', 'Key MGE-associated genes', 'Runtime')
            for label in labels:
                for tag in soup.find_all(string=re.compile(re.escape(label), re.I)):
                    node = tag
                    for _ in range(3):
                        node = node.find_parent() if node else None
                        if not node:
                            break
                        text_block = node.get_text(" ", strip=True)
                        m = re.search(rf'{re.escape(label)}\s*[:]?\s*([\d,\.]+s?)\b',
                                      text_block, re.I)
                        if m:
                            stats[label] = m.group(1)
                            break
                    if label in stats:
                        break

            per_sample = {}
            mge_table = None
            for table in soup.find_all('table'):
                header_text = ' '.join(th.get_text().strip() for th in table.find_all('th'))
                if 'Sample' in header_text and 'mobileOG' in header_text:
                    mge_table = table
                    break
            if mge_table:
                try:
                    df = pd.read_html(io.StringIO(str(mge_table)))[0]
                    df.columns = [str(c).strip() for c in df.columns]
                    for _, row in df.iterrows():
                        sample = self.normalize_sample_id(row.get('Sample', ''))
                        if not sample:
                            continue
                        entry = {}
                        for col in df.columns:
                            if col == 'Sample':
                                continue
                            v = row.get(col, 0)
                            try:
                                entry[col] = int(str(v).replace(',', '').strip())
                            except Exception:
                                entry[col] = v
                        per_sample[sample] = entry
                except Exception as e:
                    print(f"    Error parsing MGE table: {e}")

            note = ""
            for tag in soup.find_all(string=re.compile(r'Interpretation note', re.I)):
                parent = tag.find_parent()
                if parent:
                    raw = parent.get_text(" ", strip=True)
                    raw = re.sub(r'^Interpretation note\s*:?\s*', '', raw, flags=re.I)
                    if len(raw) > 40:
                        note = raw
                        break

            print(f"    MGE stats: {len(stats)} cards, {len(per_sample)} samples")
            return {'stats': stats, 'per_sample': per_sample, 'note': note}
        except Exception as e:
            print(f"    Error parsing MGE summary: {e}")
            return {'stats': {}, 'per_sample': {}, 'note': ''}


# =============================================================================
# ANALYZER
# =============================================================================
class EcoliAnalyzer:
    """Gene-centric tables, cross-genome patterns, and pathotype classification."""

    def __init__(self):
        self.critical_amr_genes = {
            'blaCTX-M', 'blaCTX-M-1', 'blaCTX-M-2', 'blaCTX-M-14', 'blaCTX-M-15',
            'blaCTX-M-27', 'blaCTX-M-55',
            'blaSHV', 'blaSHV-2', 'blaSHV-5', 'blaSHV-12',
            'blaTEM', 'blaTEM-3', 'blaTEM-52',
            'blaPER', 'blaVEB', 'blaBEL', 'blaGES',
            'blaKPC', 'blaNDM', 'blaOXA', 'blaOXA-48', 'blaOXA-181', 'blaOXA-232',
            'blaVIM', 'blaIMP', 'blaSIM', 'blaSPM',
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5',
            'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            'qnrA', 'qnrB', 'qnrC', 'qnrD', 'qnrS',
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'rmtG', 'rmtH', 'npmA',
            'cfr', 'optrA', 'poxtA', 'tetX',
        }
        self.critical_virulence_genes = {
            'stx1', 'stx1a', 'stx1b', 'stx2', 'stx2a', 'stx2b', 'stx2c', 'stx2d',
            'stx2e', 'stx2f', 'stx2g',
            'cnf1', 'cnf2', 'hlya', 'hlyb', 'hlyc', 'hlyd', 'hlye', 'hlyf',
            'elta', 'eltb', 'estia', 'estib', 'asta', 'east1',
            'clba', 'clbb', 'clbc', 'cdta', 'cdtb', 'cdtc',
            'eae', 'eaea', 'eaeh', 'tir', 'paa',
            'bfpa', 'aggr', 'aata', 'aaic',
            'ipah', 'ipaa', 'ipab', 'ipac',
            'fyua', 'irp1', 'irp2', 'iuca', 'iuta',
        }

    def classify_pathotype(self, virulence_genes: List[str]) -> str:
        genes = set(g.lower() for g in virulence_genes)
        stx = any(g in genes for g in ['stx1', 'stx2', 'stx1a', 'stx1b', 'stx2a', 'stx2b'])
        eae = 'eae' in genes or 'eaea' in genes
        bfp = 'bfpa' in genes
        agg = any(g in genes for g in ['aggr', 'aaic', 'aata'])
        lt = 'ltca' in genes or 'elta' in genes
        st = 'sta1' in genes or 'estia' in genes
        afa = any(g in genes for g in ['afac', 'afae'])
        ipa = 'ipah' in genes
        pathotypes = []
        if stx and not eae: pathotypes.append('STEC')
        if (eae or bfp) and not stx: pathotypes.append('EPEC')
        if stx and eae: pathotypes.append('EHEC')
        if agg: pathotypes.append('EAEC')
        if lt or st: pathotypes.append('ETEC')
        if afa: pathotypes.append('DAEC')
        if ipa: pathotypes.append('EIEC')
        return '-'.join(pathotypes) if pathotypes else 'none'

    def create_gene_centric_tables(self, integrated: Dict, total_samples: int) -> Dict[str, Any]:
        gc = {
            'amr_databases': {},
            'virulence_databases': {},
            'bacmet_databases': {},
            'plasmid_databases': {},
            'database_stats': {},
        }
        if 'amrfinder' in integrated.get('gene_frequencies', {}):
            data = integrated['gene_frequencies']['amrfinder']
            lst = [{'gene': g, 'database': 'AMRfinder',
                    'count': d.get('count', 0), 'percentage': d.get('percentage', 0),
                    'frequency_display': d.get('frequency_display', ''),
                    'genomes': d.get('genomes', [])} for g, d in data.items()]
            if lst:
                lst.sort(key=lambda x: x['count'], reverse=True)
                gc['amr_databases']['amrfinder'] = lst
        if 'abricate' in integrated.get('gene_frequencies', {}):
            for db_name, db_genes in integrated['gene_frequencies']['abricate'].items():
                lst = [{'gene': g, 'database': db_name.upper(),
                        'count': d.get('count', 0), 'percentage': d.get('percentage', 0),
                        'frequency_display': d.get('frequency_display', ''),
                        'genomes': d.get('genomes', [])} for g, d in db_genes.items()]
                if not lst:
                    continue
                lst.sort(key=lambda x: x['count'], reverse=True)
                if db_name in ('vfdb', 'ecoli_vf'):
                    gc['virulence_databases'][db_name] = lst
                elif db_name == 'bacmet2':
                    gc['bacmet_databases'][db_name] = lst
                elif db_name in ('plasmidfinder', 'ecoh'):
                    gc['plasmid_databases'][db_name] = lst
                else:
                    gc['amr_databases'][db_name] = lst
        for cat in ('amr_databases', 'virulence_databases', 'bacmet_databases', 'plasmid_databases'):
            for db_name, genes in gc.get(cat, {}).items():
                gc['database_stats'][db_name] = {
                    'total_genes': len(genes),
                    'total_occurrences': sum(g['count'] for g in genes),
                }
        return gc

    def create_cross_genome_patterns(self, integrated: Dict) -> Dict[str, Any]:
        patterns = {
            'st_distribution': Counter(),
            'serotype_distribution': Counter(),
            'phylogroup_distribution': Counter(),
            'ch_type_distribution': Counter(),
            'o_type_distribution': Counter(),
            'h_type_distribution': Counter(),
            'fumc_type_distribution': Counter(),
            'fimh_type_distribution': Counter(),
            'pathotype_distribution': Counter(),
            'st_o_combinations': defaultdict(list),
            'st_h_combinations': defaultdict(list),
            'st_oh_combinations': defaultdict(list),
            'st_phylogroup_combinations': defaultdict(list),
            'st_fumc_fimh_combinations': defaultdict(list),
            'st_oh_phylo_fumc_fimh_combinations': defaultdict(list),
            'gene_cooccurrence': defaultdict(Counter),
            'high_risk_combinations': [],
        }
        samples = integrated.get('samples', {})
        gc = integrated.get('gene_centric', {})
        sample_genes = defaultdict(list)
        for cat in ('amr_databases', 'virulence_databases', 'bacmet_databases', 'plasmid_databases'):
            for genes in gc.get(cat, {}).values():
                for g in genes:
                    for genome in g['genomes']:
                        if g['gene'] not in sample_genes[genome]:
                            sample_genes[genome].append(g['gene'])
        for sample, data in samples.items():
            st = data.get('mlst', {}).get('ST', 'ND')
            sero = data.get('serotype', {}).get('Serotype', 'ND')
            o_type = data.get('serotype', {}).get('O_Type', 'ND')
            h_type = data.get('serotype', {}).get('H_Type', 'ND')
            phylo = data.get('phylogrouping', {}).get('Clermont_Type', 'ND')
            fumc = data.get('chtyper', {}).get('FumC_Type', 'ND')
            fimh = data.get('chtyper', {}).get('FimH_Type', 'ND')
            cht = data.get('chtyper', {}).get('CH_Type', 'ND')
            patho = data.get('pathotype', 'none')
            if st != 'ND': patterns['st_distribution'][st] += 1
            if sero != 'ND': patterns['serotype_distribution'][sero] += 1
            if o_type != 'ND': patterns['o_type_distribution'][o_type] += 1
            if h_type != 'ND': patterns['h_type_distribution'][h_type] += 1
            if phylo != 'ND': patterns['phylogroup_distribution'][phylo] += 1
            if cht != 'ND': patterns['ch_type_distribution'][cht] += 1
            if fumc != 'ND': patterns['fumc_type_distribution'][fumc] += 1
            if fimh != 'ND': patterns['fimh_type_distribution'][fimh] += 1
            patterns['pathotype_distribution'][patho] += 1
            if st != 'ND' and o_type != 'ND':
                patterns['st_o_combinations'][f"{st} - {o_type}"].append(sample)
            if st != 'ND' and h_type != 'ND':
                patterns['st_h_combinations'][f"{st} - {h_type}"].append(sample)
            if st != 'ND' and sero != 'ND':
                patterns['st_oh_combinations'][f"{st} - {sero}"].append(sample)
            if st != 'ND' and phylo != 'ND':
                patterns['st_phylogroup_combinations'][f"{st} - {phylo}"].append(sample)
            if st != 'ND' and cht != 'ND':
                patterns['st_fumc_fimh_combinations'][f"{st} - {cht}"].append(sample)
            if st != 'ND' and sero != 'ND' and phylo != 'ND' and cht != 'ND':
                full = f"{st} - {sero} - {phylo} - {cht}"
                patterns['st_oh_phylo_fumc_fimh_combinations'][full].append(sample)
            genes = sample_genes.get(sample, [])
            for i, g1 in enumerate(genes):
                for g2 in genes[i+1:]:
                    patterns['gene_cooccurrence'][g1][g2] += 1
            amr_genes = data.get('amr_genes', [])
            vir_genes = data.get('virulence_genes', [])
            crit_amr = [g for g in amr_genes if any(c.lower() in g.lower() for c in self.critical_amr_genes)]
            crit_vir = [g for g in vir_genes if any(c.lower() in g.lower() for c in self.critical_virulence_genes)]
            if crit_amr and crit_vir:
                patterns['high_risk_combinations'].append({
                    'sample': sample, 'st': st, 'serotype': sero,
                    'phylogroup': phylo,
                    'critical_amr_genes': crit_amr,
                    'critical_virulence_genes': crit_vir,
                })
        return patterns

    def compute_mge_stats(self, mge_data: Dict[str, Any]) -> Dict[str, Any]:
        per_sample = mge_data.get('per_sample', {})
        if not per_sample:
            return {'category_totals': [], 'grand_total': 0}
        all_cols = set()
        for s in per_sample.values():
            all_cols.update(s.keys())
        numeric_cols = [c for c in all_cols
                        if all(isinstance(s.get(c, 0), (int, float)) for s in per_sample.values())]
        category_totals = []
        n = len(per_sample)
        grand = 0
        for col in numeric_cols:
            vals = {s: per_sample[s].get(col, 0) for s in per_sample}
            total = sum(vals.values())
            grand += total
            max_sample = max(vals, key=vals.get) if vals else 'ND'
            category_totals.append({
                'category': col, 'total': total,
                'mean': round(total / n, 2) if n else 0,
                'max': vals[max_sample] if vals else 0,
                'max_sample': max_sample,
            })
        category_totals.sort(key=lambda x: x['total'], reverse=True)
        for c in category_totals:
            c['pct_of_all'] = round((c['total'] / grand) * 100, 2) if grand else 0
        return {'category_totals': category_totals, 'grand_total': grand, 'n_samples': n}


# =============================================================================
# HTML GENERATOR
# =============================================================================
class EcoliHTMLGenerator:
    """Builds the interactive HTML report."""

    ACTIVE_COLOR = '#2e7d32'

    def __init__(self, analyzer: EcoliAnalyzer):
        self.data_analyzer = analyzer
        self._current_samples_data = {}

    # ------------------------------------------------------------------
    # Reusable fragments
    # ------------------------------------------------------------------
    def _credit_bar(self, color: str, icon: str, title: str, body: str) -> str:
        return (f'<div style="background:linear-gradient(135deg,#f8f9fa 0%,#f0f4f8 100%);'
                f'border-left:6px solid {color};margin-bottom:20px;padding:15px;border-radius:8px;">'
                f'<div style="display:flex;align-items:center;gap:12px;flex-wrap:wrap;">'
                f'<span style="font-size:1.4em;">{icon}</span>'
                f'<div><strong style="font-size:1.1em;color:{color};">{title}</strong><br>'
                f'<span style="font-size:.95em;color:#333;">{body}</span></div></div></div>')

    def _alert(self, kind: str, icon: str, body: str) -> str:
        return f'<div class="alert-box alert-{kind}"><i class="fas {icon} fa-2x"></i><div>{body}</div></div>'

    def _stat_card(self, value, label: str, color: str = '#4CAF50', icon: str = '') -> str:
        ic = f'<i class="fas {icon} fa-2x" style="opacity:.9;margin-bottom:6px;"></i>' if icon else ''
        return (f'<div class="stat-card" style="background:linear-gradient(135deg,{color} 0%,{color}dd 100%);">'
                f'{ic}<div class="stat-value">{value}</div><div class="stat-label">{label}</div></div>')

    def _collapsible(self, title: str, color: str, inner_html: str, open_by_default: bool = False) -> str:
        open_attr = ' open' if open_by_default else ''
        return (f'<details class="edu-details"{open_attr}>'
                f'<summary style="cursor:pointer;padding:12px 16px;background:#f0f7f0;'
                f'border-left:4px solid {color};border-radius:6px;margin:12px 0;'
                f'font-weight:600;color:#2c3e50;">{title}</summary>'
                f'<div style="padding:6px 0 16px 0;">{inner_html}</div></details>')

    def _grouping_dropdown(self, table_id: str) -> str:
        return f'''
        <div class="grouping-controls">
            <strong><i class="fas fa-layer-group"></i> Group genomes by:</strong>
            <select class="group-select" onchange="groupGenomesByTyping('{table_id}', this.value)">
                <option value="">— None (flat list) —</option>
                <optgroup label="Single field">
                    <option value="ST">MLST (ST)</option>
                    <option value="Serotype">Serotype</option>
                    <option value="Phylogroup">Phylogroup</option>
                    <option value="CHtype">CH type</option>
                </optgroup>
                <optgroup label="Two fields">
                    <option value="ST+Serotype">ST + Serotype</option>
                    <option value="ST+Phylogroup">ST + Phylogroup</option>
                    <option value="ST+CHtype">ST + CH type</option>
                    <option value="Serotype+Phylogroup">Serotype + Phylogroup</option>
                    <option value="Serotype+CHtype">Serotype + CH type</option>
                    <option value="Phylogroup+CHtype">Phylogroup + CH type</option>
                </optgroup>
                <optgroup label="Three fields">
                    <option value="ST+Serotype+Phylogroup">ST + Serotype + Phylogroup</option>
                    <option value="ST+Serotype+CHtype">ST + Serotype + CH type</option>
                    <option value="ST+Phylogroup+CHtype">ST + Phylogroup + CH type</option>
                </optgroup>
                <optgroup label="All four">
                    <option value="Full">Full (all four)</option>
                </optgroup>
            </select>
            <button class="group-btn" onclick="resetGenomeList('{table_id}')"><i class="fas fa-undo"></i> Reset</button>
        </div>'''

    def _search_pair(self, table_id: str, placeholder: str) -> str:
        return (f'<input type="text" class="search-box" id="search-{table_id}" '
                f'onkeyup="searchTable(\'{table_id}\',\'search-{table_id}\')" '
                f'placeholder="🔍 {placeholder}">'
                f'<input type="text" class="search-box" id="highlight-{table_id}" '
                f'onkeyup="highlightGenome(\'{table_id}\',\'highlight-{table_id}\')" '
                f'placeholder="🔍 Highlight genome tags containing text">')

    def _combo_table(self, table_id: str, title: str, combo_dict: dict) -> str:
        if not combo_dict:
            return ''
        rows = ''
        for combo, samples in sorted(combo_dict.items(), key=lambda x: len(x[1]), reverse=True):
            tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            rows += (f'<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td>'
                     f'<td><div class="genome-list">{tags}</div></td></tr>')
        return f'''
        <h3>🔗 {title}</h3>
        {self._search_pair(table_id, f'Search {title}...')}
        <div class="master-scrollable-container">
            <table id="{table_id}" class="data-table">
                <thead><tr>
                    <th data-sort="string">{title}</th>
                    <th data-sort="number">Count</th>
                    <th data-sort="string">Samples</th>
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    def _gene_table(self, table_id: str, title: str, gene_list: list,
                    total_samples: int, critical_set: set = None) -> str:
        critical_set = critical_set or set()
        rows = ''
        for g in gene_list:
            is_crit = any(c.lower() in g['gene'].lower() for c in critical_set)
            gene_disp = f"<strong>{g['gene']}</strong>" + (" ⚠️" if is_crit else "")
            pct = f"{g.get('percentage', 0):.1f}%"
            tags = ''.join(f'<span class="genome-tag">{x}</span>' for x in g.get('genomes', []))
            rows += (f'<tr><td>{gene_disp}</td><td>{g["database"]}</td>'
                     f'<td><strong>{g["count"]}</strong></td><td>{pct}</td>'
                     f'<td><div class="genome-list">{tags}</div></td></tr>')
        return f'''
        <h3><i class="fas fa-list"></i> {title}</h3>
        {self._grouping_dropdown(table_id)}
        {self._search_pair(table_id, f'Search {title}...')}
        <div class="master-scrollable-container">
            <table id="{table_id}" class="data-table">
                <thead><tr>
                    <th data-sort="string">Gene</th>
                    <th data-sort="string">Database</th>
                    <th data-sort="number">Count</th>
                    <th data-sort="number">Percentage</th>
                    <th data-sort="string">Genomes (scrollable, groupable)</th>
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    def _database_cards(self, gene_dict: dict) -> str:
        html = ('<h3 style="margin-top:30px;"><i class="fas fa-database"></i> Database Summary</h3>'
                '<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));'
                'gap:20px;margin:20px 0;">')
        for db, genes in gene_dict.items():
            label = db.upper() if db != 'amrfinder' else 'AMRfinder'
            top = ', '.join(f"{g['gene']} ({g['count']})" for g in genes[:3])
            total = sum(g['count'] for g in genes)
            html += (f'<div class="database-section"><h4>{label}</h4>'
                     f'<p><strong>{len(genes)} unique genes</strong> (total occurrences: {total})</p>'
                     f'<p>Top genes: {top}</p></div>')
        html += '</div>'
        return html

    # ------------------------------------------------------------------
    # CSS
    # ------------------------------------------------------------------
    def _get_css(self) -> str:
        return f"""
<style>
:root{{
  --active-color:{self.ACTIVE_COLOR};
}}
*{{margin:0;padding:0;box-sizing:border-box}}
body{{font-family:'Segoe UI',Tahoma,Geneva,Verdana,sans-serif;line-height:1.6;color:#333;
     background:#f5f5f5;min-width:1200px}}
.container{{max-width:none;margin:0 auto;padding:20px;width:100%;overflow-x:auto}}
.main-header{{background:linear-gradient(135deg,#1e3c72 0%,#2a5298 100%);color:white;
             padding:30px;border-radius:15px;margin-bottom:30px;text-align:center}}
.main-header h1{{font-size:2.6em;margin-bottom:10px;color:white}}
.metadata-bar{{background:rgba(255,255,255,.1);padding:15px;border-radius:10px;margin:20px 0;
              display:flex;justify-content:space-around;flex-wrap:wrap;gap:15px}}
.metadata-item{{display:flex;align-items:center;gap:8px;font-size:.95em}}
.dashboard-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));
                gap:20px;margin-bottom:30px}}
.dashboard-card{{background:white;padding:20px;border-radius:12px;
                box-shadow:0 5px 20px rgba(0,0,0,.1);text-align:center;
                transition:all .3s ease;cursor:pointer;border-left:5px solid}}
.dashboard-card:hover{{transform:translateY(-8px);box-shadow:0 15px 30px rgba(0,0,0,.18)}}
.card-number{{font-size:2.5em;font-weight:bold;margin:10px 0;
             background:linear-gradient(90deg,#1e3c72,#2a5298);
             -webkit-background-clip:text;-webkit-text-fill-color:transparent}}
.card-label{{font-size:.9em;color:#555;font-weight:600}}
.tab-navigation{{display:flex;gap:5px;margin-bottom:20px;flex-wrap:wrap;
                background:white;padding:15px;border-radius:12px;position:sticky;
                top:10px;z-index:100;box-shadow:0 5px 20px rgba(0,0,0,.08)}}
.tab-button{{padding:10px 16px;background:#f5f5f5;border:none;border-radius:8px;
            cursor:pointer;font-weight:600;color:#666;transition:all .3s ease;
            display:flex;align-items:center;gap:6px;font-size:.88em}}
.tab-button:hover{{background:#e8f5e9;color:var(--active-color)}}
.tab-button.active{{background:var(--active-color);color:white}}
.tab-content{{display:none;background:white;padding:30px;border-radius:15px;
             margin-bottom:30px;width:100%;overflow-x:auto}}
.tab-content.active{{display:block}}
.section-header{{color:#2c3e50;margin-bottom:25px;padding-bottom:15px;
                border-bottom:3px solid var(--active-color);font-size:1.8em;
                display:flex;justify-content:space-between;align-items:center}}
.data-table{{width:100%;border-collapse:separate;border-spacing:0;margin:20px 0;
            font-size:.95em;box-shadow:0 2px 10px rgba(0,0,0,.1)}}
.data-table th{{background:#2c3e50;color:white;padding:12px;text-align:left;
               cursor:pointer;white-space:nowrap;position:sticky;top:0;z-index:6}}
.data-table th:hover{{background:#1a252f}}
.data-table td{{padding:10px;border-bottom:1px solid #e0e0e0;vertical-align:top;
               word-wrap:break-word}}
.data-table tr:hover{{background:#f8f9fa}}
.scrollable-table,.master-scrollable-container{{max-height:640px;overflow:auto;
                  border:1px solid #e0e0e0;border-radius:8px;margin:20px 0;width:100%}}
.genome-list{{display:flex;flex-wrap:wrap;gap:5px;max-height:260px;overflow:auto;
             padding:6px;background:#f8f9fa;border-radius:5px}}
.genome-tag{{background:#e0f2f1;color:#00695c;padding:3px 10px;border-radius:12px;
            font-size:.85em;border:1px solid #b2dfdb;white-space:nowrap;margin:2px}}
.genome-tag.highlight{{background:#ffff99 !important;color:#000 !important;
                      border:1px solid #ffc107}}
.genome-group{{margin-bottom:10px;width:100%}}
.genome-group-header{{font-weight:bold;background:#e0e0e0;padding:4px 8px;
                     border-radius:4px;margin:5px 0;font-size:.85em;
                     display:inline-block;word-break:break-word}}
.genome-group-tags{{display:flex;flex-wrap:wrap;gap:5px;margin-left:10px}}
.grouping-controls{{background:#f0f7f0;padding:12px;border-radius:8px;margin:15px 0;
                   display:flex;flex-wrap:wrap;gap:10px;align-items:center;
                   border-left:4px solid var(--active-color)}}
.group-select{{padding:8px 14px;border:2px solid var(--active-color);border-radius:6px;
              background:white;color:#333;font-size:.9em;font-weight:600;
              cursor:pointer;min-width:240px}}
.group-btn{{background:white;border:1px solid var(--active-color);color:var(--active-color);
           padding:6px 12px;border-radius:20px;cursor:pointer;font-size:.85em;
           transition:all .2s}}
.group-btn:hover{{background:var(--active-color);color:white}}
.search-box{{width:100%;padding:12px;margin-bottom:10px;border:2px solid #e0e0e0;
            border-radius:8px;font-size:1em}}
.search-box:focus{{outline:none;border-color:var(--active-color);
                  box-shadow:0 0 0 3px rgba(46,125,50,.1)}}
.alert-box{{padding:20px;border-radius:10px;margin:20px 0;display:flex;
           align-items:flex-start;gap:20px;border-left:5px solid}}
.alert-info{{background:#d1ecf1;color:#0c5460;border-left-color:#17a2b8}}
.alert-success{{background:#d4edda;color:#155724;border-left-color:#28a745}}
.alert-warning{{background:#fff3cd;color:#856404;border-left-color:#ffc107}}
.alert-danger{{background:#f8d7da;color:#721c24;border-left-color:#dc3545}}
.action-buttons{{display:flex;gap:10px;margin:20px 0;flex-wrap:wrap;align-items:center}}
.action-btn{{padding:9px 18px;border:none;border-radius:8px;cursor:pointer;
            font-weight:600;display:inline-flex;align-items:center;gap:8px;
            transition:all .3s ease;text-decoration:none;font-size:.9em}}
.action-btn:hover{{transform:translateY(-2px);box-shadow:0 5px 15px rgba(0,0,0,.15)}}
.btn-primary{{background:var(--active-color);color:white}}
.btn-success{{background:#28a745;color:white}}
.btn-danger{{background:#dc3545;color:white}}
.btn-warning{{background:#ffc107;color:black}}
.btn-info{{background:#17a2b8;color:white}}
.btn-secondary{{background:#6c757d;color:white}}
.btn-light{{background:#f8f9fa;color:#212529;border:1px solid #dee2e6}}
.database-section{{margin:30px 0;padding:22px;border-radius:12px;background:#f8f9fa;
                  box-shadow:0 3px 15px rgba(0,0,0,.08)}}
.print-section-btn{{background:var(--active-color);color:white;border:none;border-radius:5px;
                   padding:8px 15px;cursor:pointer;display:flex;align-items:center;
                   gap:5px;font-size:.9em}}
.footer{{text-align:center;padding:30px;color:white;margin-top:40px;
        border-radius:15px;background:linear-gradient(135deg,#2c3e50 0%,#34495e 100%)}}
.stat-card{{color:white;padding:18px;border-radius:10px;text-align:center;
           box-shadow:0 4px 15px rgba(0,0,0,.15);transition:transform .2s}}
.stat-card:hover{{transform:translateY(-3px)}}
.stat-card .stat-value{{font-size:1.9em;font-weight:bold;margin-bottom:4px}}
.stat-card .stat-label{{font-size:.85em;opacity:.95;text-transform:uppercase;letter-spacing:.5px}}
.stats-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));
            gap:15px;margin:20px 0}}
.badge{{display:inline-block;padding:5px 15px;border-radius:20px;font-size:.85em;
       font-weight:600;margin:2px}}
.info-text{{background:#f8f9fa;padding:15px;border-radius:8px;margin:15px 0;
           border-left:4px solid #17a2b8}}
.col-toggle-wrap{{position:relative}}
.col-panel{{position:absolute;right:0;top:110%;background:#fff;border:1px solid #dee2e6;
           border-radius:10px;padding:10px 14px;box-shadow:0 10px 30px rgba(0,0,0,.15);
           z-index:30;display:none;max-height:320px;overflow-y:auto;min-width:200px}}
.col-panel.open{{display:block}}
.col-panel label{{display:block;font-size:13px;color:#1e3c72;padding:4px 0;cursor:pointer}}
.sel-counter{{font-size:12px;font-weight:700;color:#1e3c72;background:#e3f2fd;
             padding:6px 12px;border-radius:20px}}
.col-check{{width:34px;text-align:center}}
.row-selected td{{background:#e3f2fd !important}}
.mode-toggle{{display:flex;gap:8px;margin-bottom:18px}}
.mode-btn{{padding:10px 22px;border:2px solid #cbd5e1;background:white;color:#475569;
          border-radius:8px;font-weight:600;cursor:pointer;font-size:.95em;
          display:flex;align-items:center;gap:8px;transition:all .2s}}
.mode-btn:hover{{border-color:var(--active-color);color:var(--active-color)}}
.mode-btn.active{{background:var(--active-color);border-color:var(--active-color);color:white}}
.compare-mode{{display:none}}
.compare-mode.active{{display:block}}
.compare-panel{{background:#f8f9fa;border-radius:12px;padding:20px;margin-bottom:24px}}
.compare-select-row{{display:flex;gap:12px;align-items:flex-end;flex-wrap:wrap;margin-bottom:16px}}
.compare-select-row label{{display:flex;flex-direction:column;gap:4px}}
.compare-select-row label span{{font-size:.85em;font-weight:600;color:#475569}}
.compare-select-row select{{padding:10px 14px;border:2px solid #cbd5e1;border-radius:8px;
                           background:white;font-size:.95em;min-width:220px}}
.compare-visuals-grid{{display:grid;grid-template-columns:180px 1fr;gap:16px;
                      margin:20px 0;align-items:center}}
.visual-card{{background:white;border-radius:12px;padding:20px;
             box-shadow:0 2px 8px rgba(0,0,0,.06);border:1px solid #e2e8f0}}
.gauge-card{{display:flex;justify-content:center;align-items:center;padding:12px}}
.similarity-gauge{{width:160px;height:160px;display:block}}
.metrics-card{{display:flex;flex-direction:column;gap:18px}}
.metric-block{{display:flex;flex-direction:column;gap:6px}}
.metric-label{{font-size:.85em;font-weight:700;color:#475569;text-transform:uppercase;letter-spacing:.5px}}
.metric-bar-track{{background:#f1f5f9;height:22px;border-radius:11px;overflow:hidden}}
.metric-bar-fill{{height:100%;border-radius:11px;transition:width .8s ease;min-width:4px}}
.metric-value{{font-size:.88em;color:#334155;font-weight:600}}
.typing-pill{{display:inline-block;padding:5px 12px;border-radius:12px;
             font-size:.85em;font-weight:600;border:1.5px solid;
             font-family:'Consolas',monospace}}
.gene-proportion-bar{{display:flex;height:12px;border-radius:6px;overflow:hidden;
                     margin-bottom:14px;background:#f1f5f9}}
.gene-prop-segment{{height:100%;transition:width .6s ease;min-width:2px}}
.gene-prop-segment.shared{{background:#10b981}}
.gene-prop-segment.only-a{{background:#3b82f6}}
.gene-prop-segment.only-b{{background:#ec4899}}
.cluster-card{{background:white;border-radius:10px;padding:16px 20px;
              margin-bottom:12px;border-left:5px solid;
              box-shadow:0 2px 8px rgba(0,0,0,.06)}}
.cluster-header{{display:flex;align-items:center;gap:14px;margin-bottom:12px;flex-wrap:wrap}}
.cluster-badge{{color:white;padding:4px 14px;border-radius:14px;font-size:.85em;font-weight:700}}
.cluster-stats{{font-size:.9em;color:#64748b;font-weight:600}}
.cluster-samples{{display:flex;flex-wrap:wrap;gap:6px}}
.cluster-sample-pill{{background:#f1f5f9;color:#1e293b;padding:4px 12px;border-radius:12px;
                     font-size:.85em;font-weight:600;font-family:'Consolas',monospace;
                     border:1px solid #cbd5e1}}
#cluster-clusters{{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));
                  gap:12px}}
#cluster-clusters:empty{{display:none}}
#cluster-clusters .cluster-card{{margin-bottom:0}}
.heatmap-wrapper{{background:white;border-radius:12px;padding:20px;margin-top:24px;
                 box-shadow:0 2px 8px rgba(0,0,0,.06);overflow-x:auto;width:100%}}
.heatmap-grid{{display:grid;gap:3px;margin-top:16px;padding:4px;width:100%;
              min-width:max-content}}
.heatmap-cell{{width:100%;height:100%;min-width:34px;min-height:34px;
              display:flex;align-items:center;justify-content:center;
              font-size:.75em;font-weight:700;border-radius:4px}}
.heatmap-col-header{{background:transparent;color:#64748b;font-size:.72em;
                    writing-mode:vertical-rl;text-orientation:mixed;
                    height:150px;min-width:44px;font-family:'Consolas',monospace;
                    display:flex;align-items:flex-end;justify-content:center;
                    padding-bottom:6px;white-space:nowrap;overflow:visible}}
.heatmap-row-header{{background:transparent;color:#64748b;font-size:.78em;
                    display:flex;align-items:center;justify-content:flex-start;
                    padding-left:6px;font-family:'Consolas',monospace;
                    white-space:nowrap;min-width:200px}}
.heatmap-corner{{background:transparent}}
.heatmap-cell[data-tip]{{position:relative;cursor:help}}
.heatmap-cell[data-tip]:hover{{filter:brightness(1.18);z-index:60}}
.heatmap-cell[data-tip]::after{{
  content:attr(data-tip);
  position:absolute;bottom:calc(100% + 8px);left:50%;transform:translateX(-50%);
  background:#1e293b;color:#fff;padding:6px 12px;border-radius:6px;
  font-size:12px;font-weight:500;white-space:nowrap;pointer-events:none;
  opacity:0;transition:opacity .15s ease;z-index:300;
  font-family:'Segoe UI',Tahoma,sans-serif;letter-spacing:.2px;
  box-shadow:0 4px 14px rgba(0,0,0,.28);
  writing-mode:horizontal-tb;text-orientation:mixed;
}}
.heatmap-cell[data-tip]::before{{
  content:'';position:absolute;bottom:100%;left:50%;transform:translateX(-50%);
  border:6px solid transparent;border-top-color:#1e293b;pointer-events:none;
  opacity:0;transition:opacity .15s ease;z-index:300;
}}
.heatmap-cell[data-tip]:hover::after,
.heatmap-cell[data-tip]:hover::before{{opacity:1}}
.heatmap-col-header[data-tip]::after{{
  bottom:auto;top:50%;left:calc(100% + 10px);transform:translateY(-50%);
}}
.heatmap-col-header[data-tip]::before{{
  bottom:auto;top:50%;left:100%;transform:translateY(-50%) rotate(-90deg);
  border-top-color:transparent;border-right-color:#1e293b;
}}
.heatmap-row-header[data-tip]::after{{
  bottom:auto;top:50%;left:calc(100% + 10px);transform:translateY(-50%);
}}
.heatmap-row-header[data-tip]::before{{
  bottom:auto;top:50%;left:100%;transform:translateY(-50%);
  border-top-color:transparent;border-right-color:#1e293b;
}}
.compare-section{{background:white;border-radius:12px;padding:20px;margin-bottom:20px;
                 box-shadow:0 2px 8px rgba(0,0,0,.06)}}
.compare-section h3{{margin-bottom:14px;color:#1e293b;font-size:1.15em;
                    padding-bottom:8px;border-bottom:2px solid #e2e8f0}}
.compare-table{{display:flex;flex-direction:column;gap:2px}}
.compare-row{{display:grid;grid-template-columns:180px 1fr 1fr;gap:12px;padding:10px 14px;
             border-radius:6px;align-items:center;font-size:.92em}}
.compare-row.header{{background:#f1f5f9;font-weight:700;border-bottom:2px solid #cbd5e1}}
.compare-row.diff{{background:#fef2f2;border-left:3px solid #dc2626}}
.compare-row.match{{background:#f0fdf4;border-left:3px solid #16a34a}}
.compare-row .row-label{{font-weight:600;color:#475569;font-size:.88em;
                        text-transform:uppercase;letter-spacing:.3px}}
.compare-row .row-val{{color:#1e293b;word-break:break-word}}
.compare-gene-block{{margin-bottom:24px;padding-bottom:20px;border-bottom:1px dashed #e2e8f0}}
.compare-gene-block:last-child{{border-bottom:none;margin-bottom:0;padding-bottom:0}}
.gene-count{{font-size:.82em;font-weight:400;color:#64748b;background:#f1f5f9;
            padding:3px 10px;border-radius:12px}}
.gene-columns{{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:14px}}
.gene-col{{background:#f8fafc;border-radius:8px;padding:12px;border:1px solid #e2e8f0}}
.gene-col-head{{font-size:.85em;font-weight:700;color:#475569;margin-bottom:8px;
               text-transform:uppercase;letter-spacing:.3px}}
.compare-gene{{display:inline-block;font-size:.82em;padding:3px 10px;border-radius:12px;
              margin:2px;font-family:'Consolas',monospace}}
.compare-gene.shared{{background:#dcfce7;color:#14532d;border:1px solid #86efac}}
.compare-gene.only-a{{background:#dbeafe;color:#1e3a8a;border:1px solid #93c5fd}}
.compare-gene.only-b{{background:#fce7f3;color:#831843;border:1px solid #f9a8d4}}
.compare-full-details{{background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
                      padding:12px 16px;margin-bottom:10px}}
.compare-full-details[open]{{background:#fff;border-color:var(--active-color)}}
.compare-full-details summary{{cursor:pointer;font-weight:700;color:var(--active-color);
                              list-style:none;padding:4px 0}}
.compare-full-profile{{display:grid;grid-template-columns:1fr 1fr;gap:16px;
                      margin-top:14px;padding-top:14px;border-top:1px dashed #e2e8f0}}
.profile-col{{background:#f8fafc;border-radius:10px;padding:14px;border:1px solid #e2e8f0}}
.profile-head{{display:flex;justify-content:space-between;align-items:center;gap:8px;
              margin-bottom:10px;padding-bottom:8px;border-bottom:2px solid #e2e8f0}}
.profile-sample-label{{font-weight:700;color:#1e293b;font-size:.92em;word-break:break-word}}
.profile-sample-count{{font-size:.78em;font-weight:600;color:var(--active-color);background:#dcfce7;
                      padding:2px 10px;border-radius:12px;white-space:nowrap}}
.empty-note{{font-size:.85em;color:#94a3b8;font-style:italic;padding:6px 0;display:inline-block}}
.compare-warning{{padding:20px;background:#fef2f2;color:#7f1d1d;border-radius:8px;
                 border-left:4px solid #dc2626}}
.compare-verdict-box{{padding:16px 20px;border-radius:10px;margin-bottom:16px;
                     display:flex;gap:12px;align-items:flex-start;font-size:.95em}}
.compare-verdict-box.danger{{background:#fef2f2;border-left:5px solid #dc2626;color:#7f1d1d}}
.compare-verdict-box.warning{{background:#fffbeb;border-left:5px solid #f59e0b;color:#78350f}}
.compare-verdict-box.info{{background:#eff6ff;border-left:5px solid #2563eb;color:#1e3a8a}}
.cluster-controls{{display:flex;gap:14px;align-items:flex-end;flex-wrap:wrap;margin-bottom:20px}}
.cluster-controls select{{padding:10px 14px;border:2px solid #cbd5e1;border-radius:8px;
                         background:white;font-size:.95em;min-width:220px}}
.no-alerts{{padding:20px;background:#eff6ff;color:#1e3a8a;border-radius:8px;
           border-left:4px solid #2563eb}}
.edu-details summary{{list-style:none}}
.edu-details summary::-webkit-details-marker{{display:none}}
.edu-details summary::before{{content:'▸ ';display:inline-block;margin-right:6px;transition:transform .2s}}
.edu-details[open] summary::before{{content:'▾ '}}
@media (max-width:900px){{
  .compare-full-profile{{grid-template-columns:1fr}}
  .compare-row{{grid-template-columns:1fr}}
  .compare-visuals-grid{{grid-template-columns:1fr}}
}}
@media print{{.tab-navigation,.dashboard-grid,.search-box,.action-buttons,.grouping-controls,
             .print-section-btn,.col-panel{{display:none!important}}}}
</style>"""

    # ------------------------------------------------------------------
    # JavaScript
    # ------------------------------------------------------------------
    def _get_js(self, typing_json: str) -> str:
        js = """
<script>
var sampleTyping = __TYPING__;
var originalGenomeLists = {};

function switchTab(name) {
  var contents = document.querySelectorAll('.tab-content');
  var buttons  = document.querySelectorAll('.tab-button');
  for (var i = 0; i < contents.length; i++) contents[i].classList.remove('active');
  for (var j = 0; j < buttons.length; j++)  buttons[j].classList.remove('active');
  var ct = document.getElementById(name + '-tab');
  var bt = document.querySelector('.tab-button.' + name);
  if (ct) ct.classList.add('active');
  if (bt) bt.classList.add('active');
}

function escapeHtml(s) {
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}

function searchTable(tableId, searchId) {
  var input = document.getElementById(searchId); if (!input) return;
  var filter = input.value.toUpperCase();
  var table = document.getElementById(tableId); if (!table) return;
  var rows = table.tBodies[0].rows;
  for (var i = 0; i < rows.length; i++) {
    var found = false;
    var cells = rows[i].getElementsByTagName('td');
    for (var j = 0; j < cells.length; j++) {
      if ((cells[j].textContent || cells[j].innerText).toUpperCase().indexOf(filter) > -1) {
        found = true; break;
      }
    }
    rows[i].style.display = found ? '' : 'none';
  }
}

function highlightGenome(tableId, searchId) {
  var el = document.getElementById(searchId); if (!el) return;
  var filter = el.value.toUpperCase().trim();
  var table = document.getElementById(tableId); if (!table) return;
  var tags = table.querySelectorAll('.genome-tag');
  for (var i = 0; i < tags.length; i++) {
    tags[i].classList.remove('highlight');
    if (filter && tags[i].textContent.toUpperCase().indexOf(filter) > -1) {
      tags[i].classList.add('highlight');
    }
  }
}

function getTypingValue(genome, groupBy) {
  var info = sampleTyping[genome]; if (!info) return "Unknown";
  var st = info.ST || "ND", sero = info.Serotype || "ND", phy = info.Phylogroup || "ND", cht = info.CHtype || "ND";
  switch (groupBy) {
    case 'ST': return st;
    case 'Serotype': return sero;
    case 'Phylogroup': return phy;
    case 'CHtype': return cht;
    case 'ST+Serotype': return st + ' - ' + sero;
    case 'ST+Phylogroup': return st + ' - ' + phy;
    case 'ST+CHtype': return st + ' - ' + cht;
    case 'Serotype+Phylogroup': return sero + ' - ' + phy;
    case 'Serotype+CHtype': return sero + ' - ' + cht;
    case 'Phylogroup+CHtype': return phy + ' - ' + cht;
    case 'ST+Serotype+Phylogroup': return st + ' - ' + sero + ' - ' + phy;
    case 'ST+Serotype+CHtype': return st + ' - ' + sero + ' - ' + cht;
    case 'ST+Phylogroup+CHtype': return st + ' - ' + phy + ' - ' + cht;
    case 'Full': return st + ' - ' + sero + ' - ' + phy + ' - ' + cht;
  }
  return "Unknown";
}

function groupRowGenomes(row, groupBy, originalList) {
  var cell = null;
  for (var i = 0; i < row.cells.length; i++) {
    if (row.cells[i].querySelector('.genome-list')) { cell = row.cells[i]; break; }
  }
  if (!cell) return;
  if (originalList.length === 0) { cell.innerHTML = '<div class="genome-list">None</div>'; return; }
  var groups = {};
  for (var k = 0; k < originalList.length; k++) {
    var key = getTypingValue(originalList[k], groupBy);
    if (!groups[key]) groups[key] = [];
    groups[key].push(originalList[k]);
  }
  var html = '<div class="genome-list">';
  for (var key in groups) {
    var tags = groups[key].map(function(g){ return '<span class="genome-tag">' + escapeHtml(g) + '</span>'; }).join('');
    html += '<div class="genome-group"><div class="genome-group-header">' + escapeHtml(key) +
            '</div><div class="genome-group-tags">' + tags + '</div></div>';
  }
  html += '</div>';
  cell.innerHTML = html;
}

function groupGenomesByTyping(tableId, groupBy) {
  var table = document.getElementById(tableId); if (!table) return;
  var tbody = table.tBodies[0]; if (!tbody) return;
  var rows = tbody.rows;
  for (var i = 0; i < rows.length; i++) {
    var row = rows[i];
    var nameCell = row.cells[0]; if (!nameCell) continue;
    var name = nameCell.textContent.trim().replace(/⚠️/g, '').trim();
    if (!originalGenomeLists[name]) {
      var cell = null;
      for (var c = 0; c < row.cells.length; c++) {
        if (row.cells[c].querySelector('.genome-list')) { cell = row.cells[c]; break; }
      }
      if (cell) {
        var tags = cell.querySelectorAll('.genome-tag');
        originalGenomeLists[name] = Array.prototype.map.call(tags, function(t){ return t.textContent.trim(); });
      } else {
        originalGenomeLists[name] = [];
      }
    }
  }
  if (!groupBy) { resetGenomeList(tableId); return; }
  for (var i2 = 0; i2 < rows.length; i2++) {
    var r = rows[i2];
    var nc = r.cells[0]; if (!nc) continue;
    var n = nc.textContent.trim().replace(/⚠️/g, '').trim();
    groupRowGenomes(r, groupBy, originalGenomeLists[n] || []);
  }
}

function resetGenomeList(tableId) {
  var table = document.getElementById(tableId); if (!table) return;
  var tbody = table.tBodies[0]; if (!tbody) return;
  var rows = tbody.rows;
  for (var i = 0; i < rows.length; i++) {
    var row = rows[i];
    var nameCell = row.cells[0]; if (!nameCell) continue;
    var name = nameCell.textContent.trim().replace(/⚠️/g, '').trim();
    var original = originalGenomeLists[name] || [];
    for (var j = 0; j < row.cells.length; j++) {
      if (row.cells[j].querySelector('.genome-list')) {
        row.cells[j].innerHTML = '<div class="genome-list">' +
          original.map(function(g){ return '<span class="genome-tag">' + escapeHtml(g) + '</span>'; }).join('') +
          '</div>';
        break;
      }
    }
  }
  var wrapper = table.closest('.tab-content');
  if (wrapper) {
    var sel = wrapper.querySelector('.group-select');
    if (sel) sel.value = '';
  }
}

function sortTable(tableId, colIndex, type) {
  var table = document.getElementById(tableId); if (!table) return;
  var tbody = table.tBodies[0];
  var rows = Array.prototype.slice.call(tbody.rows);
  var asc = table.getAttribute('data-sort-dir') !== 'asc';
  rows.sort(function(a, b) {
    var av = a.cells[colIndex].innerText.trim();
    var bv = b.cells[colIndex].innerText.trim();
    if (type === 'number') {
      av = parseFloat(av.replace(/,/g, '')) || 0;
      bv = parseFloat(bv.replace(/,/g, '')) || 0;
      return asc ? av - bv : bv - av;
    }
    return asc ? av.localeCompare(bv) : bv.localeCompare(av);
  });
  for (var i = 0; i < rows.length; i++) tbody.appendChild(rows[i]);
  table.setAttribute('data-sort-dir', asc ? 'asc' : 'desc');
}

function exportTableToCSV(tableId, filename) {
  var table = document.getElementById(tableId); if (!table) return;
  var rows = table.querySelectorAll('tr');
  var csv = [];
  for (var i = 0; i < rows.length; i++) {
    var cols = rows[i].querySelectorAll('td, th');
    var row = [];
    for (var j = 0; j < cols.length; j++) {
      row.push('"' + (cols[j].innerText || '').replace(/"/g, '""') + '"');
    }
    csv.push(row.join(','));
  }
  var blob = new Blob([csv.join('\\n')], {type: 'text/csv'});
  var a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(a.href);
}

function printSection(id) {
  var content = document.getElementById(id); if (!content) return;
  var w = window.open('', '_blank');
  var style = document.querySelector('style');
  w.document.write('<html><head><title>Print</title>');
  if (style) w.document.write('<style>' + style.textContent + '</style>');
  w.document.write('</head><body>' + content.innerHTML + '</body></html>');
  w.document.close();
  w.print();
}

function initColumnToggle() {
  var btn = document.getElementById('colToggleBtn');
  var panel = document.getElementById('colPanel');
  if (!btn || !panel) return;
  btn.addEventListener('click', function(e){ e.stopPropagation(); panel.classList.toggle('open'); });
  document.addEventListener('click', function(){ panel.classList.remove('open'); });
  var table = document.getElementById('samples-table'); if (!table) return;
  var headers = table.querySelectorAll('th');
  Array.prototype.forEach.call(headers, function(th, idx) {
    if (th.classList.contains('col-check')) return;
    var cb = document.createElement('input');
    cb.type = 'checkbox'; cb.checked = true;
    cb.addEventListener('change', function() {
      var vis = cb.checked;
      var allRows = table.querySelectorAll('tr');
      Array.prototype.forEach.call(allRows, function(row) {
        var cell = row.children[idx];
        if (cell) cell.style.display = vis ? '' : 'none';
      });
    });
    var lbl = document.createElement('label');
    lbl.appendChild(cb);
    lbl.appendChild(document.createTextNode(' ' + th.textContent.trim()));
    panel.appendChild(lbl);
  });
}

function initRowSelection() {
  var table = document.getElementById('samples-table'); if (!table) return;
  var rows = Array.prototype.slice.call(table.querySelectorAll('tbody tr'));
  var selCount = document.getElementById('selCount');
  var selectAll = document.getElementById('selectAll');
  function update() {
    var n = rows.filter(function(r){
      var cb = r.querySelector('.row-check');
      return cb && cb.checked && r.style.display !== 'none';
    }).length;
    if (selCount) selCount.textContent = n + ' selected';
  }
  rows.forEach(function(r) {
    var cb = r.querySelector('.row-check');
    if (cb) cb.addEventListener('change', function() {
      r.classList.toggle('row-selected', cb.checked);
      update();
    });
  });
  if (selectAll) selectAll.addEventListener('change', function() {
    rows.forEach(function(r) {
      if (r.style.display !== 'none') {
        var cb = r.querySelector('.row-check');
        if (cb) { cb.checked = selectAll.checked; r.classList.toggle('row-selected', cb.checked); }
      }
    });
    update();
  });
}

function setCompareMode(mode) {
  var btns = document.querySelectorAll('.mode-btn');
  Array.prototype.forEach.call(btns, function(b){ b.classList.toggle('active', b.dataset.mode === mode); });
  var modes = document.querySelectorAll('.compare-mode');
  Array.prototype.forEach.call(modes, function(m){ m.classList.remove('active'); });
  var el = document.getElementById('mode-' + mode);
  if (el) el.classList.add('active');
  if (mode === 'cluster') runCluster();
}

function similarityColor(pct) {
  if (pct >= 95) return '#059669'; if (pct >= 90) return '#10b981';
  if (pct >= 80) return '#84cc16'; if (pct >= 70) return '#eab308';
  if (pct >= 50) return '#f59e0b'; if (pct >= 30) return '#f97316';
  return '#dc2626';
}

function computeSimilarity(A, B) {
  var tkeys = Object.keys(A.typing);
  var matches = 0;
  tkeys.forEach(function(k){ if ((A.typing[k] || '') === (B.typing[k] || '')) matches++; });
  var typingPct = tkeys.length ? (matches / tkeys.length) * 100 : 0;
  var cats = ['AMR', 'Virulence', 'BACMET', 'Plasmids', 'Mutations'];
  var scores = [];
  cats.forEach(function(cat) {
    var sa = new Set(A.genes[cat] || []); var sb = new Set(B.genes[cat] || []);
    var u = new Set(); sa.forEach(function(g){ u.add(g); }); sb.forEach(function(g){ u.add(g); });
    if (u.size === 0) return;
    var shared = 0; u.forEach(function(g){ if (sa.has(g) && sb.has(g)) shared++; });
    scores.push((shared / u.size) * 100);
  });
  var genePct = scores.length ? scores.reduce(function(a, b){ return a + b; }, 0) / scores.length : 0;
  return {typing: Math.round(typingPct), genes: Math.round(genePct),
          overall: Math.round(typingPct * 0.5 + genePct * 0.5),
          typingMatches: matches, typingTotal: tkeys.length};
}

function getPrecomputedSim(i, j) {
  if (i === j) return 100;
  if (i > j) { var t = i; i = j; j = t; }
  var sim = window.ECOLITYPER_SIM;
  if (!sim || !sim.flat) return 0;
  var n = sim.n;
  return sim.flat[i * n - i * (i + 1) / 2 + (j - i - 1)];
}

function drawGauge(pct) {
  var color = similarityColor(pct);
  var r = 60, circ = 2 * Math.PI * r, off = circ * (1 - pct / 100);
  return '<svg viewBox="0 0 160 160" class="similarity-gauge">' +
    '<circle cx="80" cy="80" r="' + r + '" fill="none" stroke="#e2e8f0" stroke-width="14"/>' +
    '<circle cx="80" cy="80" r="' + r + '" fill="none" stroke="' + color + '" stroke-width="14" ' +
    'stroke-dasharray="' + circ + '" stroke-dashoffset="' + off + '" stroke-linecap="round" ' +
    'transform="rotate(-90 80 80)"/>' +
    '<text x="80" y="76" text-anchor="middle" font-size="28" font-weight="700" fill="' + color + '">' + pct + '%</text>' +
    '<text x="80" y="98" text-anchor="middle" font-size="11" fill="#64748b">overall</text></svg>';
}

function runCompare() {
  var data = window.ECOLITYPER_COMPARE_DATA || {};
  var a = document.getElementById('compare-a').value;
  var b = document.getElementById('compare-b').value;
  if (!a || !b) return;
  var output = document.getElementById('compare-output');
  var visuals = document.getElementById('compare-visuals');
  var verdict = document.getElementById('compare-verdict');
  if (a === b) {
    output.innerHTML = '<div class="compare-warning">Pick two different samples.</div>';
    visuals.innerHTML = ''; verdict.style.display = 'none'; return;
  }
  if (!data[a] || !data[b]) {
    output.innerHTML = '<div class="compare-warning">Sample data not found.</div>'; return;
  }
  var A = data[a], B = data[b], sim = computeSimilarity(A, B);
  visuals.innerHTML =
    '<div class="compare-visuals-grid">' +
      '<div class="visual-card gauge-card">' + drawGauge(sim.overall) + '</div>' +
      '<div class="visual-card metrics-card">' +
        '<div class="metric-block"><div class="metric-label">Typing match</div>' +
        '<div class="metric-bar-track"><div class="metric-bar-fill" style="width:' + sim.typing + '%;background:' + similarityColor(sim.typing) + ';"></div></div>' +
        '<div class="metric-value">' + sim.typingMatches + '/' + sim.typingTotal + ' fields (' + sim.typing + '%)</div></div>' +
        '<div class="metric-block"><div class="metric-label">Gene content match</div>' +
        '<div class="metric-bar-track"><div class="metric-bar-fill" style="width:' + sim.genes + '%;background:' + similarityColor(sim.genes) + ';"></div></div>' +
        '<div class="metric-value">' + sim.genes + '% Jaccard similarity</div></div>' +
      '</div>' +
    '</div>';
  var typingRows = '';
  Object.keys(A.typing).forEach(function(k) {
    var va = A.typing[k] || 'ND', vb = B.typing[k] || 'ND';
    var diff = va !== vb;
    var vaPill = diff
      ? '<span class="typing-pill" style="background:#fee2e2;color:#dc2626;border-color:#dc2626;">✗ ' + escapeHtml(va) + '</span>'
      : '<span class="typing-pill" style="background:#dcfce7;color:#10b981;border-color:#10b981;">✓ ' + escapeHtml(va) + '</span>';
    var vbPill = diff
      ? '<span class="typing-pill" style="background:#fee2e2;color:#dc2626;border-color:#dc2626;">✗ ' + escapeHtml(vb) + '</span>'
      : '<span class="typing-pill" style="background:#dcfce7;color:#10b981;border-color:#10b981;">✓ ' + escapeHtml(vb) + '</span>';
    typingRows += '<div class="compare-row ' + (diff ? 'diff' : 'match') + '">' +
      '<span class="row-label">' + escapeHtml(k) + '</span>' +
      '<span class="row-val">' + vaPill + '</span>' +
      '<span class="row-val">' + vbPill + '</span></div>';
  });
  var cats = ['AMR', 'Virulence', 'BACMET', 'Plasmids', 'Mutations'];
  var geneHtml = '';
  cats.forEach(function(cat) {
    var ga = A.genes[cat] || [], gb = B.genes[cat] || [];
    var sa = new Set(ga), sb = new Set(gb);
    var u = new Set(); sa.forEach(function(g){ u.add(g); }); sb.forEach(function(g){ u.add(g); });
    if (u.size === 0) return;
    var shared = [], onlyA = [], onlyB = [];
    u.forEach(function(g) {
      if (sa.has(g) && sb.has(g)) shared.push(g);
      else if (sa.has(g)) onlyA.push(g);
      else onlyB.push(g);
    });
    var t = u.size, sPct = shared.length / t * 100, aPct = onlyA.length / t * 100, bPct = onlyB.length / t * 100;
    geneHtml += '<div class="compare-gene-block"><h4>' + cat + ' <span class="gene-count">' +
      shared.length + ' shared · ' + onlyA.length + ' only A · ' + onlyB.length + ' only B (' + sPct.toFixed(0) + '% match)</span></h4>' +
      '<div class="gene-proportion-bar">' +
        '<div class="gene-prop-segment shared" style="width:' + sPct + '%"></div>' +
        '<div class="gene-prop-segment only-a" style="width:' + aPct + '%"></div>' +
        '<div class="gene-prop-segment only-b" style="width:' + bPct + '%"></div>' +
      '</div><div class="gene-columns">' +
        '<div class="gene-col"><div class="gene-col-head">Shared (' + shared.length + ')</div>' +
          (shared.map(function(g){ return '<span class="compare-gene shared">' + escapeHtml(g) + '</span>'; }).join('') || '<em>None</em>') + '</div>' +
        '<div class="gene-col"><div class="gene-col-head">Only in A (' + onlyA.length + ')</div>' +
          (onlyA.map(function(g){ return '<span class="compare-gene only-a">' + escapeHtml(g) + '</span>'; }).join('') || '<em>None</em>') + '</div>' +
        '<div class="gene-col"><div class="gene-col-head">Only in B (' + onlyB.length + ')</div>' +
          (onlyB.map(function(g){ return '<span class="compare-gene only-b">' + escapeHtml(g) + '</span>'; }).join('') || '<em>None</em>') + '</div>' +
      '</div></div>';
  });
  var fullHtml = '<div class="compare-section"><h3>📋 Full Gene Profiles</h3>' +
    '<p style="font-size:.9em;color:#64748b;margin-bottom:12px;">Click a category to expand both samples\\' complete gene lists.</p>';
  cats.forEach(function(cat) {
    var ga = A.genes[cat] || [], gb = B.genes[cat] || [];
    if (ga.length === 0 && gb.length === 0) return;
    var tA = ga.length ? ga.map(function(g){ return '<span class="genome-tag">' + escapeHtml(g) + '</span>'; }).join('') : '<em class="empty-note">None detected</em>';
    var tB = gb.length ? gb.map(function(g){ return '<span class="genome-tag">' + escapeHtml(g) + '</span>'; }).join('') : '<em class="empty-note">None detected</em>';
    fullHtml += '<details class="compare-full-details"><summary>' + cat +
      ' <span class="gene-count">A: ' + ga.length + ' · B: ' + gb.length + '</span></summary>' +
      '<div class="compare-full-profile">' +
        '<div class="profile-col"><div class="profile-head"><span class="profile-sample-label">' + escapeHtml(a) + '</span>' +
        '<span class="profile-sample-count">' + ga.length + ' gene(s)</span></div><div class="genome-list">' + tA + '</div></div>' +
        '<div class="profile-col"><div class="profile-head"><span class="profile-sample-label">' + escapeHtml(b) + '</span>' +
        '<span class="profile-sample-count">' + gb.length + ' gene(s)</span></div><div class="genome-list">' + tB + '</div></div>' +
      '</div></details>';
  });
  fullHtml += '</div>';
  var vHtml;
  if (sim.overall >= 95) {
    vHtml = '<div class="compare-verdict-box danger"><i class="fas fa-exclamation-triangle"></i> <strong>Near-identical (' + sim.overall + '%)</strong> — likely the same strain.</div>';
  } else if (sim.overall >= 85) {
    vHtml = '<div class="compare-verdict-box warning"><i class="fas fa-info-circle"></i> <strong>High similarity (' + sim.overall + '%)</strong> — very likely same clonal lineage.</div>';
  } else if (sim.overall >= 60) {
    vHtml = '<div class="compare-verdict-box warning"><i class="fas fa-info-circle"></i> <strong>Moderate similarity (' + sim.overall + '%)</strong> — related but not recent transmission.</div>';
  } else {
    vHtml = '<div class="compare-verdict-box info"><i class="fas fa-info-circle"></i> <strong>Low similarity (' + sim.overall + '%)</strong> — distinct lineages.</div>';
  }
  verdict.innerHTML = vHtml; verdict.style.display = 'block';
  output.innerHTML =
    '<div class="compare-section"><h3>🧬 Typing Profile</h3><div class="compare-table">' +
    '<div class="compare-row header"><span class="row-label">Field</span>' +
    '<span class="row-val"><strong>' + escapeHtml(a) + '</strong></span>' +
    '<span class="row-val"><strong>' + escapeHtml(b) + '</strong></span></div>' +
    typingRows + '</div></div>' +
    '<div class="compare-section"><h3>💊 Gene Content</h3>' +
    (geneHtml || '<p><em>No gene data available.</em></p>') + '</div>' + fullHtml;
}

function swapCompare() {
  var a = document.getElementById('compare-a'), b = document.getElementById('compare-b');
  var t = a.value; a.value = b.value; b.value = t; runCompare();
}

function resetCompare() {
  var sel = document.getElementById('compare-a'); if (!sel) return;
  if (sel.options.length >= 2) {
    sel.value = sel.options[0].value;
    document.getElementById('compare-b').value = sel.options[1].value;
  }
  document.getElementById('compare-output').innerHTML = '';
  document.getElementById('compare-visuals').innerHTML = '';
  document.getElementById('compare-verdict').style.display = 'none';
}

function runCluster() {
  var simData = window.ECOLITYPER_SIM;
  var data = window.ECOLITYPER_COMPARE_DATA || {};
  var samples = simData ? simData.samples : Object.keys(data).sort();
  var n = samples.length;
  if (n < 2) {
    document.getElementById('cluster-heatmap').innerHTML = '<div class="compare-warning">Need at least two samples.</div>';
    return;
  }
  var thresholdEl = document.getElementById('cluster-threshold');
  var threshold = thresholdEl ? parseInt(thresholdEl.value) : 90;

  var parent = {};
  samples.forEach(function(s){ parent[s] = s; });
  function find(x) { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }
  function union(a, b) { var ra = find(a), rb = find(b); if (ra !== rb) parent[ra] = rb; }
  if (simData && simData.flat) {
    for (var i = 0; i < n; i++) {
      for (var j = i + 1; j < n; j++) {
        if (getPrecomputedSim(i, j) >= threshold) union(samples[i], samples[j]);
      }
    }
  } else {
    for (var i = 0; i < n; i++) {
      for (var j = i + 1; j < n; j++) {
        if (computeSimilarity(data[samples[i]], data[samples[j]]).overall >= threshold) union(samples[i], samples[j]);
      }
    }
  }
  var groups = {};
  samples.forEach(function(s){ var r = find(s); (groups[r] = groups[r] || []).push(s); });
  var clusters = Object.keys(groups).map(function(k){ return groups[k]; })
                       .filter(function(g){ return g.length >= 2; })
                       .sort(function(a, b){ return b.length - a.length; });
  var clusterHtml = '';
  if (clusters.length === 0) {
    clusterHtml = '<div class="no-alerts">No clusters detected at &gt;=' + threshold + '%.</div>';
  } else {
    clusters.forEach(function(cluster, idx) {
      var sum = 0, count = 0;
      for (var a = 0; a < cluster.length; a++) {
        for (var b = a + 1; b < cluster.length; b++) {
          var ia = samples.indexOf(cluster[a]);
          var ib = samples.indexOf(cluster[b]);
          sum += getPrecomputedSim(ia, ib);
          count++;
        }
      }
      var mean = count ? (sum / count).toFixed(1) : '100';
      var color = similarityColor(parseFloat(mean));
      clusterHtml += '<div class="cluster-card" style="border-left-color:' + color + ';">' +
        '<div class="cluster-header"><span class="cluster-badge" style="background:' + color + ';">Cluster ' + (idx + 1) + '</span>' +
        '<span class="cluster-stats">' + cluster.length + ' samples &middot; mean ' + mean + '% similar</span></div>' +
        '<div class="cluster-samples">' + cluster.map(function(s){ return '<span class="cluster-sample-pill">' + escapeHtml(s) + '</span>'; }).join('') + '</div></div>';
    });
  }
  document.getElementById('cluster-clusters').innerHTML = clusterHtml;

  if (n <= 60) {
    renderHeatmapDOM(samples, n);
  } else {
    renderHeatmapCanvas(samples, n);
  }

  var s = document.getElementById('cluster-summary');
  if (s) s.textContent = clusters.length + ' cluster(s) at >=' + threshold + '% (' + n + ' samples)';
}

function renderHeatmapDOM(samples, n) {
  var rowLabelWidth = 220;
  var cellSize = 28;
  var truncate = function(s, max) { return s.length > max ? s.substring(0, max - 1) + '…' : s; };
  var html = '<div class="heatmap-wrapper">' +
    '<h3>🔷 Pairwise Similarity Matrix</h3>' +
    '<p style="font-size:0.9em;color:#64748b;margin-bottom:12px;">' +
    'Hover any cell or label for the full sample name.</p>' +
    '<div style="overflow-x:auto;padding-bottom:10px;">' +
    '<div class="heatmap-grid" style="grid-template-columns: ' +
    rowLabelWidth + 'px repeat(' + n + ', ' + cellSize + 'px);">';
  html += '<div class="heatmap-cell heatmap-corner"></div>';
  samples.forEach(function(s) {
    html += '<div class="heatmap-cell heatmap-col-header" data-tip="' + escapeHtml(s) + '">' +
            escapeHtml(s) + '</div>';
  });
  samples.forEach(function(s1, i) {
    html += '<div class="heatmap-cell heatmap-row-header" data-tip="' + escapeHtml(s1) + '" style="justify-content:flex-start;padding-left:8px;width:' +
            rowLabelWidth + 'px;">' + escapeHtml(truncate(s1, 26)) + '</div>';
    samples.forEach(function(s2, j) {
      var pct = (i === j) ? 100 : getPrecomputedSim(i, j);
      var color = (i === j) ? '#1e293b' : similarityColor(pct);
      var text = (i === j) ? '—' : Math.round(pct);
      var tip = (i === j) ? escapeHtml(s1) + ' (self)' :
                escapeHtml(s1) + ' vs ' + escapeHtml(s2) + ' → ' + Math.round(pct) + '%';
      html += '<div class="heatmap-cell" data-tip="' + tip + '" style="background:' + color +
              ';color:white;">' + text + '</div>';
    });
  });
  html += '</div></div></div>';
  document.getElementById('cluster-heatmap').innerHTML = html;
}

function renderHeatmapCanvas(samples, n) {
  var cellSize = Math.max(3, Math.min(12, Math.floor(700 / n)));
  var rowLabelW = 160;
  var headerH = 60;
  var canvasW = rowLabelW + n * cellSize + 20;
  var canvasH = headerH + n * cellSize + 20;
  var canvas = document.createElement('canvas');
  canvas.width = canvasW;
  canvas.height = canvasH;
  canvas.style.display = 'block';
  canvas.style.maxWidth = '100%';
  canvas.style.cursor = 'crosshair';
  var ctx = canvas.getContext('2d');
  ctx.font = Math.max(8, Math.min(11, cellSize - 1)) + 'px sans-serif';
  ctx.textBaseline = 'middle';
  for (var i = 0; i < n; i++) {
    for (var j = 0; j < n; j++) {
      var pct = (i === j) ? 100 : getPrecomputedSim(i, j);
      ctx.fillStyle = (i === j) ? '#1e293b' : similarityColor(pct);
      ctx.fillRect(rowLabelW + j * cellSize, headerH + i * cellSize, cellSize, cellSize);
    }
  }
  ctx.fillStyle = '#333';
  ctx.textAlign = 'right';
  var rowStep = Math.max(1, Math.floor(n / 40));
  for (var i = 0; i < n; i += rowStep) {
    var lb = samples[i].length > 20 ? samples[i].substring(0, 19) + '…' : samples[i];
    ctx.fillText(lb, rowLabelW - 5, headerH + i * cellSize + cellSize / 2);
  }
  ctx.save();
  ctx.translate(0, headerH - 5);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'left';
  var colStep = Math.max(1, Math.floor(n / 40));
  for (var j = 0; j < n; j += colStep) {
    var lb2 = samples[j].length > 15 ? samples[j].substring(0, 14) + '…' : samples[j];
    ctx.fillText(lb2, 5, rowLabelW + j * cellSize + cellSize / 2);
  }
  ctx.restore();
  canvas.addEventListener('click', function(e) {
    var rect = canvas.getBoundingClientRect();
    var x = (e.clientX - rect.left) * (canvas.width / rect.width) - rowLabelW;
    var y = (e.clientY - rect.top) * (canvas.height / rect.height) - headerH;
    var jj = Math.floor(x / cellSize);
    var ii = Math.floor(y / cellSize);
    if (ii >= 0 && ii < n && jj >= 0 && jj < n && ii !== jj) {
      document.getElementById('compare-a').value = samples[ii];
      document.getElementById('compare-b').value = samples[jj];
      setCompareMode('pair');
      runCompare();
    }
  });
  var wrapper = document.getElementById('cluster-heatmap');
  wrapper.innerHTML =
    '<div class="heatmap-wrapper"><h3>🔷 Pairwise Similarity Matrix</h3>' +
    '<p style="font-size:0.9em;color:#64748b;margin-bottom:12px;">' +
    '🎨 Canvas-rendered for ' + n + ' samples. ' +
    '<strong>Click any cell to open pair compare for those two samples.</strong></p>' +
    '<div id="canvas-heatmap-host" style="overflow-x:auto;"></div></div>';
  wrapper.querySelector('#canvas-heatmap-host').appendChild(canvas);
}

document.addEventListener('DOMContentLoaded', function() {
  var hash = window.location.hash.substring(1);
  var target = hash ? document.querySelector('.tab-button.' + hash) : null;
  if (target) {
    switchTab(hash);
  } else {
    var first = document.querySelector('.tab-button');
    if (first) {
      var firstCls = null;
      for (var i = 0; i < first.classList.length; i++) {
        var c = first.classList[i];
        if (c !== 'tab-button' && c !== 'active') { firstCls = c; break; }
      }
      if (firstCls) switchTab(firstCls);
    }
  }
  var tables = document.querySelectorAll('.data-table');
  Array.prototype.forEach.call(tables, function(table) {
    var headers = table.querySelectorAll('th');
    Array.prototype.forEach.call(headers, function(th, idx) {
      if (th.classList.contains('col-check')) return;
      var type = th.getAttribute('data-sort') || 'string';
      th.style.cursor = 'pointer';
      th.addEventListener('click', function() { sortTable(table.id, idx, type); });
    });
  });
  initColumnToggle();
  initRowSelection();
});
</script>"""
        return js.replace('__TYPING__', typing_json)

    # ==================================================================
    # Report assembly
    # ==================================================================
    def generate_main_report(self, integrated_data: Dict[str, Any], output_dir: Path) -> str:
        print("\nGenerating ULTIMATE E. coli HTML report...")
        samples_data = integrated_data.get('samples', {})
        patterns = integrated_data.get('patterns', {})
        gene_centric = integrated_data.get('gene_centric', {})
        metadata = integrated_data.get('metadata', {})
        self._current_samples_data = samples_data

        typing_js = {}
        for s, d in samples_data.items():
            typing_js[s] = {
                "ST": d.get('mlst', {}).get('ST', 'ND'),
                "Serotype": d.get('serotype', {}).get('Serotype', 'ND'),
                "Phylogroup": d.get('phylogrouping', {}).get('Clermont_Type', 'ND'),
                "CHtype": d.get('chtyper', {}).get('CH_Type', 'ND'),
            }

        css = self._get_css()
        js = self._get_js(json.dumps(typing_js))

        TABS = [
            ('summary', 'Summary', 'chart-pie'),
            ('samples', 'Sample Overview', 'list-alt'),
            ('qc', 'FASTA QC', 'chart-line'),
            ('mlst', 'MLST', 'code-branch'),
            ('serotype', 'Serotype', 'tag'),
            ('chtype', 'CH Type', 'project-diagram'),
            ('phylogroup', 'Phylogroup', 'sitemap'),
            ('pathotype', 'Pathotype', 'virus'),
            ('amr', 'AMR', 'biohazard'),
            ('virulence', 'Virulence', 'virus'),
            ('plasmid', 'Plasmids & Serotyping', 'dna'),
            ('bacmet', 'BACMET', 'flask'),
            ('mutation', 'Mutations', 'dna'),
            ('mge', 'MGE', 'mobile-alt'),
            ('compare', 'Compare & Cluster', 'balance-scale'),
            ('patterns', 'Patterns', 'project-diagram'),
            ('databases', 'Database Metrics', 'database'),
            ('aiguide', 'AI Guide', 'robot'),
            ('citation', 'Citation', 'book'),
            ('funding', 'Funding', 'coffee'),
            ('calltoaction', 'Call to Action', 'globe'),
            ('export', 'Export', 'download'),
        ]

        nav_html = ''
        for i, (tid, title, icon) in enumerate(TABS):
            active = ' active' if i == 0 else ''
            nav_html += (f'<button class="tab-button {tid}{active}" '
                         f'onclick="switchTab(\'{tid}\')">'
                         f'<i class="fas fa-{icon}"></i> {title}</button>')

        section_methods = {
            'summary': self._sec_summary,
            'samples': self._sec_samples,
            'qc': self._sec_qc,
            'mlst': self._sec_mlst,
            'serotype': self._sec_serotype,
            'chtype': self._sec_chtype,
            'phylogroup': self._sec_phylogroup,
            'pathotype': self._sec_pathotype,
            'amr': self._sec_amr,
            'virulence': self._sec_virulence,
            'plasmid': self._sec_plasmid,
            'bacmet': self._sec_bacmet,
            'mutation': self._sec_mutation,
            'mge': self._sec_mge,
            'compare': self._sec_compare,
            'patterns': self._sec_patterns,
            'databases': self._sec_databases,
            'aiguide': lambda d: self._sec_aiguide(),
            'citation': lambda d: self._sec_citation(),
            'funding': lambda d: self._sec_funding(),
            'calltoaction': lambda d: self._sec_calltoaction(),
            'export': lambda d: self._sec_export(),
        }

        tabs_html = ''
        for i, (tid, title, icon) in enumerate(TABS):
            active = ' active' if i == 0 else ''
            content = section_methods[tid](integrated_data)
            tabs_html += f'''
            <div id="{tid}-tab" class="tab-content{active}">
                <h2 class="section-header {tid}-header">
                    <span><i class="fas fa-{icon}"></i> {title}</span>
                    <button class="print-section-btn" onclick="printSection('{tid}-tab')">
                        <i class="fas fa-print"></i> Print</button>
                </h2>
                {content}
            </div>'''

        total = len(samples_data)
        n_st = len(patterns.get('st_distribution', {}))
        n_sero = len(patterns.get('serotype_distribution', {}))
        n_phylo = len(patterns.get('phylogroup_distribution', {}))
        n_patho = len([k for k in patterns.get('pathotype_distribution', {}).keys() if k != 'none'])
        n_amr = sum(len(g) for g in gene_centric.get('amr_databases', {}).values())
        n_vir = sum(len(g) for g in gene_centric.get('virulence_databases', {}).values())
        n_high = len(patterns.get('high_risk_combinations', []))

        dash = f'''
        <div class="dashboard-grid">
            <div class="dashboard-card" onclick="switchTab('summary')" style="border-left-color:#4CAF50;">
                <i class="fas fa-vial fa-2x" style="color:#4CAF50;"></i>
                <div class="card-number">{total}</div><div class="card-label">Total Samples</div></div>
            <div class="dashboard-card" onclick="switchTab('mlst')" style="border-left-color:#FF9800;">
                <i class="fas fa-code-branch fa-2x" style="color:#FF9800;"></i>
                <div class="card-number">{n_st}</div><div class="card-label">Unique STs</div></div>
            <div class="dashboard-card" onclick="switchTab('serotype')" style="border-left-color:#9C27B0;">
                <i class="fas fa-tag fa-2x" style="color:#9C27B0;"></i>
                <div class="card-number">{n_sero}</div><div class="card-label">Serotypes</div></div>
            <div class="dashboard-card" onclick="switchTab('phylogroup')" style="border-left-color:#795548;">
                <i class="fas fa-sitemap fa-2x" style="color:#795548;"></i>
                <div class="card-number">{n_phylo}</div><div class="card-label">Phylogroups</div></div>
            <div class="dashboard-card" onclick="switchTab('pathotype')" style="border-left-color:#E91E63;">
                <i class="fas fa-virus fa-2x" style="color:#E91E63;"></i>
                <div class="card-number">{n_patho}</div><div class="card-label">Pathotypes</div></div>
            <div class="dashboard-card" onclick="switchTab('amr')" style="border-left-color:#F44336;">
                <i class="fas fa-biohazard fa-2x" style="color:#F44336;"></i>
                <div class="card-number">{n_amr}</div><div class="card-label">AMR Genes</div></div>
            <div class="dashboard-card" onclick="switchTab('virulence')" style="border-left-color:#8BC34A;">
                <i class="fas fa-virus fa-2x" style="color:#8BC34A;"></i>
                <div class="card-number">{n_vir}</div><div class="card-label">Virulence Genes</div></div>
            <div class="dashboard-card" onclick="switchTab('patterns')" style="border-left-color:#FF5722;">
                <i class="fas fa-project-diagram fa-2x" style="color:#FF5722;"></i>
                <div class="card-number">{n_high}</div><div class="card-label">High-Risk Combos</div></div>
        </div>'''

        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>EcoliTyper Ultimate Report v4.0.0</title>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
{css}
{js}
</head>
<body>
<div class="container">
  <div class="main-header">
    <h1><i class="fas fa-dna"></i> EcoliTyper Ultimate Analysis Report</h1>
    <p>Gene-Centric Cross-Genome Analysis · Pathotype Classification · MGE Profiling · Cluster Detection</p>
    <div class="metadata-bar">
      <div class="metadata-item"><i class="fas fa-calendar"></i><span>Generated: {metadata.get('analysis_date','Unknown')}</span></div>
      <div class="metadata-item"><i class="fas fa-database"></i><span>Samples: {total}</span></div>
      <div class="metadata-item"><i class="fas fa-user-md"></i><span>EcoliTyper Ultimate v4.0.0</span></div>
      <div class="metadata-item"><i class="fas fa-university"></i><span>University of Ghana Medical School</span></div>
    </div>
  </div>
  {dash}
  <div class="tab-navigation">{nav_html}</div>
  {tabs_html}
  <div class="footer">
    <h3>EcoliTyper Ultimate Reporter v4.0.0</h3>
    <p>University of Ghana Medical School | Brown Beckley &lt;brownbeckley94@gmail.com&gt;</p>
    <p>If you find this tool helpful, please <a href="https://github.com/bbeckley-hub/ecolityper" target="_blank">⭐ star us on GitHub</a>.</p>
    <p>Generated on {metadata.get('analysis_date','Unknown')}</p>
  </div>
</div>
</body>
</html>'''

        output_file = output_dir / "genius_ecoli_ultimate_gene_centric_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"    HTML report saved: {output_file}")
        return str(output_file)

    # ==================================================================
    # Sections
    # ==================================================================
    def _sec_summary(self, d):
        samples = d['samples']
        patterns = d['patterns']
        gc = d['gene_centric']
        total = len(samples)
        amr_total = sum(len(g) for g in gc.get('amr_databases', {}).values())
        vir_total = sum(len(g) for g in gc.get('virulence_databases', {}).values())
        bac_total = sum(len(g) for g in gc.get('bacmet_databases', {}).values())
        pla_total = sum(len(g) for g in gc.get('plasmid_databases', {}).values())
        high_risk = len(patterns.get('high_risk_combinations', []))
        patho_counts = patterns.get('pathotype_distribution', Counter())
        mge_n = len(d.get('mge_data', {}).get('per_sample', {}))

        features = [
            ('Gene-Centric Tables', 'fa-gene',
             'Each resistance, virulence, BACMET, plasmid, or mutation is shown with <strong>all genomes</strong> that carry it.'),
            ('Dynamic Grouping by Typing', 'fa-layer-group',
             'Group by MLST, Serotype, Phylogroup, CH type, or any combination of the four.'),
            ('Pathotype Classification', 'fa-virus',
             'Computational STEC / EPEC / EHEC / EAEC / ETEC / DAEC / EIEC prediction with prevalence tables.'),
            ('MGE Profiling', 'fa-mobile-alt',
             'mobileOG-db hit counts per sample, plus aggregate profiles by typing.'),
            ('Compare &amp; Cluster Detection', 'fa-balance-scale',
             'Pair-compare with similarity gauge; cluster detection via union-find on a precomputed pairwise similarity matrix — scales to 1000+ genomes with canvas heatmaps.'),
            ('In-Box Scrolling with Sticky Headers', 'fa-arrows-alt-v',
             'Every large table scrolls inside its container while the header row stays fixed.'),
            ('On-Demand Educational Expanders', 'fa-book-open',
             'Interpretation guides and database explanations are collapsed by default and open below each table.'),
            ('Section-Specific Printing', 'fa-print',
             'Print any tab individually with one click.'),
            ('Full Data Export', 'fa-download',
             'All tables as CSV; complete dataset as JSON for AI-assisted analysis.'),
        ]
        feature_tiles = ''
        for label, icon, desc in features:
            feature_tiles += f'<div class="database-section"><h4><i class="fas {icon}"></i> {label}</h4><p>{desc}</p></div>'

        return f'''
        {self._alert('info', 'fa-info-circle',
            f'<h3>📊 Analysis Overview — Gene-Centric E. coli Report v4.0.0</h3>'
            f'<p>Comprehensive genomic analysis of <strong>{total}</strong> <em>Escherichia coli</em> genomes using a '
            f'<strong>gene-centric</strong> approach: each resistance or virulence gene is displayed with all genomes '
            f'that carry it — making outbreak tracking and co-occurrence analysis immediate and transparent.</p>'
            f'<p>Includes MLST, serotype (O:H), phylogroup, CH type, pathotype, AMR, virulence, BACMET, plasmids, '
            f'point mutations, and mobile genetic element profiles.</p>')}

        <div class="action-buttons">
          <button class="action-btn btn-danger" onclick="switchTab('amr')"><i class="fas fa-biohazard"></i> AMR</button>
          <button class="action-btn btn-success" onclick="switchTab('virulence')"><i class="fas fa-virus"></i> Virulence</button>
          <button class="action-btn btn-warning" onclick="switchTab('pathotype')"><i class="fas fa-chart-pie"></i> Pathotypes</button>
          <button class="action-btn btn-info" onclick="switchTab('compare')"><i class="fas fa-balance-scale"></i> Compare</button>
          <button class="action-btn btn-secondary" onclick="switchTab('mge')"><i class="fas fa-mobile-alt"></i> MGE</button>
        </div>

        {self._alert('success', 'fa-magic',
            '<h3>📘 How to Use Grouping</h3>'
            '<ol style="margin-left:20px;">'
            '<li>Go to any gene-centric table (AMR, Virulence, BACMET, Plasmids, Mutations).</li>'
            '<li>Use the <strong>"Group genomes by"</strong> dropdown above the table.</li>'
            '<li>Pick any single field or combination — the genome column re-organises instantly.</li>'
            '<li>Click <strong>"Reset"</strong> to return to a flat list.</li>'
            '</ol>'
            '<p style="margin-top:8px;"><strong>Why it matters:</strong> you can instantly ask <em>"Does blaCTX-M only appear in ST131?"</em> or <em>"Which phylogroups carry Shiga toxins?"</em></p>')}

        <h3><i class="fas fa-chart-bar"></i> Key Statistics</h3>
        <div class="scrollable-table"><table class="data-table">
          <thead><tr><th>Metric</th><th>Count</th><th>Details</th></tr></thead>
          <tbody>
            <tr><td>Total Samples</td><td><strong>{total}</strong></td><td>Complete E. coli genomes</td></tr>
            <tr><td>Unique STs</td><td><strong>{len(patterns.get("st_distribution", {}))}</strong></td><td>MLST sequence types (population structure)</td></tr>
            <tr><td>Unique Serotypes</td><td><strong>{len(patterns.get("serotype_distribution", {}))}</strong></td><td>O:H antigens</td></tr>
            <tr><td>Unique Phylogroups</td><td><strong>{len(patterns.get("phylogroup_distribution", {}))}</strong></td><td>Clermont scheme</td></tr>
            <tr><td>Unique CH Types</td><td><strong>{len(patterns.get("ch_type_distribution", {}))}</strong></td><td>FumC:FimH fine-scale typing</td></tr>
            <tr><td>AMR Genes</td><td><strong>{amr_total}</strong></td><td>Across AMRFinder + ABRicate databases</td></tr>
            <tr><td>Virulence Genes</td><td><strong>{vir_total}</strong></td><td>VFDB + EcoLi_VF</td></tr>
            <tr><td>BACMET Genes</td><td><strong>{bac_total}</strong></td><td>Biocide / heavy-metal resistance</td></tr>
            <tr><td>Plasmid / Serotyping Markers</td><td><strong>{pla_total}</strong></td><td>PlasmidFinder + EcoH</td></tr>
            <tr><td>Pathotypes Detected</td><td><strong>{len(patho_counts)}</strong></td><td>{', '.join(k for k in patho_counts if k != "none") or "None"}</td></tr>
            <tr><td>MGE-profiled Samples</td><td><strong>{mge_n}</strong></td><td>mobileOG-db hit counts</td></tr>
            <tr><td>High-Risk AMR + Virulence Combos</td><td><span class="badge" style="background:#8B0000;color:white;">{high_risk}</span></td><td>Critical resistance + critical virulence</td></tr>
          </tbody>
        </table></div>

        <h3 style="margin-top:30px;"><i class="fas fa-lightbulb"></i> ✨ Key Features</h3>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:20px;margin:20px 0;">
          {feature_tiles}
        </div>'''

    def _sec_samples(self, d):
        samples = d['samples']
        rows = ''
        for s, data in sorted(samples.items()):
            mlst = data.get('mlst', {})
            sero = data.get('serotype', {})
            phylo = data.get('phylogrouping', {})
            cht = data.get('chtyper', {})
            vir_genes = data.get('virulence_genes', [])
            patho = data.get('pathotype', 'none')
            vir_tags = ''.join(f'<span class="genome-tag">{g}</span>' for g in sorted(set(vir_genes)))
            vir_details = (f'<details class="vir-details"><summary>{len(vir_genes)} gene(s)</summary>'
                           f'<div class="genome-list">{vir_tags or "<em>None</em>"}</div></details>')
            patho_badge = (f'<span class="badge" style="background:#e91e63;color:white;">{patho}</span>'
                           if patho != 'none' else '<span class="badge" style="background:#e0e0e0;color:#555;">none</span>')
            rows += f'''<tr>
                <td class="col-check"><input type="checkbox" class="row-check"></td>
                <td><strong>{s}</strong></td>
                <td>{mlst.get('ST','ND')}</td>
                <td>{sero.get('Serotype','ND')}</td>
                <td>{sero.get('O_Type','ND')}</td>
                <td>{sero.get('H_Type','ND')}</td>
                <td>{sero.get('Serotype_Confidence','ND')}</td>
                <td>{phylo.get('Clermont_Type','ND')}</td>
                <td>{cht.get('FumC_Type','ND')}</td>
                <td>{cht.get('FimH_Type','ND')}</td>
                <td>{cht.get('CH_Type','ND')}</td>
                <td>{patho_badge}</td>
                <td>{vir_details}</td>
            </tr>'''
        column_guide = self._collapsible(
            '📖 What Each Column Means', '#2196F3',
            '<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:18px;">'
            '<div style="background:#fff;padding:16px;border-radius:8px;border-left:4px solid #FF9800;">'
            '<strong style="color:#FF9800;">🔬 MLST (Multi-Locus Sequence Typing)</strong>'
            '<p style="font-size:.92em;margin-top:6px;">Seven housekeeping genes define a Sequence Type (ST). '
            'The gold standard for global <em>E. coli</em> epidemiology. Examples: ST131 (pandemic UPEC/ESBL), '
            'ST10 (dominant commensal), ST11 (O157:H7).</p></div>'
            '<div style="background:#fff;padding:16px;border-radius:8px;border-left:4px solid #9C27B0;">'
            '<strong style="color:#9C27B0;">🧬 Serotype (O:H)</strong>'
            '<p style="font-size:.92em;margin-top:6px;">O = LPS antigen; H = flagellar antigen. '
            'O157:H7 (EHEC), O104:H4 (2011 outbreak), O25b:H4 (ST131). Used worldwide for food-safety surveillance.</p></div>'
            '<div style="background:#fff;padding:16px;border-radius:8px;border-left:4px solid #795548;">'
            '<strong style="color:#795548;">🌿 Phylogroup</strong>'
            '<p style="font-size:.92em;margin-top:6px;">Deep evolutionary lineages (A, B1, B2, C, D, E, F, G). '
            'B2 and D dominate extra-intestinal infections; A and B1 are mostly commensals.</p></div>'
            '<div style="background:#fff;padding:16px;border-radius:8px;border-left:4px solid #009688;">'
            '<strong style="color:#009688;">🧬 CH Type (FumC:FimH)</strong>'
            '<p style="font-size:.92em;margin-top:6px;">Higher resolution than MLST within a clonal complex. '
            'Two ST131 isolates can carry different <em>fimH</em> alleles and be epidemiologically distinct.</p></div>'
            '<div style="background:#fff;padding:16px;border-radius:8px;border-left:4px solid #E91E63;">'
            '<strong style="color:#E91E63;">💊 Pathotype</strong>'
            '<p style="font-size:.92em;margin-top:6px;">Computational prediction from VFDB gene content. '
            'STEC/EHEC are food-safety priorities; UPEC/ExPEC dominate urinary and bloodstream infections.</p></div>'
            '</div>')
        why_matters = self._collapsible(
            '💡 Why This Matters for AMR Research', '#00695c',
            '<ul style="margin-top:8px;line-height:1.7;">'
            '<li><strong>Lineage–resistance linkage.</strong> ESBL and carbapenemase genes cluster in specific clones — ST131, ST10, ST405, ST648, ST167.</li>'
            '<li><strong>Outbreak detection.</strong> Shared <strong>ST + Serotype + Phylogroup + CH type</strong> strongly suggests recent transmission.</li>'
            '<li><strong>Surveillance.</strong> Which STs and phylogroups dominate in your setting? Do they match regional trends?</li>'
            '<li><strong>Horizontal transfer signals.</strong> The same AMR gene across unrelated STs implies plasmid or transposon spread.</li>'
            '</ul>')
        return f'''
        {self._alert('info', 'fa-list-alt',
            '<h3>🧬 Sample Overview — The Population Snapshot</h3>'
            '<p>This is the <strong>master reference table</strong>. Every downstream tab breaks down one slice of the dataset — '
            'this one shows all typing layers for every isolate on a single row.</p>')}

        <div class="action-buttons" style="margin-bottom:10px;">
          <button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table','sample_overview.csv')">
            <i class="fas fa-download"></i> Export CSV</button>
          <div class="col-toggle-wrap">
            <button class="action-btn btn-secondary" id="colToggleBtn">
              <i class="fas fa-cog"></i> Columns</button>
            <div class="col-panel" id="colPanel"></div>
          </div>
          <span class="sel-counter" id="selCount">0 selected</span>
        </div>

        <input type="text" class="search-box" id="search-samples"
               onkeyup="searchTable('samples-table','search-samples')"
               placeholder="🔍 Search samples by ID, ST, serotype, phylogroup, CH type, pathotype...">

        <div class="scrollable-table">
          <table id="samples-table" class="data-table">
            <thead><tr>
              <th class="col-check"><input type="checkbox" id="selectAll"></th>
              <th data-sort="string">Sample</th>
              <th data-sort="string">ST</th>
              <th data-sort="string">Serotype</th>
              <th data-sort="string">O Type</th>
              <th data-sort="string">H Type</th>
              <th data-sort="string">Serotype Conf.</th>
              <th data-sort="string">Phylogroup</th>
              <th data-sort="string">FumC</th>
              <th data-sort="string">FimH</th>
              <th data-sort="string">CH Type</th>
              <th data-sort="string">Pathotype</th>
              <th data-sort="number">Virulence</th>
            </tr></thead>
            <tbody>{rows}</tbody>
          </table>
        </div>

        <h3 style="margin-top:30px;"><i class="fas fa-book-medical"></i> Interpretation Guide</h3>
        {column_guide}
        {why_matters}'''

    def _sec_qc(self, d):
        qc = d.get('qc_data', {})
        if not qc:
            return self._alert('warning', 'fa-exclamation-circle',
                               '<h3>No QC Data Available</h3><p>FASTA_QC_summary.html was not found.</p>')
        metrics = sorted({m for v in qc.values() for m in v.keys()})
        header = ''.join(f'<th data-sort="number">{m}</th>' for m in metrics)
        rows = ''
        for s, vals in sorted(qc.items()):
            cells = ''
            for m in metrics:
                v = vals.get(m, 'ND')
                if isinstance(v, float):
                    v = f"{v:,.0f}" if v > 1000 else f"{v:.2f}"
                cells += f'<td>{v}</td>'
            rows += f'<tr><td><strong>{s}</strong></td>{cells}</tr>'
        guide = self._collapsible(
            '📏 FASTA Quality Control Guide', '#17a2b8',
            '<ul style="margin-top:8px;">'
            '<li><strong>Number of contigs</strong> – lower is better; typical E. coli: &lt;200 (good), &lt;100 (excellent).</li>'
            '<li><strong>N50</strong> – higher is better; &gt;100 kb excellent.</li>'
            '<li><strong>GC%</strong> – E. coli is typically ~50.8%.</li>'
            '<li><strong>Total length</strong> – ~4.6 Mb (commensal) to ~5.5 Mb (EHEC).</li>'
            '</ul>')
        return f'''
        {self._credit_bar('#17a2b8', '📊', 'FASTA QC Metrics – Tools & Acknowledgements',
            '<strong>Biopython</strong> by <a href="https://biopython.org/" target="_blank" style="color:#17a2b8;font-weight:bold;">The Biopython Consortium</a> → '
            'assembly statistics (contigs, N50, GC%, total length) from FASTA files.<br>'
            '<strong>FastANI</strong> by <a href="https://github.com/ParBLiSS/FastANI" target="_blank" style="color:#17a2b8;font-weight:bold;">ParBLiSS Lab</a> → '
            'average nucleotide identity (ANI) for species-level confirmation.')}
        {self._alert('info', 'fa-chart-line',
            '<h3>FASTA Quality Control</h3>'
            '<p>Assembly quality metrics for each genome. Poor assembly quality can lead to missing genes or false positives — '
            'always inspect QC before drawing conclusions.</p>')}
        <input type="text" class="search-box" id="search-qc"
               onkeyup="searchTable('qc-table','search-qc')" placeholder="🔍 Search sample...">
        <div class="action-buttons">
          <button class="action-btn btn-primary" onclick="exportTableToCSV('qc-table','fasta_qc.csv')">
            <i class="fas fa-download"></i> Export QC</button>
        </div>
        <div class="master-scrollable-container">
          <table id="qc-table" class="data-table">
            <thead><tr><th data-sort="string">Sample</th>{header}</tr></thead>
            <tbody>{rows}</tbody>
          </table>
        </div>
        <h3 style="margin-top:30px;"><i class="fas fa-book-medical"></i> Interpretation Guide</h3>
        {guide}'''

    def _sec_mlst(self, d):
        P = d['patterns']
        credit = self._credit_bar('#FF9800', '🧬', 'MLST Typing – Acknowledgments & Licensing',
            '<strong>MLST scheme</strong> powered by '
            '<a href="https://github.com/tseemann/mlst" target="_blank" style="color:#FF9800;font-weight:bold;">Prof. Torsten Seemann\'s scripts</a> '
            'and the <a href="https://pubmlst.org/" target="_blank" style="color:#FF9800;font-weight:bold;">PubMLST database</a> '
            '(Jolley et al., <em>Wellcome Open Res</em> 2018).<br>'
            '<span style="color:#856404;"><i class="fas fa-info-circle"></i> <strong>Note:</strong> '
            'Due to recent licensing changes and the new PubMLST data policy, allele definitions used in this report '
            'are current as of <strong>2024</strong>. Check the pipeline <code>--help</code> for update instructions.</span><br>'
            '<span style="font-size:.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We are grateful to the PubMLST curators and Torsten Seemann for providing these essential open-source tools.</span>')
        st_dist = P.get('st_distribution', Counter())
        total = sum(st_dist.values())
        rows = ''
        for st, cnt in st_dist.most_common():
            pct = cnt / total * 100 if total else 0
            seros = set()
            for combo in P.get('st_oh_combinations', {}):
                if combo.startswith(f"{st} - "):
                    seros.add(combo.split(" - ", 1)[1])
            sero_str = ', '.join(sorted(seros)) if seros else 'ND'
            rows += f'<tr><td><strong>{st}</strong></td><td>{cnt}</td><td>{pct:.1f}%</td><td>{sero_str}</td></tr>'
        return f'''
        {credit}
        {self._alert('info', 'fa-code-branch',
            '<h3>🔬 MLST (Multi-Locus Sequence Typing)</h3>'
            '<p>MLST indexes internal fragments of seven housekeeping genes. Each unique combination defines a Sequence Type (ST) — '
            'the gold standard for global <em>E. coli</em> epidemiology. Closely related STs belong to the same Clonal Complex (CC).</p>'
            f'<p><strong>{len(st_dist)} unique STs</strong> identified across the dataset.</p>')}
        <h3>📊 ST Distribution</h3>
        <div class="master-scrollable-container"><table id="mlst-dist" class="data-table">
          <thead><tr><th data-sort="string">ST</th><th data-sort="number">Count</th>
                     <th data-sort="number">%</th><th data-sort="string">Associated Serotypes</th></tr></thead>
          <tbody>{rows}</tbody></table></div>
        {self._combo_table("mlst-o-table", "ST – O antigen", P.get("st_o_combinations", {}))}
        {self._combo_table("mlst-h-table", "ST – H antigen", P.get("st_h_combinations", {}))}
        {self._combo_table("mlst-oh-table", "ST – Serotype (O:H)", P.get("st_oh_combinations", {}))}
        {self._combo_table("mlst-phylo-table", "ST – Phylogroup", P.get("st_phylogroup_combinations", {}))}
        {self._combo_table("mlst-ch-table", "ST – CH type (FumC:FimH)", P.get("st_fumc_fimh_combinations", {}))}
        {self._combo_table("mlst-full-table", "ST – O:H – Phylogroup – CH type (full combination)", P.get("st_oh_phylo_fumc_fimh_combinations", {}))}'''

    def _sec_serotype(self, d):
        P = d['patterns']
        samples = d['samples']
        credit = self._credit_bar('#9C27B0', '🔬', 'Serotype Analysis – Acknowledgments',
            'Powered by the <a href="https://cge.food.dtu.dk/services/SerotypeFinder/" target="_blank" style="color:#9C27B0;font-weight:bold;">CGE SerotypeFinder</a> database '
            '(Center for Genomic Epidemiology) and <strong>EcoH</strong> (Ingle et al., 2016).<br>'
            'Please cite: Joensen KG, et al. Rapid and easy in silico serotyping of <em>Escherichia coli</em> using whole genome sequencing (WGS) data. <em>J Clin Microbiol</em>. 2015;53(8):2410-2426. '
            '<a href="https://doi.org/10.1128/JCM.00008-15" target="_blank" style="color:#9C27B0;font-weight:bold;">doi:10.1128/JCM.00008-15</a><br>'
            '<span style="font-size:.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> We thank the CGE team, Flemming Scheutz, and SSI for maintaining this essential resource.</span>')
        dist = P.get('serotype_distribution', Counter())
        total = sum(dist.values())
        rows = ''
        for sero, cnt in dist.most_common():
            if sero == 'ND': continue
            pct = cnt / total * 100 if total else 0
            rows += f'<tr><td><strong>{sero}</strong></td><td>{cnt}</td><td>{pct:.1f}%</td></tr>'
        def simple_table(distribution, label):
            out = ''
            for v, c in distribution.most_common():
                if v == 'ND': continue
                smp = [s for s, x in samples.items() if x.get('serotype', {}).get(label) == v]
                tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in smp)
                out += f'<tr><td><strong>{v}</strong></td><td>{c}</td><td><div class="genome-list">{tags}</div></td></tr>'
            return out
        return f'''
        {credit}
        {self._alert('info', 'fa-tag',
            '<h3>🔬 Serotype Analysis (O:H)</h3>'
            '<p>Serotyping identifies lipopolysaccharide (O) and flagellar (H) antigens. Specific O:H combinations are associated with '
            'clinical syndromes and are used worldwide for food-safety and outbreak surveillance.</p>'
            f'<p><strong>{len(dist)} unique serotypes</strong> identified.</p>')}
        <h3>📊 Serotype Distribution</h3>
        <div class="master-scrollable-container"><table class="data-table">
          <thead><tr><th data-sort="string">Serotype</th><th data-sort="number">Count</th>
                     <th data-sort="number">%</th></tr></thead>
          <tbody>{rows}</tbody></table></div>
        <h3>📊 O Type Distribution</h3>
        <div class="master-scrollable-container"><table class="data-table">
          <thead><tr><th>O Type</th><th>Count</th><th>Samples</th></tr></thead>
          <tbody>{simple_table(P.get("o_type_distribution", Counter()), "O_Type")}</tbody></table></div>
        <h3>📊 H Type Distribution</h3>
        <div class="master-scrollable-container"><table class="data-table">
          <thead><tr><th>H Type</th><th>Count</th><th>Samples</th></tr></thead>
          <tbody>{simple_table(P.get("h_type_distribution", Counter()), "H_Type")}</tbody></table></div>'''

    def _sec_chtype(self, d):
        P = d['patterns']
        samples = d['samples']
        credit = self._credit_bar('#009688', '🧬', 'CHTyper – Acknowledgments',
            'Powered by the <a href="https://cge.food.dtu.dk/services/CHTyper/" target="_blank" style="color:#009688;font-weight:bold;">CGE CHTyper</a> '
            'database curated by <strong>Henrik Hasman (CGE)</strong>.<br>'
            'Please cite: Roer L, et al. CHTyper, a Web Tool for Subtyping of Extraintestinal Pathogenic <em>Escherichia coli</em> '
            'Based on the <em>fumC</em> and <em>fimH</em> Alleles. <em>J Clin Microbiol</em>. 2018;56(4):e00063-18. '
            '<a href="https://doi.org/10.1128/JCM.00063-18" target="_blank" style="color:#009688;font-weight:bold;">doi:10.1128/JCM.00063-18</a><br>'
            '<span style="font-size:.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> We thank the CGE team and Henrik Hasman for maintaining this resource.</span>')
        dist = P.get('ch_type_distribution', Counter())
        total = sum(dist.values())
        rows = ''
        for cht, cnt in dist.most_common():
            if cht == 'ND': continue
            pct = cnt / total * 100 if total else 0
            smp = [s for s, x in samples.items() if x.get('chtyper', {}).get('CH_Type') == cht]
            tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in smp)
            rows += f'<tr><td><strong>{cht}</strong></td><td>{cnt}</td><td>{pct:.1f}%</td><td><div class="genome-list">{tags}</div></td></tr>'
        return f'''
        {credit}
        {self._alert('info', 'fa-project-diagram',
            '<h3>🧬 CH Type (FumC:FimH)</h3>'
            '<p>CH typing combines alleles of the housekeeping gene <em>fumC</em> and the type-1 fimbrial adhesin <em>fimH</em>. '
            'It provides high resolution for population studies and has been used to define epidemic clones.</p>')}
        <div class="master-scrollable-container"><table class="data-table">
          <thead><tr><th data-sort="string">CH Type</th><th data-sort="number">Count</th>
                     <th data-sort="number">%</th><th data-sort="string">Samples</th></tr></thead>
          <tbody>{rows}</tbody></table></div>'''

    def _sec_phylogroup(self, d):
        P = d['patterns']
        samples = d['samples']
        credit = self._credit_bar('#795548', '🌿', 'Phylogrouping – Acknowledgments',
            'Powered by <a href="https://github.com/nickp60/EzClermont" target="_blank" style="color:#795548;font-weight:bold;">EzClermont</a> (nickp60) — '
            'a Python implementation of the Clermont phylotyping method.<br>'
            'Original method: Clermont O, et al. The Clermont <em>Escherichia coli</em> phylo-typing method revisited: improvement of specificity and detection of new phylo-groups. '
            '<em>Environ Microbiol Rep</em>. 2013;5(1):58-65. '
            '<a href="https://doi.org/10.1111/1758-2229.12019" target="_blank" style="color:#795548;font-weight:bold;">doi:10.1111/1758-2229.12019</a>')
        dist = P.get('phylogroup_distribution', Counter())
        total = sum(dist.values())
        rows = ''
        for ph, cnt in dist.most_common():
            if ph == 'ND': continue
            pct = cnt / total * 100 if total else 0
            smp = [s for s, x in samples.items() if x.get('phylogrouping', {}).get('Clermont_Type') == ph]
            tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in smp)
            rows += f'<tr><td><strong>{ph}</strong></td><td>{cnt}</td><td>{pct:.1f}%</td><td><div class="genome-list">{tags}</div></td></tr>'
        return f'''
        {credit}
        {self._alert('info', 'fa-sitemap',
            '<h3>🌿 Phylogroup Analysis (Clermont scheme)</h3>'
            '<p><em>E. coli</em> is divided into phylogroups (A, B1, B2, C, D, E, F, G) that reflect ecological niches and pathogenic potential. '
            'B2 and D are associated with extra-intestinal infections, while A and B1 are common in commensals and environmental strains.</p>')}
        <div class="master-scrollable-container"><table class="data-table">
          <thead><tr><th data-sort="string">Phylogroup</th><th data-sort="number">Count</th>
                     <th data-sort="number">%</th><th data-sort="string">Samples</th></tr></thead>
          <tbody>{rows}</tbody></table></div>'''

    def _sec_pathotype(self, d):
        P = d['patterns']
        samples = d['samples']
        patho_dist = P.get('pathotype_distribution', Counter())
        total = sum(patho_dist.values())
        rows = ''
        for p, c in patho_dist.most_common():
            pct = c / total * 100 if total else 0
            smp = [s for s, x in samples.items() if x.get('pathotype') == p]
            tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in smp[:200])
            rows += f'<tr><td><strong>{p}</strong></td><td>{c}</td><td>{pct:.1f}%</td><td><div class="genome-list">{tags}</div></td></tr>'
        sample_rows = ''
        for s, x in samples.items():
            sample_rows += f'''<tr>
                <td><strong>{s}</strong></td>
                <td>{x.get("pathotype", "none")}</td>
                <td>{x.get("mlst", {}).get("ST", "ND")}</td>
                <td>{x.get("serotype", {}).get("Serotype", "ND")}</td>
                <td>{x.get("phylogrouping", {}).get("Clermont_Type", "ND")}</td>
                <td>{x.get("chtyper", {}).get("CH_Type", "ND")}</td>
            </tr>'''
        threshold_notice = self._collapsible(
            '📏 Identity & Coverage Thresholds — Read This First', '#0891b2',
            '<p style="margin:6px 0;font-size:.95em;line-height:1.55;">'
            'Pathotype calls are derived from <strong>gene presence/absence</strong> in VFDB and EcoLi_VF, screened by ABRicate. '
            'This reporter uses whatever identity and coverage thresholds were used in the ABRicate run.</p>'
            '<p style="margin:10px 0 6px 0;font-size:.95em;"><strong>Default ABRicate thresholds:</strong></p>'
            '<ul style="margin:6px 0 6px 22px;font-size:.93em;">'
            '<li><strong>Minimum identity:</strong> <code>80%</code> (<code>--minid 80</code>)</li>'
            '<li><strong>Minimum coverage:</strong> <code>80%</code> (<code>--mincov 80</code>)</li>'
            '</ul>'
            '<p style="margin:10px 0 0 0;font-size:.92em;background:#fff3cd;padding:10px 14px;border-radius:4px;border-left:3px solid #ffc107;">'
            '<strong>Want stricter or looser calls?</strong> Re-run the ABRicate step with different <code>--minid</code> / '
            '<code>--mincov</code> flags, then regenerate this report.</p>')
        how_it_works = self._collapsible(
            '🔍 How Pathotype Is Called', '#1e3c72',
            '<table class="data-table" style="font-size:.9em;">'
            '<thead><tr><th style="width:100px;">Pathotype</th><th>Rule</th><th>Clinical context</th></tr></thead>'
            '<tbody>'
            '<tr><td><strong>STEC</strong></td><td><em>stx1</em> OR <em>stx2</em> present AND <em>eae</em> absent</td><td>Haemorrhagic colitis, HUS</td></tr>'
            '<tr><td><strong>EPEC</strong></td><td>(<em>eae</em> OR <em>bfpA</em>) present AND <em>stx</em> absent</td><td>Infantile diarrhoea</td></tr>'
            '<tr><td><strong>EHEC</strong></td><td><em>stx</em> present AND <em>eae</em> present</td><td>Bloody diarrhoea, HUS</td></tr>'
            '<tr><td><strong>EAEC</strong></td><td><em>aggR</em> OR <em>aaiC</em> OR <em>aatA</em> present</td><td>Persistent diarrhoea</td></tr>'
            '<tr><td><strong>ETEC</strong></td><td><em>lt</em> OR <em>sta1</em> present</td><td>Traveller\'s diarrhoea</td></tr>'
            '<tr><td><strong>DAEC</strong></td><td><em>afaC</em> OR <em>afaE</em> present</td><td>Paediatric diarrhoea</td></tr>'
            '<tr><td><strong>EIEC</strong></td><td><em>ipaH</em> present</td><td>Dysentery</td></tr>'
            '<tr><td><strong>none</strong></td><td>None of the above</td><td>No diarrhoeagenic pathotype predicted</td></tr>'
            '</tbody></table>'
            '<ul style="font-size:.9em;color:#555;margin:10px 0 0 22px;line-height:1.6;">'
            '<li>Pathotype is a <strong>computational prediction</strong> — confirm with phenotype.</li>'
            '<li>Markers are inferred at the current ABRicate thresholds (default 80/80).</li>'
            '<li>Upstream/downstream context (e.g. LEE island intactness) is not evaluated.</li>'
            '</ul>')
        return f'''
        {self._alert('info', 'fa-virus',
            '<h3>🧬 Pathotype Classification &amp; Prevalence</h3>'
            '<p>E. coli pathotypes are defined by specific combinations of virulence factors. '
            'This section shows the prevalence of each pathotype and lets you explore which samples belong to each category.</p>'
            '<p><strong>Definitions credit:</strong> '
            '<a href="https://github.com/kevinsblake/PathotypeR" target="_blank" style="color:#17a2b8;font-weight:bold;">Kevin Blake, PhD (PathotypeR)</a>.</p>')}

        <h3>📊 Pathotype Prevalence</h3>
        <div class="master-scrollable-container"><table class="data-table">
          <thead><tr><th data-sort="string">Pathotype</th><th data-sort="number">Count</th>
                     <th data-sort="number">%</th><th data-sort="string">Samples</th></tr></thead>
          <tbody>{rows}</tbody></table></div>

        <h3 style="margin-top:30px;">🧬 Sample Pathotype Details</h3>
        <div class="grouping-controls">
          <strong><i class="fas fa-filter"></i> Quick filters:</strong>
          <button class="group-btn" onclick="document.getElementById('search-patho').value='STEC'; searchTable('patho-samples','search-patho')">🔴 STEC</button>
          <button class="group-btn" onclick="document.getElementById('search-patho').value='EPEC'; searchTable('patho-samples','search-patho')">🟠 EPEC</button>
          <button class="group-btn" onclick="document.getElementById('search-patho').value='EHEC'; searchTable('patho-samples','search-patho')">🔴🔴 EHEC</button>
          <button class="group-btn" onclick="document.getElementById('search-patho').value='EAEC'; searchTable('patho-samples','search-patho')">🟣 EAEC</button>
          <button class="group-btn" onclick="document.getElementById('search-patho').value='ETEC'; searchTable('patho-samples','search-patho')">🔷 ETEC</button>
          <button class="group-btn" onclick="document.getElementById('search-patho').value='DAEC'; searchTable('patho-samples','search-patho')">🌸 DAEC</button>
          <button class="group-btn" onclick="document.getElementById('search-patho').value='EIEC'; searchTable('patho-samples','search-patho')">🟤 EIEC</button>
          <button class="group-btn" onclick="document.getElementById('search-patho').value='none'; searchTable('patho-samples','search-patho')">⚪ none</button>
          <button class="group-btn" onclick="document.getElementById('search-patho').value=''; searchTable('patho-samples','search-patho')">Clear</button>
        </div>
        <input type="text" class="search-box" id="search-patho"
               onkeyup="searchTable('patho-samples','search-patho')" placeholder="🔍 Search samples...">
        <div class="action-buttons">
          <button class="action-btn btn-primary" onclick="exportTableToCSV('patho-samples','pathotype_samples.csv')">
            <i class="fas fa-download"></i> Export CSV</button>
        </div>
        <div class="master-scrollable-container"><table id="patho-samples" class="data-table">
          <thead><tr><th>Sample</th><th>Pathotype</th><th>ST</th><th>Serotype</th>
                     <th>Phylogroup</th><th>CH Type</th></tr></thead>
          <tbody>{sample_rows}</tbody></table></div>

        <h3 style="margin-top:30px;"><i class="fas fa-book-medical"></i> Interpretation Guide</h3>
        {threshold_notice}
        {how_it_works}'''

    def _sec_amr(self, d):
        gc = d['gene_centric']
        amr_dbs = gc.get('amr_databases', {})
        all_genes = []
        for genes in amr_dbs.values():
            all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        total = len(d['samples'])

        credit = self._credit_bar('#F44336', '🧬', 'AMR Detection Tools & Databases',
            '<strong>AMRFinderPlus</strong> by <a href="https://github.com/ncbi/amr" target="_blank" style="color:#dc3545;font-weight:bold;">NCBI</a> → '
            'resistance gene and point-mutation detection (Feldgarden et al., <em>Sci Rep</em> 2021).<br>'
            '<strong>ABRicate</strong> by <a href="https://github.com/tseemann/abricate" target="_blank" style="color:#dc3545;font-weight:bold;">Prof. Torsten Seemann</a> → '
            'powered by databases: '
            '<a href="https://card.mcmaster.ca/" target="_blank" style="color:#28a745;font-weight:600;">CARD 2023</a>, '
            '<a href="https://www.mediterranee-infection.com/amr-databases/" target="_blank" style="color:#17a2b8;font-weight:600;">ARG-ANNOT</a>, '
            '<a href="https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/" target="_blank" style="color:#007bff;font-weight:600;">NCBI AMR</a>, '
            '<a href="https://megares.meglab.org/" target="_blank" style="color:#fd7e14;font-weight:600;">MEGARes 3.0</a>, '
            '<a href="https://genepi.food.dtu.dk/resfinder" target="_blank" style="color:#20c997;font-weight:600;">ResFinder</a>.<br>'
            '<span style="font-size:.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We are deeply grateful to all developers and curators for their open-source contributions.</span>')

        why_multi = self._collapsible(
            '🔬 Why Do We Use Multiple Databases?', '#00695c',
            '<p style="margin:6px 0;font-size:.95em;line-height:1.5;">EcoliTyper screens every genome against '
            '<strong>five independent AMR databases</strong> plus <strong>AMRFinderPlus</strong> — because '
            '<strong>no single database is comprehensive</strong>. Each has unique strengths, biases, and update cadences.</p>'
            '<p style="margin:10px 0 6px 0;font-size:.95em;"><strong>⚠️ The problem with choosing one database:</strong></p>'
            '<ul style="margin:6px 0 10px 20px;font-size:.93em;line-height:1.6;">'
            '<li>Some researchers pick a favourite database (often CARD or ResFinder) and only report hits from that one.</li>'
            '<li>This is <strong>fast but incomplete</strong> — a gene absent from CARD may still be present in MEGARes or AMRFinderPlus.</li>'
            '<li>Conversely, single-database hits with weak support can be <strong>false positives</strong> that a second database would have flagged.</li>'
            '<li>The result: <em>biased prevalence estimates</em> and <em>missed resistance signals</em> that matter clinically.</li>'
            '</ul>'
            '<p style="margin:10px 0 6px 0;font-size:.95em;"><strong>✅ Our approach — report everything, provenance preserved:</strong></p>'
            '<ul style="margin:6px 0 10px 20px;font-size:.93em;line-height:1.6;">'
            '<li>All hits from all databases are kept <strong>separate and unmerged</strong> — the <strong>Database column</strong> tells you exactly which source found each gene.</li>'
            '<li>You get the <strong>full picture</strong> and can filter/interpret yourself — no silent filtering, no cherry-picking.</li>'
            '<li>Cross-database agreement becomes a <strong>confidence signal</strong>.</li>'
            '</ul>'
            '<p style="margin:6px 0 0 0;font-size:.92em;background:#fff3cd;padding:8px 14px;border-radius:4px;border-left:3px solid #ffc107;">'
            '<i class="fas fa-lightbulb" style="color:#856404;"></i> <strong>Pro tip:</strong> Use the search box to filter by database name '
            '(type <code>CARD</code> or <code>ResFinder</code>). Group by typing to see which lineages carry which resistance signatures.</p>')
        confidence = self._collapsible(
            '🎯 Cross-Database Confidence Tiers', '#0891b2',
            '<p style="margin:6px 0;font-size:.93em;line-height:1.5;">A gene detected by <strong>multiple databases</strong> is far more likely to be a true positive than a single-database hit.</p>'
            '<ul style="margin:6px 0 0 20px;font-size:.93em;line-height:1.6;">'
            '<li><span style="color:#16a34a;font-weight:bold;">🟢 High confidence</span> — gene found in <strong>3 or more databases</strong>.</li>'
            '<li><span style="color:#f59e0b;font-weight:bold;">🟡 Moderate confidence</span> — gene found in <strong>2 databases</strong>.</li>'
            '<li><span style="color:#dc2626;font-weight:bold;">🔴 Low confidence / investigate</span> — gene found in <strong>only 1 database</strong>.</li>'
            '</ul>')
        acquired = self._collapsible(
            '🧬 Acquired vs Intrinsic Resistance — Both Matter', '#6f42c1',
            '<ul style="margin:6px 0 10px 20px;font-size:.93em;line-height:1.6;">'
            '<li><strong style="color:#7c3aed;">Intrinsic genes</strong> — part of the species baseline. Present in almost all <em>E. coli</em>. '
            'Examples: <em>acrAB-tolC</em>, <em>ampC</em>, chromosomal efflux pumps. They tell you what the bacterium is naturally tolerant to.</li>'
            '<li><strong style="color:#e11d48;">Acquired genes</strong> — gained by horizontal transfer (plasmids, transposons, integrons). '
            'Examples: <em>blaCTX-M</em>, <em>mcr-1</em>, <em>tet</em>, <em>sul</em>. They tell you what selective pressure the population has been under.</li>'
            '</ul>'
            '<p style="margin:6px 0 0 0;font-size:.92em;"><i class="fas fa-lightbulb" style="color:#6f42c1;"></i> '
            '<strong>Why both:</strong> A population with only intrinsic efflux genes still responds to standard therapy; '
            'a population with <em>blaCTX-M + sul1 + tetA</em> is multidrug-resistant and needs alternative agents.</p>')
        caveat = self._collapsible(
            '⚠️ Gene Presence ≠ Phenotypic Resistance', '#f59e0b',
            '<ul style="margin:6px 0 12px 20px;font-size:.92em;line-height:1.6;">'
            '<li><strong>Silent / truncated genes</strong> — a partial <em>blaCTX-M</em> fragment produces no enzyme.</li>'
            '<li><strong>Expression regulation</strong> — <em>ampC</em> is inducible; a mutation in the promoter can silence it.</li>'
            '<li><strong>Mechanism matters</strong> — efflux (<em>tetA</em>) and ribosomal protection (<em>tetM</em>) produce different phenotypes.</li>'
            '<li><strong>Naming ambiguity</strong> — ResFinder\'s <code>blaCTX-M-15_1</code> = CARD\'s <code>blaCTX-M-15</code>.</li>'
            '<li><strong>Dose and route</strong> — low-level efflux may be overcome in vivo.</li>'
            '</ul>'
            '<p style="margin:8px 0 0 0;font-size:.92em;background:#fff8e1;padding:8px 12px;border-radius:4px;border-left:3px solid #f59e0b;">'
            '<strong>Clinical bottom line:</strong> <strong>AST remains the gold standard</strong> for treatment decisions. '
            'Genomic prediction is a powerful triage tool and surveillance instrument — not a replacement for culture.</p>')
        gene_info = self._collapsible(
            '🧪 Role of Each AMR Gene Family', '#F44336',
            '<div style="font-size:.9em;line-height:1.8;">'
            '• <strong>blaCTX-M / blaTEM / blaSHV</strong> – ESBLs; hydrolyse extended-spectrum cephalosporins<br>'
            '• <strong>blaKPC / blaNDM / blaOXA / blaVIM / blaIMP</strong> – Carbapenemases; destroy last-line carbapenems<br>'
            '• <strong>mcr</strong> – Plasmid-borne colistin resistance; threatens last-resort drug<br>'
            '• <strong>qnr</strong> – Plasmid-mediated quinolone resistance (PMQR)<br>'
            '• <strong>aac / aph / ant / armA / rmt</strong> – Aminoglycoside resistance mechanisms<br>'
            '• <strong>tet</strong> – Tetracycline efflux or ribosomal protection<br>'
            '• <strong>sul / dfr</strong> – Sulfonamide / trimethoprim resistance<br>'
            '• <strong>cat / floR</strong> – Chloramphenicol / florfenicol resistance<br>'
            '</div>')
        filter_buttons = '''
        <div class="action-buttons">
          <button class="action-btn btn-danger" onclick="document.getElementById('search-amr-table').value='blaCTX-M'; searchTable('amr-table','search-amr-table')">ESBL (CTX-M)</button>
          <button class="action-btn btn-danger" onclick="document.getElementById('search-amr-table').value='blaTEM'; searchTable('amr-table','search-amr-table')">ESBL (TEM)</button>
          <button class="action-btn btn-danger" onclick="document.getElementById('search-amr-table').value='blaSHV'; searchTable('amr-table','search-amr-table')">ESBL (SHV)</button>
          <button class="action-btn btn-danger" onclick="document.getElementById('search-amr-table').value='blaNDM'; searchTable('amr-table','search-amr-table')">Carbapenemase (NDM)</button>
          <button class="action-btn btn-danger" onclick="document.getElementById('search-amr-table').value='blaKPC'; searchTable('amr-table','search-amr-table')">Carbapenemase (KPC)</button>
          <button class="action-btn btn-danger" onclick="document.getElementById('search-amr-table').value='blaOXA'; searchTable('amr-table','search-amr-table')">Carbapenemase (OXA)</button>
          <button class="action-btn btn-warning" onclick="document.getElementById('search-amr-table').value='mcr'; searchTable('amr-table','search-amr-table')">Colistin (mcr)</button>
          <button class="action-btn btn-info" onclick="document.getElementById('search-amr-table').value='qnr'; searchTable('amr-table','search-amr-table')">Quinolone (qnr)</button>
          <button class="action-btn btn-info" onclick="document.getElementById('search-amr-table').value='aac'; searchTable('amr-table','search-amr-table')">Aminoglycoside (aac)</button>
          <button class="action-btn btn-info" onclick="document.getElementById('search-amr-table').value='armA'; searchTable('amr-table','search-amr-table')">High-level AG (armA)</button>
          <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr-table').value='tet'; searchTable('amr-table','search-amr-table')">Tetracycline (tet)</button>
          <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr-table').value='sul'; searchTable('amr-table','search-amr-table')">Sulfonamide (sul)</button>
          <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr-table').value='dfr'; searchTable('amr-table','search-amr-table')">Trimethoprim (dfr)</button>
          <button class="action-btn btn-light" onclick="document.getElementById('search-amr-table').value=''; searchTable('amr-table','search-amr-table')">Clear</button>
        </div>'''
        return f'''
        {credit}
        {self._alert('info', 'fa-biohazard',
            '<h3>🧬 AMR Genes — Gene-Centric View</h3>'
            '<p>Each AMR gene is shown with <strong>all genomes that carry it</strong>. Use the filter buttons to focus on '
            'important resistance classes, the genome search box to highlight specific isolates, and the grouping dropdown '
            'to see which clonal backgrounds carry each gene. Critical resistance genes are flagged with ⚠️.</p>')}
        {filter_buttons}
        {self._gene_table('amr-table', 'All AMR Genes Across Databases', all_genes, total, self.data_analyzer.critical_amr_genes)}
        {self._database_cards(amr_dbs)}
        <h3 style="margin-top:30px;"><i class="fas fa-book-medical"></i> Interpretation Guide</h3>
        {why_multi}
        {confidence}
        {acquired}
        {caveat}
        {gene_info}'''

    def _sec_virulence(self, d):
        gc = d['gene_centric']
        vir_dbs = gc.get('virulence_databases', {})
        all_genes = []
        for genes in vir_dbs.values():
            all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        total = len(d['samples'])

        credit = self._credit_bar('#8BC34A', '🦠', 'Virulence Detection Tools & Databases',
            '<strong>ABRicate</strong> by <a href="https://github.com/tseemann/abricate" target="_blank" style="color:#8BC34A;font-weight:bold;">Prof. Torsten Seemann</a> → '
            'powered by the <a href="http://www.mgc.ac.cn/VFs/" target="_blank" style="color:#8BC34A;font-weight:bold;">VFDB 2025 (Virulence Factor Database)</a> '
            '(Zhou et al., <em>Nucleic Acids Res</em> 2025).<br>'
            '<span style="font-size:.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We thank the VFDB curators and Prof. Seemann for their invaluable contributions.</span>')

        why_vir = self._collapsible(
            '🔬 Virulence — Why It Matters Clinically', '#00695c',
            '<p>Virulence gene content distinguishes harmless commensals from pathogens. Combined with typing, you can answer:</p>'
            '<ul style="margin-left:20px;">'
            '<li>Which STs carry <em>stx1/stx2</em>? (STEC outbreak risk)</li>'
            '<li>Do <em>eae</em>-positive strains cluster in specific serotypes? (EPEC/EHEC)</li>'
            '<li>Which lineages carry <em>cnf1</em>, <em>hlyA</em>, <em>iroN</em>? (UPEC/ExPEC)</li>'
            '<li>Do <em>aggR</em>+ strains belong to ST10, ST31, or ST40? (EAEC)</li>'
            '</ul>')
        gene_info = self._collapsible(
            '🧪 Role of Each Virulence Family', '#8BC34A',
            '<div style="font-size:.9em;line-height:1.8;">'
            '• <strong>stx1/stx2</strong> – Shiga toxins (STEC/EHEC)<br>'
            '• <strong>eae / tir / paa</strong> – Intimin and T3SS effectors (EPEC/EHEC)<br>'
            '• <strong>bfpA</strong> – Bundle-forming pilus (EPEC)<br>'
            '• <strong>aggR / aatA</strong> – EAEC master regulator + dispersin<br>'
            '• <strong>lt / sta1</strong> – Heat-labile/stable enterotoxins (ETEC)<br>'
            '• <strong>ipaH</strong> – Invasion plasmid antigen (EIEC)<br>'
            '• <strong>cnf1/2</strong> – Cytotoxic necrotising factors (ExPEC)<br>'
            '• <strong>hlyA</strong> – α-hemolysin (UPEC/ExPEC)<br>'
            '• <strong>kps / neu</strong> – Capsule (ExPEC, NMEC)<br>'
            '• <strong>iuc / iro / fyuA</strong> – Iron acquisition (siderophores)<br>'
            '• <strong>clb</strong> – Colibactin genotoxin (colorectal cancer association)<br>'
            '</div>')
        filter_buttons = '''
        <div class="action-buttons">
          <button class="action-btn btn-danger" onclick="document.getElementById('search-vir-table').value='stx'; searchTable('vir-table','search-vir-table')">Shiga toxins (stx)</button>
          <button class="action-btn btn-danger" onclick="document.getElementById('search-vir-table').value='eae'; searchTable('vir-table','search-vir-table')">Intimin (eae)</button>
          <button class="action-btn btn-danger" onclick="document.getElementById('search-vir-table').value='bfp'; searchTable('vir-table','search-vir-table')">Bundle-forming pilus (bfp)</button>
          <button class="action-btn btn-warning" onclick="document.getElementById('search-vir-table').value='aggR'; searchTable('vir-table','search-vir-table')">EAEC (aggR)</button>
          <button class="action-btn btn-warning" onclick="document.getElementById('search-vir-table').value='aat'; searchTable('vir-table','search-vir-table')">EAEC (aat)</button>
          <button class="action-btn btn-warning" onclick="document.getElementById('search-vir-table').value='ipaH'; searchTable('vir-table','search-vir-table')">EIEC (ipaH)</button>
          <button class="action-btn btn-info" onclick="document.getElementById('search-vir-table').value='cnf'; searchTable('vir-table','search-vir-table')">CNF toxins</button>
          <button class="action-btn btn-info" onclick="document.getElementById('search-vir-table').value='hly'; searchTable('vir-table','search-vir-table')">Hemolysins</button>
          <button class="action-btn btn-info" onclick="document.getElementById('search-vir-table').value='iuc'; searchTable('vir-table','search-vir-table')">Aerobactin (iuc)</button>
          <button class="action-btn btn-info" onclick="document.getElementById('search-vir-table').value='iro'; searchTable('vir-table','search-vir-table')">Salmochelin (iro)</button>
          <button class="action-btn btn-secondary" onclick="document.getElementById('search-vir-table').value='clb'; searchTable('vir-table','search-vir-table')">Colibactin</button>
          <button class="action-btn btn-secondary" onclick="document.getElementById('search-vir-table').value='kps'; searchTable('vir-table','search-vir-table')">Capsule (kps)</button>
          <button class="action-btn btn-secondary" onclick="document.getElementById('search-vir-table').value='fim'; searchTable('vir-table','search-vir-table')">Fimbriae (fim)</button>
          <button class="action-btn btn-light" onclick="document.getElementById('search-vir-table').value=''; searchTable('vir-table','search-vir-table')">Clear</button>
        </div>'''
        return f'''
        {credit}
        {self._alert('info', 'fa-virus',
            '<h3>🧬 Virulence Factors — Gene-Centric View</h3>'
            '<p>Each virulence gene is shown with <strong>all genomes that carry it</strong>. Use filter buttons to focus on key '
            'pathotype determinants; the genome search box highlights specific isolates. Critical virulence genes are flagged ⚠️.</p>')}
        {filter_buttons}
        {self._gene_table('vir-table', 'All Virulence Genes', all_genes, total, self.data_analyzer.critical_virulence_genes)}
        {self._database_cards(vir_dbs)}
        <h3 style="margin-top:30px;"><i class="fas fa-book-medical"></i> Interpretation Guide</h3>
        {why_vir}
        {gene_info}'''

    def _sec_plasmid(self, d):
        gc = d['gene_centric']
        pla_dbs = gc.get('plasmid_databases', {})
        all_genes = []
        for genes in pla_dbs.values():
            all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        total = len(d['samples'])
        if not all_genes:
            return self._alert('warning', 'fa-plug', '<h3>No Plasmid Data</h3>')
        credit = self._credit_bar('#3F51B5', '🧬', 'Plasmid & Serotyping Markers – Acknowledgments',
            'Powered by <a href="https://cge.food.dtu.dk/services/PlasmidFinder/" target="_blank" style="color:#3F51B5;font-weight:bold;">PlasmidFinder</a> '
            '(Center for Genomic Epidemiology) and <a href="https://github.com/djingle/EcoH" target="_blank" style="color:#3F51B5;font-weight:bold;">EcoH</a> (Ingle et al.), '
            'both screened via <a href="https://github.com/tseemann/abricate" target="_blank" style="color:#3F51B5;font-weight:bold;">ABRicate</a>.<br>'
            'Please cite: Carattoli A, et al. <em>In silico</em> detection and typing of plasmids using PlasmidFinder. '
            '<em>Antimicrob Agents Chemother</em>. 2014;58(7):3895-903. '
            '<a href="https://doi.org/10.1128/AAC.02412-14" target="_blank" style="color:#3F51B5;font-weight:bold;">doi:10.1128/AAC.02412-14</a>')
        guide = self._collapsible(
            '🧬 Plasmid Replicons & Serotyping Markers', '#3F51B5',
            '<ul>'
            '<li><strong>PlasmidFinder:</strong> Identifies plasmid replicon markers (e.g. IncF, IncI, Col), which are vehicles for '
            'horizontal gene transfer of AMR and virulence genes.</li>'
            '<li><strong>EcoH:</strong> Predicts O and H antigens for classical serotype assignment.</li>'
            '</ul>')
        return f'''
        {credit}
        {self._alert('info', 'fa-plug',
            '<h3>🧬 Plasmid Replicons &amp; Serotyping Markers</h3>'
            '<p>This table combines two analyses: <strong>PlasmidFinder</strong> (replicon markers) and '
            '<strong>EcoH</strong> (O and H antigen prediction).</p>')}
        {self._gene_table('plasmid-table', 'All Plasmid / Serotyping Markers', all_genes, total)}
        {self._database_cards(pla_dbs)}
        <h3 style="margin-top:30px;"><i class="fas fa-book-medical"></i> Interpretation Guide</h3>
        {guide}'''

    def _sec_bacmet(self, d):
        gc = d['gene_centric']
        bac_dbs = gc.get('bacmet_databases', {})
        all_genes = []
        for genes in bac_dbs.values():
            all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        total = len(d['samples'])
        if not all_genes:
            return self._alert('warning', 'fa-flask', '<h3>No BACMET Data</h3>')
        credit = self._credit_bar('#607D8B', '🧪', 'BACMET – Biocide &amp; Heavy Metal Resistance',
            'Powered by the <a href="https://bacmet.biomedicine.gu.se/" target="_blank" style="color:#607D8B;font-weight:bold;">BacMet database</a> '
            'and screened via <a href="https://github.com/tseemann/abricate" target="_blank" style="color:#607D8B;font-weight:bold;">ABRicate</a>.<br>'
            'Please cite: Pal C, et al. BacMet: antibacterial biocide and metal resistance genes database. <em>Nucleic Acids Res</em>. 2014;42(Database issue):D737-43. '
            '<a href="https://doi.org/10.1093/nar/gkt1252" target="_blank" style="color:#607D8B;font-weight:bold;">doi:10.1093/nar/gkt1252</a>')
        guide = self._collapsible(
            '🧪 BACMET — Biocide & Heavy-Metal Resistance', '#607D8B',
            '<p>Environmental co-selection markers. When a plasmid carries both an antibiotic resistance gene and a '
            'metal/biocide resistance gene, exposure to the metal or biocide can select for antibiotic resistance '
            'without any antibiotic pressure — a major One Health concern.</p>'
            '<ul>'
            '<li><strong>qac / cep / form</strong> — hospital disinfectant resistance (QACs, chlorhexidine, formaldehyde).</li>'
            '<li><strong>mer / ars / pbr / chr</strong> — mercury, arsenic, lead, chromate resistance.</li>'
            '<li><strong>sil / cop / pco</strong> — silver and copper resistance (wound dressings, feed additives).</li>'
            '<li><strong>czc / cad / znt / nik</strong> — cadmium, zinc, cobalt, nickel efflux.</li>'
            '<li><strong>soxR / cpxR / baeR</strong> — stress regulators that also upregulate multidrug efflux pumps.</li>'
            '</ul>')
        return f'''
        {credit}
        {self._alert('info', 'fa-flask',
            '<h3>🧪 BACMET — Biocide &amp; Heavy-Metal Resistance</h3>'
            '<p>Environmental co-selection markers. Use the grouping dropdown to see which clones carry them.</p>')}
        {self._gene_table('bac-table', 'All BACMET Genes', all_genes, total)}
        {self._database_cards(bac_dbs)}
        <h3 style="margin-top:30px;"><i class="fas fa-book-medical"></i> Interpretation Guide</h3>
        {guide}'''

    def _sec_mutation(self, d):
        mut = d.get('mutation_data', {})
        mutations = mut.get('mutations', [])
        if not mutations:
            return self._alert('warning', 'fa-dna', '<h3>No Mutation Data</h3><p>mutation_summary.html not found.</p>')
        total = len(d['samples'])
        credit = self._credit_bar('#00BCD4', '🧬', 'Point Mutations – AMRFinderPlus',
            'Point mutations detected by <a href="https://github.com/ncbi/amr" target="_blank" style="color:#00BCD4;font-weight:bold;">NCBI AMRFinderPlus</a>.<br>'
            'Please cite: Feldgarden M, et al. AMRFinderPlus and the Reference Gene Catalog. <em>Sci Rep</em>. 2021;11(1):12728. '
            '<a href="https://doi.org/10.1038/s41598-021-91456-0" target="_blank" style="color:#00BCD4;font-weight:bold;">doi:10.1038/s41598-021-91456-0</a>')
        rows = ''
        for m in mutations:
            pct = (m['count'] / total * 100) if total else 0
            tags = ''.join(f'<span class="genome-tag">{g}</span>' for g in m['genomes'])
            rows += f'''<tr>
                <td class="mut-gene"><strong>{m["gene"]}</strong></td>
                <td class="mut-class">{m["class"]}</td>
                <td class="mut-count"><strong>{m["count"]}</strong> <span style="color:#666;">({pct:.1f}%)</span></td>
                <td class="mut-genomes"><div class="genome-list">{tags}</div></td>
            </tr>'''
        return f'''
        {credit}
        {self._alert('info', 'fa-dna',
            '<h3>🧬 Point Mutations — Gene-Centric View with Grouping</h3>'
            '<p>Chromosomal mutations that confer resistance when no acquired gene is present. '
            'Key examples: <em>gyrA</em>/<em>parC</em> (quinolones), <em>rpoB</em> (rifampin), '
            '<em>pmrAB</em>/<em>arnT</em> (colistin), <em>ampC</em> promoter (cephalosporins).</p>'
            '<p>Each unique mutation is shown with <strong>all genomes</strong> that carry it.</p>')}
        {self._grouping_dropdown('mutation-table')}
        <input type="text" class="search-box" id="search-mut"
               onkeyup="searchTable('mutation-table','search-mut')" placeholder="🔍 Search mutations...">
        <input type="text" class="search-box" id="highlight-mutation-table"
               onkeyup="highlightGenome('mutation-table','highlight-mutation-table')"
               placeholder="🔍 Highlight genome tags">
        <div class="action-buttons">
          <button class="action-btn btn-primary" onclick="exportTableToCSV('mutation-table','mutations.csv')">
            <i class="fas fa-download"></i> Export CSV</button>
        </div>
        <style>
          #mutation-table {{ table-layout: fixed; width: 100%; }}
          #mutation-table th.mut-gene, #mutation-table td.mut-gene {{ width: 22%; word-break: break-word; }}
          #mutation-table th.mut-class, #mutation-table td.mut-class {{ width: 20%; word-break: break-word; }}
          #mutation-table th.mut-count, #mutation-table td.mut-count {{ width: 12%; text-align: center; white-space: nowrap; }}
          #mutation-table th.mut-genomes, #mutation-table td.mut-genomes {{ width: 46%; }}
          #mutation-table td.mut-genomes .genome-list {{ max-height: 220px; }}
        </style>
        <div class="master-scrollable-container">
          <table id="mutation-table" class="data-table">
            <thead><tr>
              <th class="mut-gene" data-sort="string">Gene</th>
              <th class="mut-class" data-sort="string">Class</th>
              <th class="mut-count" data-sort="number">Count</th>
              <th class="mut-genomes" data-sort="string">Genomes</th>
            </tr></thead>
            <tbody>{rows}</tbody>
          </table>
        </div>'''

    def _mge_category_cards(self) -> str:
        cards = [
            ('Integr.', 'Integrase', '#ef4444', 'fa-cut',
             'Enzymes that catalyse integration and excision of mobile elements.',
             'Active chromosomal integration of foreign DNA.'),
            ('Transfer', 'Transfer / Conjugation', '#f97316', 'fa-share-alt',
             'Type IV secretion systems (T4SS), relaxases, coupling proteins — the machinery of horizontal gene transfer.',
             'Primary route by which AMR spreads.'),
            ('Stab.', 'Stability / Maintenance', '#eab308', 'fa-anchor',
             'Partitioning systems (parAB), toxin-antitoxin cassettes, addiction modules.',
             'A plasmid with strong stability persists even without antibiotic selection.'),
            ('Phage', 'Bacteriophage', '#22c55e', 'fa-virus',
             'Complete phage structural and lysogeny proteins.',
             'Phages carry virulence genes (e.g. stx) and mobilise chromosomal DNA.'),
            ('Repl.', 'Replication', '#14b8a6', 'fa-sync',
             'Replication initiators, primases, helicases.',
             'Only elements with functional replication persist as extrachromosomal DNA.'),
            ('IS-assoc.', 'Insertion Sequences', '#0891b2', 'fa-exchange-alt',
             'Transposases and IS elements.',
             'Frequently disrupt genes or move resistance genes into new contexts.'),
            ('ICE-assoc.', 'Integrative Conjugative Elements', '#3b82f6', 'fa-project-diagram',
             'Self-transmissible elements that integrate into the chromosome.',
             'Move large multi-gene cassettes between unrelated bacteria.'),
            ('Plasmid-assoc.', 'Plasmid', '#8b5cf6', 'fa-plug',
             'Plasmid backbone, replication, partitioning, conjugal transfer proteins.',
             'Single most important vehicle for acquired antibiotic resistance.'),
            ('Phage-assoc.', 'Phage-associated (lysogeny)', '#d946ef', 'fa-dna',
             'Prophage and lysogeny-associated proteins distinct from complete phage particles.',
             'Prophage carriage correlates with virulence gene acquisition.'),
            ('Key MGE-assoc.', 'Key MGE signatures', '#dc2626', 'fa-star',
             'Curated high-confidence mobileOG entries.',
             'Cleanest single indicator of MGE activity in a genome.'),
        ]
        html = '<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:15px;">'
        for key, label, color, icon, definition, importance in cards:
            html += f'''
            <div style="background:#fff;padding:16px;border-radius:10px;border-left:4px solid {color};box-shadow:0 2px 8px rgba(0,0,0,.06);">
                <div style="display:flex;align-items:center;gap:10px;margin-bottom:8px;">
                    <i class="fas {icon}" style="color:{color};font-size:1.4em;"></i>
                    <div><strong style="color:{color};font-size:1.05em;">{label}</strong>
                    <div style="font-size:.8em;color:#888;font-family:monospace;">{key}</div></div>
                </div>
                <p style="font-size:.88em;color:#333;margin:0 0 8px 0;">{definition}</p>
                <p style="font-size:.85em;color:#555;margin:0;padding-top:8px;border-top:1px dashed #e0e0e0;">
                    <i class="fas fa-bullseye" style="color:{color};"></i> <em>{importance}</em>
                </p>
            </div>'''
        html += '</div>'
        return html

    def _build_mge_aggregate(self, per_sample):
        samples_data = getattr(self, '_current_samples_data', {})
        if not samples_data or not per_sample:
            return ''
        fields = [('ST', 'MLST'), ('Serotype', 'Serotype'),
                  ('Phylogroup', 'Phylogroup'), ('CH_Type', 'CH type')]
        sample_cols = list(next(iter(per_sample.values())).keys())
        blocks = ''
        for field_key, field_label in fields:
            groups = defaultdict(list)
            for sample, vals in per_sample.items():
                sd = samples_data.get(sample, {})
                if field_key == 'ST':
                    typ = sd.get('mlst', {}).get('ST', 'ND')
                elif field_key == 'Serotype':
                    typ = sd.get('serotype', {}).get('Serotype', 'ND')
                elif field_key == 'Phylogroup':
                    typ = sd.get('phylogrouping', {}).get('Clermont_Type', 'ND')
                else:
                    typ = sd.get('chtyper', {}).get('CH_Type', 'ND')
                if typ in ('ND', ''):
                    continue
                groups[typ].append(vals)
            if not groups:
                continue
            rows = ''
            for grp, records in sorted(groups.items(), key=lambda x: -len(x[1])):
                means = {}
                for c in sample_cols:
                    vals_ = [r.get(c, 0) for r in records if isinstance(r.get(c, 0), (int, float))]
                    means[c] = round(sum(vals_) / len(vals_), 1) if vals_ else 0
                cells = ''.join(f'<td>{means[c]}</td>' for c in sample_cols)
                rows += f'<tr><td><strong>{grp}</strong></td><td>{len(records)}</td>{cells}</tr>'
            header = ''.join(f'<th>{c}</th>' for c in sample_cols)
            blocks += f'''
            <h4 style="margin-top:20px;">Mean MGE profile by {field_label}</h4>
            <div class="master-scrollable-container">
              <table class="data-table">
                <thead><tr>
                  <th data-sort="string">{field_label}</th>
                  <th data-sort="number">N samples</th>
                  {header}
                </tr></thead>
                <tbody>{rows}</tbody>
              </table>
            </div>'''
        return blocks

    def _sec_mge(self, d):
        mge = d.get('mge_data', {})
        stats = mge.get('stats', {})
        per_sample = mge.get('per_sample', {})
        analysis = mge.get('analysis', {})
        note = mge.get('note', '')
        credit = self._credit_bar('#16A085', '🧬', 'Mobile Genetic Elements – mobileOG-db',
            '<strong>mobileOG-db</strong> developed by '
            '<a href="https://github.com/clb21565/mobileOG-db" target="_blank" style="color:#16A085;font-weight:bold;">Brown, Mullet and colleagues</a> — '
            'a manually curated database of protein families mediating the life cycle of bacterial mobile genetic elements.<br>'
            'Please cite: Brown CL, et al. mobileOG-db: a Manually Curated Database of Protein Families Mediating the Life Cycle of Bacterial Mobile Genetic Elements. '
            '<em>Appl Environ Microbiol</em>. 2022;88(18):e00991-22. '
            '<a href="https://doi.org/10.1128/aem.00991-22" target="_blank" style="color:#16A085;font-weight:bold;">doi:10.1128/aem.00991-22</a>')
        why_mge = self._collapsible(
            '🧬 Why Mobile Genetic Elements Matter for AMR', '#16A085',
            '<p>Mobile genetic elements (MGEs) are the <strong>physical vehicles of horizontal gene transfer</strong> — '
            'the mechanism by which bacteria acquire resistance and virulence genes from other bacteria.</p>'
            '<ul style="margin-top:8px;">'
            '<li><strong>Plasmids</strong> carry <em>blaCTX-M</em>, <em>mcr-1</em>, <em>tet</em>, <em>sul</em> between lineages.</li>'
            '<li><strong>Integrons</strong> capture and express cassettes of resistance genes.</li>'
            '<li><strong>Transposons and IS elements</strong> move genes within and between replicons.</li>'
            '<li><strong>Prophages</strong> carry Shiga toxin genes (<em>stx1/stx2</em>) and other virulence modules.</li>'
            '</ul>'
            '<p style="margin-top:8px;"><strong>Interpreting this tab:</strong> higher mobileOG hit counts suggest a genome that has '
            'accumulated more mobile-element machinery — potentially more permissive for acquiring new resistance genes.</p>')
        category_cards = self._collapsible(
            '📚 Understanding the Categories', '#16A085', self._mge_category_cards())
        if not per_sample:
            return f'''{credit}{why_mge}{category_cards}
            {self._alert('warning', 'fa-mobile-alt', '<h3>No MGE Data Available</h3>')}'''
        cards = ''
        cards += self._stat_card(stats.get('Total Samples', len(per_sample)), 'Total Samples', '#4CAF50', 'fa-vial')
        cards += self._stat_card(stats.get('With mobileOG hits', '—'), 'With mobileOG hits', '#2196F3', 'fa-check')
        cards += self._stat_card(stats.get('Without mobileOG hits', '—'), 'Without mobileOG hits', '#9E9E9E', 'fa-times')
        cards += self._stat_card(stats.get('Total mobileOG hits', analysis.get('grand_total', '—')), 'Total mobileOG hits', '#FF9800', 'fa-database')
        cards += self._stat_card(stats.get('Key MGE-associated genes', '—'), 'Key MGE-associated genes', '#8B5CF6', 'fa-star')
        cards += self._stat_card(stats.get('Runtime', '—'), 'Runtime', '#E53935', 'fa-clock')
        cat_rows = ''
        for c in analysis.get('category_totals', []):
            cat_rows += (f'<tr><td><strong>{c["category"]}</strong></td>'
                         f'<td>{c["total"]:,}</td><td>{c["pct_of_all"]}%</td>'
                         f'<td>{c["mean"]}</td><td>{c["max"]}</td>'
                         f'<td>{c["max_sample"]}</td></tr>')
        all_cols = ['Sample']
        seen = set()
        for s in per_sample.values():
            for k in s.keys():
                if k not in seen:
                    seen.add(k); all_cols.append(k)
        header = ''.join(f'<th data-sort="{"string" if c=="Sample" else "number"}">{c}</th>' for c in all_cols)
        sample_rows = ''
        for s, vals in sorted(per_sample.items()):
            cells = f'<td><strong>{s}</strong></td>'
            for c in all_cols[1:]:
                cells += f'<td>{vals.get(c, 0)}</td>'
            sample_rows += f'<tr>{cells}</tr>'
        agg_html = self._build_mge_aggregate(per_sample)
        note_html = ''
        if note:
            note_html = (f'<div class="alert-box alert-warning" style="font-size:.9em;">'
                         f'<strong>Interpretation note:</strong> {note}</div>')
        return f'''
        {credit}
        {self._alert('info', 'fa-mobile-alt',
            '<h3>🧬 Mobile Genetic Elements — mobileOG Profile</h3>'
            '<p>Per-sample counts of mobileOG-db protein family hits across multiple functional categories. '
            'Individual matches represent MGE-associated protein signatures and should not be interpreted as '
            'independent mobile genetic elements — multiple proteins from the same element may generate multiple hits, '
            'and some families (particularly replication, recombination, and repair proteins) also occur in the bacterial chromosome. '
            'Stronger element-level inference requires genomic context and co-localization of compatible MGE-associated genes.</p>')}
        <div class="stats-grid">{cards}</div>
        {note_html}

        <h3 style="margin-top:30px;">📊 Category Totals (aggregated across all samples)</h3>
        <div class="scrollable-table"><table class="data-table">
          <thead><tr><th>Category</th><th>Total hits</th><th>% of all hits</th>
                     <th>Mean/sample</th><th>Max/sample</th><th>Sample with max</th></tr></thead>
          <tbody>{cat_rows}</tbody></table></div>

        <h3 style="margin-top:30px;">🧬 Per-sample mobileOG Profile</h3>
        <input type="text" class="search-box" id="search-mge"
               onkeyup="searchTable('mge-table','search-mge')" placeholder="🔍 Search sample...">
        <div class="action-buttons">
          <button class="action-btn btn-primary" onclick="exportTableToCSV('mge-table','mge_profile.csv')">
            <i class="fas fa-download"></i> Export MGE CSV</button>
        </div>
        <div class="master-scrollable-container">
          <table id="mge-table" class="data-table">
            <thead><tr>{header}</tr></thead>
            <tbody>{sample_rows}</tbody>
          </table>
        </div>

        <h3 style="margin-top:30px;">📈 Aggregate MGE Profile by Typing</h3>
        {agg_html}

        <h3 style="margin-top:30px;"><i class="fas fa-book-medical"></i> Interpretation Guide</h3>
        {why_mge}
        {category_cards}'''

    def _precompute_similarity(self, compare_data: Dict[str, Dict]) -> Dict[str, Any]:
        samples = sorted(compare_data.keys())
        n = len(samples)
        if n < 2:
            return {'samples': samples, 'flat': [], 'n': n}
        cats = ('AMR', 'Virulence', 'BACMET', 'Plasmids', 'Mutations')
        gene_sets, typing_data = {}, {}
        for s in samples:
            d = compare_data[s]
            gene_sets[s] = {cat: set(d['genes'].get(cat, [])) for cat in cats}
            typing_data[s] = d['typing']
        typing_keys = list(typing_data[samples[0]].keys()) if samples else []
        n_keys = len(typing_keys)
        flat = []
        for i in range(n):
            s1 = samples[i]
            for j in range(i + 1, n):
                s2 = samples[j]
                matches = sum(1 for k in typing_keys if typing_data[s1].get(k) == typing_data[s2].get(k))
                typing_pct = (matches / n_keys * 100) if n_keys else 0
                cat_scores = []
                for cat in cats:
                    sa, sb = gene_sets[s1][cat], gene_sets[s2][cat]
                    union = sa | sb
                    if not union:
                        continue
                    cat_scores.append(len(sa & sb) / len(union) * 100)
                gene_pct = (sum(cat_scores) / len(cat_scores)) if cat_scores else 0
                overall = round(typing_pct * 0.5 + gene_pct * 0.5, 1)
                flat.append(overall)
        return {'samples': samples, 'flat': flat, 'n': n}

    def _sec_compare(self, d):
        samples_data = d.get('samples', {})
        mutation_data = d.get('mutation_data', {})
        mutations_list = mutation_data.get('mutations', [])
        if len(samples_data) < 2:
            return self._alert('warning', 'fa-balance-scale',
                f'<h3>Not Enough Samples</h3><p>Compare requires ≥2 samples. Dataset has {len(samples_data)}.</p>')
        sample_mutations = defaultdict(list)
        for m in mutations_list:
            label = f"{m.get('gene','')} {m.get('mutation','')}".strip()
            for g in m.get('genomes', []):
                sample_mutations[g].append(label)
        compare_data = {}
        def db_genes(sample, db):
            return list(samples_data.get(sample, {}).get('abricate_databases', {}).get(db, []))
        for s, sd in samples_data.items():
            amr_genes = sorted(set(
                sd.get('amrfinder', {}).get('all_genes', [])
                + db_genes(s, 'card') + db_genes(s, 'resfinder')
                + db_genes(s, 'argannot') + db_genes(s, 'megares')
                + db_genes(s, 'ncbi')))
            vir_genes = sorted(set(db_genes(s, 'vfdb') + db_genes(s, 'ecoli_vf')))
            bac_genes = sorted(set(db_genes(s, 'bacmet2')))
            pla_genes = sorted(set(db_genes(s, 'plasmidfinder') + db_genes(s, 'ecoh')))
            muts = sorted(set(sample_mutations.get(s, [])))
            compare_data[s] = {
                'typing': {
                    'ST': sd.get('mlst', {}).get('ST', 'ND'),
                    'Serotype': sd.get('serotype', {}).get('Serotype', 'ND'),
                    'Phylogroup': sd.get('phylogrouping', {}).get('Clermont_Type', 'ND'),
                    'CH Type': sd.get('chtyper', {}).get('CH_Type', 'ND'),
                },
                'genes': {
                    'AMR': amr_genes, 'Virulence': vir_genes,
                    'BACMET': bac_genes, 'Plasmids': pla_genes,
                    'Mutations': muts,
                },
            }
        sorted_samples = sorted(compare_data.keys())
        options_a = ''.join(f'<option value="{s}"{" selected" if i==0 else ""}>{s}</option>' for i, s in enumerate(sorted_samples))
        options_b = ''.join(f'<option value="{s}"{" selected" if i==1 else ""}>{s}</option>' for i, s in enumerate(sorted_samples))
        blob_json = json.dumps(compare_data, default=str)
        precomp = self._precompute_similarity(compare_data)
        precomp_json = json.dumps(precomp, default=str)
        return f'''
        <div class="alert-box alert-info">
          <i class="fas fa-balance-scale fa-2x"></i>
          <div>
            <h3>⚖️ Compare &amp; Cluster Analysis</h3>
            <p>Two modes: <strong>Pair Compare</strong> for a deep side-by-side of two isolates, and
            <strong>Cluster Detection</strong> to find groups of highly similar samples across the whole dataset.</p>
            <ul style="margin-top:6px;">
              <li><strong>Outbreak confirmation</strong> — do two isolates share ST + Serotype + Phylogroup + CH type?</li>
              <li><strong>Discordance analysis</strong> — what changed between an ancestor and a resistant descendant?</li>
              <li><strong>Cluster detection</strong> — which samples form tight transmission clusters (≥95% similarity)?</li>
            </ul>
          </div>
        </div>
        <div class="mode-toggle">
          <button class="mode-btn active" data-mode="pair" onclick="setCompareMode('pair')">
            <i class="fas fa-balance-scale"></i> Pair Compare</button>
          <button class="mode-btn" data-mode="cluster" onclick="setCompareMode('cluster')">
            <i class="fas fa-project-diagram"></i> Cluster Detection</button>
        </div>
        <div id="mode-pair" class="compare-mode active">
          <div class="compare-panel">
            <div class="compare-select-row">
              <label for="compare-a"><span>Sample A</span>
                <select id="compare-a">{options_a}</select></label>
              <button class="action-btn btn-primary" onclick="runCompare()">
                <i class="fas fa-balance-scale"></i> Compare</button>
              <button class="action-btn btn-light" onclick="swapCompare()">
                <i class="fas fa-exchange-alt"></i> Swap</button>
              <button class="action-btn btn-light" onclick="resetCompare()">
                <i class="fas fa-sync"></i> Reset</button>
              <label for="compare-b"><span>Sample B</span>
                <select id="compare-b">{options_b}</select></label>
            </div>
            <div id="compare-verdict" style="display:none;"></div>
            <div id="compare-visuals"></div>
            <div id="compare-output"></div>
          </div>
        </div>
        <div id="mode-cluster" class="compare-mode">
          <div class="compare-panel">
            <div class="cluster-controls">
              <label><span>Similarity threshold:</span>
                <select id="cluster-threshold" onchange="runCluster()">
                  <option value="95">Tight clusters (≥95%)</option>
                  <option value="90" selected>Moderate clusters (≥90%)</option>
                  <option value="85">Loose clusters (≥85%)</option>
                  <option value="80">Broad clusters (≥80%)</option>
                </select></label>
              <button class="action-btn btn-primary" onclick="runCluster()">
                <i class="fas fa-project-diagram"></i> Detect Clusters</button>
              <span id="cluster-summary" class="sel-counter"></span>
            </div>
            <div id="cluster-clusters"></div>
            <div id="cluster-heatmap"></div>
          </div>
        </div>
        <script>
        window.ECOLITYPER_COMPARE_DATA = {blob_json};
        window.ECOLITYPER_SIM = {precomp_json};
        </script>'''

    def _sec_patterns(self, d):
        P = d['patterns']
        html = self._alert('info', 'fa-project-diagram',
            '<h3>🔍 Cross-Genome Pattern Discovery</h3>'
            '<p>Samples that carry both critical AMR genes (ESBL, carbapenemase, mcr) AND critical virulence genes '
            '(stx, eae, cnf1, etc.) are flagged as high-risk. These combinations pose a serious threat to public health '
            'because they combine difficult-to-treat resistance with pathogenic potential.</p>'
            '<p>The co-occurrence table shows which pairs of genes are frequently found together across your samples.</p>')
        high = P.get('high_risk_combinations', [])
        if high:
            html += '<div class="alert-box alert-danger"><i class="fas fa-exclamation-triangle"></i><div><h3>⚠️ High-Risk Samples Detected</h3><p>These samples carry both critical AMR and critical virulence genes.</p></div></div>'
            html += '<input type="text" class="search-box" id="search-highrisk" onkeyup="searchTable(\'highrisk-table\',\'search-highrisk\')" placeholder="🔍 Search high-risk combinations...">'
            rows = ''
            for c in high:
                rows += (f'<tr><td><strong>{c["sample"]}</strong></td>'
                         f'<td>{c["st"]}</td><td>{c["serotype"]}</td>'
                         f'<td>{c.get("phylogroup","ND")}</td>'
                         f'<td>{", ".join(c["critical_amr_genes"])}</td>'
                         f'<td>{", ".join(c["critical_virulence_genes"])}</td></tr>')
            html += f'''<h3>⚠️ High-Risk Combinations ({len(high)})</h3>
            <div class="master-scrollable-container"><table id="highrisk-table" class="data-table">
              <thead><tr><th>Sample</th><th>ST</th><th>Serotype</th><th>Phylogroup</th>
              <th>Critical AMR</th><th>Critical Virulence</th></tr></thead>
              <tbody>{rows}</tbody></table></div>'''
        else:
            html += '<p>No high-risk combinations (critical AMR + critical virulence) were detected.</p>'
        cooc = P.get('gene_cooccurrence', {})
        if cooc:
            cooc_list = []
            for g1, partners in cooc.items():
                for g2, cnt in partners.items():
                    pair = tuple(sorted([g1, g2]))
                    cooc_list.append((pair[0], pair[1], cnt))
            unique_pairs = {}
            for g1, g2, cnt in cooc_list:
                key = (g1, g2) if g1 <= g2 else (g2, g1)
                unique_pairs[key] = unique_pairs.get(key, 0) + cnt
            cooc_list = [(g1, g2, cnt) for (g1, g2), cnt in unique_pairs.items()]
            cooc_list.sort(key=lambda x: x[2], reverse=True)
            html += '<h3>📈 Gene Co-occurrence (top 500 pairs)</h3>'
            html += '<input type="text" class="search-box" id="search-cooc" onkeyup="searchTable(\'cooc-table\',\'search-cooc\')" placeholder="🔍 Search gene pairs...">'
            rows = ''.join(f'<tr><td>{g1}</td><td>{g2}</td><td>{c}</td></tr>' for g1, g2, c in cooc_list[:500])
            html += f'''<div class="master-scrollable-container"><table id="cooc-table" class="data-table">
              <thead><tr><th>Gene 1</th><th>Gene 2</th><th>Co-occurrence</th></tr></thead>
              <tbody>{rows}</tbody></table></div>'''
        return html

    def _sec_databases(self, d):
        stats = d['gene_centric'].get('database_stats', {})
        if not stats:
            return self._alert('warning', 'fa-database', '<h3>No database metrics available.</h3>')
        rows = ''
        for db, v in stats.items():
            rows += f'<tr><td><strong>{db.upper()}</strong></td><td>{v["total_genes"]}</td><td>{v["total_occurrences"]}</td></tr>'
        return f'''
        {self._alert('info', 'fa-database',
            '<h3>Database Metrics</h3>'
            '<p>Number of unique genes and total occurrences per database. This helps assess the completeness of your '
            'screening and the diversity of resistance/virulence genes detected.</p>')}
        <div class="master-scrollable-container"><table class="data-table">
          <thead><tr><th>Database</th><th>Unique Genes</th><th>Total Occurrences</th></tr></thead>
          <tbody>{rows}</tbody></table></div>'''

    def _sec_aiguide(self):
        return '''
        <div class="alert-box alert-info">
          <i class="fas fa-robot fa-2x"></i>
          <div>
            <h3>🤖 AI Assistant Guide — Your Co-Pilot for Science, Not a Substitute</h3>
            <p>You've generated a rich dataset. Now let AI help you turn numbers into narratives — <strong>ethically, responsibly, and with you in the driver's seat</strong>.</p>
          </div>
        </div>

        <div class="info-text">
          <i class="fas fa-lightbulb"></i> <strong>The Golden Rule of AI in Science:</strong><br>
          AI is like a hyper-caffeinated grad student who never sleeps — great at summarising, spotting patterns, and drafting, but <strong>it can hallucinate, cite fake papers, and misunderstand context</strong>. Always verify, always think critically, and never let AI replace your expert judgment.
        </div>

        <h4>🚀 How to Use This Report with AI (ChatGPT, Claude, Gemini, etc.)</h4>
        <ol style="margin-bottom:20px;margin-left:20px;">
          <li><strong>Export the full data as JSON</strong> (click the <em>Export</em> tab → "Download JSON").</li>
          <li><strong>Upload the JSON file</strong> to your preferred AI chat interface (e.g. ChatGPT Plus with file upload, Claude, or Gemini Advanced).</li>
          <li><strong>Ask specific, contextual questions</strong> — the AI will read your actual data and answer based on it.</li>
          <li><strong>Use the generated insights</strong> to draft sections of your manuscript, thesis, or public health report — but <strong>always fact-check against the original tables</strong>.</li>
        </ol>

        <h4>📝 Using AI to Draft Manuscripts &amp; Theses — A Step-by-Step Workflow</h4>
        <div class="info-text" style="background:#f0f8ff;">
          <p><strong>Scenario:</strong> You've analysed 200 E. coli genomes and want to write a paper on "Genomic epidemiology of ESBL-producing E. coli in Ghana".</p>
          <p><strong>Step 1 – Let AI generate an outline:</strong><br>
          <em>"Based on the JSON data, suggest a logical structure for a scientific paper covering: population structure (MLST, serotypes), AMR gene prevalence, pathotype distribution, and high-risk combinations."</em></p>
          <p><strong>Step 2 – Draft the Methods section:</strong><br>
          <em>"Write a methods paragraph describing how we used AMRFinderPlus, ABRicate (VFDB, PlasmidFinder, EcoH, BacMet2), and the custom EcoliTyper pipeline to classify pathotypes."</em></p>
          <p><strong>Step 3 – Summarise key results in plain English:</strong><br>
          <em>"From the JSON, what are the top 5 most common AMR genes and their percentages? Format as bullet points."</em></p>
          <p><strong>Step 4 – Generate a draft discussion:</strong><br>
          <em>"Based on the prevalence of blaCTX-M (45%), mcr-1 (3%), and the identification of ST131 (30% of isolates), draft a discussion paragraph on the public health implications."</em></p>
          <p><strong>Step 5 – Polish and cross-reference:</strong><br>
          Ask AI to rephrase, shorten, or expand specific sentences, but <strong>manually verify every number and citation</strong>.</p>
        </div>

        <h4>🎓 For Thesis Writers – Let AI Be Your Writing Coach</h4>
        <ul style="margin-left:20px;">
          <li><strong>Generate figure/table captions:</strong> <em>"Write a descriptive caption for the pathotype prevalence table (STEC 20%, EPEC 15%, …)."</em></li>
          <li><strong>Rewrite passive voice to active:</strong> <em>"Change this sentence to active voice: 'The blaCTX-M gene was detected in 45% of isolates.'"</em></li>
          <li><strong>Create a graphical abstract concept:</strong> <em>"Describe a graphical abstract that visualises the main findings: ST131 dominance, high ESBL rate, and presence of mcr-1."</em></li>
          <li><strong>Summarise limitations:</strong> <em>"Based on the QC metrics (e.g., N50 values), what limitations should I mention in the discussion?"</em></li>
        </ul>

        <h4>⚠️ Ethical Guardrails</h4>
        <div class="info-text" style="background:#fff3cd;border-left-color:#ffc107;">
          <ul style="margin:5px 0 5px 20px;">
            <li><strong>❌ Never paste unreviewed AI output directly into a manuscript</strong> – treat it as a first draft, not a final product.</li>
            <li><strong>✅ Always cite the AI tool if you use verbatim text</strong> (check journal policies – many now require disclosure).</li>
            <li><strong>🔍 Fact-check every number, gene name, and reference</strong> – AI can invent plausible-sounding nonsense.</li>
            <li><strong>👥 Keep human experts in the loop</strong> – AI is a co-pilot, not the captain.</li>
            <li><strong>🌍 Use AI to democratise science</strong> – non-native English speakers can use it to improve clarity, but always acknowledge its use.</li>
          </ul>
        </div>

        <div class="info-text" style="margin-top:20px;">
          <i class="fas fa-globe-africa"></i> <strong>Our Philosophy:</strong> We are in the era of AI-augmented science. Used ethically, AI can accelerate public health research, reduce language barriers, and empower researchers in low-resource settings. But it must never replace human curiosity, critical thinking, or domain expertise. <strong>You are the expert – AI is your enthusiastic, occasionally forgetful assistant.</strong>
        </div>'''

    def _sec_citation(self):
        entries = [
            ('#2196F3', '🧬 EcoliTyper (Main Citation)',
             'Beckley B, Amarh V. EcoliTyper: a species-optimized computational pipeline for comprehensive genotyping and surveillance of <em>Escherichia coli</em>. <em>BMC Bioinformatics</em>. 2026.',
             'https://doi.org/10.1186/s12859-026-06529-6'),
            ('#e91e63', '🔬 CHTyper',
             'Roer L, et al. CHTyper, a Web Tool for Subtyping of Extraintestinal Pathogenic <em>Escherichia coli</em> Based on the <em>fumC</em> and <em>fimH</em> Alleles. <em>J Clin Microbiol</em>. 2018;56(4):e00063-18.',
             'https://doi.org/10.1128/JCM.00063-18'),
            ('#4caf50', '🧬 AMRFinderPlus',
             'Feldgarden M, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. <em>Sci Rep</em>. 2021;11(1):12728.',
             'https://doi.org/10.1038/s41598-021-91456-0'),
            ('#ff9800', '🛡️ CARD 2023',
             'Alcock BP, et al. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. <em>Nucleic Acids Res</em>. 2023;51(D1):D690-D699.',
             'https://doi.org/10.1093/nar/gkac920'),
            ('#9c27b0', '🧬 ARG-ANNOT',
             'Gupta SK, et al. ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. <em>Antimicrob Agents Chemother</em>. 2014;58(1):212-20.',
             'https://doi.org/10.1128/AAC.01310-13'),
            ('#00bcd4', '🌊 MEGARes 3.0',
             'Bonin N, et al. MEGARes and AMR++, v3.0: an updated comprehensive database of antimicrobial resistance determinants and an improved software pipeline for classification using high-throughput sequencing. <em>Nucleic Acids Res</em>. 2023;51(D1):D744-D752.',
             'https://doi.org/10.1093/nar/gkac1047'),
            ('#f44336', '🦠 VFDB 2025',
             'Zhou S, et al. VFDB 2025: an integrated resource for exploring anti-virulence compounds. <em>Nucleic Acids Res</em>. 2025;53(D1):D871-D877.',
             'https://doi.org/10.1093/nar/gkae968'),
            ('#3f51b5', '🧬 PlasmidFinder',
             'Carattoli A, et al. <em>In silico</em> detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. <em>Antimicrob Agents Chemother</em>. 2014;58(7):3895-903.',
             'https://doi.org/10.1128/AAC.02412-14'),
            ('#ffc107', '🔬 EcoH / Serotyping',
             'Ingle DJ, et al. In silico serotyping of <em>E. coli</em> from short read data identifies limited novel O-loci but extensive diversity of O:H serotype combinations within and between pathogenic lineages. <em>Microb Genom</em>. 2016;2(7):e000064.',
             'https://doi.org/10.1099/mgen.0.000064'),
            ('#009688', '🧪 BacMet',
             'Pal C, et al. BacMet: antibacterial biocide and metal resistance genes database. <em>Nucleic Acids Res</em>. 2014;42(D1):D737-43.',
             'https://doi.org/10.1093/nar/gkt1252'),
            ('#ff5722', '🐍 Biopython',
             'Cock PJ, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. <em>Bioinformatics</em>. 2009;25(11):1422-3.',
             'https://doi.org/10.1093/bioinformatics/btp163'),
            ('#673ab7', '🌐 PubMLST / BIGSdb',
             'Jolley KA, et al. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. <em>Wellcome Open Res</em>. 2018;3:124.',
             'https://doi.org/10.12688/wellcomeopenres.14826.1'),
            ('#cddc39', '🔍 ResFinder',
             'Florensa AF, et al. ResFinder – an open online resource for identification of antimicrobial resistance genes in next-generation sequencing data and prediction of phenotypes from genotypes. <em>Microb Genom</em>. 2022;8(1):000748.',
             'https://doi.org/10.1099/mgen.0.000748'),
            ('#8bc34a', '📦 ABRicate',
             'Seemann T. ABRicate: mass screening of contigs for antibiotic resistance genes. <a href="https://github.com/tseemann/abricate" target="_blank">github.com/tseemann/abricate</a>.',
             ''),
            ('#26c6da', '🔧 MLST',
             'Seemann T. MLST: Scan contig files against PubMLST typing schemes. <a href="https://github.com/tseemann/mlst" target="_blank">github.com/tseemann/mlst</a>.',
             ''),
            ('#f06292', '🧬 mobileOG-db',
             'Brown CL, et al. mobileOG-db: a Manually Curated Database of Protein Families Mediating the Life Cycle of Bacterial Mobile Genetic Elements. <em>Appl Environ Microbiol</em>. 2022;88(18):e00991-22.',
             'https://doi.org/10.1128/aem.00991-22'),
            ('#00897b', '🌿 Clermont phylotyping',
             'Clermont O, et al. The Clermont <em>Escherichia coli</em> phylo-typing method revisited: improvement of specificity and detection of new phylo-groups. <em>Environ Microbiol Rep</em>. 2013;5(1):58-65.',
             'https://doi.org/10.1111/1758-2229.12019'),
            ('#ffb300', '🔬 Joensen serotyping',
             'Joensen KG, et al. Rapid and easy <em>in silico</em> serotyping of <em>Escherichia coli</em> using whole genome sequencing (WGS) data. <em>J Clin Microbiol</em>. 2015;53(8):2410-2426.',
             'https://doi.org/10.1128/JCM.00008-15'),
            ('#5c6bc0', '🔬 KMA',
             'Clausen PTLC, et al. Rapid and precise alignment of raw reads against redundant databases with KMA. <em>BMC Bioinformatics</em>. 2018;19:307.',
             'https://doi.org/10.1186/s12859-018-2336-6'),
            ('#e91e63', '⚡ BLAST+',
             'Camacho C, et al. BLAST+: architecture and applications. <em>BMC Bioinformatics</em>. 2009;10:421.',
             'https://doi.org/10.1186/1471-2105-10-421'),
        ]
        cards = ''
        for color, name, text, url in entries:
            link = f'<a href="{url}" target="_blank" style="color:{color};font-weight:bold;">{url}</a>' if url else ''
            cards += (f'<div style="background:{color}1a;padding:20px;border-radius:12px;margin-bottom:15px;'
                      f'border-left:5px solid {color};">'
                      f'<strong style="color:{color};font-size:1.02em;">{name}</strong><br>'
                      f'<span style="font-size:.95em;">{text} {link}</span></div>')
        return f'''
        <div class="alert-box alert-info" style="background:linear-gradient(135deg,#667eea 0%,#764ba2 100%);color:white;border-left:none;">
          <i class="fas fa-quote-right fa-2x"></i>
          <div>
            <h3 style="color:white;">📚 How to Cite EcoliTyper and Its Dependencies</h3>
            <p style="color:white;">If you use this pipeline in your research, please cite the main tool and the relevant third-party tools and databases below.</p>
          </div>
        </div>
        <div style="margin:20px 0;">{cards}</div>'''

    def _sec_funding(self):
        return '''
        <div style="background:linear-gradient(135deg,#43a047 0%,#1b5e20 100%);padding:25px;
                    border-radius:15px;margin-bottom:25px;color:white;text-align:center;">
          <i class="fas fa-coffee fa-3x" style="margin-bottom:10px;"></i>
          <h2 style="color:white;margin:0;">☕ No Grants, Just Passion (and a lot of caffeine)</h2>
          <p style="margin-top:10px;font-size:1.1em;"><strong>EcoliTyper</strong> is an independent, unfunded project born out of late nights, curiosity, and the urgent need for accessible genomic surveillance of <em>E. coli</em> in low-resource settings.</p>
          <p>No sponsors, no institutional backing – just a laptop, a never-empty coffee mug, and a burning desire to help researchers fight AMR.</p>
        </div>

        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:20px;margin:25px 0;">
          <div style="background:#fff3e0;padding:20px;border-radius:12px;border-left:5px solid #ff9800;">
            <i class="fas fa-star fa-2x" style="color:#ff9800;"></i>
            <h3 style="color:#e65100;margin:10px 0;">⭐ Star us on GitHub</h3>
            <p>One click = a million smiles. It costs nothing but makes us feel like rockstars.</p></div>
          <div style="background:#e8f5e9;padding:20px;border-radius:12px;border-left:5px solid #4caf50;">
            <i class="fas fa-bug fa-2x" style="color:#4caf50;"></i>
            <h3 style="color:#1b5e20;margin:10px 0;">🐛 Report bugs</h3>
            <p>If something breaks, tell us. We'll fix it with joy (and maybe a little swearing).</p></div>
          <div style="background:#e3f2fd;padding:20px;border-radius:12px;border-left:5px solid #2196f3;">
            <i class="fas fa-lightbulb fa-2x" style="color:#2196f3;"></i>
            <h3 style="color:#0d47a1;margin:10px 0;">💡 Suggest features</h3>
            <p>Have an idea? We're all ears – and we actually implement them.</p></div>
          <div style="background:#fce4ec;padding:20px;border-radius:12px;border-left:5px solid #e91e63;">
            <i class="fas fa-dna fa-2x" style="color:#e91e63;"></i>
            <h3 style="color:#880e4f;margin:10px 0;">🧬 Share your data</h3>
            <p>Used EcoliTyper? We'd love to hear your story and collaborate on cool science.</p></div>
          <div style="background:#f1f8e9;padding:20px;border-radius:12px;border-left:5px solid #8bc34a;">
            <i class="fas fa-share-alt fa-2x" style="color:#8bc34a;"></i>
            <h3 style="color:#33691e;margin:10px 0;">📢 Spread the word</h3>
            <p>Tweet, blog, or tell your lab mate. The more people use it, the better it gets.</p></div>
          <div style="background:#ede7f6;padding:20px;border-radius:12px;border-left:5px solid #673ab7;">
            <i class="fas fa-envelope fa-2x" style="color:#673ab7;"></i>
            <h3 style="color:#311b92;margin:10px 0;">👋 Say hello</h3>
            <p>Drop an email to <strong>brownbeckley94@gmail.com</strong>. It literally makes our day.</p></div>
        </div>

        <div style="background:#ffebee;padding:20px;border-radius:12px;margin:20px 0;border-left:5px solid #f44336;">
          <i class="fas fa-microbe fa-2x" style="color:#f44336;float:left;margin-right:15px;"></i>
          <p style="margin:0;"><strong>Fun fact:</strong> This project runs on <strong>0% grant money, 100% volunteer tears</strong>. But we're not bitter – we're just <strong>caffeinated</strong>. ☕</p>
          <div style="clear:both;"></div>
        </div>

        <div style="background:linear-gradient(135deg,#0d47a1 0%,#1565c0 100%);padding:25px;border-radius:15px;color:white;text-align:center;">
          <i class="fas fa-hand-holding-heart fa-3x"></i>
          <h3 style="color:white;">🤝 Join the Fight Against AMR</h3>
          <p>We believe open-source tools are the key to beating antimicrobial resistance. <strong>EcoliTyper</strong> is just one part – but we welcome <strong>your contributions</strong>:</p>
          <p>🐍 Code improvements • 📚 Documentation • 🌍 New databases • 🤔 Ideas</p>
          <p style="margin-top:15px;">
            <a href="https://github.com/bbeckley-hub/EcoliTyper" target="_blank" style="background:white;color:#0d47a1;padding:10px 20px;border-radius:25px;text-decoration:none;font-weight:bold;display:inline-block;margin:5px;">
              <i class="fab fa-github"></i> Star on GitHub</a>
            <a href="mailto:brownbeckley94@gmail.com" target="_blank" style="background:#ffc107;color:#0d47a1;padding:10px 20px;border-radius:25px;text-decoration:none;font-weight:bold;display:inline-block;margin:5px;">
              <i class="fas fa-envelope"></i> Contact me</a>
          </p>
          <p style="margin-top:15px;font-size:.9em;"><i class="fas fa-laugh-beam"></i> <strong>P.S.</strong> If you ever meet Brown in person, buy him a coffee – he'll talk your ear off about O157:H7 and ST131, but it's totally worth it.</p>
        </div>'''

    def _sec_calltoaction(self):
        return '''
        <div class="alert-box alert-info">
          <i class="fas fa-globe fa-2x"></i>
          <div>
            <h3>The Global Burden of AMR and Our Call to Action</h3>
            <p>Antimicrobial resistance (AMR) is one of the top global public health threats, with an estimated <strong>1.27 million direct deaths annually</strong>. <em>Escherichia coli</em> is a major contributor, causing urinary tract infections, bloodstream infections, and foodborne outbreaks. Tracking AMR and virulence determinants is essential to inform treatment guidelines and infection control.</p>
            <p>We developed <strong>ECOLITYPER</strong> to empower researchers and clinicians – especially in low-resource settings – to analyse their own sequencing data without extensive bioinformatics expertise.</p>
          </div>
        </div>

        <div style="background:#e8f5e9;padding:20px;border-radius:12px;margin:20px 0;">
          <h3><i class="fas fa-bacterium"></i> ESCAPE AMR – Our Ongoing Project (ESKAPE Pathogens)</h3>
          <p><strong>ECOLITYPER</strong> is one of the first modules of a larger initiative called <strong>ESCAPE AMR</strong> (formerly ESKAPE). We target the notorious <strong>ESKAPE pathogens</strong>:</p>
          <ul>
            <li><strong>E</strong>nterococcus faecium</li>
            <li><strong>S</strong>taphylococcus aureus</li>
            <li><strong>K</strong>lebsiella pneumoniae</li>
            <li><strong>A</strong>cinetobacter baumannii</li>
            <li><strong>P</strong>seudomonas aeruginosa</li>
            <li><strong>E</strong>nterobacter species (including E. coli)</li>
          </ul>
          <p>These bacteria "escape" the effects of antibiotics – hence the name. But we believe the name is also a global call to action:</p>
          <div style="background:#fff3e0;padding:15px;border-radius:8px;margin:15px 0;">
            <p><strong>🔹 E</strong>veryone must join forces – researchers, clinicians, policymakers, and citizens.<br>
            <strong>🔹 S</strong>mart surveillance is our first line of defence. No more guessing – we need genomic data.<br>
            <strong>🔹 K</strong>nowledge must be shared openly. No paywalls, no closed silos.<br>
            <strong>🔹 A</strong>frica bears a heavy AMR burden, but African solutions are already emerging.<br>
            <strong>🔹 P</strong>revention is cheaper than cure. Let's stop resistant infections before they spread.<br>
            <strong>🔹 E</strong>very day we delay, more lives are at stake. The time to act is now, not tomorrow.</p>
          </div>
          <p><i class="fas fa-laugh-squint"></i> <strong>"We didn't choose the name ESKAPE because it sounds cool (though it does). We chose it because it reminds us every single day: we must ESCAPE the AMR crisis – together, urgently, and with the best science we have."</strong><br>
          — Brown Beckley, lead developer</p>
        </div>

        <div style="text-align:center;margin:40px 0;">
          <i class="fas fa-star" style="font-size:3em;color:#ffc107;"></i>
          <h3>🤝 We Invite You to Contribute!</h3>
          <p><strong>If you find this tool useful, please:</strong></p>
          <div class="action-buttons" style="justify-content:center;">
            <a href="https://github.com/bbeckley-hub/ecolityper" target="_blank" class="action-btn btn-primary" style="text-decoration:none;"><i class="fab fa-github"></i> ⭐ Star us on GitHub</a>
            <a href="mailto:brownbeckley94@gmail.com" class="action-btn btn-success" style="text-decoration:none;"><i class="fas fa-envelope"></i> Contact the team</a>
            <a href="https://github.com/bbeckley-hub/ecolityper/issues" target="_blank" class="action-btn btn-warning" style="text-decoration:none;"><i class="fas fa-bug"></i> Report issues</a>
          </div>
          <p style="margin-top:20px;"><i class="fas fa-chalkboard-user"></i> <strong>We welcome collaborations</strong> to adapt this tool for other pathogens and to improve AMR surveillance in Africa and beyond.</p>
          <p><i class="fas fa-hand-holding-heart"></i> If you are a funder or organisation interested in supporting the <strong>ESCAPE AMR</strong> project, please reach out. Together we can build a free, open-source ecosystem for genomic surveillance of all ESKAPE pathogens.</p>
        </div>

        <div style="background:#f8f9fa;padding:15px;border-radius:8px;">
          <i class="fas fa-quote-left"></i> "AMR is a silent pandemic, but we have the tools to fight it – if we share them, if we teach each other, and if we act with urgency. Let's escape the era of untreatable infections."<br>
          <strong>— The ESCAPE AMR Team, University of Ghana Medical School</strong>
        </div>'''

    def _sec_export(self):
        return '''
        <div class="alert-box alert-info">
          <i class="fas fa-download"></i>
          <div>
            <h3>📥 Export Data</h3>
            <p>Download any table as CSV, or the complete dataset as JSON for AI-assisted analysis.</p>
          </div>
        </div>
        <div class="action-buttons">
          <button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table','sample_overview.csv')">
            <i class="fas fa-table"></i> Sample Overview CSV</button>
          <button class="action-btn btn-danger" onclick="exportTableToCSV('amr-table','amr_genes.csv')">
            <i class="fas fa-biohazard"></i> AMR Genes CSV</button>
          <button class="action-btn btn-success" onclick="exportTableToCSV('vir-table','virulence_genes.csv')">
            <i class="fas fa-virus"></i> Virulence CSV</button>
          <button class="action-btn btn-info" onclick="exportTableToCSV('plasmid-table','plasmid_markers.csv')">
            <i class="fas fa-plug"></i> Plasmid CSV</button>
          <button class="action-btn btn-secondary" onclick="exportTableToCSV('bac-table','bacmet_genes.csv')">
            <i class="fas fa-flask"></i> BACMET CSV</button>
          <button class="action-btn btn-warning" onclick="exportTableToCSV('mutation-table','mutations.csv')">
            <i class="fas fa-dna"></i> Mutations CSV</button>
          <button class="action-btn btn-primary" onclick="exportTableToCSV('mge-table','mge_profile.csv')">
            <i class="fas fa-mobile-alt"></i> MGE Profile CSV</button>
          <button class="action-btn btn-success" onclick="location.href='genius_ecoli_ultimate_gene_centric_report.json'">
            <i class="fas fa-file-code"></i> Complete JSON</button>
        </div>'''


# =============================================================================
# ORCHESTRATOR
# =============================================================================
class EcoliReporter:
    """Coordinates file discovery, parsing, analysis, and report generation."""

    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "GENIUS_ECOLI_ULTIMATE_GENE_CENTRIC_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = EcoliParser()
        self.analyzer = EcoliAnalyzer()
        self.html_generator = EcoliHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "EcoliTyper Ultimate Reporter",
            "version": "4.3.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir),
        }

    def find_files(self) -> Dict[str, Any]:
        print("Searching for input files...")
        found = {
            'master_tsv': None, 'qc_html': None, 'amrfinder_html': None,
            'abricate_htmls': [], 'mutation_html': None, 'mge_html': None,
        }
        for cand in ('ecolityper_comprehensive_report.tsv', 'ecoli_comprehensive_report.tsv'):
            p = self.input_dir / cand
            if p.exists():
                found['master_tsv'] = p
                print(f"    Master TSV: {p.name}")
                break
        for cand in ('FASTA_QC_summary.html', 'fasta_qc_summary.html'):
            p = self.input_dir / cand
            if p.exists():
                found['qc_html'] = p
                print(f"    QC HTML: {p.name}")
                break
        for f in self.input_dir.glob("**/*.html"):
            if 'amrfinder_summary' in f.name.lower() or 'amr_summary' in f.name.lower():
                found['amrfinder_html'] = f
                print(f"    AMRfinder: {f.name}")
                break
        for f in self.input_dir.glob("**/*.html"):
            if any(f'{db}_summary_report' in f.name.lower() for db in self.parser.abricate_databases):
                found['abricate_htmls'].append(f)
        if found['abricate_htmls']:
            print(f"    ABRicate: {len(found['abricate_htmls'])} file(s)")
        for f in self.input_dir.glob("**/mutation_summary.html"):
            found['mutation_html'] = f
            print(f"    Mutation: {f.name}")
            break
        for f in self.input_dir.glob("**/*.html"):
            if 'mge' in f.name.lower() and 'summary' in f.name.lower():
                found['mge_html'] = f
                print(f"    MGE: {f.name}")
                break
        return found

    def integrate(self, files: Dict[str, Any]) -> Dict[str, Any]:
        print("\nIntegrating data...")
        data = {
            'metadata': self.metadata, 'samples': {}, 'patterns': {},
            'gene_centric': {}, 'qc_data': {}, 'mutation_data': {}, 'mge_data': {},
        }
        typing_data = self.parser.load_master_tsv(files['master_tsv']) if files['master_tsv'] else {}
        if files['qc_html']:
            data['qc_data'] = self.parser.parse_qc_report(files['qc_html'])
        if files['mutation_html']:
            data['mutation_data'] = self.parser.parse_mutation_summary_html(files['mutation_html'])
        total_samples = len(typing_data) if typing_data else 0
        amr_by_sample, amr_freq = {}, {}
        if files['amrfinder_html']:
            amr_by_sample, amr_freq = self.parser.parse_amrfinder_report(files['amrfinder_html'], total_samples)
        abricate_by_sample = defaultdict(dict)
        abricate_freq = {}
        for f in files['abricate_htmls']:
            db, by_sample, freq = self.parser.parse_abricate_database_report(f, total_samples)
            if db != 'unknown':
                for s, g in by_sample.items():
                    abricate_by_sample[s][db] = g
                abricate_freq[db] = freq
        if files['mge_html']:
            mge_raw = self.parser.parse_mge_summary_html(files['mge_html'])
            mge_raw['analysis'] = self.analyzer.compute_mge_stats(mge_raw)
            data['mge_data'] = mge_raw
        all_samples = set(typing_data) | set(amr_by_sample) | set(abricate_by_sample) \
                      | set(data['qc_data']) | set(data['mge_data'].get('per_sample', {}))
        all_samples = sorted(all_samples)
        if not all_samples:
            print("No samples found in any report")
            return {}
        print(f"Found {len(all_samples)} unique samples")
        for sample in all_samples:
            t = typing_data.get(sample, {})
            vir_genes = []
            for db in ('vfdb', 'ecoli_vf'):
                vir_genes.extend(abricate_by_sample.get(sample, {}).get(db, []))
            vir_genes = sorted(set(vir_genes))
            patho = self.analyzer.classify_pathotype(vir_genes)
            mge_rec = data['mge_data'].get('per_sample', {}).get(sample, {})
            mge_hits = 0
            if isinstance(mge_rec, dict):
                for k in ('mobileOG Hits', 'Total mobileOG hits', 'Total hits'):
                    v = mge_rec.get(k)
                    if isinstance(v, (int, float)):
                        mge_hits = v
                        break
            data['samples'][sample] = {
                'mlst': {'ST': t.get('ST', 'ND')},
                'serotype': {
                    'Serotype': t.get('Serotype', 'ND'),
                    'O_Type': t.get('O_Type', 'ND'),
                    'H_Type': t.get('H_Type', 'ND'),
                    'Serotype_Confidence': t.get('Serotype_Confidence', 'ND'),
                },
                'phylogrouping': {'Clermont_Type': t.get('Phylogroup', 'ND')},
                'chtyper': {
                    'FumC_Type': t.get('FumC_Type', 'ND'),
                    'FimH_Type': t.get('FimH_Type', 'ND'),
                    'CH_Type': t.get('CH_Type', 'ND'),
                },
                'amrfinder': {'all_genes': amr_by_sample.get(sample, [])},
                'abricate_databases': abricate_by_sample.get(sample, {}),
                'amr_genes': list(amr_by_sample.get(sample, []))
                            + list(abricate_by_sample.get(sample, {}).get('card', []))
                            + list(abricate_by_sample.get(sample, {}).get('resfinder', []))
                            + list(abricate_by_sample.get(sample, {}).get('argannot', []))
                            + list(abricate_by_sample.get(sample, {}).get('megares', []))
                            + list(abricate_by_sample.get(sample, {}).get('ncbi', [])),
                'virulence_genes': vir_genes,
                'pathotype': patho,
                'mge_hits': mge_hits,
            }
        data['gene_frequencies'] = {'amrfinder': amr_freq, 'abricate': abricate_freq}
        print("Building gene-centric + pattern tables...")
        data['gene_centric'] = self.analyzer.create_gene_centric_tables(data, len(all_samples))
        data['patterns'] = self.analyzer.create_cross_genome_patterns(data)
        self.html_generator._current_samples_data = data['samples']
        return data

    def write_json(self, data: Dict[str, Any]) -> Path:
        print("\nWriting JSON report...")
        out = self.output_dir / "genius_ecoli_ultimate_gene_centric_report.json"

        def serialise(obj):
            if obj is None or isinstance(obj, (str, int, float, bool)):
                return obj
            if isinstance(obj, (list, tuple, set)):
                return [serialise(x) for x in obj]
            if isinstance(obj, (Counter, defaultdict)):
                return {str(k): serialise(v) for k, v in obj.items()}
            if isinstance(obj, dict):
                return {str(k): serialise(v) for k, v in obj.items()}
            if isinstance(obj, Path):
                return str(obj)
            if hasattr(obj, 'isoformat'):
                return obj.isoformat()
            return str(obj)

        with open(out, 'w', encoding='utf-8') as f:
            json.dump(serialise(data), f, indent=2, ensure_ascii=False)
        print(f"    JSON saved: {out}")
        return out

    def write_csvs(self, data: Dict[str, Any]) -> None:
        print("\nWriting CSV reports...")
        rows = []
        for s, d in data['samples'].items():
            rows.append({
                'Sample': s,
                'ST': d['mlst']['ST'],
                'Serotype': d['serotype']['Serotype'],
                'O_Type': d['serotype']['O_Type'],
                'H_Type': d['serotype']['H_Type'],
                'Phylogroup': d['phylogrouping']['Clermont_Type'],
                'CH_Type': d['chtyper']['CH_Type'],
                'Pathotype': d['pathotype'],
                'AMR_Count': len(d['amr_genes']),
                'Virulence_Count': len(d['virulence_genes']),
                'MGE_Hits': d['mge_hits'],
            })
        pd.DataFrame(rows).to_csv(self.output_dir / "sample_overview.csv", index=False)
        gc = data.get('gene_centric', {})
        total = len(data['samples']) or 1
        for cat, fname in [
            ('amr_databases', 'amr_genes.csv'),
            ('virulence_databases', 'virulence_genes.csv'),
            ('bacmet_databases', 'bacmet_genes.csv'),
            ('plasmid_databases', 'plasmid_markers.csv'),
        ]:
            out = []
            for db, genes in gc.get(cat, {}).items():
                for g in genes:
                    out.append({
                        'Gene': g['gene'], 'Database': g['database'],
                        'Count': g['count'],
                        'Percentage': f"{(g['count'] / total) * 100:.1f}%",
                        'Genomes': ';'.join(g.get('genomes', [])),
                    })
            if out:
                pd.DataFrame(out).to_csv(self.output_dir / fname, index=False)
        muts = data.get('mutation_data', {}).get('mutations', [])
        if muts:
            pd.DataFrame([{
                'Gene': m['gene'], 'Mutation': m['mutation'],
                'Class': m['class'], 'Subclass': m['subclass'],
                'Count': m['count'], 'Genomes': ';'.join(m['genomes']),
            } for m in muts]).to_csv(self.output_dir / "mutations.csv", index=False)
        mge_ps = data.get('mge_data', {}).get('per_sample', {})
        if mge_ps:
            pd.DataFrame([{'Sample': s, **v} for s, v in mge_ps.items()]) \
                .to_csv(self.output_dir / "mge_profile.csv", index=False)
        if data.get('qc_data'):
            pd.DataFrame([{'Sample': s, **m} for s, m in data['qc_data'].items()]) \
                .to_csv(self.output_dir / "fasta_qc.csv", index=False)
        print("    CSVs written")

    def run(self) -> bool:
        print("=" * 80)
        print("GENIUS E. COLI ULTIMATE REPORTER v4.0.0")
        print("=" * 80)
        print(f"Input:  {self.input_dir}")
        print(f"Output: {self.output_dir}")
        files = self.find_files()
        if not any([files['master_tsv'], files['qc_html'], files['amrfinder_html'],
                    files['abricate_htmls'], files['mutation_html'], files['mge_html']]):
            print("No recognized input files found")
            return False
        data = self.integrate(files)
        if not data:
            return False
        self.write_json(data)
        self.write_csvs(data)
        self.html_generator.generate_main_report(data, self.output_dir)
        n = len(data['samples'])
        patterns = data['patterns']
        gc = data['gene_centric']
        total_amr = sum(len(g) for g in gc.get('amr_databases', {}).values())
        total_vir = sum(len(g) for g in gc.get('virulence_databases', {}).values())
        high = len(patterns.get('high_risk_combinations', []))
        patho = len([k for k in patterns.get('pathotype_distribution', {}) if k != 'none'])
        print("\n" + "=" * 80)
        print("ANALYSIS COMPLETE")
        print("=" * 80)
        print(f"   Samples:         {n}")
        print(f"   AMR genes:       {total_amr}")
        print(f"   Virulence genes: {total_vir}")
        print(f"   Pathotypes:      {patho}")
        print(f"   High-risk:       {high}")
        print(f"   Output:          {self.output_dir}")
        print("=" * 80)
        return True


# =============================================================================
# MAIN
# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        description='EcoliTyper Ultimate Reporter v4.0.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python genius_ultimate_gene_centric_reporter.py -i /path/to/reports

Author: Brown Beckley <brownbeckley94@gmail.com>
        """)
    parser.add_argument('-i', '--input-dir', required=True,
                        help='Directory containing EcoliTyper reports')
    parser.add_argument('-o', '--output-dir', help='Custom output directory')
    args = parser.parse_args()
    inp = Path(args.input_dir)
    if not inp.exists():
        print(f"Input directory not found: {inp}")
        sys.exit(1)
    reporter = EcoliReporter(inp)
    if args.output_dir:
        reporter.output_dir = Path(args.output_dir)
        reporter.output_dir.mkdir(parents=True, exist_ok=True)
    if not reporter.run():
        sys.exit(1)


if __name__ == "__main__":
    main()