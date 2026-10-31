#!/usr/bin/env python3
"""
EcoliTyper Hybrid Reporter
==========================
Interactive HTML report for Escherichia coli genomic analysis combining
gene-centric typing tabs with sample-centric lazy-loaded isolate boxes.

Report structure
----------------
Gene-centric tabs:
    MLST, Serotype, CHTyper, Phylogroup, Pathotype, Patterns.
Sample-centric tabs (lazy-loaded isolate boxes):
    AMR, Virulence, Plasmid, Bacmet, Mutations, MGE.
Support tabs:
    Summary, Sample Overview, FASTA QC, Compare & Cluster, AI Guide,
    Citation, Funding, Call to Action, Export.

Inputs
------
Primary TSV sources:
    ecolityper_comprehensive_report.tsv
    ecoli_amrfinder_summary.tsv
    ecoli_<db>_abricate_summary.tsv
    mutation_summary.tsv
    ecolityper_mge_all_hits.tsv
Secondary HTML source:
    FASTA_QC_summary.html

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School, Department of Medical Biochemistry
License: MIT
"""

import json
import re
import argparse
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

import pandas as pd

from bs4 import BeautifulSoup


NOT_ASSIGNED = "Not Assigned"


# =============================================================================
# Parser
# =============================================================================
class EcoliHTMLParser:
    """Parses TSV summaries and HTML reports produced by EcoliTyper modules."""

    def __init__(self):
        self.abricate_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ecoh', 'ecoli_vf', 'bacmet2',
        ]
        self.abricate_tsv_files = {
            'vfdb': 'ecoli_vfdb_abricate_summary.tsv',
            'ecoli_vf': 'ecoli_ecoli_vf_abricate_summary.tsv',
            'plasmidfinder': 'ecoli_plasmidfinder_abricate_summary.tsv',
            'bacmet2': 'ecoli_bacmet2_abricate_summary.tsv',
            'resfinder': 'ecoli_resfinder_abricate_summary.tsv',
            'card': 'ecoli_card_abricate_summary.tsv',
            'argannot': 'ecoli_argannot_abricate_summary.tsv',
            'megares': 'ecoli_megares_abricate_summary.tsv',
            'ncbi': 'ecoli_ncbi_abricate_summary.tsv',
            'ecoh': 'ecoli_ecoh_abricate_summary.tsv',
        }

    def normalize_sample_id(self, sample_id: str) -> str:
        """Strip extension and directory prefix from a sample identifier."""
        s = str(sample_id).strip()
        for ext in ('.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff',
                    '.txt', '.tsv', '.csv'):
            if s.endswith(ext):
                s = s[:-len(ext)]
        if '/' in s or '\\' in s:
            s = Path(s).name
        return s.strip()

    def clean_st(self, st_str: str) -> str:
        """Normalise an MLST sequence type to the form 'STxxx'."""
        if not st_str or str(st_str).strip() in ('ND', 'nan', 'None', ''):
            return 'ND'
        s = re.sub(r'^(ST)+', 'ST', str(st_str).strip(), flags=re.I)
        if not s.upper().startswith('ST'):
            s = 'ST' + s
        return s

    def _clean_cell(self, v) -> str:
        """Return a stripped string, or an empty string for missing values."""
        if pd.isna(v):
            return ''
        return str(v)

    def _row_to_dict(self, row) -> Dict[str, Any]:
        """Convert a pandas row to a plain dict, replacing NaN with None."""
        return {k: (None if pd.isna(v) else v) for k, v in row.to_dict().items()}

    def load_master_typing_tsv(self, input_dir: Path) -> Dict[str, Dict]:
        """Load typing columns from ecolityper_comprehensive_report.tsv."""
        path = input_dir / 'ecolityper_comprehensive_report.tsv'
        if not path.exists():
            print("    Master TSV not found; typing columns will be 'ND'")
            return {}
        df = pd.read_csv(path, sep='\t', dtype=str).fillna('ND')
        typing = {}
        for _, row in df.iterrows():
            sample = self.normalize_sample_id(row.get('Sample_ID', ''))
            if not sample:
                continue
            typing[sample] = {
                'ST': self.clean_st(row.get('MLST_ST', 'ND')),
                'Allele_Profile': str(row.get('MLST_Allele_Profile', 'ND')).strip(),
                'Serotype': str(row.get('Serotype', 'ND')).strip(),
                'O_Type': str(row.get('O_Type', 'ND')).strip(),
                'H_Type': str(row.get('H_Type', 'ND')).strip(),
                'Serotype_Confidence': str(row.get('Serotype_Confidence', 'ND')).strip(),
                'Clermont_Type': str(row.get('Phylogroup', 'ND')).strip(),
                'FumC_Type': str(row.get('FumC_Type', 'ND')).strip(),
                'FimH_Type': str(row.get('FimH_Type', 'ND')).strip(),
                'CH_Type': str(row.get('CH_Type', 'ND')).strip(),
                'FumC_Identity': str(row.get('FumC_Identity', 'ND')).strip(),
                'FimH_Identity': str(row.get('FimH_Identity', 'ND')).strip(),
                'FumC_Coverage': str(row.get('FumC_Coverage', 'ND')).strip(),
                'FimH_Coverage': str(row.get('FimH_Coverage', 'ND')).strip(),
                'Status': str(row.get('Status', 'ND')).strip(),
            }
        print(f"    Loaded {len(typing)} samples from master TSV")
        return typing

    def load_amrfinder_from_tsv(self, input_dir: Path) -> Tuple[Dict[str, List[Dict]], Dict[str, int]]:
        """Load per-sample AMRFinderPlus gene details from TSV."""
        path = input_dir / 'ecoli_amrfinder_summary.tsv'
        if not path.exists():
            return {}, {}
        df = pd.read_csv(path, sep='\t')
        per_sample = defaultdict(list)
        gene_counts: Counter = Counter()
        for _, row in df.iterrows():
            sample = self.normalize_sample_id(row.get('Genome', ''))
            if not sample:
                continue
            gene_dict = self._row_to_dict(row)
            for key in ('Element symbol', 'Element_symbol', 'Gene', 'gene_symbol'):
                if key in gene_dict and gene_dict[key]:
                    gene_dict['gene'] = gene_dict.pop(key)
                    break
            if not gene_dict.get('gene'):
                continue
            per_sample[sample].append(gene_dict)
            gene_counts[gene_dict['gene']] += 1
        if per_sample:
            print(f"    Loaded AMRFinderPlus details for {len(per_sample)} samples")
        return dict(per_sample), dict(gene_counts)

    def load_abricate_from_tsv(self, input_dir: Path) -> Tuple[Dict[str, Dict[str, List[Dict]]], Dict[str, Dict[str, int]]]:
        """Load per-database, per-sample ABRicate gene details from TSVs."""
        per_sample_db = defaultdict(lambda: defaultdict(list))
        db_gene_counts: Dict[str, Counter] = defaultdict(Counter)
        loaded = 0
        for db, fname in self.abricate_tsv_files.items():
            path = input_dir / fname
            if not path.exists():
                continue
            df = pd.read_csv(path, sep='\t')
            sample_col = next((c for c in ('genome', 'file', 'sample') if c in df.columns), None)
            if not sample_col:
                continue
            for _, row in df.iterrows():
                sample = self.normalize_sample_id(row[sample_col])
                if not sample:
                    continue
                gene_dict = self._row_to_dict(row)
                if not gene_dict.get('gene'):
                    continue
                per_sample_db[sample][db].append(gene_dict)
                db_gene_counts[db][gene_dict['gene']] += 1
            loaded += 1
        if loaded:
            print(f"    Loaded ABRicate details for {len(per_sample_db)} samples across {loaded} database(s)")
        return dict(per_sample_db), {k: dict(v) for k, v in db_gene_counts.items()}

    def load_mutations_from_tsv(self, input_dir: Path) -> Dict[str, List[Dict]]:
        """Load per-sample point mutations from TSV."""
        path = input_dir / 'mutation_summary.tsv'
        if not path.exists():
            return {}
        df = pd.read_csv(path, sep='\t')
        per_sample = defaultdict(list)
        for _, row in df.iterrows():
            sample = self.normalize_sample_id(row.get('genome', ''))
            if not sample:
                continue
            per_sample[sample].append({
                'gene': self._clean_cell(row.get('gene_symbol', row.get('gene', ''))),
                'mutation': self._clean_cell(row.get('element_name', '')),
                'class': self._clean_cell(row.get('class', '')),
                'subclass': self._clean_cell(row.get('subclass', '')),
                'contig': self._clean_cell(row.get('contig_id', '')),
                'start': self._clean_cell(row.get('start', '')),
                'stop': self._clean_cell(row.get('stop', '')),
                'strand': self._clean_cell(row.get('strand', '')),
                'coverage': self._clean_cell(row.get('coverage', '')),
                'identity': self._clean_cell(row.get('identity', '')),
                'accession': self._clean_cell(row.get('accession', '')),
            })
        if per_sample:
            print(f"    Loaded mutations for {len(per_sample)} samples")
        return dict(per_sample)

    def load_mge_from_tsv(self, input_dir: Path) -> Dict[str, List[Dict]]:
        """Load per-sample mobileOG-db hits from TSV."""
        path = input_dir / 'ecolityper_mge_all_hits.tsv'
        if not path.exists():
            path = input_dir / 'ecoli_mge_all_hits.tsv'
        if not path.exists():
            return {}
        df = pd.read_csv(path, sep='\t')
        per_sample = defaultdict(list)
        for _, row in df.iterrows():
            sample = self.normalize_sample_id(row.get('Sample', ''))
            if not sample:
                continue
            per_sample[sample].append({
                'rank': self._clean_cell(row.get('Rank', '')),
                'mobileog_id': self._clean_cell(row.get('mobileOG_ID', '')),
                'gene': self._clean_cell(row.get('Gene', '')),
                'contig': self._clean_cell(row.get('Contig', '')),
                'start': self._clean_cell(row.get('Start', '')),
                'stop': self._clean_cell(row.get('Stop', '')),
                'strand': self._clean_cell(row.get('Strand', '')),
                'major_category': self._clean_cell(row.get('Major_Category', '')),
                'minor_category': self._clean_cell(row.get('Minor_Category', '')),
                'percent_identity': self._clean_cell(row.get('Percent_Identity', '')),
                'align_over_qlen': self._clean_cell(row.get('Align_Over_Qlen', '')),
                'coverage_percent': self._clean_cell(row.get('Coverage_Percent', '')),
            })
        if per_sample:
            print(f"    Loaded MGE hits for {len(per_sample)} samples")
        return dict(per_sample)

    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
        """Extract the Nth HTML table in a string as a DataFrame."""
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

    def parse_qc_report(self, file_path: Path) -> Dict[str, Dict]:
        """Parse the FASTA QC summary HTML."""
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


# =============================================================================
# Analyzer
# =============================================================================
class EcoliDataAnalyzer:
    """Builds gene-centric tables, cross-genome patterns, and pathotype calls."""

    def __init__(self):
        self.critical_amr_genes = {
            'blaCTX-M', 'blaSHV', 'blaTEM', 'blaPER', 'blaVEB', 'blaBEL', 'blaGES',
            'blaKPC', 'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP', 'blaSIM', 'blaSPM',
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6',
            'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
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
        """Return a pathotype label inferred from VFDB gene content."""
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
        """Aggregate per-sample gene lists into gene-centric tables."""
        gc = {
            'amr_databases': {},
            'virulence_databases': {},
            'bacmet_databases': {},
            'plasmid_databases': {},
            'database_stats': {},
        }
        samples_data = integrated.get('samples', {})
        buckets = {
            'amr': Counter(), 'virulence': Counter(),
            'bacmet': Counter(), 'plasmid': Counter(),
        }
        genomes = {k: defaultdict(list) for k in buckets}

        for sample, data in samples_data.items():
            for gene in data.get('amr_genes', []):
                buckets['amr'][gene] += 1
                if sample not in genomes['amr'][gene]:
                    genomes['amr'][gene].append(sample)
            for gene in data.get('virulence_genes', []):
                buckets['virulence'][gene] += 1
                if sample not in genomes['virulence'][gene]:
                    genomes['virulence'][gene].append(sample)
            for gene in data.get('bacmet_genes', []):
                buckets['bacmet'][gene] += 1
                if sample not in genomes['bacmet'][gene]:
                    genomes['bacmet'][gene].append(sample)
            for gene in data.get('plasmid_genes', []):
                buckets['plasmid'][gene] += 1
                if sample not in genomes['plasmid'][gene]:
                    genomes['plasmid'][gene].append(sample)

        mapping = [
            ('amr', 'amr_databases', 'amrfinder', 'AMRfinder'),
            ('virulence', 'virulence_databases', 'vfdb', 'VFDB'),
            ('bacmet', 'bacmet_databases', 'bacmet2', 'BACMET2'),
            ('plasmid', 'plasmid_databases', 'plasmidfinder', 'PLASMIDFINDER'),
        ]
        for bucket_key, cat, db_key, db_label in mapping:
            if not buckets[bucket_key]:
                continue
            glist = []
            for gene, count in buckets[bucket_key].most_common():
                pct = (count / total_samples) * 100 if total_samples else 0
                glist.append({
                    'gene': gene, 'database': db_label,
                    'count': count, 'percentage': round(pct, 2),
                    'frequency_display': f"{count} ({pct:.1f}%)",
                    'genomes': genomes[bucket_key][gene],
                })
            gc[cat][db_key] = glist

        for cat in ('amr_databases', 'virulence_databases',
                    'bacmet_databases', 'plasmid_databases'):
            for db_name, genes in gc.get(cat, {}).items():
                gc['database_stats'][db_name] = {
                    'total_genes': len(genes),
                    'total_occurrences': sum(g['count'] for g in genes),
                }
        return gc

    def create_cross_genome_patterns(self, integrated: Dict) -> Dict[str, Any]:
        """Build every distribution, combination, and co-occurrence table."""
        P = {
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
            'st_fimh_combinations': defaultdict(list),
            'st_fumc_combinations': defaultdict(list),
            'st_fumc_fimh_combinations': defaultdict(list),
            'st_oh_phylo_fumc_fimh_combinations': defaultdict(list),
            'gene_cooccurrence': defaultdict(Counter),
            'high_risk_combinations': [],
        }
        samples = integrated.get('samples', {})
        gc = integrated.get('gene_centric', {})
        sample_genes = defaultdict(list)
        for cat in ('amr_databases', 'virulence_databases',
                    'bacmet_databases', 'plasmid_databases'):
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
            phy = data.get('phylogrouping', {}).get('Clermont_Type', 'ND')
            fumc = data.get('chtyper', {}).get('FumC_Type', 'ND')
            fimh = data.get('chtyper', {}).get('FimH_Type', 'ND')
            cht = data.get('chtyper', {}).get('CH_Type', 'ND')
            patho = data.get('pathotype', 'none')

            if st != 'ND': P['st_distribution'][st] += 1
            if sero != 'ND': P['serotype_distribution'][sero] += 1
            if o_type != 'ND': P['o_type_distribution'][o_type] += 1
            if h_type != 'ND': P['h_type_distribution'][h_type] += 1
            if phy != 'ND': P['phylogroup_distribution'][phy] += 1
            if cht != 'ND': P['ch_type_distribution'][cht] += 1
            if fumc != 'ND': P['fumc_type_distribution'][fumc] += 1
            if fimh != 'ND': P['fimh_type_distribution'][fimh] += 1
            P['pathotype_distribution'][patho] += 1

            if st != 'ND' and o_type != 'ND':
                P['st_o_combinations'][f"{st} - {o_type}"].append(sample)
            if st != 'ND' and h_type != 'ND':
                P['st_h_combinations'][f"{st} - {h_type}"].append(sample)
            if st != 'ND' and sero != 'ND':
                P['st_oh_combinations'][f"{st} - {sero}"].append(sample)
            if st != 'ND' and phy != 'ND':
                P['st_phylogroup_combinations'][f"{st} - {phy}"].append(sample)
            if st != 'ND' and fimh != 'ND':
                P['st_fimh_combinations'][f"{st} - {fimh}"].append(sample)
            if st != 'ND' and fumc != 'ND':
                P['st_fumc_combinations'][f"{st} - {fumc}"].append(sample)
            if st != 'ND' and cht != 'ND':
                P['st_fumc_fimh_combinations'][f"{st} - {cht}"].append(sample)
            if st != 'ND' and sero != 'ND' and phy != 'ND' and cht != 'ND':
                P['st_oh_phylo_fumc_fimh_combinations'][
                    f"{st} - {sero} - {phy} - {cht}"
                ].append(sample)

            genes = sample_genes.get(sample, [])
            for i, g1 in enumerate(genes):
                for g2 in genes[i+1:]:
                    P['gene_cooccurrence'][g1][g2] += 1

            amr_genes = data.get('amr_genes', [])
            vir_genes = data.get('virulence_genes', [])
            crit_amr = [g for g in amr_genes
                        if any(c.lower() in g.lower() for c in self.critical_amr_genes)]
            crit_vir = [g for g in vir_genes
                        if any(c.lower() in g.lower() for c in self.critical_virulence_genes)]
            if crit_amr and crit_vir:
                P['high_risk_combinations'].append({
                    'sample': sample,
                    'st': st, 'serotype': sero, 'phylogroup': phy,
                    'critical_amr_genes': crit_amr,
                    'critical_virulence_genes': crit_vir,
                })
        return P


# =============================================================================
# HTML generator
# =============================================================================
class EcoliHTMLGenerator:
    """Builds the interactive hybrid HTML report."""

    def __init__(self, analyzer: EcoliDataAnalyzer):
        self.data_analyzer = analyzer
        self._current_samples_data = {}

    # -------------------------------------------------------------------------
    # Reusable fragments
    # -------------------------------------------------------------------------
    def _credit_bar(self, color: str, icon: str, title: str, body: str) -> str:
        """Coloured acknowledgement strip for tool-driven tabs."""
        return (f'<div style="background:linear-gradient(135deg,#f8f9fa 0%,#f0f4f8 100%);'
                f'border-left:6px solid {color};margin-bottom:20px;padding:15px;border-radius:8px;">'
                f'<div style="display:flex;align-items:center;gap:12px;flex-wrap:wrap;">'
                f'<span style="font-size:1.4em;">{icon}</span>'
                f'<div><strong style="font-size:1.1em;color:{color};">{title}</strong><br>'
                f'<span style="font-size:.95em;color:#333;">{body}</span></div></div></div>')

    def _alert(self, kind: str, icon: str, body: str) -> str:
        """Standard alert box for info, success, warning, danger."""
        return (f'<div class="alert-box alert-{kind}">'
                f'<i class="fas {icon} fa-2x"></i><div>{body}</div></div>')

    def _collapsible(self, title: str, color: str, inner_html: str) -> str:
        """Collapsed expander for on-demand educational content."""
        return (
            f'<details class="edu-details" style="border-left-color:{color};">'
            f'<summary style="color:{color};">{title}</summary>'
            f'<div class="edu-body">{inner_html}</div>'
            f'</details>')

    def _grouping_dropdown(self, table_id: str) -> str:
        """Grouping dropdown for gene-centric tables."""
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
            <button class="group-btn" onclick="resetGenomeList('{table_id}')">
                <i class="fas fa-undo"></i> Reset</button>
        </div>'''

    def _search_pair(self, table_id: str, placeholder: str) -> str:
        """Two search boxes: filter rows and highlight genome tags."""
        return (f'<input type="text" class="search-box" id="search-{table_id}" '
                f'onkeyup="searchTable(\'{table_id}\',\'search-{table_id}\')" '
                f'placeholder="🔍 {placeholder}">'
                f'<input type="text" class="search-box" id="highlight-{table_id}" '
                f'onkeyup="highlightGenome(\'{table_id}\',\'highlight-{table_id}\')" '
                f'placeholder="🔍 Highlight genome tags containing text">')

    def _combo_table(self, table_id: str, title: str, combo_dict: dict) -> str:
        """A 'Combination | Count | Samples' table."""
        if not combo_dict:
            return ''
        rows = ''
        for combo, samples in sorted(combo_dict.items(), key=lambda x: -len(x[1])):
            tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            rows += (f'<tr><td><strong>{combo}</strong></td>'
                     f'<td>{len(samples)}</td>'
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

    def _gene_family_info(self, color: str, title: str, items: list) -> str:
        """Explanatory box listing biological roles of gene families."""
        html = (f'<div style="margin:10px 0 20px 0;background:#f8f9fa;padding:15px;'
                f'border-radius:8px;font-size:.9em;border-left:4px solid {color};">'
                f'<strong><i class="fas fa-info-circle"></i> {title}</strong><br>')
        for name, desc in items:
            html += f'• <strong>{name}</strong> – {desc}<br>'
        html += '</div>'
        return html

    def _gene_table(self, table_id: str, title: str, gene_list: list,
                    total_samples: int, critical_set: set = None) -> str:
        """Gene-centric table with genome tags and grouping."""
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

    def _multi_db_education(self) -> str:
        """Explains why multiple AMR databases are used side-by-side."""
        return '''
        <div style="background:#e8f5e9;border-radius:8px;padding:18px 22px;">
          <p style="margin:6px 0;font-size:.95em;line-height:1.5;">EcoliTyper screens every genome against <strong>independent AMR databases</strong> because <strong>no single database is comprehensive</strong>. Each has unique strengths, biases, and update cadences.</p>
          <p style="margin:10px 0 6px 0;font-size:.95em;line-height:1.5;"><strong>⚠️ The problem with choosing one database:</strong></p>
          <ul style="margin:6px 0 10px 20px;font-size:.93em;line-height:1.6;">
            <li>Some researchers pick a favourite database (often CARD or ResFinder) and only report hits from that one.</li>
            <li>This is <strong>fast but incomplete</strong> — a gene absent from CARD may still be present in MEGARes or AMRFinderPlus.</li>
            <li>Single-database hits with weak support can be <strong>false positives</strong> that a second database would have flagged.</li>
            <li>The result: <em>biased prevalence estimates</em> and <em>missed resistance signals</em> that matter clinically.</li>
          </ul>
          <p style="margin:10px 0 6px 0;font-size:.95em;line-height:1.5;"><strong>✅ Our approach — report everything, keep provenance:</strong></p>
          <ul style="margin:6px 0 10px 20px;font-size:.93em;line-height:1.6;">
            <li>All hits from all databases are kept <strong>separate and unmerged</strong>. Each box labels the source database for every gene.</li>
            <li>You see the <strong>full picture</strong>. Nothing is silently filtered or cherry-picked.</li>
            <li>Cross-database agreement becomes a <strong>confidence signal</strong>.</li>
          </ul>
        </div>'''

    def _confidence_tiers(self) -> str:
        """Explains how cross-database agreement maps to confidence tiers."""
        return '''
        <div style="background:#e0f2fe;border-radius:8px;padding:16px 20px;">
          <p style="margin:6px 0;font-size:.93em;line-height:1.5;">A gene detected by <strong>multiple databases</strong> is far more likely to be a true positive than a single-database hit.</p>
          <ul style="margin:6px 0 0 20px;font-size:.93em;line-height:1.6;">
            <li><span style="color:#16a34a;font-weight:bold;">🟢 High confidence</span> — gene found in <strong>3 or more databases</strong>.</li>
            <li><span style="color:#f59e0b;font-weight:bold;">🟡 Moderate confidence</span> — gene found in <strong>2 databases</strong>.</li>
            <li><span style="color:#dc2626;font-weight:bold;">🔴 Low confidence / investigate</span> — gene found in <strong>only 1 database</strong>.</li>
          </ul>
        </div>'''

    def _acquired_intrinsic(self) -> str:
        """Distinguishes acquired from intrinsic resistance."""
        return '''
        <div style="background:#f3e8ff;border-radius:8px;padding:16px 20px;">
          <p style="margin:6px 0;font-size:.93em;line-height:1.5;">The AMR story is <strong>more than acquired genes</strong>. We report <strong>both</strong> because they answer different questions:</p>
          <ul style="margin:6px 0 10px 20px;font-size:.93em;line-height:1.6;">
            <li><strong style="color:#7c3aed;">Intrinsic genes</strong> — part of the species baseline. Examples: <em>acrAB-tolC</em>, <em>ampC</em>, chromosomal efflux pumps.</li>
            <li><strong style="color:#e11d48;">Acquired genes</strong> — gained by horizontal transfer. Examples: <em>blaCTX-M</em>, <em>mcr-1</em>, <em>tet</em>, <em>sul</em>.</li>
          </ul>
        </div>'''

    def _genotype_phenotype_caveat(self) -> str:
        """Reminder that gene presence is necessary but not sufficient."""
        return '''
        <div style="background:#fff8e1;border-radius:8px;padding:16px 20px;">
          <p style="margin:6px 0;font-size:.93em;line-height:1.6;">Detecting an AMR gene is <strong>necessary evidence</strong> but not sufficient to declare phenotypic resistance. Several mechanisms can break the link:</p>
          <ul style="margin:6px 0 12px 20px;font-size:.92em;line-height:1.6;">
            <li><strong>Silent or truncated genes</strong> — a partial <em>blaCTX-M</em> fragment produces no enzyme.</li>
            <li><strong>Expression regulation</strong> — <em>ampC</em> is inducible; a mutation in the promoter can silence it.</li>
            <li><strong>Mechanism matters</strong> — efflux (<em>tetA</em>) and ribosomal protection (<em>tetM</em>) produce different phenotypes.</li>
            <li><strong>Naming ambiguity</strong> — ResFinder's <code>blaCTX-M-15_1</code> equals CARD's <code>blaCTX-M-15</code>.</li>
            <li><strong>Dose and route</strong> — low-level efflux may be overcome in vivo.</li>
          </ul>
          <p style="margin:8px 0 0 0;font-size:.92em;background:#ffffff;padding:8px 12px;border-radius:4px;border-left:3px solid #f59e0b;">
            <strong>Clinical bottom line:</strong> <strong>Antimicrobial Susceptibility Testing (AST)</strong> remains the gold standard for treatment decisions.
          </p>
        </div>'''

    def _db_roles_amr(self) -> str:
        """Per-database strengths and weaknesses for AMR interpretation."""
        return '''
        <p style="color:#666;margin-bottom:15px;">Each database is optimised for a different purpose. Knowing what each is good at helps you interpret single-database hits.</p>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:20px;">
          <div style="background:#f8f9fa;padding:16px;border-radius:10px;border-left:4px solid #28a745;">
            <h4 style="color:#28a745;margin:0 0 8px 0;">🟢 CARD</h4>
            <p style="font-size:.9em;margin:0;line-height:1.5;"><strong>Strengths:</strong> Expert-curated, strict SNP-based cutoffs, precise allele assignments.<br><strong>Weaknesses:</strong> Strict thresholds may miss novel or divergent variants.<br><strong>Best for:</strong> Confident allele-level calls.</p>
          </div>
          <div style="background:#f8f9fa;padding:16px;border-radius:10px;border-left:4px solid #17a2b8;">
            <h4 style="color:#17a2b8;margin:0 0 8px 0;">🔵 ResFinder</h4>
            <p style="font-size:.9em;margin:0;line-height:1.5;"><strong>Strengths:</strong> Highly sensitive, frequently updated.<br><strong>Weaknesses:</strong> Appends <code>_1</code> to primary alleles (<code>blaCTX-M-15_1</code> = <code>blaCTX-M-15</code>).<br><strong>Best for:</strong> Broad sensitivity.</p>
          </div>
          <div style="background:#f8f9fa;padding:16px;border-radius:10px;border-left:4px solid #007bff;">
            <h4 style="color:#007bff;margin:0 0 8px 0;">🔷 NCBI AMR</h4>
            <p style="font-size:.9em;margin:0;line-height:1.5;"><strong>Strengths:</strong> Curated by NCBI, tightly linked to the Reference Gene Catalog used by AMRFinderPlus.<br><strong>Weaknesses:</strong> Overlaps substantially with AMRFinderPlus.<br><strong>Best for:</strong> Cross-checking AMRFinderPlus calls.</p>
          </div>
          <div style="background:#f8f9fa;padding:16px;border-radius:10px;border-left:4px solid #fd7e14;">
            <h4 style="color:#fd7e14;margin:0 0 8px 0;">🟠 MEGARes</h4>
            <p style="font-size:.9em;margin:0;line-height:1.5;"><strong>Strengths:</strong> Hierarchy of gene families, includes biocide and metal resistance.<br><strong>Weaknesses:</strong> Broader families can produce class-level hits.<br><strong>Best for:</strong> Environmental co-selection markers.</p>
          </div>
          <div style="background:#f8f9fa;padding:16px;border-radius:10px;border-left:4px solid #6f42c1;">
            <h4 style="color:#6f42c1;margin:0 0 8px 0;">🟣 ARG-ANNOT</h4>
            <p style="font-size:.9em;margin:0;line-height:1.5;"><strong>Strengths:</strong> Historic, well-curated ARG catalogue.<br><strong>Weaknesses:</strong> Updated less frequently.<br><strong>Best for:</strong> Historical comparisons.</p>
          </div>
          <div style="background:#f8f9fa;padding:16px;border-radius:10px;border-left:4px solid #dc3545;">
            <h4 style="color:#dc3545;margin:0 0 8px 0;">🔴 AMRFinderPlus</h4>
            <p style="font-size:.9em;margin:0;line-height:1.5;"><strong>Strengths:</strong> NCBI gold standard, includes <strong>point mutations</strong> alongside gene presence.<br><strong>Weaknesses:</strong> Conservative on some efflux families.<br><strong>Best for:</strong> Clinical-grade gene calls.</p>
          </div>
        </div>'''

    # -------------------------------------------------------------------------
    # CSS
    # -------------------------------------------------------------------------
    def _get_css(self) -> str:
        """Global stylesheet for the hybrid report."""
        return """
<style>
:root{
  --summary-color:#4CAF50;--samples-color:#2196F3;--qc-color:#17a2b8;
  --mlst-color:#FF9800;--serotype-color:#9C27B0;--chtype-color:#009688;
  --phylogroup-color:#795548;--pathotype-color:#E91E63;--amr-color:#F44336;
  --virulence-color:#8BC34A;--plasmid-color:#3F51B5;--bacmet-color:#607D8B;
  --mutation-color:#00BCD4;--mge-color:#16A085;--compare-color:#0d9488;
  --patterns-color:#FF5722;--aiguide-color:#3F51B5;
  --citation-color:#5E35B1;--funding-color:#FFB300;--calltoaction-color:#2E7D32;
  --export-color:#3F51B5;
}
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:'Segoe UI',Tahoma,Geneva,Verdana,sans-serif;line-height:1.6;color:#333;
     background:#f5f5f5;min-width:1200px}
.container{max-width:none;margin:0 auto;padding:20px;width:100%;overflow-x:auto}
.main-header{background:linear-gradient(135deg,#1e3c72 0%,#2a5298 100%);color:white;
             padding:30px;border-radius:15px;margin-bottom:30px;text-align:center}
.main-header h1{font-size:2.6em;margin-bottom:10px;color:white}
.main-header p{font-size:1.1em;opacity:.95}
.metadata-bar{background:rgba(255,255,255,.1);padding:15px;border-radius:10px;margin:20px 0;
              display:flex;justify-content:space-around;flex-wrap:wrap;gap:15px}
.metadata-item{display:flex;align-items:center;gap:8px;font-size:.95em}
.dashboard-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));
                gap:20px;margin-bottom:30px}
.dashboard-card{background:white;padding:20px;border-radius:12px;
                box-shadow:0 5px 20px rgba(0,0,0,.1);text-align:center;
                transition:all .3s ease;cursor:pointer;border-left:5px solid}
.dashboard-card:hover{transform:translateY(-8px);box-shadow:0 15px 30px rgba(0,0,0,.18)}
.card-number{font-size:2.5em;font-weight:bold;margin:10px 0;
             background:linear-gradient(90deg,#1e3c72,#2a5298);
             -webkit-background-clip:text;-webkit-text-fill-color:transparent}
.card-label{font-size:.9em;color:#555;font-weight:600}
.tab-navigation{display:flex;gap:5px;margin-bottom:20px;flex-wrap:wrap;
                background:white;padding:15px;border-radius:12px;position:sticky;
                top:10px;z-index:100;box-shadow:0 5px 20px rgba(0,0,0,.08)}
.tab-button{padding:10px 16px;background:#f5f5f5;border:none;border-radius:8px;
            cursor:pointer;font-weight:600;color:#666;transition:all .3s ease;
            display:flex;align-items:center;gap:6px;font-size:.88em}
.tab-button:hover:not(.active){background:#e3f2fd;color:#1e3c72}
.tab-button.active{background:#1e3c72;color:white}
.tab-content{display:none;background:white;padding:30px;border-radius:15px;
             margin-bottom:30px;width:100%;overflow-x:auto}
.tab-content.active{display:block}
.section-header{color:#2c3e50;margin-bottom:25px;padding-bottom:15px;
                border-bottom:3px solid;font-size:1.8em;display:flex;
                justify-content:space-between;align-items:center}
.data-table{width:100%;border-collapse:collapse;margin:20px 0;font-size:.95em;
            box-shadow:0 2px 10px rgba(0,0,0,.1)}
.data-table th{background:#2c3e50;color:white;padding:12px;text-align:left;
               cursor:pointer;white-space:nowrap;position:sticky;top:0;z-index:6}
.data-table th:hover{background:#1a252f}
.data-table td{padding:10px;border-bottom:1px solid #e0e0e0;vertical-align:top;
               word-wrap:break-word}
.data-table tr:hover{background:#f8f9fa}
.scrollable-table{max-height:640px;overflow:auto;border:1px solid #e0e0e0;
                  border-radius:8px;margin:20px 0;width:100%}
.master-scrollable-container{width:100%;overflow:auto;border:1px solid #e0e0e0;
                              border-radius:8px;margin:20px 0;max-height:640px}
.genome-list{display:flex;flex-wrap:wrap;gap:5px;max-height:260px;overflow:auto;
             padding:6px;background:#f8f9fa;border-radius:5px}
.genome-tag{background:#e0f2f1;color:#00695c;padding:3px 10px;border-radius:12px;
            font-size:.85em;border:1px solid #b2dfdb;white-space:nowrap;margin:2px}
.genome-tag.highlight{background:#ffff99 !important;color:#000 !important;
                      border:1px solid #ffc107}
.genome-group{margin-bottom:10px;width:100%}
.genome-group-header{font-weight:bold;background:#e0e0e0;padding:4px 8px;
                     border-radius:4px;margin:5px 0;font-size:.85em;
                     display:inline-block;word-break:break-word}
.genome-group-tags{display:flex;flex-wrap:wrap;gap:5px;margin-left:10px}
.grouping-controls{background:#f0f7f0;padding:12px;border-radius:8px;margin:15px 0;
                   display:flex;flex-wrap:wrap;gap:10px;align-items:center;
                   border-left:4px solid #006400}
.group-select{padding:8px 14px;border:2px solid #006400;border-radius:6px;
              background:white;color:#333;font-size:.9em;font-weight:600;
              cursor:pointer;min-width:240px}
.group-btn{background:white;border:1px solid #006400;color:#006400;
           padding:6px 12px;border-radius:20px;cursor:pointer;font-size:.85em;
           transition:all .2s}
.group-btn:hover{background:#006400;color:white}
.search-box{width:100%;padding:12px;margin-bottom:10px;border:2px solid #e0e0e0;
            border-radius:8px;font-size:1em}
.search-box:focus{outline:none;border-color:#1e3c72;box-shadow:0 0 0 3px rgba(30,60,114,.1)}
.alert-box{padding:20px;border-radius:10px;margin:20px 0;display:flex;
           align-items:flex-start;gap:20px;border-left:5px solid}
.alert-info{background:#d1ecf1;color:#0c5460;border-left-color:#17a2b8}
.alert-success{background:#d4edda;color:#155724;border-left-color:#28a745}
.alert-warning{background:#fff3cd;color:#856404;border-left-color:#ffc107}
.alert-danger{background:#f8d7da;color:#721c24;border-left-color:#dc3545}
.action-buttons{display:flex;gap:10px;margin:20px 0;flex-wrap:wrap;align-items:center}
.action-btn{padding:9px 18px;border:none;border-radius:8px;cursor:pointer;
            font-weight:600;display:inline-flex;align-items:center;gap:8px;
            transition:all .3s ease;text-decoration:none;font-size:.9em}
.action-btn:hover{transform:translateY(-2px);box-shadow:0 5px 15px rgba(0,0,0,.15)}
.btn-primary{background:#1e3c72;color:white}
.btn-success{background:#28a745;color:white}
.btn-danger{background:#dc3545;color:white}
.btn-warning{background:#ffc107;color:black}
.btn-info{background:#17a2b8;color:white}
.btn-secondary{background:#6c757d;color:white}
.btn-light{background:#f8f9fa;color:#212529;border:1px solid #dee2e6}
.database-section{margin:30px 0;padding:22px;border-radius:12px;background:#f8f9fa;
                  box-shadow:0 3px 15px rgba(0,0,0,.08)}
.print-section-btn{background:#1e3c72;color:white;border:none;border-radius:5px;
                   padding:8px 15px;cursor:pointer;display:flex;align-items:center;
                   gap:5px;font-size:.9em}
.footer{text-align:center;padding:30px;color:white;margin-top:40px;
        border-radius:15px;background:linear-gradient(135deg,#2c3e50 0%,#34495e 100%)}
.stat-card{color:white;padding:18px;border-radius:10px;text-align:center;
           box-shadow:0 4px 15px rgba(0,0,0,.15);transition:transform .2s}
.stat-card:hover{transform:translateY(-3px)}
.stat-card .stat-value{font-size:1.9em;font-weight:bold;margin-bottom:4px}
.stat-card .stat-label{font-size:.85em;opacity:.95;text-transform:uppercase;letter-spacing:.5px}
.stats-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));
            gap:15px;margin:20px 0}
.badge{display:inline-block;padding:5px 15px;border-radius:20px;font-size:.85em;
       font-weight:600;margin:2px}
.info-text{background:#f8f9fa;padding:15px;border-radius:8px;margin:15px 0;
           border-left:4px solid #17a2b8}
.col-toggle-wrap{position:relative}
.col-panel{position:absolute;right:0;top:110%;background:#fff;border:1px solid #dee2e6;
           border-radius:10px;padding:10px 14px;box-shadow:0 10px 30px rgba(0,0,0,.15);
           z-index:30;display:none;max-height:320px;overflow-y:auto;min-width:200px}
.col-panel.open{display:block}
.col-panel label{display:block;font-size:13px;color:#1e3c72;padding:4px 0;cursor:pointer}
.sel-counter{font-size:12px;font-weight:700;color:#1e3c72;background:#e3f2fd;
             padding:6px 12px;border-radius:20px}
.col-check{width:34px;text-align:center}
.row-selected td{background:#e3f2fd !important}
.mode-toggle{display:flex;gap:8px;margin-bottom:18px}
.mode-btn{padding:10px 22px;border:2px solid #cbd5e1;background:white;color:#475569;
          border-radius:8px;font-weight:600;cursor:pointer;font-size:.95em;
          display:flex;align-items:center;gap:8px;transition:all .2s}
.mode-btn:hover{border-color:#1e3c72;color:#1e3c72}
.mode-btn.active{background:#1e3c72;border-color:#1e3c72;color:white}
.compare-mode{display:none}
.compare-mode.active{display:block}
.compare-panel{background:#f8f9fa;border-radius:12px;padding:20px;margin-bottom:24px}
.compare-select-row{display:flex;gap:12px;align-items:flex-end;flex-wrap:wrap;margin-bottom:16px}
.compare-select-row label{display:flex;flex-direction:column;gap:4px}
.compare-select-row label span{font-size:.85em;font-weight:600;color:#475569}
.compare-select-row select{padding:10px 14px;border:2px solid #cbd5e1;border-radius:8px;
                           background:white;font-size:.95em;min-width:220px}
.compare-visuals-grid{display:grid;grid-template-columns:180px 1fr;gap:16px;
                      margin:20px 0;align-items:center}
.visual-card{background:white;border-radius:12px;padding:20px;
             box-shadow:0 2px 8px rgba(0,0,0,.06);border:1px solid #e2e8f0}
.gauge-card{display:flex;justify-content:center;align-items:center;padding:12px}
.similarity-gauge{width:160px;height:160px;display:block}
.metrics-card{display:flex;flex-direction:column;gap:18px}
.metric-block{display:flex;flex-direction:column;gap:6px}
.metric-label{font-size:.85em;font-weight:700;color:#475569;text-transform:uppercase;letter-spacing:.5px}
.metric-bar-track{background:#f1f5f9;height:22px;border-radius:11px;overflow:hidden}
.metric-bar-fill{height:100%;border-radius:11px;transition:width .8s ease;min-width:4px}
.metric-value{font-size:.88em;color:#334155;font-weight:600}
.typing-pill{display:inline-block;padding:5px 12px;border-radius:12px;
             font-size:.85em;font-weight:600;border:1.5px solid;
             font-family:'Consolas',monospace}
.gene-proportion-bar{display:flex;height:12px;border-radius:6px;overflow:hidden;
                     margin-bottom:14px;background:#f1f5f9}
.gene-prop-segment{height:100%;transition:width .6s ease;min-width:2px}
.gene-prop-segment.shared{background:#10b981}
.gene-prop-segment.only-a{background:#3b82f6}
.gene-prop-segment.only-b{background:#ec4899}
.cluster-card{background:white;border-radius:10px;padding:16px 20px;
              margin-bottom:12px;border-left:5px solid;
              box-shadow:0 2px 8px rgba(0,0,0,.06)}
.cluster-header{display:flex;align-items:center;gap:14px;margin-bottom:12px;flex-wrap:wrap}
.cluster-badge{color:white;padding:4px 14px;border-radius:14px;font-size:.85em;font-weight:700}
.cluster-stats{font-size:.9em;color:#64748b;font-weight:600}
.cluster-samples{display:flex;flex-wrap:wrap;gap:6px}
.cluster-sample-pill{background:#f1f5f9;color:#1e293b;padding:4px 12px;border-radius:12px;
                     font-size:.85em;font-weight:600;font-family:'Consolas',monospace;
                     border:1px solid #cbd5e1}
.heatmap-wrapper{background:white;border-radius:12px;padding:20px;margin-top:24px;
                 box-shadow:0 2px 8px rgba(0,0,0,.06);overflow-x:auto}
.heatmap-grid{display:grid;gap:2px;margin-top:16px;padding:4px}
.heatmap-cell{width:32px;height:32px;display:flex;align-items:center;justify-content:center;
              font-size:.7em;font-weight:700;border-radius:3px}
.heatmap-col-header{background:transparent;color:#64748b;font-size:.75em;
                    writing-mode:vertical-rl;text-orientation:mixed;height:70px;
                    width:32px;font-family:'Consolas',monospace}
.heatmap-row-header{background:transparent;color:#64748b;font-size:.75em;
                    justify-content:flex-start;padding-left:6px;
                    font-family:'Consolas',monospace}
.heatmap-corner{background:transparent}
.compare-section{background:white;border-radius:12px;padding:20px;margin-bottom:20px;
                 box-shadow:0 2px 8px rgba(0,0,0,.06)}
.compare-section h3{margin-bottom:14px;color:#1e293b;font-size:1.15em;
                    padding-bottom:8px;border-bottom:2px solid #e2e8f0}
.compare-table{display:flex;flex-direction:column;gap:2px}
.compare-row{display:grid;grid-template-columns:180px 1fr 1fr;gap:12px;padding:10px 14px;
             border-radius:6px;align-items:center;font-size:.92em}
.compare-row.header{background:#f1f5f9;font-weight:700;border-bottom:2px solid #cbd5e1}
.compare-row.diff{background:#fef2f2;border-left:3px solid #dc2626}
.compare-row.match{background:#f0fdf4;border-left:3px solid #16a34a}
.compare-row .row-label{font-weight:600;color:#475569;font-size:.88em;
                        text-transform:uppercase;letter-spacing:.3px}
.compare-row .row-val{color:#1e293b;word-break:break-word}
.compare-gene-block{margin-bottom:24px;padding-bottom:20px;border-bottom:1px dashed #e2e8f0}
.compare-gene-block:last-child{border-bottom:none;margin-bottom:0;padding-bottom:0}
.gene-count{font-size:.82em;font-weight:400;color:#64748b;background:#f1f5f9;
            padding:3px 10px;border-radius:12px}
.gene-columns{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:14px}
.gene-col{background:#f8fafc;border-radius:8px;padding:12px;border:1px solid #e2e8f0}
.gene-col-head{font-size:.85em;font-weight:700;color:#475569;margin-bottom:8px;
               text-transform:uppercase;letter-spacing:.3px}
.compare-gene{display:inline-block;font-size:.82em;padding:3px 10px;border-radius:12px;
              margin:2px;font-family:'Consolas',monospace}
.compare-gene.shared{background:#dcfce7;color:#14532d;border:1px solid #86efac}
.compare-gene.only-a{background:#dbeafe;color:#1e3a8a;border:1px solid #93c5fd}
.compare-gene.only-b{background:#fce7f3;color:#831843;border:1px solid #f9a8d4}
.compare-full-details{background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
                      padding:12px 16px;margin-bottom:10px}
.compare-full-details[open]{background:#fff;border-color:#1e3c72}
.compare-full-details summary{cursor:pointer;font-weight:700;color:#1e3c72;
                              list-style:none;padding:4px 0}
.compare-full-profile{display:grid;grid-template-columns:1fr 1fr;gap:16px;
                      margin-top:14px;padding-top:14px;border-top:1px dashed #e2e8f0}
.profile-col{background:#f8fafc;border-radius:10px;padding:14px;border:1px solid #e2e8f0}
.profile-head{display:flex;justify-content:space-between;align-items:center;gap:8px;
              margin-bottom:10px;padding-bottom:8px;border-bottom:2px solid #e2e8f0}
.profile-sample-label{font-weight:700;color:#1e293b;font-size:.92em;word-break:break-word}
.profile-sample-count{font-size:.78em;font-weight:600;color:#1e3c72;background:#dbeafe;
                      padding:2px 10px;border-radius:12px;white-space:nowrap}
.empty-note{font-size:.85em;color:#94a3b8;font-style:italic;padding:6px 0;display:inline-block}
.compare-warning{padding:20px;background:#fef2f2;color:#7f1d1d;border-radius:8px;
                 border-left:4px solid #dc2626}
.compare-verdict-box{padding:16px 20px;border-radius:10px;margin-bottom:16px;
                     display:flex;gap:12px;align-items:flex-start;font-size:.95em}
.compare-verdict-box.danger{background:#fef2f2;border-left:5px solid #dc2626;color:#7f1d1d}
.compare-verdict-box.warning{background:#fffbeb;border-left:5px solid #f59e0b;color:#78350f}
.compare-verdict-box.info{background:#eff6ff;border-left:5px solid #2563eb;color:#1e3a8a}
.cluster-controls{display:flex;gap:14px;align-items:flex-end;flex-wrap:wrap;margin-bottom:20px}
.cluster-controls select{padding:10px 14px;border:2px solid #cbd5e1;border-radius:8px;
                         background:white;font-size:.95em;min-width:220px}
.no-alerts{padding:20px;background:#eff6ff;color:#1e3a8a;border-radius:8px;
           border-left:4px solid #2563eb}
.isolate-box{border:1px solid #ddd;border-radius:12px;margin-bottom:20px;padding:20px;
             background:#fafafa;box-shadow:0 2px 8px rgba(0,0,0,.06);transition:box-shadow .2s}
.isolate-box:hover{box-shadow:0 4px 14px rgba(0,0,0,.10)}
.isolate-box .sample-header{display:flex;align-items:center;gap:15px;flex-wrap:wrap}
.isolate-box .sample-header h3{font-size:1.4em;margin:0}
.isolate-box .sample-header .total-badge{background:#1e3c72;color:white;padding:4px 16px;
                                         border-radius:20px;font-weight:bold;font-size:.9em}
.isolate-box .sample-header .typing-info{display:flex;flex-wrap:wrap;gap:8px}
.isolate-box .toggle-btn{background:#1e3c72;color:white;border:none;padding:8px 16px;
                         border-radius:8px;cursor:pointer;font-weight:600;font-size:.85em;
                         display:inline-flex;align-items:center;gap:6px;margin-left:auto;
                         transition:background .2s}
.isolate-box .toggle-btn:hover{background:#2a5298}
.isolate-box .box-details{margin-top:15px}
.database-table-wrapper{margin:15px 0;overflow-x:auto;border:1px solid #e0e0e0;
                        border-radius:8px}
.database-table-wrapper table{width:100%;border-collapse:collapse;font-size:.85em;
                              min-width:800px}
.database-table-wrapper table th{background:#2c3e50;color:white;padding:8px 12px;
                                 text-align:left;white-space:nowrap;position:sticky;
                                 top:0;z-index:6}
.database-table-wrapper table td{padding:8px 12px;border-bottom:1px solid #e0e0e0;
                                 white-space:nowrap}
.database-table-wrapper table tr:hover{background:#f1f1f1}
.db-title{font-weight:bold;color:#1e3c72;margin:10px 0 5px 0;font-size:1.05em;
          border-left:4px solid #1e3c72;padding-left:10px}
.filter-controls{display:flex;flex-wrap:wrap;gap:10px;align-items:center;
                 background:#f8f9fa;padding:15px;border-radius:8px;margin-bottom:20px}
.filter-controls select{padding:10px;border-radius:8px;border:2px solid #ddd;
                        background:white;min-width:150px}
.results-counter{font-size:.9em;color:#555;font-weight:600;padding:8px 12px;
                 background:#f0f0f0;border-radius:6px;white-space:nowrap}
.typing-badge{display:inline-block;padding:3px 10px;border-radius:12px;font-size:.8em;
              font-weight:600;background:#e0e0e0;color:#333;border:1px solid #ccc}
.typing-badge.badge-pathotype{background:#E91E63;color:white;border-color:#E91E63}
.typing-badge.badge-pathotype-none{background:#e0e0e0;color:#555;border-color:#ccc}
.edu-details{background:#f8f9fa;border:1px solid #e0e0e0;border-left:4px solid #1e3c72;
             border-radius:8px;margin:12px 0;overflow:hidden}
.edu-details summary{cursor:pointer;padding:14px 18px;font-weight:600;
                     background:#eef4fb;list-style:none;display:flex;
                     align-items:center;gap:8px;color:#1e3c72}
.edu-details summary:hover{background:#e0eaf7}
.edu-details summary::-webkit-details-marker{display:none}
.edu-details summary::before{content:'▸';display:inline-block;margin-right:6px;
                              transition:transform .2s}
.edu-details[open] summary::before{content:'▾'}
.edu-details[open]{background:#fff}
.edu-details .edu-body{padding:16px 20px;border-top:1px solid #e0e0e0}
@media (max-width:900px){
  .compare-full-profile{grid-template-columns:1fr}
  .compare-row{grid-template-columns:1fr}
  .compare-visuals-grid{grid-template-columns:1fr}
}
@media print{
  .tab-navigation,.dashboard-grid,.search-box,.action-buttons,.grouping-controls,
  .print-section-btn,.col-panel,.filter-controls,.toggle-btn{display:none!important}
  .tab-content.active{display:block}
  .isolate-box .box-details{display:block!important}
  .edu-details[open] .edu-body{display:block!important}
}
</style>"""

    # -------------------------------------------------------------------------
    # JavaScript
    # -------------------------------------------------------------------------
    def _get_js(self, typing_json: str) -> str:
        """Client-side scripting for tabs, search, sort, grouping, exports, compare, boxes."""
        js = """
<script>
var sampleTyping = __TYPING__;
var originalGenomeLists = {};
window.ECOLITYPER_BOX_DATA = window.ECOLITYPER_BOX_DATA || {};
window.ECOLITYPER_COMPARE_DATA = window.ECOLITYPER_COMPARE_DATA || {};

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
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;')
                  .replace(/>/g,'&gt;').replace(/"/g,'&quot;');
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
  var st = info.ST || "ND", sero = info.Serotype || "ND";
  var phy = info.Phylogroup || "ND", cht = info.CHtype || "ND";
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
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(a.href);
}

function printSection(id) {
  var content = document.getElementById(id); if (!content) return;
  content.querySelectorAll('.isolate-box').forEach(function(box) {
    var tabId = box.getAttribute('data-tab');
    var sample = box.getAttribute('data-sample');
    var details = box.querySelector('.box-details');
    if (details && !details.innerHTML.trim()) {
      details.innerHTML = buildBoxDetails(tabId, sample);
    }
  });
  var w = window.open('', '_blank');
  var style = document.querySelector('style');
  w.document.write('<html><head><title>Print</title>');
  if (style) w.document.write('<style>' + style.textContent + '</style>');
  w.document.write('</head><body>' + content.innerHTML + '</body></html>');
  w.document.close();
  w.print();
}

function buildBoxDetails(tabId, sample) {
  if (tabId === 'mutation') return buildMutationTable(sample);
  if (tabId === 'mge')      return buildMgeTable(sample);
  return buildBoxContent(tabId, sample);
}

function toggleBoxDetails(tabId, sample) {
  var box = document.querySelector('#' + tabId + '-tab .isolate-box[data-sample="' + CSS.escape(sample) + '"]');
  if (!box) return;
  var details = box.querySelector('.box-details');
  var btn = box.querySelector('.toggle-btn');
  var isOpen = details.style.display === 'block';
  if (isOpen) {
    details.style.display = 'none';
    if (btn) btn.innerHTML = '<i class="fas fa-chevron-down"></i> Show Details';
    return;
  }
  if (!details.innerHTML.trim()) {
    details.innerHTML = buildBoxDetails(tabId, sample);
  }
  details.style.display = 'block';
  if (btn) btn.innerHTML = '<i class="fas fa-chevron-up"></i> Hide Details';
}

function buildBoxContent(tabId, sample) {
  var root = window.ECOLITYPER_BOX_DATA || {};
  var byTab = root[tabId] || {};
  var entry = byTab[sample];
  if (!entry) return '<p>No data available</p>';
  var html = '';
  if (entry.amrfinder && entry.amrfinder.length) {
    html += buildDBTable('AMRfinder', entry.amrfinder, 'amrfinder');
  }
  if (entry.abricate) {
    for (var db in entry.abricate) {
      var genes = entry.abricate[db];
      if (genes && genes.length) {
        html += buildDBTable(db.toUpperCase(), genes, db);
      }
    }
  }
  return html || '<p>No gene details</p>';
}

function buildDBTable(dbName, genes, dbKey) {
  if (!genes || !genes.length) return '';
  var keys = Object.keys(genes[0]);
  var priority = ['gene', 'product', '%coverage', 'coverage_percent',
                  'coverage', '%identity', 'identity_percent', 'identity',
                  'accession', 'contig', 'start', 'stop', 'class', 'subclass',
                  'scope', 'resistance', 'database'];
  var ordered = [];
  for (var i = 0; i < priority.length; i++) {
    if (keys.indexOf(priority[i]) !== -1) ordered.push(priority[i]);
  }
  for (var k = 0; k < keys.length; k++) {
    if (ordered.indexOf(keys[k]) === -1) ordered.push(keys[k]);
  }
  var html = '<div class="database-table-wrapper" data-db="' + dbKey + '">';
  html += '<div class="db-title">' + escapeHtml(dbName) + '</div>';
  html += '<table><thead><tr>';
  for (var h = 0; h < ordered.length; h++) {
    var display = ordered[h].replace(/_/g, ' ').replace(/\\b\\w/g, function(x){ return x.toUpperCase(); });
    html += '<th>' + escapeHtml(display) + '</th>';
  }
  html += '</tr></thead><tbody>';
  for (var g = 0; g < genes.length; g++) {
    html += '<tr>';
    for (var c = 0; c < ordered.length; c++) {
      var v = (genes[g][ordered[c]] == null) ? '' : genes[g][ordered[c]];
      html += '<td>' + escapeHtml(v) + '</td>';
    }
    html += '</tr>';
  }
  html += '</tbody></table></div>';
  return html;
}

function buildMutationTable(sample) {
  var root = window.ECOLITYPER_BOX_DATA || {};
  var byTab = root['mutation'] || {};
  var entry = byTab[sample];
  if (!entry || !entry.mutations || !entry.mutations.length) return '<p>No mutations</p>';
  var cols = ['gene', 'mutation', 'class', 'subclass', 'contig', 'start',
              'stop', 'strand', 'coverage', 'identity', 'accession'];
  var html = '<div class="database-table-wrapper" data-db="mutations">';
  html += '<div class="db-title">Mutations</div><table><thead><tr>';
  for (var c = 0; c < cols.length; c++) {
    var display = cols[c].replace(/_/g, ' ').replace(/\\b\\w/g, function(x){ return x.toUpperCase(); });
    html += '<th>' + escapeHtml(display) + '</th>';
  }
  html += '</tr></thead><tbody>';
  for (var i = 0; i < entry.mutations.length; i++) {
    html += '<tr>';
    for (var j = 0; j < cols.length; j++) {
      var v = (entry.mutations[i][cols[j]] == null) ? '' : entry.mutations[i][cols[j]];
      html += '<td>' + escapeHtml(v) + '</td>';
    }
    html += '</tr>';
  }
  html += '</tbody></table></div>';
  return html;
}

function buildMgeTable(sample) {
  var root = window.ECOLITYPER_BOX_DATA || {};
  var byTab = root['mge'] || {};
  var entry = byTab[sample];
  if (!entry || !entry.mge || !entry.mge.length) return '<p>No MGE hits</p>';
  var cols = ['rank', 'mobileog_id', 'gene', 'contig', 'start', 'stop', 'strand',
              'major_category', 'minor_category', 'percent_identity',
              'align_over_qlen', 'coverage_percent'];
  var headers = ['Rank', 'mobileOG ID', 'Gene', 'Contig', 'Start', 'Stop', 'Strand',
                 'Major Category', 'Minor Category', '% Identity',
                 'Align / Qlen', 'Coverage %'];
  var html = '<div class="database-table-wrapper" data-db="mge">';
  html += '<div class="db-title">MGE (mobileOG-db)</div>';
  html += '<table><thead><tr>';
  for (var h = 0; h < headers.length; h++) {
    html += '<th>' + escapeHtml(headers[h]) + '</th>';
  }
  html += '</tr></thead><tbody>';
  for (var i = 0; i < entry.mge.length; i++) {
    html += '<tr>';
    for (var j = 0; j < cols.length; j++) {
      var v = (entry.mge[i][cols[j]] == null) ? '' : entry.mge[i][cols[j]];
      html += '<td>' + escapeHtml(v) + '</td>';
    }
    html += '</tr>';
  }
  html += '</tbody></table></div>';
  return html;
}

function expandAllBoxes(tabId) {
  var boxes = document.querySelectorAll('#' + tabId + '-tab .isolate-box');
  var visible = [];
  boxes.forEach(function(b) { if (b.style.display !== 'none') visible.push(b); });
  if (visible.length > 100) {
    if (!confirm('You have ' + visible.length + ' visible samples. Expanding all may take a moment. Continue?')) return;
  }
  var i = 0, batch = 20;
  function step() {
    var end = Math.min(i + batch, visible.length);
    for (; i < end; i++) {
      var box = visible[i];
      var sample = box.getAttribute('data-sample');
      var details = box.querySelector('.box-details');
      var btn = box.querySelector('.toggle-btn');
      if (details && details.style.display !== 'block') {
        if (!details.innerHTML.trim()) {
          details.innerHTML = buildBoxDetails(tabId, sample);
        }
        details.style.display = 'block';
        if (btn) btn.innerHTML = '<i class="fas fa-chevron-up"></i> Hide Details';
      }
    }
    if (i < visible.length) setTimeout(step, 30);
  }
  step();
}

function collapseAllBoxes(tabId) {
  document.querySelectorAll('#' + tabId + '-tab .isolate-box').forEach(function(box) {
    var details = box.querySelector('.box-details');
    if (details) details.style.display = 'none';
    var btn = box.querySelector('.toggle-btn');
    if (btn) btn.innerHTML = '<i class="fas fa-chevron-down"></i> Show Details';
  });
}

function filterBoxes(tabId) {
  var searchEl = document.getElementById('search-' + tabId);
  var dbEl = document.getElementById('dbFilter-' + tabId);
  var search = (searchEl ? searchEl.value : '').toUpperCase();
  var dbFilter = dbEl ? dbEl.value : 'all';
  var boxes = document.querySelectorAll('#' + tabId + '-tab .isolate-box');
  var visible = 0;
  boxes.forEach(function(box) {
    var sample = box.getAttribute('data-sample') || '';
    var show = (!search || sample.toUpperCase().indexOf(search) !== -1);
    box.style.display = show ? '' : 'none';
    if (show) visible++;
    if (show && dbEl) {
      var wrappers = box.querySelectorAll('.database-table-wrapper');
      wrappers.forEach(function(w) {
        var dbName = w.getAttribute('data-db') || '';
        var keep = (dbFilter === 'all' || dbName === dbFilter);
        w.style.display = keep ? '' : 'none';
      });
    }
  });
  var counter = document.getElementById('counter-' + tabId);
  if (counter) counter.textContent = visible + ' shown';
}

function resetBoxFilters(tabId) {
  var s = document.getElementById('search-' + tabId);
  var d = document.getElementById('dbFilter-' + tabId);
  if (s) s.value = '';
  if (d) d.value = 'all';
  filterBoxes(tabId);
}

function setCompareMode(mode) {
  document.querySelectorAll('.mode-btn').forEach(function(b){
    b.classList.toggle('active', b.dataset.mode === mode);
  });
  document.querySelectorAll('.compare-mode').forEach(function(m){ m.classList.remove('active'); });
  var el = document.getElementById('mode-' + mode);
  if (el) el.classList.add('active');
  if (mode === 'cluster') runCluster();
}

function similarityColor(pct) {
  if (pct >= 95) return '#059669';
  if (pct >= 90) return '#10b981';
  if (pct >= 80) return '#84cc16';
  if (pct >= 70) return '#eab308';
  if (pct >= 50) return '#f59e0b';
  if (pct >= 30) return '#f97316';
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
    var sa = new Set(A.genes[cat] || []);
    var sb = new Set(B.genes[cat] || []);
    var u = new Set(); sa.forEach(function(g){ u.add(g); }); sb.forEach(function(g){ u.add(g); });
    if (u.size === 0) return;
    var shared = 0;
    u.forEach(function(g){ if (sa.has(g) && sb.has(g)) shared++; });
    scores.push((shared / u.size) * 100);
  });
  var genePct = scores.length ? scores.reduce(function(a, b){ return a + b; }, 0) / scores.length : 0;
  return {
    typing: Math.round(typingPct),
    genes: Math.round(genePct),
    overall: Math.round(typingPct * 0.5 + genePct * 0.5),
    typingMatches: matches,
    typingTotal: tkeys.length,
  };
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
        '<div class="metric-bar-track"><div class="metric-bar-fill" style="width:' + sim.typing +
        '%;background:' + similarityColor(sim.typing) + ';"></div></div>' +
        '<div class="metric-value">' + sim.typingMatches + '/' + sim.typingTotal +
        ' fields (' + sim.typing + '%)</div></div>' +
        '<div class="metric-block"><div class="metric-label">Gene content match</div>' +
        '<div class="metric-bar-track"><div class="metric-bar-fill" style="width:' + sim.genes +
        '%;background:' + similarityColor(sim.genes) + ';"></div></div>' +
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
    geneHtml += '<div class="compare-gene-block"><h4>' + cat +
      ' <span class="gene-count">' + shared.length + ' shared · ' + onlyA.length + ' only A · ' +
      onlyB.length + ' only B (' + sPct.toFixed(0) + '% match)</span></h4>' +
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
    vHtml = '<div class="compare-verdict-box danger"><i class="fas fa-exclamation-triangle"></i> <strong>Near-identical (' + sim.overall + '%)</strong> — likely the same strain. Investigate as transmission pair or duplicate.</div>';
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
  var data = window.ECOLITYPER_COMPARE_DATA || {};
  var samples = Object.keys(data).sort();
  if (samples.length < 2) {
    document.getElementById('cluster-heatmap').innerHTML = '<div class="compare-warning">Need at least two samples.</div>';
    return;
  }
  var thresholdEl = document.getElementById('cluster-threshold');
  var threshold = thresholdEl ? parseInt(thresholdEl.value) : 90;
  var matrix = {};
  samples.forEach(function(s){ matrix[s] = {}; });
  samples.forEach(function(s1) {
    samples.forEach(function(s2) {
      matrix[s1][s2] = s1 === s2 ? 100 : computeSimilarity(data[s1], data[s2]).overall;
    });
  });
  var parent = {};
  samples.forEach(function(s){ parent[s] = s; });
  function find(x) { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }
  function union(a, b) { var ra = find(a), rb = find(b); if (ra !== rb) parent[ra] = rb; }
  for (var i = 0; i < samples.length; i++)
    for (var j = i + 1; j < samples.length; j++)
      if (matrix[samples[i]][samples[j]] >= threshold) union(samples[i], samples[j]);
  var groups = {};
  samples.forEach(function(s){ var r = find(s); (groups[r] = groups[r] || []).push(s); });
  var clusters = Object.keys(groups).map(function(k){ return groups[k]; })
                       .filter(function(g){ return g.length >= 2; })
                       .sort(function(a, b){ return b.length - a.length; });
  var clusterHtml = '';
  if (clusters.length === 0) {
    clusterHtml = '<div class="no-alerts">No clusters detected at ≥' + threshold + '%.</div>';
  } else {
    clusters.forEach(function(cluster, idx) {
      var sum = 0, count = 0;
      for (var i = 0; i < cluster.length; i++)
        for (var j = i + 1; j < cluster.length; j++) { sum += matrix[cluster[i]][cluster[j]]; count++; }
      var mean = count ? (sum / count).toFixed(1) : '100';
      var color = similarityColor(parseFloat(mean));
      clusterHtml += '<div class="cluster-card" style="border-left-color:' + color + ';">' +
        '<div class="cluster-header"><span class="cluster-badge" style="background:' + color + ';">Cluster ' + (idx + 1) + '</span>' +
        '<span class="cluster-stats">' + cluster.length + ' samples · mean ' + mean + '% similar</span></div>' +
        '<div class="cluster-samples">' + cluster.map(function(s){ return '<span class="cluster-sample-pill">' + escapeHtml(s) + '</span>'; }).join('') + '</div></div>';
    });
  }
  var heat = '<div class="heatmap-wrapper"><h3>🔷 Pairwise Similarity Matrix</h3>' +
    '<div class="heatmap-grid" style="grid-template-columns:180px repeat(' + samples.length + ',32px);">';
  heat += '<div class="heatmap-cell heatmap-corner"></div>';
  samples.forEach(function(s){ heat += '<div class="heatmap-cell heatmap-col-header">' + escapeHtml(s.slice(-8)) + '</div>'; });
  samples.forEach(function(s1) {
    heat += '<div class="heatmap-cell heatmap-row-header">' + escapeHtml(s1) + '</div>';
    samples.forEach(function(s2) {
      var pct = matrix[s1][s2];
      var color = (s1 === s2) ? '#1e293b' : similarityColor(pct);
      heat += '<div class="heatmap-cell" style="background:' + color + ';color:white;">' + (s1 === s2 ? '—' : Math.round(pct)) + '</div>';
    });
  });
  heat += '</div></div>';
  document.getElementById('cluster-clusters').innerHTML = clusterHtml;
  document.getElementById('cluster-heatmap').innerHTML = heat;
  var s = document.getElementById('cluster-summary');
  if (s) s.textContent = clusters.length + ' cluster(s) at ≥' + threshold + '%';
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
});
</script>"""
        return js.replace('__TYPING__', typing_json)

    # =====================================================================
    # Section: Sample-centric box generator (shared)
    # =====================================================================
    def _generate_sample_centric_boxes(self, integrated_data: Dict, tab_id: str,
                                       title: str, db_list: List[str]) -> str:
        """Build a sample-centric tab with lazy-loaded isolate boxes."""
        samples_data = integrated_data.get('samples', {})
        amr_details = integrated_data.get('amrfinder_details', {})
        abricate_details = integrated_data.get('abricate_details', {})

        relevant = []
        for sample in samples_data:
            has_data = False
            if 'amrfinder' in db_list and amr_details.get(sample):
                has_data = True
            for db in db_list:
                if db != 'amrfinder' and abricate_details.get(sample, {}).get(db):
                    has_data = True
            if has_data:
                relevant.append(sample)
        relevant.sort()

        if not relevant:
            return self._alert('warning', 'fa-exclamation-circle',
                f'<h3>No {title} Data Available</h3>'
                f'<p>No TSV data found for {title}. Ensure the corresponding '
                f'TSV files are present in the input directory.</p>')

        blob = {}
        for sample in relevant:
            entry = {'amrfinder': [], 'abricate': {}}
            if 'amrfinder' in db_list:
                entry['amrfinder'] = amr_details.get(sample, [])
            for db in db_list:
                if db != 'amrfinder':
                    entry['abricate'][db] = abricate_details.get(sample, {}).get(db, [])
            blob[sample] = entry
        blob_json = json.dumps(blob, default=str).replace('</', '<\\/')

        db_options = '<option value="all">All Databases</option>'
        for db in db_list:
            display = 'AMRfinder' if db == 'amrfinder' else db.upper()
            db_options += f'<option value="{db}">{display}</option>'

        html = ''
        html += self._alert('info', 'fa-info-circle',
            f'<h3>🧬 {title} – Interactive Isolate Boxes</h3>'
            f'<p>Each box represents one isolate. Click <strong>Show Details</strong> '
            f'to load its gene tables on demand. The report opens instantly even '
            f'with thousands of samples because tables are only built when requested.</p>'
            f'<p>Use the search box and the database dropdown below to filter boxes '
            f'and gene tables live.</p>')

        html += f'''
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-{tab_id}"
                   onkeyup="filterBoxes('{tab_id}')"
                   placeholder="🔍 Search sample by name..." style="max-width:320px;">
            <select id="dbFilter-{tab_id}" onchange="filterBoxes('{tab_id}')">{db_options}</select>
            <button class="action-btn btn-info" onclick="expandAllBoxes('{tab_id}')">
                <i class="fas fa-expand-alt"></i> Expand All (visible)</button>
            <button class="action-btn btn-light" onclick="collapseAllBoxes('{tab_id}')">
                <i class="fas fa-compress-alt"></i> Collapse All</button>
            <button class="action-btn btn-success" onclick="resetBoxFilters('{tab_id}')">
                <i class="fas fa-sync"></i> Clear Filters</button>
            <span class="results-counter" id="counter-{tab_id}">{len(relevant)} shown</span>
        </div>
        <div id="box-container-{tab_id}">'''

        for sample in relevant:
            sd = samples_data.get(sample, {})
            mlst = sd.get('mlst', {}).get('ST', 'ND')
            sero = sd.get('serotype', {}).get('Serotype', 'ND')
            phy = sd.get('phylogrouping', {}).get('Clermont_Type', 'ND')
            cht = sd.get('chtyper', {}).get('CH_Type', 'ND')
            patho = sd.get('pathotype', 'none')

            total = 0
            if 'amrfinder' in db_list:
                total += len(amr_details.get(sample, []))
            for db in db_list:
                if db != 'amrfinder':
                    total += len(abricate_details.get(sample, {}).get(db, []))

            patho_class = 'badge-pathotype' if patho != 'none' else 'badge-pathotype-none'

            html += f'''
            <div class="isolate-box" data-sample="{sample}" data-tab="{tab_id}">
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {sample}</h3>
                    <span class="total-badge">Total Genes: {total}</span>
                    <div class="typing-info">
                        <span class="typing-badge">ST: {mlst}</span>
                        <span class="typing-badge">Serotype: {sero}</span>
                        <span class="typing-badge">Phylogroup: {phy}</span>
                        <span class="typing-badge">CHtype: {cht}</span>
                        <span class="typing-badge {patho_class}">Pathotype: {patho}</span>
                    </div>
                    <button class="toggle-btn" onclick="toggleBoxDetails('{tab_id}', '{sample}')">
                        <i class="fas fa-chevron-down"></i> Show Details
                    </button>
                </div>
                <div class="box-details" style="display:none;"></div>
            </div>'''

        html += '</div>'
        html += f'''
        <script>
        window.ECOLITYPER_BOX_DATA = window.ECOLITYPER_BOX_DATA || {{}};
        window.ECOLITYPER_BOX_DATA["{tab_id}"] = {blob_json};
        </script>'''

        if 'amrfinder' in db_list:
            html += '<h3 style="margin-top:30px;"><i class="fas fa-book-medical"></i> Interpretation Guide</h3>'
            html += '<p style="color:#666;font-size:.9em;margin-bottom:10px;">Expand any section below for AMR interpretation guidance.</p>'
            html += self._collapsible('🔬 Why Do We Use Multiple Databases?', '#00695c', self._multi_db_education())
            html += self._collapsible('🎯 Cross-Database Confidence Tiers', '#0891b2', self._confidence_tiers())
            html += self._collapsible('🧬 Acquired vs Intrinsic Resistance', '#6f42c1', self._acquired_intrinsic())
            html += self._collapsible('⚠️ Gene Presence ≠ Phenotypic Resistance', '#f59e0b', self._genotype_phenotype_caveat())
            html += self._collapsible('🗄️ Roles & Strengths of Each AMR Database', '#3b82f6', self._db_roles_amr())

        return html

    # =====================================================================
    # Section methods
    # =====================================================================
    def _sec_summary(self, d: Dict) -> str:
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
        patho_detected = [k for k in patho_counts if k != 'none']
        mge_samples = len(d.get('mge_details', {}))
        mge_hits = sum(len(h) for h in d.get('mge_details', {}).values())

        features = [
            ('Gene-Centric Typing Tabs', 'fa-layer-group',
             'MLST, Serotype, CHTyper, and Phylogroup tables present every marker with all genomes that carry it, with combination tables for cross-tabulation.'),
            ('Sample-Centric Isolate Boxes', 'fa-microbe',
             'AMR, Virulence, Plasmid, BACMET, Mutations, and MGE open as collapsed boxes per isolate. Each box loads its tables on demand, keeping the report responsive on large datasets.'),
            ('Live Filtering', 'fa-filter',
             'Search boxes, database dropdowns, and grouping controls work across every table and every box, without reloading the page.'),
            ('Pathotype Classification', 'fa-virus',
             'Computational STEC / EPEC / EHEC / EAEC / ETEC / DAEC / EIEC prediction based on VFDB gene content.'),
            ('MGE Profiling', 'fa-mobile-alt',
             'mobileOG-db per-isolate hits with full positional and coverage detail, lazy-loaded inside each isolate box.'),
            ('Cross-Database Confidence', 'fa-shield-alt',
             'AMR hits are reported from all databases side-by-side; genes detected by multiple databases carry higher confidence.'),
            ('Compare and Cluster', 'fa-balance-scale',
             'Pair-compare two isolates with a similarity gauge, or run cluster detection across the entire dataset with a colour heatmap.'),
            ('Section-Specific Printing', 'fa-print',
             'Every tab has a print button that expands isolate boxes on the fly and prints the current view.'),
            ('Full Data Export', 'fa-download',
             'Every table exports to CSV. The complete integrated dataset exports to JSON for downstream analysis.'),
        ]
        tiles = ''
        for label, icon, desc in features:
            tiles += (f'<div class="database-section">'
                      f'<h4><i class="fas {icon}"></i> {label}</h4>'
                      f'<p>{desc}</p></div>')

        return f'''
        {self._alert('info', 'fa-info-circle',
            f'<h3>📊 Analysis Overview</h3>'
            f'<p>Comprehensive genomic analysis of <strong>{total}</strong> '
            f'<em>Escherichia coli</em> genomes. This hybrid report combines '
            f'gene-centric views for typing with sample-centric interactive boxes '
            f'for resistance, virulence, plasmid, biocide, mutation, and mobile '
            f'element detail.</p>')}

        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="switchTab('samples')">
                <i class="fas fa-list-alt"></i> Sample Overview</button>
            <button class="action-btn btn-danger" onclick="switchTab('amr')">
                <i class="fas fa-biohazard"></i> AMR</button>
            <button class="action-btn btn-success" onclick="switchTab('virulence')">
                <i class="fas fa-virus"></i> Virulence</button>
            <button class="action-btn btn-warning" onclick="switchTab('pathotype')">
                <i class="fas fa-chart-pie"></i> Pathotypes</button>
            <button class="action-btn btn-info" onclick="switchTab('mge')">
                <i class="fas fa-mobile-alt"></i> MGE</button>
            <button class="action-btn btn-secondary" onclick="switchTab('compare')">
                <i class="fas fa-balance-scale"></i> Compare</button>
        </div>

        <h3><i class="fas fa-chart-bar"></i> Key Statistics</h3>
        <div class="scrollable-table"><table class="data-table">
          <thead><tr><th>Metric</th><th>Count</th><th>Details</th></tr></thead>
          <tbody>
            <tr><td>Total Samples</td><td><strong>{total}</strong></td><td>Complete <em>E. coli</em> genomes</td></tr>
            <tr><td>Unique STs</td><td><strong>{len(patterns.get("st_distribution", {}))}</strong></td><td>MLST sequence types</td></tr>
            <tr><td>Unique Serotypes</td><td><strong>{len(patterns.get("serotype_distribution", {}))}</strong></td><td>O:H antigens</td></tr>
            <tr><td>Unique Phylogroups</td><td><strong>{len(patterns.get("phylogroup_distribution", {}))}</strong></td><td>Clermont scheme</td></tr>
            <tr><td>Unique CH Types</td><td><strong>{len(patterns.get("ch_type_distribution", {}))}</strong></td><td>FumC:FimH</td></tr>
            <tr><td>AMR Genes</td><td><strong>{amr_total}</strong></td><td>Across AMRFinderPlus + ABRicate databases</td></tr>
            <tr><td>Virulence Genes</td><td><strong>{vir_total}</strong></td><td>VFDB + EcoLi_VF</td></tr>
            <tr><td>BACMET Genes</td><td><strong>{bac_total}</strong></td><td>Biocide / heavy-metal resistance</td></tr>
            <tr><td>Plasmid / Serotyping Markers</td><td><strong>{pla_total}</strong></td><td>PlasmidFinder + EcoH</td></tr>
            <tr><td>MGE Hits</td><td><strong>{mge_hits}</strong></td><td>mobileOG-db hits across {mge_samples} sample(s)</td></tr>
            <tr><td>Pathotypes Detected</td><td><strong>{len(patho_detected)}</strong></td>
                <td>{', '.join(patho_detected) or 'None'}</td></tr>
            <tr><td>High-Risk Combos</td>
                <td><span class="badge" style="background:#8B0000;color:white;">{high_risk}</span></td>
                <td>Critical AMR + critical virulence</td></tr>
          </tbody>
        </table></div>

        <h3 style="margin-top:30px;"><i class="fas fa-lightbulb"></i> What This Report Provides</h3>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:20px;margin:20px 0;">
          {tiles}
        </div>'''

    def _sec_samples(self, d: Dict) -> str:
        samples = d['samples']
        rows = ''
        for s, data in sorted(samples.items()):
            mlst = data.get('mlst', {})
            sero = data.get('serotype', {})
            phy = data.get('phylogrouping', {})
            cht = data.get('chtyper', {})
            vir_genes = data.get('virulence_genes', [])
            patho = data.get('pathotype', 'none')

            vir_tags = ''.join(f'<span class="genome-tag">{g}</span>'
                               for g in sorted(set(vir_genes)))
            vir_details = (f'<details><summary>{len(vir_genes)} gene(s)</summary>'
                           f'<div class="genome-list">{vir_tags or "<em>None</em>"}</div>'
                           f'</details>')
            patho_badge = (f'<span class="badge typing-badge badge-pathotype">{patho}</span>'
                           if patho != 'none' else
                           '<span class="badge typing-badge badge-pathotype-none">none</span>')

            rows += f'''<tr>
                <td><strong>{s}</strong></td>
                <td>{mlst.get('ST','ND')}</td>
                <td>{sero.get('Serotype','ND')}</td>
                <td>{sero.get('O_Type','ND')}</td>
                <td>{sero.get('H_Type','ND')}</td>
                <td>{sero.get('Serotype_Confidence','ND')}</td>
                <td>{phy.get('Clermont_Type','ND')}</td>
                <td>{cht.get('FumC_Type','ND')}</td>
                <td>{cht.get('FimH_Type','ND')}</td>
                <td>{cht.get('CH_Type','ND')}</td>
                <td>{patho_badge}</td>
                <td>{vir_details}</td>
            </tr>'''

        return f'''
        {self._alert('info', 'fa-list-alt',
            '<h3>🧬 Sample Overview — The Population Snapshot</h3>'
            '<p>Every isolate on a single row, sourced from the master typing TSV.'
            ' Sort any column, '
            'search across all columns, and expand the Virulence column to reveal '
            'full VFDB gene lists inline.</p>')}

        <h3><i class="fas fa-book-medical"></i> What Each Column Means</h3>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:18px;margin:20px 0;">

          <div style="background:#fff;padding:18px;border-radius:10px;border-left:5px solid #FF9800;">
            <strong style="color:#FF9800;font-size:1.05em;">🔬 MLST (Multi-Locus Sequence Typing)</strong>
            <p style="font-size:.93em;margin-top:8px;line-height:1.55;">Sequences seven conserved housekeeping genes (<em>adk, fumC, gyrB, icd, mdh, purA, recA</em>). Each unique allele combination defines a <strong>Sequence Type (ST)</strong>. STs are the lingua franca of global <em>E. coli</em> epidemiology — ST131 (pandemic UPEC), ST10 (dominant commensal), ST11 (O157:H7), ST405, ST648 (emerging ESBL), ST95/ST73 (ExPEC).</p>
          </div>

          <div style="background:#fff;padding:18px;border-radius:10px;border-left:5px solid #9C27B0;">
            <strong style="color:#9C27B0;font-size:1.05em;">🧬 Serotype (O:H)</strong>
            <p style="font-size:.93em;margin-top:8px;line-height:1.55;">Traditional surface-antigen typing. <strong>O</strong> = lipopolysaccharide (LPS) antigen; <strong>H</strong> = flagellar antigen. Specific O:H combinations mark notorious pathotypes: O157:H7 (EHEC), O104:H4 (2011 German outbreak), O25b:H4 (ST131 UPEC), O26:H11 (STEC). <strong>Confidence</strong> is High / Medium / Low based on marker completeness.</p>
          </div>

          <div style="background:#fff;padding:18px;border-radius:10px;border-left:5px solid #795548;">
            <strong style="color:#795548;font-size:1.05em;">🌿 Phylogroup (Clermont scheme)</strong>
            <p style="font-size:.93em;margin-top:8px;line-height:1.55;">Deep evolutionary lineages: <strong>A, B1, B2, C, D, E, F, G</strong>. B2 and D carry most virulence factors and dominate extra-intestinal infections (UTIs, bloodstream); A and B1 are mostly commensals or environmental strains. Phylogrouping is the fastest way to place a new isolate in the global <em>E. coli</em> tree.</p>
          </div>

          <div style="background:#fff;padding:18px;border-radius:10px;border-left:5px solid #009688;">
            <strong style="color:#009688;font-size:1.05em;">🧬 CH Type (FumC:FimH)</strong>
            <p style="font-size:.93em;margin-top:8px;line-height:1.55;">Combines alleles of the housekeeping gene <em>fumC</em> and the type-1 fimbrial adhesin <em>fimH</em>. Provides <strong>higher resolution than MLST</strong> within a single clonal complex — two ST131 isolates can carry different <em>fimH</em> alleles and be epidemiologically distinct. Some <em>fimH</em> variants correlate with host tropism (human vs. avian).</p>
          </div>

          <div style="background:#fff;padding:18px;border-radius:10px;border-left:5px solid #E91E63;">
            <strong style="color:#E91E63;font-size:1.05em;">💊 Pathotype</strong>
            <p style="font-size:.93em;margin-top:8px;line-height:1.55;">Computational prediction from VFDB gene content. Not a clinical diagnosis — a hypothesis that should be confirmed by phenotype. STEC/EHEC are food-safety priorities; EAEC/ETEC are endemic diarrhoea drivers; UPEC/ExPEC dominate urinary and bloodstream infections.</p>
          </div>

          <div style="background:#fff;padding:18px;border-radius:10px;border-left:5px solid #F44336;">
            <strong style="color:#F44336;font-size:1.05em;">🧪 Virulence (expandable)</strong>
            <p style="font-size:.93em;margin-top:8px;line-height:1.55;">Click the count to reveal the full VFDB gene list inline. This is a fast way to see whether a suspicious isolate carries Shiga toxin, intimin, iron-capture systems, or capsule genes.</p>
          </div>
        </div>

        <div class="action-buttons" style="margin-bottom:10px;">
          <button class="action-btn btn-primary"
                  onclick="exportTableToCSV('samples-table','sample_overview.csv')">
              <i class="fas fa-download"></i> Export CSV</button>
        </div>

        <input type="text" class="search-box" id="search-samples"
               onkeyup="searchTable('samples-table','search-samples')"
               placeholder="🔍 Search samples by ID, ST, serotype, phylogroup, CH type, pathotype...">

        <div class="scrollable-table">
          <table id="samples-table" class="data-table">
            <thead><tr>
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
        </div>'''

    def _sec_qc(self, d: Dict) -> str:
        qc = d.get('qc_data', {})
        if not qc:
            return self._alert('warning', 'fa-exclamation-circle',
                '<h3>No QC Data Available</h3>'
                '<p>FASTA_QC_summary.html was not found in the input directory.</p>')

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

        return f'''
        {self._credit_bar('#17a2b8', '📊', 'FASTA QC Metrics – Tools & Acknowledgements',
            '<strong>Biopython</strong> by '
            '<a href="https://biopython.org/" target="_blank" style="color:#17a2b8;font-weight:bold;">The Biopython Consortium</a> '
            '→ computes assembly statistics (contigs, N50, GC%, total length) from FASTA files.<br>'
            '<strong>fastANI</strong> by '
            '<a href="https://github.com/ParBLiSS/FastANI" target="_blank" style="color:#17a2b8;font-weight:bold;">ParBLiSS Lab</a> '
            '→ provides average nucleotide identity (ANI) for species-level confirmation.<br>'
            '<span style="font-size:.9em;color:#6c757d;">'
            '<i class="fas fa-book-open"></i> Cite: Jain C, et al. <em>Nat Commun</em>. 2018;9(1):5114. '
            '<a href="https://doi.org/10.1038/s41467-018-07641-9" target="_blank" style="color:#17a2b8;font-weight:bold;">DOI</a></span>')}
        {self._alert('info', 'fa-chart-line',
            '<h3>📏 FASTA Quality Control</h3>'
            '<p>Assembly quality metrics for each genome. Poor assembly quality '
            'can lead to missing genes or false positives — always inspect QC '
            'before drawing biological conclusions.</p>'
            '<ul style="margin-top:8px;">'
            '<li><strong>Number of contigs</strong> – lower is better; typical '
            '<em>E. coli</em>: &lt;200 (good), &lt;100 (excellent).</li>'
            '<li><strong>N50</strong> – higher is better; &gt;100 kb excellent.</li>'
            '<li><strong>GC%</strong> – <em>E. coli</em> is typically ~50.8%.</li>'
            '<li><strong>Total length</strong> – ~4.6 Mb (commensal) to ~5.5 Mb (EHEC).</li>'
            '<li><strong>ANI ≥ 95%</strong> against a reference confirms species identity.</li>'
            '</ul>')}
        <input type="text" class="search-box" id="search-qc"
               onkeyup="searchTable('qc-table','search-qc')"
               placeholder="🔍 Search sample...">
        <div class="action-buttons">
          <button class="action-btn btn-primary"
                  onclick="exportTableToCSV('qc-table','fasta_qc.csv')">
              <i class="fas fa-download"></i> Export QC</button>
        </div>
        <div class="master-scrollable-container">
          <table id="qc-table" class="data-table">
            <thead><tr><th data-sort="string">Sample</th>{header}</tr></thead>
            <tbody>{rows}</tbody>
          </table>
        </div>'''

    def _sec_mlst(self, d: Dict) -> str:
        P = d['patterns']
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
            rows += (f'<tr><td><strong>{st}</strong></td><td>{cnt}</td>'
                     f'<td>{pct:.1f}%</td><td>{sero_str}</td></tr>')

        return f'''
        {self._credit_bar('#FF9800', '🧬', 'MLST Typing – Acknowledgments & Licensing',
            '<strong>MLST scheme</strong> powered by '
            '<a href="https://github.com/tseemann/mlst" target="_blank" style="color:#FF9800;font-weight:bold;">Prof. Torsten Seemann\'s scripts</a> '
            'and the '
            '<a href="https://pubmlst.org/" target="_blank" style="color:#FF9800;font-weight:bold;">PubMLST database</a> '
            '(Jolley et al., <em>Wellcome Open Res</em> 2018).<br>'
            '<span style="color:#856404;"><i class="fas fa-info-circle"></i> '
            '<strong>Note:</strong> Allele definitions current as of '
            '<strong>2024</strong>. Future updates may require manual downloads '
            'from the PubMLST website.</span>')}
        {self._alert('info', 'fa-code-branch',
            '<h3>🔬 MLST Analysis</h3>'
            '<p>MLST indexes internal fragments of seven housekeeping genes. '
            'Each unique combination defines a Sequence Type (ST) — the gold '
            'standard for global <em>E. coli</em> epidemiology. Closely related '
            'STs belong to the same Clonal Complex (CC).</p>'
            f'<p><strong>{len(st_dist)} unique STs</strong> identified.</p>')}
        <h3>📊 ST Distribution</h3>
        <div class="scrollable-table"><table id="mlst-dist" class="data-table">
          <thead><tr>
            <th data-sort="string">ST</th><th data-sort="number">Count</th>
            <th data-sort="number">%</th><th data-sort="string">Associated Serotypes</th>
          </tr></thead>
          <tbody>{rows}</tbody>
        </table></div>
        {self._combo_table("mlst-o-table", "ST – O antigen", P.get("st_o_combinations", {}))}
        {self._combo_table("mlst-h-table", "ST – H antigen", P.get("st_h_combinations", {}))}
        {self._combo_table("mlst-oh-table", "ST – Serotype (O:H)", P.get("st_oh_combinations", {}))}
        {self._combo_table("mlst-phylo-table", "ST – Phylogroup", P.get("st_phylogroup_combinations", {}))}
        {self._combo_table("mlst-fumc-table", "ST – FumC allele", P.get("st_fumc_combinations", {}))}
        {self._combo_table("mlst-fimh-table", "ST – FimH allele", P.get("st_fimh_combinations", {}))}
        {self._combo_table("mlst-ch-table", "ST – CH type (FumC:FimH)", P.get("st_fumc_fimh_combinations", {}))}
        {self._combo_table("mlst-full-table", "ST – O:H – Phylogroup – CH type", P.get("st_oh_phylo_fumc_fimh_combinations", {}))}'''

    def _sec_serotype(self, d: Dict) -> str:
        P = d['patterns']
        samples = d['samples']
        dist = P.get('serotype_distribution', Counter())
        total = sum(dist.values())
        rows = ''
        for sero, cnt in dist.most_common():
            if sero == 'ND': continue
            pct = cnt / total * 100 if total else 0
            sts = set()
            for combo in P.get('st_oh_combinations', {}):
                if combo.endswith(f" - {sero}"):
                    sts.add(combo.split(" - ")[0])
            st_list = ', '.join(sorted(sts)) if sts else 'Various'
            rows += (f'<tr><td><strong>{sero}</strong></td><td>{cnt}</td>'
                     f'<td>{pct:.1f}%</td><td>{st_list}</td></tr>')

        def _simple_table(distribution, label):
            out = ''
            for v, c in distribution.most_common():
                if v == 'ND': continue
                smp = [s for s, x in samples.items()
                       if x.get('serotype', {}).get(label) == v]
                tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in smp)
                out += (f'<tr><td><strong>{v}</strong></td><td>{c}</td>'
                        f'<td><div class="genome-list">{tags}</div></td></tr>')
            return out

        return f'''
        {self._credit_bar('#9C27B0', '🔬', 'Serotype Analysis – Acknowledgments',
            'Powered by the '
            '<a href="https://cge.food.dtu.dk/services/SerotypeFinder/" target="_blank" style="color:#9C27B0;font-weight:bold;">CGE SerotypeFinder</a> '
            'database (Center for Genomic Epidemiology) and '
            '<strong>EcoH</strong> (Ingle et al., 2016).<br>'
            'Please cite: Joensen KG, et al. Rapid and easy in silico serotyping '
            'of <em>Escherichia coli</em> using whole genome sequencing (WGS) data. '
            '<em>J Clin Microbiol</em>. 2015;53(8):2410-2426. '
            '<a href="https://doi.org/10.1128/JCM.00008-15" target="_blank" style="color:#9C27B0;font-weight:bold;">doi:10.1128/JCM.00008-15</a><br>'
            '<span style="font-size:.9em;color:#6c757d;">'
            '<i class="fas fa-gratipay"></i> We thank the CGE team, Flemming Scheutz, '
            'and SSI for maintaining this essential resource.</span>')}
        {self._alert('info', 'fa-tag',
            '<h3>🔬 Serotype Analysis (O:H)</h3>'
            '<p>Serotyping identifies lipopolysaccharide (O) and flagellar (H) '
            'antigens. Specific O:H combinations are associated with clinical '
            'syndromes and are used worldwide for food-safety and outbreak '
            'surveillance. The O antigen is a major vaccine target; the H antigen '
            'is a classic epidemiological marker.</p>'
            f'<p><strong>{len(dist)} unique serotypes</strong> identified.</p>')}
        <h3>📊 Serotype Distribution</h3>
        <div class="scrollable-table"><table class="data-table">
          <thead><tr>
            <th data-sort="string">Serotype</th>
            <th data-sort="number">Count</th>
            <th data-sort="number">%</th>
            <th data-sort="string">Common STs</th>
          </tr></thead>
          <tbody>{rows}</tbody>
        </table></div>
        <h3>📊 O Type Distribution</h3>
        <div class="scrollable-table"><table class="data-table">
          <thead><tr><th>O Type</th><th>Count</th><th>Samples</th></tr></thead>
          <tbody>{_simple_table(P.get("o_type_distribution", Counter()), "O_Type")}</tbody>
        </table></div>
        <h3>📊 H Type Distribution</h3>
        <div class="scrollable-table"><table class="data-table">
          <thead><tr><th>H Type</th><th>Count</th><th>Samples</th></tr></thead>
          <tbody>{_simple_table(P.get("h_type_distribution", Counter()), "H_Type")}</tbody>
        </table></div>'''

    def _sec_chtype(self, d: Dict) -> str:
        P = d['patterns']
        samples = d['samples']
        dist = P.get('ch_type_distribution', Counter())
        total = sum(dist.values())
        rows = ''
        for cht, cnt in dist.most_common():
            if cht == 'ND': continue
            pct = cnt / total * 100 if total else 0
            smp = [s for s, x in samples.items()
                   if x.get('chtyper', {}).get('CH_Type') == cht]
            tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in smp)
            rows += (f'<tr><td><strong>{cht}</strong></td><td>{cnt}</td>'
                     f'<td>{pct:.1f}%</td>'
                     f'<td><div class="genome-list">{tags}</div></td></tr>')

        fumc_dist = P.get('fumc_type_distribution', Counter())
        fumc_rows = ''
        for v, c in fumc_dist.most_common():
            if v == 'ND': continue
            smp = [s for s, x in samples.items()
                   if x.get('chtyper', {}).get('FumC_Type') == v]
            tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in smp)
            fumc_rows += (f'<tr><td><strong>{v}</strong></td><td>{c}</td>'
                          f'<td><div class="genome-list">{tags}</div></td></tr>')

        fimh_dist = P.get('fimh_type_distribution', Counter())
        fimh_rows = ''
        for v, c in fimh_dist.most_common():
            if v == 'ND': continue
            smp = [s for s, x in samples.items()
                   if x.get('chtyper', {}).get('FimH_Type') == v]
            tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in smp)
            fimh_rows += (f'<tr><td><strong>{v}</strong></td><td>{c}</td>'
                          f'<td><div class="genome-list">{tags}</div></td></tr>')

        return f'''
        {self._credit_bar('#009688', '🧬', 'CHTyper – Acknowledgments',
            'Powered by the '
            '<a href="https://cge.food.dtu.dk/services/CHTyper/" target="_blank" style="color:#009688;font-weight:bold;">CGE CHTyper</a> '
            'database curated by <strong>Henrik Hasman (CGE)</strong>.<br>'
            'Please cite: Roer L, et al. CHTyper, a Web Tool for Subtyping of '
            'Extraintestinal Pathogenic <em>Escherichia coli</em> Based on the '
            '<em>fumC</em> and <em>fimH</em> Alleles. <em>J Clin Microbiol</em>. '
            '2018;56(4):e00063-18. '
            '<a href="https://doi.org/10.1128/JCM.00063-18" target="_blank" style="color:#009688;font-weight:bold;">doi:10.1128/JCM.00063-18</a><br>'
            '<span style="font-size:.9em;color:#6c757d;">'
            '<i class="fas fa-gratipay"></i> We thank the CGE team and Henrik '
            'Hasman for maintaining this resource.</span>')}
        {self._alert('info', 'fa-project-diagram',
            '<h3>🧬 CH Type (FumC:FimH)</h3>'
            '<p>CH typing combines alleles of the housekeeping gene <em>fumC</em> '
            'and the type-1 fimbrial adhesin <em>fimH</em>. It provides high '
            'resolution for population studies and has been used to define '
            'epidemic clones. Higher resolution than MLST within a single clonal '
            'complex.</p>')}
        <h3>📊 CH Type Distribution</h3>
        <div class="scrollable-table"><table class="data-table">
          <thead><tr>
            <th data-sort="string">CH Type</th>
            <th data-sort="number">Count</th>
            <th data-sort="number">%</th>
            <th data-sort="string">Samples</th>
          </tr></thead>
          <tbody>{rows}</tbody>
        </table></div>
        <h3>📊 FumC Allele Distribution</h3>
        <div class="scrollable-table"><table class="data-table">
          <thead><tr><th>FumC Type</th><th>Count</th><th>Samples</th></tr></thead>
          <tbody>{fumc_rows}</tbody>
        </table></div>
        <h3>📊 FimH Allele Distribution</h3>
        <div class="scrollable-table"><table class="data-table">
          <thead><tr><th>FimH Type</th><th>Count</th><th>Samples</th></tr></thead>
          <tbody>{fimh_rows}</tbody>
        </table></div>'''

    def _sec_phylogroup(self, d: Dict) -> str:
        P = d['patterns']
        samples = d['samples']
        dist = P.get('phylogroup_distribution', Counter())
        total = sum(dist.values())
        rows = ''
        for ph, cnt in dist.most_common():
            if ph == 'ND': continue
            pct = cnt / total * 100 if total else 0
            smp = [s for s, x in samples.items()
                   if x.get('phylogrouping', {}).get('Clermont_Type') == ph]
            tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in smp)
            rows += (f'<tr><td><strong>{ph}</strong></td><td>{cnt}</td>'
                     f'<td>{pct:.1f}%</td>'
                     f'<td><div class="genome-list">{tags}</div></td></tr>')

        return f'''
        {self._credit_bar('#795548', '🌿', 'Phylogrouping – Acknowledgments',
            'Powered by '
            '<a href="https://github.com/nickp60/EzClermont" target="_blank" style="color:#795548;font-weight:bold;">EzClermont</a> '
            '(nickp60) — a Python implementation of the Clermont phylotyping method.<br>'
            'Original method: Clermont O, et al. The Clermont <em>Escherichia coli</em> '
            'phylo-typing method revisited: improvement of specificity and detection '
            'of new phylo-groups. <em>Environ Microbiol Rep</em>. 2013;5(1):58-65. '
            '<a href="https://doi.org/10.1111/1758-2229.12019" target="_blank" style="color:#795548;font-weight:bold;">doi:10.1111/1758-2229.12019</a>')}
        {self._alert('info', 'fa-sitemap',
            '<h3>🌿 Phylogroup Analysis (Clermont scheme)</h3>'
            '<p><em>E. coli</em> is divided into phylogroups (A, B1, B2, C, D, E, F, G) '
            'that reflect ecological niches and pathogenic potential. B2 and D are '
            'associated with extra-intestinal infections, while A and B1 are common '
            'in commensals and environmental strains.</p>')}
        <h3>📊 Phylogroup Distribution</h3>
        <div class="scrollable-table"><table class="data-table">
          <thead><tr>
            <th data-sort="string">Phylogroup</th>
            <th data-sort="number">Count</th>
            <th data-sort="number">%</th>
            <th data-sort="string">Samples</th>
          </tr></thead>
          <tbody>{rows}</tbody>
        </table></div>'''

    def _sec_pathotype(self, d: Dict) -> str:
        P = d['patterns']
        samples = d['samples']
        dist = P.get('pathotype_distribution', Counter())
        total = sum(dist.values())
        rows = ''
        for p, c in dist.most_common():
            pct = c / total * 100 if total else 0
            smp = [s for s, x in samples.items() if x.get('pathotype') == p]
            tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in smp[:200])
            rows += (f'<tr><td><strong>{p}</strong></td><td>{c}</td>'
                     f'<td>{pct:.1f}%</td>'
                     f'<td><div class="genome-list">{tags}</div></td></tr>')

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

        threshold_notice_inner = '''
        <p style="margin:6px 0;font-size:.95em;line-height:1.55;">
          Pathotype calls in this report are derived from <strong>gene presence/absence</strong>
          in the VFDB and EcoLi_VF databases, which are screened by <strong>ABRicate</strong>.
          This reporter does <strong>not</strong> apply its own identity or coverage cutoffs —
          it uses whatever thresholds were used when the ABRicate run was executed.
        </p>
        <p style="margin:10px 0 6px 0;font-size:.95em;line-height:1.55;">
          <strong>Default ABRicate thresholds:</strong>
        </p>
        <ul style="margin:6px 0 6px 22px;font-size:.93em;line-height:1.7;">
          <li><strong>Minimum identity:</strong> <code>80%</code> (<code>--minid 80</code>)</li>
          <li><strong>Minimum coverage:</strong> <code>80%</code> (<code>--mincov 80</code>)</li>
        </ul>
        <p style="margin:10px 0 6px 0;font-size:.95em;line-height:1.55;">
          A gene hit is only reported if it meets <strong>both</strong> thresholds.
          Lowering them increases sensitivity (more hits, more false positives);
          raising them increases specificity (fewer hits, more confidence).
        </p>
        <p style="margin:10px 0 0 0;font-size:.92em;line-height:1.5;background:#fff3cd;padding:10px 14px;border-radius:4px;border-left:3px solid #ffc107;">
          <i class="fas fa-lightbulb" style="color:#856404;"></i>
          <strong>Want stricter or looser calls?</strong> Re-run the ABRicate step
          of EcoliTyper with different <code>--minid</code> / <code>--mincov</code>
          flags, then regenerate this report. Whatever made it into the ABRicate
          HTML/TSV files is what this report uses.
        </p>'''

        how_it_works_inner = '''
        <p style="font-size:.93em;color:#555;margin-bottom:10px;">
          Pathotype is determined by simple marker presence in the VFDB and
          EcoLi_VF gene lists. A genome can receive <strong>multiple labels</strong>,
          joined by a dash (for example <code>EHEC-EAEC</code>).
        </p>
        <div class="scrollable-table">
          <table class="data-table" style="font-size:.9em;">
            <thead>
              <tr><th style="width:100px;">Pathotype</th><th>Rule</th><th>Clinical context</th></tr>
            </thead>
            <tbody>
              <tr><td><strong>STEC</strong></td><td><em>stx1</em> OR <em>stx2</em> present AND <em>eae</em> absent</td><td>Shiga toxin-producing — haemorrhagic colitis, HUS</td></tr>
              <tr><td><strong>EPEC</strong></td><td>(<em>eae</em> OR <em>bfpA</em>) present AND <em>stx</em> absent</td><td>Enteropathogenic — infantile diarrhoea</td></tr>
              <tr><td><strong>EHEC</strong></td><td><em>stx</em> present AND <em>eae</em> present</td><td>Enterohaemorrhagic — bloody diarrhoea, HUS</td></tr>
              <tr><td><strong>EAEC</strong></td><td><em>aggR</em> OR <em>aaiC</em> OR <em>aatA</em> present</td><td>Enteroaggregative — persistent diarrhoea</td></tr>
              <tr><td><strong>ETEC</strong></td><td><em>lt</em> OR <em>sta1</em> present</td><td>Enterotoxigenic — traveller\'s diarrhoea</td></tr>
              <tr><td><strong>DAEC</strong></td><td><em>afaC</em> OR <em>afaE</em> present</td><td>Diffusely adherent — paediatric diarrhoea</td></tr>
              <tr><td><strong>EIEC</strong></td><td><em>ipaH</em> present</td><td>Enteroinvasive — dysentery</td></tr>
              <tr><td><strong>none</strong></td><td>None of the above markers detected</td><td>No diarrhoeagenic pathotype predicted</td></tr>
            </tbody>
          </table>
        </div>
        <p style="font-size:.9em;color:#666;margin-top:12px;line-height:1.55;">
          <strong>Interpretation notes:</strong>
        </p>
        <ul style="font-size:.9em;color:#555;margin:6px 0 0 22px;line-height:1.6;">
          <li>Pathotype is a <strong>computational prediction</strong> — not a clinical diagnosis. Confirm with phenotype.</li>
          <li>Markers are inferred from ABRicate hits at the current identity/coverage thresholds (default 80/80).</li>
          <li>A gene that is present but below threshold will <strong>not</strong> be counted.</li>
          <li>Upstream/downstream context (for example LEE island intactness) is not evaluated — only marker gene presence.</li>
          <li>STEC/EHEC are food-safety priorities in most jurisdictions and are typically reportable.</li>
        </ul>'''

        return f'''
        {self._alert('info', 'fa-virus',
            '<h3>🧬 Pathotype Classification &amp; Prevalence</h3>'
            '<p><em>E. coli</em> pathotypes are defined by specific combinations of '
            'virulence factors. This section shows the prevalence of each pathotype '
            'in your dataset and lets you explore which samples belong to which category.</p>'
            '<p><strong>Definitions credit:</strong> '
            '<a href="https://github.com/kevinsblake/PathotypeR" target="_blank" style="color:#17a2b8;font-weight:bold;">Kevin Blake, PhD (PathotypeR)</a>. '
            'We are grateful for this open-source reference implementation of pathotype rules.</p>'
            '<ul style="margin-top:10px;">'
            '<li><strong>STEC</strong> (Shiga toxin-producing): <em>stx1/stx2</em>, <em>eae</em> absent → haemorrhagic colitis and HUS.</li>'
            '<li><strong>EPEC</strong> (Enteropathogenic): <em>eae/bfpA</em>, no <em>stx</em> → infant diarrhoea.</li>'
            '<li><strong>EHEC</strong> (Enterohaemorrhagic): <em>stx</em> + <em>eae</em> → severe bloody diarrhoea, HUS.</li>'
            '<li><strong>EAEC</strong> (Enteroaggregative): <em>aggR/aaiC/aatA</em> → persistent diarrhoea.</li>'
            '<li><strong>ETEC</strong> (Enterotoxigenic): LT/ST toxins → traveller\'s diarrhoea.</li>'
            '<li><strong>DAEC</strong> (Diffusely adherent): Afa/Dr adhesins → diarrhoea in children.</li>'
            '<li><strong>EIEC</strong> (Enteroinvasive): <em>ipaH</em> → dysentery-like illness.</li>'
            '<li><strong>none</strong>: No markers for the above pathotypes.</li>'
            '</ul>'
            '<p style="margin-top:10px;"><strong>Surveillance importance:</strong> Knowing '
            'which pathotypes circulate helps guide clinical management, vaccine '
            'development, and outbreak response. STEC/EHEC are particularly important '
            'for food safety.</p>')}

        <h3>📊 Pathotype Prevalence</h3>
        <div class="scrollable-table"><table class="data-table">
          <thead><tr>
            <th data-sort="string">Pathotype</th>
            <th data-sort="number">Count</th>
            <th data-sort="number">%</th>
            <th data-sort="string">Samples</th>
          </tr></thead>
          <tbody>{rows}</tbody>
        </table></div>

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
               onkeyup="searchTable('patho-samples','search-patho')"
               placeholder="🔍 Search by sample name, pathotype, ST, or serotype...">
        <div class="action-buttons">
          <button class="action-btn btn-primary"
                  onclick="exportTableToCSV('patho-samples','pathotype_samples.csv')">
              <i class="fas fa-download"></i> Export CSV</button>
        </div>
        <div class="scrollable-table"><table id="patho-samples" class="data-table">
          <thead><tr>
            <th>Sample</th><th>Pathotype</th><th>ST</th><th>Serotype</th>
            <th>Phylogroup</th><th>CH Type</th>
          </tr></thead>
          <tbody>{sample_rows}</tbody>
        </table></div>

        <h3 style="margin-top:30px;"><i class="fas fa-book-medical"></i> Interpretation Guide</h3>
        <p style="color:#666;font-size:.9em;margin-bottom:10px;">Expand the sections below for pathotype thresholds and calling rules.</p>
        {self._collapsible('📏 Identity &amp; Coverage Thresholds — Read This First', '#0891b2', threshold_notice_inner)}
        {self._collapsible('🔍 How Pathotype Is Called (exact rules)', '#1e3c72', how_it_works_inner)}'''

    def _sec_amr(self, d: Dict) -> str:
        return self._generate_sample_centric_boxes(
            d, 'amr', 'AMR',
            ['amrfinder', 'resfinder', 'card', 'argannot', 'megares', 'ncbi'])

    def _sec_virulence(self, d: Dict) -> str:
        return self._generate_sample_centric_boxes(
            d, 'virulence', 'Virulence', ['vfdb', 'ecoli_vf'])

    def _sec_plasmid(self, d: Dict) -> str:
        return self._generate_sample_centric_boxes(
            d, 'plasmid', 'Plasmid', ['plasmidfinder', 'ecoh'])

    def _sec_bacmet(self, d: Dict) -> str:
        return self._generate_sample_centric_boxes(
            d, 'bacmet', 'BACMET', ['bacmet2'])

    def _sec_mutation(self, d: Dict) -> str:
        """Mutation tab using sample-centric lazy boxes."""
        samples_data = d.get('samples', {})
        mutation_details = d.get('mutation_details', {})
        relevant = [s for s in samples_data if mutation_details.get(s)]
        relevant.sort()

        if not relevant:
            return self._alert('warning', 'fa-dna',
                '<h3>No Mutation Data Available</h3>'
                '<p>mutation_summary.tsv was not found or contains no rows.</p>')

        blob = {s: {'mutations': mutation_details.get(s, [])} for s in relevant}
        blob_json = json.dumps(blob, default=str).replace('</', '<\\/')

        credit = self._credit_bar('#00BCD4', '🧬', 'Point Mutations – AMRFinderPlus',
            'Point mutations detected by '
            '<a href="https://github.com/ncbi/amr" target="_blank" style="color:#00BCD4;font-weight:bold;">NCBI AMRFinderPlus</a>.<br>'
            'Please cite: Feldgarden M, et al. AMRFinderPlus and the Reference '
            'Gene Catalog. <em>Sci Rep</em>. 2021;11(1):12728. '
            '<a href="https://doi.org/10.1038/s41598-021-91456-0" target="_blank" style="color:#00BCD4;font-weight:bold;">doi:10.1038/s41598-021-91456-0</a>')

        info_inner = self._gene_family_info('#00BCD4',
            'Clinical relevance of key mutations:', [
                ('23S rRNA (linezolid)', 'Mutations such as G2576T and T2500A confer linezolid resistance.'),
                ('gyrA / parC (quinolones)', 'QRDR mutations reduce susceptibility to fluoroquinolones.'),
                ('rpoB (rifampin)', 'High-level rifampin resistance; combination therapy consideration.'),
                ('mprF (daptomycin)', 'Daptomycin non-susceptibility.'),
                ('rplC / rplD (linezolid)', 'Ribosomal protein mutations confer linezolid resistance.'),
                ('pmrAB / arnT (colistin)', 'Chromosomal mutations that confer colistin resistance.'),
                ('ampC promoter (cephalosporins)', 'Promoter mutations lead to ampC overexpression.'),
            ])

        html = f'''
        {credit}
        {self._alert('info', 'fa-dna',
            '<h3>🧬 Point Mutations – Sample-Centric View</h3>'
            '<p>Chromosomal mutations that confer resistance when no acquired gene '
            'is present. Each box represents one isolate; click '
            '<strong>Show Details</strong> to load its mutation table on demand.</p>')}
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-mutation"
                   onkeyup="filterBoxes('mutation')"
                   placeholder="🔍 Search sample by name..." style="max-width:320px;">
            <button class="action-btn btn-info" onclick="expandAllBoxes('mutation')">
                <i class="fas fa-expand-alt"></i> Expand All (visible)</button>
            <button class="action-btn btn-light" onclick="collapseAllBoxes('mutation')">
                <i class="fas fa-compress-alt"></i> Collapse All</button>
            <button class="action-btn btn-success" onclick="resetBoxFilters('mutation')">
                <i class="fas fa-sync"></i> Clear Filters</button>
            <span class="results-counter" id="counter-mutation">{len(relevant)} shown</span>
        </div>
        <div id="box-container-mutation">'''

        for sample in relevant:
            sd = samples_data.get(sample, {})
            mlst = sd.get('mlst', {}).get('ST', 'ND')
            sero = sd.get('serotype', {}).get('Serotype', 'ND')
            phy = sd.get('phylogrouping', {}).get('Clermont_Type', 'ND')
            cht = sd.get('chtyper', {}).get('CH_Type', 'ND')
            patho = sd.get('pathotype', 'none')
            n_mut = len(mutation_details.get(sample, []))
            patho_class = 'badge-pathotype' if patho != 'none' else 'badge-pathotype-none'

            html += f'''
            <div class="isolate-box" data-sample="{sample}" data-tab="mutation">
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {sample}</h3>
                    <span class="total-badge">Total Mutations: {n_mut}</span>
                    <div class="typing-info">
                        <span class="typing-badge">ST: {mlst}</span>
                        <span class="typing-badge">Serotype: {sero}</span>
                        <span class="typing-badge">Phylogroup: {phy}</span>
                        <span class="typing-badge">CHtype: {cht}</span>
                        <span class="typing-badge {patho_class}">Pathotype: {patho}</span>
                    </div>
                    <button class="toggle-btn" onclick="toggleBoxDetails('mutation', '{sample}')">
                        <i class="fas fa-chevron-down"></i> Show Details
                    </button>
                </div>
                <div class="box-details" style="display:none;"></div>
            </div>'''

        html += '</div>'
        html += f'''
        <script>
        window.ECOLITYPER_BOX_DATA = window.ECOLITYPER_BOX_DATA || {{}};
        window.ECOLITYPER_BOX_DATA["mutation"] = {blob_json};
        </script>'''
        html += '<h3 style="margin-top:30px;"><i class="fas fa-book-medical"></i> Interpretation Guide</h3>'
        html += '<p style="color:#666;font-size:.9em;margin-bottom:10px;">Expand the section below for mutation interpretation guidance.</p>'
        html += self._collapsible('🧪 Clinical relevance of key mutations', '#00BCD4', info_inner)
        return html

    def _sec_mge(self, d: Dict) -> str:
        """MGE tab using sample-centric lazy boxes."""
        samples_data = d.get('samples', {})
        mge_details = d.get('mge_details', {})
        relevant = [s for s in samples_data if mge_details.get(s)]
        relevant.sort()

        if not relevant:
            return self._alert('warning', 'fa-mobile-alt',
                '<h3>No MGE Data Available</h3>'
                '<p>ecolityper_mge_all_hits.tsv was not found or contains no rows.</p>')

        blob = {s: {'mge': mge_details.get(s, [])} for s in relevant}
        blob_json = json.dumps(blob, default=str).replace('</', '<\\/')

        credit = self._credit_bar('#16A085', '🧬', 'Mobile Genetic Elements – mobileOG-db',
            '<strong>mobileOG-db</strong> by '
            '<a href="https://github.com/clb21565/mobileOG-db" target="_blank" '
            'style="color:#16A085;font-weight:bold;">Brown, Mullet and colleagues</a>.<br>'
            'Please cite: Brown CL, et al. mobileOG-db: a Manually Curated Database of '
            'Protein Families Mediating the Life Cycle of Bacterial Mobile Genetic Elements. '
            '<em>Appl Environ Microbiol</em>. 2022;88(18):e00991-22. '
            '<a href="https://doi.org/10.1128/aem.00991-22" target="_blank" '
            'style="color:#16A085;font-weight:bold;">doi:10.1128/aem.00991-22</a>')

        info_inner = self._gene_family_info('#16A085',
            'MGE categories and why they matter:', [
                ('replication/recombination/repair', 'Core machinery for element replication and host DNA integration.'),
                ('transfer', 'Conjugal transfer systems (T4SS) — the primary route for horizontal gene transfer of AMR.'),
                ('integration/excision', 'Integrases and transposases that move genetic material between replicons.'),
                ('stability/defense', 'Partitioning and toxin-antitoxin systems ensuring persistence even without selection.'),
                ('phage', 'Bacteriophage structural and lysogeny proteins — can carry virulence genes such as Shiga toxin.'),
                ('plasmid', 'Plasmid backbone proteins — the most important vehicle for acquired antibiotic resistance.'),
                ('IS', 'Insertion sequence transposases — move resistance genes within and between genomes.'),
                ('ICE', 'Integrative conjugative elements — move large multi-gene cassettes between unrelated bacteria.'),
            ])

        html = f'''
        {credit}
        {self._alert('info', 'fa-mobile-alt',
            '<h3>🧬 Mobile Genetic Elements – Sample-Centric View</h3>'
            '<p>mobileOG-db protein family hits for each isolate. Individual matches '
            'represent MGE-associated protein signatures and should not be interpreted '
            'as independent mobile genetic elements. Click <strong>Show Details</strong> '
            'on any box to load its MGE table on demand.</p>')}
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-mge"
                   onkeyup="filterBoxes('mge')"
                   placeholder="🔍 Search sample by name..." style="max-width:320px;">
            <button class="action-btn btn-info" onclick="expandAllBoxes('mge')">
                <i class="fas fa-expand-alt"></i> Expand All (visible)</button>
            <button class="action-btn btn-light" onclick="collapseAllBoxes('mge')">
                <i class="fas fa-compress-alt"></i> Collapse All</button>
            <button class="action-btn btn-success" onclick="resetBoxFilters('mge')">
                <i class="fas fa-sync"></i> Clear Filters</button>
            <span class="results-counter" id="counter-mge">{len(relevant)} shown</span>
        </div>
        <div id="box-container-mge">'''

        for sample in relevant:
            sd = samples_data.get(sample, {})
            mlst = sd.get('mlst', {}).get('ST', 'ND')
            sero = sd.get('serotype', {}).get('Serotype', 'ND')
            phy = sd.get('phylogrouping', {}).get('Clermont_Type', 'ND')
            cht = sd.get('chtyper', {}).get('CH_Type', 'ND')
            patho = sd.get('pathotype', 'none')
            n_hits = len(mge_details.get(sample, []))
            patho_class = 'badge-pathotype' if patho != 'none' else 'badge-pathotype-none'

            html += f'''
            <div class="isolate-box" data-sample="{sample}" data-tab="mge">
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {sample}</h3>
                    <span class="total-badge">Total MGE hits: {n_hits}</span>
                    <div class="typing-info">
                        <span class="typing-badge">ST: {mlst}</span>
                        <span class="typing-badge">Serotype: {sero}</span>
                        <span class="typing-badge">Phylogroup: {phy}</span>
                        <span class="typing-badge">CHtype: {cht}</span>
                        <span class="typing-badge {patho_class}">Pathotype: {patho}</span>
                    </div>
                    <button class="toggle-btn" onclick="toggleBoxDetails('mge', '{sample}')">
                        <i class="fas fa-chevron-down"></i> Show Details
                    </button>
                </div>
                <div class="box-details" style="display:none;"></div>
            </div>'''

        html += '</div>'
        html += f'''
        <script>
        window.ECOLITYPER_BOX_DATA = window.ECOLITYPER_BOX_DATA || {{}};
        window.ECOLITYPER_BOX_DATA["mge"] = {blob_json};
        </script>'''
        html += '<h3 style="margin-top:30px;"><i class="fas fa-book-medical"></i> Interpretation Guide</h3>'
        html += '<p style="color:#666;font-size:.9em;margin-bottom:10px;">Expand the section below for MGE category descriptions.</p>'
        html += self._collapsible('🧪 MGE categories and why they matter', '#16A085', info_inner)
        return html

    def _sec_patterns(self, d: Dict) -> str:
        P = d['patterns']
        html = self._alert('info', 'fa-project-diagram',
            '<h3>🔍 Cross-Genome Pattern Discovery</h3>'
            '<p>Samples that carry both critical AMR genes (ESBL, carbapenemase, '
            '<em>mcr</em>) and critical virulence genes (<em>stx</em>, <em>eae</em>, '
            '<em>cnf1</em>) are flagged as high-risk. These combinations combine '
            'difficult-to-treat resistance with pathogenic potential.</p>'
            '<p>The co-occurrence table shows which pairs of genes are frequently '
            'found together across your samples.</p>')

        high = P.get('high_risk_combinations', [])
        if high:
            rows = ''
            for c in high:
                rows += (f'<tr><td><strong>{c["sample"]}</strong></td>'
                         f'<td>{c["st"]}</td><td>{c["serotype"]}</td>'
                         f'<td>{c.get("phylogroup","ND")}</td>'
                         f'<td>{", ".join(c["critical_amr_genes"])}</td>'
                         f'<td>{", ".join(c["critical_virulence_genes"])}</td></tr>')
            html += f'''
            <h3>⚠️ High-Risk Combinations ({len(high)})</h3>
            <div class="alert-box alert-danger">
              <i class="fas fa-radiation fa-2x"></i>
              <div><strong>{len(high)}</strong> sample(s) carry both critical AMR
              and critical virulence genes.</div>
            </div>
            <input type="text" class="search-box" id="search-highrisk"
                   onkeyup="searchTable('highrisk-table','search-highrisk')"
                   placeholder="🔍 Search high-risk samples...">
            <div class="master-scrollable-container">
              <table id="highrisk-table" class="data-table">
                <thead><tr>
                  <th>Sample</th><th>ST</th><th>Serotype</th><th>Phylogroup</th>
                  <th>Critical AMR</th><th>Critical Virulence</th>
                </tr></thead>
                <tbody>{rows}</tbody>
              </table>
            </div>'''

        cooc = P.get('gene_cooccurrence', {})
        if cooc:
            pairs = []
            for g1, partners in cooc.items():
                for g2, cnt in partners.items():
                    pairs.append((g1, g2, cnt))
            pairs.sort(key=lambda x: x[2], reverse=True)
            rows = ''.join(f'<tr><td>{g1}</td><td>{g2}</td><td>{c}</td></tr>'
                           for g1, g2, c in pairs[:500])
            html += f'''
            <h3>📈 Gene Co-occurrence (top 500)</h3>
            <input type="text" class="search-box" id="search-cooc"
                   onkeyup="searchTable('cooc-table','search-cooc')"
                   placeholder="🔍 Search gene pairs...">
            <div class="master-scrollable-container">
              <table id="cooc-table" class="data-table">
                <thead><tr><th>Gene 1</th><th>Gene 2</th><th>Co-occurrence</th></tr></thead>
                <tbody>{rows}</tbody>
              </table>
            </div>'''
        return html

    def _sec_compare(self, d: Dict) -> str:
        samples_data = d.get('samples', {})
        mutation_details = d.get('mutation_details', {})
        amr_details = d.get('amrfinder_details', {})
        abricate_details = d.get('abricate_details', {})

        if len(samples_data) < 2:
            return self._alert('warning', 'fa-balance-scale',
                f'<h3>Not Enough Samples</h3>'
                f'<p>Compare requires at least two samples. Dataset has '
                f'{len(samples_data)}.</p>')

        sample_mutations = defaultdict(list)
        for sample, muts in mutation_details.items():
            for m in muts:
                sample_mutations[sample].append(
                    f"{m.get('gene','')} {m.get('mutation','')}".strip())

        compare_data = {}
        for s, sd in samples_data.items():
            amr_genes = sorted(set(
                [g.get('gene', '') for g in amr_details.get(s, [])]
                + [g.get('gene', '') for g in abricate_details.get(s, {}).get('card', [])]
                + [g.get('gene', '') for g in abricate_details.get(s, {}).get('resfinder', [])]
                + [g.get('gene', '') for g in abricate_details.get(s, {}).get('argannot', [])]
                + [g.get('gene', '') for g in abricate_details.get(s, {}).get('megares', [])]
                + [g.get('gene', '') for g in abricate_details.get(s, {}).get('ncbi', [])]
            ))
            vir_genes = sorted(set(
                [g.get('gene', '') for g in abricate_details.get(s, {}).get('vfdb', [])]
                + [g.get('gene', '') for g in abricate_details.get(s, {}).get('ecoli_vf', [])]
            ))
            bac_genes = sorted(set(
                [g.get('gene', '') for g in abricate_details.get(s, {}).get('bacmet2', [])]
            ))
            pla_genes = sorted(set(
                [g.get('gene', '') for g in abricate_details.get(s, {}).get('plasmidfinder', [])]
                + [g.get('gene', '') for g in abricate_details.get(s, {}).get('ecoh', [])]
            ))
            muts = sorted(set(sample_mutations.get(s, [])))

            compare_data[s] = {
                'typing': {
                    'ST': sd.get('mlst', {}).get('ST', 'ND'),
                    'Serotype': sd.get('serotype', {}).get('Serotype', 'ND'),
                    'Phylogroup': sd.get('phylogrouping', {}).get('Clermont_Type', 'ND'),
                    'CH Type': sd.get('chtyper', {}).get('CH_Type', 'ND'),
                    'Pathotype': sd.get('pathotype', 'none'),
                },
                'genes': {
                    'AMR': [g for g in amr_genes if g],
                    'Virulence': [g for g in vir_genes if g],
                    'BACMET': [g for g in bac_genes if g],
                    'Plasmids': [g for g in pla_genes if g],
                    'Mutations': muts,
                },
            }

        sorted_samples = sorted(compare_data.keys())
        options_a = ''.join(
            f'<option value="{s}"{" selected" if i == 0 else ""}>{s}</option>'
            for i, s in enumerate(sorted_samples))
        options_b = ''.join(
            f'<option value="{s}"{" selected" if i == 1 else ""}>{s}</option>'
            for i, s in enumerate(sorted_samples))

        blob_json = json.dumps(compare_data, default=str).replace('</', '<\\/')

        return f'''
        <div class="alert-box alert-info">
          <i class="fas fa-balance-scale fa-2x"></i>
          <div>
            <h3>⚖️ Compare &amp; Cluster Analysis</h3>
            <p>Two modes: <strong>Pair Compare</strong> for a deep side-by-side of
            two isolates, and <strong>Cluster Detection</strong> to find groups of
            highly similar samples across the whole dataset.</p>
            <ul style="margin-top:6px;">
              <li><strong>Outbreak confirmation</strong> — do two isolates share ST + Serotype + Phylogroup + CH type?</li>
              <li><strong>Discordance analysis</strong> — what changed between an ancestor and a resistant descendant?</li>
              <li><strong>Cluster detection</strong> — which samples form tight transmission clusters?</li>
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
        </script>'''

    def _sec_aiguide(self) -> str:
        return '''
        <div class="alert-box alert-info">
          <i class="fas fa-robot fa-2x"></i>
          <div>
            <h3>🤖 AI Assistant Guide — Your Co-Pilot for Science, Not a Substitute</h3>
            <p>You have generated a rich dataset. Now let AI help you turn numbers into
            narratives — <strong>ethically, responsibly, and with you in the driver\'s seat</strong>.</p>
          </div>
        </div>

        <div class="info-text">
          <i class="fas fa-lightbulb"></i> <strong>The Golden Rule of AI in Science:</strong><br>
          AI is like a hyper-caffeinated grad student who never sleeps — great at summarising,
          spotting patterns, and drafting, but <strong>it can hallucinate, cite fake papers,
          and misunderstand context</strong>. Always verify, always think critically, and never
          let AI replace your expert judgment.
        </div>

        <h4>🚀 How to Use This Report with AI</h4>
        <ol style="margin-bottom:20px;margin-left:20px;">
          <li><strong>Export the full data as JSON</strong> (click the <em>Export</em> tab).</li>
          <li><strong>Upload the JSON file</strong> to your preferred AI chat interface (ChatGPT, Claude, or Gemini).</li>
          <li><strong>Ask specific, contextual questions</strong> — the AI will read your actual data and answer based on it.</li>
          <li><strong>Use the generated insights</strong> to draft sections of your manuscript or thesis — but <strong>always fact-check against the original tables</strong>.</li>
        </ol>

        <h4>📝 Using AI to Draft Manuscripts — A Workflow</h4>
        <div class="info-text" style="background:#f0f8ff;">
          <p><strong>Scenario:</strong> You analysed 200 E. coli genomes and want to write a paper on "Genomic epidemiology of ESBL-producing E. coli in Ghana".</p>
          <p><strong>Step 1 — Outline:</strong> <em>"Based on the JSON data, suggest a logical structure for a scientific paper covering population structure (MLST, serotypes), AMR gene prevalence, pathotype distribution, and high-risk combinations."</em></p>
          <p><strong>Step 2 — Methods:</strong> <em>"Write a methods paragraph describing how we used AMRFinderPlus, ABRicate (VFDB, PlasmidFinder, EcoH, BacMet2), and the EcoliTyper pipeline to classify pathotypes."</em></p>
          <p><strong>Step 3 — Results:</strong> <em>"From the JSON, what are the top 5 most common AMR genes and their percentages? Format as bullet points."</em></p>
          <p><strong>Step 4 — Discussion:</strong> <em>"Based on the prevalence of blaCTX-M (45%), mcr-1 (3%), and the identification of ST131 (30% of isolates), draft a discussion paragraph on the public health implications."</em></p>
          <p><strong>Step 5 — Polish:</strong> Ask the AI to rephrase or expand specific sentences, but <strong>manually verify every number and citation</strong>.</p>
        </div>

        <h4>🎓 For Thesis Writers</h4>
        <ul style="margin-left:20px;">
          <li><strong>Generate figure or table captions.</strong></li>
          <li><strong>Rewrite passive voice to active.</strong></li>
          <li><strong>Create a graphical abstract concept.</strong></li>
          <li><strong>Summarise limitations from the QC metrics.</strong></li>
        </ul>

        <h4>⚠️ Ethical Guardrails</h4>
        <div class="info-text" style="background:#fff3cd;border-left-color:#ffc107;">
          <ul style="margin:5px 0 5px 20px;">
            <li><strong>Never paste unreviewed AI output directly into a manuscript.</strong></li>
            <li><strong>Cite the AI tool if you use verbatim text</strong> (check journal policies).</li>
            <li><strong>Fact-check every number, gene name, and reference.</strong></li>
            <li><strong>Keep human experts in the loop.</strong> AI is a co-pilot, not the captain.</li>
            <li><strong>Use AI to democratise science</strong> — but always acknowledge its use.</li>
          </ul>
          <p style="margin-top:10px;"><i class="fas fa-laugh-beam"></i> <strong>Humour break:</strong> AI once told a colleague that "the plasmid IncF was discovered by Galileo in 1610". <strong>Don\'t be that person.</strong></p>
        </div>

        <div class="info-text" style="margin-top:20px;">
          <i class="fas fa-globe-africa"></i> <strong>Philosophy:</strong> Used ethically, AI can
          accelerate public health research, reduce language barriers, and empower researchers
          in low-resource settings. But it must never replace human curiosity, critical thinking,
          or domain expertise. <strong>You are the expert — AI is your enthusiastic, occasionally
          forgetful assistant.</strong>
        </div>'''

    def _sec_citation(self) -> str:
        entries = [
            ('#2196F3', 'EcoliTyper (Main Citation)',
             'Beckley B, Amarh V. EcoliTyper: a species-optimized computational pipeline for comprehensive genotyping and surveillance of <em>Escherichia coli</em>. <em>BMC Bioinformatics</em>. 2026.',
             'https://doi.org/10.1186/s12859-026-06529-6'),
            ('#e91e63', 'CHTyper',
             'Roer L, et al. CHTyper, a Web Tool for Subtyping of Extraintestinal Pathogenic <em>Escherichia coli</em> Based on the <em>fumC</em> and <em>fimH</em> Alleles. <em>J Clin Microbiol</em>. 2018;56(4):e00063-18.',
             'https://doi.org/10.1128/JCM.00063-18'),
            ('#4caf50', 'AMRFinderPlus',
             'Feldgarden M, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. <em>Sci Rep</em>. 2021;11(1):12728.',
             'https://doi.org/10.1038/s41598-021-91456-0'),
            ('#ff9800', 'CARD',
             'McArthur AG, et al. The comprehensive antibiotic resistance database. <em>Antimicrob Agents Chemother</em>. 2013;57(7):3348-57.',
             'https://doi.org/10.1128/AAC.00419-13'),
            ('#9c27b0', 'ARG-ANNOT',
             'Gupta SK, et al. ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. <em>Antimicrob Agents Chemother</em>. 2014;58(1):212-20.',
             'https://doi.org/10.1128/AAC.01310-13'),
            ('#00bcd4', 'MEGARes 2.0',
             'Doster E, et al. MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data. <em>Nucleic Acids Res</em>. 2020;48(D1):D561-D569.',
             'https://doi.org/10.1093/nar/gkz1010'),
            ('#f44336', 'VFDB',
             'Chen L, et al. VFDB 2012 update: toward the genetic diversity and molecular evolution of bacterial virulence factors. <em>Nucleic Acids Res</em>. 2012;40(Database issue):D641-5.',
             'https://doi.org/10.1093/nar/gkr989'),
            ('#3f51b5', 'PlasmidFinder',
             'Carattoli A, et al. <em>In silico</em> detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. <em>Antimicrob Agents Chemother</em>. 2014;58(7):3895-903.',
             'https://doi.org/10.1128/AAC.02412-14'),
            ('#ffc107', 'EcoH / Serotyping',
             'Ingle DJ, et al. In silico serotyping of <em>E. coli</em> from short read data identifies limited novel O-loci but extensive diversity of O:H serotype combinations. <em>Microb Genom</em>. 2016;2(7):e000064.',
             'https://doi.org/10.1099/mgen.0.000064'),
            ('#009688', 'BacMet',
             'Pal C, et al. BacMet: antibacterial biocide and metal resistance genes database. <em>Nucleic Acids Res</em>. 2014;42(Database issue):D737-43.',
             'https://doi.org/10.1093/nar/gkt1252'),
            ('#16a085', 'mobileOG-db',
             'Brown CL, et al. mobileOG-db: a Manually Curated Database of Protein Families Mediating the Life Cycle of Bacterial Mobile Genetic Elements. <em>Appl Environ Microbiol</em>. 2022;88(18):e00991-22.',
             'https://doi.org/10.1128/aem.00991-22'),
            ('#ff5722', 'Biopython',
             'Cock PJ, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. <em>Bioinformatics</em>. 2009;25(11):1422-3.',
             'https://doi.org/10.1093/bioinformatics/btp163'),
            ('#673ab7', 'PubMLST / BIGSdb',
             'Jolley KA, et al. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. <em>Wellcome Open Res</em>. 2018;3:124.',
             'https://doi.org/10.12688/wellcomeopenres.14826.1'),
            ('#cddc39', 'ResFinder',
             'Florensa AF, et al. ResFinder — an open online resource for identification of antimicrobial resistance genes in next-generation sequencing data. <em>Microb Genom</em>. 2022;8(1):000748.',
             'https://doi.org/10.1099/mgen.0.000748'),
            ('#8bc34a', 'ABRicate',
             'Seemann T. ABRicate: mass screening of contigs for antibiotic resistance genes. <a href="https://github.com/tseemann/abricate" target="_blank">github.com/tseemann/abricate</a>.',
             ''),
            ('#26c6da', 'MLST',
             'Seemann T. MLST: Scan contig files against PubMLST typing schemes. <a href="https://github.com/tseemann/mlst" target="_blank">github.com/tseemann/mlst</a>.',
             ''),
            ('#00897b', 'Clermont phylotyping',
             'Clermont O, et al. The Clermont <em>Escherichia coli</em> phylo-typing method revisited. <em>Environ Microbiol Rep</em>. 2013;5(1):58-65.',
             'https://doi.org/10.1111/1758-2229.12019'),
            ('#ffb300', 'Joensen serotyping',
             'Joensen KG, et al. Rapid and easy <em>in silico</em> serotyping of <em>Escherichia coli</em> using whole genome sequencing (WGS) data. <em>J Clin Microbiol</em>. 2015;53(8):2410-2426.',
             'https://doi.org/10.1128/JCM.00008-15'),
            ('#5c6bc0', 'KMA',
             'Clausen PTLC, et al. Rapid and precise alignment of raw reads against redundant databases with KMA. <em>BMC Bioinformatics</em>. 2018;19:307.',
             'https://doi.org/10.1186/s12859-018-2336-6'),
            ('#e91e63', 'BLAST+',
             'Camacho C, et al. BLAST+: architecture and applications. <em>BMC Bioinformatics</em>. 2009;10:421.',
             'https://doi.org/10.1186/1471-2105-10-421'),
            ('#7c3aed', 'PathotypeR',
             'Blake KS. PathotypeR: R package for classifying <em>E. coli</em> pathotypes. <a href="https://github.com/kevinsblake/PathotypeR" target="_blank">github.com/kevinsblake/PathotypeR</a>.',
             ''),
        ]

        cards = ''
        for color, name, text, url in entries:
            link = (f'<a href="{url}" target="_blank" style="color:{color};font-weight:bold;">'
                    f'{url}</a>') if url else ''
            cards += (f'<div style="background:{color}1a;padding:20px;border-radius:12px;'
                      f'margin-bottom:15px;border-left:5px solid {color};">'
                      f'<strong style="color:{color};font-size:1.02em;">{name}</strong><br>'
                      f'<span style="font-size:.95em;">{text} {link}</span></div>')

        return f'''
        <div class="alert-box alert-info" style="background:linear-gradient(135deg,#1e3c72 0%,#2a5298 100%);color:white;border-left:none;">
          <i class="fas fa-quote-right fa-2x"></i>
          <div>
            <h3 style="color:white;">📚 How to Cite EcoliTyper and Its Dependencies</h3>
            <p style="color:white;">If you use this pipeline in your research, please cite
            the main tool and the relevant third-party tools and databases below.</p>
          </div>
        </div>
        <div style="margin:20px 0;">{cards}</div>'''

    def _sec_funding(self) -> str:
        return '''
        <div style="background:linear-gradient(135deg,#43a047 0%,#1b5e20 100%);padding:25px;
                    border-radius:15px;margin-bottom:25px;color:white;text-align:center;">
          <i class="fas fa-coffee fa-3x" style="margin-bottom:10px;"></i>
          <h2 style="color:white;margin:0;">☕ No Grants, Just Passion (and a lot of caffeine)</h2>
          <p style="margin-top:10px;font-size:1.1em;"><strong>EcoliTyper</strong> is an independent,
          unfunded project born out of late nights, curiosity, and the urgent need for accessible
          genomic surveillance of <em>E. coli</em> in low-resource settings.</p>
          <p>No sponsors, no institutional backing — just a laptop, a coffee mug, and a burning
          desire to help researchers fight AMR.</p>
        </div>

        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:20px;margin:25px 0;">
          <div style="background:#fff3e0;padding:20px;border-radius:12px;border-left:5px solid #ff9800;">
            <i class="fas fa-star fa-2x" style="color:#ff9800;"></i>
            <h3 style="color:#e65100;margin:10px 0;">⭐ Star us on GitHub</h3>
            <p>One click — it costs nothing but makes us feel like rockstars.</p></div>
          <div style="background:#e8f5e9;padding:20px;border-radius:12px;border-left:5px solid #4caf50;">
            <i class="fas fa-bug fa-2x" style="color:#4caf50;"></i>
            <h3 style="color:#1b5e20;margin:10px 0;">🐛 Report bugs</h3>
            <p>If something breaks, tell us. We fix it with joy.</p></div>
          <div style="background:#e3f2fd;padding:20px;border-radius:12px;border-left:5px solid #2196f3;">
            <i class="fas fa-lightbulb fa-2x" style="color:#2196f3;"></i>
            <h3 style="color:#0d47a1;margin:10px 0;">💡 Suggest features</h3>
            <p>We are all ears, and we actually implement them.</p></div>
          <div style="background:#fce4ec;padding:20px;border-radius:12px;border-left:5px solid #e91e63;">
            <i class="fas fa-dna fa-2x" style="color:#e91e63;"></i>
            <h3 style="color:#880e4f;margin:10px 0;">🧬 Share your data</h3>
            <p>Used EcoliTyper? We would love to hear your story.</p></div>
          <div style="background:#f1f8e9;padding:20px;border-radius:12px;border-left:5px solid #8bc34a;">
            <i class="fas fa-share-alt fa-2x" style="color:#8bc34a;"></i>
            <h3 style="color:#33691e;margin:10px 0;">📢 Spread the word</h3>
            <p>Tweet, blog, or tell your lab mate.</p></div>
          <div style="background:#ede7f6;padding:20px;border-radius:12px;border-left:5px solid #673ab7;">
            <i class="fas fa-envelope fa-2x" style="color:#673ab7;"></i>
            <h3 style="color:#311b92;margin:10px 0;">👋 Say hello</h3>
            <p><strong>brownbeckley94@gmail.com</strong></p></div>
        </div>

        <div style="background:#ffebee;padding:20px;border-radius:12px;margin:20px 0;border-left:5px solid #f44336;">
          <i class="fas fa-microbe fa-2x" style="color:#f44336;float:left;margin-right:15px;"></i>
          <p style="margin:0;"><strong>Fun fact:</strong> This project runs on <strong>0% grant money,
          100% volunteer tears</strong>. But we are not bitter — we are just <strong>caffeinated</strong>.</p>
          <div style="clear:both;"></div>
        </div>

        <div style="background:linear-gradient(135deg,#0d47a1 0%,#1565c0 100%);padding:25px;
                    border-radius:15px;color:white;text-align:center;">
          <i class="fas fa-hand-holding-heart fa-3x"></i>
          <h3 style="color:white;">🤝 Join the Fight Against AMR</h3>
          <p>Open-source tools are the key to beating antimicrobial resistance.
          <strong>EcoliTyper</strong> is one part — we welcome <strong>your contributions</strong>:</p>
          <p>🐍 Code • 📚 Documentation • 🌍 New databases • 🤔 Ideas</p>
          <p style="margin-top:15px;">
            <a href="https://github.com/bbeckley-hub/ecolityper" target="_blank"
               style="background:white;color:#0d47a1;padding:10px 20px;border-radius:25px;
                      text-decoration:none;font-weight:bold;display:inline-block;margin:5px;">
              <i class="fab fa-github"></i> Star on GitHub</a>
            <a href="mailto:brownbeckley94@gmail.com"
               style="background:#ffc107;color:#0d47a1;padding:10px 20px;border-radius:25px;
                      text-decoration:none;font-weight:bold;display:inline-block;margin:5px;">
              <i class="fas fa-envelope"></i> Contact</a>
          </p>
        </div>'''

    def _sec_calltoaction(self) -> str:
        return '''
        <div class="alert-box alert-info">
          <i class="fas fa-globe fa-2x"></i>
          <div>
            <h3>The Global Burden of AMR and Our Call to Action</h3>
            <p>Antimicrobial resistance (AMR) is one of the top global public health threats,
            with an estimated <strong>1.27 million direct deaths annually</strong>.
            <em>Escherichia coli</em> is a major contributor, causing urinary tract infections,
            bloodstream infections, and foodborne outbreaks. Tracking AMR and virulence
            determinants is essential to inform treatment guidelines and infection control.</p>
            <p>We developed <strong>ECOLITYPER</strong> to empower researchers and clinicians —
            especially in low-resource settings — to analyse their own sequencing data without
            extensive bioinformatics expertise.</p>
          </div>
        </div>

        <div style="background:#e8f5e9;padding:20px;border-radius:12px;margin:20px 0;">
          <h3><i class="fas fa-bacterium"></i> ESCAPE AMR — Our Ongoing Project</h3>
          <p><strong>ECOLITYPER</strong> is one of the first modules of a larger initiative
          called <strong>ESCAPE AMR</strong> (formerly ESKAPE). We target the notorious
          <strong>ESKAPE pathogens</strong>:</p>
          <ul>
            <li><strong>E</strong>nterococcus faecium</li>
            <li><strong>S</strong>taphylococcus aureus</li>
            <li><strong>K</strong>lebsiella pneumoniae</li>
            <li><strong>A</strong>cinetobacter baumannii</li>
            <li><strong>P</strong>seudomonas aeruginosa</li>
            <li><strong>E</strong>nterobacter species (including <em>E. coli</em>)</li>
          </ul>
          <p>These bacteria "escape" the effects of antibiotics — hence the name. But the name
          is also a global call to action:</p>
          <div style="background:#fff3e0;padding:15px;border-radius:8px;margin:15px 0;">
            <p><strong>🔹 E</strong>veryone must join forces — researchers, clinicians, policymakers, citizens.<br>
            <strong>🔹 S</strong>mart surveillance is our first line of defence.<br>
            <strong>🔹 K</strong>nowledge must be shared openly. No paywalls, no closed silos.<br>
            <strong>🔹 A</strong>frica bears a heavy AMR burden, but African solutions are emerging.<br>
            <strong>🔹 P</strong>revention is cheaper than cure.<br>
            <strong>🔹 E</strong>very day we delay, more lives are at stake. Act now, not tomorrow.</p>
          </div>
        </div>

        <div style="text-align:center;margin:40px 0;">
          <i class="fas fa-star" style="font-size:3em;color:#ffc107;"></i>
          <h3>🤝 We Invite You to Contribute</h3>
          <p><strong>If you find this tool useful, please:</strong></p>
          <div class="action-buttons" style="justify-content:center;">
            <a href="https://github.com/bbeckley-hub/ecolityper" target="_blank"
               class="action-btn btn-primary" style="text-decoration:none;">
              <i class="fab fa-github"></i> ⭐ Star us on GitHub</a>
            <a href="mailto:brownbeckley94@gmail.com" class="action-btn btn-success"
               style="text-decoration:none;">
              <i class="fas fa-envelope"></i> Contact</a>
            <a href="https://github.com/bbeckley-hub/ecolityper/issues" target="_blank"
               class="action-btn btn-warning" style="text-decoration:none;">
              <i class="fas fa-bug"></i> Report issues</a>
          </div>
          <p style="margin-top:20px;"><i class="fas fa-chalkboard-user"></i> 
          <strong>We welcome collaborations</strong> to adapt this tool for other pathogens
          and improve AMR surveillance in Africa and beyond.</p>
        </div>

        <div class="info-text" style="background:#f8f9fa;">
          <i class="fas fa-quote-left"></i> "AMR is a silent pandemic, but we have the tools
          to fight it — if we share them, if we teach each other, and if we act with urgency.
          Let us escape the era of untreatable infections."<br>
          <strong>— The ESCAPE AMR Team, University of Ghana Medical School</strong>
        </div>'''

    def _sec_export(self) -> str:
        return '''
        <div class="alert-box alert-info">
          <i class="fas fa-download fa-2x"></i>
          <div>
            <h3>📥 Export Data</h3>
            <p>Download any table as CSV, or the complete dataset as JSON for
            AI-assisted analysis.</p>
          </div>
        </div>
        <div class="action-buttons">
          <button class="action-btn btn-primary"
                  onclick="exportTableToCSV('samples-table','sample_overview.csv')">
            <i class="fas fa-table"></i> Sample Overview CSV</button>
          <button class="action-btn btn-info"
                  onclick="exportTableToCSV('qc-table','fasta_qc.csv')">
            <i class="fas fa-chart-line"></i> FASTA QC CSV</button>
          <button class="action-btn btn-warning"
                  onclick="exportTableToCSV('patho-samples','pathotype_samples.csv')">
            <i class="fas fa-virus"></i> Pathotype Samples CSV</button>
          <button class="action-btn btn-success"
                  onclick="location.href='genius_ecoli_ultimate_gene_centric_report.json'">
            <i class="fas fa-file-code"></i> Complete JSON</button>
        </div>
        <p style="color:#666;font-size:.92em;">
          <i class="fas fa-info-circle"></i>
          AMR, Virulence, Plasmid, BACMET, Mutation, and MGE tabs are lazy-loaded.
          Click <strong>Show Details</strong> on each isolate box to reveal its
          tables, or open the JSON file for the full dataset.
        </p>'''

    # =====================================================================
    # Main report assembly
    # =====================================================================
    def generate_main_report(self, integrated_data: Dict[str, Any], output_dir: Path) -> str:
        """Assemble the full HTML report and write it to disk."""
        print("\nGenerating ULTIMATE E. coli HTML report...")
        samples_data = integrated_data.get('samples', {})
        patterns = integrated_data.get('patterns', {})
        gene_centric = integrated_data.get('gene_centric', {})
        metadata = integrated_data.get('metadata', {})

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
            'aiguide': lambda x: self._sec_aiguide(),
            'citation': lambda x: self._sec_citation(),
            'funding': lambda x: self._sec_funding(),
            'calltoaction': lambda x: self._sec_calltoaction(),
            'export': lambda x: self._sec_export(),
        }

        tabs_html = ''
        for i, (tid, title, icon) in enumerate(TABS):
            active = ' active' if i == 0 else ''
            content = section_methods[tid](integrated_data)
            tabs_html += f'''
            <div id="{tid}-tab" class="tab-content{active}">
                <h2 class="section-header {tid}-header" style="border-color:#1e3c72;">
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
        n_patho = len([k for k in patterns.get('pathotype_distribution', {}).keys()
                       if k != 'none'])
        n_amr = sum(len(g) for g in gene_centric.get('amr_databases', {}).values())
        n_vir = sum(len(g) for g in gene_centric.get('virulence_databases', {}).values())
        n_high = len(patterns.get('high_risk_combinations', []))
        mge_hits_total = sum(len(h) for h in integrated_data.get('mge_details', {}).values())

        dash = f'''
        <div class="dashboard-grid">
            <div class="dashboard-card" onclick="switchTab('summary')" style="border-left-color:#1e3c72;">
                <i class="fas fa-vial fa-2x" style="color:#1e3c72;"></i>
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
            <div class="dashboard-card" onclick="switchTab('mge')" style="border-left-color:#16A085;">
                <i class="fas fa-mobile-alt fa-2x" style="color:#16A085;"></i>
                <div class="card-number">{mge_hits_total}</div><div class="card-label">MGE Hits</div></div>
            <div class="dashboard-card" onclick="switchTab('patterns')" style="border-left-color:#FF5722;">
                <i class="fas fa-project-diagram fa-2x" style="color:#FF5722;"></i>
                <div class="card-number">{n_high}</div><div class="card-label">High-Risk Combos</div></div>
        </div>'''

        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>EcoliTyper Ultimate Hybrid Report v2.0.0</title>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
{css}
{js}
</head>
<body>
<div class="container">
  <div class="main-header">
    <h1><i class="fas fa-dna"></i> EcoliTyper Ultimate Analysis Report</h1>
    <p>Gene-Centric Typing · Interactive Sample-Centric Isolate Boxes · Pathotype Classification · MGE Profiling · Cluster Detection</p>
    <div class="metadata-bar">
      <div class="metadata-item"><i class="fas fa-calendar"></i><span>Generated: {metadata.get('analysis_date','Unknown')}</span></div>
      <div class="metadata-item"><i class="fas fa-database"></i><span>Samples: {total}</span></div>
      <div class="metadata-item"><i class="fas fa-user-md"></i><span>EcoliTyper Ultimate v2.0.0</span></div>
      <div class="metadata-item"><i class="fas fa-university"></i><span>University of Ghana Medical School</span></div>
    </div>
  </div>
  {dash}
  <div class="tab-navigation">{nav_html}</div>
  {tabs_html}
  <div class="footer">
    <h3>EcoliTyper Ultimate Hybrid Reporter v2.0.0</h3>
    <p>University of Ghana Medical School | Brown Beckley &lt;brownbeckley94@gmail.com&gt;</p>
    <p>If you find this tool helpful, please <a href="https://github.com/bbeckley-hub/ecolityper" target="_blank">⭐ star us on GitHub</a>.</p>
    <p>Tracked: ESBLs · Carbapenemases · Colistin resistance · Pathotype determinants · Biofilm · Plasmids · MGE</p>
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


# =============================================================================
# Orchestrator
# =============================================================================
class EcoliReporter:
    """Coordinates file discovery, parsing, analysis, and report generation."""

    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "GENIUS_ECOLI_ULTIMATE_GENE_CENTRIC_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = EcoliHTMLParser()
        self.analyzer = EcoliDataAnalyzer()
        self.generator = EcoliHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "EcoliTyper Ultimate Hybrid Reporter",
            "version": "2.1.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir),
        }

    def find_files(self) -> Dict[str, Any]:
        """Locate all recognised input files in the input directory."""
        print("Searching for input files...")
        found = {
            'master_tsv': None, 'qc_html': None,
        }
        for cand in ('ecolityper_comprehensive_report.tsv',
                     'ecoli_comprehensive_report.tsv'):
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
        return found

    def integrate(self, files: Dict[str, Any]) -> Dict[str, Any]:
        """Load every source, merge per-sample data, and build derived tables."""
        print("\nIntegrating data...")
        data = {
            'metadata': self.metadata,
            'samples': {},
            'patterns': {},
            'gene_centric': {},
            'qc_data': {},
            'amrfinder_details': {},
            'abricate_details': {},
            'mutation_details': {},
            'mge_details': {},
        }

        typing_data = (self.parser.load_master_typing_tsv(self.input_dir)
                       if files['master_tsv'] else {})
        if files['qc_html']:
            data['qc_data'] = self.parser.parse_qc_report(files['qc_html'])

        amr_details, amr_freq = self.parser.load_amrfinder_from_tsv(self.input_dir)
        abricate_details, abricate_freq = self.parser.load_abricate_from_tsv(self.input_dir)
        mutation_details = self.parser.load_mutations_from_tsv(self.input_dir)
        mge_details = self.parser.load_mge_from_tsv(self.input_dir)

        data['amrfinder_details'] = amr_details
        data['abricate_details'] = abricate_details
        data['mutation_details'] = mutation_details
        data['mge_details'] = mge_details

        all_samples = (set(typing_data) | set(amr_details) | set(abricate_details)
                       | set(data['qc_data']) | set(mutation_details) | set(mge_details))
        all_samples = sorted(all_samples)
        if not all_samples:
            print("No samples found in any source.")
            return {}
        print(f"Found {len(all_samples)} unique samples")

        for sample in all_samples:
            t = typing_data.get(sample, {})
            amr_genes = []
            for g in amr_details.get(sample, []):
                gene = g.get('gene')
                if gene:
                    amr_genes.append(gene)
            vir_genes = []
            bac_genes = []
            pla_genes = []
            for db, genes in abricate_details.get(sample, {}).items():
                gene_names = [g.get('gene', '') for g in genes if g.get('gene')]
                if db in ('vfdb', 'ecoli_vf'):
                    vir_genes.extend(gene_names)
                elif db == 'bacmet2':
                    bac_genes.extend(gene_names)
                elif db in ('plasmidfinder', 'ecoh'):
                    pla_genes.extend(gene_names)

            amr_genes = sorted(set(amr_genes))
            vir_genes = sorted(set(vir_genes))
            bac_genes = sorted(set(bac_genes))
            pla_genes = sorted(set(pla_genes))

            patho = self.analyzer.classify_pathotype(vir_genes)

            data['samples'][sample] = {
                'mlst': {
                    'ST': t.get('ST', 'ND'),
                    'Allele_Profile': t.get('Allele_Profile', 'ND'),
                },
                'serotype': {
                    'Serotype': t.get('Serotype', 'ND'),
                    'O_Type': t.get('O_Type', 'ND'),
                    'H_Type': t.get('H_Type', 'ND'),
                    'Serotype_Confidence': t.get('Serotype_Confidence', 'ND'),
                },
                'phylogrouping': {
                    'Clermont_Type': t.get('Clermont_Type', 'ND'),
                },
                'chtyper': {
                    'FumC_Type': t.get('FumC_Type', 'ND'),
                    'FimH_Type': t.get('FimH_Type', 'ND'),
                    'CH_Type': t.get('CH_Type', 'ND'),
                },
                'amr_genes': amr_genes,
                'virulence_genes': vir_genes,
                'bacmet_genes': bac_genes,
                'plasmid_genes': pla_genes,
                'pathotype': patho,
            }

        data['gene_frequencies'] = {
            'amrfinder': amr_freq,
            'abricate': abricate_freq,
        }
        print("\nBuilding gene-centric + pattern tables...")
        data['gene_centric'] = self.analyzer.create_gene_centric_tables(
            data, len(all_samples))
        data['patterns'] = self.analyzer.create_cross_genome_patterns(data)
        return data

    def write_json(self, data: Dict[str, Any]) -> Path:
        """Write the full integrated dataset as JSON."""
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
        """Write flat CSV exports next to the HTML report."""
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
                'BACMET_Count': len(d['bacmet_genes']),
                'Plasmid_Count': len(d['plasmid_genes']),
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
                        'Gene': g['gene'],
                        'Database': g['database'],
                        'Count': g['count'],
                        'Percentage': f"{(g['count'] / total) * 100:.1f}%",
                        'Genomes': ';'.join(g.get('genomes', [])),
                    })
            if out:
                pd.DataFrame(out).to_csv(self.output_dir / fname, index=False)

        muts = data.get('mutation_details', {})
        if muts:
            mut_rows = []
            for sample, mlist in muts.items():
                for m in mlist:
                    mut_rows.append({'Sample': sample, **m})
            if mut_rows:
                pd.DataFrame(mut_rows).to_csv(
                    self.output_dir / "mutations.csv", index=False)

        mge = data.get('mge_details', {})
        if mge:
            mge_rows = []
            for sample, hits in mge.items():
                for h in hits:
                    mge_rows.append({'Sample': sample, **h})
            if mge_rows:
                pd.DataFrame(mge_rows).to_csv(
                    self.output_dir / "mge_hits.csv", index=False)

        if data.get('qc_data'):
            pd.DataFrame(
                [{'Sample': s, **m} for s, m in data['qc_data'].items()]
            ).to_csv(self.output_dir / "fasta_qc.csv", index=False)

        print("    CSVs written")

    def run(self) -> bool:
        """Execute the full pipeline end-to-end."""
        print("=" * 80)
        print("ECOLITYPER ULTIMATE HYBRID REPORTER v2.0.0")
        print("=" * 80)
        print(f"Input:  {self.input_dir}")
        print(f"Output: {self.output_dir}")

        files = self.find_files()
        data = self.integrate(files)
        if not data:
            print("No data to process.")
            return False

        self.write_json(data)
        self.write_csvs(data)
        self.generator.generate_main_report(data, self.output_dir)

        n = len(data['samples'])
        patterns = data['patterns']
        gc = data['gene_centric']
        total_amr = sum(len(g) for g in gc.get('amr_databases', {}).values())
        total_vir = sum(len(g) for g in gc.get('virulence_databases', {}).values())
        high = len(patterns.get('high_risk_combinations', []))
        patho = len([k for k in patterns.get('pathotype_distribution', {})
                     if k != 'none'])
        mge_n = len(data.get('mge_details', {}))
        mge_hits = sum(len(h) for h in data.get('mge_details', {}).values())

        print("\n" + "=" * 80)
        print("ANALYSIS COMPLETE")
        print("=" * 80)
        print(f"   Samples:         {n}")
        print(f"   AMR genes:       {total_amr}")
        print(f"   Virulence genes: {total_vir}")
        print(f"   Pathotypes:      {patho}")
        print(f"   MGE hits:        {mge_hits} (across {mge_n} samples)")
        print(f"   High-risk:       {high}")
        print(f"   Output:          {self.output_dir}")
        print("=" * 80)
        return True


# =============================================================================
# Main
# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        description='EcoliTyper Ultimate Hybrid Reporter v2.0.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python genius_ultimate_gene_centric_reporter.py -i /path/to/reports

Author: Brown Beckley <brownbeckley94@gmail.com>
        """)
    parser.add_argument('-i', '--input-dir', required=True,
                        help='Directory containing EcoliTyper TSV summaries and QC HTML')
    parser.add_argument('-o', '--output-dir',
                        help='Custom output directory')
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