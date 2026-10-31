#!/usr/bin/env python3
"""
ECOLITYPER VISUALIZATION MODULE v2.0.0
========================================
Single-file interactive Plotly dashboard + publication-quality static exports
for Escherichia coli genomics data.

Reads the CSV/TSV outputs that EcoliTyper already produces:
  • ecolityper_comprehensive_report.tsv   master typing (Sample_ID, Phylogroup,
    Serotype, O_Type, H_Type, MLST_ST, FumC/FimH, CH_Type, etc.)
  • amr_genes.csv / virulence_genes.csv / bacmet_genes.csv / plasmid_replicons.csv
  • mutations.csv
  • mge_profile.csv  (wide)
  • fasta_qc.csv

Falls back to raw module TSVs if the reporter's CSVs are absent:
  • ecoli_amrfinder_summary.tsv
  • ecoli_<db>_abricate_summary.tsv
  • mutation_summary.tsv
  • ecolityper_mge_all_hits.tsv   (long format — aggregated)
  • FASTA_QC_summary.tsv

Produces:
  • ecolityper_dashboard.html            — flagship interactive Plotly dashboard
  • PNG/ PDF/ SVG/                       — publication-quality static exports
  • DATA/                                — flat CSVs of every plot's underlying data
  • ecolityper_visualization_report.txt  — plain-text summary
  • ecolityper_visualizations_bundle.zip — one-click export bundle

Dashboard tabs:
  Overview · Typing · Pathotype · QC · AMR · Virulence · Plasmids ·
  BACMET · Mutations · MGE · Alerts · Compare

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 2.0.0
Date: 2026-10-26
MIT License
"""

import os
import re
import sys
import json
import argparse
import zipfile
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Set, Tuple, Any, Optional
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np

# ------------------------------------------------------------------------------
# Matplotlib (static publication exports)
# ------------------------------------------------------------------------------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

plt.style.use('seaborn-v0_8-whitegrid')
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['figure.max_open_warning'] = 50

# ------------------------------------------------------------------------------
# Plotly (interactive dashboard)
# ------------------------------------------------------------------------------
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False


# ==============================================================================
# CONSTANTS
# ==============================================================================
TOOL_VERSION = "2.0.0"

# E. coli-specific palette (fresh blues/purples, distinct from  green)
PALETTE = {
    'st':        '#3b82f6',   # blue
    'serotype':  '#8b5cf6',   # violet
    'phylogroup':'#10b981',   # emerald
    'chtype':    '#f59e0b',   # amber
    'pathotype': '#ec4899',   # pink
    'amr':       '#dc2626',   # red
    'virulence': '#f97316',   # orange
    'bacmet':    '#64748b',   # slate
    'plasmid':   '#7c3aed',   # violet
    'mutation':  '#0891b2',   # cyan
    'mge':       '#14b8a6',   # teal
    'primary':   '#2563eb',   # royal blue
    'danger':    '#dc2626',
    'warning':   '#f59e0b',
    'success':   '#10b981',
}

CHART_COLORS = [
    '#3b82f6', '#8b5cf6', '#10b981', '#f59e0b', '#ec4899',
    '#06b6d4', '#84cc16', '#f97316', '#a855f7', '#14b8a6',
    '#ef4444', '#6366f1', '#22c55e', '#eab308', '#e879f9',
    '#0ea5e9', '#a3e635', '#fb923c', '#c084fc', '#2dd4bf',
    '#f43f5e', '#818cf8', '#4ade80', '#facc15', '#f472b6',
]

# ------------------------------------------------------------------------------
# Pathotype detection
# ------------------------------------------------------------------------------
PATHOTYPE_MARKERS = {
    'STEC': {'stx1', 'stx1a', 'stx1b', 'stx2', 'stx2a', 'stx2b',
             'stx2c', 'stx2d', 'stx2e', 'stx2f', 'stx2g'},
    'EPEC': {'eae', 'eaea', 'eaeh', 'bfpa'},
    'EHEC': {'eae', 'eaea', 'eaeh', 'bfpa', 'tir'},
    'EAEC': {'aggr', 'aata', 'aaic'},
    'ETEC': {'ltca', 'elta', 'eltb', 'sta1', 'estia'},
    'DAEC': {'afac', 'afae'},
    'EIEC': {'ipah'},
}

PATHOTYPE_DESCRIPTIONS = {
    'STEC': 'Shiga toxin-producing — haemorrhagic colitis, HUS',
    'EPEC': 'Enteropathogenic — infantile diarrhoea',
    'EHEC': 'Enterohaemorrhagic — bloody diarrhoea, HUS',
    'EAEC': 'Enteroaggregative — persistent diarrhoea',
    'ETEC': 'Enterotoxigenic — traveller\'s diarrhoea',
    'DAEC': 'Diffusely adherent — paediatric diarrhoea',
    'EIEC': 'Enteroinvasive — dysentery',
    'none': 'No diarrhoeagenic markers detected',
}

# ------------------------------------------------------------------------------
# Critical AMR classes for alert engine
# ------------------------------------------------------------------------------
AMR_CLASS_PATTERNS = {
    'beta-lactam':  ['bla', 'ctx-m', 'tem-', 'shv-', 'kpc', 'ndm',
                     'oxa-', 'vim-', 'imp-', 'ampc'],
    'aminoglycoside': ['aac', 'aph', 'ant(', 'arma', 'rmt'],
    'tetracycline': ['tet'],
    'macrolide':    ['erm', 'mph', 'mef', 'msr', 'ere'],
    'fluoroquinolone': ['qnr', 'gyra', 'parc'],
    'sulfonamide':  ['sul'],
    'trimethoprim': ['dfr'],
    'phenicol':     ['cat', 'flor', 'cmla'],
    'colistin':     ['mcr-'],
    'vancomycin':   ['van'],
    'fosfomycin':   ['fos'],
    'rifampin':     ['rpoB'],
}

CARBAPENEMASE_PATTERNS = ['kpc', 'ndm', 'oxa-48', 'oxa-181', 'oxa-232',
                          'vim-', 'imp-', 'spm', 'sim', 'gim']
ESBL_PATTERNS = ['ctx-m', 'shv-2', 'shv-5', 'shv-12', 'tem-3', 'tem-52',
                 'tem-10', 'tem-26', 'per-', 'veb', 'ges', 'bel']

ALERT_RULES = {
    'stec_detected': {
        'severity': 'critical', 'icon': '🚨',
        'title': 'STEC / EHEC Detected',
        'description': 'Shiga toxin-producing E. coli — food-safety priority; '
                       'often reportable.',
    },
    'eiec_detected': {
        'severity': 'high', 'icon': '💀',
        'title': 'EIEC Detected',
        'description': 'Enteroinvasive E. coli — dysentery-like illness, '
                       'reportable in many jurisdictions.',
    },
    'carbapenemase': {
        'severity': 'critical', 'icon': '☣️',
        'title': 'Carbapenemase producer',
        'description': 'Resistance to last-resort carbapenems.',
    },
    'mcr_stec': {
        'severity': 'critical', 'icon': '🚫',
        'title': 'Colistin-resistant STEC',
        'description': 'mcr-positive STEC — dual threat, public-health emergency.',
    },
    'esbl_stec': {
        'severity': 'high', 'icon': '⚠️',
        'title': 'ESBL-producing STEC',
        'description': 'Extended-spectrum β-lactamase in a Shiga-toxin background.',
    },
    'mdr': {
        'severity': 'medium', 'icon': '💊',
        'title': 'Multi-drug resistant (MDR)',
        'description': 'Resistance to 3 or more distinct drug classes.',
    },
}

# ------------------------------------------------------------------------------
# Gene deduplication
# ------------------------------------------------------------------------------
GENE_ALIASES = {
    'meca': 'meca', 'mec a': 'meca', 'mec-a': 'meca',
    'blactx-m': 'blactx-m', 'blactxm': 'blactx-m',
    'blactx-m-15': 'blactx-m-15', 'blactxm-15': 'blactx-m-15',
    'blatem': 'blatem', 'blatem-1': 'blatem-1',
    'blashv': 'blashv',
    'blakpc': 'blakpc', 'blakpc-2': 'blakpc-2',
    'blandm': 'blandm', 'blandm-1': 'blandm-1',
    'blaoxa': 'blaoxa', 'blaoxa-48': 'blaoxa-48',
    'mcr-1': 'mcr-1', 'mcr1': 'mcr-1',
    'mcr-1.1': 'mcr-1',
    'aaca-aphd': 'aaca-aphd', "aac(6')-aph(2'')": 'aaca-aphd',
    'aph3-iiia': 'aph3-iiia', "aph(3')-iiia": 'aph3-iiia',
    'ant4-ia': 'ant4-ia', "ant(4')-ia": 'ant4-ia',
    'sul1': 'sul1', 'sul2': 'sul2', 'sul3': 'sul3',
    'dfra1': 'dfra1', 'dfra12': 'dfra12',
    'tet(a)': 'teta', 'teta': 'teta',
    'tet(b)': 'tetb', 'tetb': 'tetb',
    'tet(m)': 'tetm', 'tetm': 'tetm',
    'cat': 'cat', 'cat1': 'cat',
    'flor': 'flor', 'flor1': 'flor',
    'ermb': 'ermb', 'erm(b)': 'ermb',
    'ermc': 'ermc', 'erm(c)': 'ermc',
    'stx1': 'stx1', 'stx1a': 'stx1', 'stx1b': 'stx1',
    'stx2': 'stx2', 'stx2a': 'stx2', 'stx2b': 'stx2',
    'stx2c': 'stx2', 'stx2d': 'stx2', 'stx2e': 'stx2',
    'stx2f': 'stx2', 'stx2g': 'stx2',
    'eae': 'eae', 'eaea': 'eae',
    'ipah': 'ipah', 'ipaH': 'ipah',
    'cnf1': 'cnf1', 'cnf2': 'cnf2',
    'hlya': 'hlya',
    'fima': 'fima', 'fimh': 'fimh',
    'papc': 'papc', 'papa': 'papa',
    'iuc': 'iuc', 'iucd': 'iuc',
    'iro': 'iro', 'iron': 'iro',
    'sitabcd': 'sitabcd',
}


# ==============================================================================
# HELPERS
# ==============================================================================
_EXTENSIONS = ('.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff')


def normalize_sample(name: str) -> str:
    """Strip paths and FASTA extensions from a sample identifier."""
    s = str(name).strip()
    for ext in _EXTENSIONS:
        if s.endswith(ext):
            s = s[: -len(ext)]
    if '/' in s or '\\' in s:
        s = Path(s).name
    for ext in _EXTENSIONS:
        if s.endswith(ext):
            s = s[: -len(ext)]
    return s.strip()


def clean_st(st: str) -> str:
    """Normalize an MLST sequence type to 'ST<n>' form."""
    s = str(st or '').strip()
    if not s or s in ('ND', 'nan', 'None', ''):
        return 'ND'
    s = re.sub(r'^(ST)+', 'ST', s, flags=re.I)
    if not s.upper().startswith('ST'):
        s = 'ST' + s
    return s


def esc(value: Any) -> str:
    """Escape a value for safe HTML embedding."""
    if value is None:
        return ""
    return (str(value).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def truncate_label(label: str, max_len: int = 35) -> str:
    """Truncate a label for axis rendering with an ellipsis."""
    s = str(label)
    if len(s) <= max_len:
        return s
    return s[: max_len - 1] + '…'


def auto_tickangle(n_categories: int) -> int:
    """Return a Plotly tickangle appropriate for the number of categories."""
    if n_categories <= 4:
        return 0
    if n_categories <= 8:
        return 25
    if n_categories <= 14:
        return 45
    if n_categories <= 22:
        return 60
    return 75


def auto_bar_height(n_categories: int, per_row: int = 22) -> int:
    """Return a Plotly figure height that scales with the number of categories."""
    return max(380, min(1400, 60 + per_row * n_categories))


def canonical_gene_key(name: str) -> str:
    """Return a canonical key for gene deduplication."""
    key = str(name).strip().lower()
    key = re.sub(r'\s+', ' ', key)
    return GENE_ALIASES.get(key, key)


def dedupe_gene_table(df: pd.DataFrame) -> pd.DataFrame:
    """Merge gene rows that differ only by case or alias."""
    if df.empty:
        return df
    df = df.copy()
    for col in ('Gene', 'Database', 'Count', 'Genomes'):
        if col not in df.columns:
            df[col] = 0 if col == 'Count' else ''
    df['_key'] = df['Gene'].astype(str).apply(canonical_gene_key)

    merged = []
    for key, group in df.groupby('_key', sort=False):
        genomes: Set[str] = set()
        for g in group['Genomes'].astype(str):
            for s in g.split(';'):
                s = s.strip()
                if s:
                    genomes.add(s)
        databases = sorted({str(d).strip() for d in group['Database']
                            if str(d).strip()})
        merged.append({
            'Gene': group.iloc[0]['Gene'],
            'Database': ','.join(databases),
            'Count': len(genomes) if genomes else int(group['Count'].max()),
            'Genomes': ';'.join(sorted(genomes)),
        })
    out = pd.DataFrame(merged).sort_values('Count', ascending=False).reset_index(drop=True)
    return out


def classify_pathotype(virulence_genes: List[str]) -> str:
    """Classify E. coli pathotype from VFDB gene content."""
    genes = {str(g).strip().lower() for g in virulence_genes if g}
    stx = any(g.startswith('stx') for g in genes)
    eae = 'eae' in genes or 'eaea' in genes or 'eaeh' in genes
    bfp = 'bfpa' in genes
    agg = any(g in genes for g in ('aggr', 'aata', 'aaic'))
    lt = any(g in genes for g in ('ltca', 'elta', 'eltb'))
    st = any(g in genes for g in ('sta1', 'estia'))
    afa = any(g in genes for g in ('afac', 'afae'))
    ipa = 'ipah' in genes

    pathotypes = []
    if stx and eae:
        pathotypes.append('EHEC')
    elif stx:
        pathotypes.append('STEC')
    if (eae or bfp) and not stx:
        pathotypes.append('EPEC')
    if agg:
        pathotypes.append('EAEC')
    if lt or st:
        pathotypes.append('ETEC')
    if afa:
        pathotypes.append('DAEC')
    if ipa:
        pathotypes.append('EIEC')
    return '-'.join(pathotypes) if pathotypes else 'none'


def match_amr_classes(gene_name: str) -> List[str]:
    """Return drug classes a gene belongs to (heuristic pattern matching)."""
    g = str(gene_name).lower()
    classes = []
    for cls, patterns in AMR_CLASS_PATTERNS.items():
        if any(p in g for p in patterns):
            classes.append(cls)
    return classes


# ==============================================================================
# DATA LOADER
# ==============================================================================
class EcoliDataLoader:
    """Loads typing, gene, mutation, MGE, and QC data from EcoliTyper outputs.

    Prefers the reporter's pre-aggregated CSVs; falls back to raw module TSVs.
    Handles both naming conventions ("ecolityper_*" and legacy "ecoli_*") so the
    visualizer works regardless of orchestrator naming.
    """

    MASTER_TSV_VARIANTS = [
        'ecolityper_comprehensive_report.tsv',
        'ecoli_comprehensive_report.tsv',
        '_comprehensive_report.tsv',  # legacy
    ]
    GENE_CSV_VARIANTS = {
        'amr':       ['amr_genes.csv'],
        'virulence': ['virulence_genes.csv'],
        'bacmet':    ['bacmet_genes.csv'],
        'plasmid':   ['plasmid_replicons.csv', 'plasmid_markers.csv'],
    }
    GENE_RAW_VARIANTS = {
        'amr': [
            'ecoli_amrfinder_summary.tsv',
            'ecolityper_amrfinder_summary.tsv',
        ],
        'virulence': [
            'ecoli_vfdb_abricate_summary.tsv',
            'ecoli_ecoli_vf_abricate_summary.tsv',
        ],
        'bacmet': ['ecoli_bacmet2_abricate_summary.tsv'],
        'plasmid': [
            'ecoli_plasmidfinder_abricate_summary.tsv',
            'ecoli_ecoh_abricate_summary.tsv',
        ],
    }
    MUTATION_VARIANTS = [
        'mutations.csv',
        'mutation_summary.tsv',
    ]
    MGE_VARIANTS = [
        'mge_profile.csv',
        'ecolityper_mge_all_hits.tsv',
        'ecoli_mge_all_hits.tsv',
    ]
    QC_VARIANTS = [
        'fasta_qc.csv',
        'FASTA_QC_summary.csv',
        'FASTA_QC_summary.tsv',
    ]

    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)

    # ------------------------------------------------------------------
    def _find(self, variants: List[str]) -> Optional[Path]:
        for name in variants:
            p = self.input_dir / name
            if p.exists():
                return p
        return None

    # ------------------------------------------------------------------
    def load_typing(self) -> pd.DataFrame:
        path = self._find(self.MASTER_TSV_VARIANTS)
        if not path:
            print("  ⚠️ Master typing TSV not found")
            return pd.DataFrame()

        df = pd.read_csv(path, sep='\t', dtype=str).fillna('')
        # Map all known column spellings to internal canonical names
        rename = {
            'Sample_ID': 'Sample',
            'Sample': 'Sample',
            'MLST_ST': 'ST',
            'MLST': 'ST',
            'Phylogroup': 'Phylogroup',
            'Clermont_Type': 'Phylogroup',
            'Serotype': 'Serotype',
            'O_Type': 'O_Type',
            'H_Type': 'H_Type',
            'CH_Type': 'CH_Type',
            'FumC_Type': 'FumC_Type',
            'FimH_Type': 'FimH_Type',
            'MRSA_Status': 'MRSA_Status',
        }
        df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})

        if 'Sample' not in df.columns:
            print(f"  ⚠️ No Sample column in {path.name}")
            return pd.DataFrame()

        df['Sample'] = df['Sample'].apply(normalize_sample)
        if 'ST' in df.columns:
            df['ST'] = df['ST'].apply(clean_st)
        for col in ('Phylogroup', 'Serotype', 'O_Type', 'H_Type', 'CH_Type'):
            if col in df.columns:
                df[col] = df[col].astype(str).str.strip()

        print(f"  ✅ Typing: {len(df)} samples from {path.name}")
        return df

    # ------------------------------------------------------------------
    def load_gene_tables(self) -> Dict[str, pd.DataFrame]:
        out: Dict[str, pd.DataFrame] = {}
        for cat, variants in self.GENE_CSV_VARIANTS.items():
            path = self._find(variants)
            if path:
                try:
                    df = pd.read_csv(path, dtype=str).fillna('')
                    df = dedupe_gene_table(df)
                    if not df.empty:
                        out[cat] = df
                        print(f"  ✅ {cat.capitalize()} genes: {len(df)} unique "
                              f"({path.name})")
                        continue
                except Exception as e:
                    print(f"  ⚠️ {path.name}: {e}")
            # Fallback to raw module TSVs
            raw_df = self._load_gene_from_raw(cat)
            if raw_df is not None and not raw_df.empty:
                out[cat] = raw_df
                print(f"  ✅ {cat.capitalize()} genes (rebuilt from raw): "
                      f"{len(raw_df)} unique")
        return out

    # ------------------------------------------------------------------
    def _load_gene_from_raw(self, category: str) -> Optional[pd.DataFrame]:
        """Rebuild gene table from raw ABRicate/AMRFinder TSV."""
        variants = self.GENE_RAW_VARIANTS.get(category, [])
        gene_col = 'Gene_Symbol' if category == 'amr' else 'gene'
        per_gene: Dict[str, Set[str]] = defaultdict(set)
        db_label = 'unknown'

        for name in variants:
            path = self.input_dir / name
            if not path.exists():
                continue
            try:
                df = pd.read_csv(path, sep='\t', dtype=str).fillna('')
            except Exception:
                continue
            if df.empty:
                continue
            sample_col = next((c for c in ('Genome', 'genome', 'file', 'sample')
                               if c in df.columns), None)
            if category == 'amr' and 'Gene_Symbol' in df.columns:
                effective_gene_col = 'Gene_Symbol'
            elif 'gene' in df.columns:
                effective_gene_col = 'gene'
            else:
                effective_gene_col = gene_col
            if not sample_col or effective_gene_col not in df.columns:
                continue
            db_label = 'AMRFinder' if category == 'amr' else \
                       name.replace('ecoli_', '').replace('_abricate_summary.tsv', '').upper()
            for _, row in df.iterrows():
                gene = str(row.get(effective_gene_col, '')).strip()
                sample = normalize_sample(str(row.get(sample_col, '')))
                if gene and sample:
                    per_gene[gene].add(sample)

        if not per_gene:
            return None

        rows = [{'Gene': g, 'Database': db_label,
                 'Count': len(samples),
                 'Genomes': ';'.join(sorted(samples))}
                for g, samples in per_gene.items()]
        return dedupe_gene_table(pd.DataFrame(rows))

    # ------------------------------------------------------------------
    def load_mutations(self) -> pd.DataFrame:
        path = self._find(self.MUTATION_VARIANTS)
        if not path:
            return pd.DataFrame()
        try:
            if path.suffix == '.csv':
                df = pd.read_csv(path, dtype=str).fillna('')
            else:
                df = pd.read_csv(path, sep='\t', dtype=str).fillna('')
        except Exception as e:
            print(f"  ⚠️ Mutations load failed: {e}")
            return pd.DataFrame()

        # Rebuild gene-centric view if we got raw per-sample rows
        if 'Count' not in df.columns and 'gene' in df.columns or 'gene_symbol' in df.columns:
            return self._aggregate_raw_mutations(df)

        if 'Count' in df.columns:
            df['Count'] = pd.to_numeric(df['Count'], errors='coerce').fillna(0).astype(int)
        df = df.sort_values('Count', ascending=False).reset_index(drop=True)
        print(f"  ✅ Mutations: {len(df)} entries ({path.name})")
        return df

    def _aggregate_raw_mutations(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert per-sample mutation rows to gene-centric (with Genomes)."""
        sample_col = next((c for c in ('genome', 'Genome', 'Sample', 'sample')
                          if c in df.columns), None)
        gene_col = next((c for c in ('gene_symbol', 'gene', 'Gene')
                        if c in df.columns), None)
        mut_col = next((c for c in ('element_name', 'mutation', 'Mutation')
                       if c in df.columns), None)
        if not sample_col or not gene_col:
            return pd.DataFrame()
        per_mut: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
        for _, row in df.iterrows():
            gene = str(row.get(gene_col, '')).strip()
            mut = str(row.get(mut_col, '')).strip() if mut_col else ''
            sample = normalize_sample(str(row.get(sample_col, '')))
            if gene and sample:
                per_mut[(gene, mut)].add(sample)
        rows = []
        for (gene, mut), samples in per_mut.items():
            rows.append({'Gene': gene, 'Mutation': mut,
                         'Count': len(samples),
                         'Genomes': ';'.join(sorted(samples))})
        out = pd.DataFrame(rows).sort_values('Count', ascending=False)
        print(f"  ✅ Mutations (rebuilt): {len(out)} entries")
        return out.reset_index(drop=True)

    # ------------------------------------------------------------------
    def load_mge(self) -> pd.DataFrame:
        """Return MGE data as a wide per-sample DataFrame with category counts."""
        for name in self.MGE_VARIANTS:
            path = self.input_dir / name
            if not path.exists():
                continue
            try:
                if name.endswith('.csv'):
                    df = pd.read_csv(path, dtype=str).fillna('0')
                    # Already-wide format (Sample + category columns)
                    sample_col = next((c for c in df.columns
                                      if c.lower() == 'sample'), None)
                    if sample_col:
                        df['Sample'] = df[sample_col].apply(normalize_sample)
                        for col in df.columns:
                            if col in ('Sample', sample_col):
                                continue
                            df[col] = pd.to_numeric(
                                df[col].astype(str).str.replace(',', ''),
                                errors='coerce'
                            ).fillna(0).astype(int)
                        print(f"  ✅ MGE: {len(df)} samples (wide, {path.name})")
                        return df
                else:
                    # Long format — aggregate by Sample × Major_Category
                    df = pd.read_csv(path, sep='\t', dtype=str).fillna('')
                    if 'Sample' not in df.columns or 'Major_Category' not in df.columns:
                        continue
                    df['Sample'] = df['Sample'].apply(normalize_sample)
                    wide = df.groupby(['Sample', 'Major_Category']).size().unstack(fill_value=0)
                    wide['mobileOG Hits'] = wide.sum(axis=1)
                    wide = wide.reset_index()
                    print(f"  ✅ MGE: {len(wide)} samples "
                          f"(aggregated from {path.name})")
                    return wide
            except Exception as e:
                print(f"  ⚠️ MGE load failed ({name}): {e}")
        return pd.DataFrame()

    # ------------------------------------------------------------------
    def load_qc(self) -> pd.DataFrame:
        path = self._find(self.QC_VARIANTS)
        if not path:
            return pd.DataFrame()
        try:
            sep = '\t' if path.suffix == '.tsv' else ','
            df = pd.read_csv(path, sep=sep, dtype=str).fillna('')
        except Exception as e:
            print(f"  ⚠️ QC load failed: {e}")
            return pd.DataFrame()

        sample_col = next((c for c in df.columns if c.lower() == 'sample'), None)
        if not sample_col:
            return pd.DataFrame()
        df['Sample'] = df[sample_col].apply(normalize_sample)

        numeric_keywords = ('n50', 'n75', 'n90', 'length', 'sequences',
                            'content', 'bases', 'run', 'homopolymer',
                            'duplicate', 'size', 'ani', 'contigs', 'gc')
        for col in df.columns:
            if col in ('Sample', sample_col):
                continue
            if any(k in col.lower() for k in numeric_keywords):
                df[col] = pd.to_numeric(
                    df[col].astype(str)
                          .str.replace('%', '').str.replace(',', '').str.strip(),
                    errors='coerce'
                )
        print(f"  ✅ QC: {len(df)} samples ({path.name})")
        return df


# ==============================================================================
# GENE MATRIX (samples × genes)
# ==============================================================================
def build_gene_matrix(gene_tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Build a samples × genes 0/1 presence matrix from the Genomes column."""
    rows: List[Tuple[str, str]] = []
    for df in (gene_tables or {}).values():
        if df is None or df.empty:
            continue
        for _, r in df.iterrows():
            gene = str(r.get('Gene', '')).strip()
            if not gene:
                continue
            for s in str(r.get('Genomes', '')).split(';'):
                s = s.strip()
                if s:
                    rows.append((s, gene))
    if not rows:
        return pd.DataFrame()
    long = pd.DataFrame(rows, columns=['Sample', 'Gene']).drop_duplicates()
    return long.assign(Value=1).pivot_table(
        index='Sample', columns='Gene', values='Value', fill_value=0
    )


# ==============================================================================
# ALERT ENGINE
# ==============================================================================
class AlertEngine:
    """Evaluates E. coli-specific clinical alert rules."""

    def __init__(self):
        self.alerts: List[Dict[str, Any]] = []

    def run(self, typing_df: pd.DataFrame,
            gene_tables: Dict[str, pd.DataFrame]) -> List[Dict[str, Any]]:
        # Build per-sample gene set (lowercase)
        sample_genes: Dict[str, Set[str]] = defaultdict(set)
        for df in (gene_tables or {}).values():
            if df is None or df.empty:
                continue
            for _, r in df.iterrows():
                gene = str(r.get('Gene', '')).strip().lower()
                for s in str(r.get('Genomes', '')).split(';'):
                    s = s.strip()
                    if s and gene:
                        sample_genes[s].add(gene)

        # Rule 1 — STEC / EHEC detection
        for sample, genes in sample_genes.items():
            if any(g.startswith('stx') for g in genes):
                has_eae = 'eae' in genes or 'eaea' in genes
                patho = 'EHEC' if has_eae else 'STEC'
                self.alerts.append({
                    'id': 'stec_detected', 'sample': sample,
                    'detail': f"{patho} — Shiga toxin gene present",
                })

        # Rule 2 — EIEC detection
        for sample, genes in sample_genes.items():
            if 'ipah' in genes:
                self.alerts.append({
                    'id': 'eiec_detected', 'sample': sample,
                    'detail': 'EIEC — ipaH invasion gene present',
                })

        # Rule 3 — Carbapenemase producer
        for sample, genes in sample_genes.items():
            hits = [g for g in genes
                    if any(p in g for p in CARBAPENEMASE_PATTERNS)]
            if hits:
                self.alerts.append({
                    'id': 'carbapenemase', 'sample': sample,
                    'detail': f"Carbapenemase: {', '.join(sorted(hits))}",
                })

        # Rule 4 — mcr + STEC
        for sample, genes in sample_genes.items():
            mcr = [g for g in genes if 'mcr-' in g]
            stx = [g for g in genes if g.startswith('stx')]
            if mcr and stx:
                self.alerts.append({
                    'id': 'mcr_stec', 'sample': sample,
                    'detail': f"mcr: {', '.join(sorted(mcr))} · "
                              f"stx: {', '.join(sorted(stx))}",
                })

        # Rule 5 — ESBL + STEC
        for sample, genes in sample_genes.items():
            esbl = [g for g in genes
                    if any(p in g for p in ESBL_PATTERNS)]
            stx = [g for g in genes if g.startswith('stx')]
            if esbl and stx and not any('mcr-' in g for g in genes):
                self.alerts.append({
                    'id': 'esbl_stec', 'sample': sample,
                    'detail': f"ESBL: {', '.join(sorted(esbl))} · "
                              f"stx present",
                })

        # Rule 6 — MDR (≥3 distinct drug classes)
        for sample, genes in sample_genes.items():
            classes_hit: Set[str] = set()
            for g in genes:
                for cls in match_amr_classes(g):
                    classes_hit.add(cls)
            if len(classes_hit) >= 3:
                self.alerts.append({
                    'id': 'mdr', 'sample': sample,
                    'detail': f"{len(classes_hit)} classes: "
                              f"{', '.join(sorted(classes_hit))}",
                })

        order = {'critical': 0, 'high': 1, 'medium': 2}
        self.alerts.sort(key=lambda a: order.get(
            ALERT_RULES[a['id']]['severity'], 3))
        return self.alerts


# ==============================================================================
# INTERACTIVE DASHBOARD BUILDER
# ==============================================================================
class DashboardBuilder:
    """Builds the single-file Plotly dashboard."""

    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.charts: Dict[str, str] = {}
        self.alerts: List[Dict[str, Any]] = []
        self.summary_facts: List[str] = []
        self.meta: Dict[str, Any] = {}
        self.sample_table_html: str = ''
        self.color_maps: Dict[str, Dict[str, str]] = {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def build(self, typing_df: pd.DataFrame,
              gene_tables: Dict[str, pd.DataFrame],
              mutations_df: pd.DataFrame,
              mge_df: pd.DataFrame,
              qc_df: pd.DataFrame) -> Optional[Path]:
        if typing_df.empty:
            print("  ⚠️ Dashboard skipped — no typing data")
            return None

        print("🎨 Building interactive dashboard...")
        gene_matrix = build_gene_matrix(gene_tables)

        # Pre-compute per-category color maps so categories keep the same
        # color across every chart
        self._precompute_colors(typing_df)

        # Add pathotype column
        if not typing_df.empty:
            typing_df = typing_df.copy()
            pathotypes = {}
            for sample in typing_df['Sample']:
                genes_for_sample = self._genes_for_sample(sample, gene_tables)
                pathotypes[sample] = classify_pathotype(genes_for_sample)
            typing_df['Pathotype'] = typing_df['Sample'].map(pathotypes)

        # Render view variants for typing charts
        for view in ('all', 'stec', 'esbl'):
            view_df = self._slice_view(typing_df, view, gene_tables)
            if view_df.empty and view != 'all':
                continue
            self._render_pathotype_donut(view, view_df)
            self._render_phylogroup_donut(view, view_df)
            self._render_typing_bar(view, view_df, 'ST', 'MLST Sequence Types')
            self._render_typing_bar(view, view_df, 'Serotype', 'Serotypes')
            self._render_typing_bar(view, view_df, 'Phylogroup', 'Phylogroups')
            self._render_typing_bar(view, view_df, 'CH_Type', 'CH Types')
            self._render_typing_bar(view, view_df, 'Pathotype', 'Pathotypes')
            self._render_sunburst(view, view_df)

        # View-independent charts
        self._render_amr_top(gene_tables)
        self._render_amr_matrix(gene_matrix)
        self._render_amr_by_phylogroup(typing_df, gene_tables)
        self._render_virulence_top(gene_tables)
        self._render_virulence_by_pathotype(typing_df, gene_tables)
        self._render_plasmid_top(gene_tables)
        self._render_bacmet_top(gene_tables)
        self._render_mutations_top(mutations_df)
        self._render_mutations_matrix(mutations_df)
        self._render_mge_breakdown(mge_df)
        self._render_qc_boxplots(qc_df)
        self._render_ani_chart(qc_df)

        # Data tables & alerts
        self.sample_table_html = self._build_sample_table(typing_df)
        engine = AlertEngine()
        self.alerts = engine.run(typing_df, gene_tables)
        self._narrate(typing_df, gene_tables, qc_df)

        self.meta = {
            'generated': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'n_samples': len(typing_df),
            'n_pathotypes': int((typing_df.get('Pathotype', '') != 'none').sum()),
            'n_stec': int((typing_df.get('Pathotype', '').str.contains('STEC', na=False)).sum()),
            'n_ehec': int((typing_df.get('Pathotype', '').str.contains('EHEC', na=False)).sum()),
            'n_alerts': len(self.alerts),
            'tool_version': TOOL_VERSION,
        }

        html = self._assemble_html(typing_df, qc_df)
        out = self.output_dir / 'ecolityper_dashboard.html'
        out.write_text(html, encoding='utf-8')
        size_mb = out.stat().st_size / (1024 * 1024)
        print(f"  ✅ Dashboard: {out.name} ({size_mb:.1f} MB)")
        return out

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    def _precompute_colors(self, typing_df: pd.DataFrame):
        for col in ('ST', 'Serotype', 'Phylogroup', 'CH_Type', 'Pathotype'):
            if col not in typing_df.columns:
                continue
            values = sorted({str(v) for v in typing_df[col].dropna().unique()
                             if v and str(v) not in ('Not Assigned', 'ND', '', 'nan')})
            self.color_maps[col] = {v: CHART_COLORS[i % len(CHART_COLORS)]
                                    for i, v in enumerate(values)}

    @staticmethod
    def _genes_for_sample(sample: str,
                          gene_tables: Dict[str, pd.DataFrame]) -> List[str]:
        genes = []
        for df in (gene_tables or {}).values():
            if df is None or df.empty:
                continue
            for _, r in df.iterrows():
                if sample in str(r.get('Genomes', '')).split(';'):
                    genes.append(str(r.get('Gene', '')).strip())
        return genes

    @staticmethod
    def _slice_view(df: pd.DataFrame, view: str,
                    gene_tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        if view == 'all':
            return df
        if 'Pathotype' in df.columns and view == 'stec':
            return df[df['Pathotype'].str.contains('STEC|EHEC', na=False)]
        if view == 'esbl':
            # Build per-sample ESBL presence
            esbl_samples: Set[str] = set()
            for df_gene in (gene_tables or {}).values():
                if df_gene is None or df_gene.empty:
                    continue
                for _, r in df_gene.iterrows():
                    g = str(r.get('Gene', '')).lower()
                    if any(p in g for p in ESBL_PATTERNS):
                        for s in str(r.get('Genomes', '')).split(';'):
                            s = s.strip()
                            if s:
                                esbl_samples.add(s)
            return df[df['Sample'].isin(esbl_samples)]
        return df

    # ------------------------------------------------------------------
    # Chart storage
    # ------------------------------------------------------------------
    def _store(self, key: str, fig, view: str = 'all'):
        html = fig.to_html(
            include_plotlyjs=False, full_html=False,
            div_id=f"chart-{key}-{view}",
            config={
                'displayModeBar': True, 'displaylogo': False,
                'modeBarButtonsToRemove': ['lasso2d', 'select2d', 'autoScale2d'],
                'toImageButtonOptions': {
                    'format': 'png',
                    'filename': f"ecolityper_{key}_{view}",
                    'height': 900, 'width': 1400, 'scale': 2,
                },
            },
        )
        self.charts[f"{key}_{view}"] = html

    def _chart_block(self, key: str, title: str, subtitle: str = '') -> str:
        variants = []
        for v in ('all', 'stec', 'esbl'):
            chart_key = f"{key}_{v}"
            if chart_key not in self.charts:
                continue
            display = 'block' if v == 'all' else 'none'
            variants.append(
                f'<div class="chart-view" data-view="{v}" style="display:{display};">'
                f'{self.charts[chart_key]}</div>'
            )
        if not variants:
            return ''
        sub_html = f'<p class="chart-subtitle">{esc(subtitle)}</p>' if subtitle else ''
        return f'''<div class="chart-card">
            <h3 class="chart-title">{title}</h3>
            {sub_html}
            {''.join(variants)}
        </div>'''

    # ------------------------------------------------------------------
    # Donut renderers
    # ------------------------------------------------------------------
    def _render_pathotype_donut(self, view: str, df: pd.DataFrame):
        if 'Pathotype' not in df.columns:
            return
        counts = df['Pathotype'].value_counts()
        if counts.empty:
            return
        colors = [self.color_maps.get('Pathotype', {}).get(l, PALETTE['pathotype'])
                  for l in counts.index]
        fig = go.Figure(go.Pie(
            labels=counts.index.tolist(),
            values=counts.values.tolist(),
            hole=0.55,
            marker=dict(colors=colors, line=dict(color='#0a0e1a', width=2)),
            textinfo='label+percent+value',
            textfont=dict(size=12),
            hovertemplate=('<b>%{label}</b><br>%{value} samples (%{percent})'
                          '<extra></extra>'),
        ))
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=400, showlegend=False,
            margin=dict(t=10, l=10, r=10, b=10),
            annotations=[dict(
                text=f"<b>{counts.sum()}</b><br>"
                     f"<span style='font-size:12px'>samples</span>",
                x=0.5, y=0.5,
                font=dict(size=24, color='#e2e8f0'),
                showarrow=False,
            )],
            font=dict(family='Segoe UI, sans-serif', size=13),
        )
        self._store('pathotype_donut', fig, view)

    def _render_phylogroup_donut(self, view: str, df: pd.DataFrame):
        if 'Phylogroup' not in df.columns:
            return
        valid = df[~df['Phylogroup'].isin(['Not Assigned', 'ND', '', 'nan'])]
        counts = valid['Phylogroup'].value_counts()
        if counts.empty:
            return
        colors = [self.color_maps.get('Phylogroup', {}).get(l, PALETTE['phylogroup'])
                  for l in counts.index]
        fig = go.Figure(go.Pie(
            labels=counts.index.tolist(),
            values=counts.values.tolist(),
            hole=0.55,
            marker=dict(colors=colors, line=dict(color='#0a0e1a', width=2)),
            textinfo='label+percent+value',
            textfont=dict(size=12),
            hovertemplate=('<b>%{label}</b><br>%{value} samples (%{percent})'
                          '<extra></extra>'),
        ))
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=400, showlegend=False,
            margin=dict(t=10, l=10, r=10, b=10),
            annotations=[dict(
                text=f"<b>{counts.sum()}</b><br>"
                     f"<span style='font-size:12px'>samples</span>",
                x=0.5, y=0.5,
                font=dict(size=24, color='#e2e8f0'),
                showarrow=False,
            )],
            font=dict(family='Segoe UI, sans-serif', size=13),
        )
        self._store('phylogroup_donut', fig, view)

    # ------------------------------------------------------------------
    # Typing bar
    # ------------------------------------------------------------------
    def _render_typing_bar(self, view: str, df: pd.DataFrame,
                           column: str, title: str, top_n: int = 25):
        if column not in df.columns:
            return
        valid = df[~df[column].isin(['Not Assigned', 'ND', '', 'nan'])]
        counts = valid[column].value_counts()
        if counts.empty:
            return
        top = counts.head(top_n)
        n = len(top)
        max_label_len = 45 if n <= 8 else 35 if n <= 15 else 28
        labels = [truncate_label(str(i), max_label_len) for i in top.index]
        colors = [self.color_maps.get(column, {}).get(i,
                                                       CHART_COLORS[idx % len(CHART_COLORS)])
                  for idx, i in enumerate(top.index)]
        sample_lists = [';'.join(
            valid[valid[column] == label]['Sample'].tolist())
            for label in top.index]

        total = len(df)
        pct_text = [f"{v} ({v/total*100:.0f}%)" for v in top.values]

        fig = go.Figure(go.Bar(
            x=top.values,
            y=labels,
            orientation='h',
            marker=dict(color=colors, line=dict(color='#0a0e1a', width=1)),
            text=pct_text,
            textposition='outside',
            customdata=sample_lists,
            hovertemplate=('<b>%{y}</b><br>%{x} samples'
                          '<br><i>Click to reveal sample list</i><extra></extra>'),
        ))
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=auto_bar_height(n),
            showlegend=False,
            margin=dict(t=20, l=180, r=80, b=20),
            yaxis=dict(autorange='reversed',
                      tickfont=dict(size=12), automargin=True),
            xaxis=dict(gridcolor='rgba(148,163,184,0.12)', title=None,
                      tickfont=dict(size=11)),
            font=dict(family='Segoe UI, sans-serif', size=12),
        )
        self._store(f'bar-{column}', fig, view)

    # ------------------------------------------------------------------
    # Sunburst — Phylogroup → ST → Serotype
    # ------------------------------------------------------------------
    def _render_sunburst(self, view: str, df: pd.DataFrame):
        if not all(c in df.columns for c in ('Phylogroup', 'ST', 'Serotype')):
            return
        df2 = df[~df['Phylogroup'].isin(['Not Assigned', 'ND', '', 'nan'])]
        df2 = df2[~df2['ST'].isin(['Not Assigned', 'ND', '', 'nan'])]
        if df2.empty:
            return

        ids, labels, parents, colors = [], [], [], []
        seen = set()

        def add(node_id, label, parent, color):
            if node_id in seen:
                return
            seen.add(node_id)
            ids.append(node_id)
            labels.append(label)
            parents.append(parent)
            colors.append(color)

        for pg in df2['Phylogroup'].unique():
            pg_color = self.color_maps.get('Phylogroup', {}).get(pg, PALETTE['phylogroup'])
            add(f'pg:{pg}', pg, '', pg_color)
            pg_df = df2[df2['Phylogroup'] == pg]
            for st in pg_df['ST'].unique():
                st_color = self.color_maps.get('ST', {}).get(st, PALETTE['st'])
                st_id = f'pg:{pg}|st:{st}'
                add(st_id, st, f'pg:{pg}', st_color)
                st_df = pg_df[pg_df['ST'] == st]
                for sero in st_df['Serotype'].unique():
                    if sero in ('Not Assigned', '', 'ND', 'nan'):
                        continue
                    sero_id = f'{st_id}|sero:{sero}'
                    sero_color = self.color_maps.get('Serotype', {}).get(sero, PALETTE['serotype'])
                    add(sero_id, sero, st_id, sero_color)

        if not ids:
            return

        # Compute values from leaf counts
        values = []
        for node_id in ids:
            parts = node_id.split('|')
            mask = pd.Series(True, index=df2.index)
            for part in parts:
                k, v = part.split(':', 1)
                if k == 'pg':
                    mask &= (df2['Phylogroup'] == v)
                elif k == 'st':
                    mask &= (df2['ST'] == v)
                elif k == 'sero':
                    mask &= (df2['Serotype'] == v)
            values.append(int(mask.sum()))

        fig = go.Figure(go.Sunburst(
            ids=ids, labels=labels, parents=parents, values=values,
            branchvalues='remainder',
            marker=dict(colors=colors, line=dict(color='#0a0e1a', width=2)),
            hovertemplate='<b>%{label}</b><br>Samples: %{value}<extra></extra>',
            textinfo='label+value',
            insidetextorientation='radial',
        ))
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=620,
            margin=dict(t=10, l=10, r=10, b=10),
            font=dict(family='Segoe UI, sans-serif', size=13),
        )
        self._store('sunburst', fig, view)

    # ------------------------------------------------------------------
    # AMR top genes bar
    # ------------------------------------------------------------------
    def _render_amr_top(self, gene_tables: Dict[str, pd.DataFrame]):
        df = (gene_tables or {}).get('amr')
        if df is None or df.empty:
            return
        top = df.nlargest(30, 'Count')
        n = len(top)
        labels = [truncate_label(g, 40) for g in top['Gene'].astype(str)]
        sample_lists = top['Genomes'].astype(str).str.replace(';', ', ')
        fig = go.Figure(go.Bar(
            x=top['Count'], y=labels,
            orientation='h',
            marker=dict(color=PALETTE['amr'],
                       line=dict(color='#0a0e1a', width=1)),
            text=top['Count'], textposition='outside',
            customdata=top['Genomes'],
            hovertemplate=('<b>%{y}</b><br>%{x} genomes'
                          '<br><i>Click to reveal genomes</i><extra></extra>'),
        ))
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=auto_bar_height(n, 24),
            showlegend=False,
            margin=dict(t=20, l=200, r=80, b=20),
            yaxis=dict(autorange='reversed',
                      tickfont=dict(size=11), automargin=True),
            xaxis=dict(gridcolor='rgba(148,163,184,0.12)', title=None),
            font=dict(family='Segoe UI, sans-serif', size=12),
        )
        self._store('amr_top', fig, 'all')

    # ------------------------------------------------------------------
    # AMR presence matrix
    # ------------------------------------------------------------------
    def _render_amr_matrix(self, gene_matrix: pd.DataFrame):
        if gene_matrix.empty:
            return
        top_genes = gene_matrix.sum(axis=0).sort_values(ascending=False).head(30).index.tolist()
        mat = gene_matrix[top_genes]
        mat = mat.loc[mat.sum(axis=1).sort_values(ascending=False).index]
        if mat.empty:
            return
        x_labels = [truncate_label(g, 22) for g in mat.columns]
        fig = go.Figure(go.Heatmap(
            z=mat.values,
            x=x_labels,
            y=mat.index.tolist(),
            colorscale=[[0.0, '#1e293b'], [1.0, PALETTE['amr']]],
            showscale=False, xgap=1, ygap=1,
            hovertemplate='<b>%{y}</b><br>%{x}: %{z}<extra></extra>',
        ))
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=max(420, 20 * len(mat)),
            margin=dict(t=20, l=10, r=10, b=10),
            xaxis=dict(tickangle=auto_tickangle(len(x_labels)),
                      tickfont=dict(size=10), automargin=True),
            yaxis=dict(tickfont=dict(size=10), automargin=True),
            font=dict(family='Segoe UI, sans-serif', size=11),
        )
        self._store('amr_matrix', fig, 'all')

    # ------------------------------------------------------------------
    # AMR by phylogroup
    # ------------------------------------------------------------------
    def _render_amr_by_phylogroup(self, typing_df: pd.DataFrame,
                                  gene_tables: Dict[str, pd.DataFrame]):
        if 'Phylogroup' not in typing_df.columns:
            return
        amr = (gene_tables or {}).get('amr')
        if amr is None or amr.empty:
            return

        top_genes = amr.nlargest(15, 'Count')['Gene'].tolist()
        phylogroups = [p for p in typing_df['Phylogroup'].unique()
                       if p not in ('Not Assigned', 'ND', '', 'nan')]
        phylogroups = sorted(phylogroups)

        # For each phylogroup × gene: fraction of phylogroup samples with the gene
        matrix = []
        for pg in phylogroups:
            pg_samples = set(typing_df[typing_df['Phylogroup'] == pg]['Sample'])
            pg_total = len(pg_samples)
            row = []
            for gene in top_genes:
                gene_row = amr[amr['Gene'] == gene]
                if gene_row.empty:
                    row.append(0.0)
                    continue
                gene_samples = set(str(gene_row.iloc[0]['Genomes']).split(';'))
                hits = len(pg_samples & gene_samples)
                row.append(hits / pg_total if pg_total else 0.0)
            matrix.append(row)

        x_labels = [truncate_label(g, 24) for g in top_genes]
        fig = go.Figure(go.Heatmap(
            z=matrix,
            x=x_labels,
            y=phylogroups,
            colorscale='Reds',
            showscale=True,
            colorbar=dict(title='Fraction', tickfont=dict(size=10)),
            hovertemplate='<b>%{y}</b><br>%{x}: %{z:.1%}<extra></extra>',
        ))
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=max(360, 60 + 50 * len(phylogroups)),
            margin=dict(t=20, l=80, r=20, b=80),
            xaxis=dict(tickangle=45, tickfont=dict(size=10), automargin=True),
            yaxis=dict(tickfont=dict(size=11), automargin=True),
            font=dict(family='Segoe UI, sans-serif', size=11),
        )
        self._store('amr_by_phylo', fig, 'all')

    # ------------------------------------------------------------------
    # Virulence top genes
    # ------------------------------------------------------------------
    def _render_virulence_top(self, gene_tables: Dict[str, pd.DataFrame]):
        df = (gene_tables or {}).get('virulence')
        if df is None or df.empty:
            return
        top = df.nlargest(25, 'Count')
        n = len(top)
        labels = [truncate_label(g, 35) for g in top['Gene'].astype(str)]
        fig = go.Figure(go.Bar(
            x=top['Count'], y=labels,
            orientation='h',
            marker=dict(color=PALETTE['virulence'],
                       line=dict(color='#0a0e1a', width=1)),
            text=top['Count'], textposition='outside',
            customdata=top['Genomes'],
            hovertemplate=('<b>%{y}</b><br>%{x} genomes'
                          '<br><i>Click to reveal genomes</i><extra></extra>'),
        ))
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=auto_bar_height(n, 24),
            showlegend=False,
            margin=dict(t=20, l=200, r=80, b=20),
            yaxis=dict(autorange='reversed',
                      tickfont=dict(size=11), automargin=True),
            xaxis=dict(gridcolor='rgba(148,163,184,0.12)', title=None),
            font=dict(family='Segoe UI, sans-serif', size=12),
        )
        self._store('virulence_top', fig, 'all')

    # ------------------------------------------------------------------
    # Virulence by pathotype
    # ------------------------------------------------------------------
    def _render_virulence_by_pathotype(self, typing_df: pd.DataFrame,
                                       gene_tables: Dict[str, pd.DataFrame]):
        if 'Pathotype' not in typing_df.columns:
            return
        vir = (gene_tables or {}).get('virulence')
        if vir is None or vir.empty:
            return
        top_genes = vir.nlargest(20, 'Count')['Gene'].tolist()
        pathotypes = [p for p in typing_df['Pathotype'].unique()
                      if p not in ('Not Assigned', 'ND', '', 'nan')]
        pathotypes = sorted(pathotypes)

        matrix = []
        for pt in pathotypes:
            pt_samples = set(typing_df[typing_df['Pathotype'] == pt]['Sample'])
            pt_total = len(pt_samples) if pt_samples else 1
            row = []
            for gene in top_genes:
                gene_row = vir[vir['Gene'] == gene]
                if gene_row.empty:
                    row.append(0.0)
                    continue
                gene_samples = set(str(gene_row.iloc[0]['Genomes']).split(';'))
                row.append(len(pt_samples & gene_samples) / pt_total)
            matrix.append(row)

        x_labels = [truncate_label(g, 22) for g in top_genes]
        fig = go.Figure(go.Heatmap(
            z=matrix,
            x=x_labels,
            y=pathotypes,
            colorscale='Oranges',
            showscale=True,
            colorbar=dict(title='Fraction', tickfont=dict(size=10)),
            hovertemplate='<b>%{y}</b><br>%{x}: %{z:.1%}<extra></extra>',
        ))
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=max(360, 80 + 50 * len(pathotypes)),
            margin=dict(t=20, l=100, r=20, b=100),
            xaxis=dict(tickangle=45, tickfont=dict(size=10), automargin=True),
            yaxis=dict(tickfont=dict(size=11), automargin=True),
            font=dict(family='Segoe UI, sans-serif', size=11),
        )
        self._store('virulence_by_patho', fig, 'all')

    # ------------------------------------------------------------------
    # Plasmid top
    # ------------------------------------------------------------------
    def _render_plasmid_top(self, gene_tables: Dict[str, pd.DataFrame]):
        df = (gene_tables or {}).get('plasmid')
        if df is None or df.empty:
            return
        top = df.nlargest(20, 'Count')
        n = len(top)
        labels = [truncate_label(g, 35) for g in top['Gene'].astype(str)]
        fig = go.Figure(go.Bar(
            x=top['Count'], y=labels,
            orientation='h',
            marker=dict(color=PALETTE['plasmid'],
                       line=dict(color='#0a0e1a', width=1)),
            text=top['Count'], textposition='outside',
            customdata=top['Genomes'],
            hovertemplate=('<b>%{y}</b><br>%{x} genomes'
                          '<br><i>Click to reveal genomes</i><extra></extra>'),
        ))
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=auto_bar_height(n, 24),
            showlegend=False,
            margin=dict(t=20, l=200, r=80, b=20),
            yaxis=dict(autorange='reversed',
                      tickfont=dict(size=11), automargin=True),
            xaxis=dict(gridcolor='rgba(148,163,184,0.12)', title=None),
            font=dict(family='Segoe UI, sans-serif', size=12),
        )
        self._store('plasmid_top', fig, 'all')

    # ------------------------------------------------------------------
    # BACMET top
    # ------------------------------------------------------------------
    def _render_bacmet_top(self, gene_tables: Dict[str, pd.DataFrame]):
        df = (gene_tables or {}).get('bacmet')
        if df is None or df.empty:
            return
        top = df.nlargest(25, 'Count')
        n = len(top)
        labels = [truncate_label(g, 40) for g in top['Gene'].astype(str)]
        fig = go.Figure(go.Bar(
            x=top['Count'], y=labels,
            orientation='h',
            marker=dict(color=PALETTE['bacmet'],
                       line=dict(color='#0a0e1a', width=1)),
            text=top['Count'], textposition='outside',
            customdata=top['Genomes'],
            hovertemplate=('<b>%{y}</b><br>%{x} genomes'
                          '<br><i>Click to reveal genomes</i><extra></extra>'),
        ))
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=auto_bar_height(n, 24),
            showlegend=False,
            margin=dict(t=20, l=220, r=80, b=20),
            yaxis=dict(autorange='reversed',
                      tickfont=dict(size=11), automargin=True),
            xaxis=dict(gridcolor='rgba(148,163,184,0.12)', title=None),
            font=dict(family='Segoe UI, sans-serif', size=12),
        )
        self._store('bacmet_top', fig, 'all')

    # ------------------------------------------------------------------
    # Mutations top
    # ------------------------------------------------------------------
    def _render_mutations_top(self, mutations_df: pd.DataFrame):
        if mutations_df is None or mutations_df.empty:
            return
        top = mutations_df.head(30).copy()
        if 'Mutation' in top.columns:
            labels = [truncate_label(f"{r['Gene']} {r['Mutation']}", 40)
                     for _, r in top.iterrows()]
        else:
            labels = [truncate_label(g, 40) for g in top['Gene']]
        n = len(top)
        fig = go.Figure(go.Bar(
            x=top['Count'], y=labels,
            orientation='h',
            marker=dict(color=PALETTE['mutation'],
                       line=dict(color='#0a0e1a', width=1)),
            text=top['Count'], textposition='outside',
            customdata=top['Genomes'] if 'Genomes' in top.columns else None,
            hovertemplate=('<b>%{y}</b><br>%{x} genomes'
                          '<br><i>Click to reveal genomes</i><extra></extra>'),
        ))
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=auto_bar_height(n, 24),
            showlegend=False,
            margin=dict(t=20, l=220, r=80, b=20),
            yaxis=dict(autorange='reversed',
                      tickfont=dict(size=11), automargin=True),
            xaxis=dict(gridcolor='rgba(148,163,184,0.12)', title=None),
            font=dict(family='Segoe UI, sans-serif', size=12),
        )
        self._store('mutations_top', fig, 'all')

    # ------------------------------------------------------------------
    # Mutations matrix
    # ------------------------------------------------------------------
    def _render_mutations_matrix(self, mutations_df: pd.DataFrame):
        if mutations_df is None or mutations_df.empty or 'Genomes' not in mutations_df.columns:
            return
        top = mutations_df.head(25)
        # Build sample → mutation map
        sample_muts: Dict[str, Set[str]] = defaultdict(set)
        for _, r in top.iterrows():
            mut_label = f"{r.get('Gene', '')} {r.get('Mutation', '')}".strip()
            for s in str(r['Genomes']).split(';'):
                s = s.strip()
                if s:
                    sample_muts[s].add(mut_label)
        if not sample_muts:
            return
        samples = sorted(sample_muts.keys())
        mutations = sorted({m for muts in sample_muts.values() for m in muts})
        if not mutations or not samples:
            return
        matrix = []
        for sample in samples:
            matrix.append([1 if m in sample_muts[sample] else 0 for m in mutations])

        x_labels = [truncate_label(m, 24) for m in mutations]
        fig = go.Figure(go.Heatmap(
            z=matrix,
            x=x_labels,
            y=samples,
            colorscale=[[0, '#1e293b'], [1, PALETTE['mutation']]],
            showscale=False, xgap=1, ygap=1,
            hovertemplate='<b>%{y}</b><br>%{x}: %{z}<extra></extra>',
        ))
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=max(400, 20 * len(samples)),
            margin=dict(t=20, l=10, r=10, b=10),
            xaxis=dict(tickangle=45, tickfont=dict(size=10), automargin=True),
            yaxis=dict(tickfont=dict(size=9), automargin=True),
            font=dict(family='Segoe UI, sans-serif', size=10),
        )
        self._store('mutations_matrix', fig, 'all')

    # ------------------------------------------------------------------
    # MGE breakdown
    # ------------------------------------------------------------------
    def _render_mge_breakdown(self, mge_df: pd.DataFrame):
        if mge_df is None or mge_df.empty:
            return
        skip = {'Sample', 'mobileOG Hits'}
        cats = [c for c in mge_df.columns if c not in skip
                and pd.api.types.is_numeric_dtype(mge_df[c])]
        if not cats:
            return
        means = {c: float(mge_df[c].mean()) for c in cats}
        cats_sorted = sorted(means, key=lambda c: -means[c])
        n = len(cats_sorted)

        fig = go.Figure(go.Bar(
            x=cats_sorted, y=[means[c] for c in cats_sorted],
            marker=dict(color=CHART_COLORS[:n],
                       line=dict(color='#0a0e1a', width=1)),
            text=[f"{means[c]:.1f}" for c in cats_sorted],
            textposition='outside',
            hovertemplate='<b>%{x}</b><br>Mean hits/sample: %{y:.1f}<extra></extra>',
        ))
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=440, showlegend=False,
            margin=dict(t=20, l=20, r=20, b=100),
            xaxis=dict(tickangle=auto_tickangle(n), automargin=True,
                      tickfont=dict(size=11)),
            yaxis=dict(title='Mean mobileOG hits',
                      gridcolor='rgba(148,163,184,0.12)',
                      title_font=dict(size=12)),
            font=dict(family='Segoe UI, sans-serif', size=12),
        )
        self._store('mge_breakdown', fig, 'all')

    # ------------------------------------------------------------------
    # QC box plots
    # ------------------------------------------------------------------
    def _render_qc_boxplots(self, qc_df: pd.DataFrame):
        if qc_df is None or qc_df.empty:
            return
        metrics = []
        for keyword, label, color in [
            ('N50', 'N50 (bp)', PALETTE['st']),
            ('Total Length', 'Assembly size (bp)', PALETTE['phylogroup']),
            ('GC Content', 'GC content (%)', PALETTE['mge']),
            ('Total Sequences', 'Number of contigs', PALETTE['plasmid']),
            ('Contigs', 'Contigs', PALETTE['plasmid']),
        ]:
            matching = [c for c in qc_df.columns if keyword.lower() in c.lower()]
            for col in matching:
                if pd.api.types.is_numeric_dtype(qc_df[col]) and qc_df[col].notna().any():
                    metrics.append((col, label, color))
                    break

        if not metrics:
            return

        n = len(metrics)
        cols = min(2, n)
        rows = (n + cols - 1) // cols
        fig = make_subplots(
            rows=rows, cols=cols,
            subplot_titles=[m[1] for m in metrics],
            vertical_spacing=0.18,
        )
        for i, (col, label, color) in enumerate(metrics):
            r = i // cols + 1
            c = i % cols + 1
            series = qc_df[col].dropna()
            if series.empty:
                continue
            fig.add_trace(
                go.Box(
                    y=series, name=label, marker_color=color,
                    boxmean='sd', boxpoints='all',
                    jitter=0.4, pointpos=0,
                    hovertemplate=f'<b>{label}</b><br>%{{y}}<extra></extra>',
                ),
                row=r, col=c,
            )
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=max(400, 320 * rows),
            showlegend=False,
            margin=dict(t=50, l=60, r=20, b=40),
            font=dict(family='Segoe UI, sans-serif', size=12),
        )
        self._store('qc_boxplots', fig, 'all')

    # ------------------------------------------------------------------
    # ANI chart
    # ------------------------------------------------------------------
    def _render_ani_chart(self, qc_df: pd.DataFrame):
        if qc_df is None or qc_df.empty:
            return
        ani_col = next((c for c in qc_df.columns if 'ani' in c.lower()
                       and pd.api.types.is_numeric_dtype(qc_df[c])), None)
        if not ani_col:
            return
        df = qc_df.dropna(subset=[ani_col])
        if df.empty:
            return

        fig = go.Figure()
        fig.add_trace(go.Histogram(
            x=df[ani_col],
            nbinsx=20,
            marker=dict(color=PALETTE['primary'],
                       line=dict(color='#0a0e1a', width=1)),
            hovertemplate=f'{ani_col} %{{x}}<br>%{{y}} samples<extra></extra>',
            name='ANI',
        ))
        fig.add_vline(
            x=95, line_dash='dash', line_color=PALETTE['danger'],
            annotation_text='95% species threshold',
            annotation_position='top right',
            annotation_font=dict(size=11, color=PALETTE['danger']),
        )
        fig.update_layout(
            template='plotly_dark',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=400, showlegend=False,
            margin=dict(t=50, l=60, r=20, b=50),
            xaxis=dict(title=ani_col, gridcolor='rgba(148,163,184,0.12)',
                      title_font=dict(size=12)),
            yaxis=dict(title='Samples', gridcolor='rgba(148,163,184,0.12)',
                      title_font=dict(size=12)),
            font=dict(family='Segoe UI, sans-serif', size=12),
        )
        self._store('ani_chart', fig, 'all')

    # ------------------------------------------------------------------
    # Summary narrative (compact top-of-page)
    # ------------------------------------------------------------------
    def _narrate(self, typing_df: pd.DataFrame,
                 gene_tables: Dict[str, pd.DataFrame],
                 qc_df: pd.DataFrame):
        n = len(typing_df)
        if n == 0:
            return

        top_st = '—'
        if 'ST' in typing_df.columns:
            vc = typing_df['ST'][~typing_df['ST'].isin(['Not Assigned', 'ND', ''])]
            if not vc.empty:
                top_st = f"{vc.value_counts().idxmax()} (n={vc.value_counts().max()})"

        top_patho = 'none'
        if 'Pathotype' in typing_df.columns:
            vc = typing_df['Pathotype'].value_counts()
            if not vc.empty:
                top_patho = f"{vc.idxmax()} (n={vc.max()})"

        n_amr = len(gene_tables.get('amr', pd.DataFrame())) if gene_tables else 0
        n_vir = len(gene_tables.get('virulence', pd.DataFrame())) if gene_tables else 0
        stec_count = int(typing_df.get('Pathotype', '').str.contains('STEC|EHEC',
                                                                      na=False).sum()) \
            if 'Pathotype' in typing_df.columns else 0

        facts = [
            f"<b>{n}</b> E. coli genomes analysed.",
            f"Top lineage: <b>{esc(top_st)}</b>.",
            f"Dominant pathotype: <b>{esc(top_patho)}</b>.",
            f"<b>{stec_count}</b> STEC/EHEC isolates detected.",
            f"<b>{n_amr}</b> AMR genes · <b>{n_vir}</b> virulence genes.",
        ]
        if qc_df is not None and not qc_df.empty:
            ani_col = next((c for c in qc_df.columns
                          if 'ani' in c.lower()), None)
            if ani_col:
                ani = qc_df[ani_col].dropna()
                if not ani.empty:
                    facts.append(
                        f"Mean ANI: <b>{ani.mean():.2f}%</b> "
                        f"(min {ani.min():.2f}%)."
                    )
        self.summary_facts = facts

    # ------------------------------------------------------------------
    # Sample table
    # ------------------------------------------------------------------
    def _build_sample_table(self, typing_df: pd.DataFrame) -> str:
        cols = ['Sample'] + [c for c in
                             ('ST', 'Serotype', 'Phylogroup', 'Pathotype',
                              'CH_Type', 'O_Type', 'H_Type')
                             if c in typing_df.columns]
        rows = ''
        for _, r in typing_df.iterrows():
            cells = ''
            for c in cols:
                val = r.get(c, '')
                if c == 'Pathotype' and val and val != 'none':
                    cells += (f'<td><span class="pathotype-pill">{esc(val)}</span></td>')
                else:
                    cells += f'<td>{esc(val)}</td>'
            rows += f'<tr data-sample="{esc(r["Sample"])}">{cells}</tr>'
        header = ''.join(f'<th>{esc(c)}</th>' for c in cols)
        return f'''<table class="sample-table" id="sample-table">
            <thead><tr>{header}</tr></thead>
            <tbody>{rows}</tbody>
        </table>'''

    # ------------------------------------------------------------------
    # HTML assembly
    # ------------------------------------------------------------------
    def _assemble_html(self, typing_df: pd.DataFrame,
                       qc_df: pd.DataFrame) -> str:
        has_qc = qc_df is not None and not qc_df.empty
        return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>EcoliTyper Interactive Dashboard</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js" charset="utf-8"></script>
<style>{self._css()}</style>
</head>
<body>
{self._header()}
{self._narrative_block()}
{self._tabs(has_qc)}
<main class="dashboard-main">
{self._tab_overview()}
{self._tab_typing()}
{self._tab_pathotype()}
{self._tab_qc() if has_qc else ''}
{self._tab_amr()}
{self._tab_virulence()}
{self._tab_plasmids()}
{self._tab_bacmet()}
{self._tab_mutations()}
{self._tab_mge()}
{self._tab_alerts()}
{self._tab_compare(typing_df)}
</main>
<footer class="dashboard-footer">
    <p><strong>EcoliTyper Interactive Dashboard v{TOOL_VERSION}</strong> · Generated {self.meta.get("generated", "")}</p>
    <p>Author: Brown Beckley · University of Ghana Medical School · brownbeckley94@gmail.com</p>
    <p>⭐ Star us on GitHub if you find this useful!</p>
</footer>
<div id="filter-popup" class="filter-popup" style="display:none;">
    <div class="filter-popup-header">
        <span id="filter-title">Samples</span>
        <button id="filter-close" class="filter-close">×</button>
    </div>
    <div id="filter-body" class="filter-popup-body"></div>
</div>
<script>{self._js()}</script>
</body>
</html>'''

    def _header(self) -> str:
        m = self.meta
        return f'''<header class="dashboard-header">
    <div class="header-content">
        <h1>🧬 EcoliTyper Interactive Dashboard</h1>
        <p class="subtitle">Modern visual analytics for <em>Escherichia coli</em> genomics</p>
        <div class="header-stats">
            <div class="header-stat"><div class="num">{m.get('n_samples', 0)}</div><div class="lbl">Samples</div></div>
            <div class="header-stat"><div class="num">{m.get('n_pathotypes', 0)}</div><div class="lbl">With Pathotype</div></div>
            <div class="header-stat"><div class="num">{m.get('n_stec', 0) + m.get('n_ehec', 0)}</div><div class="lbl">STEC/EHEC</div></div>
            <div class="header-stat"><div class="num">{m.get('n_alerts', 0)}</div><div class="lbl">Alerts</div></div>
        </div>
    </div>
</header>'''

    def _narrative_block(self) -> str:
        if not self.summary_facts:
            return ''
        items = ''.join(f'<li>{n}</li>' for n in self.summary_facts)
        return f'''<section class="narrative-bar">
    <h3>📖 What This Dataset Shows</h3>
    <ul>{items}</ul>
</section>'''

    def _tabs(self, has_qc: bool) -> str:
        tabs = [('overview', '📊 Overview'), ('typing', '🧬 Typing'),
                ('pathotype', '🦠 Pathotype')]
        if has_qc:
            tabs.append(('qc', '📏 QC'))
        tabs += [
            ('amr', '💊 AMR'), ('virulence', '⚔️ Virulence'),
            ('plasmids', '🧬 Plasmids'), ('bacmet', '🛡️ BACMET'),
            ('mutations', '🔬 Mutations'), ('mge', '📱 MGE'),
            ('alerts', '🚨 Alerts'), ('compare', '⚖️ Compare'),
        ]
        buttons = ''.join(
            f'<button class="tab-btn{" active" if i == 0 else ""}" data-tab="{t}">{l}</button>'
            for i, (t, l) in enumerate(tabs)
        )
        return f'''<nav class="tab-nav">{buttons}</nav>
<div class="view-filter">
    <label>Filter typing charts:</label>
    <button class="view-btn active" data-view="all">All samples</button>
    <button class="view-btn" data-view="stec">STEC / EHEC only</button>
    <button class="view-btn" data-view="esbl">ESBL producers only</button>
</div>'''

    def _tab_overview(self) -> str:
        return f'''<section class="tab-content active" data-tab-content="overview">
    <div class="grid-2">
        {self._chart_block('pathotype_donut', 'Pathotype Distribution')}
        {self._chart_block('phylogroup_donut', 'Phylogroup Distribution')}
    </div>
    <div class="grid-full">
        {self._chart_block('sunburst', 'Hierarchical Lineage — Phylogroup → ST → Serotype')}
    </div>
</section>'''

    def _tab_typing(self) -> str:
        return f'''<section class="tab-content" data-tab-content="typing">
    <div class="grid-2">
        {self._chart_block('bar-ST', 'MLST Sequence Types (Top 25)')}
        {self._chart_block('bar-Serotype', 'Serotypes (Top 25)')}
        {self._chart_block('bar-Phylogroup', 'Phylogroups')}
        {self._chart_block('bar-CH_Type', 'CH Types (Top 25)')}
    </div>
</section>'''

    def _tab_pathotype(self) -> str:
        return f'''<section class="tab-content" data-tab-content="pathotype">
    <div class="grid-2">
        {self._chart_block('bar-Pathotype', 'Pathotype Distribution')}
        {self._chart_block('virulence_by_patho', 'Virulence Genes by Pathotype')}
    </div>
</section>'''

    def _tab_qc(self) -> str:
        return f'''<section class="tab-content" data-tab-content="qc">
    <div class="grid-2">
        {self._chart_block('qc_boxplots', 'Assembly Quality Metrics')}
        {self._chart_block('ani_chart', 'Species Confirmation (fastANI)')}
    </div>
</section>'''

    def _tab_amr(self) -> str:
        return f'''<section class="tab-content" data-tab-content="amr">
    <div class="grid-full">
        {self._chart_block('amr_top', 'Top 30 AMR Genes')}
    </div>
    <div class="grid-full">
        {self._chart_block('amr_matrix', 'AMR Gene Presence Matrix')}
    </div>
    <div class="grid-full">
        {self._chart_block('amr_by_phylo', 'AMR Gene Prevalence by Phylogroup')}
    </div>
</section>'''

    def _tab_virulence(self) -> str:
        return f'''<section class="tab-content" data-tab-content="virulence">
    <div class="grid-full">
        {self._chart_block('virulence_top', 'Top 25 Virulence Genes')}
    </div>
</section>'''

    def _tab_plasmids(self) -> str:
        return f'''<section class="tab-content" data-tab-content="plasmids">
    <div class="grid-full">
        {self._chart_block('plasmid_top', 'Top 20 Plasmid Markers')}
    </div>
</section>'''

    def _tab_bacmet(self) -> str:
        return f'''<section class="tab-content" data-tab-content="bacmet">
    <div class="grid-full">
        {self._chart_block('bacmet_top', 'Top 25 BACMET Genes (Biocide / Heavy-Metal Resistance)')}
    </div>
</section>'''

    def _tab_mutations(self) -> str:
        return f'''<section class="tab-content" data-tab-content="mutations">
    <div class="grid-full">
        {self._chart_block('mutations_top', 'Top 30 Point Mutations')}
    </div>
    <div class="grid-full">
        {self._chart_block('mutations_matrix', 'Mutation Presence Matrix')}
    </div>
</section>'''

    def _tab_mge(self) -> str:
        return f'''<section class="tab-content" data-tab-content="mge">
    <div class="grid-full">
        {self._chart_block('mge_breakdown', 'Mean mobileOG Hits per Category')}
    </div>
</section>'''

    def _tab_alerts(self) -> str:
        return f'''<section class="tab-content" data-tab-content="alerts">
    {self._render_alerts()}
</section>'''

    def _render_alerts(self) -> str:
        if not self.alerts:
            return '<div class="no-alerts">✅ No alerts raised — the dataset is clean.</div>'
        groups = defaultdict(list)
        for a in self.alerts:
            groups[a['id']].append(a)
        html = '<div class="alert-panel">'
        for rule_id, items in groups.items():
            rule = ALERT_RULES[rule_id]
            lines = '<br>'.join(
                f'<b>{esc(a["sample"])}</b> — {esc(a["detail"])}'
                for a in items[:15]
            )
            if len(items) > 15:
                lines += f'<br><em>… and {len(items) - 15} more</em>'
            html += f'''<div class="alert-card {rule['severity']}">
                <div class="alert-icon">{rule['icon']}</div>
                <div class="alert-body">
                    <h4>{rule['title']} <span class="severity-tag {rule['severity']}">{rule['severity']}</span></h4>
                    <p class="desc">{rule['description']}</p>
                    <div class="detail">{lines}</div>
                </div>
            </div>'''
        html += '</div>'
        return html

    def _tab_compare(self, typing_df: pd.DataFrame) -> str:
        options = ''.join(
            f'<option value="{esc(s)}">{esc(s)}</option>'
            for s in sorted(typing_df['Sample'].tolist())
        ) if 'Sample' in typing_df.columns else ''
        return f'''<section class="tab-content" data-tab-content="compare">
    <div class="compare-controls">
        <label>Sample A:
            <select id="compare-a">{options}</select>
        </label>
        <label>Sample B:
            <select id="compare-b">{options}</select>
        </label>
        <button id="compare-run" class="btn-primary">Compare</button>
    </div>
    <div id="compare-output"></div>
    <details class="sample-table-details" open>
        <summary>📋 All Samples</summary>
        <input type="text" id="sample-search" class="search-input"
               placeholder="🔍 Search by sample name, ST, pathotype, serotype...">
        {self.sample_table_html}
    </details>
</section>'''

    # ------------------------------------------------------------------
    # CSS
    # ------------------------------------------------------------------
    def _css(self) -> str:
        return """
:root{
  --bg:#0a0e1a; --bg-card:#1a1f2e; --bg-hover:#252b3d; --bg-soft:#141922;
  --text:#e2e8f0; --text-muted:#94a3b8; --border:#2d3748;
  --primary:#3b82f6; --accent:#8b5cf6;
  --danger:#dc2626; --warning:#f59e0b; --success:#10b981;
  --pathotype:#ec4899;
}
*{margin:0;padding:0;box-sizing:border-box}
html,body{background:var(--bg);color:var(--text);
  font-family:'Segoe UI',Tahoma,Geneva,Verdana,sans-serif;line-height:1.55}
body{padding-bottom:40px}

.dashboard-header{background:linear-gradient(135deg,#1e3a8a 0%,#3b82f6 50%,#8b5cf6 100%);
  color:white;padding:40px 30px;border-radius:0 0 24px 24px;
  box-shadow:0 10px 40px rgba(0,0,0,.5);margin-bottom:24px}
.header-content{max-width:1600px;margin:0 auto}
.dashboard-header h1{font-size:2.5em;margin-bottom:8px;letter-spacing:-.5px;
  text-shadow:0 2px 8px rgba(0,0,0,.3)}
.dashboard-header .subtitle{font-size:1.08em;opacity:.94;margin-bottom:20px}
.header-stats{display:flex;gap:16px;flex-wrap:wrap;margin-top:20px}
.header-stat{background:rgba(255,255,255,.15);padding:14px 22px;border-radius:12px;
  backdrop-filter:blur(10px);border:1px solid rgba(255,255,255,.2);min-width:120px}
.header-stat .num{font-size:1.9em;font-weight:700;line-height:1}
.header-stat .lbl{font-size:.78em;opacity:.9;text-transform:uppercase;
  letter-spacing:.5px;margin-top:4px}

.narrative-bar{max-width:1600px;margin:0 auto 24px;padding:20px 30px;
  background:var(--bg-card);border-radius:14px;border-left:5px solid var(--primary)}
.narrative-bar h3{font-size:1em;color:var(--primary);margin-bottom:10px;
  text-transform:uppercase;letter-spacing:.5px;font-weight:700}
.narrative-bar ul{list-style:none;padding:0;display:flex;flex-wrap:wrap;gap:14px 28px}
.narrative-bar li{font-size:1.02em;color:var(--text)}
.narrative-bar li b{color:var(--primary);font-weight:700}

.tab-nav{position:sticky;top:0;z-index:100;background:rgba(10,14,26,.96);
  backdrop-filter:blur(12px);padding:12px 30px;display:flex;gap:8px;flex-wrap:wrap;
  border-bottom:1px solid var(--border);margin-bottom:16px}
.tab-btn{padding:10px 20px;background:transparent;border:1px solid var(--border);
  color:var(--text-muted);border-radius:8px;cursor:pointer;font-weight:600;
  font-size:.9em;transition:all .2s}
.tab-btn:hover{background:var(--bg-hover);color:var(--text)}
.tab-btn.active{background:var(--primary);color:white;border-color:var(--primary);
  box-shadow:0 0 20px rgba(59,130,246,.4)}

.view-filter{max-width:1600px;margin:0 auto 24px;padding:0 30px;display:flex;
  gap:10px;align-items:center;flex-wrap:wrap}
.view-filter label{color:var(--text-muted);font-size:.9em;font-weight:600}
.view-btn{padding:8px 18px;background:var(--bg-card);border:1px solid var(--border);
  color:var(--text-muted);border-radius:20px;cursor:pointer;font-weight:600;
  font-size:.85em;transition:all .2s}
.view-btn:hover{color:var(--text);border-color:var(--primary)}
.view-btn.active{background:var(--primary);border-color:var(--primary);color:white}

.dashboard-main{max-width:1600px;margin:0 auto;padding:0 30px}
.tab-content{display:none;animation:fadeIn .3s}
.tab-content.active{display:block}
@keyframes fadeIn{from{opacity:0;transform:translateY(8px)}
                  to{opacity:1;transform:translateY(0)}}

.grid-2{display:grid;grid-template-columns:repeat(auto-fit,minmax(520px,1fr));
  gap:20px;margin-bottom:20px}
.grid-full{margin-bottom:20px}
.chart-card{background:var(--bg-card);border:1px solid var(--border);
  border-radius:14px;padding:20px;transition:box-shadow .2s}
.chart-card:hover{box-shadow:0 8px 30px rgba(0,0,0,.35)}
.chart-title{font-size:1.05em;margin-bottom:6px;font-weight:600;
  display:flex;align-items:center;gap:8px;color:var(--text)}
.chart-title::before{content:'';width:4px;height:18px;background:var(--primary);
  border-radius:2px}
.chart-subtitle{font-size:.85em;color:var(--text-muted);margin-bottom:12px}

.alert-panel{display:flex;flex-direction:column;gap:12px}
.alert-card{background:var(--bg-card);border-radius:12px;padding:18px 22px;
  border-left:5px solid;display:flex;gap:16px;align-items:flex-start;
  transition:transform .15s}
.alert-card:hover{transform:translateX(4px)}
.alert-card.critical{border-left-color:var(--danger)}
.alert-card.high{border-left-color:var(--warning)}
.alert-card.medium{border-left-color:#eab308}
.alert-icon{font-size:1.8em;line-height:1}
.alert-body{flex:1}
.alert-body h4{font-size:1.02em;margin-bottom:4px}
.alert-body .desc{font-size:.88em;color:var(--text-muted);margin-bottom:8px}
.alert-body .detail{font-size:.85em;font-family:'Consolas',monospace;
  background:var(--bg);padding:10px 14px;border-radius:6px;line-height:1.7}
.severity-tag{font-size:.72em;padding:3px 10px;border-radius:12px;
  text-transform:uppercase;letter-spacing:.5px;font-weight:700;margin-left:8px}
.severity-tag.critical{background:rgba(220,38,38,.2);color:#fca5a5}
.severity-tag.high{background:rgba(245,158,11,.2);color:#fcd34d}
.severity-tag.medium{background:rgba(234,179,8,.2);color:#fde047}
.no-alerts{background:rgba(16,185,129,.1);border-left:5px solid var(--success);
  padding:40px;border-radius:12px;text-align:center;font-size:1.1em;color:var(--success)}

.compare-controls{background:var(--bg-card);border-radius:14px;padding:20px;
  display:flex;gap:16px;align-items:center;flex-wrap:wrap;margin-bottom:20px}
.compare-controls label{display:flex;gap:8px;align-items:center;
  color:var(--text-muted);font-weight:600;font-size:.92em}
.compare-controls select{background:var(--bg);color:var(--text);
  border:1px solid var(--border);border-radius:8px;padding:8px 12px;
  font-size:.92em;min-width:240px}
.btn-primary{background:var(--primary);color:white;border:none;
  padding:10px 22px;border-radius:8px;font-weight:700;cursor:pointer;
  transition:transform .2s}
.btn-primary:hover{transform:translateY(-2px);box-shadow:0 4px 14px rgba(59,130,246,.4)}

#compare-output{display:grid;grid-template-columns:1fr 1fr;gap:20px;margin-bottom:30px}
.compare-col{background:var(--bg-card);border-radius:14px;padding:20px;
  border:1px solid var(--border)}
.compare-col h4{color:var(--primary);margin-bottom:12px;font-size:1.05em}
.compare-row{display:flex;justify-content:space-between;padding:8px 0;
  border-bottom:1px dashed var(--border);font-size:.92em}
.compare-row.diff{background:rgba(220,38,38,.08);border-left:3px solid var(--danger);
  padding-left:8px;margin-left:-8px;border-radius:4px}
.compare-verdict{grid-column:1/-1;padding:16px 20px;border-radius:10px;
  margin-bottom:12px;display:flex;gap:12px;align-items:center;font-size:.95em}
.compare-verdict.danger{background:rgba(220,38,38,.15);border-left:5px solid var(--danger)}
.compare-verdict.warning{background:rgba(245,158,11,.15);border-left:5px solid var(--warning)}
.compare-verdict.info{background:rgba(59,130,246,.15);border-left:5px solid var(--primary)}

.sample-table-details{margin-top:20px;background:var(--bg-card);
  border-radius:14px;padding:20px}
.sample-table-details summary{cursor:pointer;font-weight:700;
  color:var(--primary);font-size:1.02em;padding:6px 0}
.search-input{width:100%;padding:10px 14px;margin:12px 0;background:var(--bg);
  border:1px solid var(--border);border-radius:8px;color:var(--text);
  font-size:.92em}
.search-input:focus{outline:none;border-color:var(--primary);
  box-shadow:0 0 0 3px rgba(59,130,246,.15)}
.sample-table{width:100%;border-collapse:collapse;margin-top:12px;font-size:.88em}
.sample-table th{background:var(--bg);color:var(--text);padding:10px 12px;
  text-align:left;position:sticky;top:0;font-weight:600;
  border-bottom:2px solid var(--primary)}
.sample-table td{padding:8px 12px;border-bottom:1px solid var(--border)}
.sample-table tr:hover{background:var(--bg-hover)}
.pathotype-pill{background:rgba(236,72,153,.2);color:#f9a8d4;
  padding:2px 10px;border-radius:10px;font-weight:600;font-size:.85em}

.filter-popup{position:fixed;bottom:20px;right:20px;width:440px;max-height:60vh;
  background:var(--bg-card);border:1px solid var(--border);border-radius:12px;
  box-shadow:0 20px 60px rgba(0,0,0,.7);z-index:1000;overflow:hidden;
  display:flex;flex-direction:column}
.filter-popup-header{background:var(--primary);color:white;padding:12px 16px;
  display:flex;justify-content:space-between;align-items:center;font-weight:700}
.filter-close{background:transparent;border:none;color:white;font-size:1.5em;
  cursor:pointer;line-height:1}
.filter-popup-body{padding:14px 16px;overflow-y:auto;
  font-size:.85em;font-family:'Consolas',monospace;line-height:1.7}
.filter-popup-body span{display:inline-block;background:var(--bg);
  padding:3px 10px;border-radius:12px;margin:2px;border:1px solid var(--border)}

.dashboard-footer{max-width:1600px;margin:40px auto 0;padding:30px;
  text-align:center;color:var(--text-muted);font-size:.88em;
  border-top:1px solid var(--border)}
.dashboard-footer p{margin:4px 0}

@media (max-width:900px){
  .dashboard-header h1{font-size:1.7em}
  .grid-2{grid-template-columns:1fr}
  #compare-output{grid-template-columns:1fr}
  .dashboard-main,.tab-nav,.view-filter,.narrative-bar{
    padding-left:15px;padding-right:15px}
}
"""

    # ------------------------------------------------------------------
    # JS
    # ------------------------------------------------------------------
    def _js(self) -> str:
        return r"""
// ---- Tab switching ----
document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => {
        const tab = btn.dataset.tab;
        document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
        document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
        btn.classList.add('active');
        const content = document.querySelector(`[data-tab-content="${tab}"]`);
        if (content) content.classList.add('active');
        setTimeout(() => {
            document.querySelectorAll('.js-plotly-plot').forEach(div => {
                if (div.offsetWidth > 0) Plotly.Plots.resize(div);
            });
        }, 60);
    });
});

// ---- View filter (typing charts only) ----
document.querySelectorAll('.view-btn').forEach(btn => {
    btn.addEventListener('click', () => {
        const view = btn.dataset.view;
        document.querySelectorAll('.view-btn').forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        document.querySelectorAll('.chart-view').forEach(v => {
            v.style.display = v.dataset.view === view ? 'block' : 'none';
        });
        setTimeout(() => {
            document.querySelectorAll('.chart-view:not([style*="display: none"]) .js-plotly-plot')
                .forEach(div => { Plotly.Plots.resize(div); });
        }, 60);
    });
});

// ---- Cross-filter popup ----
const popup = document.getElementById('filter-popup');
const popupTitle = document.getElementById('filter-title');
const popupBody = document.getElementById('filter-body');
if (popup) {
    document.getElementById('filter-close').addEventListener('click', () => {
        popup.style.display = 'none';
    });
}

document.addEventListener('DOMContentLoaded', () => {
    setTimeout(() => {
        document.querySelectorAll('.js-plotly-plot').forEach(div => {
            if (typeof div.on !== 'function') return;
            div.on('plotly_click', (data) => {
                const pt = data.points[0];
                let samples = [];
                if (typeof pt.customdata === 'string') {
                    samples = pt.customdata.split(';').filter(s => s);
                } else if (Array.isArray(pt.customdata) && pt.customdata.length) {
                    samples = [String(pt.customdata[0])];
                }
                if (!samples.length) return;
                const label = pt.label || pt.y || pt.x || 'Samples';
                popupTitle.textContent = `${label} · ${samples.length} sample(s)`;
                popupBody.innerHTML = samples
                    .map(s => `<span>${s}</span>`).join('');
                popup.style.display = 'flex';
            });
        });
    }, 500);
});

// ---- Compare tool ----
const compareBtn = document.getElementById('compare-run');
if (compareBtn) {
    compareBtn.addEventListener('click', () => {
        const a = document.getElementById('compare-a').value;
        const b = document.getElementById('compare-b').value;
        if (!a || !b) return;
        if (a === b) {
            document.getElementById('compare-output').innerHTML =
                '<div class="compare-verdict warning" style="grid-column:1/-1;">Pick two different samples.</div>';
            return;
        }
        const rows = Array.from(document.querySelectorAll('#sample-table tr'));
        const rowA = rows.find(r => r.dataset.sample === a);
        const rowB = rows.find(r => r.dataset.sample === b);
        if (!rowA || !rowB) return;
        const headers = Array.from(
            rowA.parentElement.parentElement.querySelectorAll('thead th')
        ).map(th => th.textContent);
        const aVals = Array.from(rowA.querySelectorAll('td')).map(td => td.textContent);
        const bVals = Array.from(rowB.querySelectorAll('td')).map(td => td.textContent);

        let leftHtml = `<h4>${a}</h4>`;
        let rightHtml = `<h4>${b}</h4>`;
        let diffCount = 0;
        for (let i = 1; i < headers.length; i++) {
            const diff = aVals[i] !== bVals[i];
            if (diff) diffCount++;
            leftHtml += `<div class="compare-row ${diff ? 'diff' : ''}">` +
                        `<span>${headers[i]}</span><span>${aVals[i]}</span></div>`;
            rightHtml += `<div class="compare-row ${diff ? 'diff' : ''}">` +
                         `<span>${headers[i]}</span><span>${bVals[i]}</span></div>`;
        }

        let verdict = '';
        if (diffCount === 0) {
            verdict = `<div class="compare-verdict danger">
                🚨 <strong>Identical typing profile</strong> — possible transmission
                pair or duplicate sample. Investigate further.
            </div>`;
        } else if (diffCount <= 2) {
            verdict = `<div class="compare-verdict warning">
                ⚠️ <strong>${diffCount} typing differences</strong> — likely
                closely related but epidemiologically distinct.
            </div>`;
        } else {
            verdict = `<div class="compare-verdict info">
                🔍 <strong>${diffCount} typing differences</strong> — distinct lineages.
            </div>`;
        }

        document.getElementById('compare-output').innerHTML =
            verdict +
            `<div class="compare-col">${leftHtml}</div>` +
            `<div class="compare-col">${rightHtml}</div>`;
    });
}

// ---- Sample table search ----
const searchInput = document.getElementById('sample-search');
if (searchInput) {
    searchInput.addEventListener('input', () => {
        const filter = searchInput.value.toUpperCase();
        document.querySelectorAll('#sample-table tbody tr').forEach(row => {
            const text = row.textContent.toUpperCase();
            row.style.display = text.includes(filter) ? '' : 'none';
        });
    });
}
"""


# ==============================================================================
# STATIC PLOTTER (publication-quality matplotlib — bug-fixed)
# ==============================================================================
class StaticPlotter:
    """Publication-quality matplotlib figures for download/sharing."""

    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.subdirs = {
            'png': self.output_dir / 'PNG',
            'pdf': self.output_dir / 'PDF',
            'svg': self.output_dir / 'SVG',
            'data': self.output_dir / 'DATA',
        }
        for d in self.subdirs.values():
            d.mkdir(parents=True, exist_ok=True)

    def _save(self, fig, name: str):
        for fmt in ('png', 'pdf', 'svg'):
            try:
                fig.savefig(self.subdirs[fmt] / f"{name}.{fmt}",
                            dpi=300 if fmt == 'png' else None,
                            bbox_inches='tight')
            except Exception as e:
                print(f"    ⚠️ Save {name}.{fmt} failed: {e}")
        plt.close(fig)

    @staticmethod
    def _wrap_labels(labels, max_len=30):
        return [l if len(str(l)) <= max_len
                else str(l)[:max_len - 1] + '…' for l in labels]

    def typing_summary(self, typing_df: pd.DataFrame):
        """4-panel summary: pathotype donut, ST bar, phylogroup bar, serotype bar."""
        if typing_df.empty:
            return
        fig, axes = plt.subplots(2, 2, figsize=(18, 14))

        # Pathotype donut
        ax = axes[0, 0]
        if 'Pathotype' in typing_df.columns:
            counts = typing_df['Pathotype'].value_counts()
            if not counts.empty:
                colors = [CHART_COLORS[i % len(CHART_COLORS)]
                          for i in range(len(counts))]
                ax.pie(counts.values, labels=counts.index, colors=colors,
                       autopct='%1.1f%%', startangle=90,
                       wedgeprops=dict(width=0.4, edgecolor='white'))
                ax.set_title('Pathotype Distribution', fontweight='bold', fontsize=13)
                pd.DataFrame({'Pathotype': counts.index, 'Count': counts.values}) \
                    .to_csv(self.subdirs['data'] / 'pathotype_distribution.csv',
                           index=False)

        # ST bar
        ax = axes[0, 1]
        if 'ST' in typing_df.columns:
            valid = typing_df[~typing_df['ST'].isin(['Not Assigned', 'ND', ''])]
            vc = valid['ST'].value_counts().head(15)
            if not vc.empty:
                labels = self._wrap_labels(vc.index.tolist())
                ax.barh(labels, vc.values, color=PALETTE['st'])
                ax.invert_yaxis()
                ax.set_xlabel('Samples')
                ax.set_title(f'MLST ST — Top {len(vc)}', fontweight='bold', fontsize=13)
                ax.grid(alpha=0.3, axis='x')
                ax.tick_params(axis='y', labelsize=9)

        # Phylogroup
        ax = axes[1, 0]
        if 'Phylogroup' in typing_df.columns:
            valid = typing_df[~typing_df['Phylogroup'].isin(
                ['Not Assigned', 'ND', '', 'nan'])]
            if not valid.empty:
                vc = valid['Phylogroup'].value_counts()
                ax.bar(vc.index.astype(str), vc.values,
                       color=PALETTE['phylogroup'])
                ax.set_ylabel('Samples')
                ax.set_title('Phylogroups', fontweight='bold', fontsize=13)
                ax.grid(alpha=0.3, axis='y')

        # Serotype bar
        ax = axes[1, 1]
        if 'Serotype' in typing_df.columns:
            valid = typing_df[~typing_df['Serotype'].isin(
                ['Not Assigned', 'ND', '', 'nan'])]
            vc = valid['Serotype'].value_counts().head(15)
            if not vc.empty:
                labels = self._wrap_labels(vc.index.tolist(), max_len=24)
                ax.barh(labels, vc.values, color=PALETTE['serotype'])
                ax.invert_yaxis()
                ax.set_xlabel('Samples')
                ax.set_title(f'Serotypes — Top {len(vc)}', fontweight='bold', fontsize=13)
                ax.grid(alpha=0.3, axis='x')
                ax.tick_params(axis='y', labelsize=9)

        plt.suptitle('E. coli Typing Summary', fontsize=16, fontweight='bold')
        plt.tight_layout()
        self._save(fig, 'typing_summary')
        print("    ✓ typing_summary.[png/pdf/svg]")

    def gene_top_plot(self, gene_df: pd.DataFrame, title: str,
                      filename: str, color: str):
        if gene_df is None or gene_df.empty:
            return
        top = gene_df.nlargest(20, 'Count')
        if top.empty:
            return
        labels = self._wrap_labels(top['Gene'].astype(str), max_len=40)
        height = max(6, 0.4 * len(top))
        fig, ax = plt.subplots(figsize=(11, height))
        ax.barh(labels, top['Count'], color=color)
        ax.invert_yaxis()
        ax.set_xlabel('Number of genomes')
        ax.set_title(title, fontweight='bold', fontsize=13)
        ax.grid(alpha=0.3, axis='x')
        ax.tick_params(axis='y', labelsize=9)
        plt.tight_layout()
        self._save(fig, filename)
        print(f"    ✓ {filename}.[png/pdf/svg]")

    def qc_boxplots(self, qc_df: pd.DataFrame):
        """QC box plots — bug-fixed version that works with matplotlib 3.9+."""
        if qc_df is None or qc_df.empty:
            return
        metrics = [c for c in ('N50', 'Total Length', 'GC Content (%)',
                               'Total Sequences', 'Ambiguous Bases (%)')
                   if c in qc_df.columns
                   and pd.api.types.is_numeric_dtype(qc_df[c])]
        if not metrics:
            return
        n = len(metrics)
        cols = min(3, n)
        rows = (n + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 5 * rows))
        if n == 1:
            axes = np.array([axes])
        axes = axes.flatten()

        for i, col in enumerate(metrics):
            data = qc_df[col].dropna()
            if data.empty:
                continue
            ax = axes[i]
            # FIX: no `labels=` kwarg — set ticks manually. Works in all
            # matplotlib versions (3.5 → 4.0).
            bp = ax.boxplot(
                data.values,
                vert=True,
                patch_artist=True,
                boxprops=dict(facecolor=PALETTE['primary'], alpha=0.6),
                medianprops=dict(color='black', linewidth=2),
            )
            ax.scatter(np.random.normal(1, 0.04, len(data)), data.values,
                       alpha=0.6, s=20, color='#334155')
            ax.set_ylabel(col, fontsize=11)
            ax.set_title(col, fontweight='bold', fontsize=12)
            ax.grid(alpha=0.3, axis='y')

        for j in range(len(metrics), len(axes)):
            axes[j].axis('off')

        plt.suptitle('Assembly QC Metrics', fontsize=15, fontweight='bold')
        plt.tight_layout()
        self._save(fig, 'qc_boxplots')
        print("    ✓ qc_boxplots.[png/pdf/svg]")


# ==============================================================================
# ORCHESTRATOR
# ==============================================================================
class EcoliVisualizer:
    """Coordinates data loading, dashboard assembly, and static exports."""

    def __init__(self, input_dir: Path, output_dir: Optional[Path] = None):
        self.input_dir = Path(input_dir)
        self.output_dir = (Path(output_dir) if output_dir
                          else self.input_dir / 'ECOLI_VISUALIZATIONS')
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.loader = EcoliDataLoader(self.input_dir)
        self.plotter = StaticPlotter(self.output_dir)
        self.builder = DashboardBuilder(self.output_dir)

        self.typing_df = pd.DataFrame()
        self.gene_tables: Dict[str, pd.DataFrame] = {}
        self.mutations_df = pd.DataFrame()
        self.mge_df = pd.DataFrame()
        self.qc_df = pd.DataFrame()

    def run(self):
        print("=" * 70)
        print(f"🧬 ECOLITYPER VISUALIZER v{TOOL_VERSION}")
        print("=" * 70)
        print(f"Input:  {self.input_dir}")
        print(f"Output: {self.output_dir}")
        print("-" * 70)

        start = datetime.now()

        print("\n📥 Loading data...")
        self.typing_df = self.loader.load_typing()
        self.gene_tables = self.loader.load_gene_tables()
        self.mutations_df = self.loader.load_mutations()
        self.mge_df = self.loader.load_mge()
        self.qc_df = self.loader.load_qc()

        if self.typing_df.empty and not self.gene_tables:
            print("❌ No data found. Nothing to visualize.")
            return

        print("\n📊 Generating static publication plots...")
        try:
            self.plotter.typing_summary(self.typing_df)
        except Exception as e:
            print(f"  ⚠️ typing_summary failed: {e}")
        try:
            self.plotter.gene_top_plot(
                self.gene_tables.get('amr'), 'Top 20 AMR Genes',
                'amr_top_genes', PALETTE['amr'])
        except Exception as e:
            print(f"  ⚠️ amr_top_genes failed: {e}")
        try:
            self.plotter.gene_top_plot(
                self.gene_tables.get('virulence'), 'Top 20 Virulence Genes',
                'virulence_top_genes', PALETTE['virulence'])
        except Exception as e:
            print(f"  ⚠️ virulence_top_genes failed: {e}")
        try:
            self.plotter.qc_boxplots(self.qc_df)
        except Exception as e:
            print(f"  ⚠️ qc_boxplots failed: {e}")

        if PLOTLY_AVAILABLE:
            print("\n🌐 Building interactive dashboard...")
            try:
                self.builder.build(
                    self.typing_df, self.gene_tables,
                    self.mutations_df, self.mge_df, self.qc_df)
            except Exception as e:
                print(f"  ⚠️ Dashboard failed: {e}")
                import traceback
                traceback.print_exc()
        else:
            print("\n⚠️ Plotly not available — skipping interactive dashboard")
            print("   Install with: pip install 'plotly>=5.0'")

        self._write_summary_report(datetime.now() - start)
        zip_path = self._build_export_bundle()
        if zip_path:
            print(f"\n📦 Export bundle: {zip_path.name}")

        print("\n" + "=" * 70)
        print("✅ VISUALIZATION PIPELINE COMPLETE")
        print("=" * 70)

    def _write_summary_report(self, duration):
        report = self.output_dir / 'ecolityper_visualization_report.txt'
        with open(report, 'w', encoding='utf-8') as f:
            f.write("=" * 70 + "\n")
            f.write("ECOLITYPER VISUALIZATION REPORT\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Generated: {datetime.now().isoformat()}\n")
            f.write(f"Duration:  {duration}\n")
            f.write(f"Input:     {self.input_dir}\n")
            f.write(f"Output:    {self.output_dir}\n\n")

            f.write("DATA SUMMARY\n" + "-" * 40 + "\n")
            f.write(f"Samples:  {len(self.typing_df)}\n")
            for cat, df in self.gene_tables.items():
                f.write(f"{cat.capitalize()}: {len(df)} unique genes\n")
            if not self.mutations_df.empty:
                f.write(f"Mutations: {len(self.mutations_df)} unique entries\n")
            if not self.mge_df.empty:
                f.write(f"MGE: {len(self.mge_df)} samples profiled\n")
            if not self.qc_df.empty:
                f.write(f"QC: {len(self.qc_df)} samples with metrics\n")

            if self.builder.alerts:
                f.write("\nALERTS\n" + "-" * 40 + "\n")
                for a in self.builder.alerts[:50]:
                    severity = ALERT_RULES[a['id']]['severity'].upper()
                    f.write(f"  [{severity}] {a['sample']}: {a['detail']}\n")

            f.write("\nOUTPUT FILES\n" + "-" * 40 + "\n")
            for d in ('PNG', 'PDF', 'SVG', 'DATA'):
                sub = self.output_dir / d
                if sub.exists():
                    files = list(sub.glob("*"))
                    f.write(f"\n{d} ({len(files)} files):\n")
                    for fp in sorted(files):
                        f.write(f"  • {fp.name}\n")

            dashboard = self.output_dir / 'ecolityper_dashboard.html'
            if dashboard.exists():
                f.write("\nINTERACTIVE DASHBOARD\n" + "-" * 40 + "\n")
                f.write(f"  {dashboard.name}\n")
                f.write(f"  Open in browser: file://{dashboard.resolve()}\n")

        print(f"\n📋 Report: {report.name}")

    def _build_export_bundle(self) -> Optional[Path]:
        try:
            zip_path = self.output_dir.parent / 'ecolityper_visualizations_bundle.zip'
            with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
                for file in self.output_dir.rglob('*'):
                    if file.is_file() and file != zip_path:
                        zf.write(file, file.relative_to(self.output_dir.parent))
            return zip_path
        except Exception as e:
            print(f"  ⚠️ Export bundle failed: {e}")
            return None


# ==============================================================================
# CLI
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="EcoliTyper Visualizer — unified interactive + static report",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python visualization_reporter.py
  python visualization_reporter.py --input /path/to/results --output /path/to/vis
  python visualization_reporter.py --no-static      # dashboard only
  python visualization_reporter.py --no-dashboard   # static plots only
        """,
    )
    parser.add_argument('--input', '-i', type=str, default='.',
                        help='Directory containing EcoliTyper outputs (default: cwd)')
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='Output directory (default: ECOLI_VISUALIZATIONS)')
    parser.add_argument('--no-static', action='store_true',
                        help='Skip static matplotlib exports')
    parser.add_argument('--no-dashboard', action='store_true',
                        help='Skip interactive dashboard')
    args = parser.parse_args()

    viz = EcoliVisualizer(
        input_dir=Path(args.input),
        output_dir=Path(args.output) if args.output else None,
    )
    if args.no_static:
        viz.plotter.typing_summary = lambda *a, **k: None
        viz.plotter.gene_top_plot = lambda *a, **k: None
        viz.plotter.qc_boxplots = lambda *a, **k: None
    if args.no_dashboard:
        viz.builder.build = lambda *a, **k: None

    try:
        viz.run()
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()