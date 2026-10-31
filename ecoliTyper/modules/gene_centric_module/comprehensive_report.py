#!/usr/bin/env python3
"""
EcoliTyper Comprehensive Report Generator
==========================================
Creates a unified, interactive report from CHTyper, serotype, MLST and phylogroup results.

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-10-05 (Robust: any subset of input files is accepted)

Recognized input files (in target directory):
  • chtyper_results.tsv
  • serotype_analysis_report.tsv
  • mlst_summary.tsv
  • phylogrouping_results.tsv
"""

import os
import sys
import json
import random
import argparse
import pandas as pd
from pathlib import Path
from collections import Counter
from datetime import datetime


# =========================================================================
# Configuration
# =========================================================================
SCIENCE_QUOTES = [
    {"text": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
    {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
    {"text": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
    {"text": "In science, the credit goes to the man who convinces the world, not the man to whom the idea first occurs.", "author": "Francis Darwin"},
    {"text": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"},
    {"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
    {"text": "EcoliTyper turns raw E. coli genomes into structured, decision-ready typing information.", "author": "Brown Beckley"},
    {"text": "Every sequence is a story. EcoliTyper helps you read the chapters on virulence, resistance and lineage.", "author": "Brown Beckley"},
    {"text": "Thank you for using EcoliTyper — please share any feature suggestions you have in mind.", "author": "Brown Beckley"},
    {"text": "Behind every ST, every serotype, every CH type, there is an epidemiological story worth telling.", "author": "Brown Beckley"},
]

ASCII_ART = r"""
███████╗ ██████╗ ██████╗ ██╗     ██╗████████╗██╗   ██╗██████╗ ███████╗██████╗ 
██╔════╝██╔════╝██╔═══██╗██║     ██║╚══██╔══╝╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
█████╗  ██║     ██║   ██║██║     ██║   ██║    ╚████╔╝ ██████╔╝█████╗  ██████╔╝
██╔══╝  ██║     ██║   ██║██║     ██║   ██║     ╚██╔╝  ██╔═══╝ ██╔══╝  ██╔══██╗
███████╗╚██████╗╚██████╔╝███████╗██║   ██║      ██║   ██║     ███████╗██║  ██║
╚══════╝ ╚═════╝ ╚═════╝ ╚══════╝╚═╝   ╚═╝      ╚═╝   ╚═╝     ╚══════╝╚═╝  ╚═╝
"""

NULL_VALUES = {'-', 'ND', '', 'UNKNOWN', 'Unknown', 'Not Assigned', 'None',
               'nan', 'NaN', 'NA', 'N/A', 'null', 'NULL'}

INPUT_FILES = {
    'chtyper':    'chtyper_results.tsv',
    'serotype':   'serotype_analysis_report.tsv',
    'mlst':       'mlst_summary.tsv',
    'phylogroup': 'phylogrouping_results.tsv',
}


# =========================================================================
# Helpers
# =========================================================================
def normalize_sample_name(name):
    """Strip FASTA extensions for cross-module matching."""
    if not name:
        return ""
    s = str(name).strip().strip('\t \n\r')
    for ext in ('.fna', '.fasta', '.fa', '.fsa', '.fn'):
        if s.lower().endswith(ext):
            return s[:-len(ext)]
    return s


def _clean(value, default="Not Assigned"):
    if value is None:
        return default
    v = str(value).strip()
    if not v or v in NULL_VALUES:
        return default
    return v


def _safe_read_tsv(path):
    try:
        return pd.read_csv(path, sep='\t', dtype=str).fillna('')
    except Exception as e:
        print(f"   ✗ Error reading {path}: {e}")
        return pd.DataFrame()


def _to_float(value, default=0.0):
    try:
        return float(str(value).replace('%', '').strip())
    except Exception:
        return default


def _coverage_pct(cov):
    try:
        if '/' in str(cov):
            num, denom = str(cov).split('/')
            d = float(denom)
            return (float(num) / d * 100) if d else 0.0
    except Exception:
        pass
    return 0.0


# =========================================================================
# Loaders — each returns {} if file missing / malformed
# =========================================================================
def load_chtyper_data(path):
    print(f"📊 Loading CHTyper data from: {path}")
    if not Path(path).exists():
        print("   ⚠️  File not found — CH type columns will be 'Not Assigned'")
        return {}
    df = _safe_read_tsv(path)
    if df.empty:
        return {}
    out = {}
    for _, row in df.iterrows():
        sample = normalize_sample_name(row.get('Sample_ID', row.get('sample', '')))
        if not sample:
            continue
        out[sample] = {
            'fumc_type':     _clean(row.get('FumC_Type')),
            'fumc_identity': _clean(row.get('FumC_Identity'), "0.00"),
            'fumc_coverage': _clean(row.get('FumC_Coverage'), "0/0"),
            'fimh_type':     _clean(row.get('FimH_Type')),
            'fimh_identity': _clean(row.get('FimH_Identity'), "0.00"),
            'fimh_coverage': _clean(row.get('FimH_Coverage'), "0/0"),
            'chtyper_status': _clean(row.get('Status'), "Unknown"),
        }
    print(f"   ✓ Loaded {len(out)} samples")
    return out


def load_serotype_data(path):
    print(f"📊 Loading serotype data from: {path}")
    if not Path(path).exists():
        print("   ⚠️  File not found — serotype columns will be 'Not Assigned'")
        return {}
    df = _safe_read_tsv(path)
    if df.empty:
        return {}
    out = {}
    for _, row in df.iterrows():
        sample = normalize_sample_name(row.get('Sample_ID', row.get('sample', '')))
        if not sample:
            continue
        out[sample] = {
            'serotype':   _clean(row.get('Serotype')),
            'o_type':     _clean(row.get('O_Type')),
            'h_type':     _clean(row.get('H_Type')),
            'sero_genes': _clean(row.get('Genes_Found')),
            'sero_confidence': _clean(row.get('Confidence'), "Unknown"),
            'sero_status': _clean(row.get('Status'), "Unknown"),
        }
    print(f"   ✓ Loaded {len(out)} samples")
    return out


def load_mlst_data(path):
    print(f"📊 Loading MLST data from: {path}")
    if not Path(path).exists():
        print("   ⚠️  File not found — MLST columns will be 'Not Assigned'")
        return {}
    df = _safe_read_tsv(path)
    if df.empty:
        return {}
    out = {}
    for _, row in df.iterrows():
        sample = normalize_sample_name(row.get('Sample', row.get('Sample_ID', '')))
        if not sample:
            continue
        st = _clean(row.get('ST'))
        out[sample] = {
            'mlst_st': _clean(f"ST{st}" if st != "Not Assigned" else "Not Assigned"),
            'mlst_profile': _clean(row.get('Allele_Profile')),
        }
    print(f"   ✓ Loaded {len(out)} samples")
    return out


def load_phylogroup_data(path):
    print(f"📊 Loading phylogroup data from: {path}")
    if not Path(path).exists():
        print("   ⚠️  File not found — phylogroup column will be 'Not Assigned'")
        return {}
    df = _safe_read_tsv(path)
    if df.empty:
        return {}
    out = {}
    for _, row in df.iterrows():
        sample = normalize_sample_name(row.get('Sample_ID', row.get('sample', '')))
        if not sample:
            continue
        out[sample] = _clean(row.get('Clermont_Type'))
    print(f"   ✓ Loaded {len(out)} samples")
    return out


# =========================================================================
# Combine
# =========================================================================
def combine_data(chtyper, serotype, mlst, phylogroup):
    print("\n🔗 Combining data from all sources...")
    all_samples = set()
    for d in (chtyper, serotype, mlst, phylogroup):
        all_samples.update(d.keys())
    if not all_samples:
        print("   ⚠️  No samples found in any source")
        return []

    combined = []
    for sample in sorted(all_samples):
        ch = chtyper.get(sample, {}) or {}
        se = serotype.get(sample, {}) or {}
        ml = mlst.get(sample, {}) or {}

        fumc = ch.get('fumc_type', "Not Assigned")
        fimh = ch.get('fimh_type', "Not Assigned")
        if fumc != "Not Assigned" and fimh != "Not Assigned":
            ch_type = f"{fumc}-{fimh}"
        else:
            ch_type = "Not Assigned"

        # Overall sample status: prefer CHTyper status, else serotype, else Unknown
        status = ch.get('chtyper_status') or se.get('sero_status') or "Unknown"

        combined.append({
            'sample': sample,
            'phylogroup': phylogroup.get(sample, "Not Assigned"),
            'serotype': se.get('serotype', "Not Assigned"),
            'o_type': se.get('o_type', "Not Assigned"),
            'h_type': se.get('h_type', "Not Assigned"),
            'sero_confidence': se.get('sero_confidence', "Unknown"),
            'mlst_st': ml.get('mlst_st', "Not Assigned"),
            'mlst_profile': ml.get('mlst_profile', "Not Assigned"),
            'fumc_type': fumc,
            'fimh_type': fimh,
            'ch_type': ch_type,
            'fumc_identity': ch.get('fumc_identity', "0.00"),
            'fumc_coverage': ch.get('fumc_coverage', "0/0"),
            'fimh_identity': ch.get('fimh_identity', "0.00"),
            'fimh_coverage': ch.get('fimh_coverage', "0/0"),
            'status': status,
        })

    print(f"   ✓ Combined data for {len(combined)} samples")
    return combined


# =========================================================================
# Distributions
# =========================================================================
DIST_FIELDS = {
    'phylogroup': ('Phylogroup (Clermont)', '#667eea', '#4c1d95'),
    'serotype':   ('Serotype',              '#8b5cf6', '#6d28d9'),
    'o_type':     ('O Type',                '#06b6d4', '#0e7490'),
    'h_type':     ('H Type',                '#14b8a6', '#0f766e'),
    'mlst_st':    ('MLST Sequence Type',    '#3b82f6', '#1e40af'),
    'fumc_type':  ('FumC Type',             '#f59e0b', '#b45309'),
    'fimh_type':  ('FimH Type',             '#ec4899', '#be185d'),
    'ch_type':    ('CH Type (FumC-FimH)',   '#a855f7', '#7e22ce'),
}


def compute_distributions(data):
    distributions = {}
    for key, (label, top, bottom) in DIST_FIELDS.items():
        counter = Counter()
        not_typed = 0
        for row in data:
            v = row.get(key, "Not Assigned")
            if v in ("Not Assigned", "Unknown", "None"):
                not_typed += 1
            else:
                counter[v] += 1
        ordered = sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))
        distributions[key] = {
            'label': label,
            'unique_count': len(ordered),
            'counts': [{'value': v, 'count': c} for v, c in ordered],
            'not_typed': not_typed,
            'color_top': top,
            'color_bottom': bottom,
        }
    return distributions


def compute_quality_metrics(data):
    fumc_id, fimh_id, fumc_cov, fimh_cov = [], [], [], []
    high = medium = low = 0

    for row in data:
        f_id = _to_float(row['fumc_identity'])
        h_id = _to_float(row['fimh_identity'])
        f_cov = _coverage_pct(row['fumc_coverage'])
        h_cov = _coverage_pct(row['fimh_coverage'])

        if row['fumc_type'] != "Not Assigned":
            fumc_id.append(f_id); fumc_cov.append(f_cov)
        if row['fimh_type'] != "Not Assigned":
            fimh_id.append(h_id); fimh_cov.append(h_cov)

        if f_id >= 99 and h_id >= 99:
            high += 1
        elif f_id >= 90 and h_id >= 90:
            medium += 1
        else:
            low += 1

    avg = lambda x: round(sum(x) / len(x), 2) if x else 0.0

    return {
        'avg_fumc_identity': avg(fumc_id),
        'avg_fimh_identity': avg(fimh_id),
        'avg_fumc_coverage': avg(fumc_cov),
        'avg_fimh_coverage': avg(fimh_cov),
        'high_confidence':   high,
        'medium_confidence': medium,
        'low_confidence':    low,
        'n_fumc_typed': len(fumc_id),
        'n_fimh_typed': len(fimh_id),
    }


# =========================================================================
# TSV
# =========================================================================
TSV_COLUMNS = [
    ('sample',         'Sample_ID'),
    ('phylogroup',     'Phylogroup'),
    ('serotype',       'Serotype'),
    ('o_type',         'O_Type'),
    ('h_type',         'H_Type'),
    ('sero_confidence','Serotype_Confidence'),
    ('mlst_st',        'MLST_ST'),
    ('mlst_profile',   'MLST_Allele_Profile'),
    ('fumc_type',      'FumC_Type'),
    ('fumc_identity',  'FumC_Identity'),
    ('fumc_coverage',  'FumC_Coverage'),
    ('fimh_type',      'FimH_Type'),
    ('fimh_identity',  'FimH_Identity'),
    ('fimh_coverage',  'FimH_Coverage'),
    ('ch_type',        'CH_Type'),
    ('status',         'Status'),
]


def generate_tsv_report(data, output_file):
    print(f"\n📄 Generating TSV report: {output_file}")
    try:
        df = pd.DataFrame(data)
        df = df[[k for k, _ in TSV_COLUMNS]]
        df.columns = [v for _, v in TSV_COLUMNS]
        df.to_csv(output_file, sep='\t', index=False)
        print("   ✓ TSV report generated")
        return True
    except Exception as e:
        print(f"   ✗ Error generating TSV: {e}")
        return False


# =========================================================================
# Excel (optional)
# =========================================================================
def generate_xlsx_report(data, distributions, quality, output_file):
    print(f"📄 Generating Excel report: {output_file}")
    try:
        # Merged samples
        df = pd.DataFrame(data)[[k for k, _ in TSV_COLUMNS]]
        df.columns = [v for _, v in TSV_COLUMNS]

        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Samples', index=False)

            # Distributions sheet
            rows = []
            for key, d in distributions.items():
                for c in d['counts']:
                    rows.append({'Field': d['label'], 'Value': c['value'], 'Count': c['count']})
            pd.DataFrame(rows).to_excel(writer, sheet_name='Distributions', index=False)

            # Quality sheet
            pd.DataFrame([quality]).T.reset_index().rename(
                columns={'index': 'Metric', 0: 'Value'}
            ).to_excel(writer, sheet_name='Quality', index=False)

        print("   ✓ Excel report generated")
        return True
    except ImportError:
        print("   ⚠️  openpyxl not installed — skipping Excel export (pip install openpyxl)")
        return False
    except Exception as e:
        print(f"   ✗ Error generating Excel: {e}")
        return False


# =========================================================================
# JSON
# =========================================================================
def generate_json_report(data, distributions, quality, missing_sources, output_file):
    print(f"📄 Generating JSON report: {output_file}")
    try:
        total = len(data)
        completed = sum(1 for d in data if d['status'] == 'Completed')

        report = {
            'metadata': {
                'generated_date': datetime.now().isoformat(),
                'total_samples': total,
                'report_type': 'EcoliTyper Comprehensive Report',
                'version': '1.0.0',
                'author': 'Brown Beckley',
                'email': 'brownbeckley94@gmail.com',
                'affiliation': 'University of Ghana Medical School - Department of Medical Biochemistry',
                'missing_sources': missing_sources,
            },
            'analysis_summary': {
                'total_samples': total,
                'completed': completed,
                'failed': total - completed,
                'success_rate': round(completed / total * 100, 2) if total else 0,
                'unique_phylogroups': distributions['phylogroup']['unique_count'],
                'unique_serotypes': distributions['serotype']['unique_count'],
                'unique_mlst_sts': distributions['mlst_st']['unique_count'],
                'unique_ch_types': distributions['ch_type']['unique_count'],
                'most_common_phylogroup': distributions['phylogroup']['counts'][0]['value']
                    if distributions['phylogroup']['counts'] else 'None',
                'most_common_serotype': distributions['serotype']['counts'][0]['value']
                    if distributions['serotype']['counts'] else 'None',
                'most_common_st': distributions['mlst_st']['counts'][0]['value']
                    if distributions['mlst_st']['counts'] else 'None',
                'most_common_ch_type': distributions['ch_type']['counts'][0]['value']
                    if distributions['ch_type']['counts'] else 'None',
            },
            'quality_metrics': quality,
            'distributions': {k: {kk: vv for kk, vv in v.items() if kk not in ('color_top', 'color_bottom')}
                              for k, v in distributions.items()},
            'methodology': {
                'description': 'Unified E. coli typing report — CHTyper (FumC/FimH), Serotype, MLST, Phylogroup',
                'modules': ['CHTyper', 'SerotypeFinder', 'MLST', 'ClermontTyping'],
                'purpose': 'Comprehensive E. coli characterisation for epidemiological and clinical surveillance',
            },
            'samples': data,
        }

        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        print("   ✓ JSON report generated")
        return True
    except Exception as e:
        print(f"   ✗ Error generating JSON: {e}")
        return False


# =========================================================================
# HTML
# =========================================================================
def _build_distribution_cards(distributions):
    cards = []
    for key in ('phylogroup', 'serotype', 'o_type', 'h_type',
                'mlst_st', 'fumc_type', 'fimh_type', 'ch_type'):
        d = distributions.get(key)
        if not d:
            continue
        top = d['color_top']
        bottom = d['color_bottom']
        total_typed = sum(c['count'] for c in d['counts'])
        max_count = max((c['count'] for c in d['counts']), default=1)

        if not d['counts']:
            bars = '<div class="dist-empty">No typed samples</div>'
        else:
            bars = ''
            for c in d['counts']:
                pct = (c['count'] / total_typed * 100) if total_typed else 0
                width = (c['count'] / max_count * 100) if max_count else 0
                bars += f'''
                <div class="dist-row">
                    <div class="dist-label" title="{c['value']}">{c['value']}</div>
                    <div class="dist-bar-track">
                        <div class="dist-bar-fill"
                             style="width:{width:.1f}%;background:linear-gradient(90deg,{top} 0%,{bottom} 100%);">
                            <span class="dist-bar-count">{c['count']} ({pct:.1f}%)</span>
                        </div>
                    </div>
                </div>'''

        footnote = ''
        if d['not_typed'] > 0:
            footnote = f'<div class="dist-footnote">⚠️ Not Assigned / unknown: <strong>{d["not_typed"]}</strong></div>'

        cards.append(f'''
        <div class="dist-card">
            <div class="dist-header" style="background:linear-gradient(135deg,{top} 0%,{bottom} 100%);">
                <div class="dist-title">{d['label']}</div>
                <div class="dist-unique">
                    <span class="dist-unique-num">{d['unique_count']}</span>
                    <span class="dist-unique-lbl">unique</span>
                </div>
            </div>
            <div class="dist-body">
                {bars}
                {footnote}
            </div>
        </div>''')

    return '<div class="dist-grid">' + ''.join(cards) + '</div>'


def generate_html_report(data, distributions, quality, missing_sources, output_file):
    print(f"📄 Generating HTML report: {output_file}")
    try:
        random_quote = random.choice(SCIENCE_QUOTES)
        total = len(data)
        completed = sum(1 for d in data if d['status'] == 'Completed')
        failed = total - completed

        # Table rows
        table_rows = ""
        for row in data:
            status = row['status']
            status_cls = 'status-ok' if status == 'Completed' else 'status-err'

            table_rows += f'''
            <tr data-phylogroup="{row['phylogroup']}"
                data-serotype="{row['serotype']}"
                data-st="{row['mlst_st']}"
                data-chtype="{row['ch_type']}"
                data-status="{status}">
                <td class="col-check"><input type="checkbox" class="row-check"></td>
                <td class="col-sample"><strong>{row['sample']}</strong></td>
                <td><span class="pill pill-phylo">{row['phylogroup']}</span></td>
                <td>{row['serotype']}</td>
                <td>{row['o_type']}</td>
                <td>{row['h_type']}</td>
                <td><span class="pill pill-st">{row['mlst_st']}</span></td>
                <td>{row['fumc_type']}</td>
                <td>{row['fimh_type']}</td>
                <td><strong class="ch-type">{row['ch_type']}</strong></td>
                <td>{row['fumc_identity']}</td>
                <td>{row['fimh_identity']}</td>
                <td class="{status_cls}">{status}</td>
            </tr>'''

        dist_cards_html = _build_distribution_cards(distributions)

        # Overview stat cards
        overview_stats = [
            ("Total Samples",  total,                                              ""),
            ("Completed",      completed,                                          "stat-ok"),
            ("Failed",         failed,                                             "stat-err"),
            ("Phylogroups",    distributions['phylogroup']['unique_count'],        ""),
            ("Serotypes",      distributions['serotype']['unique_count'],          ""),
            ("MLST STs",       distributions['mlst_st']['unique_count'],           ""),
            ("CH Types",       distributions['ch_type']['unique_count'],           ""),
            ("FumC Types",     distributions['fumc_type']['unique_count'],         ""),
            ("FimH Types",     distributions['fimh_type']['unique_count'],         ""),
        ]
        overview_cards = ''.join(
            f'<div class="stat-card {cls}"><div class="stat-value">{val}</div>'
            f'<div class="stat-label">{lbl}</div></div>'
            for lbl, val, cls in overview_stats
        )

        # Quality cards
        q = quality
        quality_cards = f'''
        <div class="quality-grid">
            <div class="quality-card"><div class="q-val">{q['avg_fumc_identity']}%</div><div class="q-lbl">Avg FumC Identity</div></div>
            <div class="quality-card"><div class="q-val">{q['avg_fimh_identity']}%</div><div class="q-lbl">Avg FimH Identity</div></div>
            <div class="quality-card"><div class="q-val">{q['avg_fumc_coverage']}%</div><div class="q-lbl">Avg FumC Coverage</div></div>
            <div class="quality-card"><div class="q-val">{q['avg_fimh_coverage']}%</div><div class="q-lbl">Avg FimH Coverage</div></div>
            <div class="quality-card q-high"><div class="q-val">{q['high_confidence']}</div><div class="q-lbl">High Confidence (≥99%)</div></div>
            <div class="quality-card q-med"><div class="q-val">{q['medium_confidence']}</div><div class="q-lbl">Medium (90-99%)</div></div>
            <div class="quality-card q-low"><div class="q-val">{q['low_confidence']}</div><div class="q-lbl">Low (&lt;90%)</div></div>
        </div>'''

        # Missing-sources banner
        missing_banner = ''
        if missing_sources:
            items = ''.join(f'<li><code>{fname}</code> — <em>{key}</em></li>'
                            for key, fname in missing_sources)
            missing_banner = f'''
            <div class="missing-banner">
                <div class="missing-banner-title">⚠️ Some input files were not found — those columns show <code>Not Assigned</code></div>
                <ul class="missing-banner-list">{items}</ul>
                <div class="missing-banner-hint">Run the corresponding module then re-run this report.</div>
            </div>'''

        # Unique phylogroups / STs / CH types for filter pills
        phylos = [c['value'] for c in distributions['phylogroup']['counts']]
        phylo_pills = '<button class="filter-pill active" data-filter="all">All</button>' + ''.join(
            f'<button class="filter-pill" data-filter="{p}">{p}</button>' for p in phylos[:8]
        )

        sero_opts = ''.join(f'<option value="{c["value"]}">{c["value"]} ({c["count"]})</option>'
                            for c in distributions['serotype']['counts'])
        st_opts = ''.join(f'<option value="{c["value"]}">{c["value"]} ({c["count"]})</option>'
                          for c in distributions['mlst_st']['counts'])
        ch_opts = ''.join(f'<option value="{c["value"]}">{c["value"]} ({c["count"]})</option>'
                          for c in distributions['ch_type']['counts'])

        # ----- Build JS block (kept separate to avoid brace escaping headaches)
        quotes_json = json.dumps(SCIENCE_QUOTES)
        all_data_json = json.dumps(data)
        headers_json = json.dumps([v for _, v in TSV_COLUMNS])

        js = """
const QUOTES = __QUOTES__;
const ALL_DATA = __ALL_DATA__;
const COLUMN_LABELS = __HEADERS__;

// ---------- Rotating quote ----------
(function(){
  let i = 0;
  const qText   = document.getElementById('quoteText');
  const qAuthor = document.getElementById('quoteAuthor');
  function rotate() {
    const q = QUOTES[i % QUOTES.length];
    qText.textContent   = '"' + q.text + '"';
    qAuthor.textContent = '— ' + q.author;
    i++;
  }
  rotate();
  setInterval(rotate, 10000);
})();

// ---------- Tabs ----------
document.querySelectorAll('.tab-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
    document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
    btn.classList.add('active');
    document.getElementById(btn.dataset.tab).classList.add('active');
  });
});

// ---------- Table refs ----------
const table        = document.getElementById('mainTable');
const tbody        = table.querySelector('tbody');
const allRows      = Array.from(tbody.querySelectorAll('tr'));
const searchInput  = document.getElementById('tableSearch');
const counter      = document.getElementById('resultsCounter');
const selCount     = document.getElementById('selCount');
const phyloPills   = document.querySelectorAll('.filter-pill[data-filter]');
const seroSelect   = document.getElementById('seroFilter');
const stSelect     = document.getElementById('stFilter');
const chSelect     = document.getElementById('chFilter');
const statusSelect = document.getElementById('statusFilter');
const selectAllCb  = document.getElementById('selectAll');

let activePhylo = 'all';

function applyFilters() {
  const term = (searchInput.value || '').toLowerCase().trim();
  const sero = seroSelect.value;
  const st   = stSelect.value;
  const ch   = chSelect.value;
  const stat = statusSelect.value;

  let visible = 0;
  allRows.forEach(r => {
    const text = r.textContent.toLowerCase();
    const okText  = !term || text.includes(term);
    const okPhylo = activePhylo === 'all' || r.dataset.phylogroup === activePhylo;
    const okSero  = !sero || r.dataset.serotype === sero;
    const okSt    = !st   || r.dataset.st === st;
    const okCh    = !ch   || r.dataset.chtype === ch;
    const okStat  = !stat || r.dataset.status === stat;
    const show = okText && okPhylo && okSero && okSt && okCh && okStat;
    r.style.display = show ? '' : 'none';
    if (show) visible++;
  });
  counter.textContent = visible + ' of ' + allRows.length + ' shown';
  updateSelCount();
}

function updateSelCount() {
  const n = allRows.filter(r => r.querySelector('.row-check').checked && r.style.display !== 'none').length;
  selCount.textContent = n + ' selected';
}

searchInput.addEventListener('input', applyFilters);
[seroSelect, stSelect, chSelect, statusSelect].forEach(s => s.addEventListener('change', applyFilters));

phyloPills.forEach(p => {
  p.addEventListener('click', () => {
    phyloPills.forEach(x => x.classList.remove('active'));
    p.classList.add('active');
    activePhylo = p.dataset.filter;
    applyFilters();
  });
});

// ---------- Row selection ----------
allRows.forEach(r => r.querySelector('.row-check').addEventListener('change', updateSelCount));
selectAllCb.addEventListener('change', () => {
  allRows.forEach(r => { if (r.style.display !== 'none') r.querySelector('.row-check').checked = selectAllCb.checked; });
  updateSelCount();
});

// ---------- Sorting ----------
const headers = table.querySelectorAll('th');
headers.forEach((th, idx) => {
  if (th.classList.contains('col-check') || th.classList.contains('col-nosort')) return;
  th.addEventListener('click', () => {
    const asc = !th.classList.contains('sorted-asc');
    headers.forEach(h => h.classList.remove('sorted-asc', 'sorted-desc'));
    const rows = allRows.slice().sort((a, b) => {
      const aT = a.children[idx].textContent.trim();
      const bT = b.children[idx].textContent.trim();
      const aN = parseFloat(aT.replace(/[^0-9.\\-]/g, ''));
      const bN = parseFloat(bT.replace(/[^0-9.\\-]/g, ''));
      if (!isNaN(aN) && !isNaN(bN)) return asc ? aN - bN : bN - aN;
      return asc ? aT.localeCompare(bT) : bT.localeCompare(aT);
    });
    rows.forEach(r => tbody.appendChild(r));
    th.classList.add(asc ? 'sorted-asc' : 'sorted-desc');
  });
});

// ---------- Column visibility ----------
const colToggleBtn = document.getElementById('colToggleBtn');
const colPanel     = document.getElementById('colPanel');
colToggleBtn.addEventListener('click', e => { e.stopPropagation(); colPanel.classList.toggle('open'); });
document.addEventListener('click', () => colPanel.classList.remove('open'));

headers.forEach((th, idx) => {
  if (th.classList.contains('col-check')) return;
  const cb = document.createElement('input');
  cb.type = 'checkbox'; cb.checked = true; cb.dataset.colIdx = idx;
  const lbl = document.createElement('label');
  lbl.appendChild(cb); lbl.appendChild(document.createTextNode(' ' + th.textContent.trim()));
  cb.addEventListener('change', () => {
    const vis = cb.checked;
    table.querySelectorAll('tr').forEach(row => {
      const cell = row.children[idx];
      if (cell) cell.style.display = vis ? '' : 'none';
    });
  });
  colPanel.appendChild(lbl);
});

// ---------- Export helpers ----------
function getVisibleRows() {
  return allRows.filter(r => r.style.display !== 'none');
}
function getSelectedRows() {
  return allRows.filter(r => r.querySelector('.row-check').checked && r.style.display !== 'none');
}
function rowToObj(row) {
  const cells = Array.from(row.children);
  const obj = {};
  COLUMN_LABELS.forEach((lbl, i) => {
    const cell = cells[i + 1];               // skip checkbox column
    obj[lbl] = cell ? cell.textContent.trim() : '';
  });
  return obj;
}
function download(filename, content, mime) {
  const blob = new Blob([content], {type: mime + ';charset=utf-8;'});
  const url  = URL.createObjectURL(blob);
  const a    = document.createElement('a');
  a.href = url; a.download = filename; a.click();
  URL.revokeObjectURL(url);
}
function toCSV(objs) {
  if (!objs.length) return '';
  const headers = Object.keys(objs[0]);
  const esc = v => '"' + String(v == null ? '' : v).replace(/"/g, '""') + '"';
  return [headers.map(esc).join(',')]
    .concat(objs.map(o => headers.map(h => esc(o[h])).join(',')))
    .join('\\n');
}
function toTSV(objs) {
  if (!objs.length) return '';
  const headers = Object.keys(objs[0]);
  return [headers.join('\\t')]
    .concat(objs.map(o => headers.map(h => String(o[h] == null ? '' : o[h]).replace(/\\t/g, ' ')).join('\\t')))
    .join('\\n');
}

document.getElementById('exportCSV').addEventListener('click', () => {
  const rows = document.getElementById('exportScope').value === 'selected' ? getSelectedRows() : getVisibleRows();
  download('ecolityper_samples.csv', toCSV(rows.map(rowToObj)), 'text/csv');
});
document.getElementById('exportTSV').addEventListener('click', () => {
  const rows = document.getElementById('exportScope').value === 'selected' ? getSelectedRows() : getVisibleRows();
  download('ecolityper_samples.tsv', toTSV(rows.map(rowToObj)), 'text/tab-separated-values');
});
document.getElementById('exportJSON').addEventListener('click', () => {
  const rows = document.getElementById('exportScope').value === 'selected' ? getSelectedRows() : getVisibleRows();
  download('ecolityper_samples.json', JSON.stringify(rows.map(rowToObj), null, 2), 'application/json');
});
document.getElementById('exportXLSX').addEventListener('click', () => {
  if (typeof XLSX === 'undefined') { alert('Excel export requires internet connection to load SheetJS.'); return; }
  const rows = document.getElementById('exportScope').value === 'selected' ? getSelectedRows() : getVisibleRows();
  const ws = XLSX.utils.json_to_sheet(rows.map(rowToObj));
  const wb = XLSX.utils.book_new();
  XLSX.utils.book_append_sheet(wb, ws, 'EcoliTyper');
  XLSX.writeFile(wb, 'ecolityper_samples.xlsx');
});
document.getElementById('copyBtn').addEventListener('click', () => {
  const rows = document.getElementById('exportScope').value === 'selected' ? getSelectedRows() : getVisibleRows();
  navigator.clipboard.writeText(toTSV(rows.map(rowToObj))).then(() => {
    const b = document.getElementById('copyBtn');
    const old = b.textContent; b.textContent = '✓ Copied!';
    setTimeout(() => b.textContent = old, 1500);
  });
});

// ---------- SheetJS lazy-load ----------
(function loadSheetJS() {
  const s = document.createElement('script');
  s.src = 'https://cdn.sheetjs.com/xlsx-0.20.2/package/dist/xlsx.full.min.js';
  document.head.appendChild(s);
})();
""".replace("__QUOTES__", quotes_json) \
   .replace("__ALL_DATA__", all_data_json) \
   .replace("__HEADERS__", headers_json)

        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>EcoliTyper Comprehensive Report</title>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{
  font-family:'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  min-height:100vh; color:#1f2937; padding:20px;
}}
.container {{ max-width:1500px; margin:0 auto; }}

/* Header */
.ascii-container {{
  background: rgba(0,0,0,0.75); padding:18px; border-radius:15px;
  margin-bottom:18px; box-shadow: 0 8px 32px rgba(0,0,0,0.4);
  border: 2px solid rgba(139,92,246,0.35); text-align:center; overflow-x:auto;
}}
.ascii-art {{
  font-family:'Courier New', monospace; font-size:10px; line-height:1.1;
  white-space:pre; color:#00ff00;
  text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
}}
.title-block {{ text-align:center; margin-bottom:20px; color:#fff; }}
.title-block h1 {{ font-size:2.2em; text-shadow: 0 2px 10px rgba(0,0,0,0.3); }}
.title-block p  {{ opacity:0.9; margin-top:6px; font-size:1.05em; }}

/* Quote */
.quote-container {{
  background: rgba(255,255,255,0.12); backdrop-filter: blur(10px);
  padding:18px; border-radius:12px; text-align:center; margin-bottom:20px;
  color:#fff; border:1px solid rgba(255,255,255,0.2);
  min-height:90px; display:flex; flex-direction:column; justify-content:center;
  transition: opacity 0.5s;
}}
.quote-text {{ font-style:italic; font-size:16px; margin-bottom:6px; }}
.quote-author {{ color:#fde68a; font-weight:bold; font-size:13px; }}

/* Section shell */
.section {{
  background: rgba(255,255,255,0.97); border-radius:12px; padding:24px;
  margin-bottom:20px; box-shadow: 0 4px 20px rgba(0,0,0,0.1);
}}
.section h2 {{
  color:#4c1d95; border-bottom:3px solid #8b5cf6; padding-bottom:10px;
  margin-bottom:18px; font-size:22px;
}}
.section h3 {{ color:#5b21b6; margin:16px 0 10px; font-size:17px; }}

/* Missing banner */
.missing-banner {{
  background:#fef3c7; border-left:5px solid #f59e0b; border-radius:8px;
  padding:14px 18px; margin-bottom:20px; color:#78350f; font-size:13px;
}}
.missing-banner-title {{ font-weight:700; margin-bottom:6px; }}
.missing-banner-list {{ margin:6px 0 6px 22px; }}
.missing-banner-list code {{ background:rgba(120,53,15,0.15); padding:1px 6px; border-radius:4px; }}

/* Tabs */
.tabs {{ display:flex; gap:6px; margin-bottom:18px; flex-wrap:wrap; }}
.tab-btn {{
  background: rgba(255,255,255,0.35); color:#4c1d95; border:2px solid rgba(255,255,255,0.55);
  padding:10px 20px; border-radius:24px; cursor:pointer; font-weight:700; font-size:13px;
  transition: all 0.2s;
}}
.tab-btn:hover {{ background: rgba(255,255,255,0.6); }}
.tab-btn.active {{ background:#fff; color:#4c1d95; border-color:#fff; box-shadow: 0 4px 14px rgba(0,0,0,0.15); }}
.tab-panel {{ display:none; }}
.tab-panel.active {{ display:block; }}

/* Stats */
.stats-grid {{
  display:grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
  gap:14px; margin-bottom:18px;
}}
.stat-card {{
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  color:#fff; padding:16px; border-radius:10px; text-align:center;
  box-shadow: 0 4px 15px rgba(0,0,0,0.18); transition: transform 0.2s;
}}
.stat-card:hover {{ transform: translateY(-3px); }}
.stat-card.stat-ok  {{ background: linear-gradient(135deg, #22c55e 0%, #15803d 100%); }}
.stat-card.stat-err {{ background: linear-gradient(135deg, #ef4444 0%, #991b1b 100%); }}
.stat-value {{ font-size:26px; font-weight:800; }}
.stat-label {{ font-size:11px; opacity:0.92; text-transform:uppercase; letter-spacing:0.5px; margin-top:4px; }}

/* Toolbar */
.toolbar {{
  display:flex; flex-wrap:wrap; gap:10px; align-items:center;
  background:#f5f3ff; border:1px solid #ddd6fe; padding:14px;
  border-radius:10px; margin-bottom:14px;
}}
.search-wrap {{ position:relative; flex:1; min-width:240px; }}
.search-wrap input {{
  width:100%; padding:10px 14px 10px 38px; border:2px solid #c4b5fd;
  border-radius:8px; font-size:14px; outline:none;
}}
.search-wrap input:focus {{ border-color:#8b5cf6; box-shadow: 0 0 0 3px rgba(139,92,246,0.15); }}
.search-wrap::before {{ content:'🔍'; position:absolute; left:12px; top:50%; transform:translateY(-50%); }}

select.filter-select {{
  padding:9px 12px; border:2px solid #c4b5fd; border-radius:8px;
  background:#fff; font-size:13px; color:#4c1d95; cursor:pointer; font-weight:600;
}}

.filter-pills {{ display:flex; gap:6px; flex-wrap:wrap; }}
.filter-pill {{
  padding:7px 15px; border-radius:20px; border:2px solid #c4b5fd;
  background:#fff; color:#5b21b6; cursor:pointer; font-size:12.5px;
  font-weight:700; transition: all 0.2s;
}}
.filter-pill:hover {{ background:#ede9fe; }}
.filter-pill.active {{
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  color:#fff; border-color:transparent;
}}

/* Action buttons */
.action-btn {{
  padding:9px 15px; border:none; border-radius:8px; cursor:pointer;
  font-weight:700; font-size:13px; transition: all 0.2s; color:#fff;
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  box-shadow: 0 3px 10px rgba(102,126,234,0.35);
}}
.action-btn:hover {{ transform: translateY(-1px); box-shadow: 0 5px 14px rgba(102,126,234,0.45); }}
.action-btn.secondary {{ background: linear-gradient(135deg, #a78bfa 0%, #8b5cf6 100%); }}
.action-btn.success   {{ background: linear-gradient(135deg, #22c55e 0%, #15803d 100%); }}
.action-btn.info      {{ background: linear-gradient(135deg, #06b6d4 0%, #0e7490 100%); }}

/* Column toggle panel */
.col-toggle-wrap {{ position:relative; }}
.col-panel {{
  position:absolute; right:0; top:110%; background:#fff; border:1px solid #ddd6fe;
  border-radius:10px; padding:10px 14px; box-shadow: 0 10px 30px rgba(0,0,0,0.15);
  z-index:30; display:none; max-height:320px; overflow-y:auto; min-width:200px;
}}
.col-panel.open {{ display:block; }}
.col-panel label {{ display:block; font-size:13px; color:#4c1d95; padding:4px 0; cursor:pointer; }}

/* Export bar */
.export-bar {{ display:flex; gap:8px; align-items:center; flex-wrap:wrap; margin-left:auto; }}
.sel-counter {{
  font-size:12px; font-weight:700; color:#6d28d9;
  background:#ede9fe; padding:6px 12px; border-radius:20px;
}}

/* Table */
.table-container {{
  overflow:auto; max-height:640px; border-radius:10px;
  box-shadow: 0 2px 12px rgba(0,0,0,0.06); border:1px solid #e5e7eb;
}}
table {{ width:100%; border-collapse:collapse; font-size:13px; }}
th {{
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  color:#fff; padding:12px 10px; text-align:left; position:sticky; top:0; z-index:2;
  cursor:pointer; user-select:none; white-space:nowrap;
}}
th.col-check {{ cursor:default; }}
th:not(.col-check)::after {{ content:' ⇅'; opacity:0.5; font-size:11px; }}
th.sorted-asc::after  {{ content:' ▲'; opacity:1; }}
th.sorted-desc::after {{ content:' ▼'; opacity:1; }}
td {{ padding:10px; border-bottom:1px solid #f1f5f9; vertical-align:middle; }}
tr:nth-child(even) td {{ background:#faf5ff; }}
tr:hover td {{ background:#ede9fe; }}
tr.row-selected td {{ background:#ddd6fe !important; }}
.col-check {{ width:34px; text-align:center; }}
.col-sample {{ white-space:nowrap; }}
.ch-type {{ color:#6d28d9; }}

.pill {{
  display:inline-block; padding:3px 9px; border-radius:12px;
  font-size:11.5px; font-weight:700;
}}
.pill-phylo {{ background:#ede9fe; color:#6d28d9; }}
.pill-st    {{ background:#dbeafe; color:#1e40af; }}
.status-ok  {{ color:#15803d; font-weight:700; }}
.status-err {{ color:#b91c1c; font-weight:700; }}

/* Distribution cards */
.dist-grid {{
  display:grid; grid-template-columns: repeat(auto-fit, minmax(340px, 1fr)); gap:16px;
}}
.dist-card {{
  background:#fff; border-radius:10px; overflow:hidden; border:1px solid #e5e7eb;
  box-shadow: 0 4px 15px rgba(0,0,0,0.07); transition: transform 0.2s;
}}
.dist-card:hover {{ transform: translateY(-3px); box-shadow: 0 8px 24px rgba(0,0,0,0.12); }}
.dist-header {{
  display:flex; justify-content:space-between; align-items:center;
  padding:13px 16px; color:#fff;
}}
.dist-title {{ font-size:14px; font-weight:700; }}
.dist-unique {{ background: rgba(255,255,255,0.25); padding:3px 10px; border-radius:20px;
                display:flex; gap:4px; align-items:baseline; font-size:11px; }}
.dist-unique-num {{ font-size:15px; font-weight:800; }}
.dist-body {{ padding:14px 16px; }}
.dist-row {{ margin-bottom:9px; }}
.dist-label {{
  font-size:12px; font-weight:600; color:#334155; margin-bottom:4px;
  white-space:nowrap; overflow:hidden; text-overflow:ellipsis;
}}
.dist-bar-track {{ background:#e2e8f0; height:22px; border-radius:6px; overflow:hidden; }}
.dist-bar-fill {{
  height:100%; border-radius:6px; display:flex; align-items:center;
  justify-content:flex-end; padding-right:8px; min-width:44px; transition: width 0.7s ease;
}}
.dist-bar-count {{ font-size:11px; font-weight:700; color:#fff; text-shadow: 0 1px 2px rgba(0,0,0,0.35); }}
.dist-empty {{ color:#94a3b8; font-style:italic; font-size:13px; text-align:center; padding:12px 0; }}
.dist-footnote {{
  margin-top:10px; font-size:11px; color:#b45309; background:#fef3c7;
  padding:6px 10px; border-radius:6px; border-left:3px solid #f59e0b;
}}

/* Quality */
.quality-grid {{
  display:grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap:14px;
}}
.quality-card {{
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  color:#fff; padding:20px; border-radius:12px; text-align:center;
  box-shadow: 0 4px 15px rgba(0,0,0,0.18);
}}
.quality-card.q-high {{ background: linear-gradient(135deg, #22c55e 0%, #15803d 100%); }}
.quality-card.q-med  {{ background: linear-gradient(135deg, #f59e0b 0%, #b45309 100%); }}
.quality-card.q-low  {{ background: linear-gradient(135deg, #ef4444 0%, #991b1b 100%); }}
.q-val {{ font-size:28px; font-weight:800; }}
.q-lbl {{ font-size:12px; opacity:0.92; margin-top:4px; }}

/* Footer */
.footer {{
  text-align:center; padding:22px; background: rgba(0,0,0,0.4); border-radius:12px;
  color:#fff; font-size:13px; line-height:1.7; margin-top:20px;
}}
.footer strong {{ color:#fde68a; }}
.footer a {{ color:#c4b5fd; text-decoration:none; }}
.footer a:hover {{ text-decoration:underline; }}

@media (max-width:768px) {{
  .ascii-art {{ font-size:6px; }}
  .title-block h1 {{ font-size:1.5em; }}
  .toolbar {{ flex-direction:column; align-items:stretch; }}
  .export-bar {{ margin-left:0; }}
}}
</style>
</head>
<body>
<div class="container">

  <!-- Header -->
  <div class="ascii-container"><div class="ascii-art">{ASCII_ART}</div></div>
  <div class="title-block">
    <h1>🔬 EcoliTyper Comprehensive Report</h1>
    <p>Unified CHTyper · Serotype · MLST · Phylogroup analysis</p>
  </div>

  <!-- Quote -->
  <div class="quote-container">
    <div class="quote-text" id="quoteText"></div>
    <div class="quote-author" id="quoteAuthor"></div>
  </div>

  {missing_banner}

  <!-- Tabs -->
  <div class="tabs">
    <button class="tab-btn active" data-tab="tab-overview">📊 Overview</button>
    <button class="tab-btn" data-tab="tab-samples">🧬 Sample Data</button>
    <button class="tab-btn" data-tab="tab-dist">📈 Distributions</button>
    <button class="tab-btn" data-tab="tab-quality">✅ Quality Metrics</button>
  </div>

  <!-- Overview -->
  <div class="tab-panel active" id="tab-overview">
    <div class="section">
      <h2>📊 Analysis Overview</h2>
      <div class="stats-grid">{overview_cards}</div>
      <p style="color:#6b7280; font-size:13px; margin-top:10px;">
        <strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ·
        <strong>Report version:</strong> 1.0.0 ·
        <strong>Author:</strong> Brown Beckley
      </p>
    </div>
  </div>

  <!-- Samples -->
  <div class="tab-panel" id="tab-samples">
    <div class="section">
      <h2>🧬 Sample Typing Data ({total} samples)</h2>

      <div class="toolbar">
        <div class="search-wrap">
          <input type="text" id="tableSearch" placeholder="Search across all columns...">
        </div>

        <div class="filter-pills">
          {phylo_pills}
        </div>

        <select class="filter-select" id="seroFilter">
          <option value="">Serotype — all</option>
          {sero_opts}
        </select>
        <select class="filter-select" id="stFilter">
          <option value="">ST — all</option>
          {st_opts}
        </select>
        <select class="filter-select" id="chFilter">
          <option value="">CH type — all</option>
          {ch_opts}
        </select>
        <select class="filter-select" id="statusFilter">
          <option value="">Status — all</option>
          <option value="Completed">Completed</option>
        </select>

        <div class="col-toggle-wrap">
          <button class="action-btn secondary" id="colToggleBtn">⚙ Columns</button>
          <div class="col-panel" id="colPanel"></div>
        </div>

        <span class="sel-counter" id="resultsCounter">{total} of {total} shown</span>
        <span class="sel-counter" id="selCount">0 selected</span>
      </div>

      <div class="toolbar">
        <label style="font-weight:600; color:#4c1d95; font-size:13px;">
          Export scope:
          <select class="filter-select" id="exportScope" style="margin-left:6px;">
            <option value="visible">Visible rows</option>
            <option value="selected">Selected rows</option>
          </select>
        </label>

        <div class="export-bar">
          <button class="action-btn success" id="exportCSV">⬇ CSV</button>
          <button class="action-btn success" id="exportTSV">⬇ TSV</button>
          <button class="action-btn success" id="exportXLSX">⬇ Excel</button>
          <button class="action-btn info"    id="exportJSON">⬇ JSON</button>
          <button class="action-btn secondary" id="copyBtn">📋 Copy</button>
        </div>
      </div>

      <div class="table-container">
        <table id="mainTable">
          <thead>
            <tr>
              <th class="col-check"><input type="checkbox" id="selectAll"></th>
              <th>Sample_ID</th>
              <th>Phylogroup</th>
              <th>Serotype</th>
              <th>O Type</th>
              <th>H Type</th>
              <th>MLST ST</th>
              <th>FumC</th>
              <th>FimH</th>
              <th>CH Type</th>
              <th>FumC Id%</th>
              <th>FimH Id%</th>
              <th>Status</th>
            </tr>
          </thead>
          <tbody>{table_rows}</tbody>
        </table>
      </div>
    </div>
  </div>

  <!-- Distributions -->
  <div class="tab-panel" id="tab-dist">
    <div class="section">
      <h2>📈 Type Distributions &amp; Frequencies</h2>
      <p style="color:#6b7280; font-size:13px; margin-bottom:14px;">
        Each card shows unique values and frequency (% of typed samples).
      </p>
      {dist_cards_html}
    </div>
  </div>

  <!-- Quality -->
  <div class="tab-panel" id="tab-quality">
    <div class="section">
      <h2>✅ Quality Metrics</h2>
      {quality_cards}
      <p style="color:#6b7280; font-size:13px; margin-top:14px;">
        Identity = BLAST sequence identity against reference alleles.
        Coverage = fraction of the reference gene covered by the query.
        High-confidence cutoff: ≥99% identity on both loci.
      </p>
    </div>
  </div>

  <!-- Footer -->
  <div class="footer">
    <p><strong>EcoliTyper</strong> — Comprehensive Analysis Report</p>
    <p><strong>Author:</strong> Brown Beckley · <strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">bbeckley-hub</a></p>
    <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
    <p><strong>Affiliation:</strong> University of Ghana Medical School — Department of Medical Biochemistry</p>
    <p style="margin-top:8px; opacity:0.8;"><em>Turning E. coli sequence data into structured epidemiological insight 🧬</em></p>
  </div>

</div>

<script>
{js}
</script>
</body>
</html>'''

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        print("   ✓ HTML report generated")
        return True
    except Exception as e:
        print(f"   ✗ Error generating HTML: {e}")
        import traceback; traceback.print_exc()
        return False


# =========================================================================
# Main
# =========================================================================
def main():
    parser = argparse.ArgumentParser(
        description='EcoliTyper Comprehensive Report Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
All input files are OPTIONAL — the report is generated from whatever is
present in the target directory. Missing files simply yield 'Not Assigned'
values in the corresponding columns.

Recognized files:
  • chtyper_results.tsv
  • serotype_analysis_report.tsv
  • mlst_summary.tsv
  • phylogrouping_results.tsv
        """
    )
    parser.add_argument('--dir', '-d', default=None,
                        help='Directory containing summary TSVs (default: current directory)')
    parser.add_argument('--no-xlsx', action='store_true',
                        help='Skip Excel (.xlsx) export')
    args = parser.parse_args()

    print("\n" + "=" * 64)
    print("ECOLITYPER COMPREHENSIVE REPORT GENERATOR")
    print("=" * 64)

    target_dir = Path(args.dir).resolve() if args.dir else Path.cwd()
    if not target_dir.exists() or not target_dir.is_dir():
        print(f"❌ Error: directory does not exist: {target_dir}")
        sys.exit(1)
    print(f"📂 Working directory: {target_dir}")

    resolved = {k: (target_dir / fname) for k, fname in INPUT_FILES.items()}
    present  = {k: p for k, p in resolved.items() if p.exists()}
    missing  = {k: p for k, p in resolved.items() if not p.exists()}

    print("\n🔎 Input file check:")
    for key, fname in INPUT_FILES.items():
        mark   = "✓" if key in present else "•"
        status = "found" if key in present else "missing (column will be 'Not Assigned')"
        print(f"   {mark} {fname:<40s} — {status}")

    if not present:
        print("\n❌ No recognized input files found — nothing to report.")
        print("   Expected at least one of:")
        for fname in INPUT_FILES.values():
            print(f"     • {fname}")
        sys.exit(1)

    # Load
    chtyper    = load_chtyper_data(present['chtyper'])    if 'chtyper'    in present else {}
    serotype   = load_serotype_data(present['serotype'])  if 'serotype'   in present else {}
    mlst       = load_mlst_data(present['mlst'])          if 'mlst'       in present else {}
    phylogroup = load_phylogroup_data(present['phylogroup']) if 'phylogroup' in present else {}

    combined = combine_data(chtyper, serotype, mlst, phylogroup)
    if not combined:
        print("\n❌ No samples found in any of the available files.")
        sys.exit(1)

    distributions = compute_distributions(combined)
    quality       = compute_quality_metrics(combined)
    missing_sources = [(k, INPUT_FILES[k]) for k in missing]

    base = "ecolityper_comprehensive_report"
    tsv_ok  = generate_tsv_report(combined, target_dir / f"{base}.tsv")
    json_ok = generate_json_report(combined, distributions, quality, missing_sources,
                                   target_dir / f"{base}.json")
    html_ok = generate_html_report(combined, distributions, quality, missing_sources,
                                   target_dir / f"{base}.html")
    xlsx_ok = True
    if not args.no_xlsx:
        xlsx_ok = generate_xlsx_report(combined, distributions, quality,
                                       target_dir / f"{base}.xlsx")

    print("\n" + "=" * 64)
    if tsv_ok and json_ok and html_ok:
        print("✅ COMPREHENSIVE REPORT GENERATED SUCCESSFULLY!")
        print("=" * 64)
        print(f"📁 Reports saved in: {target_dir}")
        print(f"   • {base}.html   (interactive HTML report)")
        print(f"   • {base}.json   (structured JSON)")
        print(f"   • {base}.tsv    (tab-separated data)")
        if xlsx_ok:
            print(f"   • {base}.xlsx   (Excel workbook)")

        print(f"\n📊 Summary: {len(combined)} samples")
        print(f"🧬 Unique — Phylogroups: {distributions['phylogroup']['unique_count']}, "
              f"Serotypes: {distributions['serotype']['unique_count']}, "
              f"STs: {distributions['mlst_st']['unique_count']}, "
              f"CH types: {distributions['ch_type']['unique_count']}")
        print(f"✅ Quality — avg FumC identity: {quality['avg_fumc_identity']}%, "
              f"avg FimH identity: {quality['avg_fimh_identity']}%")

        if missing_sources:
            print("\n⚠️  Missing sources (columns shown as 'Not Assigned'):")
            for key, fname in missing_sources:
                print(f"     • {fname}")

        print(f"\n💡 {random.choice(SCIENCE_QUOTES)['text']}")
    else:
        print("❌ Some reports failed to generate")
    print("=" * 64)


if __name__ == "__main__":
    main()
