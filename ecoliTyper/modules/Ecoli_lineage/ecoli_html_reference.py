#!/usr/bin/env python3
"""
EcoliDB Comprehensive Reference - Full Database HTML Generator
Captures ALL fields from the comprehensive E. coli database
Includes 2025-2026 updates: new lineages, serotypes, hybrid pathotypes
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-05-15
"""

import os
import json
from datetime import datetime
from ecoli_lineage_database import (
    LINEAGE_DATABASE, SEROTYPE_DATABASE, PHYLOGROUP_DATABASE,
    PATHOTYPE_DATABASE, SPECIALIZED_PROFILES, COMPREHENSIVE_REFERENCES,
    CARBAPENEMASE_PRODUCERS
)

def safe_join(lst, max_items=None):
    """Join list items safely, optionally truncate. If max_items is None, return all."""
    if not lst:
        return "None"
    if max_items is None or len(lst) <= max_items:
        return ", ".join(str(x) for x in lst)
    return ", ".join(str(x) for x in lst[:max_items]) + f" (+{len(lst)-max_items} more)"

def format_resistance(profile):
    """Format resistance profile dictionary into HTML"""
    if not profile:
        return "No data"
    parts = []
    for key, value in profile.items():
        if isinstance(value, list):
            parts.append(f"<strong>{key.replace('_', ' ').title()}:</strong> {safe_join(value, 8)}")
        elif isinstance(value, str):
            parts.append(f"<strong>{key.replace('_', ' ').title()}:</strong> {value}")
    return "<br>".join(parts)

def generate_html():
    """Generate comprehensive HTML reference"""
    
    # Calculate statistics
    stats = {
        'lineages': len(LINEAGE_DATABASE),
        'serotypes': len(SEROTYPE_DATABASE),
        'phylogroups': len(PHYLOGROUP_DATABASE),
        'pathotypes': len(PATHOTYPE_DATABASE),
        'specialized': len(SPECIALIZED_PROFILES),
        'carbapenemase': len(CARBAPENEMASE_PRODUCERS),
        'references_pubmed': sum(len(v) for v in COMPREHENSIVE_REFERENCES.get("PUBMED_REFERENCES", {}).values()),
        'references_doi': sum(len(v) for v in COMPREHENSIVE_REFERENCES.get("DOI_REFERENCES", {}).values()),
        'hybrid_pathotypes': sum(1 for pt in PATHOTYPE_DATABASE.values() if pt.get('category') == 'Hybrid')
    }
    
    # Category breakdowns
    lineage_cats = {}
    for info in LINEAGE_DATABASE.values():
        cat = info.get('category', 'Unknown')
        lineage_cats[cat] = lineage_cats.get(cat, 0) + 1
    
    pathotype_cats = {}
    for info in PATHOTYPE_DATABASE.values():
        cat = info.get('category', 'Unknown')
        pathotype_cats[cat] = pathotype_cats.get(cat, 0) + 1
    
    # Start building HTML
    html_parts = []
    
    html_parts.append(f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliDB v2.0 – Complete E. coli Reference (2025-2026 Update)</title>
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <style>
        :root {{
            --primary: #1a365d;
            --primary-light: #2d3748;
            --secondary: #2b6cb0;
            --accent: #dd6b20;
            --success: #38a169;
            --warning: #d69e2e;
            --danger: #e53e3e;
            --info: #3182ce;
            --new: #9b59b6;
            --hybrid: #e67e22;
            --gray-50: #f9fafb;
            --gray-100: #f3f4f6;
            --gray-200: #e5e7eb;
            --gray-300: #d1d5db;
            --gray-400: #9ca3af;
            --gray-500: #6b7280;
            --gray-600: #4b5563;
            --gray-700: #374151;
            --gray-800: #1f2937;
            --gray-900: #111827;
        }}
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: linear-gradient(135deg, var(--primary) 0%, var(--gray-900) 100%);
            line-height: 1.6;
            min-height: 100vh;
        }}
        .app-container {{ max-width: 1400px; margin: 0 auto; padding: 1rem; }}
        .header {{ text-align: center; margin-bottom: 2rem; color: white; }}
        .logo {{ font-size: 2rem; font-weight: bold; display: flex; align-items: center; justify-content: center; gap: 1rem; }}
        .version-badge {{ background: var(--new); padding: 0.2rem 0.8rem; border-radius: 20px; font-size: 0.8rem; }}
        .update-banner {{ background: linear-gradient(135deg, var(--new), var(--hybrid)); padding: 0.8rem; border-radius: 12px; margin-bottom: 1.5rem; }}
        .stats-overview {{
            display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 1rem;
            margin-bottom: 2rem;
        }}
        .stat-card {{ background: rgba(255,255,255,0.95); padding: 1rem; border-radius: 12px; text-align: center; }}
        .stat-number {{ font-size: 2rem; font-weight: bold; color: var(--primary); }}
        .main-nav {{ background: white; border-radius: 16px; overflow: hidden; margin-bottom: 2rem; }}
        .nav-tabs {{
            display: flex; flex-wrap: wrap; background: var(--gray-50); border-bottom: 1px solid var(--gray-200);
            padding: 0 1rem;
        }}
        .nav-tab {{
            padding: 1rem 1.5rem; background: none; border: none; color: var(--gray-600);
            font-weight: 500; cursor: pointer; border-bottom: 3px solid transparent;
            display: flex; align-items: center; gap: 0.5rem;
        }}
        .nav-tab.active {{ color: var(--primary); border-bottom-color: var(--accent); }}
        .content-section {{ display: none; padding: 2rem; background: white; border-radius: 16px; margin-bottom: 2rem; }}
        .content-section.active {{ display: block; animation: fadeIn 0.3s; }}
        @keyframes fadeIn {{ from {{ opacity: 0; }} to {{ opacity: 1; }} }}
        .section-title {{ color: var(--primary); border-bottom: 2px solid var(--gray-200); padding-bottom: 0.5rem; margin-bottom: 1.5rem; font-size: 1.8rem; }}
        .search-section {{ background: var(--gray-50); padding: 1rem; border-radius: 12px; margin-bottom: 1.5rem; }}
        .search-input {{ width: 100%; padding: 0.8rem; border: 1px solid var(--gray-300); border-radius: 8px; }}
        .cards-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(420px, 1fr)); gap: 1.5rem; }}
        .data-card {{ background: white; border: 1px solid var(--gray-200); border-radius: 12px; overflow: hidden; transition: 0.2s; }}
        .data-card:hover {{ transform: translateY(-4px); box-shadow: 0 10px 20px rgba(0,0,0,0.1); }}
        .card-header {{
            padding: 1rem; background: linear-gradient(135deg, var(--primary), var(--primary-light)); color: white;
        }}
        .card-title {{ font-size: 1.3rem; font-weight: bold; }}
        .card-subtitle {{ opacity: 0.9; font-size: 0.85rem; }}
        .badge {{
            display: inline-block; padding: 0.2rem 0.6rem; border-radius: 20px; font-size: 0.7rem; font-weight: 600;
            margin-right: 0.3rem; margin-top: 0.3rem;
        }}
        .badge-risk-very-high {{ background: #c53030; }}
        .badge-risk-high {{ background: var(--danger); }}
        .badge-risk-moderate {{ background: var(--warning); }}
        .badge-risk-low {{ background: var(--success); }}
        .badge-category {{ background: rgba(255,255,255,0.2); }}
        .badge-new {{ background: var(--new); }}
        .badge-hybrid {{ background: var(--hybrid); }}
        .card-content {{ padding: 1rem; }}
        .info-group {{ margin-bottom: 0.8rem; }}
        .info-label {{ font-weight: 600; color: var(--gray-700); font-size: 0.85rem; margin-bottom: 0.2rem; }}
        .info-value {{ color: var(--gray-600); font-size: 0.85rem; }}
        .gene-tag {{
            display: inline-block; background: var(--gray-100); color: var(--gray-700);
            padding: 0.15rem 0.4rem; border-radius: 4px; font-size: 0.7rem; font-family: monospace;
            margin: 0.1rem;
        }}
        .detailed-section {{ background: var(--gray-50); padding: 0.8rem; border-radius: 8px; margin-top: 0.8rem; }}
        .subsection {{ margin-bottom: 0.8rem; }}
        .subsection-title {{ font-weight: 600; color: var(--primary); margin-bottom: 0.3rem; }}
        .footer {{ text-align: center; margin-top: 2rem; padding: 1.5rem; background: rgba(255,255,255,0.1); border-radius: 16px; color: white; }}
        .custom-link {{ color: var(--accent) !important; text-decoration: none; }}
        .custom-link:hover {{ text-decoration: underline; }}
        .amr-heading {{ color: #e53e3e !important; }}
        @media (max-width: 768px) {{ .cards-grid {{ grid-template-columns: 1fr; }} }}
    </style>
</head>
<body>
<div class="app-container">
    <div class="header">
        <div class="logo"><i class="fas fa-dna"></i> EcoliDB Comprehensive Reference <span class="version-badge">v2.0 2026</span></div>
        <div class="update-banner"><i class="fas fa-star-of-life"></i> NEW: 6 emerging lineages • 4 serotypes • 4 hybrid pathotypes • Updated carbapenemase profiles</div>
        <div class="stats-overview">
            <div class="stat-card"><div class="stat-number">{stats['lineages']}</div><div>Lineages</div></div>
            <div class="stat-card"><div class="stat-number">{stats['pathotypes']}</div><div>Pathotypes</div></div>
            <div class="stat-card"><div class="stat-number">{stats['serotypes']}</div><div>Serotypes</div></div>
            <div class="stat-card"><div class="stat-number">{stats['phylogroups']}</div><div>Phylogroups</div></div>
            <div class="stat-card"><div class="stat-number">{stats['carbapenemase']}</div><div>Carbapenemase Types</div></div>
            <div class="stat-card"><div class="stat-number">{stats['references_pubmed']+stats['references_doi']}</div><div>References</div></div>
        </div>
    </div>
    <div class="main-nav">
        <div class="nav-tabs">
            <button class="nav-tab active" onclick="switchTab('overview')"><i class="fas fa-home"></i> Overview</button>
            <button class="nav-tab" onclick="switchTab('lineages')"><i class="fas fa-dna"></i> Lineages</button>
            <button class="nav-tab" onclick="switchTab('pathotypes')"><i class="fas fa-biohazard"></i> Pathotypes</button>
            <button class="nav-tab" onclick="switchTab('serotypes')"><i class="fas fa-tag"></i> Serotypes</button>
            <button class="nav-tab" onclick="switchTab('phylogroups')"><i class="fas fa-project-diagram"></i> Phylogroups</button>
            <button class="nav-tab" onclick="switchTab('carbapenemase')"><i class="fas fa-shield-virus"></i> Carbapenemase</button>
            <button class="nav-tab" onclick="switchTab('specialized')"><i class="fas fa-star"></i> Specialized</button>
            <button class="nav-tab" onclick="switchTab('references')"><i class="fas fa-book"></i> References</button>
        </div>
    </div>
    
    <!-- ========== OVERVIEW TAB (FULL RICH CONTENT) ========== -->
    <div id="overview" class="content-section active">
        <h2 class="section-title"><i class="fas fa-database"></i> Database Overview</h2>
        
        <!-- Welcome -->
        <div class="data-card" style="margin-bottom:1.5rem;">
            <div class="card-header"><div class="card-title">Welcome to EcoliDB Comprehensive Reference</div></div>
            <div class="card-content">
                <p>This database represents our ongoing effort to compile and organize comprehensive information on <strong>Escherichia coli</strong> lineages, pathotypes, serotypes, phylogroups, specialized profiles, and carbapenemase producers for global research and diagnostic applications.</p>
            </div>
        </div>
        
        <div class="cards-grid" style="margin-bottom:1.5rem;">
            <!-- AMR Section -->
            <div class="data-card" style="background: linear-gradient(135deg, var(--primary-light), var(--secondary)); color: white;">
                <div class="card-header" style="background: transparent;"><div class="card-title amr-heading"><i class="fas fa-hands-helping"></i> Join the Fight Against Antimicrobial Resistance</div></div>
                <div class="card-content" style="color: white;">
                    <p>Antimicrobial resistance (AMR) represents one of the most significant global health threats of our time. We invite researchers, clinicians, and public health professionals to collaborate with us in:</p>
                    <ul style="margin: 0.5rem 0 0 1.5rem;">
                        <li>Expanding and validating our E. coli database</li>
                        <li>Sharing regional epidemiological data</li>
                        <li>Developing standardized typing methodologies</li>
                        <li>Advancing AMR surveillance and intervention strategies</li>
                    </ul>
                    <p style="margin-top:0.5rem;"><strong>Together, we can enhance global AMR monitoring and develop more effective treatment strategies.</strong></p>
                </div>
            </div>
            
            <!-- AI Section -->
            <div class="data-card" style="background: var(--info); color: white;">
                <div class="card-header" style="background: transparent;"><div class="card-title"><i class="fas fa-robot"></i> Next Generation: AI-Powered E. coli Prediction</div></div>
                <div class="card-content" style="color: white;">
                    <p>We are currently developing <strong>machine learning and AI approaches</strong> to integrate results from EcoliTyper and predict complete pattern combinations for rapid E. coli characterization.</p>
                    <p><strong>Follow our GitHub repository for upcoming releases and contribute to this open-source initiative:</strong></p>
                    <div style="background: rgba(255,255,255,0.2); padding: 0.5rem; border-radius: 6px; margin: 0.5rem 0;">
                        <i class="fab fa-github"></i> <strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub/EcoliTyper" class="custom-link" target="_blank">https://github.com/bbeckley-hub/EcoliTyper</a>
                    </div>
                    <p>Stay tuned for predictive models that will revolutionize E. coli typing and resistance profiling!</p>
                </div>
            </div>
        </div>
        
        <!-- Feedback Section -->
        <div class="data-card" style="margin-bottom:1.5rem;">
            <div class="card-header"><div class="card-title"><i class="fas fa-comments"></i> We Value Your Input</div></div>
            <div class="card-content">
                <p><strong>Feature Suggestions & Technical Issues:</strong> We welcome feedback to improve this resource. Please contact us with:</p>
                <ul style="margin: 0.5rem 0 0 1.5rem;">
                    <li>Additional E. coli lineages or pathotypes for inclusion</li>
                    <li>Updated epidemiological data from your region</li>
                    <li>Technical issues or data discrepancies</li>
                    <li>Feature requests for future versions</li>
                    <li>Collaboration opportunities in AMR research</li>
                </ul>
                <p style="margin-top:0.5rem;"><strong>Follow our GitHub for the latest developments in AI-powered E. coli prediction models!</strong></p>
            </div>
        </div>
        
        <!-- Side-by-side categories -->
        <div class="cards-grid" style="margin-bottom:1.5rem;">
            <!-- Database Information Card -->
            <div class="data-card">
                <div class="card-header"><div class="card-title"><i class="fas fa-info-circle"></i> Database Information</div></div>
                <div class="card-content">
                    <div class="info-group"><div class="info-label">Last Updated</div><div class="info-value">{datetime.now().strftime('%Y-%m-%d')}</div></div>
                    <div class="info-group"><div class="info-label">Total Data Points</div><div class="info-value">{stats['lineages']+stats['pathotypes']+stats['serotypes']+stats['phylogroups']+stats['carbapenemase']}</div></div>
                    <div class="info-group"><div class="info-label">Coverage</div><div class="info-value">Global</div></div>
                </div>
            </div>
            
            <!-- Lineage Categories -->
            <div class="data-card">
                <div class="card-header"><div class="card-title"><i class="fas fa-sitemap"></i> Lineage Categories</div></div>
                <div class="card-content">''')
    for cat, cnt in lineage_cats.items():
        html_parts.append(f'<div class="info-group"><div class="info-label">{cat}</div><div class="info-value">{cnt} lineages</div></div>')
    html_parts.append('</div></div>')
    
    # Pathotype Categories
    html_parts.append('''
            <div class="data-card">
                <div class="card-header"><div class="card-title"><i class="fas fa-biohazard"></i> Pathotype Categories</div></div>
                <div class="card-content">''')
    for cat, cnt in pathotype_cats.items():
        html_parts.append(f'<div class="info-group"><div class="info-label">{cat}</div><div class="info-value">{cnt} pathotypes</div></div>')
    html_parts.append('</div></div></div>')
    
    # Scientific Context
    html_parts.append('''
        <div class="data-card">
            <div class="card-header"><div class="card-title"><i class="fas fa-flask"></i> Scientific Context</div></div>
            <div class="card-content">
                <p><strong>This reference database captures current understanding of E. coli molecular epidemiology. However, bacterial evolution and horizontal gene transfer continuously generate new variants. Users should supplement this information with recent publications, local surveillance data, and confirmatory laboratory testing for clinical decision-making.</strong></p>
            </div>
        </div>
    </div>
    
    <!-- ========== LINEAGES TAB ========== -->
    <div id="lineages" class="content-section">
        <h2 class="section-title"><i class="fas fa-dna"></i> Lineage Database</h2>
        <div class="search-section"><input type="text" id="lineageSearch" class="search-input" placeholder="Search by ST, name, serotype..."></div>
        <div class="cards-grid" id="lineagesGrid">''')
    
    for st, info in sorted(LINEAGE_DATABASE.items()):
        first_year = info.get('epidemiology', {}).get('first_identified', 0)
        is_new = (isinstance(first_year, int) and first_year >= 2020)
        risk = info.get('risk_level', 'MODERATE').lower().replace(' ', '-')
        html_parts.append(f'''
            <div class="data-card" data-name="{st} {info.get('primary_name','')} {info.get('serotype','')}">
                <div class="card-header">
                    <div class="card-title">{st} {info.get('primary_name','')}''')
        if is_new:
            html_parts.append(' <span class="badge badge-new">NEW 2025-2026</span>')
        html_parts.append(f'''</div>
                    <div class="card-subtitle">Serotype: {info.get('serotype','?')} | Phylogroup: {info.get('phylogroup','?')} | Risk: {info.get('risk_level','?')}</div>
                    <div><span class="badge badge-category">{info.get('category','')}</span>
                    <span class="badge badge-risk-{risk}">{info.get('risk_level','')}</span></div>
                </div>
                <div class="card-content">
                    <div class="info-group"><div class="info-label">Sublineages</div><div class="info-value">{safe_join(info.get('sublineages', []), 5)}</div></div>
                    <div class="info-group"><div class="info-label">Molecular typing</div><div class="info-value">fumC: {info.get('fumC','?')} | fimH: {info.get('fimH','?')} | Clermont: {info.get('clermont_complex','?')}</div></div>
                    <div class="info-group"><div class="info-label">Pathotypes</div><div class="info-value">{safe_join(info.get('pathotypes', []))}</div></div>
                    <div class="info-group"><div class="info-label">Key virulence genes</div><div class="gene-list">''')
        for gene in info.get('key_virulence_genes', [])[:12]:
            html_parts.append(f'<span class="gene-tag">{gene}</span>')
        html_parts.append('</div></div>')
        if info.get('resistance_profile'):
            html_parts.append(f'''<div class="detailed-section">
                <div class="subsection-title"><i class="fas fa-shield-alt"></i> Resistance Profile</div>
                <div class="info-value">{format_resistance(info['resistance_profile'])}</div></div>''')
        epi = info.get('epidemiology', {})
        if epi:
            html_parts.append(f'''<div class="detailed-section">
                <div class="subsection-title"><i class="fas fa-globe"></i> Epidemiology</div>
                <div class="info-value"><strong>First identified:</strong> {epi.get('first_identified','?')}<br>
                <strong>Global distribution:</strong> {epi.get('global_distribution', epi.get('distribution','?'))}<br>
                <strong>Transmission:</strong> {epi.get('transmission','?')}<br>
                <strong>Reservoir:</strong> {epi.get('reservoir','?')}<br>''')
            geo = epi.get('geographical_distribution', {})
            if geo:
                html_parts.append(f'<strong>High prevalence:</strong> {safe_join(geo.get("high_prevalence",[]))}<br>')
                html_parts.append(f'<strong>Medium prevalence:</strong> {safe_join(geo.get("medium_prevalence",[]))}')
            html_parts.append('</div></div>')
        clin = info.get('clinical_significance', {})
        if clin:
            html_parts.append(f'''<div class="detailed-section">
                <div class="subsection-title"><i class="fas fa-stethoscope"></i> Clinical Significance</div>
                <div class="info-value"><strong>Primary infections:</strong> {safe_join(clin.get('primary_infections',[]))}<br>''')
            if clin.get('mortality'):
                html_parts.append(f'<strong>Mortality:</strong> {clin["mortality"]}<br>')
            if clin.get('treatment_challenges'):
                html_parts.append(f'<strong>Treatment challenges:</strong> {clin["treatment_challenges"]}')
            html_parts.append('</div></div>')
        refs = info.get('key_references', [])
        if refs:
            html_parts.append(f'<div class="info-group"><div class="info-label">References</div><div class="info-value">{safe_join(refs, None)}</div></div>')
        html_parts.append('</div></div>')
    
    html_parts.append('''</div></div>
    
    <!-- ========== PATHOTYPES TAB ========== -->
    <div id="pathotypes" class="content-section">
        <h2 class="section-title"><i class="fas fa-biohazard"></i> Pathotype Database</h2>
        <div class="cards-grid">''')
    
    for pt, info in sorted(PATHOTYPE_DATABASE.items()):
        is_hybrid = info.get('category') == 'Hybrid'
        risk = info.get('risk_level', 'MODERATE').lower().replace(' ', '-')
        html_parts.append(f'''
            <div class="data-card">
                <div class="card-header">
                    <div class="card-title">{pt} {info.get('primary_name','')}''')
        if is_hybrid:
            html_parts.append(' <span class="badge badge-hybrid">HYBRID</span>')
        html_parts.append(f'''</div>
                    <div><span class="badge badge-category">{info.get('category','')}</span>
                    <span class="badge badge-risk-{risk}">{info.get('risk_level','')}</span></div>
                </div>
                <div class="card-content">
                    <div class="info-group"><div class="info-label">Subtypes</div><div class="info-value">{safe_join(info.get('subtypes', []))}</div></div>
                    <div class="info-group"><div class="info-label">Key virulence genes</div><div class="gene-list">''')
        for gene in info.get('key_virulence_genes', [])[:12]:
            html_parts.append(f'<span class="gene-tag">{gene}</span>')
        html_parts.append('</div></div>')
        path = info.get('pathogenesis', {})
        if path:
            html_parts.append(f'<div class="info-group"><div class="info-label">Pathogenesis</div><div class="info-value">{path.get("mechanism","")}</div></div>')
        clin = info.get('clinical_manifestations', {})
        if clin:
            html_parts.append(f'<div class="info-group"><div class="info-label">Clinical</div><div class="info-value">{clin.get("primary","")}<br>')
            if clin.get('complications'):
                html_parts.append(f'<strong>Complications:</strong> {clin["complications"]}')
            html_parts.append('</div></div>')
        if info.get('resistance_profile'):
            html_parts.append(f'<div class="detailed-section"><div class="subsection-title">Resistance</div><div class="info-value">{format_resistance(info["resistance_profile"])}</div></div>')
        if 'serotypes' in info:
            sero = info['serotypes']
            html_parts.append(f'<div class="info-group"><div class="info-label">Common serotypes</div><div class="info-value">{safe_join(sero.get("common", []), 6)}</div></div>')
        html_parts.append(f'<div class="info-group"><div class="info-label">Outbreak potential</div><div class="info-value">{info.get("outbreak_potential","?")}</div></div>')
        refs = info.get('key_references', [])
        if refs:
            html_parts.append(f'<div class="info-group"><div class="info-label">References</div><div class="info-value">{safe_join(refs, None)}</div></div>')
        html_parts.append('</div></div>')
    
    html_parts.append('''</div></div>
    
    <!-- ========== SEROTYPES TAB ========== -->
    <div id="serotypes" class="content-section">
        <h2 class="section-title"><i class="fas fa-tag"></i> Serotype Database</h2>
        <div class="cards-grid">''')
    
    for sero, info in sorted(SEROTYPE_DATABASE.items()):
        risk = info.get('h_us_risk', 'MODERATE').lower().replace(' ', '-')
        html_parts.append(f'''
            <div class="data-card">
                <div class="card-header"><div class="card-title">{sero}</div>
                <div class="card-subtitle">{info.get('primary_pathotype','')}</div>
                <div><span class="badge badge-risk-{risk}">HUS Risk: {info.get('h_us_risk','?')}</span></div></div>
                <div class="card-content">
                    <div class="info-group"><div class="info-label">Sequence Types</div><div class="info-value">{safe_join([f'ST{s}' for s in info.get('st', []) if s], 8)}</div></div>
                    <div class="info-group"><div class="info-label">Key virulence</div><div class="gene-list">''')
        for gene in info.get('key_virulence', [])[:8]:
            html_parts.append(f'<span class="gene-tag">{gene}</span>')
        html_parts.append('</div></div>')
        if 'shiga_toxin_profile' in info:
            stx = info['shiga_toxin_profile']
            html_parts.append(f'<div class="info-group"><div class="info-label">Shiga toxin profile</div><div class="info-value">{stx if isinstance(stx,str) else str(stx)}</div></div>')
        html_parts.append(f'<div class="info-group"><div class="info-label">Outbreak association</div><div class="info-value">{info.get("outbreak_association","?")}</div></div>')
        geo = info.get('geographical_distribution', {})
        if geo:
            html_parts.append(f'<div class="detailed-section"><div class="subsection-title">Geography</div><div class="info-value">')
            for k,v in geo.items():
                if isinstance(v, list):
                    html_parts.append(f'<strong>{k}:</strong> {safe_join(v,5)}<br>')
                elif isinstance(v, dict):
                    html_parts.append(f'<strong>{k}:</strong> {str(v)[:100]}<br>')
                else:
                    html_parts.append(f'<strong>{k}:</strong> {v}<br>')
            html_parts.append('</div></div>')
        refs = info.get('references', [])
        if refs:
            html_parts.append(f'<div class="info-group"><div class="info-label">References</div><div class="info-value">{safe_join(refs, None)}</div></div>')
        html_parts.append('</div></div>')
    
    html_parts.append('''</div></div>
    
    <!-- ========== PHYLOGROUPS TAB ========== -->
    <div id="phylogroups" class="content-section">
        <h2 class="section-title"><i class="fas fa-project-diagram"></i> Phylogroups</h2>
        <div class="cards-grid">''')
    for pg, info in sorted(PHYLOGROUP_DATABASE.items()):
        html_parts.append(f'''
            <div class="data-card">
                <div class="card-header"><div class="card-title">Phylogroup {pg}</div><div class="card-subtitle">{info.get('characteristics','')}</div></div>
                <div class="card-content">
                    <div class="info-group"><div class="info-label">Pathogenic potential</div><div class="info-value">{info.get('pathogenic_potential','')}</div></div>
                    <div class="info-group"><div class="info-label">Common STs</div><div class="info-value">{safe_join([f'ST{s}' for s in info.get('common_st',[])], 10)}</div></div>
                    <div class="info-group"><div class="info-label">Common serotypes</div><div class="info-value">{safe_join(info.get('serotypes',[]), 6)}</div></div>
                    <div class="info-group"><div class="info-label">Virulence genes</div><div class="gene-list">''')
        for gene in info.get('virulence_genes', [])[:8]:
            html_parts.append(f'<span class="gene-tag">{gene}</span>')
        html_parts.append('</div></div>')
        if info.get('notes'):
            html_parts.append(f'<div class="info-group"><div class="info-label">Notes</div><div class="info-value">{info["notes"]}</div></div>')
        html_parts.append('</div></div>')
    
    html_parts.append('''</div></div>
    
    <!-- ========== CARBAPENEMASE PRODUCERS TAB ========== -->
    <div id="carbapenemase" class="content-section">
        <h2 class="section-title"><i class="fas fa-shield-virus"></i> Carbapenemase-Producing E. coli</h2>
        <div class="cards-grid">''')
    for cp, data in CARBAPENEMASE_PRODUCERS.items():
        html_parts.append(f'''
            <div class="data-card">
                <div class="card-header"><div class="card-title">{cp.replace('_',' ').title()}</div></div>
                <div class="card-content">
                    <div class="info-group"><div class="info-label">Associated STs</div><div class="info-value">{safe_join([f'ST{s}' for s in data.get('st',[])], 12)}</div></div>
                    <div class="info-group"><div class="info-label">Carbapenemase genes</div><div class="gene-list">''')
        for gene in data.get('carbapenemase', []):
            html_parts.append(f'<span class="gene-tag">{gene}</span>')
        html_parts.append('</div></div>')
        html_parts.append(f'<div class="info-group"><div class="info-label">Enzyme class</div><div class="info-value">{data.get("enzyme_class","")}</div></div>')
        html_parts.append(f'<div class="info-group"><div class="info-label">Hydrolysis spectrum</div><div class="info-value">{safe_join(data.get("hydrolysis_spectrum",[]))}</div></div>')
        inhib = data.get('inhibitor_profile', {})
        if inhib:
            html_parts.append(f'<div class="info-group"><div class="info-label">Inhibitors</div><div class="info-value"><strong>Inhibited by:</strong> {safe_join(inhib.get("inhibited_by",[]))}<br><strong>Resistant to:</strong> {safe_join(inhib.get("resistant_to",[]))}</div></div>')
        treat = data.get('treatment_options', {})
        if treat:
            html_parts.append(f'<div class="detailed-section"><div class="subsection-title">Treatment</div><div class="info-value"><strong>First line:</strong> {safe_join(treat.get("first_line",[]))}<br><strong>Alternatives:</strong> {safe_join(treat.get("alternative",[]),5)}<br><strong>Notes:</strong> {treat.get("important_notes","")}</div></div>')
        geo_cp = data.get('geographical_distribution', {})
        if geo_cp:
            html_parts.append(f'<div class="info-group"><div class="info-label">Geography</div><div class="info-value"><strong>Endemic regions:</strong> {safe_join(geo_cp.get("endemic_regions",[]))}</div></div>')
        refs = data.get('references', [])
        if refs:
            html_parts.append(f'<div class="info-group"><div class="info-label">References</div><div class="info-value">{safe_join(refs, None)}</div></div>')
        html_parts.append('</div></div>')
    
    html_parts.append('''</div></div>
    
    <!-- ========== SPECIALIZED PROFILES TAB ========== -->
    <div id="specialized" class="content-section">
        <h2 class="section-title"><i class="fas fa-star"></i> Specialized Profiles</h2>''')
    for sp_type, profiles in SPECIALIZED_PROFILES.items():
        html_parts.append(f'<h3 style="margin-top:1rem;">{sp_type.replace("_"," ").title()}</h3><div class="cards-grid">')
        for name, data in profiles.items():
            html_parts.append(f'''
                <div class="data-card"><div class="card-header"><div class="card-title">{name}</div></div>
                <div class="card-content">''')
            for k,v in data.items():
                if isinstance(v, list):
                    html_parts.append(f'<div class="info-group"><div class="info-label">{k.replace("_"," ").title()}</div><div class="gene-list">')
                    for item in v[:6]:
                        html_parts.append(f'<span class="gene-tag">{item}</span>')
                    html_parts.append('</div></div>')
                else:
                    html_parts.append(f'<div class="info-group"><div class="info-label">{k.replace("_"," ").title()}</div><div class="info-value">{v}</div></div>')
            html_parts.append('</div></div>')
        html_parts.append('</div>')
    
    html_parts.append('''</div>
    
    <!-- ========== REFERENCES TAB ========== -->
    <div id="references" class="content-section">
        <h2 class="section-title"><i class="fas fa-book"></i> References</h2>
        <div class="cards-grid">
            <div class="data-card"><div class="card-header"><div class="card-title">PubMed References</div></div>
            <div class="card-content">''')
    for cat, refs in COMPREHENSIVE_REFERENCES.get("PUBMED_REFERENCES", {}).items():
        html_parts.append(f'<div class="info-group"><div class="info-label">{cat.replace("_"," ").title()}</div><div class="info-value">{safe_join(refs, None)}</div></div>')
    html_parts.append('</div></div>')
    html_parts.append(f'''<div class="data-card"><div class="card-header"><div class="card-title">DOI References</div></div>
            <div class="card-content">''')
    for cat, refs in COMPREHENSIVE_REFERENCES.get("DOI_REFERENCES", {}).items():
        html_parts.append(f'<div class="info-group"><div class="info-label">{cat.replace("_"," ").title()}</div><div class="info-value">{safe_join(refs, None)}</div></div>')
    html_parts.append('</div></div></div></div>')
    
    # Footer
    html_parts.append(f'''
    <div class="footer">
        <p><strong>Author:</strong> Brown Beckley | University of Ghana Medical School | <a href="mailto:brownbeckley94@gmail.com" class="custom-link">brownbeckley94@gmail.com</a></p>
        <p>Version 2.0 (2025-2026 update) | Generated {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p><i class="fab fa-github"></i> <a href="https://github.com/bbeckley-hub/EcoliTyper" class="custom-link" target="_blank">GitHub Repository</a></p>
    </div>
</div>
<script>
function switchTab(tabId) {{
    document.querySelectorAll('.content-section').forEach(s => s.classList.remove('active'));
    document.querySelectorAll('.nav-tab').forEach(t => t.classList.remove('active'));
    document.getElementById(tabId).classList.add('active');
    event.target.classList.add('active');
}}
document.getElementById('lineageSearch').addEventListener('input', function() {{
    let term = this.value.toLowerCase();
    document.querySelectorAll('#lineagesGrid .data-card').forEach(card => {{
        let text = card.getAttribute('data-name') || card.innerText.toLowerCase();
        card.style.display = text.includes(term) ? 'block' : 'none';
    }});
}});
</script>
</body>
</html>''')
    
    output = ''.join(html_parts)
    out_file = "ecoli_comprehensive_reference.html"
    with open(out_file, 'w', encoding='utf-8') as f:
        f.write(output)
    print(f"✅ Comprehensive HTML reference generated: {out_file}")
    print(f"📊 Statistics: {stats['lineages']} lineages, {stats['pathotypes']} pathotypes, {stats['serotypes']} serotypes, {stats['carbapenemase']} carbapenemase profiles")
    print(f"🔗 Open {out_file} in your browser")
    return out_file

if __name__ == "__main__":
    generate_html()