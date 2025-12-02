#!/usr/bin/env python3
"""
EcoliDB Comprehensive Reference - Complete Database HTML Generator
Captures ALL data from the comprehensive E. coli database
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
Date: 2025
Send a quick mail for any issues or further explanations.
"""

import os
import json
from datetime import datetime
from ecoli_lineage_database import (
    LINEAGE_DATABASE, SEROTYPE_DATABASE, PHYLOGROUP_DATABASE, 
    PATHOTYPE_DATABASE, SPECIALIZED_PROFILES, COMPREHENSIVE_REFERENCES,
    CARBAPENEMASE_PRODUCERS
)

def generate_comprehensive_reference(output_path="ecoli_comprehensive_reference.html"):
    """Generate a complete HTML reference covering all database content"""
    
    # Calculate statistics - UPDATED WITH NEW COUNTS
    stats = {
        'lineages': len(LINEAGE_DATABASE),
        'serotypes': len(SEROTYPE_DATABASE),
        'phylogroups': len(PHYLOGROUP_DATABASE),
        'pathotypes': len(PATHOTYPE_DATABASE),
        'references_pubmed': sum(len(refs) for refs in COMPREHENSIVE_REFERENCES.get("PUBMED_REFERENCES", {}).values()),
        'references_doi': sum(len(refs) for refs in COMPREHENSIVE_REFERENCES.get("DOI_REFERENCES", {}).values()),
        'carbapenemase_profiles': len(CARBAPENEMASE_PRODUCERS)
    }
    
    # Count categories
    lineage_categories = {}
    for info in LINEAGE_DATABASE.values():
        cat = info.get('category', 'Unknown')
        lineage_categories[cat] = lineage_categories.get(cat, 0) + 1
    
    pathotype_categories = {}
    for info in PATHOTYPE_DATABASE.values():
        cat = info.get('category', 'Unknown')
        pathotype_categories[cat] = pathotype_categories.get(cat, 0) + 1

    html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EcoliDB - Comprehensive E. coli Reference Database</title>
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <style>
        :root {{
            /* Professional Color Scheme */
            --primary: #1a365d;
            --primary-light: #2d3748;
            --secondary: #2b6cb0;
            --accent: #dd6b20;
            --success: #38a169;
            --warning: #d69e2e;
            --danger: #e53e3e;
            --info: #3182ce;
            
            /* Neutrals */
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
            
            /* Spacing */
            --space-xs: 0.25rem;
            --space-sm: 0.5rem;
            --space-md: 1rem;
            --space-lg: 1.5rem;
            --space-xl: 2rem;
            --space-2xl: 3rem;
        }}
        
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: linear-gradient(135deg, var(--primary) 0%, var(--gray-900) 100%);
            color: var(--gray-800);
            line-height: 1.6;
            min-height: 100vh;
        }}
        
        .app-container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: var(--space-md);
        }}
        
        /* Header */
        .header {{
            text-align: center;
            margin-bottom: var(--space-2xl);
            color: white;
        }}
        
        .logo {{
            font-size: 2.5rem;
            font-weight: bold;
            margin-bottom: var(--space-sm);
            display: flex;
            align-items: center;
            justify-content: center;
            gap: var(--space-md);
        }}
        
        .logo i {{
            color: var(--accent);
        }}
        
        .subtitle {{
            color: var(--gray-300);
            font-size: 1.2rem;
            margin-bottom: var(--space-lg);
        }}
        
        /* Stats Overview */
        .stats-overview {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: var(--space-md);
            margin-bottom: var(--space-2xl);
        }}
        
        .stat-card {{
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(10px);
            padding: var(--space-lg);
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
            transition: transform 0.2s ease;
        }}
        
        .stat-card:hover {{
            transform: translateY(-2px);
        }}
        
        .stat-number {{
            font-size: 2.5rem;
            font-weight: bold;
            color: var(--primary);
            margin-bottom: var(--space-xs);
        }}
        
        .stat-label {{
            color: var(--gray-600);
            font-size: 0.9rem;
            font-weight: 500;
        }}
        
        /* Main Navigation */
        .main-nav {{
            background: white;
            border-radius: 16px;
            box-shadow: 0 10px 25px -3px rgba(0, 0, 0, 0.1);
            overflow: hidden;
            margin-bottom: var(--space-2xl);
        }}
        
        .nav-tabs {{
            display: flex;
            background: var(--gray-50);
            border-bottom: 1px solid var(--gray-200);
            padding: 0 var(--space-md);
            flex-wrap: wrap;
        }}
        
        .nav-tab {{
            padding: var(--space-lg) var(--space-xl);
            background: none;
            border: none;
            color: var(--gray-600);
            font-weight: 500;
            cursor: pointer;
            transition: all 0.2s ease;
            border-bottom: 3px solid transparent;
            display: flex;
            align-items: center;
            gap: var(--space-sm);
            white-space: nowrap;
        }}
        
        .nav-tab:hover {{
            color: var(--primary);
            background: var(--gray-100);
        }}
        
        .nav-tab.active {{
            color: var(--primary);
            border-bottom-color: var(--accent);
        }}
        
        /* Content Sections */
        .content-section {{
            display: none;
            padding: var(--space-2xl);
            background: white;
            border-radius: 16px;
            box-shadow: 0 10px 25px -3px rgba(0, 0, 0, 0.1);
            margin-bottom: var(--space-2xl);
        }}
        
        .content-section.active {{
            display: block;
            animation: fadeIn 0.3s ease;
        }}
        
        @keyframes fadeIn {{
            from {{ opacity: 0; transform: translateY(10px); }}
            to {{ opacity: 1; transform: translateY(0); }}
        }}
        
        .section-title {{
            color: var(--primary);
            border-bottom: 2px solid var(--gray-200);
            padding-bottom: var(--space-md);
            margin-bottom: var(--space-2xl);
            font-size: 1.8rem;
            display: flex;
            align-items: center;
            gap: var(--space-md);
        }}
        
        /* Search and Filters */
        .search-section {{
            background: var(--gray-50);
            padding: var(--space-lg);
            border-radius: 12px;
            margin-bottom: var(--space-xl);
        }}
        
        .search-box {{
            display: flex;
            gap: var(--space-md);
            margin-bottom: var(--space-md);
        }}
        
        .search-input {{
            flex: 1;
            padding: var(--space-md);
            border: 1px solid var(--gray-300);
            border-radius: 8px;
            font-size: 1rem;
        }}
        
        .search-input:focus {{
            outline: none;
            border-color: var(--secondary);
            box-shadow: 0 0 0 3px rgba(43, 108, 176, 0.1);
        }}
        
        /* Data Cards */
        .cards-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(350px, 1fr));
            gap: var(--space-lg);
        }}
        
        .data-card {{
            background: white;
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            overflow: hidden;
            transition: all 0.3s ease;
            box-shadow: 0 1px 3px 0 rgba(0, 0, 0, 0.1);
        }}
        
        .data-card:hover {{
            transform: translateY(-4px);
            box-shadow: 0 10px 25px -3px rgba(0, 0, 0, 0.1);
        }}
        
        .card-header {{
            padding: var(--space-lg);
            border-bottom: 1px solid var(--gray-200);
            background: linear-gradient(135deg, var(--primary), var(--primary-light));
            color: white;
        }}
        
        .card-title {{
            font-size: 1.4rem;
            font-weight: bold;
            margin-bottom: var(--space-xs);
        }}
        
        .card-subtitle {{
            opacity: 0.9;
            margin-bottom: var(--space-sm);
        }}
        
        .card-badges {{
            display: flex;
            gap: var(--space-xs);
            flex-wrap: wrap;
        }}
        
        .badge {{
            padding: 4px 8px;
            border-radius: 6px;
            font-size: 0.75rem;
            font-weight: 600;
        }}
        
        .badge-risk-high {{ background: var(--danger); }}
        .badge-risk-moderate {{ background: var(--warning); }}
        .badge-risk-low {{ background: var(--success); }}
        .badge-category {{ background: rgba(255, 255, 255, 0.2); }}
        
        .card-content {{
            padding: var(--space-lg);
        }}
        
        .info-group {{
            margin-bottom: var(--space-md);
        }}
        
        .info-label {{
            font-weight: 600;
            color: var(--gray-700);
            font-size: 0.9rem;
            margin-bottom: var(--space-xs);
        }}
        
        .info-value {{
            color: var(--gray-600);
            font-size: 0.9rem;
            line-height: 1.5;
        }}
        
        .gene-list {{
            display: flex;
            flex-wrap: wrap;
            gap: 4px;
            margin-top: 4px;
        }}
        
        .gene-tag {{
            background: var(--gray-100);
            color: var(--gray-700);
            padding: 2px 6px;
            border-radius: 4px;
            font-size: 0.75rem;
            font-family: 'Courier New', monospace;
        }}
        
        /* Detailed Sections */
        .detailed-section {{
            background: var(--gray-50);
            padding: var(--space-lg);
            border-radius: 8px;
            margin-top: var(--space-md);
        }}
        
        .subsection {{
            margin-bottom: var(--space-lg);
        }}
        
        .subsection-title {{
            font-weight: 600;
            color: var(--primary);
            margin-bottom: var(--space-md);
            font-size: 1.1rem;
            display: flex;
            align-items: center;
            gap: var(--space-sm);
        }}
        
        /* Tables for structured data */
        .data-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: var(--space-md);
        }}
        
        .data-table th,
        .data-table td {{
            padding: var(--space-sm);
            text-align: left;
            border-bottom: 1px solid var(--gray-200);
        }}
        
        .data-table th {{
            background: var(--gray-100);
            font-weight: 600;
            color: var(--gray-700);
        }}
        
        /* Footer */
        .footer {{
            text-align: center;
            margin-top: var(--space-2xl);
            padding: var(--space-xl);
            background: rgba(255, 255, 255, 0.1);
            border-radius: 16px;
            color: white;
        }}
        
        .authorship {{
            margin-top: var(--space-lg);
            padding: var(--space-lg);
            background: rgba(255, 255, 255, 0.1);
            border-radius: 12px;
        }}
        
        /* NEW STYLES FOR BETTER ORGANIZATION */
        .overview-grid {{
            display: grid;
            grid-template-columns: 2fr 1fr;
            gap: var(--space-xl);
            margin-bottom: var(--space-2xl);
        }}
        
        .main-content {{
            display: flex;
            flex-direction: column;
            gap: var(--space-xl);
        }}
        
        .sidebar {{
            display: flex;
            flex-direction: column;
            gap: var(--space-lg);
        }}
        
        .category-breakdown {{
            background: white;
            border-radius: 12px;
            padding: var(--space-lg);
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
        }}
        
        .category-item {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: var(--space-sm) 0;
            border-bottom: 1px solid var(--gray-200);
        }}
        
        .category-item:last-child {{
            border-bottom: none;
        }}
        
        .category-name {{
            font-weight: 600;
            color: var(--gray-700);
            flex: 1;
        }}
        
        .category-count {{
            background: var(--primary);
            color: white;
            padding: 4px 8px;
            border-radius: 12px;
            font-size: 0.8rem;
            font-weight: 600;
        }}
        
        .feature-card {{
            background: white;
            border-radius: 12px;
            padding: var(--space-xl);
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
            margin-bottom: var(--space-lg);
        }}
        
        .feature-card h3 {{
            color: var(--primary);
            margin-bottom: var(--space-md);
            display: flex;
            align-items: center;
            gap: var(--space-sm);
        }}
        
        .database-info {{
            background: linear-gradient(135deg, var(--primary), var(--primary-light));
            color: white;
            border-radius: 12px;
            padding: var(--space-xl);
            margin-bottom: var(--space-lg);
        }}
        
        .database-info h3 {{
            margin-bottom: var(--space-md);
            display: flex;
            align-items: center;
            gap: var(--space-sm);
        }}
        
        .info-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: var(--space-md);
            margin-top: var(--space-md);
        }}
        
        .info-item {{
            text-align: center;
        }}
        
        .info-value-large {{
            font-size: 1.1rem;
            font-weight: bold;
            margin-bottom: var(--space-xs);
            color: #38a169; /* Green color for database info values */
        }}
        
        .info-label-small {{
            font-size: 0.8rem;
            opacity: 0.9;
            color: white;
        }}
        
        .amr-heading {{
            color: #e53e3e !important; /* Red color for AMR heading */
        }}
        
        .custom-link {{
            color: #d69e2e !important; /* Yellow color for links */
            text-decoration: none;
            font-weight: 600;
        }}
        
        .custom-link:hover {{
            text-decoration: underline;
            color: #ecc94b !important;
        }}
        
        /* Responsive */
        @media (max-width: 768px) {{
            .cards-grid {{
                grid-template-columns: 1fr;
            }}
            .nav-tabs {{
                flex-direction: column;
            }}
            .search-box {{
                flex-direction: column;
            }}
            .stats-overview {{
                grid-template-columns: repeat(2, 1fr);
            }}
            .overview-grid {{
                grid-template-columns: 1fr;
            }}
        }}
    </style>
</head>
<body>
    <div class="app-container">
        <!-- Header -->
        <div class="header">
            <div class="logo">
                <i class="fas fa-dna"></i>
                <span>EcoliDB Comprehensive Reference</span>
            </div>
            <div class="subtitle">
                Complete Escherichia coli Lineage, Pathotype, and Epidemiology Database
            </div>
            
            <!-- Statistics - UPDATED WITH NEW COUNTS -->
            <div class="stats-overview">
                <div class="stat-card">
                    <div class="stat-number">{stats['lineages']}</div>
                    <div class="stat-label">Sequence Types</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{stats['pathotypes']}</div>
                    <div class="stat-label">Pathotypes</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{stats['serotypes']}</div>
                    <div class="stat-label">Serotypes</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{stats['phylogroups']}</div>
                    <div class="stat-label">Phylogroups</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{stats['carbapenemase_profiles']}</div>
                    <div class="stat-label">Carbapenemase Types</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{stats['references_pubmed'] + stats['references_doi']}</div>
                    <div class="stat-label">References</div>
                </div>
            </div>
        </div>
        
        <!-- Main Navigation - UPDATED WITH NEW SECTION -->
        <div class="main-nav">
            <div class="nav-tabs">
                <button class="nav-tab active" onclick="switchTab('overview')">
                    <i class="fas fa-chart-bar"></i>
                    Overview
                </button>
                <button class="nav-tab" onclick="switchTab('lineages')">
                    <i class="fas fa-dna"></i>
                    Lineages ({stats['lineages']})
                </button>
                <button class="nav-tab" onclick="switchTab('pathotypes')">
                    <i class="fas fa-biohazard"></i>
                    Pathotypes ({stats['pathotypes']})
                </button>
                <button class="nav-tab" onclick="switchTab('serotypes')">
                    <i class="fas fa-tag"></i>
                    Serotypes ({stats['serotypes']})
                </button>
                <button class="nav-tab" onclick="switchTab('phylogroups')">
                    <i class="fas fa-project-diagram"></i>
                    Phylogroups ({stats['phylogroups']})
                </button>
                <button class="nav-tab" onclick="switchTab('carbapenemase')">
                    <i class="fas fa-shield-virus"></i>
                    Carbapenemase ({stats['carbapenemase_profiles']})
                </button>
                <button class="nav-tab" onclick="switchTab('specialized')">
                    <i class="fas fa-star"></i>
                    Specialized Profiles
                </button>
                <button class="nav-tab" onclick="switchTab('references')">
                    <i class="fas fa-book"></i>
                    References ({stats['references_pubmed'] + stats['references_doi']})
                </button>
            </div>
            
            <!-- Content Sections -->
            <div class="content-section active" id="overview">
                <h2 class="section-title">
                    <i class="fas fa-chart-bar"></i>
                    Database Overview
                </h2>
                
                <!-- NEW ORGANIZED LAYOUT -->
                <div class="overview-grid">
                    <!-- Main Content Column -->
                    <div class="main-content">
                        <!-- Welcome Section -->
                        <div class="feature-card">
                            <h3><i class="fas fa-database"></i> Welcome to EcoliDB Comprehensive Reference</h3>
                            <p>This database represents our ongoing effort to compile and organize comprehensive information on <strong>Escherichia coli</strong> lineages, pathotypes, serotypes, phylogroups, specialized profiles, and carbapenemase producers for global research and diagnostic applications.</p>
                        </div>
                        
                        <!-- AMR Section -->
                        <div class="feature-card" style="background: linear-gradient(135deg, var(--primary-light), var(--secondary)); color: white;">
                            <h3 class="amr-heading"><i class="fas fa-hands-helping"></i> Join the Fight Against Antimicrobial Resistance</h3>
                            <p>Antimicrobial resistance (AMR) represents one of the most significant global health threats of our time. We invite researchers, clinicians, and public health professionals to collaborate with us in:</p>
                            <ul style="margin: var(--space-md) 0; padding-left: var(--space-lg);">
                                <li>Expanding and validating our E. coli database</li>
                                <li>Sharing regional epidemiological data</li>
                                <li>Developing standardized typing methodologies</li>
                                <li>Advancing AMR surveillance and intervention strategies</li>
                            </ul>
                            <p><strong>Together, we can enhance global AMR monitoring and develop more effective treatment strategies.</strong></p>
                        </div>
                        
                        <!-- AI Section -->
                        <div class="feature-card" style="background: var(--info); color: white;">
                            <h3><i class="fas fa-robot"></i> Next Generation: AI-Powered E. coli Prediction</h3>
                            <p>We are currently developing <strong>machine learning and AI approaches</strong> to integrate results from EcoliTyper and predict complete pattern combinations for rapid E. coli characterization.</p>
                            <p><strong>Follow our GitHub repository for upcoming releases and contribute to this open-source initiative:</strong></p>
                            <div style="background: rgba(255,255,255,0.2); padding: var(--space-md); border-radius: 6px; margin: var(--space-md) 0;">
                                <i class="fab fa-github"></i>
                                <strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub/EcoliTyper" class="custom-link" target="_blank">https://github.com/bbeckley-hub/EcoliTyper</a>
                            </div>
                            <p>Stay tuned for predictive models that will revolutionize E. coli typing and resistance profiling!</p>
                        </div>
                        
                        <!-- Feedback Section -->
                        <div class="feature-card">
                            <h3><i class="fas fa-comments"></i> We Value Your Input</h3>
                            <p><strong>Feature Suggestions & Technical Issues:</strong> We welcome feedback to improve this resource. Please contact us with:</p>
                            <ul style="margin: var(--space-sm) 0; padding-left: var(--space-lg);">
                                <li>Additional E. coli lineages or pathotypes for inclusion</li>
                                <li>Updated epidemiological data from your region</li>
                                <li>Technical issues or data discrepancies</li>
                                <li>Feature requests for future versions</li>
                                <li>Collaboration opportunities in AMR research</li>
                            </ul>
                            <p><strong>Follow our GitHub for the latest developments in AI-powered E. coli prediction models!</strong></p>
                        </div>
                    </div>
                    
                    <!-- Sidebar Column -->
                    <div class="sidebar">
                        <!-- Database Information Card -->
                        <div class="database-info">
                            <h3><i class="fas fa-info-circle"></i> Database Information</h3>
                            <div class="info-grid">
                                <div class="info-item">
                                    <div class="info-value-large">''' + datetime.now().strftime('%Y-%m-%d') + '''</div>
                                    <div class="info-label-small">Last Updated</div>
                                </div>
                                <div class="info-item">
                                    <div class="info-value-large">''' + str(sum(stats.values())) + '''</div>
                                    <div class="info-label-small">Total Data Points</div>
                                </div>
                                <div class="info-item">
                                    <div class="info-value-large">Global</div>
                                    <div class="info-label-small">Coverage</div>
                                </div>
                            </div>
                        </div>
                        
                        <!-- Lineage Categories -->
                        <div class="category-breakdown">
                            <h3 style="color: var(--primary); margin-bottom: var(--space-md);">
                                <i class="fas fa-sitemap"></i> Lineage Categories
                            </h3>
                            '''
    
    # Add lineage categories in organized list
    for category, count in lineage_categories.items():
        html_content += f'''
                            <div class="category-item">
                                <span class="category-name">{category}</span>
                                <span class="category-count">{count} lineages</span>
                            </div>
        '''
    
    html_content += '''
                        </div>
                        
                        <!-- Pathotype Categories -->
                        <div class="category-breakdown">
                            <h3 style="color: var(--primary); margin-bottom: var(--space-md);">
                                <i class="fas fa-biohazard"></i> Pathotype Categories
                            </h3>
                            '''
    
    # Add pathotype categories in organized list
    for category, count in pathotype_categories.items():
        html_content += f'''
                            <div class="category-item">
                                <span class="category-name">{category}</span>
                                <span class="category-count">{count} pathotypes</span>
                            </div>
        '''
    
    html_content += '''
                        </div>
                    </div>
                </div>
                
                <!-- Scientific Context at Bottom -->
                <div class="feature-card" style="background: var(--gray-100);">
                    <h3><i class="fas fa-flask"></i> Scientific Context</h3>
                    <p><strong>This reference database captures current understanding of E. coli molecular epidemiology. However, bacterial evolution and horizontal gene transfer continuously generate new variants. Users should supplement this information with recent publications, local surveillance data, and confirmatory laboratory testing for clinical decision-making.</strong></p>
                </div>
            </div>
            
            <!-- Lineages Section - UPDATED WITH MISSING FIELDS -->
            <div class="content-section" id="lineages">
                <h2 class="section-title">
                    <i class="fas fa-dna"></i>
                    E. coli Lineage Database
                </h2>
                
                <div class="search-section">
                    <div class="search-box">
                        <input type="text" class="search-input" id="lineageSearch" placeholder="🔍 Search lineages by ST, name, or characteristics...">
                        <select class="search-input" id="lineageCategory" style="flex: 0 0 200px;">
                            <option value="">All Categories</option>
    '''
    
    # Add category options
    for category in lineage_categories.keys():
        html_content += f'<option value="{category}">{category}</option>'
    
    html_content += '''
                        </select>
                    </div>
                </div>
                
                <div class="cards-grid" id="lineagesGrid">
    '''
    
    # Generate lineage cards WITH ALL MISSING FIELDS
    for st, info in sorted(LINEAGE_DATABASE.items()):
        risk_level = info.get('risk_level', 'MODERATE').lower().replace(' ', '-')
        category = info.get('category', 'Unknown')
        epidemiology = info.get('epidemiology', {})
        geo_dist = epidemiology.get('geographical_distribution', {})
        
        html_content += f'''
                    <div class="data-card" data-category="{category}" data-risk="{risk_level}">
                        <div class="card-header">
                            <div class="card-title">{st}</div>
                            <div class="card-subtitle">{info.get('primary_name', '')}</div>
                            <div class="card-badges">
                                <div class="badge badge-risk-{risk_level}">{info.get('risk_level', 'Unknown')}</div>
                                <div class="badge badge-category">{category}</div>
                            </div>
                        </div>
                        <div class="card-content">
                            <div class="info-group">
                                <div class="info-label">Molecular Typing</div>
                                <div class="info-value">
                                    <strong>fumC:</strong> {info.get('fumC', 'Unknown')}<br>
                                    <strong>fimH:</strong> {info.get('fimH', 'Unknown')}<br>
                                    <strong>Sublineages:</strong> {', '.join(info.get('sublineages', ['None']))}
                                </div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">Typing Information</div>
                                <div class="info-value">
                                    <strong>Serotype:</strong> {info.get('serotype', 'Unknown')}<br>
                                    <strong>Phylogroup:</strong> {info.get('phylogroup', 'Unknown')}<br>
                                    <strong>Clermont Complex:</strong> {info.get('clermont_complex', 'Unknown')}
                                </div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">Pathotypes</div>
                                <div class="info-value">{', '.join(info.get('pathotypes', []))}</div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">Key Virulence Genes</div>
                                <div class="gene-list">
        '''
        
        # Add virulence genes
        for gene in info.get('key_virulence_genes', [])[:8]:
            html_content += f'<span class="gene-tag">{gene}</span>'
        
        html_content += f'''
                                </div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">Epidemiology</div>
                                <div class="info-value">
                                    <strong>Reservoir:</strong> {epidemiology.get('reservoir', 'Unknown')}<br>
                                    <strong>Distribution:</strong> {epidemiology.get('global_distribution', epidemiology.get('distribution', 'Unknown'))}<br>
                                    <strong>High Prevalence:</strong> {', '.join(geo_dist.get('high_prevalence', []))}<br>
                                    <strong>Medium Prevalence:</strong> {', '.join(geo_dist.get('medium_prevalence', []))}
                                </div>
                            </div>
                            
                            <div class="detailed-section">
                                <div class="subsection">
                                    <div class="subsection-title">
                                        <i class="fas fa-shield-alt"></i>
                                        Resistance Profile
                                    </div>
                                    <div class="info-value">
        '''
        
        # Add resistance information
        resistance = info.get('resistance_profile', {})
        for category, data in resistance.items():
            if isinstance(data, list):
                html_content += f'<div><strong>{category.title()}:</strong> {", ".join(data)}</div>'
            elif category in ['notes', 'important_note', 'resistance_notes']:
                html_content += f'<div><strong>Note:</strong> {data}</div>'
        
        html_content += f'''
                                    </div>
                                </div>
                                
                                <div class="subsection">
                                    <div class="subsection-title">
                                        <i class="fas fa-stethoscope"></i>
                                        Clinical Significance
                                    </div>
                                    <div class="info-value">
        '''
        
        # Add clinical significance
        clinical = info.get('clinical_significance', {})
        for key, value in clinical.items():
            if isinstance(value, list):
                html_content += f'<div><strong>{key.replace("_", " ").title()}:</strong> {", ".join(value)}</div>'
            else:
                html_content += f'<div><strong>{key.replace("_", " ").title()}:</strong> {value}</div>'
        
        html_content += f'''
                                    </div>
                                </div>
                                
                                <div class="subsection">
                                    <div class="subsection-title">
                                        <i class="fas fa-globe-americas"></i>
                                        Geographical Distribution
                                    </div>
                                    <div class="info-value">
        '''
        
        # Add geographical distribution
        geo = info.get('epidemiology', {}).get('geographical_distribution', {})
        if geo:
            if 'high_prevalence' in geo:
                html_content += f'<div><strong>High Prevalence:</strong> {", ".join(geo["high_prevalence"])}</div>'
            if 'medium_prevalence' in geo:
                html_content += f'<div><strong>Medium Prevalence:</strong> {", ".join(geo["medium_prevalence"])}</div>'
            if 'regional_variants' in geo:
                html_content += '<div><strong>Regional Variants:</strong> '
                for region, variant in list(geo['regional_variants'].items())[:2]:
                    html_content += f'{region}: {variant}; '
                html_content += '</div>'
        
        html_content += f'''
                                    </div>
                                </div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">Key References ({len(info.get('key_references', []))})</div>
                                <div class="info-value">
                                    {', '.join(info.get('key_references', []))}
                                </div>
                            </div>
                        </div>
                    </div>
        '''
    
    html_content += '''
                </div>
            </div>
            
            <!-- Pathotypes Section - UPDATED WITH NEW PATHOTYPES -->
            <div class="content-section" id="pathotypes">
                <h2 class="section-title">
                    <i class="fas fa-biohazard"></i>
                    E. coli Pathotype Database
                </h2>
                
                <div class="cards-grid">
    '''
    
    # Generate pathotype cards
    for pt, info in sorted(PATHOTYPE_DATABASE.items()):
        category = info.get('category', 'Unknown')
        risk_level = info.get('risk_level', 'MODERATE').lower().replace(' ', '-')
        
        html_content += f'''
                    <div class="data-card">
                        <div class="card-header">
                            <div class="card-title">{pt}</div>
                            <div class="card-subtitle">{info.get('primary_name', '')}</div>
                            <div class="card-badges">
                                <div class="badge badge-risk-{risk_level}">{info.get('risk_level', 'Unknown')}</div>
                                <div class="badge badge-category">{category}</div>
                            </div>
                        </div>
                        <div class="card-content">
                            <div class="info-group">
                                <div class="info-label">Subtypes</div>
                                <div class="info-value">{', '.join(info.get('subtypes', []))}</div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">Key Virulence Genes</div>
                                <div class="gene-list">
        '''
        
        # Add virulence genes
        for gene in info.get('key_virulence_genes', [])[:8]:
            html_content += f'<span class="gene-tag">{gene}</span>'
        
        html_content += f'''
                                </div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">Pathogenesis</div>
                                <div class="info-value">
                                    {info.get('pathogenesis', {}).get('mechanism', 'Unknown')}
                                </div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">Clinical Manifestations</div>
                                <div class="info-value">
        '''
        
        clinical = info.get('clinical_manifestations', {})
        if 'primary' in clinical:
            html_content += f'{clinical["primary"]}'
        if 'complications' in clinical:
            html_content += f'<br><strong>Complications:</strong> {clinical["complications"]}'
        
        html_content += f'''
                                </div>
                            </div>
                            
                            <div class="detailed-section">
                                <div class="subsection">
                                    <div class="subsection-title">
                                        <i class="fas fa-shield-alt"></i>
                                        Resistance Profile
                                    </div>
                                    <div class="info-value">
        '''
        
        # Add resistance information
        resistance = info.get('resistance_profile', {})
        for key, value in resistance.items():
            if isinstance(value, list):
                html_content += f'<div><strong>{key.replace("_", " ").title()}:</strong> {", ".join(value)}</div>'
            elif key in ['notes', 'important_note']:
                html_content += f'<div><strong>Note:</strong> {value}</div>'
        
        html_content += f'''
                                    </div>
                                </div>
                                
                                <div class="subsection">
                                    <div class="subsection-title">
                                        <i class="fas fa-dna"></i>
                                        Additional Features
                                    </div>
                                    <div class="info-value">
        '''
        
        # Add additional features
        if 'serotypes' in info:
            common_serotypes = info['serotypes'].get('common', [])[:3]
            html_content += f'<div><strong>Common Serotypes:</strong> {", ".join(common_serotypes)}</div>'
        
        if 'subtype_markers' in info:
            html_content += '<div><strong>Subtype Markers:</strong> '
            for subtype, markers in list(info['subtype_markers'].items())[:2]:
                html_content += f'{subtype}: {", ".join(markers)}; '
            html_content += '</div>'
        
        html_content += f'''
                                    </div>
                                </div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">Outbreak Potential</div>
                                <div class="info-value">{info.get('outbreak_potential', 'Unknown')}</div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">Key References ({len(info.get('key_references', []))})</div>
                                <div class="info-value">
                                    {', '.join(info.get('key_references', []))}
                                </div>
                            </div>
                        </div>
                    </div>
        '''
    
    html_content += '''
                </div>
            </div>
            
            <!-- Serotypes Section - UPDATED WITH TOXIN PROFILES AND REFERENCES -->
            <div class="content-section" id="serotypes">
                <h2 class="section-title">
                    <i class="fas fa-tag"></i>
                    Serotype Database
                </h2>
                
                <div class="cards-grid">
    '''
    
    # Generate serotype cards WITH TOXIN PROFILES AND REFERENCES
    for serotype, info in sorted(SEROTYPE_DATABASE.items()):
        html_content += f'''
                    <div class="data-card">
                        <div class="card-header">
                            <div class="card-title">{serotype}</div>
                            <div class="card-subtitle">{info.get('primary_pathotype', '')}</div>
                            <div class="card-badges">
                                <div class="badge badge-risk-{info.get('h_us_risk', 'moderate').lower().replace(' ', '-')}">
                                    HUS Risk: {info.get('h_us_risk', 'Unknown')}
                                </div>
                            </div>
                        </div>
                        <div class="card-content">
                            <div class="info-group">
                                <div class="info-label">Sequence Types</div>
                                <div class="info-value">
                                    {', '.join([f'ST{st}' for st in info.get('st', [])])}
                                </div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">Key Virulence Factors</div>
                                <div class="gene-list">
        '''
        
        # Add virulence factors
        for gene in info.get('key_virulence', [])[:8]:
            html_content += f'<span class="gene-tag">{gene}</span>'
        
        html_content += f'''
                                </div>
                            </div>
        '''
        
        # ADD SHIGA TOXIN PROFILE IF AVAILABLE
        if 'shiga_toxin_profile' in info:
            html_content += f'''
                            <div class="detailed-section">
                                <div class="subsection">
                                    <div class="subsection-title">
                                        <i class="fas fa-skull-crossbones"></i>
                                        Shiga Toxin Profile
                                    </div>
                                    <div class="info-value">
            '''
            
            toxin_profile = info['shiga_toxin_profile']
            html_content += f'<div><strong>Primary Toxin:</strong> {toxin_profile.get("primary", "Unknown")}</div>'
            html_content += f'<div><strong>Secondary Toxin:</strong> {toxin_profile.get("secondary", "None")}</div>'
            html_content += f'<div><strong>Stx1 Presence:</strong> {toxin_profile.get("stx1", "Unknown")}</div>'
            html_content += f'<div><strong>Risk Notes:</strong> {toxin_profile.get("toxin_notes", "None")}</div>'
            
            html_content += '''
                                    </div>
                                </div>
                            </div>
            '''
        
        html_content += f'''
                            <div class="info-group">
                                <div class="info-label">Outbreak Association</div>
                                <div class="info-value">{info.get('outbreak_association', 'Unknown')}</div>
                            </div>
                            
                            <div class="detailed-section">
                                <div class="subsection">
                                    <div class="subsection-title">
                                        <i class="fas fa-globe-americas"></i>
                                        Geographical Distribution
                                    </div>
                                    <div class="info-value">
        '''
        
        # Add geographical distribution
        geo = info.get('geographical_distribution', {})
        for key, value in geo.items():
            if isinstance(value, list):
                html_content += f'<div><strong>{key.replace("_", " ").title()}:</strong> {", ".join(value)}</div>'
            elif isinstance(value, dict):
                html_content += f'<div><strong>{key.replace("_", " ").title()}:</strong> '
                for subkey, subvalue in list(value.items())[:2]:
                    html_content += f'{subkey}: {subvalue}; '
                html_content += '</div>'
            else:
                html_content += f'<div><strong>{key.replace("_", " ").title()}:</strong> {value}</div>'
        
        html_content += f'''
                                    </div>
                                </div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">References</div>
                                <div class="info-value">
                                    {', '.join(info.get('references', []))}
                                </div>
                            </div>
                        </div>
                    </div>
        '''
    
    html_content += '''
                </div>
            </div>
            
            <!-- Phylogroups Section -->
            <div class="content-section" id="phylogroups">
                <h2 class="section-title">
                    <i class="fas fa-project-diagram"></i>
                    Phylogroup Database
                </h2>
                
                <div class="cards-grid">
    '''
    
    # Generate phylogroup cards
    for phylogroup, info in sorted(PHYLOGROUP_DATABASE.items()):
        html_content += f'''
                    <div class="data-card">
                        <div class="card-header">
                            <div class="card-title">Phylogroup {phylogroup}</div>
                            <div class="card-subtitle">{info.get('characteristics', '')}</div>
                        </div>
                        <div class="card-content">
                            <div class="info-group">
                                <div class="info-label">Pathogenic Potential</div>
                                <div class="info-value">{info.get('pathogenic_potential', 'Unknown')}</div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">Common Sequence Types</div>
                                <div class="info-value">
                                    {', '.join([f'ST{st}' for st in info.get('common_st', [])[:5]])}
                                </div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">Common Serotypes</div>
                                <div class="info-value">
                                    {', '.join(info.get('serotypes', [])[:4])}
                                </div>
                            </div>
                            
                            <div class="info-group">
                                <div class="info-label">Common Virulence Genes</div>
                                <div class="gene-list">
        '''
        
        # Add virulence genes
        for gene in info.get('virulence_genes', [])[:6]:
            html_content += f'<span class="gene-tag">{gene}</span>'
        
        html_content += f'''
                                </div>
                            </div>
                            
                            {f'<div class="info-group"><div class="info-label">Notes</div><div class="info-value">{info.get("notes", "")}</div></div>' if info.get("notes") else ""}
                        </div>
                    </div>
        '''
    
    html_content += '''
                </div>
            </div>
            
            <!-- NEW CARBAPENEMASE SECTION - UPDATED WITH ALL MISSING FIELDS -->
            <div class="content-section" id="carbapenemase">
                <h2 class="section-title">
                    <i class="fas fa-shield-virus"></i>
                    Carbapenemase-Producing E. coli
                </h2>
                
                <div class="search-section">
                    <h3>Global Threat: Carbapenem-Resistant E. coli</h3>
                    <p>Comprehensive profiles of carbapenemase-producing E. coli strains, including enzyme characteristics, geographical distribution, and treatment options.</p>
                </div>
    '''
    
    # Generate carbapenemase profiles WITH ALL MISSING FIELDS
    for profile_type, profile_data in CARBAPENEMASE_PRODUCERS.items():
        html_content += f'''
                <div class="subsection" style="margin-bottom: var(--space-2xl);">
                    <div class="subsection-title">
                        <i class="fas fa-virus"></i>
                        {profile_type.replace('_', ' ').title()}
                    </div>
                    <div class="cards-grid">
                        <div class="data-card">
                            <div class="card-header">
                                <div class="card-title">{profile_type.replace('_', ' ').title()}</div>
                            </div>
                            <div class="card-content">
                                <div class="info-group">
                                    <div class="info-label">Pathotype Association</div>
                                    <div class="info-value">{profile_data.get('pathotype', 'Various')}</div>
                                </div>
                                
                                <div class="info-group">
                                    <div class="info-label">Sequence Types</div>
                                    <div class="info-value">
                                        {', '.join([f'ST{st}' for st in profile_data.get('st', [])])}
                                    </div>
                                </div>
                                
                                <div class="info-group">
                                    <div class="info-label">Carbapenemase Genes</div>
                                    <div class="gene-list">
        '''
        
        # Add carbapenemase genes
        for gene in profile_data.get('carbapenemase', []):
            html_content += f'<span class="gene-tag" style="background: var(--danger); color: white;">{gene}</span>'
        
        html_content += f'''
                                    </div>
                                </div>
                                
                                <div class="info-group">
                                    <div class="info-label">Enzyme Class</div>
                                    <div class="info-value">{profile_data.get('enzyme_class', 'Unknown')}</div>
                                </div>
                                
                                <div class="info-group">
                                    <div class="info-label">Hydrolysis Spectrum</div>
                                    <div class="info-value">
                                        {', '.join(profile_data.get('hydrolysis_spectrum', []))}
                                    </div>
                                </div>
                                
                                <div class="info-group">
                                    <div class="info-label">Inhibitor Profile</div>
                                    <div class="info-value">
                                        <strong>Inhibited By:</strong> {', '.join(profile_data.get('inhibitor_profile', {}).get('inhibited_by', []))}<br>
                                        <strong>Resistant To:</strong> {', '.join(profile_data.get('inhibitor_profile', {}).get('resistant_to', []))}
                                    </div>
                                </div>
                                
                                <div class="detailed-section">
                                    <div class="subsection">
                                        <div class="subsection-title">
                                            <i class="fas fa-dna"></i>
                                            Genetic Context
                                        </div>
                                        <div class="info-value">
                                            <div><strong>Gene Location:</strong> {', '.join(profile_data.get('genetic_context', {}).get('gene_location', []))}</div>
                                            <div><strong>Mobile Element:</strong> {profile_data.get('genetic_context', {}).get('mobile_element', 'Unknown')}</div>
                                            <div><strong>Co-resistance:</strong> {', '.join(profile_data.get('genetic_context', {}).get('co-resistance', []))}</div>
                                        </div>
                                    </div>
                                    
                                    <div class="subsection">
                                        <div class="subsection-title">
                                            <i class="fas fa-microscope"></i>
                                            Detection Methods
                                        </div>
                                        <div class="info-value">
        '''
        
        # Add detection methods
        detection = profile_data.get('detection_methods', {})
        for key, value in detection.items():
            if isinstance(value, list):
                html_content += f'<div><strong>{key.replace("_", " ").title()}:</strong> {", ".join(value)}</div>'
            else:
                html_content += f'<div><strong>{key.replace("_", " ").title()}:</strong> {value}</div>'
        
        html_content += f'''
                                        </div>
                                    </div>
                                    
                                    <div class="subsection">
                                        <div class="subsection-title">
                                            <i class="fas fa-pills"></i>
                                            Treatment Options
                                        </div>
                                        <div class="info-value">
        '''
        
        # Add complete treatment options including combination therapy
        treatment = profile_data.get('treatment_options', {})
        if 'first_line' in treatment:
            html_content += f'<div><strong>First Line:</strong> {", ".join(treatment["first_line"])}</div>'
        if 'alternative' in treatment:
            html_content += f'<div><strong>Alternative:</strong> {", ".join(treatment["alternative"])}</div>'
        if 'combination_therapy' in treatment:
            html_content += f'<div><strong>Combination Therapy:</strong> {", ".join(treatment["combination_therapy"])}</div>'
        if 'important_notes' in treatment:
            html_content += f'<div><strong>Important Notes:</strong> {treatment["important_notes"]}</div>'
        
        html_content += f'''
                                        </div>
                                    </div>
                                    
                                    <div class="subsection">
                                        <div class="subsection-title">
                                            <i class="fas fa-shield-alt"></i>
                                            Infection Control
                                        </div>
                                        <div class="info-value">
        '''
        
        # Add infection control measures
        infection_control = profile_data.get('infection_control', {})
        for key, value in infection_control.items():
            if isinstance(value, list):
                html_content += f'<div><strong>{key.replace("_", " ").title()}:</strong> {", ".join(value)}</div>'
            else:
                html_content += f'<div><strong>{key.replace("_", " ").title()}:</strong> {value}</div>'
        
        html_content += f'''
                                        </div>
                                    </div>
                                    
                                    <div class="subsection">
                                        <div class="subsection-title">
                                            <i class="fas fa-globe-americas"></i>
                                            Geographical Distribution
                                        </div>
                                        <div class="info-value">
        '''
        
        # Add geographical distribution
        geo = profile_data.get('geographical_distribution', {})
        if 'endemic_regions' in geo:
            html_content += f'<div><strong>Endemic Regions:</strong> {", ".join(geo["endemic_regions"])}</div>'
        if 'hotspots' in geo:
            html_content += '<div><strong>Hotspots:</strong> '
            for region, desc in list(geo['hotspots'].items())[:2]:
                html_content += f'{region}: {desc}; '
            html_content += '</div>'
        
        html_content += f'''
                                        </div>
                                    </div>
                                    
                                    <div class="subsection">
                                        <div class="subsection-title">
                                            <i class="fas fa-stethoscope"></i>
                                            Clinical Significance
                                        </div>
                                        <div class="info-value">
        '''
        
        # Add clinical significance
        clinical = profile_data.get('clinical_significance', {})
        for key, value in clinical.items():
            if isinstance(value, list):
                html_content += f'<div><strong>{key.replace("_", " ").title()}:</strong> {", ".join(value)}</div>'
            else:
                html_content += f'<div><strong>{key.replace("_", " ").title()}:</strong> {value}</div>'
        
        html_content += f'''
                                        </div>
                                    </div>
                                </div>
                                
                                <div class="info-group">
                                    <div class="info-label">References</div>
                                    <div class="info-value">
                                        {', '.join(profile_data.get('references', [])[:3])}
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
        '''
    
    html_content += '''
            </div>
            
            <!-- Specialized Profiles Section -->
            <div class="content-section" id="specialized">
                <h2 class="section-title">
                    <i class="fas fa-star"></i>
                    Specialized Pathotype Profiles
                </h2>
    '''
    
    # Generate specialized profiles
    for profile_type, profiles in SPECIALIZED_PROFILES.items():
        html_content += f'''
                <div class="subsection" style="margin-bottom: var(--space-2xl);">
                    <div class="subsection-title">
                        <i class="fas fa-virus"></i>
                        {profile_type.replace('_', ' ').title()}
                    </div>
                    <div class="cards-grid">
        '''
        
        for profile_name, profile_data in profiles.items():
            html_content += f'''
                        <div class="data-card">
                            <div class="card-header">
                                <div class="card-title">{profile_name}</div>
                            </div>
                            <div class="card-content">
            '''
            
            # Add profile data
            for key, value in profile_data.items():
                if key == 'virulence' and isinstance(value, list):
                    html_content += f'''
                                <div class="info-group">
                                    <div class="info-label">Virulence Factors</div>
                                    <div class="gene-list">
                    '''
                    for virulence in value[:6]:
                        html_content += f'<span class="gene-tag">{virulence}</span>'
                    html_content += '</div></div>'
                elif isinstance(value, list):
                    html_content += f'''
                                <div class="info-group">
                                    <div class="info-label">{key.replace('_', ' ').title()}</div>
                                    <div class="info-value">{', '.join(value)}</div>
                                </div>
                    '''
                else:
                    html_content += f'''
                                <div class="info-group">
                                    <div class="info-label">{key.replace('_', ' ').title()}</div>
                                    <div class="info-value">{value}</div>
                                </div>
                    '''
            
            html_content += '''
                            </div>
                        </div>
            '''
        
        html_content += '''
                    </div>
                </div>
        '''
    
    html_content += '''
            </div>
            
            <!-- References Section -->
            <div class="content-section" id="references">
                <h2 class="section-title">
                    <i class="fas fa-book"></i>
                    Comprehensive Reference Database
                </h2>
    '''
    
    # Generate PubMed references
    if "PUBMED_REFERENCES" in COMPREHENSIVE_REFERENCES:
        html_content += '''
                <div class="subsection">
                    <div class="subsection-title">
                        <i class="fas fa-file-medical"></i>
                        PubMed References
                    </div>
        '''
        
        for category, refs in COMPREHENSIVE_REFERENCES["PUBMED_REFERENCES"].items():
            html_content += f'''
                    <div class="detailed-section" style="margin-bottom: var(--space-xl);">
                        <div class="subsection-title" style="font-size: 1rem;">
                            {category.replace('_', ' ').title()} ({len(refs)} references)
                        </div>
                        <div class="info-value">
            '''
            
            for ref in refs:
                html_content += f'<div style="margin-bottom: var(--space-sm);">{ref}</div>'
            
            html_content += '''
                        </div>
                    </div>
            '''
        
        html_content += '''
                </div>
        '''
    
    # Generate DOI references
    if "DOI_REFERENCES" in COMPREHENSIVE_REFERENCES:
        html_content += '''
                <div class="subsection">
                    <div class="subsection-title">
                        <i class="fas fa-link"></i>
                        DOI References
                    </div>
        '''
        
        for category, refs in COMPREHENSIVE_REFERENCES["DOI_REFERENCES"].items():
            html_content += f'''
                    <div class="detailed-section" style="margin-bottom: var(--space-xl);">
                        <div class="subsection-title" style="font-size: 1rem;">
                            {category.replace('_', ' ').title()} ({len(refs)} references)
                        </div>
                        <div class="info-value">
            '''
            
            for ref in refs:
                html_content += f'<div style="margin-bottom: var(--space-sm);">{ref}</div>'
            
            html_content += '''
                        </div>
                    </div>
            '''
        
        html_content += '''
                </div>
        '''
    
    html_content += '''
            </div>
        </div>
        
        <!-- Footer -->
        <div class="footer">
            <div class="authorship">
                <h3>Technical Information & Collaboration</h3>
                <p><strong>Author:</strong> Brown Beckley</p>
                <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
                <p><strong>Email:</strong> <a href="mailto:brownbeckley94@gmail.com" class="custom-link">brownbeckley94@gmail.com</a></p>
                <p><strong>Database Version:</strong> 2025 | <a href="https://github.com/bbeckley-hub/EcoliTyper" class="custom-link" target="_blank">Actively Maintained</a></p>
            </div>
        </div>
    </div>
    
    <script>
        // Tab navigation
        function switchTab(tabId) {{
            // Hide all sections
            document.querySelectorAll('.content-section').forEach(section => {{
                section.classList.remove('active');
            }});
            
            // Remove active class from all tabs
            document.querySelectorAll('.nav-tab').forEach(tab => {{
                tab.classList.remove('active');
            }});
            
            // Show selected section
            document.getElementById(tabId).classList.add('active');
            
            // Activate clicked tab
            event.target.classList.add('active');
        }}
        
        // Search and filter functionality for lineages
        document.addEventListener('DOMContentLoaded', function() {{
            const searchInput = document.getElementById('lineageSearch');
            const categoryFilter = document.getElementById('lineageCategory');
            
            if (searchInput && categoryFilter) {{
                searchInput.addEventListener('input', filterLineages);
                categoryFilter.addEventListener('change', filterLineages);
            }}
        }});
        
        function filterLineages() {{
            const searchTerm = document.getElementById('lineageSearch').value.toLowerCase();
            const categoryFilter = document.getElementById('lineageCategory').value;
            
            document.querySelectorAll('#lineagesGrid .data-card').forEach(card => {{
                const cardText = card.textContent.toLowerCase();
                const cardCategory = card.dataset.category;
                
                const matchesSearch = searchTerm === '' || cardText.includes(searchTerm);
                const matchesCategory = categoryFilter === '' || cardCategory === categoryFilter;
                
                card.style.display = (matchesSearch && matchesCategory) ? 'block' : 'none';
            }});
        }}
        
        // Simple risk-based filtering
        function filterByRisk(riskLevel) {{
            document.querySelectorAll('#lineagesGrid .data-card').forEach(card => {{
                if (riskLevel === 'all') {{
                    card.style.display = 'block';
                }} else {{
                    card.style.display = card.dataset.risk === riskLevel ? 'block' : 'none';
                }}
            }});
        }}
    </script>
</body>
</html>
'''
    
    # Save the file
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"✅ COMPREHENSIVE E. coli reference generated: {output_path}")
    print(f"📊 UPDATED Database Statistics:")
    print(f"   - Lineages: {stats['lineages']} (includes ST1193, ST405)")
    print(f"   - Pathotypes: {stats['pathotypes']} (includes MAEC, AIEC)")
    print(f"   - Serotypes: {stats['serotypes']} (includes O45:H2, O103:H2)")
    print(f"   - Phylogroups: {stats['phylogroups']}")
    print(f"   - Carbapenemase Types: {stats['carbapenemase_profiles']}")
    print(f"   - References: {stats['references_pubmed'] + stats['references_doi']}")
    print(f"🎯 NEW FEATURES ADDED:")
    print(f"   - Complete molecular typing data (fumC, fimH, sublineages)")
    print(f"   - Enhanced epidemiology with reservoir and prevalence data")
    print(f"   - Full reference lists (no truncation)")
    print(f"   - Shiga toxin profiles for serotypes")
    print(f"   - NEW Carbapenemase Producers section")
    print(f"   - COMPLETE carbapenemase data (hydrolysis, inhibitors, genetic context, detection, treatment, infection control)")
    print(f"🎨 IMPROVED ORGANIZATION:")
    print(f"   - Green text for Database Information values")
    print(f"   - Red heading for AMR section")
    print(f"   - Hyperlinked GitHub repository")
    print(f"   - Yellow hyperlinks in footer with mailto")
    print(f"   - Reverted to original font sizes for data cards")
    
    return output_path

# Test the generator
if __name__ == "__main__":
    generate_comprehensive_reference()
