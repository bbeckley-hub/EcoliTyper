#!/usr/bin/env python3
"""
EcoliLineageDB: Comprehensive Escherichia coli Lineage & Pathotype Database
Complete E. coli sequence types, serotypes, phylogroups, pathotypes and references
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
Scientific Review: Based on current literature and established typing schemes
Date: 2025
Send a quick mail for any issues or further explanations.
"""

# =============================================================================
# COMPREHENSIVE ESCHERICHIA COLI LINEAGE DATABASE 
# =============================================================================

LINEAGE_DATABASE = {
    # =========================================================================
    # HIGH-RISK CLONES & GLOBAL LINEAGES
    # =========================================================================
    
    "ST131": {
        "primary_name": "Global MDR Pandemic Clone",
        "category": "High-Risk Clone",
        "sublineages": ["ST131-H30", "ST131-H30-R", "ST131-H30-Rx"],
        "serotype": "O25b:H4",
        "phylogroup": "B2",
        "fumC": "fumC40",
        "fimH": "fimH30",
        "clermont_complex": "STc131",
        
        "pathotypes": ["UPEC", "SEPEC", "ExPEC"],
        "key_virulence_genes": ["fimH", "papA", "iutA", "fyuA", "kpsMII", "ompT", "malX"],
        
        "resistance_profile": {
            "esbl": ["blaCTX-M-15", "blaCTX-M-14", "blaCTX-M-27"],
            "fluoroquinolone": ["gyrA-S83L", "gyrA-D87N", "parC-S80I", "parC-E84V"],
            "aminoglycoside": ["aac(6')-Ib-cr", "aac(3)-II"],
            "other": ["sul1", "sul2", "tet(A)", "dfrA17"]
        },
        
        "epidemiology": {
            "first_identified": 2008,
            "global_distribution": "Worldwide pandemic",
            "geographical_distribution": {
                "high_prevalence": ["North America", "Europe", "Asia", "South America", "Australia"],
                "medium_prevalence": ["Africa", "Middle East"],
                "regional_variants": {
                    "North America": "H30-Rx dominant, CTX-M-15",
                    "Europe": "CTX-M-14/15/27 variants",
                    "Asia": "CTX-M-14/27 common, emerging in hospitals",
                    "South America": "Increasing prevalence in healthcare",
                    "Africa": "Emerging threat, limited surveillance"
                }
            },
            "transmission": "Healthcare-associated, community-acquired",
            "reservoir": "Human gastrointestinal tract"
        },
        
        "clinical_significance": {
            "primary_infections": ["UTI", "Bloodstream infections", "Prostatitis"],
            "mortality": "Increased mortality in bacteremia",
            "treatment_challenges": "Multidrug resistance limits options"
        },
        
        "genomic_features": {
            "plasmid": ["pEK499-like", "IncFII"],
            "prophage": ["pheV", "yehV"],
            "island": ["pks island (colibactin)"]
        },
        
        "risk_level": "VERY HIGH",
        "outbreak_potential": "HIGH",
        "key_references": ["PMID: 28344773", "PMID: 27025834", "PMID: 35054380", "PMID: 21081548", "PMID: 23926176", "PMID: 35927655", "PMID: 25241262", "PMID: 31299909"]
    },
    
    "ST73": {
        "primary_name": "Classic Uropathogenic Clone",
        "category": "High-Risk Clone", 
        "sublineages": ["ST73-1", "ST73-2"],
        "serotype": "O6:H1",
        "phylogroup": "B2",
        "fumC": "fumC32",
        "fimH": "fimH1",
        "clermont_complex": "STc73",
        
        "pathotypes": ["UPEC", "SEPEC"],
        "key_virulence_genes": ["fimH", "papA", "papGII", "hlyA", "cnf1", "iutA", "fyuA"],
        
        "resistance_profile": {
            "common": ["Ampicillin", "Tetracycline"],
            "emerging": ["Fluoroquinolones", "Third-gen cephalosporins"],
            "genes": ["blaTEM", "tet(A)", "tet(B)"]
        },
        
        "epidemiology": {
            "distribution": "Global, community-associated",
            "geographical_distribution": {
                "high_prevalence": ["Europe", "North America", "Australia"],
                "medium_prevalence": ["Asia", "South America"],
                "low_prevalence": ["Africa"],
                "regional_features": "Community-associated UTIs, less hospital adaptation than ST131"
            },
            "prevalence": "Common UPEC lineage",
            "age_distribution": "All age groups"
        },
        
        "clinical_significance": {
            "primary_infections": ["Cystitis", "Pyelonephritis"],
            "virulence": "High virulence in UTI models",
            "persistence": "Forms intracellular bacterial communities"
        },
        
        "risk_level": "HIGH",
        "outbreak_potential": "MODERATE",
        "key_references": ["PMID: 25729395", "PMID: 37471138", "PMID: 39303870", "PMID: 30810518", "PMID: 37471138", "PMID: 25355761", "PMID: 32120035"]
    },
    
    "ST95": {
        "primary_name": "Neonatal Meningitis/ExPEC Clone",
        "category": "High-Risk Clone",
        "sublineages": ["ST95-1", "ST95-2"],
        "serotype": "O18:H7",
        "phylogroup": "B2", 
        "fumC": "fumC28",
        "fimH": "fimH5",
        "clermont_complex": "STc95",
        
        "pathotypes": ["NMEC", "UPEC", "SEPEC"],
        "key_virulence_genes": ["fimH", "ibeA", "kpsMII-K1", "cnf1", "hlyA", "iutA"],
        
        "resistance_profile": {
            "common": ["Ampicillin", "Tetracycline"],
            "genes": ["blaTEM", "tet(A)"],
            "notes": "Generally less resistant than ST131"
        },
        
        "epidemiology": {
            "distribution": "Global",
            "geographical_distribution": {
                "high_prevalence": ["Europe", "North America", "East Asia"],
                "medium_prevalence": ["South America", "Australia"],
                "emerging_regions": ["Africa", "South Asia"],
                "neonatal_focus": "Major cause of neonatal meningitis worldwide"
            },
            "specialization": "Neonatal meningitis, UTI",
            "capsule": "K1 capsule predominant"
        },
        
        "clinical_significance": {
            "primary_infections": ["Neonatal meningitis", "UTI", "Bacteremia"],
            "mortality": "High in neonatal meningitis",
            "sequelae": "Neurological deficits common"
        },
        
        "risk_level": "VERY HIGH",
        "outbreak_potential": "MODERATE",
        "key_references": ["PMID: 38622998", "PMID: 34560248", "PMID: 19583828", "PMID: 35076267", "PMID: 25037925", "PMID: 36558824", "PMID: 38383596"]
    },
    
    "ST69": {
        "primary_name": "Uropathogenic/Diarrheagenic Hybrid",
        "category": "High-Risk Clone",
        "sublineages": ["ST69-1", "ST69-2"],
        "serotype": "O17/O77:H18", 
        "phylogroup": "D",
        "fumC": "fumC14",
        "fimH": "fimH27",
        "clermont_complex": "STc69",
        
        "pathotypes": ["UPEC", "EAEC", "ExPEC"],
        "key_virulence_genes": ["fimH", "aggR", "aap", "aatA", "iutA", "fyuA"],
        
        "resistance_profile": {
            "common": ["Ampicillin", "Tetracycline", "Sulfonamides"],
            "genes": ["blaTEM", "tet(A)", "sul1", "sul2"],
            "emerging": ["Fluoroquinolones"]
        },
        
        "epidemiology": {
            "distribution": "Global",
            "geographical_distribution": {
                "high_prevalence": ["North America", "Europe"],
                "medium_prevalence": ["Asia", "Australia"],
                "regional_variants": {
                    "North America": "Common in community UTIs",
                    "Europe": "Increasing in healthcare settings",
                    "Asia": "Emerging with hybrid virulence"
                }
            },
            "hybrid_nature": "UPEC with EAEC virulence factors",
            "transmission": "Community and healthcare"
        },
        
        "clinical_significance": {
            "primary_infections": ["UTI", "Diarrhea"],
            "virulence": "Enhanced adherence capabilities",
            "treatment": "May require combination therapy"
        },
        
        "risk_level": "HIGH",
        "outbreak_potential": "MODERATE",
        "key_references": ["PMID: 37764013", "PMID: 18199778", "PMID: 32634601", "PMID: 38318209", "PMID: 31477018"]
    },
    
    # =========================================================================
    # EMERGING RESISTANCE-ASSOCIATED LINEAGES
    # =========================================================================
    
    "ST410": {
        "primary_name": "Emerging MDR Clone",
        "category": "Resistance-Associated",
        "sublineages": ["ST410-A", "ST410-B", "ST410-C"],
        "serotype": "O8:H9", 
        "phylogroup": "C",
        "fumC": "fumC23",
        "fimH": "fimH24",
        "clermont_complex": "STc410",
        
        "pathotypes": ["UPEC", "SEPEC", "ExPEC"],
        "key_virulence_genes": ["fimH", "iutA", "fyuA", "kpsMII"],
        
        "resistance_profile": {
            "esbl": ["blaCTX-M-15", "blaCTX-M-27"],
            "carbapenemase": ["blaOXA-181", "blaNDM-5"],
            "fluoroquinolone": ["gyrA mutations", "parC mutations"],
            "other": ["aac(6')-Ib-cr", "rmtB"]
        },
        
        "epidemiology": {
            "emergence": "Rapid global spread since 2010",
            "geographical_distribution": {
                "high_prevalence": ["Europe", "Asia", "Middle East"],
                "medium_prevalence": ["North America", "Africa", "South America"],
                "hotspots": {
                    "Germany": "Early emergence, high prevalence",
                    "China": "Rapid dissemination, NDM variants",
                    "India": "High carbapenem resistance",
                    "United States": "Increasing in healthcare"
                }
            },
            "distribution": "Worldwide, increasing prevalence",
            "transmission": "Healthcare and community"
        },
        
        "clinical_significance": {
            "primary_infections": ["UTI", "Bloodstream infections"],
            "resistance": "Extensive drug resistance",
            "treatment": "Limited therapeutic options"
        },
        
        "risk_level": "VERY HIGH",
        "outbreak_potential": "HIGH",
        "key_references": ["PMID: 38193668", "PMID: 34925277", "PMID: 31176748", "PMID: 30021879", " PMID: 39324056", "PMID: 37107633"]
    },
    
    "ST648": {
        "primary_name": "Zoonotic MDR Clone", 
        "category": "Resistance-Associated",
        "sublineages": ["ST648-1", "ST648-2"],
        "serotype": "O8:H9",
        "phylogroup": "F",
        "fumC": "fumC26",
        "fimH": "fimH29",
        "clermont_complex": "STc648",
        
        "pathotypes": ["UPEC", "SEPEC", "ExPEC", "AVEC"],
        "key_virulence_genes": ["fimH", "iutA", "fyuA", "iss", "hlyF"],
        
        "resistance_profile": {
            "esbl": ["blaCTX-M-14", "blaCTX-M-15", "blaCTX-M-55"],
            "carbapenemase": ["blaNDM-1", "blaNDM-5", "blaOXA-48"],
            "fluoroquinolone": ["qnrS1", "aac(6')-Ib-cr"],
            "other": ["rmtB", "tet(X4)"]
        },
        
        "epidemiology": {
            "zoonotic_potential": "Human, poultry, environmental",
            "geographical_distribution": {
                "high_prevalence": ["Asia", "Europe", "Middle East"],
                "medium_prevalence": ["Africa", "North America"],
                "animal_reservoirs": {
                    "Poultry": "Global poultry industry",
                    "Pigs": "Increasing reports",
                    "Wild birds": "Environmental spread"
                },
                "regional_variants": {
                    "China": "NDM-5 dominant, tet(X4) emerging",
                    "Europe": "OXA-48, CTX-M-15 common",
                    "Middle East": "NDM-1, OXA-48 variants"
                }
            },
            "distribution": "Global, particularly Asia and Europe",
            "transmission": "Foodborne, environmental, healthcare"
        },
        
        "clinical_significance": {
            "primary_infections": ["UTI", "Bloodstream infections"],
            "resistance": "Pan-drug resistance emerging",
            "one_health": "Significant zoonotic concern"
        },
        
        "risk_level": "VERY HIGH",
        "outbreak_potential": "HIGH",
        "key_references": ["PMID: 30885899", "PMID: 24398338", "PMID: 40943570", "PMID: 36708770", "PMID: 37293536"]
    },
    
    # =========================================================================
    # DIARRHEAGENIC-ASSOCIATED LINEAGES  
    # =========================================================================
    
    "ST11": {
        "primary_name": "EHEC O157:H7 Clone",
        "category": "Diarrheagenic",
        "sublineages": ["ST11-1", "ST11-2"],
        "serotype": "O157:H7",
        "phylogroup": "E",
        "fumC": "fumC12",
        "fimH": "fimH9",
        "clermont_complex": "STc11",
        
        "pathotypes": ["EHEC", "STEC"],
        "key_virulence_genes": ["stx1", "stx2", "eae", "ehxA", "tccP", "espA"],
        
        "resistance_profile": {
            "common": ["Streptomycin", "Sulfonamides", "Tetracycline"],
            "genes": ["strA-strB", "sul2", "tet(B)"],
            "important_note": "Avoid antibiotics due to HUS risk"
        },
        
        "epidemiology": {
            "distribution": "Global, foodborne outbreaks",
            "geographical_distribution": {
                "high_incidence": ["North America", "Europe", "Japan", "Argentina"],
                "medium_incidence": ["Australia", "South Africa", "Chile"],
                "emerging_regions": ["Africa", "Asia"],
                "outbreak_patterns": {
                    "North America": "Beef-associated outbreaks",
                    "Europe": "Produce, dairy outbreaks", 
                    "Japan": "Multiple food vehicles",
                    "Argentina": "Beef-dominated cuisine high risk"
                }
            },
            "reservoir": "Cattle gastrointestinal tract",
            "transmission": "Undercooked beef, produce, water"
        },
        
        "clinical_significance": {
            "primary_infections": ["Hemorrhagic colitis", "HUS"],
            "complications": "Hemolytic uremic syndrome (5-15%)",
            "mortality": "Low for diarrhea, high for HUS"
        },
        
        "risk_level": "VERY HIGH",
        "outbreak_potential": "VERY HIGH",
        "key_references": [" PMID: 36144444", "PMID: 22355013", "PMID: 32946720", "PMID: 19245340", "PMID: 25819955", "PMID: 26311863", "PMID: 39083224"]
    },
    
    "ST10": {
        "primary_name": "Diverse Commensal/Pathogenic Clone",
        "category": "Commensal/Pathogenic",
        "sublineages": ["ST10-1", "ST10-2", "ST10-3"],
        "serotype": "O2:H6",
        "phylogroup": "A",
        "fumC": "fumC10", 
        "fimH": "fimH54",
        "clermont_complex": "STc10",
        
        "pathotypes": ["Commensal", "UPEC", "ETEC", "EAEC"],
        "key_virulence_genes": ["fimH", "eltA", "estA", "aggR", "aap"],
        
        "resistance_profile": {
            "common": ["Ampicillin", "Tetracycline"],
            "genes": ["blaTEM", "tet(A)"],
            "emerging": ["ESBLs", "Fluoroquinolones"]
        },
        
        "epidemiology": {
            "distribution": "Global, highly diverse",
            "geographical_distribution": {
                "ubiquitous": "Worldwide distribution in humans, animals, environment",
                "human_commensal": "All continents, all populations",
                "clinical_importance": "Carries diverse virulence factors globally",
                "regional_variants": {
                    "Asia": "Often carries ETEC virulence factors",
                    "Africa": "Common commensal, emerging pathogen",
                    "Europe/N America": "UPEC and commensal variants"
                }
            },
            "ecology": "Human and animal gastrointestinal",
            "pathogenic_potential": "Can acquire virulence plasmids"
        },
        
        "clinical_significance": {
            "spectrum": "Commensal to pathogenic",
            "infections": ["UTI", "Diarrhea when virulent"],
            "importance": "Genetic background for horizontal gene transfer"
        },
        
        "risk_level": "LOW to MODERATE",
        "outbreak_potential": "LOW",
        "key_references": [" PMID: 32985020", "PMID: 39382271", "PMID: 35000591", "PMID: 40591063", "PMID: 38072242", "PMID: 29897467", "PMID: 32985020"]
    },

    
    # =========================================================================
    # ANIMAL-ASSOCIATED LINEAGES
    # =========================================================================
    
    "ST117": {
        "primary_name": "Avian Pathogenic E. coli (APEC)",
        "category": "Animal-Associated",
        "sublineages": ["ST117-1", "ST117-2"],
        "serotype": "O78:H4",
        "phylogroup": "B1",
        "fumC": "fumC17",
        "fimH": "fimH5",
        "clermont_complex": "STc117",
        
        "pathotypes": ["APEC", "UPEC", "ExPEC"],
        "key_virulence_genes": ["iss", "iroN", "iutA", "hlyF", "ompT", "vat"],
        
        "resistance_profile": {
            "common": ["Tetracycline", "Sulfonamides", "Ampicillin"],
            "genes": ["tet(A)", "sul1", "sul2", "blaTEM"],
            "emerging": ["Fluoroquinolones", "Third-gen cephalosporins"]
        },
        
        "epidemiology": {
            "host_range": "Poultry primary, zoonotic potential",
            "geographical_distribution": {
                "poultry_global": "Worldwide in poultry production",
                "high_prevalence": ["Europe", "North America", "Asia", "Brazil"],
                "zoonotic_cases": "Sporadic human infections worldwide",
                "poultry_industry_impact": {
                    "Global": "Major economic losses",
                    "Europe": "Widespread in broilers and layers",
                    "North America": "Commercial and backyard flocks",
                    "Asia": "Rapidly expanding poultry industry"
                }
            },
            "distribution": "Global poultry industry",
            "economic_impact": "Major losses in poultry production"
        },
        
        "clinical_significance": {
            "animal_disease": "Colibacillosis, septicemia",
            "human_infections": "Rare, usually immunocompromised",
            "one_health": "Antibiotic use in agriculture drives resistance"
        },
        
        "risk_level": "MODERATE",
        "outbreak_potential": "HIGH in poultry",
        "key_references": ["PMID: 35945191", "PMID: 39235965", "PMID: 38175844", "PMID: 36755458", "PMID: 22615330", "PMID: 28049430", "PMID: 34047644", "PMID: 28049430"]
    },
    
    "ST88": {
        "primary_name": "Meningitis/ExPEC Clone",
        "category": "High-Risk Clone",
        "sublineages": ["ST88-1", "ST88-2"],
        "serotype": "O45:H7",
        "phylogroup": "B2",
        "fumC": "fumC35",
        "fimH": "fimH8",
        "clermont_complex": "STc88",
        
        "pathotypes": ["NMEC", "UPEC", "SEPEC"],
        "key_virulence_genes": ["fimH", "ibeA", "kpsMII", "cnf1", "hlyA"],
        
        "resistance_profile": {
            "common": ["Ampicillin", "Tetracycline"],
            "genes": ["blaTEM", "tet(A)"],
            "notes": "Variable resistance patterns"
        },
        
        "epidemiology": {
            "distribution": "Global, less common than ST95",
            "geographical_distribution": {
                "high_prevalence": ["Europe", "North America"],
                "medium_prevalence": ["Asia", "Africa"],
                "emerging_regions": "Increasing reports globally",
                "clinical_focus": "Neonatal infections and UTIs"
            },
            "specialization": "Neonatal infections",
            "capsule": "Often K1 capsule"
        },
        
        "risk_level": "HIGH",
        "outbreak_potential": "MODERATE",
        "key_references": ["PMID: 19307211", "PMID: 26865687", "PMID: 26803720", "PMID: 37370379", "PMID: 36672813", "PMID: 36246985", "PMID: 37578342"]
    },

# =========================================================================
# EMERGING & GLOBAL MDR LINEAGES - COMPREHENSIVE ADDITIONS
# =========================================================================

"ST1193": {
        "primary_name": "Emerging Fluoroquinolone-Resistant Uropathogenic Clone",
        "category": "High-Risk Clone",
        "sublineages": ["ST1193-1", "ST1193-2", "ST1193-3"],
        "serotype": "O75:H5",
        "phylogroup": "B2",
        "fumC": "fumC41",
        "fimH": "fimH64",
        "clermont_complex": "STc14",
        
        "pathotypes": ["UPEC", "SEPEC", "ExPEC"],
        "key_virulence_genes": ["fimH", "papA", "papGII", "iutA", "fyuA", "kpsMII", "ompT", "malX"],
        
        "resistance_profile": {
            "fluoroquinolone": ["gyrA-S83L", "gyrA-D87N", "parC-S80I", "parC-E84V"],
            "esbl": ["blaCTX-M-15", "blaCTX-M-27"],
            "aminoglycoside": ["aac(6')-Ib-cr", "aac(3)-II"],
            "other": ["sul1", "sul2", "tet(A)", "dfrA17"],
            "important_note": "Rapidly emerging fluoroquinolone resistance in community settings"
        },
        
        "epidemiology": {
            "first_identified": 2015,
            "emergence_timeline": "Rapid global spread since 2015",
            "geographical_distribution": {
                "high_prevalence": ["North America", "Europe", "Australia"],
                "medium_prevalence": ["Asia", "South America"],
                "emerging_regions": ["Africa", "Middle East"],
                "regional_variants": {
                    "North America": "Community-associated UTIs, nursing homes",
                    "Europe": "Healthcare and community settings",
                    "Australia": "Rapid dissemination in community",
                    "Asia": "Emerging in hospital settings"
                }
            },
            "transmission": "Community-associated, healthcare-associated",
            "reservoir": "Human gastrointestinal tract",
            "demographics": "Elderly, immunocompromised, recurrent UTI patients"
        },
        
        "clinical_significance": {
            "primary_infections": ["UTI", "Recurrent UTI", "Bloodstream infections"],
            "treatment_challenges": "Fluoroquinolone resistance limits oral options",
            "persistence": "High recurrence rates in urinary tract",
            "mortality": "Increased in bacteremia compared to susceptible strains"
        },
        
        "genomic_features": {
            "plasmid": ["IncFII", "IncFIB"],
            "prophage": ["yehV", "pheV"],
            "island": ["PAI-IV", "PAI-V"],
            "mutations": ["gyrA-S83L", "gyrA-D87N", "parC-S80I"]
        },
        
        "risk_level": "HIGH",
        "outbreak_potential": "MODERATE to HIGH",
        "surveillance_priority": "EMERGING THREAT",
        "key_references": ["PMID: 35658504", "PMID: 39674367", "PMID: 39507943", "PMID: 31262826", "PMID: 36136628", "PMID: 30787145"]
    },

    "ST405": {
        "primary_name": "Global Multidrug-Resistant ExPEC Clone",
        "category": "Resistance-Associated",
        "sublineages": ["ST405-1", "ST405-2", "ST405-3"],
        "serotype": "O102:H6",
        "phylogroup": "D",
        "fumC": "fumC24",
        "fimH": "fimH30",
        "clermont_complex": "STc405",
        
        "pathotypes": ["UPEC", "SEPEC", "ExPEC"],
        "key_virulence_genes": ["fimH", "iutA", "fyuA", "kpsMII", "ompT", "traT", "iss"],
        
        "resistance_profile": {
            "esbl": ["blaCTX-M-15", "blaCTX-M-14", "blaCTX-M-27"],
            "carbapenemase": ["blaOXA-48", "blaNDM-1", "blaNDM-5"],
            "fluoroquinolone": ["gyrA mutations", "parC mutations", "qnr genes"],
            "aminoglycoside": ["aac(6')-Ib-cr", "aac(3)-II", "rmtB"],
            "other": ["sul1", "sul2", "tet(A)", "tet(B)", "dfrA variants"],
            "resistance_notes": "Extensive drug resistance, often pan-drug resistant"
        },
        
        "epidemiology": {
            "first_identified": 2005,
            "global_spread": "Established global MDR clone",
            "geographical_distribution": {
                "high_prevalence": ["Europe", "Middle East", "Asia"],
                "medium_prevalence": ["North America", "South America", "Africa"],
                "regional_hotspots": {
                    "Europe": "OXA-48 predominant, healthcare settings",
                    "Middle East": "NDM variants, both hospital and community",
                    "Asia": "Rapid dissemination, multiple resistance mechanisms",
                    "North America": "Increasing reports, primarily imported cases"
                }
            },
            "transmission": "Healthcare-associated predominant, some community spread",
            "reservoir": "Human gastrointestinal tract, environmental persistence"
        },
        
        "clinical_significance": {
            "primary_infections": ["Bloodstream infections", "UTI", "Intra-abdominal infections"],
            "treatment_challenges": "Limited therapeutic options, often requires carbapenems",
            "mortality": "High in bacteremia due to delayed effective therapy",
            "special_populations": "Immunocompromised, ICU patients, transplant recipients"
        },
        
        "genomic_features": {
            "plasmid": ["IncFII", "IncN", "IncL/M"],
            "prophage": ["yehV", "pheV"],
            "island": ["PAI-I", "PAI-II"],
            "resistance_islands": ["Tn3 transposons", "Class 1 integrons"]
        },
        
        "risk_level": "VERY HIGH",
        "outbreak_potential": "HIGH",
        "surveillance_priority": "CRITICAL",
        "key_references": ["PMID: 30671039", "PMID: 30671039", "DOI: 10.1016/j.cmi.2017.01.020", "PMID: 34965471", "PMID: 29321680"]
    }
}    
# =============================================================================
# SEROTYPE-SPECIFIC DATABASE
# =============================================================================

SEROTYPE_DATABASE = {
    "O157:H7": {
        "st": [11, 16, 20, 301],
        "primary_pathotype": "EHEC",
        "key_virulence": ["stx1", "stx2", "eae", "ehxA"],
        "h_us_risk": "HIGH",
        "geographical_distribution": {
            "high_incidence": ["North America", "Europe", "Japan", "Argentina"],
            "medium_incidence": ["Australia", "Chile", "South Africa"],
            "emerging_regions": ["Africa", "Asia"],
            "outbreak_patterns": "Global foodborne outbreaks, regional variations in prevalence"
        },
        "outbreak_association": "Ground beef, produce, water",
        "references": ["PMID: 15040260", "PMID: 15829201"]
    },
    
    "O104:H4": {
        "st": [678],
        "primary_pathotype": "EAEC/EHEC", 
        "key_virulence": ["stx2", "aggR", "aap", "aatA", "aaiC", "pic"], 
        "h_us_risk": "VERY HIGH",  
        "geographical_distribution": {
            "major_outbreak": "Germany 2011 with international cases",
            "sporadic_cases": ["Europe", "North America", "Asia"],
            "limited_spread": "Contained primarily to 2011 outbreak",
            "significance": "Demonstrated potential for hybrid pathotype emergence"
        },
        "outbreak_association": "2011 German sprout outbreak (fenugreek sprouts)",
        "references": [
            "Bielaszewska M, et al. Characterisation of the Escherichia coli strain associated with an outbreak of haemolytic uraemic syndrome in Germany, 2011. Lancet Infect Dis. 2011;11(9):671-6.",
            "Frank C, et al. Epidemic profile of Shiga-toxin-producing Escherichia coli O104:H4 outbreak in Germany. N Engl J Med. 2011;365(19):1771-80."
        ]
    },
    
    "O26:H11": {
        "st": [21, 29, 396],
        "primary_pathotype": "EHEC",
        "key_virulence": ["stx1", "stx2", "eae", "ehxA"], 
        "h_us_risk": "HIGH",
        "geographical_distribution": {
            "high_prevalence": ["Europe", "North America", "Japan"],
            "increasing_global": "Rising incidence worldwide",
            "regional_variants": {
                "Europe": "Common non-O157 STEC",
                "North America": "Increasing in beef and dairy",
                "Asia": "Emerging food safety concern"
            }
        },
        "outbreak_association": "Raw milk, beef products, increasing non-O157 STEC",
        "references": [
            "Karch H, et al. The enemy within us: lessons from the 2011 European Escherichia coli O104:H4 outbreak. EMBO Mol Med. 2012;4(9):841-8.",
            "Bettelheim KA. The non-O157 shiga-toxigenic (verocytotoxigenic) Escherichia coli; under-rated pathogens. Crit Rev Microbiol. 2007;33(1):67-87."
        ]
    },
    
    "O25b:H4": {
        "st": [131],
        "primary_pathotype": "UPEC/ExPEC",
        "key_virulence": ["fimH", "papA", "iutA", "fyuA"],
        "resistance": "CTX-M-15 common",
        "clinical": "Global MDR pandemic clone",
        "h_us_risk": "HIGH",  
        "geographical_distribution": {
            "pandemic_distribution": "Worldwide in healthcare and community",
            "high_prevalence": ["North America", "Europe", "Asia", "South America"],
            "dominant_variants": {
                "Global": "ST131-associated",
                "Europe": "CTX-M-15, CTX-M-27 variants",
                "Asia": "Rapid dissemination in hospitals"
            }
        },
        "outbreak_association": "Healthcare-associated infections, UTIs",
        "references": ["PMID: 32989036"]  
    },
    
    "O78:H4": {
        "st": [117, 350],
        "primary_pathotype": "APEC",
        "key_virulence": ["iss", "iroN", "iutA", "hlyF"],
        "host": "Poultry primary",
        "zoonotic": "Limited human transmission",
        "h_us_risk": "LOW",  
        "geographical_distribution": {
            "poultry_global": "Worldwide in poultry production",
            "high_prevalence": ["Europe", "North America", "Brazil", "China"],
            "economic_impact": "Global poultry industry losses",
            "zoonotic_cases": "Sporadic worldwide, usually occupational"
        },
        "outbreak_association": "Poultry infections",
        "references": ["PMID: 28049430"]  
    },
    
    "O111:H8": {
        "st": [16, 21],
        "primary_pathotype": "EHEC",
        "key_virulence": ["stx1", "stx2", "eae", "ehxA"],
        "h_us_risk": "HIGH",
        "geographical_distribution": {
            "global_distribution": "Worldwide, important non-O157 STEC",
            "high_incidence": ["North America", "Europe", "Australia"],
            "outbreak_association": "Multiple food vehicles globally",
            "public_health_concern": "Significant HUS risk"
        },
        "outbreak_association": "Deli meats, unpasteurized juice",
        "references": [
            "Brooks JT, et al. Non-O157 Shiga toxin-producing Escherichia coli infections in the United States, 1983-2002. J Infect Dis. 2005;192(8):1422-9.",
            "Paton AW, et al. Characterization of Saa, a novel autoagglutinating adhesin produced by locus of enterocyte effacement-negative Shiga-toxigenic Escherichia coli strains. Infect Immun. 2001;69(11):6999-7009."
        ]
    },
    
    "O145:H28": {
        "st": [342, 441],
        "primary_pathotype": "EHEC",
        "key_virulence": ["stx2", "eae", "ehxA"],
        "h_us_risk": "HIGH",
        "geographical_distribution": {
            "high_prevalence": ["North America", "Europe"],
            "emerging_regions": ["Asia", "South America"],
            "outbreak_patterns": "Produce-associated outbreaks increasing",
            "public_health": "Big Six non-O157 STEC of concern"
        },
        "outbreak_association": "Lettuce, recreational water",
        "references": [
            "Luna-Gierke RE, et al. Outbreaks of non-O157 Shiga toxin-producing Escherichia coli infection: USA. Epidemiol Infect. 2014;142(11):2270-80.",
            "Boisen N, et al. The role of EHEC O145:H28 in human disease. PLoS One. 2015;10(4):e0123823."
        ]
    },
    
    "O121:H19": {
        "st": [655, 656],
        "primary_pathotype": "EHEC",
        "key_virulence": ["stx2", "eae"],
        "h_us_risk": "MODERATE",
        "geographical_distribution": {
            "distribution": "Primarily North America and Europe",
            "outbreak_emergence": "Flour-associated outbreaks highlighted new vehicles",
            "increasing_detection": "Improved diagnostics revealing true prevalence"
        },
        "outbreak_association": "Flour, raw milk",
        "references": ["PMID: 29166238"]
    },
    
    "O55:H7": {
        "st": [342, 443],
        "primary_pathotype": "aEPEC",
        "key_virulence": ["eae", "tccP"],
        "h_us_risk": "LOW",
        "geographical_distribution": {
            "global_distribution": "Worldwide, often sporadic cases",
            "evolutionary_significance": "Considered EHEC O157 precursor",
            "clinical_importance": "Mainly pediatric diarrhea in developing regions"
        },
        "outbreak_association": "Sporadic cases, potential EHEC precursor",
        "references": ["PMID: 16573685", "PMID: 22328665"]
    },
    
    "O127:H6": {
        "st": [2, 6],
        "primary_pathotype": "tEPEC",
        "key_virulence": ["eae", "bfpA", "perA"],
        "h_us_risk": "MODERATE",
        "geographical_distribution": {
            "high_incidence": ["Developing countries", "Tropical regions"],
            "low_incidence": ["Developed countries"],
            "epidemiology": "Classic infant diarrhea in resource-limited settings",
            "global_health": "Significant pediatric burden in endemic areas"
        },
        "outbreak_association": "Infant diarrhea in developing countries",
        "references": [
            "Trabulsi LR, et al. Typical and atypical enteropathogenic Escherichia coli. Emerg Infect Dis. 2002;8(5):508-13.", " PMID: 12857773", "PMID: 35456842"]
    },
    
    "O78:H11": {
        "st": [10, 23, 38],
        "primary_pathotype": "ExPEC",
        "key_virulence": ["fimH", "kpsMTII", "iutA"],
        "h_us_risk": "MODERATE",
        "geographical_distribution": {
            "global_distribution": "Worldwide in human infections",
            "clinical_association": "UTIs and sepsis across all regions",
            "pathotype_diversity": "Can express various virulence combinations"
        },
        "outbreak_association": "Urinary tract infections, sepsis",
        "references": ["PMID: 20802035", "PMID: 20457787"]
    },
# =========================================================================
# NON-O157 EHEC SEROTYPES - COMPREHENSIVE ADDITIONS
# =========================================================================

    "O45:H2": {
        "st": [10, 301, 342, 443],
        "primary_pathotype": "EHEC",
        "secondary_pathotypes": ["STEC", "aEPEC"],
        "key_virulence": ["stx2", "stx2a", "stx2c", "eae", "ehxA", "tccP", "espA", "espB", "nleB", "nleE"],
        "shiga_toxin_profile": {
            "primary": "stx2a (high HUS risk)",
            "secondary": "stx2c",
            "stx1": "Rare",
            "toxin_notes": "Stx2a predominant, associated with severe disease outcomes"
        },
        "h_us_risk": "HIGH",
        "clinical_severity": "Moderate to severe, HUS cases reported",
        "geographical_distribution": {
            "high_incidence": ["North America", "Europe"],
            "medium_incidence": ["Australia", "South America"],
            "emerging_regions": ["Asia", "Africa"],
            "regional_patterns": {
                "United States": "Increasing detection in beef and produce",
                "Europe": "Sporadic cases and small outbreaks",
                "Australia": "Emerging non-O157 STEC",
                "Global_trend": "Increasing recognition with improved diagnostics"
            }
        },
        "outbreak_association": {
            "food_vehicles": ["Leafy greens", "Beef products", "Raw milk"],
            "environmental": ["Contaminated water", "Animal contact"],
            "notable_outbreaks": [
                "2007 US lettuce outbreak (multistate)",
                "2010 US raw milk outbreak",
                "Sporadic cases with animal contact"
            ],
            "outbreak_frequency": "Occasional outbreaks, primarily foodborne"
        },
        "animal_reservoirs": {
            "primary": "Cattle",
            "secondary": ["Sheep", "Goats", "Deer"],
            "shedding": "Intermittent fecal shedding in ruminants",
            "transmission_cycle": "Bovine reservoir → environmental contamination → human infection"
        },
        "diagnostic_considerations": {
            "screening": "Often missed in routine STEC screening",
            "detection_methods": ["PCR for stx2/eae", "Immunomagnetic separation", "Whole genome sequencing"],
            "culture": "Sorbitol-fermenting, requires specific serotyping",
            "public_health": "Reportable in many jurisdictions as non-O157 STEC"
        },
        "treatment_recommendations": {
            "supportive_care": "Hydration, monitoring for HUS",
            "antibiotic_controversy": "Avoid antibiotics due to potential increased toxin release",
            "antimotility_agents": "Generally contraindicated",
            "hospitalization": "Consider for high-risk patients or severe symptoms"
        },
        "prevention_control": {
            "food_safety": "Thorough cooking of beef, washing produce",
            "environmental": "Avoid unpasteurized dairy, farm water sources",
            "high_risk_groups": "Children <5, elderly, immunocompromised",
            "surveillance": "Enhanced non-O157 STEC monitoring"
        },
        "genomic_features": {
            "locus_of_enterocyte_effacement": "Present (eae-positive)",
            "O_island": ["OI-122", "OI-57"],
            "plasmid": ["pO113-like", "pO157-like variants"],
            "phage_integration": "Stx2-encoding phage locations variable"
        },
        "public_health_significance": "One of the 'Big Six' non-O157 STEC serogroups of concern",
        "references": [
            "PMID: 32455956",  # Non-O157 STEC epidemiology
            "PMID: 34197587",  # O45 outbreaks and virulence
            "PMID: 38310676",  # Non-O157 STEC clinical spectrum
            "PMID: 39370525"
            "https://beaconbio.org/en/event/?eventid=6c1cd171-b304-4692-b53b-1b1be60cb713&"
        ]
    },

    "O103:H2": {
        "st": [17, 20, 21, 301],
        "primary_pathotype": "EHEC",
        "secondary_pathotypes": ["STEC", "aEPEC"],
        "key_virulence": ["stx1", "stx1a", "eae", "ehxA", "espA", "espB", "nleA", "nleB", "tccP"],
        "shiga_toxin_profile": {
            "primary": "stx1a",
            "secondary": "stx2 (rare variants)",
            "stx2": "Uncommon but reported",
            "toxin_notes": "Stx1 predominant, generally lower HUS risk than Stx2 producers"
        },
        "h_us_risk": "MODERATE",
        "clinical_severity": "Mild to moderate, HUS less common than O157",
        "geographical_distribution": {
            "high_incidence": ["Europe", "North America"],
            "medium_incidence": ["Japan", "Australia", "South America"],
            "established_regions": {
                "Europe": "Well-established non-O157 STEC, multiple outbreaks",
                "North America": "Increasing detection and recognition",
                "Japan": "Significant public health concern",
                "Global": "Widely distributed in cattle populations"
            }
        },
        "outbreak_association": {
            "food_vehicles": ["Fermented meats", "Dairy products", "Produce", "Water"],
            "environmental": ["Recreational water", "Animal contact settings"],
            "notable_outbreaks": [
                "2000 Japan school outbreak (fermented squid)",
                "2004 Norway fermented sausage outbreak",
                "2010 US recreational water outbreak",
                "Multiple European dairy-associated outbreaks"
            ],
            "outbreak_frequency": "Regular outbreaks, diverse vehicles"
        },
        "animal_reservoirs": {
            "primary": "Cattle",
            "secondary": ["Sheep", "Goats"],
            "persistence": "Can persist in farm environments",
            "transmission_dynamics": "Direct animal contact and contaminated food/water"
        },
        "diagnostic_considerations": {
            "detection_challenges": "Often stx1-only, may be missed by some screening assays",
            "optimal_methods": ["Multiplex PCR for stx1/stx2/eae", "Culture with serotyping"],
            "public_health_reporting": "Reportable as non-O157 STEC in many regions",
            "surveillance_importance": "Routine monitoring in food and clinical samples"
        },
        "treatment_recommendations": {
            "clinical_management": "Supportive care with hydration monitoring",
            "antibiotic_considerations": "Generally avoided, though lower HUS risk than Stx2 producers",
            "high_risk_monitoring": "Monitor for HUS in children and elderly",
            "outbreak_response": "Rapid identification and source tracking"
        },
        "prevention_control": {
            "food_processing": "Adequate cooking/fermentation of meats",
            "dairy_safety": "Pasteurization of milk products",
            "environmental": "Water treatment, farm hygiene",
            "public_education": "Awareness of non-O157 STEC risks"
        },
        "genomic_features": {
            "locus_of_enterocyte_effacement": "Present (eae-positive)",
            "O_island": ["OI-122", "OI-71"],
            "plasmid": ["pO103", "pO157-like"],
            "phage_variation": "Stx1-encoding phage integration sites"
        },
        "public_health_significance": "Major non-O157 STEC serogroup, significant outbreak potential",
        "historical_context": "One of the first non-O157 STEC serogroups recognized as human pathogen",
        "references": [
            "PMID: 26895282",  # O103 outbreaks and characterization
            "PMID: 39840652",  # O103 virulence factors
            "PMID: 29884330",  # Non-O157 STEC clinical importance
            "PMID: 24077709",  # STEC clinical spectrum
            "PMID: 39840652 ",
            "PMID: 29884330"
            "Ethelberg S, et al. Outbreak of non-O157 Shiga toxin-producing Escherichia coli infection from consumption of beef sausage. Clin Infect Dis. 2009;48(8):e78-81."
        ]
    }
}

# =============================================================================
# CLERMONT PHYLOGROUP CHARACTERISTICS
# =============================================================================

PHYLOGROUP_DATABASE = {
    "A": {
        "characteristics": "Commensal, environmental",
        "pathogenic_potential": "Low",
        "common_st": [10, 48, 58, 155, 167, 617],
        "serotypes": ["O2:H6", "O8:H9", "O9:H4"],
        "virulence_genes": ["fimH", "gad", "lpfA"]
    },
    
    "B1": {
        "characteristics": "Commensal, animal-associated", 
        "pathogenic_potential": "Low to Moderate",
        "common_st": [56, 117, 155, 350, 398, 602],
        "serotypes": ["O8:H9", "O78:H4", "O86:H18"],
        "virulence_genes": ["fimH", "iss", "iroN", "iutA"]
    },
    
    "B2": {
        "characteristics": "Extraintestinal pathogenic",
        "pathogenic_potential": "High",
        "common_st": [73, 95, 127, 131, 144, 372, 404],
        "serotypes": ["O2:H1", "O6:H1", "O18:H7", "O25b:H4"],
        "virulence_genes": ["fimH", "papA", "hlyA", "cnf1", "iutA", "fyuA"]
    },
    
    "C": {
        "characteristics": "Animal-associated, emerging human pathogens",
        "pathogenic_potential": "Moderate to High",
        "common_st": [88, 92, 317, 345, 410, 448],
        "serotypes": ["O15:H1", "O77:H4", "O88:H25", "O141:H4"],
        "virulence_genes": ["fimH", "stx2", "eae", "iutA", "iss"],
        "notes": "Includes some STEC and ExPEC strains"
    },
    
    "D": {
        "characteristics": "Mixed commensal/pathogenic",
        "pathogenic_potential": "Moderate to High", 
        "common_st": [38, 59, 62, 69, 393, 405, 648],
        "serotypes": ["O15:H1", "O17/O77:H18", "O81:H9"],
        "virulence_genes": ["fimH", "aggR", "aap", "iutA", "fyuA"]
    },
    
    "E": {
        "characteristics": "EHEC/STEC associated",
        "pathogenic_potential": "High for diarrheal disease",
        "common_st": [11, 16, 17, 301],
        "serotypes": ["O157:H7", "O55:H7", "O145:H28"],
        "virulence_genes": ["stx1", "stx2", "eae", "ehxA"]
    },
    
    "F": {
        "characteristics": "Emerging MDR lineages",
        "pathogenic_potential": "Variable",
        "common_st": [648, 1193, 1415],
        "serotypes": ["O8:H9", "O15:H1", "O25:H4"],
        "virulence_genes": ["fimH", "iutA", "fyuA", "iss"]
    },
    
    "G": {
        "characteristics": "Recently identified, diverse",
        "pathogenic_potential": "Variable",
        "common_st": [117, 155, 453, 602, 1193],
        "serotypes": ["O2:H6", "O8:H9", "O78:H4"],
        "virulence_genes": ["fimH", "iss", "iroN", "iutA"],
        "notes": "Shares characteristics with B1, includes some APEC strains"
    }
}

# =============================================================================
# COMPREHENSIVE ESCHERICHIA COLI PATHOTYPE DATABASE 
# =============================================================================

PATHOTYPE_DATABASE = {
    # =========================================================================
    # DIARRHEAGENIC E. COLI PATHOTYPES
    # =========================================================================
    
    "EPEC": {
        "primary_name": "Enteropathogenic E. coli",
        "category": "Diarrheagenic",
        "subtypes": ["Typical EPEC (tEPEC)", "Atypical EPEC (aEPEC)"],
        "key_virulence_genes": ["eae", "bfpA", "bfpB", "perA", "perB", "perC"],
        "subtype_markers": {
            "tEPEC": ["eae", "bfpA"],
            "aEPEC": ["eae", "!bfpA"]
        },
        "serotypes": {
            "common": ["O55:H6", "O86:H34", "O111:H2", "O119:H6", "O125:H21", "O126:H27", "O127:H6", "O128:H2", "O142:H6"],
            "notes": "Diverse O and H serotypes, historically classic infant diarrhea serogroups"
        },
        "pathogenesis": {
            "mechanism": "Attaching and effacing (A/E) lesions, bundle-forming pilus (BFP) mediated aggregation",
            "attachment": "Intimate attachment via intimin (eae gene), formation of pedestals",
            "effector_proteins": "Tir, Esp proteins via Type III Secretion System"
        },
        "clinical_manifestations": {
            "primary": "Infantile diarrhea, watery diarrhea in children <2 years",
            "complications": "Persistent diarrhea, growth faltering",
            "epidemiology": "Endemic in developing countries, sporadic in developed countries"
        },
        "resistance_profile": {
            "common_resistance": ["Ampicillin", "Trimethoprim-sulfamethoxazole", "Tetracycline"],
            "genes": ["blaTEM", "sul1", "sul2", "tet(A)", "tet(B)"],
            "notes": "Increasing multidrug resistance in endemic areas"
        },
        "outbreak_potential": "MODERATE",
        "risk_level": "HIGH in pediatric populations",
        "key_references": ["PMID: 37388762", "PMID: 9390560", "PMID: 1398907", "PMID: 15845459"]
    },
    
    "ETEC": {
        "primary_name": "Enterotoxigenic E. coli",
        "category": "Diarrheagenic", 
        "subtypes": ["LT-only", "ST-only", "LT/ST"],
        "key_virulence_genes": ["eltA", "eltB", "estA", "estB", "cfa", "cs1", "cs2", "cs3", "csw"],
        "toxin_profiles": {
            "LT": ["eltA", "eltB"],
            "STa": ["estA"],
            "STb": ["estB"],
            "CFA": ["cfa", "cs1", "cs2", "cs3"]
        },
        "serotypes": {
            "common": ["O6:H16", "O8:H9", "O25:H42", "O78:H11", "O148:H28", "O159:H20"],
            "notes": "Diverse serotypes associated with different colonization factors"
        },
        "pathogenesis": {
            "mechanism": "Colonization of small intestine followed by toxin production",
            "toxin_action": "LT: ADP-ribosylation of G proteins; ST: Guanylate cyclase activation",
            "colonization": "Via colonization factor antigens (CFAs)"
        },
        "clinical_manifestations": {
            "primary": "Traveler's diarrhea, watery diarrhea in all age groups",
            "characteristics": "Watery diarrhea without blood or leukocytes",
            "duration": "Self-limiting, 3-5 days typically"
        },
        "resistance_profile": {
            "common_resistance": ["Doxycycline", "Trimethoprim-sulfamethoxazole", "Fluoroquinolones"],
            "genes": ["tet(B)", "sul1", "sul2", "qnr variants", "aac(6')-Ib-cr"],
            "notes": "High resistance rates in endemic areas, prophylaxis often ineffective"
        },
        "outbreak_potential": "HIGH",
        "risk_level": "MODERATE to HIGH in travelers and endemic areas",
        "key_references": ["PMID: 35358002", "PMID: 40895478", "PMID: 36371425"]
    },
    
    "EHEC": {
        "primary_name": "Enterohemorrhagic E. coli",
        "category": "Diarrheagenic",
        "subtypes": ["STEC", "EHEC O157", "Non-O157 EHEC"],
        "key_virulence_genes": ["stx1", "stx2", "eae", "ehxA", "tccP", "espA", "espB", "espD"],
        "shiga_toxin_variants": {
            "Stx1": ["stx1a", "stx1c", "stx1d"],
            "Stx2": ["stx2a", "stx2b", "stx2c", "stx2d", "stx2e", "stx2f", "stx2g"],
            "high_risk": ["stx2a", "stx2d"]
        },
        "serotypes": {
            "major": ["O157:H7", "O157:NM"],
            "non_O157": ["O26:H11", "O45:H2", "O103:H2", "O111:NM", "O121:H19", "O145:NM"],
            "notes": "O157:H7 most common in outbreaks, non-O157 increasing in prevalence"
        },
        "pathogenesis": {
            "mechanism": "Shiga toxin production with attaching and effacing lesions",
            "toxin_action": "Ribosome inactivation leading to cell death, endothelial damage",
            "complications": "Hemolytic uremic syndrome (HUS), thrombotic microangiopathy"
        },
        "clinical_manifestations": {
            "primary": "Bloody diarrhea (hemorrhagic colitis), severe abdominal pain",
            "complications": "HUS (hemolytic anemia, thrombocytopenia, renal failure)",
            "progression": "Diarrhea → bloody diarrhea → HUS (in 5-15% of cases)"
        },
        "resistance_profile": {
            "common_resistance": ["Ampicillin", "Streptomycin", "Sulfonamides", "Tetracycline"],
            "genes": ["blaTEM", "strA-strB", "sul1", "sul2", "tet(A)", "tet(B)"],
            "important_note": "ANTIBIOTIC USE MAY INCREASE RISK OF HUS"
        },
        "outbreak_potential": "VERY HIGH",
        "risk_level": "VERY HIGH due to HUS risk",
        "key_references": ["PMID: 30137794", "PMID: 26104359", "PMID: 35436562", "PMID: 16766659"]
    },
    
    "EIEC": {
        "primary_name": "Enteroinvasive E. coli", 
        "category": "Diarrheagenic",
        "subtypes": ["Classical EIEC", "EIEC-like"],
        "key_virulence_genes": ["ipaH", "ipaB", "ipaC", "ipaD", "virA", "ial", "icsA"],
        "invasion_genes": {
            "plasmid_encoded": ["ipaH", "ipaB", "ipaC", "ipaD", "virA", "ial"],
            "chromosomal": ["icsA", "virK", "msbB2"]
        },
        "serotypes": {
            "common": ["O28ac:H-", "O29:H-", "O112ac:H-", "O124:H30", "O136:H-", "O143:H-", "O144:H-", "O152:H-", "O159:H-", "O164:H-"],
            "notes": "Often non-motile (H-), related to Shigella species"
        },
        "pathogenesis": {
            "mechanism": "Invasion of colonic epithelial cells, intracellular replication",
            "invasion_process": "Trigger mechanism via Type III Secretion System",
            "inflammatory_response": "Massive inflammation, epithelial destruction"
        },
        "clinical_manifestations": {
            "primary": "Dysentery (bloody, mucoid diarrhea with fever)",
            "symptoms": "Fever, abdominal cramps, tenesmus, fecal leukocytes",
            "differentiation": "Clinically similar to shigellosis"
        },
        "resistance_profile": {
            "common_resistance": ["Ampicillin", "Trimethoprim-sulfamethoxazole", "Chloramphenicol"],
            "genes": ["blaTEM", "catA1", "sul1", "sul2", "dfrA variants"],
            "notes": "Multidrug resistance common, similar to Shigella patterns"
        },
        "outbreak_potential": "HIGH",
        "risk_level": "HIGH in institutional settings",
        "key_references": ["PMID: 29259590", "PMID: 40617833", "https://asm.org/articles/2019/march/whats-the-significance-of-enteroinvasive-e-coli"]
    },
    
    "EAEC": {
        "primary_name": "Enteroaggregative E. coli",
        "category": "Diarrheagenic",
        "subtypes": ["Typical EAEC", "Atypical EAEC"],
        "key_virulence_genes": ["aggR", "aap", "aatA", "aaiC", "astA", "pet", "pic", "sigA", "sepA"],
        "virulence_plasmids": {
            "pAA": ["aggR", "aap", "aatA"],
            "dispersin": ["aap"],
            "toxins": ["pet", "pic", "sigA", "sepA"]
        },
        "serotypes": {
            "common": ["O3:H2", "O15:H18", "O44:H18", "O77:H18", "O86:H-", "O92:H33", "O111:H21", "O126:H27"],
            "notes": "Highly diverse serotypes, not well-defined"
        },
        "pathogenesis": {
            "mechanism": "Stacked-brick adherence to intestinal mucosa, biofilm formation",
            "adherence": "Via aggregative adherence fimbriae (AAF)",
            "toxin_effects": "Mucosal damage, fluid secretion, inflammation"
        },
        "clinical_manifestations": {
            "primary": "Acute and persistent diarrhea in all age groups",
            "characteristics": "Watery diarrhea with mucus, without blood",
            "special_groups": "Important cause of persistent diarrhea in children and HIV patients"
        },
        "resistance_profile": {
            "common_resistance": ["Ampicillin", "Trimethoprim-sulfamethoxazole", "Tetracycline"],
            "genes": ["blaTEM", "sul1", "sul2", "tet(A)", "tet(B)", "dfrA variants"],
            "notes": "Extensive multidrug resistance in endemic areas"
        },
        "outbreak_potential": "MODERATE to HIGH",
        "risk_level": "MODERATE to HIGH in endemic areas",
        "key_references": ["PMID: 20300577", "PMID: 24982324", "PMID: 39355069"]
    },
    
    "DAEC": {
        "primary_name": "Diffusely Adherent E. coli",
        "category": "Diarrheagenic",
        "subtypes": ["Afa/Dr DAEC", "Non-Afa/Dr DAEC"],
        "key_virulence_genes": ["afaA", "draA", "daaD", "daaE", "f1845"],
        "adhesin_families": {
            "Afa/Dr": ["afaA", "draA", "daaD", "daaE"],
            "F1845": ["f1845"],
            "other": ["nfaA", "afaB", "afaC"]
        },
        "serotypes": {
            "common": ["O1:H7", "O2:H1", "O21:H-", "O75:H-", "O101:H-", "O115:H-"],
            "notes": "Diverse serotypes, often associated with urinary isolates"
        },
        "pathogenesis": {
            "mechanism": "Diffuse adherence to epithelial cells, promoting inflammation",
            "receptors": "Bind to CD55 (DAF) and CEACAMs",
            "cellular_effects": "Epithelial cell activation, cytokine production"
        },
        "clinical_manifestations": {
            "primary": "Diarrhea in children >2 years, urinary tract infections",
            "age_distribution": "More common in older children and adults",
            "controversy": "Pathogenic role in diarrhea still debated"
        },
        "resistance_profile": {
            "common_resistance": ["Ampicillin", "Tetracycline"],
            "genes": ["blaTEM", "tet(A)", "tet(B)"],
            "notes": "Variable resistance patterns"
        },
        "outbreak_potential": "LOW",
        "risk_level": "LOW to MODERATE",
        "key_references": ["PMID: 23374248", "PMID: 11825986", "PMID: 9453606", "PMID: 33163200"]
    },
    
    # =========================================================================
    # EXTRACTSINAL PATHOGENIC E. COLI (ExPEC)
    # =========================================================================
    
    "UPEC": {
        "primary_name": "Uropathogenic E. coli",
        "category": "Extraintestinal",
        "subtypes": ["Cystitis", "Pyelonephritis", "Complicated UTI"],
        "key_virulence_genes": ["fimH", "papA", "papC", "papG", "sfaA", "focA", "hlyA", "cnf1", "iutA", "fyuA", "iroN"],
        "virulence_factors": {
            "adhesins": ["fimH", "papA", "papC", "papG", "sfaA", "focA"],
            "toxins": ["hlyA", "cnf1"],
            "siderophores": ["iutA", "fyuA", "iroN"]
        },
        "phylogroups": ["B2", "D"],
        "sequence_types": {
            "high_risk": ["ST131", "ST73", "ST95", "ST69", "ST127"],
            "global": ["ST131 (CTX-M-15)", "ST73", "ST95", "ST69"]
        },
        "pathogenesis": {
            "mechanism": "Ascending colonization of urinary tract, biofilm formation",
            "adherence": "P fimbriae (pyelonephritis), type 1 fimbriae (cystitis)",
            "invasion": "Intracellular bacterial communities, quiescent reservoirs"
        },
        "clinical_manifestations": {
            "spectrum": "Cystitis → pyelonephritis → bacteremia",
            "complications": "Renal scarring, urosepsis, recurrent infections"
        },
        "resistance_profile": {
            "common_resistance": ["Ampicillin", "Trimethoprim-sulfamethoxazole", "Fluoroquinolones"],
            "esbl_genes": ["blaCTX-M-15", "blaCTX-M-14", "blaCTX-M-27"],
            "high_risk_clones": ["ST131 (global MDR clone)", "ST410", "ST648"]
        },
        "outbreak_potential": "MODERATE in healthcare settings",
        "risk_level": "HIGH due to antimicrobial resistance",
        "key_references": ["PMID: 37764013", "PMID: 26443763", "PMID: 26443763", "PMID: 37445714"]
    },
    
    "NMEC": {
        "primary_name": "Neonatal Meningitis E. coli",
        "category": "Extraintestinal", 
        "subtypes": ["Neonatal", "Adult meningitis"],
        "key_virulence_genes": ["ibeA", "ibeB", "ibeC", "cnf1", "hlyA", "fimH", "sfaA", "kpsMTII", "neuC"],
        "specific_factors": {
            "invasion": ["ibeA", "ibeB", "ibeC"],
            "capsule": ["kpsMTII", "neuC"],
            "other": ["aslA", "ompA", "traJ"]
        },
        "serotypes": {
            "primary": ["O18:K1:H7", "O1:K1:H7", "O7:K1:H-", "O16:K1:H-", "O83:K1:H-"],
            "notes": "K1 capsule is major virulence determinant"
        },
        "phylogroups": ["B2", "D"],
        "pathogenesis": {
            "mechanism": "Gastrointestinal colonization → bacteremia → meningeal invasion",
            "key_steps": ["Maternal transmission", "Bacteremia", "BBB crossing", "Meningeal inflammation"]
        },
        "clinical_manifestations": {
            "primary": "Neonatal meningitis (first month of life)",
            "mortality": "High mortality (15-40%) despite treatment",
            "sequelae": "Neurological deficits in 30-50% of survivors"
        },
        "resistance_profile": {
            "common_resistance": ["Ampicillin", "Gentamicin"],
            "genes": ["blaTEM", "aac(3)-II", "aac(6')-Ib"],
            "emerging": "ESBL-producing NMEC increasing"
        },
        "outbreak_potential": "LOW",
        "risk_level": "VERY HIGH due to mortality and sequelae",
        "key_references": ["PMID: 22706051", "PMID: 30335297", "PMID: 40430805", "PMID: 26467858"]
    },
    
    "SEPEC": {
        "primary_name": "Sepsis-associated E. coli",
        "category": "Extraintestinal",
        "subtypes": ["Primary bacteremia", "Secondary bacteremia"],
        "key_virulence_genes": ["fimH", "papA", "hlyA", "cnf1", "iutA", "fyuA", "iroN", "kpsMTII"],
        "virulence_factors": {
            "shared_upec": ["fimH", "papA", "hlyA", "cnf1"],
            "siderophores": ["iutA", "fyuA", "iroN"],
            "capsule": ["kpsMTII"]
        },
        "phylogroups": ["B2", "D"],
        "sequence_types": {
            "common": ["ST131", "ST73", "ST95", "ST69", "ST127", "ST393", "ST405"],
            "notes": "Overlap with UPEC high-risk clones"
        },
        "pathogenesis": {
            "mechanism": "Translocation from gut or ascending UTI → bacteremia → sepsis",
            "risk_factors": ["UTI", "Abdominal source", "Immunocompromised", "Healthcare exposure"]
        },
        "clinical_manifestations": {
            "spectrum": "Bacteremia → severe sepsis → septic shock",
            "sources": "Urinary (most common), abdominal, unknown primary"
        },
        "resistance_profile": {
            "common_resistance": ["Ampicillin", "Fluoroquinolones", "Third-generation cephalosporins"],
            "carbapenemase": ["blaKPC", "blaNDM", "blaOXA-48"],
            "high_risk": "ESBL and carbapenemase producers increasing"
        },
        "outbreak_potential": "MODERATE in healthcare settings",
        "risk_level": "VERY HIGH due to mortality",
        "key_references": ["PMID: 22488222", "PMID: 35227247", "PMID: 26963151", "PMID: 38066783", "PMID: 38547882", "PMID: 22488222"]
    },
    
    # =========================================================================
    # OTHER PATHOTYPES
    # =========================================================================
    
    "EAHEC": {
        "primary_name": "Enteroaggregative-Hemorrhagic E. coli",
        "category": "Hybrid",
        "subtypes": ["O104:H4 German outbreak"],
        "key_virulence_genes": ["stx2", "aggR", "aap", "aatA", "terC", "aaiC", "pic"],
        "hybrid_features": {
            "EAEC_background": ["aggR", "aap", "aatA", "aaiC"],
            "EHEC_toxin": ["stx2"],
            "additional": ["terC", "pic"]
        },
        "serotypes": {
            "outbreak": ["O104:H4"],
            "notes": "Rare hybrid combining EAEC virulence with Shiga toxin"
        },
        "pathogenesis": {
            "mechanism": "EAEC adherence and biofilm + Shiga toxin production",
            "unique_aspects": "Enhanced adherence with toxin production",
            "outbreak": "2011 German outbreak with high HUS incidence"
        },
        "clinical_manifestations": {
            "primary": "Watery diarrhea progressing to bloody diarrhea",
            "complications": "Very high HUS rate (~20%), neurological complications",
            "demographics": "Primarily adults, unusual for STEC"
        },
        "resistance_profile": {
            "resistance": ["Ampicillin", "Tetracycline", "Trimethoprim-sulfamethoxazole", "Streptomycin"],
            "genes": ["blaTEM", "tet(A)", "tet(B)", "sul1", "sul2", "strA-strB"],
            "susceptible": "Carbapenems, fluoroquinolones effective"
        },
        "outbreak_potential": "VERY HIGH",
        "risk_level": "VERY HIGH due to high HUS risk",
        "key_references": ["PMID: 22895033", "PMID: 32874471", "PMID: 31262610", "PMID: 33939753", "PMID: 25887577"]
    },
    
    "ATEC": {
        "primary_name": "Avian Pathogenic E. coli",
        "category": "Animal",
        "subtypes": ["Poultry", "Wild birds"],
        "key_virulence_genes": ["iss", "iroN", "iutA", "hlyF", "ompT", "vat", "ibeA", "tsh"],
        "avian_factors": {
            "serum_resistance": ["iss"],
            "siderophores": ["iroN", "iutA"],
            "toxins": ["hlyF", "vat"],
            "other": ["ompT", "tsh"]
        },
        "serotypes": {
            "common": ["O1:K1", "O2:K1", "O18:K1", "O78:K80"],
            "notes": "Limited serotype diversity in poultry"
        },
        "pathogenesis": {
            "mechanism": "Respiratory colonization → septicemia → systemic infection",
            "diseases": "Colibacillosis, airsacculitis, pericarditis, perihepatitis"
        },
        "clinical_significance": "Major economic impact in poultry industry, zoonotic potential low",
        "resistance_profile": {
            "common_resistance": ["Tetracycline", "Sulfonamides", "Ampicillin", "Fluoroquinolones"],
            "genes": ["tet(A)", "tet(B)", "sul1", "sul2", "blaTEM", "qnr variants"],
            "notes": "Extensive antimicrobial use in poultry drives resistance"
        },
        "outbreak_potential": "HIGH in poultry operations",
        "risk_level": "LOW for humans, HIGH for poultry",
        "key_references": ["PMID: 33921518", "PMID: 39388979", "PMID: 39334984"]
    },


# =========================================================================
# SPECIALIZED EXTRACTSINAL & MUCOSAL PATHOTYPES - COMPREHENSIVE ADDITIONS
# =========================================================================

    "MAEC": {
        "primary_name": "Meningitis-Associated E. coli",
        "category": "Extraintestinal",
        "subtypes": ["Neonatal MAEC", "Adult MAEC", "Healthcare-associated MAEC"],
        "key_virulence_genes": ["ibeA", "ibeB", "ibeC", "cnf1", "hlyA", "kpsMTII", "neuC", "fimH", "sfaA", "ompA", "aslA", "traJ"],
        "specific_factors": {
            "blood_brain_barrier": ["ibeA", "ibeB", "ibeC", "ompA"],
            "cytotoxins": ["cnf1", "hlyA"],
            "capsule": ["kpsMTII", "neuC"],
            "adhesins": ["fimH", "sfaA"],
            "other": ["aslA", "traJ"]
        },
        "serotypes": {
            "common": ["O18:K1:H7", "O1:K1:H7", "O7:K1:H-", "O16:K1:H-", "O83:K1:H-", "O45:K1:H7"],
            "notes": "K1 capsule is critical virulence determinant, enables immune evasion"
        },
        "phylogroups": ["B2", "D"],
        "sequence_types": {
            "high_risk": ["ST95", "ST73", "ST131", "ST88", "ST127"],
            "neonatal_focused": ["ST95", "ST88"],
            "notes": "Overlap with UPEC/SEPEC high-risk clones"
        },
        "pathogenesis": {
            "mechanism": "Gastrointestinal colonization → bacteremia → meningeal invasion → CNS inflammation",
            "key_steps": [
                "Maternal vaginal colonization or environmental acquisition",
                "Gastrointestinal colonization in neonate",
                "Bacteremia via translocation or mucosal breach", 
                "Blood-brain barrier crossing via specific invasion mechanisms",
                "Meningeal colonization and inflammatory response"
            ],
            "blood_brain_barrier": "IbeA-mediated invasion of brain microvascular endothelial cells",
            "immune_evasion": "K1 capsule prevents complement deposition and phagocytosis"
        },
        "clinical_manifestations": {
            "primary": "Neonatal meningitis (first 28 days of life), adult meningitis in immunocompromised",
            "neonatal_presentation": "Fever, lethargy, poor feeding, seizures, bulging fontanelle",
            "adult_presentation": "Fever, headache, neck stiffness, altered mental status",
            "complications": "Brain abscess, ventriculitis, hydrocephalus, cerebral infarction",
            "mortality": "15-40% despite appropriate antibiotic therapy",
            "sequelae": "Neurological deficits in 30-50% of survivors (hearing loss, developmental delay, seizures)"
        },
        "diagnostic_considerations": {
            "csf_findings": "Elevated protein, decreased glucose, neutrophilic pleocytosis",
            "culture": "Blood and CSF culture positive",
            "pcr_detection": "K1 capsule genes, ibeA for rapid diagnosis",
            "imaging": "MRI may show meningeal enhancement, complications"
        },
        "resistance_profile": {
            "common_resistance": ["Ampicillin", "Gentamicin"],
            "genes": ["blaTEM", "aac(3)-II", "aac(6')-Ib"],
            "emerging_resistance": ["Third-generation cephalosporins", "Carbapenems"],
            "esbl_concern": "Increasing reports of ESBL-producing MAEC",
            "treatment_implications": "Ampicillin resistance common, requires third-gen cephalosporin + aminoglycoside"
        },
        "epidemiology": {
            "incidence": "Leading cause of neonatal meningitis in developed countries",
            "neonatal_risk_factors": ["Prematurity", "Low birth weight", "Maternal colonization", "Healthcare exposure"],
            "adult_risk_factors": ["Neurosurgical procedures", "CSF shunts", "Immunocompromised", "Head trauma"],
            "geographical_distribution": "Global distribution, higher incidence in settings with limited prenatal care",
            "transmission": "Vertical (mother to neonate), nosocomial, community-acquired"
        },
        "treatment_recommendations": {
            "empiric_therapy": "Cefotaxime or ceftriaxone + aminoglycoside",
            "duration": "Minimum 14-21 days, longer for complications",
            "adjunctive_therapy": "Dexamethasone controversial in neonates, considered in adults",
            "monitoring": "Repeat CSF analysis, hearing evaluation, neurodevelopmental follow-up"
        },
        "prevention": {
            "maternal_screening": "Not routinely recommended",
            "infection_control": "Hand hygiene, environmental cleaning in NICUs",
            "vaccine_development": "K1 capsule vaccines in preclinical development"
        },
        "outbreak_potential": "LOW to MODERATE in neonatal ICUs",
        "risk_level": "VERY HIGH due to mortality and neurological sequelae",
        "key_references": ["PMID: 22706051", "PMID: 40430805", "PMID: 26467858", "PMID: 19307211", "PMID: 26865687"]
    },

    "AIEC": {
        "primary_name": "Adherent-Invasive E. coli",
        "category": "Mucosal",
        "subtypes": ["Crohn's disease-associated", "Ulcerative colitis-associated", "Sporadic"],
        "key_virulence_genes": ["fimH", "lpfA", "ompC", "ibeA", "gspE", "yfgL", "pduC", "chuA", "ybtS"],
        "virulence_factors": {
            "adherence": ["fimH", "lpfA"],
            "invasion": ["ompC", "ibeA"],
            "secretion": ["gspE"],
            "metabolic_adaptation": ["pduC", "chuA", "ybtS"],
            "outer_membrane": ["yfgL"]
        },
        "in_vitro_characteristics": {
            "adherence": "Adherence to intestinal epithelial cells",
            "invasion": "Ability to invade intestinal epithelial cells",
            "intracellular_survival": "Survival and replication within macrophages",
            "biofilm_formation": "Enhanced biofilm formation on intestinal mucosa",
            "inflammatory_response": "Induction of TNF-α and other pro-inflammatory cytokines"
        },
        "serotypes": {
            "diverse": ["O6", "O22", "O83", "O150", "O2", "O25"],
            "notes": "No specific serotype association, highly diverse"
        },
        "phylogroups": ["B2", "D", "B1"],
        "sequence_types": {
            "reported": ["ST59", "ST73", "ST95", "ST131"],
            "notes": "No specific ST association, can emerge from various genetic backgrounds"
        },
        "pathogenesis": {
            "mechanism": "Mucosal adherence → epithelial invasion → macrophage survival → chronic inflammation",
            "key_steps": [
                "Colonization of ileal mucosa via FimH and LpfA",
                "Invasion of intestinal epithelial cells via OmpC and IbeA",
                "Survival and replication within macrophages (avoiding killing)",
                "Induction of TNF-α and pro-inflammatory cytokine production",
                "Disruption of autophagy processes in host cells",
                "Promotion of Th1/Th17 immune responses"
            ],
            "inflammatory_cascade": "NF-κB activation, IL-8 production, neutrophil recruitment",
            "autophagy_interference": "Impairment of autophagy in macrophages leading to bacterial persistence"
        },
        "clinical_association": {
            "primary": "Crohn's disease (particularly ileal disease)",
            "disease_correlation": "Higher prevalence in Crohn's patients vs. controls (30-40% vs. 5-10%)",
            "disease_activity": "Correlation with disease activity and postoperative recurrence",
            "therapeutic_implications": "Potential target for anti-adhesive therapies, bacteriophages"
        },
        "diagnostic_detection": {
            "culture_methods": "Mucosal biopsy culture with specific adherence/invasion assays",
            "molecular_detection": "PCR for virulence genes, whole genome sequencing",
            "functional_assays": [
                "Cell adherence assays (HEp-2, Caco-2)",
                "Invasion assays with gentamicin protection", 
                "Macrophage survival assays",
                "Biofilm formation assays"
            ],
            "challenges": "No standardized diagnostic criteria, research tool primarily"
        },
        "epidemiology": {
            "prevalence_crohns": "Detected in 30-40% of Crohn's disease patients",
            "geographical_distribution": "Global distribution, similar prevalence across regions",
            "tissue_localization": "Primarily ileal mucosa in Crohn's disease",
            "temporal_relationship": "May precede clinical disease onset"
        },
        "research_significance": {
            "pathogenesis_insights": "Model for microbiome-gut-brain axis in IBD",
            "therapeutic_targets": "FimH inhibitors, bacteriophage therapy, vaccines in development",
            "diagnostic_potential": "Potential biomarker for Crohn's disease subtypes",
            "microbiome_interactions": "Interacts with other gut microbiota members"
        },
        "treatment_considerations": {
            "current_antibiotics": "Variable susceptibility, not routinely targeted",
            "emerging_therapies": [
                "FimH antagonists (mannosides)",
                "Bacteriophage therapy",
                "Anti-adhesion vaccines",
                "Microbiome modulation"
            ],
            "clinical_trials": "Early phase trials for FimH inhibitors in IBD"
        },
        "outbreak_potential": "LOW (not typically outbreak-associated)",
        "risk_level": "MODERATE (chronic disease association)",
        "key_references": ["PMID: 20300577", "PMID: 23382754", "PMID: 29141957", "PMID: 36883813", "PMID: 25338542"]
    }
}

# =============================================================================
# SPECIALIZED PATHOTYPE PROFILES
# =============================================================================

SPECIALIZED_PROFILES = {
    "SHIGA_TOXIN_POSITIVE": {
        "EHEC_O157": {
            "serotype": "O157:H7",
            "stx_profile": ["stx1", "stx2", "stx2a", "stx2c"],
            "virulence": ["eae", "ehxA", "tccP"],
            "clinical_risk": "HIGH HUS risk",
            "outbreak_association": "Foodborne outbreaks, ground beef"
        },
        "EHEC_O104": {
            "serotype": "O104:H4", 
            "stx_profile": ["stx2a"],
            "virulence": ["aggR", "aap", "aatA"],
            "clinical_risk": "VERY HIGH HUS risk",
            "outbreak_association": "2011 German sprout outbreak"
        },
        "EHEC_O26": {
            "serotype": "O26:H11",
            "stx_profile": ["stx1", "stx2a"],
            "virulence": ["eae", "ehxA"],
            "clinical_risk": "HIGH HUS risk",
            "outbreak_association": "Increasing in prevalence"
        }
    },
    
    "ESBL_PRODUCERS": {
        "ST131_CTX-M-15": {
            "pathotype": "UPEC/SEPEC",
            "st": 131,
            "esbl": ["blaCTX-M-15"],
            "resistance": ["Fluoroquinolones", "Aminoglycosides", "Trimethoprim-sulfamethoxazole"],
            "global_spread": "Worldwide pandemic clone"
        },
        "ST410_CTX-M-15": {
            "pathotype": "UPEC/SEPEC", 
            "st": 410,
            "esbl": ["blaCTX-M-15"],
            "resistance": ["Fluoroquinolones", "Aminoglycosides"],
            "emergence": "Rapidly emerging globally"
        },
        "ST648_CTX-M": {
            "pathotype": "Various",
            "st": 648,
            "esbl": ["blaCTX-M-14", "blaCTX-M-15", "blaCTX-M-27"],
            "resistance": ["Multiple drug classes"],
            "notes": "Zoonotic potential, broad host range"
        }
    }
}


# =========================================================================
# CARBAPENEMASE-PRODUCING E. COLI PROFILES - COMPREHENSIVE ADDITION
# =========================================================================

CARBAPENEMASE_PRODUCERS = {
    "KPC_PRODUCERS": {
        "pathotype": "Various (UPEC/SEPEC/ExPEC predominant)",
        "st": [258, 131, 405, 410, 648, 101, 38, 69],
        "carbapenemase": ["blaKPC-2", "blaKPC-3", "blaKPC-4", "blaKPC-8", "blaKPC-14"],
        "enzyme_class": "Ambler Class A serine β-lactamase",
        "hydrolysis_spectrum": ["Carbapenems", "Penicillins", "Cephalosporins", "Aztreonam"],
        "inhibitor_profile": {
            "inhibited_by": ["Avibactam", "Relebactam", "Vaborbactam"],
            "resistant_to": ["Clavulanate", "Tazobactam", "Sulbactam"]
        },
        "genetic_context": {
            "gene_location": ["Plasmid-borne (IncFII, IncN)", "Chromosomal (rare)"],
            "mobile_element": "Tn4401 transposon variants",
            "co-resistance": ["ESBL genes", "Aminoglycoside resistance", "Fluoroquinolone resistance"]
        },
        "geographical_distribution": {
            "endemic_regions": ["United States (Northeast)", "Greece", "Israel", "Italy", "China", "Brazil"],
            "emerging_regions": ["Western Europe", "Latin America", "Southeast Asia"],
            "global_spread": "Originally Klebsiella-dominated, now spreading to E. coli",
            "hotspots": {
                "New York": "Early epicenter, healthcare networks",
                "Greece": "High prevalence in hospitals",
                "China": "Rapid dissemination in multiple provinces",
                "Brazil": "Increasing in urban hospitals"
            }
        },
        "clinical_significance": {
            "infections": ["Bloodstream infections", "UTI", "Pneumonia", "Intra-abdominal infections"],
            "mortality": "30-50% in bacteremia",
            "treatment_challenges": "Limited therapeutic options, often requires combination therapy",
            "risk_factors": ["Healthcare exposure", "ICU stay", "Device-associated", "Prior antibiotic use"]
        },
        "detection_methods": {
            "phenotypic": ["Modified Hodge test", "Carba NP", "mCIM"],
            "molecular": ["PCR for blaKPC", "Whole genome sequencing", "Multiplex arrays"],
            "automated_systems": ["VITEK2", "BD Phoenix", "MicroScan"],
            "challenges": "Variable expression, false negatives in some tests"
        },
        "treatment_options": {
            "first_line": ["Ceftazidime-avibactam", "Meropenem-vaborbactam", "Imipenem-relebactam"],
            "alternative": ["Plazomicin + meropenem", "Tigecycline", "Eravacycline", "Cefiderocol"],
            "combination_therapy": ["Aztreonam + ceftazidime-avibactam", "Multiple agent combinations"],
            "important_notes": "Therapy guided by MIC testing and local epidemiology"
        },
        "infection_control": {
            "isolation": "Contact precautions, single rooms",
            "screening": "Active surveillance in high-risk units",
            "environmental": "Enhanced cleaning, dedicated equipment",
            "outbreak_management": "Cohorting, antibiotic stewardship"
        },
        "references": [
            "PMID: 37440483",  # KPC epidemiology and treatment
            "PMID: 24354657",  # Global carbapenemase spread
            "PMID: 34202216",  # Plasmid-mediated resistance
            "PMID: 37887195"   # Carbapenemase detection and management
        ]
    },

    "NDM_PRODUCERS": {
        "pathotype": "Various (Broad host range)",
        "st": [410, 648, 167, 101, 131, 38, 405, 617],
        "carbapenemase": ["blaNDM-1", "blaNDM-5", "blaNDM-7", "blaNDM-9", "blaNDM-12"],
        "enzyme_class": "Ambler Class B metallo-β-lactamase",
        "hydrolysis_spectrum": ["Carbapenems", "Penicillins", "Cephalosporins"],
        "inhibitor_profile": {
            "inhibited_by": ["EDTA", "Dipicolinic acid"],
            "resistant_to": ["All β-lactamase inhibitors", "Aztreonam (unless combined)"]
        },
        "genetic_context": {
            "gene_location": ["Plasmid-borne (IncX3, IncFII, IncL/M)", "Chromosomal (rare)"],
            "mobile_element": "ISAba125, Tn125",
            "co-resistance": ["16S rRNA methylases", "Aminoglycoside resistance", "Fluoroquinolone resistance", "Colistin resistance"]
        },
        "geographical_distribution": {
            "endemic_regions": ["Indian subcontinent", "Balkans", "Middle East", "Southeast Asia"],
            "global_spread": "Rapid worldwide dissemination via travel and migration",
            "regional_variants": {
                "Indian subcontinent": "NDM-1, NDM-5 dominant, community and hospital",
                "Balkans": "NDM-1, healthcare-associated",
                "Middle East": "NDM-1, NDM-7, travel-associated",
                "China": "NDM-5 increasing, diverse plasmids",
                "United Kingdom": "Travel-imported cases, local transmission"
            }
        },
        "clinical_significance": {
            "infections": ["UTI", "Bloodstream infections", "Intra-abdominal infections", "Wound infections"],
            "community_acquisition": "Significant community prevalence in endemic areas",
            "mortality": "35-60% in bacteremia",
            "one_health_aspect": "Animal and environmental reservoirs, zoonotic potential"
        },
        "detection_methods": {
            "phenotypic": ["Carba NP", "mCIM", "eCIM for differentiation"],
            "molecular": ["PCR for blaNDM", "Multiplex assays", "Whole genome sequencing"],
            "automated_systems": ["VITEK2", "BD Phoenix", "MicroScan"],
            "challenges": "Co-production with other β-lactamases common"
        },
        "treatment_options": {
            "first_line": ["Aztreonam + ceftazidime-avibactam", "Cefiderocol"],
            "alternative": ["Plazomicin", "Tigecycline", "Eravacycline", "Fosfomycin"],
            "combination_therapy": ["Multiple agent regimens based on susceptibility"],
            "important_notes": "Aztreonam alone inactive due to co-production of other β-lactamases"
        },
        "infection_control": {
            "isolation": "Contact precautions, single rooms",
            "screening": "Travel history-based screening, high-risk units",
            "environmental": "Enhanced cleaning, water source monitoring",
            "travel_medicine": "Education for travelers to endemic areas"
        },
        "references": [
            "PMID: 24982081",  # NDM epidemiology and genetics
            "PMID: 37303773",  # Global NDM spread
            "PMID: 30745391",  # Plasmid dynamics
            "PMID: 38661186"   # NDM treatment challenges
        ]
    },

    "OXA-48-LIKE_PRODUCERS": {
        "pathotype": "Various (UPEC/SEPEC predominant)",
        "st": [38, 131, 410, 648, 69, 405, 101],
        "carbapenemase": ["blaOXA-48", "blaOXA-181", "blaOXA-232", "blaOXA-204", "blaOXA-244"],
        "enzyme_class": "Ambler Class D serine β-lactamase",
        "hydrolysis_spectrum": ["Carbapenems (weak)", "Penicillins"],
        "inhibitor_profile": {
            "inhibited_by": ["Avibactam", "Relebactam"],
            "resistant_to": ["Clavulanate", "Tazobactam"]
        },
        "genetic_context": {
            "gene_location": ["Plasmid-borne (IncL/M, IncX3)"],
            "mobile_element": "Tn1999 transposon variants",
            "co-resistance": ["ESBL genes", "Aminoglycoside resistance"]
        },
        "geographical_distribution": {
            "endemic_regions": ["Turkey", "North Africa", "Middle East", "Mediterranean"],
            "global_spread": "Increasing in Europe, sporadic worldwide",
            "regional_patterns": {
                "Turkey": "OXA-48 endemic, healthcare and community",
                "North Africa": "OXA-48, OXA-181 increasing",
                "Middle East": "OXA-48 variants, travel-associated",
                "Europe": "Increasing importation and local transmission"
            }
        },
        "clinical_significance": {
            "infections": ["UTI", "Bloodstream infections", "Intra-abdominal infections"],
            "detection_challenge": "Weak carbapenemase activity may be missed",
            "mortality": "25-40% in bacteremia",
            "epidemiology": "Often healthcare-associated with some community transmission"
        },
        "detection_methods": {
            "phenotypic": ["Carba NP", "mCIM", "Modified Hodge test"],
            "molecular": ["PCR for blaOXA-48-like", "Whole genome sequencing"],
            "automated_systems": ["VITEK2", "BD Phoenix"],
            "challenges": "Weak hydrolysis may give borderline results"
        },
        "treatment_options": {
            "first_line": ["Ceftazidime-avibactam", "Meropenem-vaborbactam"],
            "alternative": ["Plazomicin", "Tigecycline", "Eravacycline", "Cefiderocol"],
            "combination_therapy": ["Based on susceptibility testing"],
            "important_notes": "May appear susceptible to carbapenems but treatment failures reported"
        },
        "infection_control": {
            "isolation": "Contact precautions",
            "screening": "Active surveillance in endemic regions",
            "environmental": "Standard enhanced cleaning",
            "regional_cooperation": "Cross-border surveillance important"
        },
        "references": [
            "PMID: 36103096",  # OXA-48 epidemiology
            "PMID: 25340464",  # Global carbapenemase spread
            "PMID: 34202216",  # Plasmid-mediated resistance
            "PMID: 40885895"   # Detection and management
        ]
    },

    "VIM_IMP_PRODUCERS": {
        "pathotype": "Various (Less common in E. coli)",
        "st": [131, 410, 648, 38, 69],
        "carbapenemase": ["blaVIM-1", "blaVIM-2", "blaIMP-1", "blaIMP-8"],
        "enzyme_class": "Ambler Class B metallo-β-lactamase",
        "hydrolysis_spectrum": ["Carbapenems", "Penicillins", "Cephalosporins"],
        "inhibitor_profile": {
            "inhibited_by": ["EDTA", "Dipicolinic acid"],
            "resistant_to": ["All β-lactamase inhibitors"]
        },
        "genetic_context": {
            "gene_location": ["Plasmid-borne", "Integron-associated"],
            "mobile_element": ["Class 1 integrons", "Tn21-like transposons"],
            "co-resistance": ["Aminoglycoside resistance", "Fluoroquinolone resistance"]
        },
        "geographical_distribution": {
            "endemic_regions": ["Greece", "Italy", "Taiwan", "Japan"],
            "sporadic_distribution": "Rare but globally dispersed",
            "regional_variants": {
                "Greece": "VIM-1 endemic in certain regions",
                "Italy": "VIM variants in hospital settings",
                "East Asia": "IMP variants more common"
            }
        },
        "clinical_significance": {
            "infections": ["UTI", "Bloodstream infections"],
            "prevalence": "Rare in E. coli compared to other Enterobacteriaceae",
            "mortality": "Similar to other carbapenemase producers when they occur"
        },
        "references": [
            "PMID: 35847092",  # Metallo-β-lactamase epidemiology
            "PMID: 38638826",  # Resistance mechanisms
            "PMID: 30003866"   # Detection and management
        ]
    }
}

# =============================================================================
# COMPREHENSIVE REFERENCE SECTION
# =============================================================================

COMPREHENSIVE_REFERENCES = {
    # =========================================================================
    # PART 1: PMID REFERENCES
    # =========================================================================
    "PUBMED_REFERENCES": {
        "EPIDEMIOLOGY_REVIEWS": [
            "PMID: 30266330",  # Global burden of E. coli diarrhea
            "PMID: 30266330",  # ETEC epidemiology and vaccine development  
            "PMID: 40895478",  # ST131 global spread
            "PMID: 22491693",  # EIEC and Shigella evolution
            "PMID: 14562958",  # ExPEC pathogenesis and epidemiology
            "PMID: 29667573",  # EHEC outbreaks and pathogenesis
            "PMID: 22919681",  # EPEC molecular pathogenesis
            "PMID: 35187883",  # EAEC epidemiology and pathogenesis
            "PMID: 24982324",  # Global antimicrobial resistance in E. coli
            "PMID: 18455741",  # EPEC disease burden update
            "PMID: 41196327",  # ETEC vaccine progress
            "PMID: 23382754",  # MDR E. coli in healthcare settings
            "PMID: 36668830"   # STEC diagnosis and management
        ],
        
        "PATHOGENESIS_MECHANISMS": [
            "PMID: 22289607",  # EAHEC O104:H4 outbreak analysis
            "PMID: 21713444",  # EAHEC genomic analysis
            "PMID: 35921067",  # EAHEC virulence mechanisms
            "PMID: 11446654",  # EHEC virulence factors
            "PMID: 29259590",  # EIEC invasion mechanisms
            "PMID: 29084270",  # EPEC attaching/effacing lesions
            "PMID: 11555281",  # EAEC biofilm formation
            "PMID: 34832511"   # ExPEC virulence evolution
        ],
        
        "ANTIMICROBIAL_RESISTANCE": [
            "PMID: 31856254",  # ST131 resistance mechanisms
            "PMID: 31437486",  # Global AMR surveillance
            "PMID: 36838380",  # ESBL and carbapenemase epidemiology
            "PMID: 32435624",  # STEC treatment considerations
            "PMID: 34202216",  # Plasmid-mediated resistance
            "PMID: 33411745",  # Colistin resistance mechanisms
            "PMID: 40859410",  # Carbapenemase evolution
            "PMID: 39375848"   # One Health AMR perspectives
        ],
        
        "GENOMICS_EVOLUTION": [
            "PMID: 32848072",  # EHEC genome plasticity
            "PMID: 14562958",  # EIEC/Shigella evolution
            "PMID: 35187883",  # EPEC pathoadaptation
            "PMID: 6099090",   # E. coli population genetics
            "PMID: 28766584",  # EPEC lineage diversification
            "PMID: 21078854",  # ETEC genomic diversity
            "PMID: 31431529"   # Horizontal gene transfer
        ],
        
        "DIAGNOSTICS_TYPING": [
            "PMID: 36668830",  # STEC diagnostic approaches
            "PMID: 12149382",  # EAEC detection methods
            "PMID: 25143581",  # Molecular typing schemes
            "PMID: 19527295",  # EPEC pathotyping
            "PMID: 38921808",  # STEC surveillance
            "PMID: 39530852"   # Genomic epidemiology
        ]
    },
    
    # =========================================================================
    # PART 2: DOI REFERENCES  
    # =========================================================================
    "DOI_REFERENCES": {
        "COMPREHENSIVE_REVIEWS": [
            "https://doi.org/10.1038/nrmicro818",    # E. coli pathotype evolution
            "https://doi.org/10.1128/cmr.11.1.142",  # E. coli diarrheagenic pathotypes
            "https://doi.org/10.1128/cmr.00135-18",  # ExPEC pathogenesis
            "https://doi.org/10.1128/cmr.00101-19",  # E. coli metabolism and virulence
            "https://doi.org/10.1128/msystems.01700-24",  # Antimicrobial resistance evolution
            "https://doi.org/10.3390/antibiotics13070662", # One Health AMR perspectives
            "https://doi.org/10.1128/iai.00368-25"         # E. coli vaccine development
        ],
        
        "GENOME_ANALYSES": [
            "https://doi.org/10.1128/msystems.00105-25",  # Global E. coli genomic diversity
            "https://doi.org/10.1128/mbio.02162-15",      # ST131 genome evolution
            "https://doi.org/10.1128/msphere.00532-24",   # E. coli pangenome analysis
            "https://doi.org/10.3390/antibiotics14050506",# Plasmid-mediated AMR spread
            "https://doi.org/10.1128/AEM.03771-15",       # E. coli population structure
            "https://doi.org/10.1128/jcm.01309-19"        # Hybrid pathotype emergence
        ],
        
        "PATHOGENESIS_STUDIES": [
            "https://doi.org/10.1046/j.1462-5822.2003.00281.x",  # EPEC attaching/effacing mechanism
            "https://doi.org/10.1128/microbiolspec.ehec-0008-2013",    # EHEC Shiga toxin action
            "https://doi.org/10.3389/fmicb.2019.01965",  # Bacterial secretion systems
            "https://doi.org/10.1128/mbio.01070-18",       # UPEC intracellular communities
            "DOI: 10.1016/S1413-8670(11)70158-1",       # ETEC colonization factors
            "https://doi.org/10.1128/jcm.02572-12"        # EAEC biofilm regulation
        ],
        
        "RESISTANCE_MECHANISMS": [
            "https://doi.org/10.1128/msphere.00108-21",    # β-lactamase evolution
            "https://doi.org/10.1128/microbiolspec.arba-0026-2017",  # AMR gene transfer
            "https://doi.org/10.1128/aac.02201-18",  # Plasmid dynamics
            "https://doi.org/10.1186/s12866-024-03215-6",  # Carbapenemase spread
            "https://doi.org/10.1128/mbio.01191-16", # Colistin resistance
            "https://doi.org/10.1128/jcm.06002-11"  # ESBL epidemiology
        ],
        
        "METHODS_TYPING": [
            "https://doi.org/10.1128/msphere.00738-20",     # Genomic epidemiology methods
            "https://doi.org/10.1002/vms3.1101",            # Molecular typing approaches
            "https://doi.org/10.1128/msystems.01236-22",    # Phylogenetic analysis
            "https://doi.org/10.1038/s41564-022-01079-y",   # Population genomics
            "https://doi.org/10.1128/mra.01095-22",         # Hybrid assembly approaches
            "https://doi.org/10.1128/AEM.01686-06"          # Pathotype prediction algorithms
        ],
        
        "OUTBREAK_INVESTIGATIONS": [
            "https://doi.org/10.4315/0362-028X.JFP-11-452",  # EAHEC O104:H4 outbreak
            "DOI: 10.1016/j.medine.2012.10.009",             # STEC O104 clinical features
            "https://doi.org/10.1186/1297-9716-43-13",       # EHEC O157 outbreaks
            "https://doi.org/10.1128/mbio.00377-13",         # ST131 healthcare outbreaks
            "https://doi.org/10.3390/ani14172490",           # ESBL-producing E. coli
            "https://doi.org/10.1128/aac.02571-13"           # Carbapenemase outbreaks
        ]
    }
}


# =============================================================================
# COMBINED DATABASE CLASS
# =============================================================================

class EcoliLineageDB:
    """Comprehensive E. coli lineage and pathotype database for typing and surveillance"""
    
    def __init__(self):
        self.lineages = LINEAGE_DATABASE
        self.serotypes = SEROTYPE_DATABASE
        self.phylogroups = PHYLOGROUP_DATABASE
        self.pathotypes = PATHOTYPE_DATABASE
        self.specialized_profiles = SPECIALIZED_PROFILES
        self.references = COMPREHENSIVE_REFERENCES
    
    def get_lineage_by_st(self, st: int) -> dict:
        """Get lineage data by sequence type"""
        st_str = f"ST{st}"
        return self.lineages.get(st_str, {})
    
    def get_pathotype_by_name(self, pathotype: str) -> dict:
        """Get pathotype data by name"""
        return self.pathotypes.get(pathotype, {})
    
    def get_pathotypes_by_category(self, category: str) -> dict:
        """Get all pathotypes of a specific category"""
        return {pt: data for pt, data in self.pathotypes.items() if data.get("category") == category}
    
    def predict_pathotype(self, virulence_genes: list, serotype: str = None) -> dict:
        """Predict pathotype based on virulence genes and optional serotype"""
        predictions = {}
        
        for pt_name, pt_data in self.pathotypes.items():
            score = 0
            matched_genes = []
            
            # Check key virulence genes
            key_genes = pt_data.get("key_virulence_genes", [])
            for gene in key_genes:
                if gene in virulence_genes:
                    score += 1
                    matched_genes.append(gene)
            
            # Check subtype markers for EPEC
            if pt_name == "EPEC":
                if "eae" in virulence_genes:
                    score += 2
                    if "bfpA" in virulence_genes:
                        subtype = "tEPEC"
                    else:
                        subtype = "aEPEC"
                    matched_genes.append(f"subtype: {subtype}")
            
            # Check for EHEC
            if pt_name == "EHEC":
                stx_genes = [g for g in virulence_genes if g.startswith('stx')]
                if stx_genes:
                    score += 2
                    matched_genes.extend(stx_genes)
            
            # Check serotype if provided
            if serotype and "serotypes" in pt_data:
                common_serotypes = pt_data["serotypes"].get("common", [])
                if serotype in common_serotypes:
                    score += 1
                    matched_genes.append(f"serotype_match: {serotype}")
            
            if score > 0:
                predictions[pt_name] = {
                    "score": score,
                    "matched_genes": matched_genes,
                    "confidence": self._get_confidence_level(score, len(key_genes)),
                    "pathotype_data": {
                        "primary_name": pt_data["primary_name"],
                        "category": pt_data["category"],
                        "risk_level": pt_data["risk_level"]
                    }
                }
        
        return dict(sorted(predictions.items(), key=lambda x: x[1]["score"], reverse=True))
    
    def _get_confidence_level(self, score: int, total_key_genes: int) -> str:
        """Determine confidence level based on matching score"""
        if total_key_genes == 0:
            return "LOW"
        
        ratio = score / total_key_genes
        if ratio >= 0.7:
            return "HIGH"
        elif ratio >= 0.4:
            return "MEDIUM"
        else:
            return "LOW"
    
    def get_specialized_profiles(self, profile_type: str) -> dict:
        """Get specialized pathotype profiles"""
        return self.specialized_profiles.get(profile_type, {})
    
    def export_complete_database(self) -> dict:
        """Export complete database"""
        return {
            "lineages": self.lineages,
            "serotypes": self.serotypes,
            "phylogroups": self.phylogroups,
            "pathotypes": self.pathotypes,
            "specialized_profiles": self.specialized_profiles,
            "references": self.references
        }

# =============================================================================
# VALIDATION AND UTILITY FUNCTIONS
# =============================================================================

def validate_complete_database():
    """Validate the integrity of the complete database"""
    issues = []
    
    # Validate lineages
    for st, data in LINEAGE_DATABASE.items():
        required = ["primary_name", "category", "phylogroup", "serotype", "pathotypes", "risk_level"]
        for field in required:
            if field not in data:
                issues.append(f"Missing {field} in {st}")
    
    # Validate pathotypes
    for pt, data in PATHOTYPE_DATABASE.items():
        required = ["primary_name", "category", "key_virulence_genes", "clinical_manifestations", "risk_level"]
        for field in required:
            if field not in data:
                issues.append(f"Missing field '{field}' in {pt}")
    
    return issues

def export_database_summary():
    """Export summary of complete database"""
    summary = {
        "lineages": len(LINEAGE_DATABASE),
        "serotypes": len(SEROTYPE_DATABASE),
        "phylogroups": len(PHYLOGROUP_DATABASE),
        "pathotypes": len(PATHOTYPE_DATABASE),
        "specialized_profiles": len(SPECIALIZED_PROFILES)
    }
    
    # Count pathotypes by category
    pathotype_categories = {}
    for pt, data in PATHOTYPE_DATABASE.items():
        category = data["category"]
        pathotype_categories[category] = pathotype_categories.get(category, 0) + 1
    
    summary["pathotype_categories"] = pathotype_categories
    return summary

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("=== E. COLI COMPREHENSIVE DATABASE ===")
    
    # Validate database
    issues = validate_complete_database()
    if issues:
        print(f"Validation issues: {len(issues)}")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print("✓ Database validation passed")
    
    # Create instance
    db = EcoliLineageDB()
    
    # Test queries
    print(f"\n=== DATABASE SUMMARY ===")
    summary = export_database_summary()
    print(f"Lineages: {summary['lineages']}")
    print(f"Serotypes: {summary['serotypes']}") 
    print(f"Phylogroups: {summary['phylogroups']}")
    print(f"Pathotypes: {summary['pathotypes']}")
    
    print(f"\nPathotypes by category:")
    for category, count in summary['pathotype_categories'].items():
        print(f"  {category}: {count}")
    
    # Test pathotype prediction
    print(f"\n=== TESTING PATHOTYPE PREDICTION ===")
    test_genes = ["stx1", "stx2", "eae", "ehxA", "fimH"]
    predictions = db.predict_pathotype(test_genes, "O157:H7")
    print(f"Test genes: {test_genes}")
    for pt, info in list(predictions.items())[:3]:
        print(f"  {pt}: score={info['score']}, confidence={info['confidence']}")
    
    print(f"\n✅ EcoliLineageDB - Comprehensive E. coli Database Ready!")
