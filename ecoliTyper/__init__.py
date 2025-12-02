"""
EcoliTyper - Comprehensive E. coli Typing Pipeline
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "Brown Beckley"
__email__ = "brownbeckley94@gmail.com"
__affiliation__ = "University of Ghana Medical School"

# Import the main orchestrator
from .ecolityper import EcoliTyperOrchestrator, main

__all__ = ["EcoliTyperOrchestrator", "main"]
