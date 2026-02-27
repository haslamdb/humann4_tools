#\!/usr/bin/env python3
"""
HUMAnN3 Tools Differential Abundance Module

This module runs differential abundance analysis on HUMAnN3 output files.
It supports multiple statistical methods for comparing pathway or gene abundances between groups.
"""

import os
import sys
import time
import logging
import argparse
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Union

# Import internal modules
try:
    from src.humann4_tools.analysis.differential_abundance import run_differential_abundance_analysis
except ImportError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
    from src.humann4_tools.analysis.differential_abundance import run_differential_abundance_analysis

# Set up logging
logger = logging.getLogger('humann4_tools')

def setup_logger(log_file=None, log_level=logging.INFO):
    """Set up the logger with console and optional file output."""
    # Remove any existing handlers to avoid duplication
    logger.handlers = []
    logger.setLevel(log_level)

    # Format for logs
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Optional file handler
    if log_file:
        os.makedirs(os.path.dirname(os.path.abspath(log_file)), exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger

def parse_args():
    """Parse command line arguments for the Differential Abundance module."""
    parser = argparse.ArgumentParser(
        description="Run differential abundance analysis on HUMAnN3 output files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Description:
  Differential abundance module performs sophisticated differential abundance testing to identify
  pathways or genes that significantly differ between experimental groups. It uses robust
  compositional data analysis methods that account for the special properties of microbiome data.

Key Features:
  • Supports multiple differential abundance methods:
    - ALDEx2: Uses Bayesian approach with Monte-Carlo instances of the Dirichlet distribution
    - ANCOM: Analysis of composition of microbiomes, robust to compositional effects
    - ANCOM-BC: ANCOM with bias correction for uneven sampling depth
  • Handles both pathway and gene family data
  • Accounts for compositional nature of microbiome data
  • Option to filter by specific groups of interest

Output Files:
  • [method]_results.csv: Results from each differential abundance method with statistics
  • diff_abundance_summary.txt: Summary of analysis with significant feature counts
  • Visualizations of significant features (when applicable)

Common Usage:
  # Basic differential abundance analysis:
  humann4-tools diff --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv --metadata-file metadata.csv

  # For gene family data:
  humann4-tools diff --abundance-file joined_output/genefamilies_cpm_unstratified.tsv --metadata-file metadata.csv --feature-type gene

  # Specify specific groups to compare:
  humann4-tools diff --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv --metadata-file metadata.csv --filter-groups Control,Treatment
"""
    )
    
    # Required arguments
    parser.add_argument("--abundance-file", required=True, 
                      help="Path to the unstratified abundance file (pathway or gene family)")
    parser.add_argument("--metadata-file", required=True,
                      help="Path to sample metadata CSV file")
    
    # Analysis options
    parser.add_argument("--output-dir", default="./DifferentialAbundance",
                      help="Directory for output files")
    parser.add_argument("--feature-type", choices=["pathway", "gene"], default="pathway",
                      help="Type of features in the abundance file (pathway or gene)")
    parser.add_argument("--group-col", default="Group",
                      help="Column name in metadata for grouping samples")
    parser.add_argument("--sample-id-col", 
                      help="Column name in metadata for sample IDs (autodetected if not specified)")
    parser.add_argument("--methods", default="aldex2,ancom,ancom-bc",
                      help="Comma-separated list of methods to use (aldex2,ancom,ancom-bc)")
    parser.add_argument("--exclude-unmapped", action="store_true",
                      help="Exclude unmapped features from analysis")
    parser.add_argument("--filter-groups",
                      help="Comma-separated list of group names to include in the analysis. "
                            "For ALDEx2, exactly 2 groups must be specified.")
    parser.add_argument("--alpha", type=float, default=0.05,
                      help="Significance threshold for statistical tests (default: 0.05)")
    
    # Additional options
    parser.add_argument("--log-file", 
                      help="Path to log file")
    parser.add_argument("--log-level", default="INFO", 
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                      help="Logging level")
    
    return parser.parse_args()

def main():
    """Main function to run differential abundance analysis."""
    # Parse arguments
    args = parse_args()
    
    # Setup logging
    log_level = getattr(logging, args.log_level.upper())
    setup_logger(args.log_file, log_level)
    
    logger.info("Starting HUMAnN3 Tools Differential Abundance Module")
    start_time = time.time()
    
    # Process methods
    methods = [m.strip().lower() for m in args.methods.split(',')]
    logger.info(f"Using methods: {', '.join(methods)}")
    
    # Process filter groups
    filter_groups = None
    if args.filter_groups:
        filter_groups = [g.strip() for g in args.filter_groups.split(',')]
        logger.info(f"Filtering to groups: {filter_groups}")
    
    # Check if files exist
    for file_path, desc in [(args.abundance_file, "Abundance file"), (args.metadata_file, "Metadata file")]:
        if not os.path.exists(file_path):
            logger.error(f"ERROR: {desc} not found: {file_path}")
            return 1
    
    # Set denominator based on exclude_unmapped flag
    denom = "unmapped_excluded" if args.exclude_unmapped else "all"
    
    # Run differential abundance analysis
    success = run_differential_abundance_analysis(
        abundance_file=args.abundance_file,
        metadata_file=args.metadata_file,
        output_dir=args.output_dir,
        group_col=args.group_col,
        methods=methods,
        feature_type=args.feature_type,
        denom=denom,
        filter_groups=filter_groups,
        sample_id_col=args.sample_id_col,
        alpha=args.alpha
    )
    
    if not success:
        logger.error("Differential abundance analysis failed")
        return 1
    
    # Print elapsed time
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    
    logger.info(f"Total processing time: {int(minutes)}m {int(seconds)}s")
    
    # Print next steps
    logger.info("\nNext Steps:")
    logger.info("  For visualizations:")
    logger.info(f"  humann4-tools viz --abundance-file {args.abundance_file} --metadata-file {args.metadata_file}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
