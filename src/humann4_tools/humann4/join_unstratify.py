# humann4_tools/humann4_tools/join_unstratify.py
import os
import sys
import argparse
import logging
import time

from src.humann4_tools.logger import setup_logger, log_print
from src.humann4_tools.utils.sample_utils import validate_sample_key, check_input_files_exist
from src.humann4_tools.humann4.pathway_processing import process_pathway_abundance
from src.humann4_tools.humann4.gene_processing import process_gene_families

def parse_args():
    """Parse command line arguments for join_unstratify_humann_output."""
    parser = argparse.ArgumentParser(
        description="Join and unstratify HUMAnN3 gene family and pathway abundance files without running downstream analysis"
    )
    
    # Required arguments
    parser.add_argument("--sample-key", required=True, help="CSV file with columns for sample names and metadata")
    parser.add_argument("--pathway-dir", required=True, help="Directory containing raw pathway abundance files")
    parser.add_argument("--gene-dir", required=True, help="Directory containing raw gene family files")
    
    # Output options
    parser.add_argument("--output-dir", default="./HUMAnN3_Processed", help="Directory for output files")
    parser.add_argument("--output-prefix", default="ProcessedFiles", help="Prefix for output files")
    
    # Processing options
    parser.add_argument("--skip-pathway", action="store_true", help="Skip pathway processing")
    parser.add_argument("--skip-gene", action="store_true", help="Skip gene family processing")
    parser.add_argument("--units", default="cpm", choices=["cpm", "relab"], help="Units for normalization (default: cpm)")
    
    # Additional options
    parser.add_argument("--no-interactive", action="store_true", help="Non-interactive mode for sample key column selection")
    parser.add_argument("--log-file", default=None, help="Path to log file")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], 
                      help="Logging level (default: INFO)")
    
    return parser.parse_args()

def process_join_unstratify(
    sample_key, 
    pathway_dir=None, 
    gene_dir=None, 
    output_dir="./HUMAnN3_Processed", 
    output_prefix="ProcessedFiles", 
    units="cpm", 
    no_interactive=False, 
    log_file=None, 
    log_level="INFO"
):
    """
    Function to join and unstratify HUMAnN3 gene family and pathway files.
    This function only processes the HUMAnN3 output files without running 
    preprocessing or downstream analysis.
    
    Args:
        sample_key: Path to sample metadata CSV file
        pathway_dir: Directory containing pathway files (None to skip)
        gene_dir: Directory containing gene family files (None to skip)
        output_dir: Directory for output files
        output_prefix: Prefix for output files
        units: Units for normalization ("cpm" or "relab")
        no_interactive: Disable interactive prompts
        log_file: Path to log file
        log_level: Logging level
        
    Returns:
        Boolean success flag
    """
    # Setup logging
    logger = setup_logger(log_file=log_file, log_level=getattr(logging, log_level.upper()))
    log_print("Starting HUMAnN3 Join and Unstratify Pipeline", level="info")
    start_time = time.time()
    
    # Check if at least one directory is provided
    if not pathway_dir and not gene_dir:
        log_print("ERROR: At least one of pathway_dir or gene_dir must be provided", level="error")
        return False
    
    # Process sample metadata
    samples, selected_columns = validate_sample_key(sample_key, no_interactive=no_interactive)
    
    # Check input files
    valid_path_samples, valid_gene_samples = [], []
    
    if pathway_dir:
        valid_path_samples, _ = check_input_files_exist(samples, pathway_dir, gene_dir if gene_dir else pathway_dir)
    
    if gene_dir:
        _, valid_gene_samples = check_input_files_exist(samples, pathway_dir if pathway_dir else gene_dir, gene_dir)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Process pathway files
    pathway_unstrat_file = None
    skip_pathway = pathway_dir is None
    
    if not skip_pathway and valid_path_samples:
        log_print(f"Processing pathway abundance files with units: {units}...", level="info")
        pathway_unstrat_file = process_pathway_abundance(
            valid_path_samples,
            pathway_dir,
            output_dir,
            output_prefix,
            selected_columns=selected_columns,
            units=units
        )
        
        if pathway_unstrat_file:
            log_print(f"Pathway processing completed. Unstratified file: {pathway_unstrat_file}", level="info")
        else:
            log_print("Pathway processing failed to produce an unstratified file", level="warning")
    else:
        if skip_pathway:
            log_print("Skipping pathway processing as no pathway_dir provided", level="info")
        else:
            log_print("No valid pathway files found; skipping pathway processing", level="warning")
    
    # Process gene family files
    gene_unstrat_file = None
    skip_gene = gene_dir is None
    
    if not skip_gene and valid_gene_samples:
        log_print(f"Processing gene family files with units: {units}...", level="info")
        gene_unstrat_file = process_gene_families(
            valid_gene_samples,
            gene_dir,
            output_dir,
            output_prefix,
            selected_columns=selected_columns,
            units=units
        )
        
        if gene_unstrat_file:
            log_print(f"Gene family processing completed. Unstratified file: {gene_unstrat_file}", level="info")
        else:
            log_print("Gene family processing failed to produce an unstratified file", level="warning")
    else:
        if skip_gene:
            log_print("Skipping gene family processing as no gene_dir provided", level="info")
        else:
            log_print("No valid gene family files found; skipping gene family processing", level="warning")
    
    # Report summary
    elapsed = time.time() - start_time
    mm, ss = divmod(elapsed, 60)
    
    log_print(f"Pipeline completed in {int(mm)}m {int(ss)}s", level="info")
    log_print("Summary:", level="info")
    
    if pathway_unstrat_file:
        log_print(f"  - Pathway unstratified file: {pathway_unstrat_file}", level="info")
    if gene_unstrat_file:
        log_print(f"  - Gene family unstratified file: {gene_unstrat_file}", level="info")
    
    success = pathway_unstrat_file is not None or gene_unstrat_file is not None
    if not success:
        log_print("  - No output files were generated", level="warning")
    
    return success

def join_unstratify_humann_output():
    """
    Main function to join and unstratify HUMAnN3 gene family and pathway files.
    This function only processes the HUMAnN3 output files without running 
    preprocessing or downstream analysis.
    """
    args = parse_args()
    
    # Use the process function instead of duplicating code
    success = process_join_unstratify(
        sample_key=args.sample_key,
        pathway_dir=None if args.skip_pathway else args.pathway_dir,
        gene_dir=None if args.skip_gene else args.gene_dir,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        units=args.units,
        no_interactive=args.no_interactive,
        log_file=args.log_file,
        log_level=args.log_level
    )
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(join_unstratify_humann_output())