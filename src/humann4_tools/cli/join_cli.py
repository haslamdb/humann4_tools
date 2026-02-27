#!/usr/bin/env python3
"""
HUMAnN3 Tools Join Module

This module joins, normalizes, and unstratifies HUMAnN3 output files.
It can process:
- Pathway abundance files
- Pathway coverage files
- Gene family files

Examples:
  # Join and normalize pathway abundance files:
  humann4-tools join --input-dir ./PathwayAbundance --pathabundance --output-dir ./joined_output --units cpm
  
  # Join and normalize gene family files:
  humann4-tools join --input-dir ./GeneFamilies --genefamilies --output-dir ./joined_output --units cpm
  
  # Join pathway coverage files (no normalization needed):
  humann4-tools join --input-dir ./PathwayCoverage --pathcoverage --output-dir ./joined_output
  
  # With sample name updates:
  humann4-tools join --input-dir ./PathwayAbundance --pathabundance --output-dir ./joined_output --update-snames
"""

import os
import sys
import time
import argparse
import logging
import subprocess
import glob
from typing import Dict, List, Optional, Tuple, Union

# Import internal modules
try:
    from src.humann4_tools.utils.cmd_utils import run_cmd
    from src.humann4_tools.utils.resource_utils import track_peak_memory
    from src.humann4_tools.utils.file_utils import strip_suffixes_from_file_headers
except ImportError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
    from src.humann4_tools.utils.cmd_utils import run_cmd
    from src.humann4_tools.utils.resource_utils import track_peak_memory
    from src.humann4_tools.utils.file_utils import strip_suffixes_from_file_headers

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

def check_humann_util_installation(util_name: str) -> Tuple[bool, str]:
    """Check if a HUMAnN utility is installed and available."""
    try:
        result = subprocess.run([util_name, "--help"], 
                               capture_output=True, text=True, check=False)
        if result.returncode == 0:
            return True, f"{util_name} is available"
        return False, f"{util_name} command exists but returned an error"
    except FileNotFoundError:
        return False, f"{util_name} not found in PATH"

def join_normalize_tables(
    input_dir: str,
    output_dir: str,
    file_type: str,
    units: Optional[str] = None,
    output_basename: Optional[str] = None,
    update_snames: bool = False,
    file_pattern: Optional[str] = None,
    strip_headers: bool = True
) -> Optional[Dict[str, str]]:
    """
    Join, normalize, and unstratify HUMAnN3 output files.
    
    Args:
        input_dir: Directory containing HUMAnN3 output files
        output_dir: Directory for output files
        file_type: Type of files to process (pathabundance, pathcoverage, genefamilies)
        units: Units for normalization (cpm, relab)
        output_basename: Base filename for output (default derived from file type)
        update_snames: Whether to update sample names during normalization
        file_pattern: Pattern for input files
        strip_headers: Whether to strip suffixes from column headers
        
    Returns:
        Dictionary mapping output types (unstratified, stratified) to file paths,
        or None if processing fails
    """
    
    # Determine file type and pattern
    if file_type == "pathabundance":
        default_pattern = "*pathabundance.tsv"
        utility_description = "pathway abundance"
    elif file_type == "pathcoverage":
        default_pattern = "*pathcoverage.tsv"
        utility_description = "pathway coverage"
    elif file_type == "genefamilies":
        default_pattern = "*genefamilies.tsv"
        utility_description = "gene families"
    else:
        logger.error(f"Invalid file type: {file_type}")
        return None
    
    # Use provided pattern or default
    pattern = file_pattern if file_pattern else default_pattern
    
    # Find input files
    search_pattern = os.path.join(input_dir, pattern)
    input_files = glob.glob(search_pattern)
    
    if not input_files:
        logger.error(f"No files found matching pattern: {search_pattern}")
        return None
    
    logger.info(f"Found {len(input_files)} {utility_description} files")
    
    # Create normalized files directory (except for pathcoverage)
    need_normalization = file_type != "pathcoverage" and units
    
    if need_normalization:
        norm_dir = os.path.join(output_dir, "normalized")
        os.makedirs(norm_dir, exist_ok=True)
        
        # Normalize each file
        logger.info(f"Normalizing files to {units} units...")
        normalized_files = []
        
        for input_file in input_files:
            basename = os.path.basename(input_file)
            sample_name = basename.replace(f"_{file_type}.tsv", "").replace(f".{file_type}.tsv", "")
            output_norm = os.path.join(norm_dir, f"{sample_name}_{file_type}_{units}.tsv")
            
            cmd = [
                "humann_renorm_table",
                "--input", input_file,
                "--output", output_norm,
                "--units", units
            ]
            
            if update_snames:
                cmd.append("--update-snames")
            
            success = run_cmd(cmd, exit_on_error=False)
            if success:
                normalized_files.append(output_norm)
                logger.debug(f"Normalized: {os.path.basename(input_file)} → {os.path.basename(output_norm)}")
            else:
                logger.warning(f"Failed to normalize: {os.path.basename(input_file)}")
        
        if not normalized_files:
            logger.error("No files were successfully normalized")
            return None
        
        logger.info(f"Successfully normalized {len(normalized_files)} files")
        
        # Use normalized files for joining
        files_to_join = norm_dir
    else:
        # Skip normalization for pathcoverage
        files_to_join = input_dir
    
    # Set output basename if not provided
    if not output_basename:
        if need_normalization:
            output_basename = f"{file_type}_{units}"
        else:
            output_basename = file_type
    
    # Join files
    logger.info("Joining files...")
    joined_output = os.path.join(output_dir, f"{output_basename}.tsv")
    
    join_cmd = [
        "humann_join_tables",
        "-i", files_to_join,
        "-o", joined_output
    ]
    
    success = run_cmd(join_cmd, exit_on_error=False)
    if not success:
        logger.error("Failed to join files")
        return None
    
    logger.info(f"Successfully joined files: {joined_output}")
    
    # Split stratified table
    logger.info("Splitting stratified table...")
    
    split_cmd = [
        "humann_split_stratified_table",
        "-i", joined_output,
        "-o", output_dir
    ]
    
    success = run_cmd(split_cmd, exit_on_error=False)
    if not success:
        logger.error("Failed to split stratified table")
        return None
    
    # Find unstratified and stratified files
    unstrat_pattern = os.path.join(output_dir, f"*{output_basename}*_unstratified.tsv")
    strat_pattern = os.path.join(output_dir, f"*{output_basename}*_stratified.tsv")
    
    unstrat_files = glob.glob(unstrat_pattern)
    strat_files = glob.glob(strat_pattern)
    
    # Standardize file names
    output_files = {}
    
    if unstrat_files:
        unstrat_file = unstrat_files[0]
        standard_unstrat = os.path.join(output_dir, f"{output_basename}_unstratified.tsv")
        
        try:
            if os.path.abspath(unstrat_file) != os.path.abspath(standard_unstrat):
                os.rename(unstrat_file, standard_unstrat)
                logger.info(f"Renamed unstratified file to: {os.path.basename(standard_unstrat)}")
                
            output_files['unstratified'] = standard_unstrat
            
            # Strip suffixes from headers if requested
            if strip_headers:
                strip_suffixes_from_file_headers(standard_unstrat)
                logger.info("Stripped suffixes from unstratified file headers")
                
        except Exception as e:
            logger.warning(f"Error standardizing unstratified filename: {str(e)}")
            output_files['unstratified'] = unstrat_file
    else:
        logger.warning("Could not find unstratified output file")
    
    if strat_files:
        strat_file = strat_files[0]
        standard_strat = os.path.join(output_dir, f"{output_basename}_stratified.tsv")
        
        try:
            if os.path.abspath(strat_file) != os.path.abspath(standard_strat):
                os.rename(strat_file, standard_strat)
                logger.info(f"Renamed stratified file to: {os.path.basename(standard_strat)}")
                
            output_files['stratified'] = standard_strat
            
            # Strip suffixes from headers if requested
            if strip_headers:
                strip_suffixes_from_file_headers(standard_strat)
                logger.info("Stripped suffixes from stratified file headers")
                
        except Exception as e:
            logger.warning(f"Error standardizing stratified filename: {str(e)}")
            output_files['stratified'] = strat_file
    else:
        logger.warning("Could not find stratified output file")
    
    return output_files

def parse_args(args=None, parent_parser=None):
    """
    Parse command line arguments for the Join module.
    
    Args:
        args: List of arguments to parse (default: sys.argv[1:])
        parent_parser: Parent parser to add arguments to (default: create new parser)
        
    Returns:
        Parsed arguments if args is provided, otherwise the configured parser
    """
    if parent_parser:
        # Use the provided parser
        parser = parent_parser
    else:
        # Create a new parser
        parser = argparse.ArgumentParser(
            description="Join, normalize, and unstratify HUMAnN3 output files",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Description:
  Join module combines multiple HUMAnN3 output files into unified tables, applies
  normalization, and separates stratified (species-level) from unstratified (community-level)
  data. This is a critical step for preparing data for statistical testing and visualization.

Key Features:
  • Combines outputs from multiple samples into a single table
  • Normalizes abundance data to counts per million (CPM) or relative abundance
  • Generates both stratified and unstratified tables
  • Handles all three HUMAnN3 output types (pathway abundance, pathway coverage, gene families)

Output Files:
  • *_stratified.tsv: Combined table with species contributions preserved
  • *_unstratified.tsv: Combined table with only community-level totals

Common Usage:
  # For pathway abundance files with CPM normalization:
  humann4-tools join --input-dir PathwayAbundance --pathabundance --output-dir joined_output --units cpm

  # For gene families with relative abundance normalization:
  humann4-tools join --input-dir GeneFamilies --genefamilies --output-dir joined_output --units relab

  # Next step after joining:
  humann4-tools stats --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv --metadata-file metadata.csv
"""
        )
    
    # Required arguments
    input_group = parser.add_argument_group("Required Options")
    input_group.add_argument("--input-dir", required=True, 
                      help="Directory containing HUMAnN3 output files (e.g., PathwayAbundance/)")
    
    # File type options
    file_group = parser.add_argument_group("File Type Options (choose one)")
    type_group = file_group.add_mutually_exclusive_group(required=True)
    type_group.add_argument("--pathabundance", action="store_true", 
                     help="Process pathway abundance files (*pathabundance.tsv)")
    type_group.add_argument("--pathcoverage", action="store_true", 
                     help="Process pathway coverage files (*pathcoverage.tsv)")
    type_group.add_argument("--genefamilies", action="store_true", 
                     help="Process gene family files (*genefamilies.tsv)")
    
    # Output options
    output_group = parser.add_argument_group("Output Options")
    output_group.add_argument("--output-dir", default="./joined_output",
                      help="Directory for output files (default: ./joined_output)")
    output_group.add_argument("--output-basename", 
                      help="Base filename for output (default: derived from file type, e.g., 'pathway_abundance_cpm')")
    output_group.add_argument("--units", default="cpm", choices=["cpm", "relab"],
                      help="Normalization units: cpm (counts per million) or relab (relative abundance) (default: cpm)")
    
    # Additional options
    format_group = parser.add_argument_group("Format Options")
    format_group.add_argument("--update-snames", action="store_true",
                      help="Update sample names during normalization to remove file suffixes")
    format_group.add_argument("--file-pattern", 
                      help="Glob pattern for input files (default is based on file type, e.g., '*pathabundance.tsv')")
    format_group.add_argument("--no-strip-headers", action="store_true",
                      help="Don't strip suffixes from column headers (keep full sample names)")
    
    # Logging options
    log_group = parser.add_argument_group("Logging Options")
    log_group.add_argument("--log-file", 
                      help="Path to log file (if not specified, logs to console only)")
    log_group.add_argument("--log-level", default="INFO",
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                      help="Logging level verbosity (default: INFO)")
    
    # If args is provided, parse and return args
    if args is not None:
        return parser.parse_args(args)
    
    # Otherwise return the parser
    return parser

def main(args=None):
    """
    Main function to run join, normalize, and unstratify operations.
    
    Args:
        args: Command line arguments (optional)
    """
    # Parse arguments
    if isinstance(args, list):
        args = parse_args(args)
    else:
        args = parse_args()
    
    # Setup logging
    log_level = getattr(logging, args.log_level.upper())
    setup_logger(args.log_file, log_level)
    
    logger.info("Starting HUMAnN3 Tools Join Module")
    start_time = time.time()
    
    # Check required utilities
    utils_to_check = [
        "humann_join_tables",
        "humann_split_stratified_table"
    ]
    
    # Add humann_renorm_table if normalization is needed
    if (args.pathabundance or args.genefamilies) and args.units:
        utils_to_check.append("humann_renorm_table")
    
    # Check all required utilities
    for util in utils_to_check:
        util_ok, util_msg = check_humann_util_installation(util)
        if not util_ok:
            logger.error(f"{util} not properly installed: {util_msg}")
            return 1
    
    # Determine file type
    file_type = None
    if args.pathabundance:
        file_type = "pathabundance"
    elif args.pathcoverage:
        file_type = "pathcoverage"
    elif args.genefamilies:
        file_type = "genefamilies"
    
    # Process files
    results = join_normalize_tables(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        file_type=file_type,
        units=args.units if file_type != "pathcoverage" else None,
        output_basename=args.output_basename,
        update_snames=args.update_snames,
        file_pattern=args.file_pattern,
        strip_headers=not args.no_strip_headers
    )
    
    if not results:
        logger.error("Join/normalize/unstratify process failed")
        return 1
    
    # Log results
    logger.info("Join/normalize/unstratify process completed successfully")
    
    for file_type, file_path in results.items():
        logger.info(f"{file_type.capitalize()} file: {file_path}")
    
    # Print elapsed time
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    
    logger.info(f"Total processing time: {int(minutes)}m {int(seconds)}s")
    
    # Print next steps
    logger.info("\nNext Steps:")
    if 'unstratified' in results:
        logger.info("  For statistical testing:")
        logger.info(f"  humann4-tools stats --abundance-file {results['unstratified']} --metadata-file metadata.csv")
        logger.info("\n  For differential abundance analysis:")
        logger.info(f"  humann4-tools diff --abundance-file {results['unstratified']} --metadata-file metadata.csv")
        logger.info("\n  For visualizations:")
        logger.info(f"  humann4-tools viz --abundance-file {results['unstratified']} --metadata-file metadata.csv")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())