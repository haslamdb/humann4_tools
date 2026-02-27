#!/usr/bin/env python3
"""
HUMAnN3 Tools HUMAnN3 Module

This module runs HUMAnN3 on preprocessed sequence files (typically from KneadData).
It supports four input methods:
1. Direct input files
2. Sample list file
3. Metadata-driven input
4. Directory of KneadData output files

Examples:
  # Basic usage with direct input files:
  humann4-tools humann4 --input-files sample1_cleaned.fastq --threads 8 --output-dir ./humann4_output
  
  # Using a directory of KneadData outputs:
  humann4-tools humann4 --input-dir ./kneaddata_output --output-dir ./humann4_output
  
  # Using a sample list file:
  humann4-tools humann4 --samples-file samples.txt --output-dir ./humann4_output
  
  # Using metadata:
  humann4-tools humann4 --metadata-file metadata.csv --seq-dir /path/to/cleaned_sequences --output-dir ./humann4_output
  
  # Specifying database paths:
  humann4-tools humann4 --input-files sample.fastq --nucleotide-db /path/to/chocophlan --protein-db /path/to/uniref --output-dir ./humann4_output
"""

import os
import sys
import time
import argparse
import logging
import subprocess
import multiprocessing
import glob
from typing import Dict, List, Optional, Tuple, Union

# Import internal modules
try:
    from src.humann4_tools.utils.input_handler import get_input_files, find_kneaddata_output_files
    from src.humann4_tools.utils.cmd_utils import run_cmd
    from src.humann4_tools.utils.resource_utils import track_peak_memory
except ImportError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
    from src.humann4_tools.utils.input_handler import get_input_files, find_kneaddata_output_files
    from src.humann4_tools.utils.cmd_utils import run_cmd
    from src.humann4_tools.utils.resource_utils import track_peak_memory

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

def check_humann4_installation() -> Tuple[bool, str]:
    """Check if HUMAnN3 is installed and available."""
    try:
        result = subprocess.run(["humann", "--version"], 
                               capture_output=True, text=True, check=False)
        if result.returncode == 0:
            return True, result.stdout.strip()
        return False, "HUMAnN3 command exists but returned an error"
    except FileNotFoundError:
        return False, "HUMAnN3 not found in PATH"

# [... rest of the original function definitions ...]

def parse_args(args=None, parent_parser=None):
    """
    Parse command line arguments for the HUMAnN3 module.
    
    Args:
        args: List of arguments to parse (default: sys.argv[1:])
        parent_parser: Parent parser to add arguments to (default: create new parser)
        
    Returns:
        Parsed arguments if args is provided, otherwise the configured parser
    """
    # Use existing parser or create new one
    parser = parent_parser or argparse.ArgumentParser(
        description="Run HUMAnN3 on preprocessed sequence files"
    )
    
    # Input methods (choose one)
    input_group = parser.add_argument_group("Input Options (choose one)")
    input_group.add_argument("--input-dir", 
                           help="Directory containing KneadData output files")
    input_group.add_argument("--input-files", nargs="+", 
                           help="Input FASTQ file(s) for HUMAnN3")
    input_group.add_argument("--samples-file", 
                           help="Tab-delimited file with sample IDs and file paths")
    input_group.add_argument("--metadata-file", 
                           help="CSV file with sample metadata")
    
    # Metadata-driven input options
    metadata_group = parser.add_argument_group("Metadata-driven Input Options")
    metadata_group.add_argument("--seq-dir", 
                             help="Directory containing sequence files (for metadata-driven)")
    metadata_group.add_argument("--sample-col", 
                             help="Column name for sample IDs in metadata")
    metadata_group.add_argument("--r1-col", 
                             help="Column name for R1 file paths in metadata")
    metadata_group.add_argument("--r2-col", 
                             help="Column name for R2 file paths in metadata")
    metadata_group.add_argument("--file-pattern", 
                             help="Pattern for finding files (e.g., {sample}_S*_R*.fastq.gz)")
    metadata_group.add_argument("--r1-suffix", 
                             help="Suffix for R1 files (e.g., _R1.fastq.gz)")
    metadata_group.add_argument("--r2-suffix", 
                             help="Suffix for R2 files (e.g., _R2.fastq.gz)")
    
    # Output options
    output_group = parser.add_argument_group("Output Options")
    output_group.add_argument("--output-dir", default="./humann4_output", 
                           help="Directory for output files (default: ./humann4_output)")
    output_group.add_argument("--organize-outputs", action="store_true",
                           help="Organize output files into type-specific directories")
    output_group.add_argument("--log-file", 
                           help="Path to log file")
    output_group.add_argument("--log-level", default="INFO", 
                           choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                           help="Logging level")
    
    # HUMAnN3 database options
    db_group = parser.add_argument_group("HUMAnN3 Database Options")
    db_group.add_argument("--nucleotide-db", 
                        help="Path to nucleotide database (ChocoPhlAn)")
    db_group.add_argument("--protein-db", 
                        help="Path to protein database (UniRef)")
    
    # HUMAnN3 processing options
    processing_group = parser.add_argument_group("HUMAnN3 Processing Options")
    processing_group.add_argument("--threads", type=int, default=1, 
                               help="Number of threads per sample (default: 1)")
    processing_group.add_argument("--paired", action="store_true", 
                               help="Input files are paired-end reads")
    processing_group.add_argument("--use-parallel", action="store_true",
                               help="Process samples in parallel")
    processing_group.add_argument("--max-parallel", type=int,
                               help="Maximum number of samples to process in parallel")
    
    # HUMAnN3 pipeline options
    pipeline_group = parser.add_argument_group("HUMAnN3 Pipeline Options")
    pipeline_group.add_argument("--bypass-prescreen", action="store_true",
                             help="Bypass the MetaPhlAn taxonomic prescreen")
    pipeline_group.add_argument("--bypass-nucleotide-index", action="store_true",
                             help="Bypass the nucleotide index database")
    pipeline_group.add_argument("--bypass-translated-search", action="store_true",
                             help="Bypass the translated search")
    pipeline_group.add_argument("--humann4-options", nargs="+",
                             help="Additional options to pass to HUMAnN3 (format: key=value)")
    
    # If args is provided, parse and return args, otherwise parse from sys.argv
    if args is not None:
        return parser.parse_args(args)
    else:
        return parser.parse_args()
    
def main(args=None):
    """
    Main function to run HUMAnN3 processing.
    
    Args:
        args: Command line arguments (optional)
    """
    # Handle different types of args input:
    # 1. If args is None: parse from command line
    # 2. If args is a list: parse those args
    # 3. If args is an argparse.Namespace: use directly
    if args is None:
        # No args provided, parse from command line
        args = parse_args()
    elif isinstance(args, list):
        # A list of args was provided
        args = parse_args(args)
    # else: args is already an argparse.Namespace object, use as is
    
    # Setup logging
    log_level = getattr(logging, args.log_level.upper())
    setup_logger(args.log_file, log_level)
    
    start_time = time.time()
    logger.info("Starting HUMAnN3 Tools HUMAnN3 Module")
    
    # Check HUMAnN3 installation
    humann4_ok, humann4_version = check_humann4_installation()
    if not humann4_ok:
        logger.error(f"HUMAnN3 not properly installed: {humann4_version}")
        return 1
    logger.info(f"Using HUMAnN3 version: {humann4_version}")
    
    # Process input methods to get sample dictionary
    samples = {}
    
    # Method 1: Direct input files
    if args.input_files:
        logger.info(f"Processing {len(args.input_files)} input files")
        
        # Create a basic sample from the input files
        if args.paired and len(args.input_files) >= 2:
            # For paired-end data, use the first two files
            r1_file = args.input_files[0]
            r2_file = args.input_files[1]
            sample_name = os.path.basename(r1_file).split('_')[0]
            
            samples[sample_name] = {
                'files': [r1_file, r2_file],
                'metadata': {}
            }
            logger.info(f"Created sample {sample_name} from paired input files")
        elif args.input_files:
            # For single-end data, treat each file as a separate sample
            for i, file in enumerate(args.input_files):
                sample_name = os.path.basename(file).split('.')[0]
                
                # Handle potential duplicates by adding an index
                if sample_name in samples:
                    sample_name = f"{sample_name}_{i+1}"
                
                samples[sample_name] = {
                    'files': [file],
                    'metadata': {}
                }
                logger.info(f"Created sample {sample_name} from input file")
    
    # Method 2: Sample list file
    elif args.samples_file:
        logger.info(f"Reading samples from file: {args.samples_file}")
        
        # Use input handler to process sample list file
        samples = get_input_files(args, input_type="sequence")
    
    # Method 3: Metadata-driven
    elif args.metadata_file and args.seq_dir:
        logger.info(f"Reading samples from metadata: {args.metadata_file}")
        
        # Use input handler to process metadata and find sequence files
        samples = get_input_files(args, input_type="sequence")
    
    # Method 4: KneadData output directory
    elif args.input_dir:
        # Find KneadData output files in the input directory
        logger.info(f"Searching for KneadData output files in {args.input_dir}")
        kneaddata_files_by_sample = find_kneaddata_files_in_dir(args.input_dir)
        
        # Convert to standard sample format
        for sample_id, files in kneaddata_files_by_sample.items():
            samples[sample_id] = {
                'kneaddata_files': files,
                'files': [],
                'metadata': {}
            }
    
    if not samples:
        logger.error("No valid input samples found")
        return 1
    
    logger.info(f"Prepared {len(samples)} samples for HUMAnN3 processing")
    
    # Process HUMAnN3 options
    humann4_options = {}
    
    if args.bypass_prescreen:
        humann4_options["bypass-prescreen"] = True
    if args.bypass_nucleotide_index:
        humann4_options["bypass-nucleotide-index"] = True
    if args.bypass_translated_search:
        humann4_options["bypass-translated-search"] = True
    
    # Add any additional options
    if hasattr(args, 'humann4_options') and args.humann4_options:
        for option in args.humann4_options:
            if '=' in option:
                key, value = option.split('=', 1)
                humann4_options[key] = value
            else:
                humann4_options[option] = True
    
    # Run HUMAnN3
    if args.use_parallel:
        # Parallel processing
        logger.info("Using parallel processing for HUMAnN3")
        results = run_humann4_parallel(
            samples=samples,
            output_dir=args.output_dir,
            nucleotide_db=args.nucleotide_db,
            protein_db=args.protein_db,
            threads_per_sample=args.threads,
            max_parallel=args.max_parallel,
            options=humann4_options,
            paired=args.paired
        )
    else:
        # Sequential processing
        logger.info("Processing samples sequentially")
        results = {}
        
        for sample_id, sample_info in samples.items():
            # Prepare input file
            input_file = None
            
            # Check if we have KneadData files or direct input files
            if 'kneaddata_files' in sample_info and sample_info['kneaddata_files']:
                input_file = prepare_humann4_input(
                    sample_info['kneaddata_files'], 
                    sample_id, 
                    args.output_dir, 
                    args.paired
                )
            elif 'files' in sample_info and sample_info['files']:
                # If single input file, use it directly; if multiple, concatenate if paired
                if len(sample_info['files']) == 1:
                    input_file = sample_info['files'][0]
                elif args.paired and len(sample_info['files']) >= 2:
                    input_file = prepare_humann4_input(
                        sample_info['files'],
                        sample_id,
                        args.output_dir,
                        args.paired
                    )
                else:
                    logger.warning(f"Sample {sample_id} has {len(sample_info['files'])} files but paired={args.paired}. Using first file.")
                    input_file = sample_info['files'][0]
            
            if not input_file:
                logger.warning(f"Skipping sample {sample_id}: no valid input file")
                continue
                
            # Process the sample with HUMAnN3
            output_files = process_sample_humann4(
                sample_id=sample_id,
                input_file=input_file,
                output_dir=args.output_dir,
                nucleotide_db=args.nucleotide_db,
                protein_db=args.protein_db,
                threads=args.threads,
                options=humann4_options
            )
            
            if output_files:
                results[sample_id] = output_files
    
    # Organize outputs if requested
    if args.organize_outputs and results:
        logger.info("Organizing output files into type-specific directories")
        output_dirs = organize_output_files(results, args.output_dir)
        
        # Log where to find outputs
        for output_type, dir_path in output_dirs.items():
            count = sum(1 for sample in results.values() if sample.get(output_type))
            if count > 0:
                logger.info(f"{output_type.capitalize()} files ({count}): {dir_path}")
    
    # Log results summary
    if results:
        successful_samples = len(results)
        
        logger.info(f"HUMAnN3 processing completed successfully for {successful_samples}/{len(samples)} samples")
        
        # Count outputs by type
        output_counts = {
            'pathabundance': sum(1 for sample in results.values() if sample.get('pathabundance')),
            'genefamilies': sum(1 for sample in results.values() if sample.get('genefamilies')),
            'pathcoverage': sum(1 for sample in results.values() if sample.get('pathcoverage')),
            'metaphlan': sum(1 for sample in results.values() if sample.get('metaphlan'))
        }
        
        for output_type, count in output_counts.items():
            logger.info(f"  {output_type.capitalize()} files: {count}")
        
        logger.info(f"All output files are in: {os.path.abspath(args.output_dir)}")
    else:
        logger.error("HUMAnN3 processing failed for all samples")
        return 1
    
    # Print elapsed time
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    hours, minutes = divmod(minutes, 60)
    
    logger.info(f"Total processing time: {int(hours)}h {int(minutes)}m {int(seconds)}s")
    
    # Print next steps
    logger.info("\nNext Steps:")
    logger.info("  To join and normalize these files, use:")
    logger.info("  humann4-tools join --input-dir PathwayAbundance --output-dir joined_output")
    logger.info("  humann4-tools join --input-dir GeneFamilies --output-dir joined_output")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())