#!/usr/bin/env python3
"""
HUMAnN3 Tools KneadData Module

This module handles quality control and host depletion using KneadData.
It supports three input methods:
1. Direct input files
2. Sample list file
3. Metadata-driven input

Examples:
  # Basic usage with direct input files:
  humann4-kneaddata --input-files sample1_R1.fastq.gz sample1_R2.fastq.gz --paired --threads 8 --output-dir ./kneaddata_output

  # Using a sample list file:
  humann4-kneaddata --samples-file samples.txt --reference-dbs human --output-dir ./kneaddata_output

  # Using metadata:
  humann4-kneaddata --metadata-file metadata.csv --seq-dir /path/to/sequences --r1-suffix _R1.fastq.gz --r2-suffix _R2.fastq.gz --paired --output-dir ./kneaddata_output
  
  # With additional KneadData options:
  humann4-kneaddata --input-files sample.fastq.gz --reference-dbs human --output-dir ./kneaddata_output --kneaddata-options trimmomatic-options=SLIDINGWINDOW:4:20,MINLEN:50
"""

import os
import sys
import time
import argparse
import logging
import subprocess
import multiprocessing
from typing import Dict, List, Optional, Tuple

# Set up path for imports - add parent directory to path if module imports fail
try:
    from src.humann4_tools.utils.input_handler import get_input_files
    from src.humann4_tools.utils.cmd_utils import run_cmd
    try:
        from src.humann4_tools.utils.resource_utils import track_peak_memory
        TRACK_MEMORY = True
    except ImportError:
        TRACK_MEMORY = False
except ImportError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
    try:
        from src.humann4_tools.utils.input_handler import get_input_files
        from src.humann4_tools.utils.cmd_utils import run_cmd
        try:
            from src.humann4_tools.utils.resource_utils import track_peak_memory
            TRACK_MEMORY = True
        except ImportError:
            TRACK_MEMORY = False
    except ImportError:
        print("ERROR: Failed to import required modules. Please ensure the humann4_tools package is properly installed.")
        sys.exit(1)

# Create a no-op decorator to use if track_peak_memory is not available
if not TRACK_MEMORY:
    def track_peak_memory(func):
        return func

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

def check_kneaddata_installation() -> Tuple[bool, str]:
    """Check if KneadData is installed and available."""
    try:
        result = subprocess.run(["kneaddata", "--version"], 
                               capture_output=True, text=True, check=False)
        if result.returncode == 0:
            return True, result.stdout.strip()
        return False, "KneadData command exists but returned an error"
    except FileNotFoundError:
        return False, "KneadData not found in PATH"

def process_sample_kneaddata(sample_id: str, 
                            input_files: List[str], 
                            output_dir: str,
                            reference_dbs: List[str],
                            threads: int = 1,
                            paired: bool = False,
                            options: Optional[Dict] = None) -> List[str]:
    """
    Process a single sample with KneadData.
    
    Args:
        sample_id: Sample identifier
        input_files: List of input file paths (1 for single-end, 2 for paired-end)
        output_dir: Directory for output files
        reference_dbs: List of reference database paths
        threads: Number of threads to use
        paired: Whether input is paired-end
        options: Additional KneadData options
        
    Returns:
        List of output file paths
    """
    # Create sample output directory
    sample_outdir = os.path.join(output_dir, sample_id)
    os.makedirs(sample_outdir, exist_ok=True)
    
    # Build KneadData command
    cmd = ["kneaddata"]
    
    if paired and len(input_files) >= 2:
        cmd.extend(["--input1", input_files[0], "--input2", input_files[1]])
        logger.info(f"Running KneadData on paired files for sample {sample_id}")
    elif len(input_files) >= 1:
        cmd.extend(["--input", input_files[0]])
        logger.info(f"Running KneadData on single-end file for sample {sample_id}")
    else:
        logger.error(f"No input files provided for sample {sample_id}")
        return []
    
    # Add output directory
    cmd.extend(["--output", sample_outdir])
    
    # Add reference databases
    for db in reference_dbs:
        cmd.extend(["--reference-db", db])
    
    # Add threads
    cmd.extend(["--threads", str(threads)])
    
    # Add additional options
    if options:
        for key, value in options.items():
            if key == 'paired':
                continue  # Skip as it's handled differently
                
            if value is True:
                cmd.append(f"--{key}")
            elif value is not None and value != "":
                cmd.extend([f"--{key}", str(value)])
    
    # Run KneadData
    logger.info(f"Running KneadData command: {' '.join(cmd)}")
    success = run_cmd(cmd, exit_on_error=False)
    
    if not success:
        logger.error(f"KneadData failed for sample {sample_id}")
        return []
    
    # Find output files
    output_files = []
    for file in os.listdir(sample_outdir):
        if file.endswith(".fastq") and (
            "paired" in file or "kneaddata" in file or "cleaned" in file
        ):
            output_files.append(os.path.join(sample_outdir, file))
    
    if output_files:
        logger.info(f"KneadData completed for sample {sample_id} with {len(output_files)} output files")
    else:
        logger.warning(f"No output files found for sample {sample_id}")
    
    return output_files

def run_kneaddata_parallel(samples: Dict[str, Dict],
                        output_dir: str,
                        reference_dbs: List[str],
                        threads_per_sample: int = 1,
                        max_parallel: Optional[int] = None,
                        paired: bool = False,
                        options: Optional[Dict] = None) -> Dict[str, List[str]]:
    """
    Run KneadData on multiple samples in parallel.
    
    Args:
        samples: Dictionary mapping sample IDs to sample information
        output_dir: Base directory for outputs
        reference_dbs: List of reference database paths
        threads_per_sample: Number of threads per sample
        max_parallel: Maximum number of parallel samples
        paired: Whether to process as paired-end
        options: Dictionary of additional KneadData options
        
    Returns:
        Dictionary mapping sample IDs to lists of output file paths
    """
    from concurrent.futures import ProcessPoolExecutor
    
    # Set default max_parallel based on CPU count if not specified
    if max_parallel is None:
        available_cpus = multiprocessing.cpu_count()
        max_parallel = max(1, available_cpus // threads_per_sample)
    
    logger.info(f"Running KneadData in parallel: {len(samples)} samples, " 
                f"{max_parallel} parallel processes, {threads_per_sample} threads per sample")
    
    results = {}
    
    # Define a wrapper for parallel processing
    def process_sample_wrapper(sample_id, sample_info):
        if not sample_info['files']:
            logger.warning(f"Skipping sample {sample_id}: no input files")
            return sample_id, []
        
        output_files = process_sample_kneaddata(
            sample_id=sample_id, 
            input_files=sample_info['files'],
            output_dir=output_dir,
            reference_dbs=reference_dbs,
            threads=threads_per_sample,
            paired=paired,
            options=options
        )
        return sample_id, output_files
    
    # Process samples in parallel
    with ProcessPoolExecutor(max_workers=max_parallel) as executor:
        # Submit jobs
        futures = {
            executor.submit(process_sample_wrapper, sample_id, sample_info): sample_id 
            for sample_id, sample_info in samples.items()
        }
        
        # Collect results
        for future in futures:
            sample_id = futures[future]
            try:
                result_id, output_files = future.result()
                if output_files:
                    results[result_id] = output_files
                    logger.info(f"Successfully processed sample {result_id}")
                else:
                    logger.error(f"Failed to process sample {result_id}")
            except Exception as e:
                logger.error(f"Error processing sample {sample_id}: {str(e)}")
    
    return results

def parse_args():
    """Parse command line arguments for the KneadData module."""
    parser = argparse.ArgumentParser(
        description="Process FASTQ files with KneadData for quality control and host depletion"
    )
    
    # Input options group (three methods supported)
    input_group = parser.add_argument_group("Input Options (choose one)")
    # Method 1: Direct input files
    input_group.add_argument("--input-files", nargs="+", 
                           help="Input FASTQ file(s) for KneadData")
    # Method 2: Sample list file
    input_group.add_argument("--samples-file", 
                           help="Tab-delimited file with sample IDs and file paths")
    # Method 3: Metadata-driven
    input_group.add_argument("--metadata-file", 
                           help="CSV file with sample metadata")
    input_group.add_argument("--seq-dir", 
                           help="Directory containing sequence files (for metadata-driven)")
    input_group.add_argument("--sample-col", 
                           help="Column name for sample IDs in metadata")
    input_group.add_argument("--r1-col", 
                           help="Column name for R1 file paths in metadata")
    input_group.add_argument("--r2-col", 
                           help="Column name for R2 file paths in metadata")
    input_group.add_argument("--file-pattern", 
                           help="Pattern for finding files (e.g., {sample}_S*_R*.fastq.gz)")
    input_group.add_argument("--r1-suffix", 
                           help="Suffix for R1 files (e.g., _R1.fastq.gz)")
    input_group.add_argument("--r2-suffix", 
                           help="Suffix for R2 files (e.g., _R2.fastq.gz)")
    
    # Required arguments
    parser.add_argument("--reference-dbs", nargs="+", required=True,
                      help="Path(s) to KneadData reference database(s)")
    parser.add_argument("--output-dir", default="./kneaddata_output",
                      help="Directory for KneadData output")
    
    # KneadData options
    parser.add_argument("--paired", action="store_true",
                      help="Input files are paired-end reads")
    parser.add_argument("--decontaminate-pairs", default="strict", 
                      choices=["strict", "lenient", "unpaired"],
                      help="Method for decontaminating paired-end reads (default: strict)")
    
    # Performance options
    parser.add_argument("--threads", type=int, default=1,
                      help="Number of threads to use per sample")
    parser.add_argument("--use-parallel", action="store_true",
                      help="Process multiple samples in parallel")
    parser.add_argument("--max-parallel", type=int, default=None,
                      help="Maximum number of samples to process in parallel")
    
    # Logging options
    parser.add_argument("--log-file", 
                      help="Path to log file")
    parser.add_argument("--log-level", default="INFO",
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                      help="Logging level")
    
    # Additional KneadData options
    parser.add_argument("--kneaddata-options", nargs="+",
                      help="Additional options to pass to KneadData (format: key=value)")
    
    return parser.parse_args()

@track_peak_memory
def main():
    """Main function to run KneadData processing."""
    # Parse arguments
    args = parse_args()
    
    # Setup logging
    log_level = getattr(logging, args.log_level.upper())
    setup_logger(args.log_file, log_level)
    
    start_time = time.time()
    logger.info("Starting HUMAnN3 Tools KneadData Module")
    
    # Check KneadData installation
    kneaddata_ok, kneaddata_version = check_kneaddata_installation()
    if not kneaddata_ok:
        logger.error(f"KneadData not properly installed: {kneaddata_version}")
        return 1
    logger.info(f"Using KneadData version: {kneaddata_version}")
    
    # Get input files using the input handler
    samples = get_input_files(args, input_type="sequence")
    
    if not samples:
        logger.error("No valid input samples found")
        return 1
    
    logger.info(f"Prepared {len(samples)} samples for KneadData processing")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Process KneadData options
    kneaddata_options = {}
    if args.paired:
        kneaddata_options["decontaminate-pairs"] = args.decontaminate_pairs
    
    # Add any additional options
    if args.kneaddata_options:
        for option in args.kneaddata_options:
            if '=' in option:
                key, value = option.split('=', 1)
                kneaddata_options[key] = value
            else:
                kneaddata_options[option] = True
    
    # Run KneadData
    if args.use_parallel:
        # Parallel processing
        logger.info("Using parallel processing for KneadData")
        results = run_kneaddata_parallel(
            samples=samples,
            output_dir=args.output_dir,
            reference_dbs=args.reference_dbs,
            threads_per_sample=args.threads,
            max_parallel=args.max_parallel,
            paired=args.paired,
            options=kneaddata_options
        )
    else:
        # Sequential processing
        logger.info("Processing samples sequentially")
        results = {}
        for sample_id, sample_info in samples.items():
            if not sample_info['files']:
                logger.warning(f"Skipping sample {sample_id}: no input files")
                continue
                
            output_files = process_sample_kneaddata(
                sample_id=sample_id,
                input_files=sample_info['files'],
                output_dir=args.output_dir,
                reference_dbs=args.reference_dbs,
                threads=args.threads,
                paired=args.paired,
                options=kneaddata_options
            )
            
            if output_files:
                results[sample_id] = output_files
    
    # Log results summary
    if results:
        successful_samples = len(results)
        total_files = sum(len(files) for files in results.values())
        
        logger.info(f"KneadData processing completed successfully for {successful_samples}/{len(samples)} samples")
        logger.info(f"Generated a total of {total_files} output files")
        
        # Print locations of some example output files
        logger.info("Example output files:")
        for sample_id, files in list(results.items())[:2]:  # Show first 2 samples
            for i, file_path in enumerate(files[:2]):  # Show first 2 files per sample
                logger.info(f"  {sample_id} ({i+1}): {os.path.basename(file_path)}")
        
        logger.info(f"All output files are in: {os.path.abspath(args.output_dir)}")
    else:
        logger.error("KneadData processing failed for all samples")
        return 1
    
    # Print elapsed time
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    hours, minutes = divmod(minutes, 60)
    
    logger.info(f"Total processing time: {int(hours)}h {int(minutes)}m {int(seconds)}s")
    
    # Print next steps
    logger.info("\nNext Steps:")
    logger.info("  To run HUMAnN3 on these files, use:")
    logger.info("  humann4-tools humann4 --input-dir kneaddata_output --output-dir humann4_output")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
