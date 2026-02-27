#!/usr/bin/env python3
"""
HUMAnN3 Tools Preprocessing Script

This standalone script handles the preprocessing steps for HUMAnN3 analysis:
1. Quality control with KneadData (optional)
2. Running HUMAnN3 on preprocessed sequences

It can be used independently from the main humann4_tools pipeline.

Examples:
  # Basic usage with paired-end reads:
  python humann4_tools_preprocessing.py --input-fastq sample1_R1.fastq.gz sample1_R2.fastq.gz --paired --threads 8 --output-dir ./my_output
  
  # Running with multiple samples:
  python humann4_tools_preprocessing.py --input-fastq sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz --paired --threads 8
  
  # Specifying database paths:
  python humann4_tools_preprocessing.py --input-fastq sample1.fastq.gz --kneaddata-dbs /path/to/kneaddata_db --humann4-nucleotide-db /path/to/chocophlan --humann4-protein-db /path/to/uniref
  
  # Skipping KneadData and using existing files:
  python humann4_tools_preprocessing.py --input-fastq sample1_R1.fastq.gz sample1_R2.fastq.gz --skip-kneaddata --kneaddata-output-files /path/to/kneaddata/sample1_paired_1.fastq /path/to/kneaddata/sample1_paired_2.fastq

Dependencies:
  - KneadData (v0.7.0+)
  - HUMAnN3 (v3.0.0+)
  - Python 3.6+
"""

import os
import sys
import time
import argparse
import logging
import glob
from typing import List, Dict, Optional, Union, Tuple

# Set up logging
logger = logging.getLogger('humann4_preprocessing')
logger.setLevel(logging.INFO)

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

def log_print(message, level='info'):
    """Print to console and also log with the logger."""
    print(message)
    if level.lower() == 'debug':
        logger.debug(message)
    elif level.lower() == 'warning':
        logger.warning(message)
    elif level.lower() == 'error':
        logger.error(message)
    elif level.lower() == 'critical':
        logger.critical(message)
    else:
        logger.info(message)

def run_cmd(cmd, exit_on_error=True, verbose=True):
    """
    Run a shell command with subprocess.
    
    Args:
        cmd (list): Command to run as a list of strings
        exit_on_error (bool): Whether to exit the program if the command fails
        verbose (bool): Whether to print/log the command being run
        
    Returns:
        bool: True if command executed successfully, False otherwise
    """
    import subprocess
    
    if verbose:
        logger.info(f"Running: {' '.join(cmd)}")
    
    # Additional check for 'cp' commands
    if cmd[0] == "cp" and len(cmd) >= 3:
        src = cmd[1]
        if not os.path.exists(src):
            logger.error(f"ERROR: Source file does not exist: {src}")
            if exit_on_error:
                sys.exit(1)
            return False
        
        dst_dir = os.path.dirname(cmd[2])
        if not os.path.exists(dst_dir):
            logger.info(f"Creating directory: {dst_dir}")
            os.makedirs(dst_dir, exist_ok=True)
    
    try:
        process = subprocess.run(cmd, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout = process.stdout.decode('utf-8', errors='replace')
        stderr = process.stderr.decode('utf-8', errors='replace')
        if stdout.strip():
            logger.debug(f"Command stdout: {stdout}")
        if stderr.strip():
            logger.debug(f"Command stderr: {stderr}")
        return True
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr.decode('utf-8', errors='replace')
        logger.error(f"ERROR: Command failed with exit code {e.returncode}")
        logger.error(f"Error message: {error_msg}")
        logger.error(f"Failed command: {' '.join(cmd)}")
        if exit_on_error:
            sys.exit(1)
        return False

def check_kneaddata_installation():
    """Check if KneadData is installed and available."""
    import subprocess
    try:
        result = subprocess.run(["kneaddata", "--version"], 
                               capture_output=True, text=True, check=False)
        if result.returncode == 0:
            return True, result.stdout.strip()
        return False, "KneadData command exists but returned an error"
    except FileNotFoundError:
        return False, "KneadData not found in PATH"

def check_humann4_installation():
    """Check if HUMAnN3 is installed and available."""
    import subprocess
    try:
        result = subprocess.run(["humann", "--version"], 
                               capture_output=True, text=True, check=False)
        if result.returncode == 0:
            return True, result.stdout.strip()
        return False, "HUMAnN3 command exists but returned an error"
    except FileNotFoundError:
        return False, "HUMAnN3 not found in PATH"

def run_kneaddata(input_files, output_dir, threads=1, reference_dbs=None, 
                 paired=False, additional_options=None):
    """
    Run KneadData on input sequence files.
    
    Args:
        input_files: List of input FASTQ files
        output_dir: Directory for KneadData output
        threads: Number of threads to use
        reference_dbs: List of reference database paths
        paired: Whether input files are paired
        additional_options: Dict of additional KneadData options
        
    Returns:
        Dict mapping sample IDs to output FASTQ files
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Process samples individually to keep track of outputs by sample
    sample_outputs = {}
    
    # Handle paired vs single end
    if paired:
        # Ensure we have an even number of files
        if len(input_files) % 2 != 0:
            logger.error("Paired mode requires an even number of input files")
            sys.exit(1)
        
        # Process in pairs
        for i in range(0, len(input_files), 2):
            r1_file = input_files[i]
            r2_file = input_files[i+1]
            
            # Extract sample name from R1 file
            sample_name = os.path.basename(r1_file).split('_')[0]
            
            # Create sample-specific output directory
            sample_outdir = os.path.join(output_dir, sample_name)
            os.makedirs(sample_outdir, exist_ok=True)
            
            # Build command
            cmd = ["kneaddata", "--input1", r1_file, "--input2", r2_file, 
                   "--output", sample_outdir, "--threads", str(threads)]
            
            # Add reference database(s)
            if reference_dbs:
                if isinstance(reference_dbs, str):
                    cmd.extend(["--reference-db", reference_dbs])
                else:
                    for db in reference_dbs:
                        cmd.extend(["--reference-db", db])
            
            # Add additional options
            if additional_options:
                for key, value in additional_options.items():
                    if key == 'paired':
                        continue  # Skip the paired option as it's handled differently
                        
                    if value is True:
                        cmd.append(f"--{key}")
                    elif value is not None and value != "":
                        cmd.extend([f"--{key}", str(value)])
            
            # Run KneadData
            logger.info(f"Running KneadData for paired sample {sample_name}")
            success = run_cmd(cmd, exit_on_error=False)
            
            if not success:
                logger.error(f"KneadData failed for sample {sample_name}")
                continue
            
            # Find output files
            output_files = []
            for file in os.listdir(sample_outdir):
                if file.endswith(".fastq") and "paired" in file:
                    output_files.append(os.path.join(sample_outdir, file))
            
            if output_files:
                sample_outputs[sample_name] = output_files
                logger.info(f"KneadData completed for sample {sample_name} with {len(output_files)} output files")
            else:
                logger.warning(f"No valid output files found for sample {sample_name}")
    else:
        # Single-end processing
        for input_file in input_files:
            # Extract sample name
            sample_name = os.path.basename(input_file).split('.')[0]
            
            # Create sample-specific output directory
            sample_outdir = os.path.join(output_dir, sample_name)
            os.makedirs(sample_outdir, exist_ok=True)
            
            # Build command
            cmd = ["kneaddata", "--input", input_file, "--output", sample_outdir, 
                   "--threads", str(threads)]
            
            # Add reference database(s)
            if reference_dbs:
                if isinstance(reference_dbs, str):
                    cmd.extend(["--reference-db", reference_dbs])
                else:
                    for db in reference_dbs:
                        cmd.extend(["--reference-db", db])
            
            # Add additional options
            if additional_options:
                for key, value in additional_options.items():
                    if key == 'paired':
                        continue  # Skip the paired option
                        
                    if value is True:
                        cmd.append(f"--{key}")
                    elif value is not None and value != "":
                        cmd.extend([f"--{key}", str(value)])
            
            # Run KneadData
            logger.info(f"Running KneadData for single-end sample {sample_name}")
            success = run_cmd(cmd, exit_on_error=False)
            
            if not success:
                logger.error(f"KneadData failed for sample {sample_name}")
                continue
            
            # Find output files
            output_files = []
            for file in os.listdir(sample_outdir):
                if file.endswith(".fastq"):
                    output_files.append(os.path.join(sample_outdir, file))
            
            if output_files:
                sample_outputs[sample_name] = output_files
                logger.info(f"KneadData completed for sample {sample_name} with {len(output_files)} output files")
            else:
                logger.warning(f"No valid output files found for sample {sample_name}")
    
    return sample_outputs

def prepare_humann4_input(kneaddata_outputs, paired=False):
    """
    Prepare KneadData outputs for HUMAnN3 input.
    For paired data, this concatenates the paired reads into a single file.
    
    Args:
        kneaddata_outputs: Dict mapping sample IDs to lists of KneadData output files
        paired: Whether we processed paired-end data
        
    Returns:
        Dict mapping sample IDs to HUMAnN3 input files
    """
    humann4_inputs = {}
    
    for sample_id, files in kneaddata_outputs.items():
        if paired:
            # Find and sort paired files
            r1_files = [f for f in files if "_paired_1.fastq" in f or "_kneaddata_paired_1.fastq" in f]
            r2_files = [f for f in files if "_paired_2.fastq" in f or "_kneaddata_paired_2.fastq" in f]
            
            if not r1_files or not r2_files:
                logger.warning(f"Missing paired files for sample {sample_id}")
                continue
            
            # Sort to ensure consistent ordering
            r1_files.sort()
            r2_files.sort()
            
            # Take the first pair
            r1_file = r1_files[0]
            r2_file = r2_files[0]
            
            # Create concatenated file
            out_dir = os.path.dirname(r1_file)
            concat_file = os.path.join(out_dir, f"{sample_id}_paired_concat.fastq")
            
            try:
                with open(concat_file, 'w') as outfile:
                    logger.info(f"Concatenating files for {sample_id}: {os.path.basename(r1_file)} + {os.path.basename(r2_file)}")
                    for input_file in [r1_file, r2_file]:
                        with open(input_file, 'r') as infile:
                            outfile.write(infile.read())
                
                if os.path.getsize(concat_file) > 0:
                    humann4_inputs[sample_id] = concat_file
                    logger.info(f"Created concatenated file for HUMAnN3: {concat_file}")
                else:
                    logger.error(f"Created empty concatenated file for {sample_id}")
            except Exception as e:
                logger.error(f"Error concatenating files for {sample_id}: {str(e)}")
        else:
            # For single-end, just use the first file
            if files:
                # Look for files containing "_kneaddata." in their name
                kneaddata_output_files = [f for f in files if "_kneaddata." in f]
                
                if kneaddata_output_files:
                    # Use the first kneaddata output file
                    humann4_inputs[sample_id] = kneaddata_output_files[0]
                else:
                    # If no kneaddata-specific files, use the first available file
                    humann4_inputs[sample_id] = files[0]
                
                logger.info(f"Using file for HUMAnN3: {humann4_inputs[sample_id]}")
    
    return humann4_inputs

def find_existing_kneaddata_files(kneaddata_output_dir, kneaddata_output_pattern=None, 
                                 kneaddata_output_files=None, input_files=None):
    """
    Find existing KneadData output files.
    
    Args:
        kneaddata_output_dir: Directory containing KneadData output
        kneaddata_output_pattern: Pattern to match KneadData output files
        kneaddata_output_files: Explicitly provided KneadData output files
        input_files: Original input files (to extract sample names)
        
    Returns:
        Dict mapping sample IDs to lists of KneadData output files
    """
    sample_outputs = {}
    
    # Case 1: Explicitly provided files
    if kneaddata_output_files:
        logger.info(f"Using {len(kneaddata_output_files)} provided KneadData output files")
        
        for file_path in kneaddata_output_files:
            if os.path.exists(file_path) and os.path.isfile(file_path):
                filename = os.path.basename(file_path)
                
                # Try to extract sample name from filename
                sample_name = None
                if "kneaddata_paired_1.fastq" in filename:
                    sample_name = filename.split("kneaddata_paired_1.fastq")[0].rstrip("_")
                elif "kneaddata_paired_2.fastq" in filename:
                    sample_name = filename.split("kneaddata_paired_2.fastq")[0].rstrip("_")
                elif "paired_1.fastq" in filename:
                    sample_name = filename.split("paired_1.fastq")[0].rstrip("_")
                elif "paired_2.fastq" in filename:
                    sample_name = filename.split("paired_2.fastq")[0].rstrip("_")
                
                # If we could extract sample name, add file to sample's group
                if sample_name:
                    if sample_name not in sample_outputs:
                        sample_outputs[sample_name] = []
                    sample_outputs[sample_name].append(file_path)
            else:
                logger.warning(f"KneadData output file not found: {file_path}")
    
    # Case 2: Pattern matching
    elif kneaddata_output_pattern:
        logger.info(f"Looking for KneadData output files with pattern: {kneaddata_output_pattern}")
        
        # Extract potential sample names from input files if available
        sample_names = set()
        if input_files:
            for input_file in input_files:
                sample_name = os.path.basename(input_file).split("_")[0]
                sample_names.add(sample_name)
        
        # If the pattern includes {sample}, replace for each sample name
        if '{sample}' in kneaddata_output_pattern:
            for sample_name in sample_names:
                pattern = kneaddata_output_pattern.format(sample=sample_name)
                files = glob.glob(pattern)
                
                if files:
                    sample_outputs[sample_name] = files
                    logger.info(f"Found {len(files)} files for sample {sample_name}")
        else:
            # Use pattern directly and try to extract sample names
            files = glob.glob(kneaddata_output_pattern)
            logger.info(f"Found {len(files)} files matching pattern {kneaddata_output_pattern}")
            
            for file_path in files:
                filename = os.path.basename(file_path)
                
                # Try to extract sample name from filename
                sample_name = None
                if "kneaddata_paired_1.fastq" in filename:
                    sample_name = filename.split("kneaddata_paired_1.fastq")[0].rstrip("_")
                elif "kneaddata_paired_2.fastq" in filename:
                    sample_name = filename.split("kneaddata_paired_2.fastq")[0].rstrip("_")
                elif "paired_1.fastq" in filename:
                    sample_name = filename.split("paired_1.fastq")[0].rstrip("_")
                elif "paired_2.fastq" in filename:
                    sample_name = filename.split("paired_2.fastq")[0].rstrip("_")
                
                # If we could extract sample name, add file to sample's group
                if sample_name:
                    if sample_name not in sample_outputs:
                        sample_outputs[sample_name] = []
                    sample_outputs[sample_name].append(file_path)
    
    # Case 3: Search in output directory
    else:
        logger.info(f"Looking for KneadData output files in {kneaddata_output_dir}")
        
        if os.path.isdir(kneaddata_output_dir):
            # Look in subdirectories (samples are often in their own directories)
            for root, dirs, files in os.walk(kneaddata_output_dir):
                for filename in files:
                    if filename.endswith(".fastq") and ("paired_1" in filename or "paired_2" in filename):
                        file_path = os.path.join(root, filename)
                        
                        # Try to extract sample name from filename or parent directory
                        sample_name = None
                        
                        # First try from filename
                        if "kneaddata_paired_1.fastq" in filename:
                            sample_name = filename.split("kneaddata_paired_1.fastq")[0].rstrip("_")
                        elif "kneaddata_paired_2.fastq" in filename:
                            sample_name = filename.split("kneaddata_paired_2.fastq")[0].rstrip("_")
                        elif "paired_1.fastq" in filename:
                            sample_name = filename.split("paired_1.fastq")[0].rstrip("_")
                        elif "paired_2.fastq" in filename:
                            sample_name = filename.split("paired_2.fastq")[0].rstrip("_")
                        
                        # If not found, try from parent directory
                        if not sample_name:
                            parent_dir = os.path.basename(root)
                            if parent_dir != os.path.basename(kneaddata_output_dir):
                                sample_name = parent_dir
                        
                        # If we could extract sample name, add file to sample's group
                        if sample_name:
                            if sample_name not in sample_outputs:
                                sample_outputs[sample_name] = []
                            sample_outputs[sample_name].append(file_path)
        else:
            logger.warning(f"KneadData output directory not found: {kneaddata_output_dir}")
    
    return sample_outputs

def run_humann4(input_files, output_dir, threads=1, nucleotide_db=None, protein_db=None, 
               additional_options=None):
    """
    Run HUMAnN3 on input sequence files.
    
    Args:
        input_files: Dict mapping sample IDs to input FASTQ files
        output_dir: Directory for HUMAnN3 output
        threads: Number of threads to use
        nucleotide_db: Path to nucleotide database
        protein_db: Path to protein database
        additional_options: Dict of additional HUMAnN3 options
        
    Returns:
        Dict mapping sample IDs to dicts of output file paths by type
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Create subdirectories for different output types
    pathabundance_dir = os.path.join(output_dir, "PathwayAbundance")
    genefamilies_dir = os.path.join(output_dir, "GeneFamilies")
    pathcoverage_dir = os.path.join(output_dir, "PathwayCoverage")
    metaphlan_dir = os.path.join(output_dir, "MetaphlanFiles")
    
    os.makedirs(pathabundance_dir, exist_ok=True)
    os.makedirs(genefamilies_dir, exist_ok=True)
    os.makedirs(pathcoverage_dir, exist_ok=True)
    os.makedirs(metaphlan_dir, exist_ok=True)
    
    # Process each sample
    humann4_results = {}
    
    for sample_id, input_file in input_files.items():
        # Create sample output directory
        sample_outdir = os.path.join(output_dir, sample_id)
        os.makedirs(sample_outdir, exist_ok=True)
        
        # Build command
        cmd = ["humann", "--input", input_file, "--output", sample_outdir]
        
        # Add threads
        cmd.extend(["--threads", str(threads)])
        
        # Add database paths if provided
        if nucleotide_db:
            cmd.extend(["--nucleotide-database", nucleotide_db])
        if protein_db:
            cmd.extend(["--protein-database", protein_db])
        
        # Add additional options
        if additional_options:
            for key, value in additional_options.items():
                if value is True:
                    cmd.append(f"--{key}")
                elif value is not None and value != "":
                    cmd.extend([f"--{key}", str(value)])
        
        # Run HUMAnN3
        logger.info(f"Running HUMAnN3 for sample {sample_id}")
        success = run_cmd(cmd, exit_on_error=False)
        
        if not success:
            logger.error(f"HUMAnN3 failed for sample {sample_id}")
            continue
        
        # Find output files
        sample_outputs = {
            'genefamilies': None,
            'pathabundance': None,
            'pathcoverage': None,
            'metaphlan': None
        }
        
        # Look for output files in sample output directory
        for root, _, files in os.walk(sample_outdir):
            for filename in files:
                filepath = os.path.join(root, filename)
                
                if not filename.endswith(".tsv"):
                    continue
                
                if "genefamilies" in filename.lower():
                    sample_outputs['genefamilies'] = filepath
                    # Copy to genefamilies directory
                    dest_path = os.path.join(genefamilies_dir, filename)
                    run_cmd(["cp", filepath, dest_path], exit_on_error=False)
                
                elif "pathabundance" in filename.lower():
                    sample_outputs['pathabundance'] = filepath
                    # Copy to pathabundance directory
                    dest_path = os.path.join(pathabundance_dir, filename)
                    run_cmd(["cp", filepath, dest_path], exit_on_error=False)
                
                elif "pathcoverage" in filename.lower():
                    sample_outputs['pathcoverage'] = filepath
                    # Copy to pathcoverage directory
                    dest_path = os.path.join(pathcoverage_dir, filename)
                    run_cmd(["cp", filepath, dest_path], exit_on_error=False)
                
                elif "metaphlan_bugs_list" in filename.lower():
                    sample_outputs['metaphlan'] = filepath
                    # Copy to metaphlan directory
                    dest_path = os.path.join(metaphlan_dir, filename)
                    run_cmd(["cp", filepath, dest_path], exit_on_error=False)
        
        # Log which files were found
        found_files = [k for k, v in sample_outputs.items() if v is not None]
        missing_files = [k for k, v in sample_outputs.items() if v is None]
        
        logger.info(f"Found output files for {sample_id}: {', '.join(found_files)}")
        if missing_files:
            logger.warning(f"Missing output files for {sample_id}: {', '.join(missing_files)}")
        
        humann4_results[sample_id] = sample_outputs
    
    return humann4_results

def print_sample_summary_table(results):
    """
    Print a formatted summary table of processed samples.
    
    Args:
        results: Dict with preprocessing results
    """
    if not results or not results.get('humann4_results'):
        log_print("No results to display", level="warning")
        return
    
    # Calculate max width for sample names
    sample_names = list(results['humann4_results'].keys())
    max_name_width = max(len(name) for name in sample_names) if sample_names else 10
    max_name_width = max(max_name_width, 10)  # at least 10 characters wide
    
    # Print table header
    header = f"| {'Sample'.ljust(max_name_width)} | {'Pathway'.center(10)} | {'Gene Fam.'.center(10)} | {'PathCov.'.center(10)} | {'MetaPhlAn'.center(10)} |"
    separator = f"|-{'-' * max_name_width}-|-{'-' * 10}-|-{'-' * 10}-|-{'-' * 10}-|-{'-' * 10}-|"
    
    log_print("\nSample Results Summary:", level="info")
    log_print(separator, level="info")
    log_print(header, level="info")
    log_print(separator, level="info")
    
    # Print each sample's results
    for sample_id, files in results['humann4_results'].items():
        path_status = "✓" if files.get('pathabundance') else "-"
        gene_status = "✓" if files.get('genefamilies') else "-"
        cov_status = "✓" if files.get('pathcoverage') else "-"
        meta_status = "✓" if files.get('metaphlan') else "-"
        
        row = f"| {sample_id.ljust(max_name_width)} | {path_status.center(10)} | {gene_status.center(10)} | {cov_status.center(10)} | {meta_status.center(10)} |"
        log_print(row, level="info")
    
    log_print(separator, level="info")
    log_print("", level="info")  # Add a blank line after the table

def run_preprocessing_pipeline(input_files, output_dir, threads=1, 
                              kneaddata_dbs=None, nucleotide_db=None, protein_db=None,
                              paired=False, kneaddata_options=None, humann4_options=None,
                              skip_kneaddata=False, kneaddata_output_files=None,
                              kneaddata_output_pattern=None, kneaddata_output_dir=None):
    """
    Run the full preprocessing pipeline: KneadData → HUMAnN3.
    
    Args:
        input_files: List of input FASTQ files
        output_dir: Base directory for outputs
        threads: Number of threads per sample
        kneaddata_dbs: Path(s) to KneadData reference database(s)
        nucleotide_db: Path to HUMAnN3 nucleotide database
        protein_db: Path to HUMAnN3 protein database
        paired: Whether input files are paired
        kneaddata_options: Dict of additional KneadData options
        humann4_options: Dict of additional HUMAnN3 options
        skip_kneaddata: Whether to skip KneadData processing
        kneaddata_output_files: List of existing KneadData output files to use
        kneaddata_output_pattern: Pattern to find KneadData output files
        kneaddata_output_dir: Custom directory for KneadData outputs
        
    Returns:
        Dict with preprocessing results
    """
    # Create the main output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Set default output directories if not provided
    if kneaddata_output_dir is None:
        kneaddata_output_dir = os.path.join(output_dir, "kneaddata_output")
    
    humann4_output_dir = os.path.join(output_dir, "humann4_output")
    
    # Create output directories
    os.makedirs(kneaddata_output_dir, exist_ok=True)
    os.makedirs(humann4_output_dir, exist_ok=True)
    
    # 1. KneadData processing (or find existing files)
    kneaddata_results = {}
    
    if skip_kneaddata:
        logger.info("Skipping KneadData processing as requested")
        # Find existing KneadData output files
        kneaddata_results = find_existing_kneaddata_files(
            kneaddata_output_dir,
            kneaddata_output_pattern,
            kneaddata_output_files,
            input_files
        )
        
        if not kneaddata_results:
            logger.error("No KneadData output files found. Cannot proceed.")
            sys.exit(1)
        
        logger.info(f"Found KneadData output files for {len(kneaddata_results)} samples")
    else:
        # Run KneadData
        kneaddata_results = run_kneaddata(
            input_files=input_files,
            output_dir=kneaddata_output_dir,
            threads=threads,
            reference_dbs=kneaddata_dbs,
            paired=paired,
            additional_options=kneaddata_options
        )
        
        if not kneaddata_results:
            logger.error("KneadData processing failed. Cannot proceed.")
            sys.exit(1)
        
        logger.info(f"KneadData processing completed for {len(kneaddata_results)} samples")
    
    # 2. Prepare files for HUMAnN3
    humann4_input_files = prepare_humann4_input(kneaddata_results, paired=paired)
    
    if not humann4_input_files:
        logger.error("Failed to prepare input files for HUMAnN3. Cannot proceed.")
        sys.exit(1)
    
    logger.info(f"Prepared {len(humann4_input_files)} input files for HUMAnN3")
    
    # 3. Run HUMAnN3
    humann4_results = run_humann4(
        input_files=humann4_input_files,
        output_dir=humann4_output_dir,
        threads=threads,
        nucleotide_db=nucleotide_db,
        protein_db=protein_db,
        additional_options=humann4_options
    )
    
    if not humann4_results:
        logger.error("HUMAnN3 processing failed.")
        return {
            "kneaddata_results": kneaddata_results,
            "humann4_results": {}
        }
    
    logger.info(f"HUMAnN3 processing completed for {len(humann4_results)} samples")
    
    # Return the results
    return {
        "kneaddata_results": kneaddata_results,
        "humann4_results": humann4_results
    }

def main():
    """Main function to run the preprocessing pipeline."""
    # Parse arguments
    args = parse_args()
    
    # Setup logging
    log_level = getattr(logging, args.log_level.upper())
    setup_logger(log_file=args.log_file, log_level=log_level)
    
    start_time = time.time()
    log_print("Starting HUMAnN3 Tools Preprocessing Pipeline", level="info")
    
    # Check installations (skip KneadData check if we're skipping KneadData processing)
    if not args.skip_kneaddata:
        kneaddata_ok, kneaddata_version = check_kneaddata_installation()
        if not kneaddata_ok:
            log_print(f"ERROR: KneadData not properly installed: {kneaddata_version}", level="error")
            sys.exit(1)
        log_print(f"Using KneadData version: {kneaddata_version}", level="info")
    
    humann4_ok, humann4_version = check_humann4_installation()
    if not humann4_ok:
        log_print(f"ERROR: HUMAnN3 not properly installed: {humann4_version}", level="error")
        sys.exit(1)
    log_print(f"Using HUMAnN3 version: {humann4_version}", level="info")
    
    # Prepare KneadData options
    kneaddata_options = {}
    if args.paired:
        kneaddata_options["decontaminate-pairs"] = args.decontaminate_pairs
    
    # Prepare HUMAnN3 options
    humann4_options = {}
    if args.bypass_prescreen:
        humann4_options["bypass-prescreen"] = True
    if args.bypass_nucleotide_index:
        humann4_options["bypass-nucleotide-index"] = True
    if args.bypass_translated_search:
        humann4_options["bypass-translated-search"] = True
    
    # Run the preprocessing pipeline
    results = run_preprocessing_pipeline(
        input_files=args.input_fastq,
        output_dir=args.output_dir,
        threads=args.threads,
        kneaddata_dbs=args.kneaddata_dbs,
        nucleotide_db=args.humann4_nucleotide_db,
        protein_db=args.humann4_protein_db,
        paired=args.paired,
        kneaddata_options=kneaddata_options,
        humann4_options=humann4_options,
        skip_kneaddata=args.skip_kneaddata,
        kneaddata_output_files=args.kneaddata_output_files,
        kneaddata_output_pattern=args.kneaddata_output_pattern,
        kneaddata_output_dir=args.kneaddata_output_dir
    )
    
    # Generate a summary of results
    log_print("\nPreprocessing Pipeline Summary:", level="info")
    
    if args.skip_kneaddata:
        log_print(f"KneadData: Skipped (using existing files)", level="info")
    else:
        log_print(f"KneadData: Processed {len(results['kneaddata_results'])} samples", level="info")
    
    log_print(f"HUMAnN3: Processed {len(results['humann4_results'])} samples", level="info")
    
    # Print output locations
    log_print("\nOutput Directories:", level="info")
    log_print(f"  KneadData Output: {args.kneaddata_output_dir or os.path.join(args.output_dir, 'kneaddata_output')}", level="info")
    log_print(f"  HUMAnN3 Output: {args.humann4_output_dir or os.path.join(args.output_dir, 'humann4_output')}", level="info")
    
    # Print output files summary
    log_print("\nOutput Files Summary:", level="info")
    
    if results['humann4_results']:
        path_count = sum(1 for sample in results['humann4_results'].values() if sample.get('pathabundance'))
        gene_count = sum(1 for sample in results['humann4_results'].values() if sample.get('genefamilies'))
        pathcov_count = sum(1 for sample in results['humann4_results'].values() if sample.get('pathcoverage'))
        metaphlan_count = sum(1 for sample in results['humann4_results'].values() if sample.get('metaphlan'))
        
        if path_count > 0:
            pathabundance_dir = os.path.join(args.humann4_output_dir or os.path.join(args.output_dir, 'humann4_output'), "PathwayAbundance")
            log_print(f"  Pathway Abundance Files: {path_count} (in {pathabundance_dir})", level="info")
        
        if gene_count > 0:
            genefamilies_dir = os.path.join(args.humann4_output_dir or os.path.join(args.output_dir, 'humann4_output'), "GeneFamilies")
            log_print(f"  Gene Families Files: {gene_count} (in {genefamilies_dir})", level="info")
            
        if pathcov_count > 0:
            pathcoverage_dir = os.path.join(args.humann4_output_dir or os.path.join(args.output_dir, 'humann4_output'), "PathwayCoverage")
            log_print(f"  Pathway Coverage Files: {pathcov_count} (in {pathcoverage_dir})", level="info")
            
        if metaphlan_count > 0:
            metaphlan_dir = os.path.join(args.humann4_output_dir or os.path.join(args.output_dir, 'humann4_output'), "MetaphlanFiles")
            log_print(f"  MetaPhlAn Files: {metaphlan_count} (in {metaphlan_dir})", level="info")
            
        # Print detailed sample summary table
        print_sample_summary_table(results)
    
    # Print elapsed time
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    hours, minutes = divmod(minutes, 60)
    
    log_print(f"\nTotal preprocessing time: {int(hours)}h {int(minutes)}m {int(seconds)}s", level="info")
    
    # Print next steps
    log_print("\nNext Steps:", level="info")
    log_print("  To process these files with humann4_tools, you can use:", level="info")
    log_print("  humann4_tools --pathway-dir <pathabundance_dir> --gene-dir <genefamilies_dir> --sample-key <metadata.csv> --output-dir <output_dir>", level="info")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="HUMAnN3 Tools Preprocessing: Run KneadData and HUMAnN3 preprocessing steps"
    )
    
    # Required arguments
    parser.add_argument("--input-fastq", nargs="+", required=True,
                      help="Input FASTQ file(s) for preprocessing")
    
    # Output options
    parser.add_argument("--output-dir", default="./preprocessing_output",
                      help="Directory for output files")
    
    # KneadData options
    kneaddata_group = parser.add_argument_group("KneadData Options")
    kneaddata_group.add_argument("--skip-kneaddata", action="store_true",
                               help="Skip KneadData processing and use existing KneadData output files")
    kneaddata_group.add_argument("--kneaddata-dbs", nargs="+",
                               help="Path(s) to KneadData reference database(s)")
    kneaddata_group.add_argument("--kneaddata-output-dir",
                               help="Directory for KneadData output files")
    kneaddata_group.add_argument("--kneaddata-output-files", nargs="+",
                               help="Existing KneadData output files to use when --skip-kneaddata is specified")
    kneaddata_group.add_argument("--kneaddata-output-pattern",
                               help="Pattern to find KneadData output files (e.g. '/path/to/kneaddata/{sample}*paired*.fastq')")
    kneaddata_group.add_argument("--decontaminate-pairs", default="strict", 
                               choices=["strict", "lenient", "unpaired"],
                               help="Method for decontaminating paired-end reads (default: strict)")
    
    # HUMAnN3 options
    humann4_group = parser.add_argument_group("HUMAnN3 Options")
    humann4_group.add_argument("--humann4-nucleotide-db",
                             help="Path to HUMAnN3 nucleotide database (ChocoPhlAn)")
    humann4_group.add_argument("--humann4-protein-db",
                             help="Path to HUMAnN3 protein database (UniRef)")
    humann4_group.add_argument("--humann4-output-dir",
                             help="Directory for HUMAnN3 output files")
    humann4_group.add_argument("--bypass-prescreen", action="store_true",
                             help="Bypass the MetaPhlAn taxonomic prescreen")
    humann4_group.add_argument("--bypass-nucleotide-index", action="store_true",
                             help="Bypass the nucleotide index database")
    humann4_group.add_argument("--bypass-translated-search", action="store_true",
                             help="Bypass the translated search")
    
    # Shared options
    parser.add_argument("--paired", action="store_true",
                      help="Input files are paired-end reads")
    parser.add_argument("--threads", type=int, default=1,
                      help="Number of threads to use")
    parser.add_argument("--log-file",
                      help="Path to log file")
    parser.add_argument("--log-level", default="INFO",
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                      help="Logging level")
    
    return parser.parse_args()