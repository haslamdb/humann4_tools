"""
Input handler module for HUMAnN3 Tools.

This module provides standardized methods to handle different input formats:
1. Direct file list - A list of input files passed directly
2. Sample list file - A file containing sample names or paths
3. Metadata file - A CSV file with sample information

It returns a standardized sample dictionary that can be used by any step 
in the workflow.
"""

import os
import glob
import logging
import pandas as pd
from typing import List, Dict, Union, Optional, Tuple

logger = logging.getLogger('humann4_tools')

def find_sample_files(sample_id: str, 
                     search_dir: str, 
                     file_pattern: Optional[str] = None,
                     r1_suffix: Optional[str] = None, 
                     r2_suffix: Optional[str] = None,
                     paired: bool = False) -> List[str]:
    """
    Find sequence files for a sample based on patterns.
    
    Args:
        sample_id: Sample identifier
        search_dir: Directory to search for sequence files
        file_pattern: Optional pattern to use for file search
        r1_suffix: Suffix for R1 files in paired-end mode 
        r2_suffix: Suffix for R2 files in paired-end mode
        paired: Whether to look for paired-end files
        
    Returns:
        List of file paths matching the sample
    """
    files = []
    
    # If a specific file pattern is provided
    if file_pattern:
        pattern = file_pattern.replace('{sample}', sample_id)
        search_pattern = os.path.join(search_dir, pattern)
        files = sorted(glob.glob(search_pattern))
        logger.debug(f"Searching with pattern {search_pattern}, found {len(files)} files")
        if paired and len(files) >= 2:
            # Take only the first two files for paired mode
            files = files[:2]
        
    # For paired-end reads with specific suffixes
    elif paired and r1_suffix and r2_suffix:
        r1_pattern = f"{sample_id}{r1_suffix}"
        r2_pattern = f"{sample_id}{r2_suffix}"
        r1_search = os.path.join(search_dir, r1_pattern)
        r2_search = os.path.join(search_dir, r2_pattern)
        
        r1_files = sorted(glob.glob(r1_search))
        r2_files = sorted(glob.glob(r2_search))
        
        logger.debug(f"Searching for R1: {r1_search}, found {len(r1_files)} files")
        logger.debug(f"Searching for R2: {r2_search}, found {len(r2_files)} files")
        
        # Ensure we have matching pairs
        if len(r1_files) > 0 and len(r2_files) > 0:
            # Only add one pair (the first pair found)
            files = [r1_files[0], r2_files[0]]
            logger.debug(f"Added paired files: {os.path.basename(r1_files[0])} and {os.path.basename(r2_files[0])}")
        else:
            logger.warning(f"Couldn't find matching paired files for sample {sample_id}")
    
    # Generic pattern search
    else:
        # For paired-end data, look for R1/R2 pairs explicitly
        if paired:
            # Try to find R1 and R2 files
            r1_patterns = [
                f"{sample_id}_R1.fastq",
                f"{sample_id}_R1.fastq.gz",
                f"{sample_id}_1.fastq",
                f"{sample_id}_1.fastq.gz",
                f"{sample_id}.R1.fastq",
                f"{sample_id}.R1.fastq.gz"
            ]
            
            r2_patterns = [
                f"{sample_id}_R2.fastq",
                f"{sample_id}_R2.fastq.gz",
                f"{sample_id}_2.fastq", 
                f"{sample_id}_2.fastq.gz",
                f"{sample_id}.R2.fastq",
                f"{sample_id}.R2.fastq.gz"
            ]
            
            r1_file = None
            r2_file = None
            
            # Find R1 file
            for pattern in r1_patterns:
                search_pattern = os.path.join(search_dir, pattern)
                matches = glob.glob(search_pattern)
                logger.debug(f"Trying R1 pattern {search_pattern}, found {len(matches)} files")
                if matches:
                    r1_file = matches[0]
                    break
                    
            # Find R2 file
            for pattern in r2_patterns:
                search_pattern = os.path.join(search_dir, pattern)
                matches = glob.glob(search_pattern)
                logger.debug(f"Trying R2 pattern {search_pattern}, found {len(matches)} files")
                if matches:
                    r2_file = matches[0]
                    break
            
            # Add both files if found
            if r1_file and r2_file:
                files = [r1_file, r2_file]
                logger.debug(f"Found paired files: {os.path.basename(r1_file)} and {os.path.basename(r2_file)}")
            else:
                logger.warning(f"Could not find complete paired files for sample {sample_id}")
        else:
            # Single-end mode - try common patterns
            common_patterns = [
                f"{sample_id}.fastq", 
                f"{sample_id}.fq", 
                f"{sample_id}.fastq.gz", 
                f"{sample_id}.fq.gz",
                f"{sample_id}_*.fastq.gz",
                f"{sample_id}*.fastq.gz",
                f"{sample_id}_R1.fastq",
                f"{sample_id}_R1.fastq.gz",
                f"{sample_id}_*.fastq"
            ]
            
            for pattern in common_patterns:
                search_pattern = os.path.join(search_dir, pattern)
                matches = glob.glob(search_pattern)
                logger.debug(f"Trying pattern {search_pattern}, found {len(matches)} files")
                if matches:
                    # Only take the first file found for single-end mode
                    files = [matches[0]]
                    logger.debug(f"Selected file: {os.path.basename(matches[0])}")
                    break
    
    if files:
        logger.debug(f"Found {len(files)} files for sample {sample_id}: {[os.path.basename(f) for f in files]}")
    
    return files

def parse_metadata_file(metadata_path: str, 
                       sample_col: Optional[str] = None,
                       group_col: Optional[str] = None) -> Tuple[List[str], Dict, str, str]:
    """
    Parse a metadata file to extract sample IDs and optional grouping information.
    
    Args:
        metadata_path: Path to the metadata CSV file
        sample_col: Column name with sample identifiers (auto-detected if None)
        group_col: Column name for sample grouping (auto-detected if None)
        
    Returns:
        Tuple of (sample_list, metadata_dict, sample_col, group_col)
            - sample_list: List of sample IDs
            - metadata_dict: Dictionary mapping sample IDs to their metadata rows
            - sample_col: The column name used for sample IDs
            - group_col: The column name used for grouping
    """
    # Read metadata file
    try:
        metadata_df = pd.read_csv(metadata_path)
        logger.info(f"Read metadata file with {len(metadata_df)} rows and columns: {metadata_df.columns.tolist()}")
    except Exception as e:
        logger.error(f"Error reading metadata file: {e}")
        return [], {}, None, None
    
    # Identify sample ID column if not specified
    if not sample_col:
        common_sample_cols = ['SampleID', 'Sample_ID', 'SampleName', 'sample_id', 
                             'sample_name', 'Sample', 'ID']
        for col in common_sample_cols:
            if col in metadata_df.columns:
                sample_col = col
                logger.info(f"Auto-detected sample ID column: {sample_col}")
                break
        
        if not sample_col:
            sample_col = metadata_df.columns[0]
            logger.warning(f"Could not auto-detect sample ID column, using the first column: {sample_col}")
    
    # Verify sample column exists
    if sample_col not in metadata_df.columns:
        logger.error(f"Sample column '{sample_col}' not found in metadata file")
        return [], {}, None, None
    
    # Identify grouping column if not specified
    if not group_col:
        common_group_cols = ['Group', 'Treatment', 'Condition', 'group', 'treatment', 'condition']
        for col in common_group_cols:
            if col in metadata_df.columns:
                group_col = col
                logger.info(f"Auto-detected group column: {group_col}")
                break
    
    # Extract sample IDs and create metadata dictionary
    sample_list = []
    metadata_dict = {}
    
    for _, row in metadata_df.iterrows():
        sample_id = str(row[sample_col]).strip()
        
        # Skip empty sample IDs
        if not sample_id or pd.isna(sample_id):
            continue
            
        sample_list.append(sample_id)
        metadata_dict[sample_id] = row.to_dict()
    
    logger.info(f"Extracted {len(sample_list)} sample IDs from metadata file")
    return sample_list, metadata_dict, sample_col, group_col

def read_samples_file(samples_file_path: str) -> Dict[str, List[str]]:
    """
    Read a tab-delimited samples file.
    
    Format:
        sample_id    file1 file2 ...
    
    Args:
        samples_file_path: Path to the tab-delimited samples file
        
    Returns:
        Dictionary mapping sample IDs to lists of file paths
    """
    samples = {}
    
    try:
        with open(samples_file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                    
                parts = line.split('\t')
                if len(parts) < 2:
                    parts = line.split()  # Try space delimiter if tab doesn't work
                
                if len(parts) >= 2:
                    sample_id = parts[0]
                    file_paths = parts[1:]
                    
                    # Allow space-separated file paths in the second column
                    if len(file_paths) == 1 and ' ' in file_paths[0]:
                        file_paths = file_paths[0].split()
                    
                    # Check if files exist
                    valid_files = []
                    for file_path in file_paths:
                        if os.path.exists(file_path):
                            valid_files.append(file_path)
                        else:
                            logger.warning(f"File not found: {file_path}")
                    
                    if valid_files:
                        samples[sample_id] = valid_files
                    
        logger.info(f"Loaded {len(samples)} samples from {samples_file_path}")
        return samples
    except Exception as e:
        logger.error(f"Error reading samples file {samples_file_path}: {e}")
        return {}

def collect_files_from_metadata(metadata_file: str, 
                               seq_dir: str, 
                               sample_col: Optional[str] = None,
                               group_col: Optional[str] = None,
                               r1_col: Optional[str] = None, 
                               r2_col: Optional[str] = None,
                               file_pattern: Optional[str] = None,
                               r1_suffix: Optional[str] = None, 
                               r2_suffix: Optional[str] = None,
                               paired: bool = False) -> Dict[str, Dict]:
    """
    Collect sample information from metadata file and find associated files.
    
    Args:
        metadata_file: Path to metadata CSV file
        seq_dir: Directory to search for sequence files
        sample_col: Column name for sample identifiers
        group_col: Column name for grouping
        r1_col: Column name for R1 file paths
        r2_col: Column name for R2 file paths
        file_pattern: Pattern for finding files (e.g., "{sample}_S*_R*.fastq.gz")
        r1_suffix: Suffix for R1 files (e.g., "_R1.fastq.gz")
        r2_suffix: Suffix for R2 files (e.g., "_R2.fastq.gz")
        paired: Whether samples are paired-end
        
    Returns:
        Dictionary mapping sample IDs to dicts with 'files' and 'metadata'
    """
    # Parse metadata to get sample IDs
    sample_list, metadata_dict, detected_sample_col, detected_group_col = parse_metadata_file(
        metadata_file, sample_col, group_col
    )
    
    if not sample_list:
        logger.error("No samples found in metadata file")
        return {}
    
    # Use detected columns if none were specified
    if not sample_col and detected_sample_col:
        sample_col = detected_sample_col
    if not group_col and detected_group_col:
        group_col = detected_group_col
    
    # Check if there are explicit file path columns
    has_file_cols = False
    if r1_col:
        try:
            # Read metadata without parsing to check column existence
            temp_df = pd.read_csv(metadata_file)
            has_file_cols = r1_col in temp_df.columns
            if paired:
                has_file_cols = has_file_cols and r2_col in temp_df.columns
        except:
            has_file_cols = False
    
    # Collect sample information
    samples = {}
    
    for sample_id in sample_list:
        # Initialize sample entry with metadata
        samples[sample_id] = {
            'metadata': metadata_dict.get(sample_id, {}),
            'files': []
        }
        
        # Case 1: File paths are in the metadata
        if has_file_cols:
            r1_path = metadata_dict[sample_id].get(r1_col)
            
            if r1_path:
                # Check if path is relative or absolute
                if not os.path.isabs(r1_path):
                    r1_path = os.path.join(seq_dir, r1_path)
                
                if paired and r2_col:
                    r2_path = metadata_dict[sample_id].get(r2_col)
                    if r2_path:
                        if not os.path.isabs(r2_path):
                            r2_path = os.path.join(seq_dir, r2_path)
                        
                        if os.path.exists(r1_path) and os.path.exists(r2_path):
                            samples[sample_id]['files'] = [r1_path, r2_path]
                else:
                    if os.path.exists(r1_path):
                        samples[sample_id]['files'] = [r1_path]
        
        # Case 2: Need to find files based on patterns
        else:
            files = find_sample_files(
                sample_id=sample_id,
                search_dir=seq_dir,
                file_pattern=file_pattern,
                r1_suffix=r1_suffix,
                r2_suffix=r2_suffix,
                paired=paired
            )
            
            if files:
                samples[sample_id]['files'] = files
            else:
                logger.warning(f"No sequence files found for sample {sample_id}")
    
    # Log the results
    valid_samples = sum(1 for s in samples.values() if s['files'])
    logger.info(f"Found files for {valid_samples} out of {len(sample_list)} samples")
    
    return samples

def parse_input_files(input_files: List[str], paired: bool = False) -> Dict[str, Dict]:
    """
    Parse a list of input files into sample dictionaries.
    
    Args:
        input_files: List of input file paths
        paired: Whether samples are paired-end
        
    Returns:
        Dictionary mapping sample IDs to dicts with 'files'
    """
    samples = {}
    
    if paired:
        # Ensure we have an even number of files
        if len(input_files) % 2 != 0:
            logger.error("Paired mode requires an even number of input files")
            return {}
        
        # Group files in pairs
        for i in range(0, len(input_files), 2):
            r1_file = input_files[i]
            r2_file = input_files[i+1]
            
            # Extract sample name from R1 file
            sample_id = os.path.basename(r1_file).split('_')[0]
            
            # Avoid duplicates by appending an index if needed
            base_id = sample_id
            idx = 1
            while sample_id in samples:
                sample_id = f"{base_id}_{idx}"
                idx += 1
            
            samples[sample_id] = {
                'files': [r1_file, r2_file],
                'metadata': {}
            }
            logger.debug(f"Added paired files for sample {sample_id}: {os.path.basename(r1_file)} and {os.path.basename(r2_file)}")
    else:
        # Single-end data: one file per sample
        for file_path in input_files:
            # Extract sample name from filename
            sample_id = os.path.basename(file_path).split('.')[0]
            
            # Avoid duplicates by appending an index if needed
            base_id = sample_id
            idx = 1
            while sample_id in samples:
                sample_id = f"{base_id}_{idx}"
                idx += 1
            
            samples[sample_id] = {
                'files': [file_path],
                'metadata': {}
            }
            logger.debug(f"Added file for sample {sample_id}: {os.path.basename(file_path)}")
    
    logger.info(f"Parsed {len(samples)} samples from {len(input_files)} input files")
    return samples

def get_input_files(args, input_type: str = "sequence") -> Dict[str, Dict]:
    """
    Process input arguments and return standardized sample dictionary.
    
    This is the main entry point for the input handler module. It detects
    the input method (direct files, sample list, metadata) and returns a
    standardized dictionary for workflow steps.
    
    Args:
        args: Parsed argument object from argparse
        input_type: Type of input files ("sequence", "kneaddata", "humann4", etc.)
        
    Returns:
        Dictionary mapping sample IDs to dicts with 'files' and 'metadata' keys
    """
    # Method 1: Direct file list
    if hasattr(args, 'input_files') and args.input_files:
        logger.info(f"Processing {len(args.input_files)} direct input files")
        return parse_input_files(args.input_files, getattr(args, 'paired', False))
    
    # Method 2: Sample list file
    elif hasattr(args, 'samples_file') and args.samples_file:
        logger.info(f"Processing samples from file: {args.samples_file}")
        sample_dict = read_samples_file(args.samples_file)
        
        # Convert to standardized format
        samples = {}
        for sample_id, files in sample_dict.items():
            samples[sample_id] = {
                'files': files,
                'metadata': {}
            }
        return samples
    
    # Method 3: Metadata-driven
    elif hasattr(args, 'metadata_file') and args.metadata_file and hasattr(args, 'seq_dir') and args.seq_dir:
        logger.info(f"Processing samples from metadata: {args.metadata_file}")
        return collect_files_from_metadata(
            metadata_file=args.metadata_file,
            seq_dir=args.seq_dir,
            sample_col=getattr(args, 'sample_col', None),
            group_col=getattr(args, 'group_col', None),
            r1_col=getattr(args, 'r1_col', None),
            r2_col=getattr(args, 'r2_col', None),
            file_pattern=getattr(args, 'file_pattern', None),
            r1_suffix=getattr(args, 'r1_suffix', None),
            r2_suffix=getattr(args, 'r2_suffix', None),
            paired=getattr(args, 'paired', False)
        )
    
    # No valid input method specified
    else:
        logger.error("No valid input method specified. Use --input-files, --samples-file, or --metadata-file + --seq-dir")
        return {}

def find_output_files(sample_id: str, output_dir: str, suffixes: List[str], file_type: str) -> List[str]:
    """
    Find output files for a sample based on common naming patterns.
    
    Args:
        sample_id: Sample identifier
        output_dir: Directory to search for output files
        suffixes: List of file suffixes to look for
        file_type: Type of files to look for (for logging)
        
    Returns:
        List of matched file paths
    """
    matched_files = []
    
    for suffix in suffixes:
        pattern = f"{sample_id}{suffix}"
        search_pattern = os.path.join(output_dir, pattern)
        matches = glob.glob(search_pattern)
        
        if matches:
            matched_files.extend(matches)
            logger.debug(f"Found {file_type} files for {sample_id} with pattern {pattern}: {len(matches)} files")
    
    if not matched_files:
        logger.warning(f"No {file_type} files found for sample {sample_id} in {output_dir}")
    
    return matched_files

def find_humann4_output_files(samples: Dict[str, Dict], output_dir: str) -> Dict[str, Dict]:
    """
    Find HUMAnN3 output files for samples in output directory.
    
    Args:
        samples: Dictionary mapping sample IDs to sample information
        output_dir: HUMAnN3 output directory
        
    Returns:
        Updated samples dictionary with 'humann4_files' key for each sample
    """
    for sample_id, sample_info in samples.items():
        # Define common suffixes for HUMAnN3 output files
        pathabundance_suffixes = [
            "_pathabundance.tsv",
            ".pathabundance.tsv",
            "_humann_pathabundance.tsv"
        ]
        
        genefamilies_suffixes = [
            "_genefamilies.tsv",
            ".genefamilies.tsv",
            "_humann_genefamilies.tsv"
        ]
        
        pathcoverage_suffixes = [
            "_pathcoverage.tsv",
            ".pathcoverage.tsv",
            "_humann_pathcoverage.tsv"
        ]
        
        # Find files for each type
        pathabundance_files = find_output_files(sample_id, output_dir, pathabundance_suffixes, "pathabundance")
        genefamilies_files = find_output_files(sample_id, output_dir, genefamilies_suffixes, "genefamilies")
        pathcoverage_files = find_output_files(sample_id, output_dir, pathcoverage_suffixes, "pathcoverage")
        
        # Update sample info with found files
        sample_info['humann4_files'] = {
            'pathabundance': pathabundance_files[0] if pathabundance_files else None,
            'genefamilies': genefamilies_files[0] if genefamilies_files else None,
            'pathcoverage': pathcoverage_files[0] if pathcoverage_files else None
        }
    
    return samples

def find_kneaddata_output_files(samples: Dict[str, Dict], output_dir: str, paired: bool = False) -> Dict[str, Dict]:
    """
    Find KneadData output files for samples in output directory.
    
    Args:
        samples: Dictionary mapping sample IDs to sample information
        output_dir: KneadData output directory
        paired: Whether to look for paired-end output files
        
    Returns:
        Updated samples dictionary with 'kneaddata_files' key for each sample
    """
    for sample_id, sample_info in samples.items():
        if paired:
            # Define common suffixes for paired KneadData output files
            paired_suffixes = [
                "_paired_1.fastq",
                "_paired_2.fastq",
                "_kneaddata_paired_1.fastq",
                "_kneaddata_paired_2.fastq"
            ]
            
            # Find potential sample directories
            sample_dirs = [output_dir]
            potential_sample_dir = os.path.join(output_dir, sample_id)
            if os.path.isdir(potential_sample_dir):
                sample_dirs.append(potential_sample_dir)
            
            # Search in all potential directories
            paired_files = []
            for search_dir in sample_dirs:
                for suffix in paired_suffixes:
                    pattern = f"{sample_id}{suffix}"
                    search_pattern = os.path.join(search_dir, pattern)
                    matches = glob.glob(search_pattern)
                    
                    if matches:
                        paired_files.extend(matches)
                        logger.debug(f"Found KneadData files for {sample_id} with pattern {pattern}: {len(matches)} files")
            
            # Sort to ensure R1 comes before R2
            paired_files.sort()
            
            if len(paired_files) >= 2:
                # Keep only the first pair
                sample_info['kneaddata_files'] = paired_files[:2]
            else:
                logger.warning(f"Not enough paired KneadData files found for sample {sample_id}")
                sample_info['kneaddata_files'] = []
        else:
            # For single-end data
            suffixes = [
                "_kneaddata.fastq",
                ".kneaddata.fastq",
                "_cleaned.fastq"
            ]
            
            single_files = find_output_files(sample_id, output_dir, suffixes, "kneaddata")
            
            if single_files:
                sample_info['kneaddata_files'] = [single_files[0]]
            else:
                logger.warning(f"No KneadData files found for sample {sample_id}")
                sample_info['kneaddata_files'] = []
    
    return samples
