# humann4_tools/humann4_tools/utils/file_utils.py
import os
import logging

def check_file_exists(filepath, description):
    """Check if a file exists and is readable."""
    logger = logging.getLogger('humann4_analysis')
    
    if not os.path.isfile(filepath):
        logger.error(f"ERROR: {description} file does not exist: {filepath}")
        return False
    
    if not os.access(filepath, os.R_OK):
        logger.error(f"ERROR: {description} file exists but is not readable: {filepath}")
        return False
    
    return True

def check_file_exists_with_logger(filepath, description, logger):
    """
    Check if file exists & is readable (logger-based).
    """
    if not os.path.isfile(filepath):
        logger.error(f"{description} file does not exist: {filepath}")
        return False
    if not os.access(filepath, os.R_OK):
        logger.error(f"{description} file is not readable: {filepath}")
        return False
    return True

def sanitize_filename(filename):
    """Replace invalid filename characters with underscores."""
    import re
    return re.sub(r'[<>:"/\\|?*]', '_', filename)

def strip_suffix(col):
    """
    Enhanced function to remove HUMAnN3 abundance suffix from column names.
    Handles various suffix patterns more robustly.
    """
    # First, try exact matches with known suffixes
    suffixes = [
        ".paired_Abundance-CPM", 
        "_Abundance-CPM",
        ".paired_Abundance-RELAB", 
        "_Abundance-RELAB",
        "-cpm",
        "-relab",
        ".cpm",
        ".relab"
    ]
    
    for suffix in suffixes:
        if col.lower().endswith(suffix.lower()):
            return col[:-len(suffix)]
    
    # If no exact match, try pattern-based detection
    
    # Pattern 1: Sample.UNIT or Sample.something.UNIT
    parts = col.split('.')
    if len(parts) > 1:
        last_part = parts[-1].lower()
        if any(unit in last_part for unit in ['cpm', 'relab', 'abundance']):
            return '.'.join(parts[:-1])
    
    # Pattern 2: Sample_UNIT or Sample-UNIT or variations
    for marker in ['_abundance', '-abundance', '_cpm', '-cpm', '_relab', '-relab']:
        pos = col.lower().find(marker)
        if pos > 0:
            return col[:pos]
    
    # No match found
    return col

def strip_suffixes_from_file_headers(file_path, logger=None):
    """
    Remove HUMAnN3 abundance suffixes from column headers in a file.
    Includes safety checks for file format compatibility.
    
    Args:
        file_path: Path to the file with headers to strip
        logger: Optional logger for messages
        
    Returns:
        True if successful, False otherwise
    """
    if logger is None:
        import logging
        logger = logging.getLogger('humann4_analysis')
    
    try:
        # Read the file
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            logger.warning(f"Empty file: {file_path}")
            return False
        
        # Find the first non-comment line (header line)
        header_index = 0
        for i, line in enumerate(lines):
            if not line.startswith('#'):
                header_index = i
                break
        
        # Get the header line and split it
        header = lines[header_index].strip()
        cols = header.split('\t')
        
        # Check if this is a typical HUMAnN3 output file with named sample columns
        # First check if we have at least a feature column and one sample column
        if len(cols) < 2:
            logger.warning(f"File has fewer than 2 columns, cannot process: {file_path}")
            return False
        
        # Check if columns appear to be numeric IDs rather than sample names with suffixes
        numeric_cols = sum(1 for col in cols[1:] if col.isdigit())
        if numeric_cols > 0 and numeric_cols / (len(cols) - 1) > 0.5:  # If >50% are numeric
            logger.info(f"File appears to have numeric column headers, skipping suffix stripping: {file_path}")
            return True
        
        # Check if any of the columns have a recognizable suffix pattern
        has_suffix_pattern = False
        for col in cols[1:]:  # Skip first column (feature ID)
            # Check for common patterns we expect in HUMAnN3 output
            if any(suffix in col.lower() for suffix in ['abundance', 'cpm', 'relab']):
                has_suffix_pattern = True
                break
            
            # Check for patterns with separators
            for sep in ['.', '_', '-']:
                if sep in col and any(unit in col.lower().split(sep)[-1] for unit in ['cpm', 'relab']):
                    has_suffix_pattern = True
                    break
        
        if not has_suffix_pattern:
            logger.info(f"No recognizable suffix patterns found in column headers: {file_path}")
            return True
        
        # Now we've confirmed this is a suitable file for suffix stripping
        # Apply strip_suffix to each column except the first one (which is usually the feature ID)
        new_cols = [cols[0]]
        change_count = 0
        
        for col in cols[1:]:
            new_col = strip_suffix(col)
            new_cols.append(new_col)
            if new_col != col:
                change_count += 1
                logger.debug(f"Stripped suffix: '{col}' -> '{new_col}'")
        
        # Only update the file if we actually made changes
        if change_count == 0:
            logger.info(f"No columns were modified in: {file_path}")
            return True
        
        # Update the header line
        lines[header_index] = '\t'.join(new_cols) + '\n'
        
        # Write the file back
        with open(file_path, 'w') as f:
            f.writelines(lines)
        
        logger.info(f"Successfully stripped {change_count} suffixes from headers in: {file_path}")
        return True
    except Exception as e:
        if logger:
            logger.error(f"Error stripping suffixes from headers in {file_path}: {str(e)}")
        return False
