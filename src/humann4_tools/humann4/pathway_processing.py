# humann4_tools/humann4_tools/humann4/pathway_processing.py
import os
import logging

from src.humann4_tools.logger import log_print
from src.humann4_tools.utils.cmd_utils import run_cmd

def process_pathway_abundance(valid_samples, pathway_dir, output_dir, output_prefix, selected_columns=None, units="cpm"):
    """
    Runs humann_renorm_table, humann_join_tables, humann_split_stratified_table, 
    then returns the path to the final unstratified pathway file (renamed to pathway_abundance.tsv).
    
    Args:
        valid_samples: List of tuples (sample, filepath) for valid samples
        pathway_dir: Directory containing original pathway files
        output_dir: Directory where processed files will be stored
        output_prefix: Prefix for output filenames
        selected_columns: Optional dict with sample key column selections
        units: Units for normalization (default: cpm, can be "relab")
        
    Returns:
        Path to unstratified pathway file, or None if processing failed
    """
    logger = logging.getLogger('humann4_analysis')
    log_print(f"PROCESSING PATHWAY ABUNDANCE FILES (using {units} units)", level='info')
    path_abundance_out = os.path.join(output_dir, "pathways", output_prefix)
    path_abundance_norm = os.path.join(path_abundance_out, "Normalized")
    os.makedirs(path_abundance_out, exist_ok=True)
    os.makedirs(path_abundance_norm, exist_ok=True)
    
    if not valid_samples:
        log_print("No valid pathway abundance files to process", level='warning')
        return None
    
    # Optionally save column selections
    if selected_columns:
        column_info_file = os.path.join(path_abundance_out, "column_selections.txt")
        try:
            with open(column_info_file, "w") as f:
                f.write(f"sample_id_column: {selected_columns['sample_id']}\n")
                f.write("grouping_columns:\n")
                for col in selected_columns.get('grouping', {}):
                    f.write(f"  - {col}\n")
        except Exception as e:
            log_print(f"Warning: Could not save column selection info: {e}", level='warning')
    
    # Build the suffix based on the units
    units_suffix = f"-{units}"
    
    # Copy + renormalize
    for (sample, src_path) in valid_samples:
        dst = os.path.join(path_abundance_out, f"{sample}_pathabundance.tsv")
        if not run_cmd(["cp", src_path, dst], exit_on_error=False):
            continue
        
        out_norm = os.path.join(path_abundance_out, f"{sample}_pathabundance{units_suffix}.tsv")
        run_cmd([
            "humann_renorm_table",
            "--input", dst,
            "--output", out_norm,
            "--units", units,
            "--update-snames"
        ], exit_on_error=False)
        
        run_cmd(["mv", out_norm, path_abundance_norm], exit_on_error=False)
    
    # Join
    norm_files = [f for f in os.listdir(path_abundance_norm) if f.endswith(f"{units_suffix}.tsv")]
    if not norm_files:
        log_print(f"WARNING: No normalized pathway files to join (with {units_suffix} suffix)", level='warning')
        return None
    
    joined_output = os.path.join(path_abundance_out, f"{output_prefix}_pathabundance{units_suffix}.tsv")
    run_cmd([
        "humann_join_tables",
        "-i", path_abundance_norm,
        "-o", joined_output
    ])
    
    if not os.path.exists(joined_output):
        log_print("WARNING: Joined pathway file not found after humann_join_tables", level='warning')
        return None
    
    # Split stratified
    run_cmd([
        "humann_split_stratified_table",
        "-i", joined_output,
        "-o", path_abundance_out
    ])
    
    # Locate unstratified file
    unstrat_file = None
    base_name = os.path.basename(joined_output).replace('.tsv', '')
    possible_patterns = [
        f"{base_name}_unstratified.tsv",
        f"{output_prefix}_pathabundance{units_suffix}_unstratified.tsv",
        f"{output_prefix}_pathabundance-{units}_unstratified.tsv"
    ]
    
    for pattern in possible_patterns:
        test_path = os.path.join(path_abundance_out, pattern)
        if os.path.isfile(test_path):
            unstrat_file = test_path
            break
    
    if not unstrat_file:
        # Try searching subdirs or any file containing 'unstratified'
        for root, dirs, files in os.walk(path_abundance_out):
            for fname in files:
                if 'unstratified' in fname.lower():
                    unstrat_file = os.path.join(root, fname)
                    break
            if unstrat_file:
                break
    
    if not unstrat_file:
        log_print("WARNING: Could not locate unstratified pathway file", level='warning')
        return None
    
    # Rename final unstratified file
    final_name = os.path.join(path_abundance_out, f"pathway_abundance{units_suffix}_unstratified.tsv")
    try:
        os.rename(unstrat_file, final_name)
        unstrat_file = final_name
        log_print(f"Renamed unstratified file to: {unstrat_file}", level='info')
    except Exception as e:
        log_print(f"WARNING: Could not rename file: {e}", level='warning')

        # Rename final unstratified file
    final_name = os.path.join(path_abundance_out, f"pathway_abundance{units_suffix}_unstratified.tsv")
    try:
        os.rename(unstrat_file, final_name)
        unstrat_file = final_name
        log_print(f"Renamed unstratified file to: {unstrat_file}", level='info')
        
        # Strip suffixes from headers in all output files
        from src.humann4_tools.utils.file_utils import strip_suffixes_from_file_headers
        strip_suffixes_from_file_headers(unstrat_file)
        
        # Also process the stratified file if it exists
        stratified_file = os.path.join(path_abundance_out, f"pathway_abundance{units_suffix}_stratified.tsv")
        if os.path.exists(stratified_file):
            strip_suffixes_from_file_headers(stratified_file)
            
    except Exception as e:
        log_print(f"WARNING: Could not rename file: {e}", level='warning')

    return unstrat_file
    
