# humann4_tools/preprocessing/humann4_run.py
import os
import subprocess
import logging
import shutil
from src.humann4_tools.utils.cmd_utils import run_cmd
from src.humann4_tools.logger import log_print
from src.humann4_tools.utils.resource_utils import track_peak_memory

def check_humann4_installation():
    """Check if HUMAnN3 is installed and available."""
    try:
        result = subprocess.run(["humann", "--version"], 
                               capture_output=True, text=True, check=False)
        if result.returncode == 0:
            return True, result.stdout.strip()
        return False, "HUMAnN3 command exists but returned an error"
    except FileNotFoundError:
        return False, "HUMAnN3 not found in PATH"
    

def process_single_sample_humann4(input_file, sample_id=None, output_dir=None, 
                                 threads=1, nucleotide_db=None, protein_db=None, 
                                 additional_options=None, logger=None, pathabdirectory=None,
                                 genedirectory=None, pathcovdirectory=None, metadirectory=None):
    """Process a single sample with HUMAnN3."""
    if logger is None:
        logger = logging.getLogger('humann4_analysis')
    
    if sample_id is None:
        sample_id = os.path.basename(input_file).split('.')[0]

    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), "humann4_output", sample_id)

    if pathabdirectory is None:
        pathabdirectory = os.path.join(output_dir, "Pathabundance")

    if genedirectory is None:
        genedirectory = os.path.join(output_dir, "GeneFamilies")

    if pathcovdirectory is None:
        pathcovdirectory = os.path.join(output_dir, "PathCoverage")

    if metadirectory is None:
        metadirectory = os.path.join(output_dir, "MetaphlanFiles")

    
    os.makedirs(output_dir, exist_ok=True)
    
    # Build command
    cmd = ["humann", "--input", input_file, "--output", output_dir]
    
    # Add threads (per sample)
    cmd.extend(["--threads", str(threads)])
    
    # Add database paths if provided
    if nucleotide_db:
        cmd.extend(["--nucleotide-database", str(nucleotide_db)])
    if protein_db:
        cmd.extend(["--protein-database", str(protein_db)])
    
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
        logger.error(f"HUMAnN3 run failed for sample {sample_id}")
        return None
    
# Now collect output file paths in one pass over the entire output_dir
    output_files = {
        'genefamilies': None,
        'pathabundance': None,
        'pathcoverage': None,
        'metaphlan': None
    }
    
    # Find the output files
    for root, dirs, files in os.walk(output_dir):
        for f in files:
            f_lower = f.lower()
            # Must be a .tsv
            if not f_lower.endswith(".tsv"):
                continue
            full_path = os.path.join(root, f)
            
            if "genefamilies" in f_lower:
                output_files["genefamilies"] = full_path
            elif "pathabundance" in f_lower:
                output_files["pathabundance"] = full_path
            elif "pathcoverage" in f_lower:
                output_files["pathcoverage"] = full_path
            elif "metaphlan_bugs_list" in f_lower:
                output_files["metaphlan"] = full_path
    
    # Log which files were found and which are missing
    found_files = [k for k, v in output_files.items() if v is not None]
    missing_files = [k for k, v in output_files.items() if v is None]
    
    logger.info(f"Found output files: {', '.join(found_files)}")
    if missing_files:
        logger.warning(f"Missing output files: {', '.join(missing_files)}")
    
    logger.info(f"HUMAnN3 completed for sample {sample_id}")

    # Send to respective directories
    pathabundance_dir = pathabdirectory
    genefamilies_dir = genedirectory
    pathcoverage_dir = pathcovdirectory
    metaphlan_dir = metadirectory

    # Create target directories if they do not exist
    os.makedirs(pathabundance_dir, exist_ok=True)
    os.makedirs(genefamilies_dir, exist_ok=True)
    
    if pathcoverage_dir:
        os.makedirs(pathcoverage_dir, exist_ok=True)
    if metaphlan_dir:
        os.makedirs(metaphlan_dir, exist_ok=True)

    # Process the output files
    try:
        path_file = output_files.get("pathabundance")
        if path_file and os.path.isfile(path_file):
            new_location = os.path.join(pathabundance_dir, os.path.basename(path_file))
            shutil.copy(path_file, new_location)
            logger.info(f"Copied pathabundance for {sample_id} to {new_location}")

        gene_file = output_files.get("genefamilies")
        if gene_file and os.path.isfile(gene_file):
            new_location = os.path.join(genefamilies_dir, os.path.basename(gene_file))
            shutil.copy(gene_file, new_location)
            logger.info(f"Copied genefamilies for {sample_id} to {new_location}")

        path_coverage_file = output_files.get("pathcoverage")
        if pathcoverage_dir and path_coverage_file and os.path.isfile(path_coverage_file):
            new_location = os.path.join(pathcoverage_dir, os.path.basename(path_coverage_file))
            shutil.copy(path_coverage_file, new_location)
            logger.info(f"Copied pathcoverage for {sample_id} to {new_location}")

        metaphlan_file = output_files.get("metaphlan")
        if metaphlan_dir and metaphlan_file and os.path.isfile(metaphlan_file):
            new_location = os.path.join(metaphlan_dir, os.path.basename(metaphlan_file))
            shutil.copy(metaphlan_file, new_location)
            logger.info(f"Copied metaphlan for {sample_id} to {new_location}")
    except Exception as e:
        logger.error(f"Error processing output files for sample {sample_id}: {str(e)}")
        logger.debug(f"output_files content: {output_files}")

    return output_files


# Add support for parallel processing
@track_peak_memory
def run_humann4_parallel(input_files, output_dir, threads=1, max_parallel=None,
                        nucleotide_db=None, protein_db=None, additional_options=None, 
                        logger=None):
    """
    Run HUMAnN3 on multiple samples in parallel.
    
    Args:
        input_files: List of input FASTQ files from KneadData
        output_dir: Base directory for outputs
        threads: Number of threads per sample
        max_parallel: Maximum number of parallel samples (None = CPU count)
        nucleotide_db: Path to nucleotide database
        protein_db: Path to protein database
        additional_options: Dict of additional HUMAnN3 options
        logger: Logger instance
        
    Returns:
        Dict mapping sample IDs to output files
    """
    from src.humann4_tools.preprocessing.parallel import run_parallel
    
    if logger is None:
        logger = logging.getLogger('humann4_analysis')
    
    # Create sample list
    #TODO This could fail if the sample name is not the first part of the file name
    sample_list = []
    for file in input_files:
        sample_name = os.path.basename(file).split('_')[0]
        sample_list.append((sample_name, file))
    
    # Prepare common arguments for all samples
    kwargs = {
        'output_dir': output_dir,
        'threads': threads,
        'nucleotide_db': nucleotide_db,
        'protein_db': protein_db,
        'additional_options': additional_options,
        'logger': logger
    }
    
    # Run in parallel
    results = run_parallel(sample_list, process_single_sample_humann4, 
                          max_workers=max_parallel, **kwargs)
    
    # Post-process to ensure we found metaphlan files
    for sample_id, sample_outputs in results.items():
        # If we didn't find the metaphlan file in the main function, try to find it now
        if isinstance(sample_outputs, dict) and sample_outputs.get('metaphlan') is None:
            sample_dir = os.path.join(output_dir, sample_id)
            if os.path.exists(sample_dir):
                for root, dirs, files in os.walk(sample_dir):
                    for file in files:
                        if "metaphlan_bugs_list" in file.lower() and file.endswith(".tsv"):
                            sample_outputs['metaphlan'] = os.path.join(root, file)
                            logger.info(f"Post-process: Found MetaPhlAn output for {sample_id}: {file}")
                            break
    return results

def run_humann4(input_files, output_dir, threads=1, nucleotide_db=None, 
               protein_db=None, additional_options=None, logger=None,  pathabdirectory=None,
                genedirectory=None, pathcovdirectory=None, metadirectory=None):
    """
    Run HUMAnN3 on input sequence files.
    
    Args:
        input_files: List of input FASTQ files from KneadData
        output_dir: Directory for HUMAnN3 output
        threads: Number of threads to use
        nucleotide_db: Path to nucleotide database
        protein_db: Path to protein database
        additional_options: Dict of additional HUMAnN3 options
        logger: Logger instance
        
    Returns:
        Dict of output file paths by type
    """
    if logger is None:
        logger = logging.getLogger('humann4_analysis')
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each input file
    output_files = {}
    for input_file in input_files:
        # Extract sample name from file path
        sample_name = os.path.basename(input_file).split("_")[0]
        sample_output_dir = os.path.join(output_dir, sample_name)
        os.makedirs(sample_output_dir, exist_ok=True)
        
        # Basic command
        cmd = ["humann", "--input", input_file, "--output", sample_output_dir]
        
        # Add threads - ensure it's a string
        cmd.extend(["--threads", str(threads)])
        
        # Add database paths if provided
        if nucleotide_db:
            cmd.extend(["--nucleotide-database", str(nucleotide_db)])
        if protein_db:
            cmd.extend(["--protein-database", str(protein_db)])
        
        # Add any additional options
        if additional_options:
            for key, value in additional_options.items():
                if value is True:
                    cmd.append(f"--{key}")
                elif value is not None:
                    # Convert value to string to ensure it's not a complex object
                    cmd.extend([f"--{key}", str(value)])
        
        # Ensure all command elements are strings before joining
        str_cmd = [str(item) for item in cmd]
        
        # Run HUMAnN3
        logger.info(f"Running HUMAnN3 for sample {sample_name}: {' '.join(str_cmd)}")
        success = run_cmd(cmd, exit_on_error=False)
        
        if not success:
            logger.error(f"HUMAnN3 run failed for sample {sample_name}")
            continue
        
    # Now collect output file paths in one pass over the entire output_dir
        sample_outputs = {
            'genefamilies': None,
            'pathabundance': None,
            'pathcoverage': None,
            'metaphlan': None
        }
        
        # Use os.walk to find the output files
        for root, dirs, files in os.walk(output_dir):
            for f in files:
                f_lower = f.lower()
                if not f_lower.endswith(".tsv"):
                    continue
                full_path = os.path.join(root, f)
                
                if "genefamilies" in f_lower:
                    sample_outputs["genefamilies"] = full_path
                elif "pathabundance" in f_lower:
                    sample_outputs["pathabundance"] = full_path
                elif "pathcoverage" in f_lower:
                    sample_outputs["pathcoverage"] = full_path
                elif "metaphlan_bugs_list" in f_lower:
                    sample_outputs["metaphlan"] = full_path
        
        # Log which files were found and which are missing
        found_files = [k for k, v in sample_outputs.items() if v is not None]
        missing_files = [k for k, v in sample_outputs.items() if v is None]
        
        logger.info(f"Sample {sample_name} output files found: {', '.join(found_files)}")
        if missing_files:
            logger.warning(f"Sample {sample_name} missing files: {', '.join(missing_files)}")
            
        output_files[sample_name] = sample_outputs

    # Send to respective directories
        pathabundance_dir = pathabdirectory
        genefamilies_dir = genedirectory
        pathcoverage_dir = pathcovdirectory
        metaphlan_dir = metadirectory
    
        # Create target directories if they do not exist
        os.makedirs(pathabundance_dir, exist_ok=True)
        os.makedirs(genefamilies_dir, exist_ok=True)

        for sample_id, file_paths in sample_outputs.items():


            if pathabdirectory is None:
                pathabdirectory = os.path.join(output_dir, "PathwayAbundance")

            if genedirectory is None:
                genedirectory = os.path.join(output_dir, "GeneFamilies")

            if pathcovdirectory is None:
                pathcovdirectory = os.path.join(output_dir, "PathwayCoverage")

            if metadirectory is None:
                metadirectory = os.path.join(output_dir, "MetaphlanFiles")

            path_file = file_paths.get("pathabundance")
            if path_file and os.path.isfile(path_file):
                new_location = os.path.join(pathabundance_dir, os.path.basename(path_file))
                shutil.copy(path_file, new_location)
                print(f"Copied pathabundance for {sample_id} to {new_location}")

    
            gene_file = file_paths.get("genefamilies")
            if gene_file and os.path.isfile(gene_file):
                new_location = os.path.join(genefamilies_dir, os.path.basename(gene_file))
                shutil.copy(gene_file, new_location)
                print(f"Copied genefamilies for {sample_id} to {new_location}")

            path_coverage_file = file_paths.get("pathcoverage")
            if path_coverage_file and os.path.isfile(path_coverage_file):
                new_location = os.path.join(pathcoverage_dir, os.path.basename(path_coverage_file))
                shutil.copy(path_coverage_file, new_location)
                print(f"Copied pathcoverage for {sample_id} to {new_location}")

            metaphlan_file = file_paths.get("metaphlan")
            if metaphlan_file and os.path.isfile(metaphlan_file):
                new_location = os.path.join(metaphlan_dir, os.path.basename(metaphlan_file))
                shutil.copy(metaphlan_file, new_location)
                print(f"Copied metaphlan for {sample_id} to {new_location}")
        
    logger.info(f"HUMAnN3 completed for {len(output_files)} samples")
    return output_files