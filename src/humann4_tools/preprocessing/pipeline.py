# humann4_tools/preprocessing/pipeline.py
import os
import logging
from src.humann4_tools.preprocessing.kneaddata import run_kneaddata, check_kneaddata_installation, run_kneaddata_parallel
from src.humann4_tools.preprocessing.humann4_run import run_humann4, check_humann4_installation, run_humann4_parallel
from src.humann4_tools.logger import log_print
from src.humann4_tools.utils.resource_utils import (
    track_peak_memory, 
    monitor_memory_usage, 
    stop_memory_monitoring
)
from src.humann4_tools.humann4.join_unstratify import process_join_unstratify, join_unstratify_humann_output

def run_preprocessing_pipeline(
    input_files,
    output_dir,
    threads=1,
    kneaddata_dbs=None,
    nucleotide_db=None,
    protein_db=None,
    kneaddata_options=None,
    humann4_options=None,
    paired=False,
    kneaddata_output_dir=None,
    humann4_output_dir=None,
    skip_kneaddata=False,
    kneaddata_output_files=None,
    kneaddata_output_pattern=None,
    logger=None
):
    """
    Run the full preprocessing pipeline: KneadData → HUMAnN3.

    This version simplifies:
      1) Avoiding duplicates,
      2) Finding correct paired files, and
      3) Logging only the essentials.

    Args:
        input_files: List of input FASTQ files
        output_dir: Base directory for outputs
        threads: Number of threads per sample
        kneaddata_dbs: Path(s) to KneadData reference database(s)
        nucleotide_db: Path to HUMAnN3 nucleotide database
        protein_db: Path to HUMAnN3 protein database
        kneaddata_options: Dict of additional KneadData options
        humann4_options: Dict of additional HUMAnN3 options
        paired: Whether input files are paired
        kneaddata_output_dir: Custom directory for KneadData outputs
        humann4_output_dir: Custom directory for HUMAnN3 outputs
        skip_kneaddata: Whether to skip KneadData processing
        kneaddata_output_files: List of existing KneadData output files to use
        kneaddata_output_pattern: Pattern to find KneadData output files
        logger: Logger instance

    Returns:
        Dict with:
          'kneaddata_files': list of all final kneaddata files
          'humann4_results': dict of HUMAnN3 outputs per sample
    """
    if logger is None:
        logger = logging.getLogger("humann4_analysis")

    # 1. Check installations (skip KneadData check if we're skipping KneadData processing)
    if not skip_kneaddata:
        kneaddata_ok, kneaddata_version = check_kneaddata_installation()
        if not kneaddata_ok:
            logger.error(f"KneadData not properly installed: {kneaddata_version}")
            return None

    humann4_ok, humann4_version = check_humann4_installation()
    if not humann4_ok:
        logger.error(f"HUMAnN3 not properly installed: {humann4_version}")
        return None

    if skip_kneaddata:
        logger.info("Skipping KneadData processing as requested")
    else:
        logger.info(f"Starting preprocessing pipeline with {len(input_files)} input files")
        logger.info(f"KneadData version: {kneaddata_version}")
    
    logger.info(f"HUMAnN3 version: {humann4_version}")

    # 2. Create output directories
    if kneaddata_output_dir is None:
        kneaddata_output_dir = os.path.join(output_dir, "kneaddata_output")
    if humann4_output_dir is None:
        humann4_output_dir = os.path.join(output_dir, "humann4_output")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(kneaddata_output_dir, exist_ok=True)
    os.makedirs(humann4_output_dir, exist_ok=True)

    logger.info(f"Using KneadData output directory: {kneaddata_output_dir}")
    logger.info(f"Using HUMAnN3 output directory: {humann4_output_dir}")

    # 3. Handle KneadData step
    kneaddata_files = []
    
    if skip_kneaddata:
        # Skip KneadData and work with existing output files
        if kneaddata_output_files:
            # Use explicitly provided files
            logger.info(f"Using {len(kneaddata_output_files)} provided KneadData output files")
            for file in kneaddata_output_files:
                if os.path.exists(file) and os.path.isfile(file):
                    kneaddata_files.append(file)
                else:
                    logger.warning(f"KneadData output file not found: {file}")
        
        elif kneaddata_output_pattern:
            # Use files matching a pattern
            import glob
            
            logger.info(f"Looking for KneadData output files with pattern: {kneaddata_output_pattern}")
            
            # If the pattern includes {sample}, we need to replace it for each potential sample
            if "{sample}" in kneaddata_output_pattern:
                # Extract potential sample names from input files
                sample_names = set()
                for input_file in input_files:
                    sample_name = os.path.basename(input_file).split("_")[0]
                    sample_names.add(sample_name)
                
                # Find files for each sample
                for sample_name in sample_names:
                    pattern = kneaddata_output_pattern.format(sample=sample_name)
                    files = glob.glob(pattern)
                    logger.info(f"Found {len(files)} files for sample {sample_name} with pattern {pattern}")
                    kneaddata_files.extend(files)
            else:
                # Use pattern directly
                files = glob.glob(kneaddata_output_pattern)
                logger.info(f"Found {len(files)} files with pattern {kneaddata_output_pattern}")
                kneaddata_files.extend(files)
        
        else:
            # Look in the kneaddata_output_dir for paired files
            logger.info(f"Looking for KneadData output files in {kneaddata_output_dir}")
            
            for root, _, files in os.walk(kneaddata_output_dir):
                for filename in files:
                    if ("paired_1.fastq" in filename or "paired_2.fastq" in filename or
                        "kneaddata_paired_1.fastq" in filename or "kneaddata_paired_2.fastq" in filename):
                        file_path = os.path.join(root, filename)
                        kneaddata_files.append(file_path)
                        logger.debug(f"Found KneadData file: {filename}")
    else:
        # Run KneadData as usual
        logger.info(f"Starting KneadData step (paired={paired})...")

        if paired:
            # Ensure an even number of input files for paired reads
            if len(input_files) % 2 != 0:
                logger.error("Paired mode requires an even number of input files")
                return None

            # Process input files in pairs
            for i in range(0, len(input_files), 2):
                pair = [input_files[i], input_files[i + 1]]
                sample_name = os.path.basename(pair[0]).split("_")[0]
                sample_outdir = os.path.join(kneaddata_output_dir, sample_name)
                os.makedirs(sample_outdir, exist_ok=True)

                logger.info(f"KneadData on paired files for sample {sample_name}")
                results = run_kneaddata(
                    input_files=pair,
                    output_dir=sample_outdir,
                    threads=threads,
                    reference_dbs=kneaddata_dbs,
                    paired=True,
                    additional_options=kneaddata_options,
                    logger=logger,
                )
                if results:
                    kneaddata_files.extend(results)
                else:
                    logger.warning(f"No KneadData output for sample {sample_name}")

        else:
            # Single-end mode
            for f in input_files:
                sample_name = os.path.basename(f).split("_")[0]
                sample_outdir = os.path.join(kneaddata_output_dir, sample_name)
                os.makedirs(sample_outdir, exist_ok=True)

                logger.info(f"KneadData on single-end file for sample {sample_name}")
                results = run_kneaddata(
                    input_files=[f],
                    output_dir=sample_outdir,
                    threads=threads,
                    reference_dbs=kneaddata_dbs,
                    paired=False,
                    additional_options=kneaddata_options,
                    logger=logger,
                )
                if results:
                    kneaddata_files.extend(results)
                else:
                    logger.warning(f"No KneadData output for sample {sample_name}")

    if not kneaddata_files:
        logger.error("No KneadData output files found; cannot continue.")
        return None

    logger.info(f"Found {len(kneaddata_files)} total KneadData output files")

    # 4. Organize KneadData outputs by sample
    logger.info("Preparing KneadData outputs for HUMAnN3...")
    sample_files = {}

    for f in kneaddata_files:
        fname = os.path.basename(f)
        # Only process .fastq files
        if not fname.endswith(".fastq"):
            continue

        # Identify if file is R1 or R2 based on naming
        if "_paired_1.fastq" in fname:
            sname = fname.replace("_paired_1.fastq", "")
        elif "_paired_2.fastq" in fname:
            sname = fname.replace("_paired_2.fastq", "")
        elif "_kneaddata_paired_1.fastq" in fname:
            sname = fname.replace("_kneaddata_paired_1.fastq", "")
        elif "_kneaddata_paired_2.fastq" in fname:
            sname = fname.replace("_kneaddata_paired_2.fastq", "")
        else:
            logger.debug(f"Skipping non-paired file: {fname}")
            continue

        # Remove trailing "_kneaddata" if it exists
        if sname.endswith("_kneaddata"):
            sname = sname[: -len("_kneaddata")]

        # Use a set to avoid duplicates
        if sname not in sample_files:
            sample_files[sname] = set()
        sample_files[sname].add(f)

    # 5. For each sample, find exactly one _paired_1.fastq and one _paired_2.fastq,
    #    then concatenate them for HUMAnN3.
    humann4_input_files = []

    for sample_name, file_set in sample_files.items():
        flist = sorted(list(file_set))
        logger.info(f"Sample {sample_name} has {len(flist)} candidate files")

        r1_list = []
        r2_list = []
        
        for x in flist:
            basename = os.path.basename(x)
            if ("_paired_1.fastq" in basename) or ("_kneaddata_paired_1.fastq" in basename):
                r1_list.append(x)
            elif ("_paired_2.fastq" in basename) or ("_kneaddata_paired_2.fastq" in basename):
                r2_list.append(x)

        if len(r1_list) == 1 and len(r2_list) == 1:
            r1, r2 = r1_list[0], r2_list[0]
            logger.info(
                f"Concatenating R1 and R2 for sample {sample_name} -> single HUMAnN3 input"
            )
            out_dir = os.path.dirname(r1)
            concat_file = os.path.join(out_dir, f"{sample_name}_paired_concat.fastq")

            try:
                with open(concat_file, "w") as outfile:
                    for rf in (r1, r2):
                        size = os.path.getsize(rf)
                        logger.debug(f"  Adding {os.path.basename(rf)} ({size} bytes)")
                        with open(rf, "r") as infile:
                            outfile.write(infile.read())

                if os.path.getsize(concat_file) > 0:
                    logger.info(f"Created {os.path.basename(concat_file)} for {sample_name}")
                    humann4_input_files.append(concat_file)
                else:
                    logger.error(f"Concat file was empty for {sample_name}")
            except Exception as e:
                logger.error(f"Concat error for {sample_name}: {str(e)}")
        else:
            logger.warning(
                f"Sample {sample_name}: found {len(r1_list)} R1 files and {len(r2_list)} R2 files. Skipping."
            )

    logger.info(f"Prepared {len(humann4_input_files)} input files for HUMAnN3")

    if not humann4_input_files:
        logger.error("No valid input files after KneadData processing; cannot run HUMAnN3.")
        return None

    # 6. Run HUMAnN3
    logger.info("Starting HUMAnN3...")
    humann4_results = run_humann4(
        input_files=humann4_input_files,
        output_dir=humann4_output_dir,
        threads=threads,
        nucleotide_db=nucleotide_db,
        protein_db=protein_db,
        additional_options=humann4_options,
        logger=logger
    )
    
    if not humann4_results:
        logger.error("HUMAnN3 step failed.")
        return None

    logger.info(f"HUMAnN3 completed for {len(humann4_results)} samples.")

    if humann4_results:
        for sample_id, files in humann4_results.items():
            if isinstance(files, dict) and files.get('metaphlan') is None:
                logger.warning(f"MetaPhlAn output file not found for sample {sample_id}. Continuing anyway.")

    # Return combined results
    return {
        "kneaddata_files": kneaddata_files,
        "humann4_results": humann4_results,
    }

def run_preprocessing_pipeline_parallel(input_files, output_dir, threads_per_sample=1, 
                                       max_parallel=None, kneaddata_dbs=None, 
                                       nucleotide_db=None, protein_db=None,
                                       kneaddata_options=None, humann4_options=None, 
                                       paired=False, kneaddata_output_dir=None, humann4_output_dir=None,
                                       skip_kneaddata=False, kneaddata_output_files=None, 
                                       kneaddata_output_pattern="kneaddata_paired", logger=None):
    """
    Run the full preprocessing pipeline in parallel: KneadData → HUMAnN3.
    
    Args:
        input_files: List of input FASTQ files
        output_dir: Base directory for outputs
        threads_per_sample: Number of threads per sample
        max_parallel: Maximum number of parallel samples (None = CPU count)
        kneaddata_dbs: Path(s) to KneadData reference database(s). Can be a string or a list of paths.
        nucleotide_db: Path to HUMAnN3 nucleotide database
        protein_db: Path to HUMAnN3 protein database
        kneaddata_options: Dict of additional KneadData options
        humann4_options: Dict of additional HUMAnN3 options
        paired: Whether input files are paired
        kneaddata_output_dir: Custom directory for KneadData outputs
        humann4_output_dir: Custom directory for HUMAnN3 outputs
        skip_kneaddata: Whether to skip KneadData processing
        kneaddata_output_files: List of existing KneadData output files to use
        kneaddata_output_pattern: Pattern to find KneadData output files
        logger: Logger instance
        
    Returns:
        Dict of final HUMAnN3 output file paths by sample and type
    """
    if logger is None:
        logger = logging.getLogger('humann4_analysis')
    
    # Check installations (skip KneadData check if we're skipping KneadData processing)
    if not skip_kneaddata:
        kneaddata_ok, kneaddata_version = check_kneaddata_installation()
        if not kneaddata_ok:
            logger.error(f"KneadData not properly installed: {kneaddata_version}")
            return None
    
    humann4_ok, humann4_version = check_humann4_installation()
    if not humann4_ok:
        logger.error(f"HUMAnN3 not properly installed: {humann4_version}")
        return None
    
    if skip_kneaddata:
        logger.info("Skipping KneadData processing as requested")
    else:
        logger.info(f"Starting parallel preprocessing pipeline with {len(input_files)} input files")
        logger.info(f"KneadData version: {kneaddata_version}")
    
    logger.info(f"HUMAnN3 version: {humann4_version}")
    
    # Use specified output directories or create defaults
    if kneaddata_output_dir is None:
        kneaddata_output = os.path.join(output_dir, "kneaddata_output")
    else:
        kneaddata_output = kneaddata_output_dir
        
    if humann4_output_dir is None:
        humann4_output = os.path.join(output_dir, "humann4_output")
    else:
        humann4_output = humann4_output_dir
    
    logger.info(f"KneadData output dir: {kneaddata_output}")
    logger.info(f"HUMAnN3 output dir: {humann4_output}")
    
    os.makedirs(kneaddata_output, exist_ok=True)
    os.makedirs(humann4_output, exist_ok=True)
    
    # Handling for skipping KneadData
    if skip_kneaddata:
        # Prepare HUMAnN3 input by finding existing KneadData output files
        humann4_input_files = []
        
        # Case 1: KneadData output files explicitly provided
        if kneaddata_output_files:
            logger.info(f"Using {len(kneaddata_output_files)} provided KneadData output files")
            # Find and group KneadData paired output files by sample
            sample_files = {}
            
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
                        if sample_name not in sample_files:
                            sample_files[sample_name] = []
                        sample_files[sample_name].append(file_path)
                else:
                    logger.warning(f"KneadData output file not found: {file_path}")
        
        # Case 2: KneadData output pattern provided
        elif kneaddata_output_pattern:
            import glob
            
            logger.info(f"Looking for KneadData output files with pattern: {kneaddata_output_pattern}")
            
            # If the pattern includes {sample}, we need to extract sample names from available files
            if "{sample}" in kneaddata_output_pattern:
                # First, find all potential sample names from input files
                sample_names = []
                for input_file in input_files:
                    sample_name = os.path.basename(input_file).split("_")[0]
                    if sample_name not in sample_names:
                        sample_names.append(sample_name)
                
                logger.info(f"Found {len(sample_names)} potential sample names from input files")
                
                # Now find files for each sample
                sample_files = {}
                for sample_name in sample_names:
                    pattern = kneaddata_output_pattern.format(sample=sample_name)
                    files = glob.glob(pattern)
                    if files:
                        sample_files[sample_name] = files
                        logger.debug(f"Found {len(files)} files for sample {sample_name}")
            else:
                # Use the pattern directly
                files = glob.glob(kneaddata_output_pattern)
                logger.info(f"Found {len(files)} files matching pattern")
                
                # Group files by extracting sample names
                sample_files = {}
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
                        if sample_name not in sample_files:
                            sample_files[sample_name] = []
                        sample_files[sample_name].append(file_path)
        
        # Case 3: No specific files or pattern - try to infer from kneaddata_output directory
        else:
            logger.info(f"Looking for KneadData output files in {kneaddata_output}")
            
            # Just find all potential paired files in the directory
            sample_files = {}
            for root, _, files in os.walk(kneaddata_output):
                for filename in files:
                    if ("paired_1.fastq" in filename or "paired_2.fastq" in filename or
                        "kneaddata_paired_1.fastq" in filename or "kneaddata_paired_2.fastq" in filename):
                        
                        file_path = os.path.join(root, filename)
                        
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
                            if sample_name not in sample_files:
                                sample_files[sample_name] = []
                            sample_files[sample_name].append(file_path)
        
        # Now, for each sample with files, prepare for HUMAnN3
        logger.info(f"Found KneadData paired files for {len(sample_files)} samples")
        
        for sample_id, files in sample_files.items():
            logger.info(f"Processing KneadData output for sample {sample_id}")
            
            # First, strictly filter to only paired output files that belong to THIS sample
            paired_files = []
            for file in files:
                filename = os.path.basename(file)
                # Only include files that match both the pattern AND belong to this sample
                if (f"{sample_id}" in filename) and ("_paired_1.fastq" in filename or "_paired_2.fastq" in filename):
                    paired_files.append(file)
                    logger.debug(f"Found paired file for {sample_id}: {filename}")
            
            # Skip samples without exactly 2 paired files
            if len(paired_files) != 2:
                logger.warning(f"Sample {sample_id} has {len(paired_files)} paired files (expected 2). Skipping.")
                continue
            
            # Sort to ensure R1 comes before R2
            paired_files.sort()
            
            # Verify we have proper paired files
            file1 = os.path.basename(paired_files[0])
            file2 = os.path.basename(paired_files[1])
            
            if "_paired_1.fastq" not in file1 or "_paired_2.fastq" not in file2:
                logger.warning(f"Files for sample {sample_id} don't match expected pattern: {file1}, {file2}. Skipping.")
                continue
            
            # Create concatenated file for HUMAnN3
            concat_dir = os.path.dirname(paired_files[0])
            concatenated_file = os.path.join(concat_dir, f"{sample_id}_paired_concat.fastq")
            logger.info(f"Concatenating paired files for sample {sample_id} to {os.path.basename(concatenated_file)}")
            
            try:
                # Create concatenated file
                with open(concatenated_file, 'w') as outfile:
                    for file in paired_files:
                        logger.debug(f"  Adding file: {os.path.basename(file)} (Size: {os.path.getsize(file)} bytes)")
                        with open(file, 'r') as infile:
                            outfile.write(infile.read())
                
                # Verify the concatenated file
                if os.path.exists(concatenated_file) and os.path.getsize(concatenated_file) > 0:
                    logger.info(f"Successfully created concatenated file: {os.path.basename(concatenated_file)} "
                             f"(Size: {os.path.getsize(concatenated_file)} bytes)")
                    humann4_input_files.append(concatenated_file)
                else:
                    logger.error(f"Failed to create valid concatenated file for {sample_id}.")
            except Exception as e:
                logger.error(f"Error concatenating files for sample {sample_id}: {str(e)}")
    else:
        # Run KneadData normally
        logger.info("Starting KneadData step in parallel...")
        kneaddata_results = run_kneaddata_parallel(
            input_files=input_files,
            output_dir=kneaddata_output,
            threads=threads_per_sample,
            max_parallel=max_parallel,
            reference_dbs=kneaddata_dbs,
            paired=paired,
            additional_options=kneaddata_options,
            logger=logger
        )
        
        if not kneaddata_results:
            logger.error("KneadData step failed, stopping pipeline")
            return None
        
        # Prepare files for HUMAnN3
        logger.info("Preparing KneadData outputs for HUMAnN3...")
        sample_files = {}
        humann4_input_files = []

        # Process each sample's output files
        for sample_id, files in kneaddata_results.items():
            logger.info(f"Processing KneadData output for sample {sample_id}")
            
            if not isinstance(files, list) or not files:
                logger.warning(f"No KneadData output files found for sample {sample_id}")
                continue
                
            # First, strictly filter to only paired output files that belong to THIS sample
            paired_files = []
            for file in files:
                filename = os.path.basename(file)
                # Only include files that match both the pattern AND belong to this sample
                if (f"{sample_id}" in filename) and ("_paired_1.fastq" in filename or "_paired_2.fastq" in filename):
                    paired_files.append(file)
                    logger.debug(f"Found paired file for {sample_id}: {filename}")
            
            # Skip samples without exactly 2 paired files
            if len(paired_files) != 2:
                logger.warning(f"Sample {sample_id} has {len(paired_files)} paired files (expected 2). Skipping.")
                continue
            
            # Sort to ensure R1 comes before R2
            paired_files.sort()
            
            # Verify we have proper paired files
            file1 = os.path.basename(paired_files[0])
            file2 = os.path.basename(paired_files[1])
            
            if "_paired_1.fastq" not in file1 or "_paired_2.fastq" not in file2:
                logger.warning(f"Files for sample {sample_id} don't match expected pattern: {file1}, {file2}. Skipping.")
                continue
            
            # Create concatenated file for HUMAnN3
            concat_dir = os.path.dirname(paired_files[0])
            concatenated_file = os.path.join(concat_dir, f"{sample_id}_paired_concat.fastq")
            logger.info(f"Concatenating paired files for sample {sample_id} to {os.path.basename(concatenated_file)}")
            
            try:
                # Create concatenated file
                with open(concatenated_file, 'w') as outfile:
                    for file in paired_files:
                        logger.debug(f"  Adding file: {os.path.basename(file)} (Size: {os.path.getsize(file)} bytes)")
                        with open(file, 'r') as infile:
                            outfile.write(infile.read())
                
                # Verify the concatenated file
                if os.path.exists(concatenated_file) and os.path.getsize(concatenated_file) > 0:
                    logger.info(f"Successfully created concatenated file: {os.path.basename(concatenated_file)} "
                             f"(Size: {os.path.getsize(concatenated_file)} bytes)")
                    humann4_input_files.append(concatenated_file)
                else:
                    logger.error(f"Failed to create valid concatenated file for {sample_id}.")
            except Exception as e:
                logger.error(f"Error concatenating files for sample {sample_id}: {str(e)}")

    logger.info(f"Prepared {len(humann4_input_files)} input files for HUMAnN3")

    # Now check that we have valid files to run HUMAnN3 on
    if not humann4_input_files:
        logger.error("No valid input files prepared for HUMAnN3")
        return None

    # Step 2: Run HUMAnN3 in parallel only on our concatenated files
    logger.info("Starting HUMAnN3 step in parallel...")
    humann4_results = run_humann4_parallel(
        input_files=humann4_input_files,
        output_dir=humann4_output,
        threads=threads_per_sample,
        max_parallel=max_parallel,
        nucleotide_db=nucleotide_db,
        protein_db=protein_db,
        additional_options=humann4_options,
        logger=logger
    )
    
    if not humann4_results:
        logger.error("HUMAnN3 step failed")
        return None
    
    logger.info(f"HUMAnN3 completed for {len(humann4_results)} samples")
    
    # Return combined results
    return {
        'kneaddata_files': humann4_input_files,  
        'humann4_results': humann4_results
    }
