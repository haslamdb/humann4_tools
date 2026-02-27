# humann4_tools/core/humann4.py
"""
Main module for humann4_tools package.

This module provides functions to run the full analysis pipeline
or individual components as needed.
"""
import os
import logging
import time
import traceback
import pandas as pd

from src.humann4_tools.logger import setup_logger, log_print
from src.humann4_tools.utils.sample_utils import validate_sample_key, check_input_files_exist
from src.humann4_tools.utils.file_utils import check_file_exists_with_logger
from src.humann4_tools.humann4.pathway_processing import process_pathway_abundance
from src.humann4_tools.humann4.gene_processing import process_gene_families
from src.humann4_tools.analysis.metadata import read_and_process_metadata
from src.humann4_tools.analysis.visualizations import (
    read_and_process_gene_families,
    read_and_process_pathways,
)
from src.humann4_tools.analysis.statistical import run_statistical_tests
from src.humann4_tools.analysis.differential_abundance import run_differential_abundance_analysis
from src.humann4_tools.preprocessing.pipeline import run_preprocessing_pipeline


def run_full_pipeline(
    sample_key,
    pathway_dir,
    gene_dir,
    output_dir,
    output_prefix="processed_files",
    group_col="Group",
    skip_pathway=False,
    skip_gene=False,
    skip_downstream=False,
    no_interactive=False,
    log_file=None,
):
    """
    Run the full HUMAnN3 processing and analysis pipeline.
    
    Args:
        sample_key: Path to CSV file with sample metadata
        pathway_dir: Directory containing HUMAnN3 pathway files
        gene_dir: Directory containing HUMAnN3 gene family files
        output_dir: Directory for output files
        output_prefix: Prefix for output files
        group_col: Column name to use for statistical grouping
        skip_pathway: Skip pathway processing
        skip_gene: Skip gene family processing
        skip_downstream: Skip downstream analysis
        no_interactive: Disable interactive prompts
        log_file: Path to log file
        
    Returns:
        Tuple of (pathway_file, gene_file, success_flag)
    """
    # Setup logging
    logger = setup_logger(log_file=log_file)
    log_print("Starting HUMAnN3 Analysis Pipeline", level="info")
    start_time = time.time()

    # Process sample metadata
    samples, selected_columns = validate_sample_key(sample_key, no_interactive=no_interactive)

    # Check input files
    valid_path_samples, valid_gene_samples = check_input_files_exist(samples, pathway_dir, gene_dir)

    # Process pathways
    pathway_unstrat_file = None
    if not skip_pathway and valid_path_samples:
        pathway_unstrat_file = process_pathway_abundance(
            valid_path_samples,
            pathway_dir,
            output_dir,
            output_prefix,
            selected_columns=selected_columns
        )
    else:
        if skip_pathway:
            log_print("Skipping HUMAnN3 pathway processing", level="info")
        else:
            log_print("No valid pathway files; skipping pathway stage", level="warning")

    # Process gene families
    gene_unstrat_file = None
    if not skip_gene and valid_gene_samples:
        gene_unstrat_file = process_gene_families(
            valid_gene_samples,
            gene_dir,
            output_dir,
            output_prefix,
            selected_columns=selected_columns
        )
    else:
        if skip_gene:
            log_print("Skipping HUMAnN3 gene processing", level="info")
        else:
            log_print("No valid gene files; skipping gene stage", level="warning")

    # Downstream analysis
    if skip_downstream:
        log_print("Skipping downstream analysis stage", level="info")
        success = True
    else:
        if not pathway_unstrat_file and not gene_unstrat_file:
            log_print("No unstratified HUMAnN3 outputs found; cannot run downstream analysis", level="warning")
            success = False
        else:
            try:
                # Create analysis output directory
                downstream_out = os.path.join(output_dir, "DownstreamAnalysis")
                os.makedirs(downstream_out, exist_ok=True)
                logger.info(f"Downstream analysis output will be in: {downstream_out}")

                # Read sample metadata
                if not check_file_exists_with_logger(sample_key, "Sample key", logger):
                    log_print("Cannot proceed with downstream analysis; missing sample key", level="error")
                    success = False
                else:
                    sample_key_df = read_and_process_metadata(sample_key, logger)

                    # Process gene families if available
                    if gene_unstrat_file and check_file_exists_with_logger(gene_unstrat_file, "Gene families", logger):
                        read_and_process_gene_families(gene_unstrat_file, sample_key_df, downstream_out, logger)

                    # Process pathways if available
                    if pathway_unstrat_file and check_file_exists_with_logger(pathway_unstrat_file, "Pathways", logger):
                        pathways_merged = read_and_process_pathways(
                            pathway_unstrat_file,
                            sample_key_df,
                            downstream_out,
                            logger
                        )
                        # Run statistical tests
                        run_statistical_tests(pathways_merged, downstream_out, logger, group_col=group_col)

                    success = True
            except Exception as e:
                logger.error(f"Downstream analysis failed: {e}")
                logger.error(traceback.format_exc())
                success = False

    elapsed = time.time() - start_time
    hh, rr = divmod(elapsed, 3600)
    mm, ss = divmod(rr, 60)
    log_print(f"Pipeline finished in {int(hh)}h {int(mm)}m {int(ss)}s", level="info")

    return pathway_unstrat_file, gene_unstrat_file, success


# Preprocessing and Humann3 run
def run_preprocessing_and_analysis(
    input_fastq,
    sample_key,
    output_dir,
    output_prefix="ProcessedFiles",
    group_col="Group",
    threads=1,
    kneaddata_db=None,
    nucleotide_db=None,
    protein_db=None,
    paired=False,
    skip_pathway=False,
    skip_gene=False,
    skip_downstream=False,
    log_file=None
):
    """
    Run the full preprocessing and analysis pipeline:
    KneadData → HUMAnN3 → HUMAnN3 processing → Downstream analysis.
    
    Args:
        input_fastq: List of input FASTQ files
        sample_key: Path to sample metadata CSV
        output_dir: Directory for output files
        output_prefix: Prefix for output files
        group_col: Column name for statistical grouping
        threads: Number of CPU threads to use
        kneaddata_db: Path to KneadData reference database
        nucleotide_db: Path to HUMAnN3 nucleotide database
        protein_db: Path to HUMAnN3 protein database
        paired: Whether input FASTQ files are paired-end
        skip_pathway: Skip pathway processing
        skip_gene: Skip gene family processing
        skip_downstream: Skip downstream analysis
        log_file: Path to log file
        
    Returns:
        Tuple of (pathway_file, gene_file, success_flag)
    """
    # Setup logging
    logger = setup_logger(log_file=log_file)
    log_print("Starting Full Microbiome Analysis Pipeline", level="info")
    start_time = time.time()
    
    # Create preprocessing output directory
    preproc_dir = os.path.join(output_dir, "processed_files")
    os.makedirs(preproc_dir, exist_ok=True)
    
    # Step 1: Run preprocessing (KneadData + HUMAnN3)
    preprocessing_results = run_preprocessing_pipeline(
        input_files=input_fastq,
        output_dir=preproc_dir,
        threads=threads,
        kneaddata_db=kneaddata_db,
        nucleotide_db=nucleotide_db,
        protein_db=protein_db,
        paired=paired,
        logger=logger
    )
    
    if not preprocessing_results:
        log_print("Preprocessing pipeline failed", level="error")
        return None, None, False
    
    # Extract paths to HUMAnN3 output files
    humann4_results = preprocessing_results['humann4_results']
    pathway_files = []
    gene_files = []
    
    for sample, files in humann4_results.items():
        if files.get('pathabundance'):
            pathway_files.append((sample, files['pathabundance']))
        if files.get('genefamilies'):
            gene_files.append((sample, files['genefamilies']))
    
    # Step 2: Process sample metadata
    samples, selected_columns = validate_sample_key(sample_key, no_interactive=True)
    
    # Step 3: Process pathway files
    pathway_unstrat_file = None
    if not skip_pathway and pathway_files:
        pathway_unstrat_file = process_pathway_abundance(
            pathway_files,
            preproc_dir,
            output_dir,
            output_prefix,
            selected_columns=selected_columns
        )
    
    # Step 4: Process gene family files
    gene_unstrat_file = None
    if not skip_gene and gene_files:
        gene_unstrat_file = process_gene_families(
            gene_files,
            preproc_dir,
            output_dir,
            output_prefix,
            selected_columns=selected_columns
        )
    
    # Step 5: Downstream analysis
    success = True
    if not skip_downstream:
        if not pathway_unstrat_file and not gene_unstrat_file:
            log_print("No unstratified HUMAnN3 outputs found; cannot run downstream analysis", level="warning")
            success = False
        else:
            try:
                # Create downstream analysis directory
                downstream_out = os.path.join(output_dir, "DownstreamAnalysis")
                os.makedirs(downstream_out, exist_ok=True)
                
                # Read sample metadata
                sample_key_df = read_and_process_metadata(sample_key, logger)
                
                # Process gene families if available
                if gene_unstrat_file:
                    read_and_process_gene_families(gene_unstrat_file, sample_key_df, downstream_out, logger)
                
                # Process pathways if available
                if pathway_unstrat_file:
                    pathways_merged = read_and_process_pathways(pathway_unstrat_file, sample_key_df, downstream_out, logger)
                    run_statistical_tests(pathways_merged, downstream_out, logger, group_col=group_col)
            except Exception as e:
                logger.error(f"Downstream analysis failed: {e}")
                logger.error(traceback.format_exc())
                success = False
    
    elapsed = time.time() - start_time
    hh, rr = divmod(elapsed, 3600)
    mm, ss = divmod(rr, 60)
    log_print(f"Pipeline finished in {int(hh)}h {int(mm)}m {int(ss)}s", level="info")
    
    return pathway_unstrat_file, gene_unstrat_file, success


def process_humann4_files_only(
    sample_key,
    pathway_dir,
    gene_dir,
    output_dir,
    output_prefix="ProcessedFiles",
    skip_pathway=False,
    skip_gene=False,
    no_interactive=False,
    log_file=None,
):
    """
    Run only the HUMAnN3 processing stage (no downstream analysis).
    
    Args:
        Same as run_full_pipeline but without downstream parameters
        
    Returns:
        Tuple of (pathway_file, gene_file)
    """
    # Run full pipeline with skip_downstream=True
    pathway_unstrat_file, gene_unstrat_file, _ = run_full_pipeline(
        sample_key=sample_key,
        pathway_dir=pathway_dir,
        gene_dir=gene_dir,
        output_dir=output_dir,
        output_prefix=output_prefix,
        skip_pathway=skip_pathway,
        skip_gene=skip_gene,
        skip_downstream=True,
        no_interactive=no_interactive,
        log_file=log_file
    )
    
    return pathway_unstrat_file, gene_unstrat_file


def analyze_existing_humann4_files(
    pathway_file,
    gene_file,
    sample_key,
    output_dir,
    group_col="Group",
    log_file=None,
):
    """
    Run only the downstream analysis on existing HUMAnN3 unstratified files.
    
    Args:
        pathway_file: Path to unstratified pathway file
        gene_file: Path to unstratified gene family file
        sample_key: Path to sample metadata CSV
        output_dir: Directory for output files
        group_col: Column name to use for statistical grouping
        log_file: Path to log file
        
    Returns:
        Boolean success flag
    """
    # Setup logging
    logger = setup_logger(log_file=log_file)
    log_print("Starting downstream analysis of existing HUMAnN3 files", level="info")
    start_time = time.time()

    if not pathway_file and not gene_file:
        log_print("Error: At least one of pathway_file or gene_file must be provided", level="error")
        return False

    try:
        # Create analysis output directory
        downstream_out = os.path.join(output_dir, "DownstreamAnalysis")
        os.makedirs(downstream_out, exist_ok=True)
        logger.info(f"Downstream analysis output will be in: {downstream_out}")

        # Read sample metadata
        if not check_file_exists_with_logger(sample_key, "Sample key", logger):
            log_print("Cannot proceed with analysis; missing sample key", level="error")
            return False

        sample_key_df = read_and_process_metadata(sample_key, logger)

        # Process gene families if available
        if gene_file and check_file_exists_with_logger(gene_file, "Gene families", logger):
            read_and_process_gene_families(gene_file, sample_key_df, downstream_out, logger)

        # Process pathways if available
        if pathway_file and check_file_exists_with_logger(pathway_file, "Pathways", logger):
            pathways_merged = read_and_process_pathways(pathway_file, sample_key_df, downstream_out, logger)
            # Run statistical tests
            run_statistical_tests(pathways_merged, downstream_out, logger, group_col=group_col)

        elapsed = time.time() - start_time
        hh, rr = divmod(elapsed, 3600)
        mm, ss = divmod(rr, 60)
        log_print(f"Analysis finished in {int(hh)}h {int(mm)}m {int(ss)}s", level="info")

        return True

    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        logger.error(traceback.format_exc())
        return False


def run_pathway_differential_abundance(
    pathway_file,
    sample_key,
    output_dir,
    group_col="Group",
    methods=["aldex2", "ancom", "ancom-bc"],
    include_unmapped=True,
    log_file=None,
):
    """
    Run differential abundance tests on pathway data.
    
    Args:
        pathway_file: Path to unstratified pathway file
        sample_key: Path to sample metadata
        output_dir: Directory for output files
        group_col: Column name for grouping
        methods: List of methods to run
        include_unmapped: Whether to include unmapped features
        log_file: Path to log file
        
    Returns:
        Dictionary with results from each method
    """
    logger = logging.getLogger("humann4_analysis")

    if not check_file_exists_with_logger(pathway_file, "Pathway file", logger):
        logger.error("Cannot run differential abundance analysis: pathway file not found")
        return None

    if not check_file_exists_with_logger(sample_key, "Sample key", logger):
        logger.error("Cannot run differential abundance analysis: sample key not found")
        return None

    try:
        # Make output directory
        diff_abund_dir = os.path.join(output_dir, "DifferentialAbundance", "Pathways")
        os.makedirs(diff_abund_dir, exist_ok=True)

        # Read data
        pathway_df = pd.read_csv(pathway_file, sep="\t", index_col=0)
        metadata_df = pd.read_csv(sample_key, index_col=None)

        # Get sample ID column (attempt common naming)
        sample_id_col = None
        common_id_names = ["SampleName", "Sample", "SampleID", "Sample_ID", "sample_name", "sample_id"]
        for col in common_id_names:
            if col in metadata_df.columns:
                sample_id_col = col
                break

        if sample_id_col is None:
            logger.error("Could not find sample ID column in metadata")
            return None

        # Set the sample ID as index
        metadata_df = metadata_df.set_index(sample_id_col)

        # Decide how to handle unmapped features
        denom = "all" if include_unmapped else "unmapped_excluded"

        # Run differential abundance analysis
        results = run_differential_abundance_analysis(
            pathway_df,
            metadata_df,
            diff_abund_dir,
            group_col=group_col,
            methods=methods,
            denom=denom,
            logger=logger
        )

        return results

    except Exception as e:
        logger.error(f"Error in differential abundance analysis: {str(e)}")
        logger.error(traceback.format_exc())
        return None


def run_gene_differential_abundance(
    gene_file,
    sample_key,
    output_dir,
    group_col="Group",
    methods=["aldex2", "ancom", "ancom-bc"],
    include_unmapped=True,
    log_file=None,
):
    """
    Run differential abundance tests on gene family data.
    
    Args:
        gene_file: Path to unstratified gene family file
        sample_key: Path to sample metadata
        output_dir: Directory for output files
        group_col: Column name for grouping
        methods: List of methods to run
        include_unmapped: Whether to include unmapped features
        log_file: Path to log file
        
    Returns:
        Dictionary with results from each method
    """
    logger = logging.getLogger("humann4_analysis")

    if not check_file_exists_with_logger(gene_file, "Gene family file", logger):
        logger.error("Cannot run differential abundance analysis: gene family file not found")
        return None

    if not check_file_exists_with_logger(sample_key, "Sample key", logger):
        logger.error("Cannot run differential abundance analysis: sample key not found")
        return None

    try:
        # Make output directory
        diff_abund_dir = os.path.join(output_dir, "DifferentialAbundance", "Genes")
        os.makedirs(diff_abund_dir, exist_ok=True)

        # Read data
        gene_df = pd.read_csv(gene_file, sep="\t", index_col=0)
        metadata_df = pd.read_csv(sample_key, index_col=None)

        # Get sample ID column (attempt common naming)
        sample_id_col = None
        common_id_names = ["SampleName", "Sample", "SampleID", "Sample_ID", "sample_name", "sample_id"]
        for col in common_id_names:
            if col in metadata_df.columns:
                sample_id_col = col
                break

        if sample_id_col is None:
            logger.error("Could not find sample ID column in metadata")
            return None

        # Set the sample ID as index
        metadata_df = metadata_df.set_index(sample_id_col)

        # Decide how to handle unmapped features
        denom = "all" if include_unmapped else "unmapped_excluded"

        # Run differential abundance analysis
        results = run_differential_abundance_analysis(
            gene_df,
            metadata_df,
            diff_abund_dir,
            group_col=group_col,
            methods=methods,
            denom=denom,
            logger=logger
        )

        return results

    except Exception as e:
        logger.error(f"Error in differential abundance analysis: {str(e)}")
        logger.error(traceback.format_exc())
        return None