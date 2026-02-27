
import os
import sys
import time
import argparse
import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Union
from scipy.stats import kruskal

try:
    import scikit_posthocs as sp
except ImportError:
    sp = None

# Import internal modules
try:
    from src.humann4_tools.utils.resource_utils import track_peak_memory
    from src.humann4_tools.utils.file_utils import sanitize_filename
except ImportError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
    from src.humann4_tools.utils.resource_utils import track_peak_memory
    from src.humann4_tools.utils.file_utils import sanitize_filename

# Set up logging
logger = logging.getLogger('humann4_tools')

def read_and_process_data(
    abundance_file: str, 
    metadata_file: str, 
    sample_id_col: Optional[str] = None, 
    group_col: str = "Group",
    feature_type: str = "pathway"
) -> Tuple[pd.DataFrame, List[str], str, str]:
    """
    Read abundance and metadata files and merge them into a long-format DataFrame.
    
    Args:
        abundance_file: Path to abundance file (unstratified)
        metadata_file: Path to metadata CSV file
        sample_id_col: Column in metadata for sample IDs (auto-detected if None)
        group_col: Column in metadata for grouping
        feature_type: Type of features in abundance file ("pathway" or "gene")
        
    Returns:
        Tuple of (merged_long_df, groups, feature_col, sample_id_col)
    """
    logger.info(f"Reading abundance file: {abundance_file}")
    
    # Read abundance file (make sure it's tab-delimited and has row names in first column)
    try:
        abundance_df = pd.read_csv(abundance_file, sep='\t', index_col=0)
        logger.info(f"Loaded abundance data with {abundance_df.shape[0]} features and {abundance_df.shape[1]} samples")
    except Exception as e:
        logger.error(f"Error reading abundance file: {str(e)}")
        return pd.DataFrame(), [], "", ""
    
    logger.info(f"Reading metadata file: {metadata_file}")
    
    # Read metadata
    try:
        metadata_df = pd.read_csv(metadata_file)
        logger.info(f"Loaded metadata with {metadata_df.shape[0]} samples and {metadata_df.shape[1]} columns")
    except Exception as e:
        logger.error(f"Error reading metadata file: {str(e)}")
        return pd.DataFrame(), [], "", ""
    
    # Auto-detect sample ID column if not specified
    if not sample_id_col:
        common_id_cols = ["SampleName", "Sample", "SampleID", "Sample_ID", "sample_name", "sample_id"]
        for col in common_id_cols:
            if col in metadata_df.columns:
                sample_id_col = col
                logger.info(f"Auto-detected sample ID column: {sample_id_col}")
                break
        
        if not sample_id_col:
            sample_id_col = metadata_df.columns[0]
            logger.warning(f"Could not auto-detect sample ID column, using the first column: {sample_id_col}")
    
    # Verify sample ID column exists
    if sample_id_col not in metadata_df.columns:
        logger.error(f"Sample ID column '{sample_id_col}' not found in metadata")
        return pd.DataFrame(), [], "", ""
    
    # Verify group column exists
    if group_col not in metadata_df.columns:
        logger.error(f"Group column '{group_col}' not found in metadata")
        return pd.DataFrame(), [], "", ""
    
    # Convert abundance data to long format
    feature_col = "Pathway" if feature_type == "pathway" else "Gene_Family"
    long_df = abundance_df.reset_index().melt(
        id_vars=abundance_df.index.name, 
        var_name=sample_id_col, 
        value_name="Abundance"
    )
    
    # Rename feature column
    long_df = long_df.rename(columns={abundance_df.index.name: feature_col})
    
    # Merge with metadata
    merged_df = pd.merge(
        long_df, 
        metadata_df,
        on=sample_id_col,
        how="inner"
    )
    
    # Check if merge was successful
    if merged_df.empty:
        logger.error("No matching samples between abundance data and metadata")
        return pd.DataFrame(), [], "", ""
    
    logger.info(f"Successfully merged abundance data with metadata. Working with {len(merged_df)} rows.")
    
    # Get unique groups
    groups = merged_df[group_col].unique().tolist()
    logger.info(f"Found {len(groups)} groups: {groups}")
    
    return merged_df, groups, feature_col, sample_id_col

def run_statistical_tests(
    abundance_file: str,
    metadata_file: str,
    output_dir: str,
    group_col: str = "Group",
    feature_type: str = "pathway",
    sample_id_col: Optional[str] = None,
    alpha: float = 0.05
) -> bool:
    """
    Run statistical tests on HUMAnN3 output data.
    
    Args:
        abundance_file: Path to abundance file (unstratified)
        metadata_file: Path to metadata CSV file
        output_dir: Directory for output files
        group_col: Column in metadata for grouping
        feature_type: Type of features in abundance file ("pathway" or "gene")
        sample_id_col: Column in metadata for sample IDs (auto-detected if None)
        alpha: Significance threshold
        
    Returns:
        Boolean indicating success or failure
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Read and process data
    merged_df, groups, feature_col, detected_sample_id_col = read_and_process_data(
        abundance_file, metadata_file, sample_id_col, group_col, feature_type
    )
    
    if merged_df.empty:
        logger.error("Failed to process input data")
        return False
    
    # Check number of groups
    if len(groups) < 2:
        logger.error(f"Need at least 2 groups for statistical testing, found {len(groups)}")
        return False
    
    # Run Kruskal-Wallis and Dunn's test
    kw_results, dunn_results = kruskal_wallis_dunn(
        df_long=merged_df,
        group_col=group_col,
        feature_col=feature_col,
        abundance_col="Abundance",
        alpha=alpha
    )
    
    # Save results
    if kw_results.empty:
        logger.error("No valid Kruskal-Wallis results produced")
        return False
    
    # Save Kruskal-Wallis results
    kw_path = os.path.join(output_dir, "kruskal_wallis_results.csv")
    kw_results.to_csv(kw_path, index=False)
    logger.info(f"Saved Kruskal-Wallis results to {kw_path}")
    
    # Save Dunn's post-hoc test results
    if dunn_results:
        dunn_dir = os.path.join(output_dir, "dunn_posthoc_tests")
        os.makedirs(dunn_dir, exist_ok=True)
        
        for feat, pdf in dunn_results.items():
            feat_safe = sanitize_filename(feat)
            dunn_path = os.path.join(dunn_dir, f"dunn_{feat_safe}.csv")
            pdf.to_csv(dunn_path)
        
        logger.info(f"Saved Dunn's post-hoc results for {len(dunn_results)} features in {dunn_dir}")
    elif sp is None:
        logger.warning("Dunn's post-hoc tests were skipped (scikit_posthocs not available)")
    else:
        logger.warning("No Dunn's post-hoc results to save")
    
    # Create a summary file
    summary_path = os.path.join(output_dir, "statistical_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(f"Statistical Analysis Summary\n")
        f.write(f"===========================\n\n")
        f.write(f"Abundance File: {os.path.basename(abundance_file)}\n")
        f.write(f"Metadata File: {os.path.basename(metadata_file)}\n")
        f.write(f"Feature Type: {feature_type}\n")
        f.write(f"Group Column: {group_col}\n")
        f.write(f"Groups: {', '.join(groups)}\n\n")
        
        f.write(f"Kruskal-Wallis Test Results\n")
        f.write(f"-------------------------\n")
        f.write(f"Total Features Tested: {len(kw_results)}\n")
        f.write(f"Significant Features (q < {alpha}): {sum(kw_results['Reject_H0'])}\n\n")
        
        if dunn_results:
            f.write(f"Dunn's Post-hoc Tests\n")
            f.write(f"--------------------\n")
            f.write(f"Features with Post-hoc Tests: {len(dunn_results)}\n")
    
    logger.info(f"Saved statistical summary to {summary_path}")
    return True

def parse_args():
    """Parse command line arguments for the Statistical Testing module."""
    parser = argparse.ArgumentParser(
        description="Run statistical tests on HUMAnN3 output data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Description:
  Statistics module performs non-parametric statistical tests to identify significant
  differences in pathway or gene abundances between sample groups. It uses Kruskal-Wallis
  tests followed by Dunn's post-hoc tests for pairwise comparisons.

Key Features:
  • Performs Kruskal-Wallis tests to detect overall differences between groups
  • Applies Benjamini-Hochberg FDR correction for multiple testing
  • Conducts Dunn's post-hoc tests for significant features to identify specific group differences
  • Works with both pathway and gene family data

Output Files:
  • kruskal_wallis_results.csv: Contains test statistics and adjusted p-values for all features
  • dunn_posthoc_tests/: Directory with pairwise comparison results for significant features
  • statistical_summary.txt: Summary of the analysis, including counts of significant features

Common Usage:
  # Basic statistical testing on pathway data:
  humann4-tools stats --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv --metadata-file metadata.csv

  # Specify a different group column in metadata:
  humann4-tools stats --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv --metadata-file metadata.csv --group-col Treatment

  # With gene family data:
  humann4-tools stats --abundance-file joined_output/genefamilies_cpm_unstratified.tsv --metadata-file metadata.csv --feature-type gene
"""
    )
    
    # Required arguments
    parser.add_argument("--abundance-file", required=True, 
                      help="Path to the unstratified abundance file (pathway or gene family)")
    parser.add_argument("--metadata-file", required=True,
                      help="Path to sample metadata CSV file")
    
    # Analysis options
    parser.add_argument("--output-dir", default="./StatisticalTests",
                      help="Directory for output files")
    parser.add_argument("--feature-type", choices=["pathway", "gene"], default="pathway",
                      help="Type of features in the abundance file (pathway or gene)")
    parser.add_argument("--group-col", default="Group",
                      help="Column name in metadata for grouping samples")
    parser.add_argument("--sample-id-col", 
                      help="Column name in metadata for sample IDs (autodetected if not specified)")
    parser.add_argument("--alpha", type=float, default=0.05,
                      help="Significance threshold for statistical tests (default: 0.05)")
    
    # Logging options
    parser.add_argument("--log-file", 
                      help="Path to log file")
    parser.add_argument("--log-level", default="INFO", 
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                      help="Logging level")
    
    return parser.parse_args()

@track_peak_memory
def main():
    """Main function to run statistical tests."""
    # Parse arguments
    args = parse_args()
    
    # Setup logging
    log_level = getattr(logging, args.log_level.upper())
    setup_logger(args.log_file, log_level)
    
    logger.info("Starting HUMAnN3 Tools Statistical Testing Module")
    start_time = time.time()
    
    # Check if required Python packages are available
    if sp is None:
        logger.warning("scikit_posthocs is not installed. Dunn's post-hoc tests will be skipped.")
    
    # Check if files exist
    for file_path, desc in [(args.abundance_file, "Abundance file"), (args.metadata_file, "Metadata file")]:
        if not os.path.exists(file_path):
            logger.error(f"ERROR: {desc} not found: {file_path}")
            return 1
    
    # Run statistical tests
    success = run_statistical_tests(
        abundance_file=args.abundance_file,
        metadata_file=args.metadata_file,
        output_dir=args.output_dir,
        group_col=args.group_col,
        feature_type=args.feature_type,
        sample_id_col=args.sample_id_col,
        alpha=args.alpha
    )
    
    if not success:
        logger.error("Statistical testing failed")
        return 1
    
    # Print elapsed time
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    
    logger.info(f"Total processing time: {int(minutes)}m {int(seconds)}s")
    
    # Print next steps
    logger.info("\nNext Steps:")
    logger.info("  For differential abundance analysis:")
    logger.info(f"  humann4-tools diff --abundance-file {args.abundance_file} --metadata-file {args.metadata_file}")
    logger.info("\n  For visualizations:")
    logger.info(f"  humann4-tools viz --abundance-file {args.abundance_file} --metadata-file {args.metadata_file}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())#!/usr/bin/env python3
"""
HUMAnN3 Tools Statistical Testing Module

This module performs statistical tests on HUMAnN3 output data, primarily
Kruskal-Wallis tests followed by Dunn's post-hoc tests for multiple groups.

Examples:
  # Basic usage:
  humann4-tools stats --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv \
                    --metadata-file metadata.csv \
                    --output-dir ./statistical_results
  
  # Specify group column:
  humann4-tools stats --abundance-file joined_output/gene_families_cpm_unstratified.tsv \
                    --metadata-file metadata.csv \
                    --group-col Treatment \
                    --output-dir ./statistical_results
                    
  # Specify feature type:
  humann4-tools stats --abundance-file joined_output/gene_families_cpm_unstratified.tsv \
                    --metadata-file metadata.csv \
                    --feature-type gene \
                    --output-dir ./statistical_results
"""

import os
import sys
import time
import argparse
import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Union
from scipy.stats import kruskal

try:
    import scikit_posthocs as sp
except ImportError:
    sp = None

# Import internal modules
try:
    from src.humann4_tools.utils.resource_utils import track_peak_memory
    from src.humann4_tools.utils.file_utils import sanitize_filename
except ImportError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
    from src.humann4_tools.utils.resource_utils import track_peak_memory
    from src.humann4_tools.utils.file_utils import sanitize_filename

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

def kruskal_wallis_dunn(
    df_long: pd.DataFrame, 
    group_col: str = "Group", 
    feature_col: str = "Pathway", 
    abundance_col: str = "Abundance", 
    alpha: float = 0.05
) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """
    Perform Kruskal-Wallis tests followed by Dunn's post-hoc tests.
    
    Steps:
    1. Kruskal-Wallis across multiple groups
    2. Adjust p-values (Benjamini–Hochberg)
    3. Dunn's post-hoc for significant features
    
    Args:
        df_long: Long-format DataFrame with samples, features, and abundances
        group_col: Column name for grouping variable
        feature_col: Column name for feature (pathway or gene)
        abundance_col: Column name for abundance values
        alpha: Significance threshold
        
    Returns:
        Tuple of (kw_results_df, dict_of_posthoc_dfs)
    """
    # Check if scikit_posthocs is available for Dunn's test
    if sp is None:
        logger.warning("scikit_posthocs not available; will skip Dunn's post-hoc tests")
    
    # Get unique features
    features = df_long[feature_col].unique()
    logger.info(f"Running Kruskal-Wallis tests on {len(features)} features")
    
    # Run Kruskal-Wallis test for each feature
    results = []
    
    from statsmodels.stats.multitest import multipletests
    
    for i, feat in enumerate(features):
        if i % 100 == 0:
            logger.debug(f"Processing feature {i+1}/{len(features)}")
            
        sub = df_long[df_long[feature_col] == feat]
        groups = sub[group_col].unique()
        
        # Skip features with too few groups
        if len(groups) < 2:
            logger.debug(f"Skipping feature {feat}: only {len(groups)} group(s)")
            continue
        
        # Get data for each group
        group_data = [sub.loc[sub[group_col] == g, abundance_col].values for g in groups]
        
        # Skip if any group has too few samples
        if any(len(x) < 2 for x in group_data):
            groups_too_small = [groups[i] for i, x in enumerate(group_data) if len(x) < 2]
            logger.debug(f"Skipping feature {feat}: groups {groups_too_small} have fewer than 2 samples")
            continue
        
        try:
            # Run Kruskal-Wallis test
            stat, pval = kruskal(*group_data)
            
            results.append({
                feature_col: feat,
                "KW_stat": stat,
                "KW_pvalue": pval,
                "Group_count": len(groups),
                "Groups": ",".join(groups)
            })
        except Exception as e:
            logger.warning(f"Error in Kruskal-Wallis test for feature {feat}: {str(e)}")
    
    # Check if we have any results
    if not results:
        logger.warning("No valid Kruskal-Wallis test results")
        return pd.DataFrame(), {}
    
    # Create results DataFrame
    kw_df = pd.DataFrame(results)
    
    # Apply multiple testing correction
    reject, pvals_corrected, _, _ = multipletests(kw_df["KW_pvalue"], alpha=alpha, method="fdr_bh")
    kw_df["KW_padj"] = pvals_corrected
    kw_df["Reject_H0"] = reject
    
    # Count significant features
    sig_count = sum(kw_df["Reject_H0"])
    logger.info(f"Found {sig_count} significant features after FDR correction")
    
    # Run Dunn's post-hoc tests for significant features
    posthoc_results = {}
    
    if sp is not None and sig_count > 0:
        logger.info("Running Dunn's post-hoc tests for significant features")
        sig_features = kw_df[kw_df["Reject_H0"]][feature_col].tolist()
        
        for i, feat in enumerate(sig_features):
            if i % 20 == 0:
                logger.debug(f"Running Dunn's test for feature {i+1}/{len(sig_features)}")
                
            sub = df_long[df_long[feature_col] == feat]
            
            try:
                # Run Dunn's test
                posthoc_df = sp.posthoc_dunn(sub, val_col=abundance_col, group_col=group_col, p_adjust="holm")
                posthoc_results[feat] = posthoc_df
            except Exception as e:
                logger.warning(f"Error in Dunn's post-hoc test for feature {feat}: {str(e)}")
    
    return kw_df, posthoc_results