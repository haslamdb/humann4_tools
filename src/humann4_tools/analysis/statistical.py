# humann4_tools/humann4_tools/analysis/statistical.py
import pandas as pd
import numpy as np
import logging
import traceback
import os

from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests
import scikit_posthocs as sp

from src.humann4_tools.utils.file_utils import sanitize_filename

def kruskal_wallis_dunn(df_long, group_col="Group", feature_col="Pathway", 
                       abundance_col="Abundance", alpha=0.05, logger=None):
    """
    1) Kruskal-Wallis across multiple groups
    2) Adjust p-values (Benjamini–Hochberg)
    3) Dunn's post-hoc for significant features
    
    Args:
        df_long: Long-format DataFrame with samples, features, and abundances
        group_col: Column name for grouping variable
        feature_col: Column name for feature (pathway or gene)
        abundance_col: Column name for abundance values
        alpha: Significance threshold
        logger: Logger instance for logging
        
    Returns:
        Tuple of (kw_results_df, dict_of_posthoc_dfs)
    """
    if logger:
        logger.info(f"Running Kruskal-Wallis and Dunn's (group={group_col}, feature={feature_col})")
    features = df_long[feature_col].unique()
    results = []
    for i, feat in enumerate(features):
        sub = df_long[df_long[feature_col] == feat]
        groups = sub[group_col].unique()
        group_data = [sub.loc[sub[group_col] == g, abundance_col] for g in groups]
        if any(len(x) < 2 for x in group_data):
            continue  # skip
        try:
            stat, pval = kruskal(*group_data)
            results.append({
                feature_col: feat,
                "KW_stat": stat,
                "KW_pvalue": pval
            })
        except Exception as e:
            if logger:
                logger.warning(f"Error Kruskal-Wallis on {feat}: {str(e)}")
    
    if not results:
        return pd.DataFrame(), {}
    kw_df = pd.DataFrame(results)
    reject, pvals_corrected, _, _ = multipletests(kw_df["KW_pvalue"], alpha=alpha, method="fdr_bh")
    kw_df["KW_padj"] = pvals_corrected
    kw_df["Reject_H0"] = reject
    
    # Dunn's post-hoc for those with Reject_H0 = True
    posthoc_results = {}
    sig_features = kw_df[kw_df["Reject_H0"]][feature_col].tolist()
    for feat in sig_features:
        sub = df_long[df_long[feature_col] == feat]
        try:
            posthoc_df = sp.posthoc_dunn(sub, val_col=abundance_col, group_col=group_col, p_adjust="holm")
            posthoc_results[feat] = posthoc_df
        except Exception as e:
            if logger:
                logger.warning(f"Error Dunn's on {feat}: {str(e)}")
    return kw_df, posthoc_results


def run_statistical_tests(pathways_merged, output_dir, logger, group_col="Group"):
    """
    Run statistical tests on pathway data and save results.
    
    Args:
        pathways_merged: Merged pathway data in long format
        output_dir: Directory to save results
        logger: Logger instance
        group_col: Column name for grouping variable
    """
    logger.info(f"Running statistical tests on pathways data (Kruskal-Wallis + Dunn) with grouping variable '{group_col}'.")
    try:
        kw_results, dunn_results = kruskal_wallis_dunn(
            pathways_merged,
            group_col=group_col,  # Use the user-specified grouping column
            feature_col="Pathway",
            abundance_col="Abundance",
            alpha=0.05,
            logger=logger
        )
        if kw_results.empty:
            logger.warning("No valid KW results. No file saved.")
            return
        kw_path = os.path.join(output_dir, "kruskal_wallis_pathways.csv")
        kw_results.to_csv(kw_path, index=False)
        logger.info(f"Saved Kruskal-Wallis results: {kw_path}")
        
        sig_count = sum(kw_results["Reject_H0"])
        logger.info(f"{sig_count} significant pathways found after FDR correction")
        
        if dunn_results:
            dunn_dir = os.path.join(output_dir, "dunn_posthoc_tests")
            os.makedirs(dunn_dir, exist_ok=True)
            for feat, pdf in dunn_results.items():
                safe_feat = sanitize_filename(feat)
                path = os.path.join(dunn_dir, f"dunn_{safe_feat}.csv")
                pdf.to_csv(path)
            logger.info(f"Saved Dunn's post-hoc results for {len(dunn_results)} features.")
    except Exception as e:
        logger.error(f"Error in statistical tests: {str(e)}")
        logger.error(traceback.format_exc())