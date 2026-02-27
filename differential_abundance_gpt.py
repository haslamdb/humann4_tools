#!/usr/bin/env python3
"""
Differential abundance analysis module providing ALDEx2‐like, ANCOM, and ANCOM-BC methods.
Refactored for maintainability, reproducibility, and CLI use.
"""

import os
import argparse
import logging
from typing import Optional, List, Dict, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests


def clr_transform(
    data_matrix: pd.DataFrame
) -> pd.DataFrame:
    """
    Centered log‐ratio (CLR) transform of a pandas DataFrame.
    
    Parameters
    ----------
    data_matrix
        rows = features, columns = samples, all values > 0
    
    Returns
    -------
    clr_data : same shape as data_matrix
    """
    log_data = np.log(data_matrix)
    gm = log_data.mean(axis=1)
    clr_data = log_data.subtract(gm, axis=0)
    return clr_data


def _prepare(
    abundance_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    group_col: str,
    denom: str,
    filter_groups: Optional[List[str]],
    pseudocount: Optional[float]
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series, np.ndarray]:
    """
    Common preprocessing: match samples, drop UNMAPPED, add pseudocount, filter groups.
    
    Returns
    -------
    abundance, metadata, groups, unique_groups
    """
    # intersect samples
    shared = list(set(abundance_df.columns) & set(metadata_df.index))
    if not shared:
        raise ValueError("No shared samples between abundance and metadata.")
    abund = abundance_df[shared].copy()
    meta = metadata_df.loc[shared].copy()
    
    # drop UNMAPPED if requested
    if denom == "unmapped_excluded" and "UNMAPPED" in abund.index:
        abund = abund.drop("UNMAPPED", axis=0)
    
    # determine pseudocount if needed
    if pseudocount is None:
        # half minimum non-zero across entire matrix
        min_nonzero = abund.values[abund.values > 0].min()
        pseudocount = float(min_nonzero) / 2.0
    abund = abund.replace(0, pseudocount)
    
    # group series
    groups = meta[group_col]
    unique_groups = groups.unique()
    
    # apply group filter
    if filter_groups is not None:
        sel = groups.isin(filter_groups)
        if not sel.any():
            raise ValueError(f"No samples found for filter_groups={filter_groups}")
        abund = abund.loc[:, sel]
        meta = meta.loc[sel]
        groups = meta[group_col]
        unique_groups = groups.unique()
        logging.getLogger(__name__).info(f"Filtered to groups: {unique_groups}")
    
    return abund, meta, groups, unique_groups


def aldex2_like(
    abundance_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    group_col: str,
    mc_samples: int = 128,
    denom: str = "all",
    filter_groups: Optional[List[str]] = None,
    pseudocount: Optional[float] = None,
    random_state: Optional[int] = None
) -> pd.DataFrame:
    """
    ALDEx2-like differential abundance via Monte Carlo Dirichlet + CLR + Welch's t-test.
    
    Returns sorted DataFrame with columns:
      - feature, effect_size, p_value, q_value,
      - mean_abundance_<group1>, mean_abundance_<group2>
    """
    logger = logging.getLogger(__name__)
    abund, meta, groups, unique_groups = _prepare(
        abundance_df, metadata_df, group_col, denom, filter_groups, pseudocount
    )
    if unique_groups.size != 2:
        raise ValueError(f"ALDEx2 requires exactly 2 groups, found {unique_groups}")
    g1, g2 = unique_groups
    samples_g1 = groups[groups == g1].index
    samples_g2 = groups[groups == g2].index
    
    rng = np.random.default_rng(random_state)
    n_feats, n_samps = abund.shape
    feature_list = abund.index.tolist()
    
    # Collect MC CLR instances: shape (mc_samples, n_feats, n_samps)
    clr_stack = np.empty((mc_samples, n_feats, n_samps))
    for j, sample in enumerate(abund.columns):
        # Dirichlet draws: shape (mc_samples, n_feats)
        alphas = abund[sample].values
        draws = rng.dirichlet(alphas, size=mc_samples) * alphas.sum()
        # CLR on each draw
        log_draws = np.log(draws + 0.5)
        gm = log_draws.mean(axis=1, keepdims=True)
        clr_stack[:, :, j] = log_draws - gm
    
    # Pre-allocate results arrays
    effects = np.zeros((mc_samples, n_feats))
    pvals = np.zeros((mc_samples, n_feats))
    
    # compute per-MC effect sizes and p-values
    for i in range(mc_samples):
        mat = pd.DataFrame(
            clr_stack[i, :, :],
            index=feature_list,
            columns=abund.columns
        )
        vals1 = mat.loc[:, samples_g1].values
        vals2 = mat.loc[:, samples_g2].values
        # mean difference
        effects[i, :] = vals1.mean(axis=1) - vals2.mean(axis=1)
        # Welch's t-test per feature
        for f_idx in range(n_feats):
            _, p = stats.ttest_ind(
                vals1[f_idx], vals2[f_idx], equal_var=False
            )
            pvals[i, f_idx] = p
    
    # median across MC
    median_effect = np.median(effects, axis=0)
    median_pval = np.median(pvals, axis=0)
    
    results = pd.DataFrame({
        'feature': feature_list,
        'effect_size': median_effect,
        'p_value': median_pval
    }).set_index('feature')
    
    # multiple testing correction
    results['q_value'] = multipletests(results['p_value'], method='fdr_bh')[1]
    
    # mean abundances
    results[f'mean_abundance_{g1}'] = abund[samples_g1].mean(axis=1)
    results[f'mean_abundance_{g2}'] = abund[samples_g2].mean(axis=1)
    
    return results.sort_values('q_value')


def ancom(
    abundance_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    group_col: str,
    alpha: float = 0.05,
    denom: str = "all",
    filter_groups: Optional[List[str]] = None,
    pseudocount: Optional[float] = None
) -> pd.DataFrame:
    """
    ANCOM: feature‐by‐feature log‐ratio tests.
    
    Returns DataFrame with columns:
      - feature, W, W_ratio, significant, mean_abundance_<group>
    """
    logger = logging.getLogger(__name__)
    abund, meta, groups, unique_groups = _prepare(
        abundance_df, metadata_df, group_col, denom, filter_groups, pseudocount
    )
    n_feats = abund.shape[0]
    feature_list = abund.index.tolist()
    
    W: Dict[str, int] = {f: 0 for f in feature_list}
    for i, fi in enumerate(feature_list):
        if i % 50 == 0:
            logger.info(f"ANCOM processing {i+1}/{n_feats}")
        for fj in feature_list:
            if fi == fj:
                continue
            lr = np.log(abund.loc[fi] / abund.loc[fj])
            groups_values = [lr[groups == g] for g in unique_groups]
            if unique_groups.size == 2:
                _, p = stats.ttest_ind(
                    groups_values[0], groups_values[1], equal_var=False
                )
            else:
                _, p = stats.f_oneway(*groups_values)
            if p < alpha:
                W[fi] += 1
    
    df = pd.DataFrame({
        'feature': feature_list,
        'W': [W[f] for f in feature_list]
    }).set_index('feature')
    df['W_ratio'] = df['W'] / (n_feats - 1)
    df['significant'] = df['W_ratio'] > 0.7
    
    for g in unique_groups:
        idx = groups[groups == g].index
        df[f'mean_abundance_{g}'] = abund[idx].mean(axis=1)
    
    return df.sort_values('W', ascending=False)


def ancom_bc(
    abundance_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    group_col: str,
    formula: Optional[str] = None,
    denom: str = "all",
    filter_groups: Optional[List[str]] = None,
    pseudocount: Optional[float] = None
) -> pd.DataFrame:
    """
    ANCOM-BC simplified: CLR + OLS per feature.
    
    Returns DataFrame with columns:
      - feature, p_value, effect_size, q_value, mean_abundance_<group>
    """
    logger = logging.getLogger(__name__)
    try:
        from statsmodels.formula.api import ols
    except ImportError:
        raise ImportError("statsmodels is required for ANCOM-BC")
    
    abund, meta, groups, unique_groups = _prepare(
        abundance_df, metadata_df, group_col, denom, filter_groups, pseudocount
    )
    
    # log + center (CLR)
    log_ab = np.log(abund)
    gm = log_ab.mean(axis=0)
    clr_ab = log_ab.subtract(gm, axis=1)
    clr_t = clr_ab.T  # samples x features
    
    if formula is None:
        formula = f"value ~ C({group_col})"
    pvals = []
    effects = []
    for i, feat in enumerate(clr_t.columns):
        if i % 100 == 0:
            logger.info(f"ANCOM-BC processing {i+1}/{clr_t.shape[1]}")
        df_tmp = pd.DataFrame({
            'value': clr_t[feat],
            **meta
        })
        model = ols(formula, data=df_tmp).fit()
        for term, p in model.pvalues.items():
            if term.startswith(f"C({group_col})"):
                pvals.append((feat, p))
                effects.append((feat, model.params[term]))
                break
    
    res = pd.DataFrame(pvals, columns=['feature', 'p_value']).set_index('feature')
    res['effect_size'] = pd.Series(dict(effects))
    res['q_value'] = multipletests(res['p_value'], method='fdr_bh')[1]
    
    for g in unique_groups:
        idx = groups[groups == g].index
        res[f'mean_abundance_{g}'] = abund[idx].mean(axis=1)
    
    return res.sort_values('q_value')


def _plot_aldex2_volcano(df_a: pd.DataFrame, output_dir: str):
    """Generate and save a volcano plot for ALDEx2 results."""
    plt.figure(figsize=(8, 6))
    plt.scatter(df_a['effect_size'], -np.log10(df_a['p_value']), alpha=0.6)
    sig = df_a['q_value'] < 0.05
    plt.scatter(df_a.loc[sig, 'effect_size'], -np.log10(df_a.loc[sig, 'p_value']),
                alpha=0.6, label='q<0.05')
    plt.xlabel('Effect size')
    plt.ylabel('-log10(p-value)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "aldex2_volcano.png"), dpi=300)
    plt.close()


def _plot_ancom_barplot(df_b: pd.DataFrame, output_dir: str):
    """Generate and save a bar plot for ANCOM results."""
    top20 = df_b.head(20)
    plt.figure(figsize=(8, 6))
    plt.barh(top20.index, top20['W_ratio'])
    plt.axvline(0.7, linestyle='--')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "ancom_top20.png"), dpi=300)
    plt.close()


def run_differential_abundance_analysis(
    abundance_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    output_dir: str,
    group_col: str = "Group",
    methods: Optional[List[str]] = None,
    denom: str = "all",
    filter_groups: Optional[List[str]] = None,
    mc_samples: int = 128,
    pseudocount: Optional[float] = None,
    random_state: Optional[int] = None,
    logger: Optional[logging.Logger] = None
) -> Dict[str, pd.DataFrame]:
    """
    Run selected differential abundance methods and save results/plots to output_dir.
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    os.makedirs(output_dir, exist_ok=True)

    if methods is None:
        methods = ["aldex2", "ancom", "ancom-bc"]
    logger.info(f"Methods to run: {methods}")

    results: Dict[str, pd.DataFrame] = {}

    if "aldex2" in methods:
        logger.info("Running ALDEx2-like analysis")
        df_a = aldex2_like(
            abundance_df, metadata_df, group_col,
            mc_samples=mc_samples,
            denom=denom,
            filter_groups=filter_groups,
            pseudocount=pseudocount,
            random_state=random_state
        )
        df_a.to_csv(os.path.join(output_dir, "aldex2_results.csv"))
        _plot_aldex2_volcano(df_a, output_dir)
        results['aldex2'] = df_a

    if "ancom" in methods:
        logger.info("Running ANCOM analysis")
        df_b = ancom(abundance_df, metadata_df, group_col, denom=denom,
                     filter_groups=filter_groups, pseudocount=pseudocount)
        df_b.to_csv(os.path.join(output_dir, "ancom_results.csv"))
        _plot_ancom_barplot(df_b, output_dir)
        results['ancom'] = df_b

    if "ancom-bc" in methods:
        logger.info("Running ANCOM-BC analysis")
        df_c = ancom_bc(abundance_df, metadata_df, group_col, denom=denom,
                        filter_groups=filter_groups, pseudocount=pseudocount)
        df_c.to_csv(os.path.join(output_dir, "ancom_bc_results.csv"))
        results['ancom-bc'] = df_c
    
    # Compare overlaps if multiple
    if len(results) > 1:
        sig_sets = {
            m: set(df[df.columns.intersection(['q_value', 'significant']).\
                     where(df[df.columns.intersection(['q_value', 'significant'])] < 0.05).any(axis=1)].index)
            for m, df in results.items()
        }
        comp_lines = ["Overlap of significant features:"]
        names = list(sig_sets.keys())
        for i in range(len(names)):
            for j in range(i+1, len(names)):
                o = len(sig_sets[names[i]].intersection(sig_sets[names[j]]))
                comp_lines.append(f"{names[i]} ∩ {names[j]} = {o}")
        with open(os.path.join(output_dir, "method_comparison.txt"), "w") as fw:
            fw.write("\n".join(comp_lines))
    
    logger.info("Analysis complete")
    return results


def main():
    p = argparse.ArgumentParser(
        description="Differential abundance: ALDEx2-like, ANCOM, ANCOM-BC"
    )
    p.add_argument("--abundance", required=True,
                   help="Feature table CSV (index=features, columns=samples)")
    p.add_argument("--metadata", required=True,
                   help="Metadata CSV (index=samples)")
    p.add_argument("--group-col", default="Group",
                   help="Column in metadata for grouping")
    p.add_argument("--methods", default="aldex2,ancom,ancom-bc",
                   help="Comma-separated methods")
    p.add_argument("--denom", choices=["all", "unmapped_excluded"],
                   default="all", help="Denominator choice")
    p.add_argument("--filter-groups",
                   help="Comma-separated groups to keep")
    p.add_argument("--mc-samples", type=int, default=128,
                   help="Monte Carlo samples for ALDEx2")
    p.add_argument("--pseudocount", type=float,
                   help="Pseudocount for zeros (None=half min nonzero)")
    p.add_argument("--random-state", type=int,
                   help="Random seed for reproducibility")
    p.add_argument("--output-dir", required=True,
                   help="Directory to save results")
    p.add_argument("--verbose", action="store_true",
                   help="Verbose logging")
    args = p.parse_args()
    
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s:%(name)s: %(message)s"
    )
    logger = logging.getLogger(__name__)
    
    abund = pd.read_csv(args.abundance, index_col=0)
    meta = pd.read_csv(args.metadata, index_col=0)
    methods = [m.strip() for m in args.methods.split(",")]
    filter_groups = (args.filter_groups.split(",")
                     if args.filter_groups else None)
    
    run_differential_abundance_analysis(
        abund, meta,
        output_dir=args.output_dir,
        group_col=args.group_col,
        methods=methods,
        denom=args.denom,
        filter_groups=filter_groups,
        mc_samples=args.mc_samples,
        pseudocount=args.pseudocount,
        random_state=args.random_state,
        logger=logger
    )


if __name__ == "__main__":
    main()
