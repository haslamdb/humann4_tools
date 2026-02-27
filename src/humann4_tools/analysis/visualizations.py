# humann4_tools/humann4_tools/analysis/visualization.py
import os
import traceback
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from src.humann4_tools.utils.file_utils import strip_suffix

def read_and_process_gene_families(unstrat_genefam, sample_key_df, output_dir, logger):
    """
    Read, process, and visualize gene family data.
    
    Args:
        unstrat_genefam: Path to unstratified gene family file
        sample_key_df: DataFrame with sample metadata
        output_dir: Directory to save outputs
        logger: Logger instance
        
    Returns:
        DataFrame with processed gene family data in long format
    """
    try:
        df = pd.read_csv(unstrat_genefam, sep="\t")
        logger.info(f"Loaded gene families: {df.shape}")
        # Clean columns
        cols = df.columns.tolist()
        new_cols = []
        for c in cols:
            if c.startswith("#"):
                c = c.lstrip("#").strip().replace(" ", "_")
            else:
                c = strip_suffix(c)
            new_cols.append(c)
        df.columns = new_cols
        
        # Melt
        if "Gene_Family" not in df.columns:
            # Try to rename if found a typical first column
            first_col = df.columns[0]
            logger.warning(f"Renaming first column '{first_col}' to 'Gene_Family'")
            df.rename(columns={first_col: "Gene_Family"}, inplace=True)
        
        long_df = df.melt(
            id_vars="Gene_Family",
            var_name="SampleName",
            value_name="Abundance"
        )
        merged = pd.merge(long_df, sample_key_df, on="SampleName", how="inner")
        if merged.empty:
            raise ValueError("No matching samples after merging gene families with sample key.")
        
        # Plot example bar
        grouped = merged.groupby(["Group", "Gene_Family"])["Abundance"].mean().reset_index()
        plt.figure(figsize=(8,4))
        sns.barplot(data=grouped.head(20), x="Group", y="Abundance", hue="Gene_Family")
        plt.title("Mean Abundance of First 20 Gene Families by Group")
        plt.xticks(rotation=45)
        plt.tight_layout()
        bar_path = os.path.join(output_dir, "gene_families_bar.svg")
        plt.savefig(bar_path, format="svg", dpi=300)
        plt.close()
        logger.info(f"Saved gene families bar plot: {bar_path}")
        
        # PCA
        pivoted = merged.pivot_table(index="SampleName", columns="Gene_Family", values="Abundance", aggfunc="sum")
        pivoted_log = np.log10(pivoted+1)
        sc = StandardScaler()
        scaled = sc.fit_transform(pivoted_log)
        pca = PCA(n_components=2)
        pca_scores = pca.fit_transform(scaled)
        pca_df = pd.DataFrame(pca_scores, columns=["PC1","PC2"], index=pivoted.index).reset_index()
        pca_merged = pd.merge(pca_df, sample_key_df, on="SampleName", how="left")
        logger.info(f"Gene families PCA variance ratio: {pca.explained_variance_ratio_[:2]}")
        
        plt.figure(figsize=(8, 5))  # Increased width to accommodate the legend
        sns.scatterplot(data=pca_merged, x="PC1", y="PC2", hue="Group", style="BMTStatus")
        plt.title("PCA on Gene Families")

        # Move legend to the left side outside the plot
        plt.legend(bbox_to_anchor=(-0.3, 0.5), loc='center right', borderaxespad=0)

        plt.tight_layout()  # This will adjust the plot to make room for the legend
        pca_path = os.path.join(output_dir, "gene_families_pca.svg")
        plt.savefig(pca_path, format="svg", dpi=300, bbox_inches='tight')  # Added bbox_inches='tight' to ensure the legend is included
        plt.close()
        logger.info(f"Saved gene families PCA plot: {pca_path}")
        
        return merged
    except Exception as e:
        logger.error(f"Error reading gene families: {str(e)}")
        logger.error(traceback.format_exc())
        raise

def read_and_process_pathways(unstrat_pathways, sample_key_df, output_dir, logger):
    """
    Read, process, and visualize pathway data.
    
    Args:
        unstrat_pathways: Path to unstratified pathway file
        sample_key_df: DataFrame with sample metadata
        output_dir: Directory to save outputs
        logger: Logger instance
        
    Returns:
        DataFrame with processed pathway data in long format
    """
    try:
        df = pd.read_csv(unstrat_pathways, sep="\t")
        logger.info(f"Loaded pathways: {df.shape}")
        
        # Clean columns
        cols = df.columns.tolist()
        new_cols = []
        if "Pathway" not in cols:
            logger.warning("Attempting to rename first column to 'Pathway'")
        for i, c in enumerate(cols):
            if i == 0 and "Pathway" not in cols:
                new_cols.append("Pathway")
            else:
                new_cols.append(strip_suffix(c))
        df.columns = new_cols
        
        # Melt
        pathways_long = df.melt(
            id_vars="Pathway",
            var_name="SampleName",
            value_name="Abundance"
        )
        merged = pd.merge(pathways_long, sample_key_df, on="SampleName", how="inner")
        if merged.empty:
            raise ValueError("No matching samples after merging pathways with sample key.")
        
        # Bar plot
        grouped = merged.groupby(["Group", "Pathway"])["Abundance"].mean().reset_index()
        plt.figure(figsize=(8,4))
        sns.barplot(data=grouped.head(20), x="Group", y="Abundance", hue="Pathway")
        plt.title("Mean Abundance of First 20 Pathways by Group")
        plt.xticks(rotation=45)
        plt.tight_layout()
        bar_path = os.path.join(output_dir, "pathways_bar.svg")
        plt.savefig(bar_path, format="svg", dpi=300)
        plt.close()
        logger.info(f"Saved pathways bar plot: {bar_path}")
        
        # PCA
        pivoted = merged.pivot_table(index="SampleName", columns="Pathway", values="Abundance", aggfunc="sum")
        pivoted_log = np.log10(pivoted+1)
        sc = StandardScaler()
        scaled = sc.fit_transform(pivoted_log)
        pca = PCA(n_components=2)
        pca_scores = pca.fit_transform(scaled)
        pca_df = pd.DataFrame(pca_scores, columns=["PC1","PC2"], index=pivoted.index).reset_index()
        pca_merged = pd.merge(pca_df, sample_key_df, on="SampleName", how="left")
        logger.info(f"Pathways PCA variance ratio: {pca.explained_variance_ratio_[:2]}")
        
        plt.figure(figsize=(8,5))
        sns.scatterplot(data=pca_merged, x="PC1", y="PC2", hue="Group", style="BMTStatus")
        plt.title("PCA on Pathways")
        plt.legend(bbox_to_anchor=(-0.3, 0.5), loc='center right', borderaxespad=0)
        plt.tight_layout()
        pca_path = os.path.join(output_dir, "pathways_pca.svg")
        plt.savefig(pca_path, format="svg", dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved pathways PCA plot: {pca_path}")
        
        return merged
    except Exception as e:
        logger.error(f"Error reading pathways: {str(e)}")
        logger.error(traceback.format_exc())
        raise