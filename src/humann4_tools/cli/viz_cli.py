#!/usr/bin/env python3
"""
HUMAnN3 Tools Visualization Module

This module creates various visualizations from HUMAnN3 output files, including:
- PCA plots
- Heatmaps
- Barplots
- Abundance histograms
- Feature boxplots

Examples:
  # Basic usage (creates PCA plot by default):
  humann4-tools viz --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv \
                  --metadata-file metadata.csv \
                  --output-dir ./visualizations
  
  # Generate multiple plot types:
  humann4-tools viz --abundance-file joined_output/gene_families_cpm_unstratified.tsv \
                  --metadata-file metadata.csv \
                  --pca --heatmap --barplot \
                  --output-dir ./visualizations
                    
  # Generate a boxplot for a specific feature:
  humann4-tools viz --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv \
                  --metadata-file metadata.csv \
                  --feature "PWY-5484: glycolysis I" \
                  --output-dir ./visualizations
                    
  # Specify group column:
  humann4-tools viz --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv \
                  --metadata-file metadata.csv \
                  --group-col Treatment \
                  --output-dir ./visualizations
"""

import os
import sys
import time
import argparse
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional, Tuple, Union
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Import internal modules
try:
    from src.humann4_tools.utils.resource_utils import track_peak_memory
except ImportError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
    from src.humann4_tools.utils.resource_utils import track_peak_memory

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

def read_and_process_data(
    abundance_file: str, 
    metadata_file: str, 
    sample_id_col: Optional[str] = None, 
    group_col: str = "Group",
    feature_type: str = "pathway",
    log_transform: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, List[str], str, str]:
    """
    Read abundance and metadata files and prepare data for visualization.
    
    Args:
        abundance_file: Path to abundance file (unstratified)
        metadata_file: Path to metadata CSV file
        sample_id_col: Column in metadata for sample IDs (auto-detected if None)
        group_col: Column in metadata for grouping
        feature_type: Type of features in abundance file ("pathway" or "gene")
        log_transform: Whether to apply log10(x+1) transformation
        
    Returns:
        Tuple of (abundance_df, abundance_transformed, merged_long_df, groups, feature_col, sample_id_col)
    """
    logger.info(f"Reading abundance file: {abundance_file}")
    
    # Read abundance file
    try:
        abundance_df = pd.read_csv(abundance_file, sep='\t', index_col=0)
        logger.info(f"Loaded abundance data with {abundance_df.shape[0]} features and {abundance_df.shape[1]} samples")
    except Exception as e:
        logger.error(f"Error reading abundance file: {str(e)}")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), [], "", ""
    
    logger.info(f"Reading metadata file: {metadata_file}")
    
    # Read metadata
    try:
        metadata_df = pd.read_csv(metadata_file)
        logger.info(f"Loaded metadata with {metadata_df.shape[0]} samples and {metadata_df.shape[1]} columns")
    except Exception as e:
        logger.error(f"Error reading metadata file: {str(e)}")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), [], "", ""
    
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
    
    # Verify columns exist
    if sample_id_col not in metadata_df.columns:
        logger.error(f"Sample ID column '{sample_id_col}' not found in metadata")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), [], "", ""
    
    if group_col not in metadata_df.columns:
        logger.error(f"Group column '{group_col}' not found in metadata")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), [], "", ""
    
    # Get intersection of samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df[sample_id_col]))
    if not shared_samples:
        logger.error("No matching samples between abundance data and metadata")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), [], "", ""
    
    # Filter abundance data to shared samples
    abundance_filtered = abundance_df[shared_samples]
    
    # Apply log transformation if requested
    if log_transform:
        logger.info("Applying log10(x+1) transformation")
        abundance_transformed = np.log10(abundance_filtered + 1)
    else:
        abundance_transformed = abundance_filtered
    
    # Convert to long format for certain plots
    feature_col = "Pathway" if feature_type == "pathway" else "Gene_Family"
    long_df = abundance_filtered.reset_index().melt(
        id_vars=abundance_filtered.index.name, 
        var_name=sample_id_col, 
        value_name="Abundance"
    )
    
    # Rename feature column
    long_df = long_df.rename(columns={abundance_filtered.index.name: feature_col})
    
    # Merge with metadata
    merged_df = pd.merge(
        long_df, 
        metadata_df,
        on=sample_id_col,
        how="inner"
    )
    
    # Get unique groups
    groups = metadata_df[group_col].unique().tolist()
    
    logger.info(f"Successfully processed data with {len(shared_samples)} samples and {len(groups)} groups")
    
    return abundance_filtered, abundance_transformed, merged_df, groups, feature_col, sample_id_col

def generate_pca_plot(
    abundance_transformed: pd.DataFrame,
    metadata_df: pd.DataFrame,
    sample_id_col: str,
    group_col: str,
    shape_col: Optional[str] = None,
    output_dir: str = "./visualizations",
    output_format: str = "svg",
    dpi: int = 300,
    feature_type: str = "pathway"
) -> str:
    """
    Generate PCA plot from abundance data.
    
    Args:
        abundance_transformed: Transformed abundance data
        metadata_df: Metadata DataFrame
        sample_id_col: Column in metadata for sample IDs
        group_col: Column in metadata for grouping/coloring
        shape_col: Optional column in metadata for point shapes
        output_dir: Directory for output files
        output_format: Output file format (svg, png, pdf)
        dpi: DPI for raster formats
        feature_type: Type of features ("pathway" or "gene")
        
    Returns:
        Path to output file
    """
    logger.info("Generating PCA plot...")
    
    # Get shared samples
    shared_samples = list(set(abundance_transformed.columns) & set(metadata_df[sample_id_col]))
    
    # Filter data to shared samples
    abundance_filtered = abundance_transformed[shared_samples]
    
    # Scale data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(abundance_filtered.T)
    
    # Run PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_data)
    
    # Create DataFrame with PCA results
    pca_df = pd.DataFrame(
        data=pca_result,
        columns=['PC1', 'PC2'],
        index=shared_samples
    )
    
    # Add metadata
    pca_df = pd.merge(
        pca_df.reset_index(),
        metadata_df,
        left_on='index',
        right_on=sample_id_col,
        how='inner'
    )
    
    # Create plot
    plt.figure(figsize=(10, 8))
    
    if shape_col and shape_col in metadata_df.columns:
        # Use both color and shape for grouping
        ax = sns.scatterplot(
            data=pca_df,
            x='PC1',
            y='PC2',
            hue=group_col,
            style=shape_col,
            s=100,
            alpha=0.8
        )
    else:
        # Use only color for grouping
        ax = sns.scatterplot(
            data=pca_df,
            x='PC1',
            y='PC2',
            hue=group_col,
            s=100,
            alpha=0.8
        )
    
    # Add variance explained
    variance_explained = pca.explained_variance_ratio_ * 100
    plt.xlabel(f'PC1 ({variance_explained[0]:.1f}%)')
    plt.ylabel(f'PC2 ({variance_explained[1]:.1f}%)')
    
    # Add title
    feature_type_plural = "Pathways" if feature_type == "pathway" else "Gene Families"
    plt.title(f'PCA of {feature_type_plural}')
    
    # Move legend outside the plot if needed
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    # Save plot
    output_file = os.path.join(output_dir, f"pca_{feature_type}.{output_format}")
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"PCA plot saved to {output_file}")
    return output_file

def generate_heatmap(
    abundance_transformed: pd.DataFrame,
    metadata_df: pd.DataFrame,
    sample_id_col: str,
    group_col: str,
    output_dir: str = "./visualizations",
    output_format: str = "svg",
    dpi: int = 300,
    feature_type: str = "pathway",
    top_n: int = 25
) -> str:
    """
    Generate heatmap of top features.
    
    Args:
        abundance_transformed: Transformed abundance data
        metadata_df: Metadata DataFrame
        sample_id_col: Column in metadata for sample IDs
        group_col: Column in metadata for grouping/coloring
        output_dir: Directory for output files
        output_format: Output file format (svg, png, pdf)
        dpi: DPI for raster formats
        feature_type: Type of features ("pathway" or "gene")
        top_n: Number of top features to include
        
    Returns:
        Path to output file
    """
    logger.info(f"Generating heatmap of top {top_n} features...")
    
    # Get shared samples
    shared_samples = list(set(abundance_transformed.columns) & set(metadata_df[sample_id_col]))
    
    # Filter data to shared samples
    abundance_filtered = abundance_transformed[shared_samples]
    
    # Calculate feature means and select top features
    feature_means = abundance_filtered.mean(axis=1)
    top_features = feature_means.sort_values(ascending=False).head(top_n).index
    
    # Filter to top features
    top_data = abundance_filtered.loc[top_features]
    
    # Get sample grouping
    sample_groups = metadata_df.set_index(sample_id_col).loc[shared_samples, group_col]
    
    # Sort samples by group
    sorted_samples = sample_groups.sort_values().index
    
    # Prepare data for heatmap
    heatmap_data = top_data[sorted_samples]
    
    # Create row colors based on groups
    group_colors = sns.color_palette("husl", n_colors=len(sample_groups.unique()))
    group_color_map = dict(zip(sample_groups.unique(), group_colors))
    col_colors = pd.Series(sample_groups).map(group_color_map)
    
    # Create heatmap plot
    plt.figure(figsize=(14, 10))
    g = sns.clustermap(
        heatmap_data,
        col_cluster=False,
        yticklabels=True,
        cmap="viridis",
        col_colors=col_colors,
        figsize=(14, 10)
    )
    
    # Create legend for groups
    for group, color in group_color_map.items():
        g.ax_heatmap.bar(0, 0, color=color, label=group, alpha=0.8)
    g.ax_heatmap.legend(title=group_col, loc="center left", bbox_to_anchor=(1, 0.5))
    
    # Set title
    feature_type_plural = "Pathways" if feature_type == "pathway" else "Gene Families"
    plt.suptitle(f'Heatmap of Top {top_n} {feature_type_plural}', y=1.02)
    
    # Save plot
    output_file = os.path.join(output_dir, f"heatmap_top{top_n}_{feature_type}.{output_format}")
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Heatmap saved to {output_file}")
    return output_file

def generate_barplot(
    abundance_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    sample_id_col: str,
    group_col: str,
    output_dir: str = "./visualizations",
    output_format: str = "svg",
    dpi: int = 300,
    feature_type: str = "pathway",
    top_n: int = 25
) -> str:
    """
    Generate barplot of top features by group.
    
    Args:
        abundance_df: Original abundance data
        metadata_df: Metadata DataFrame
        sample_id_col: Column in metadata for sample IDs
        group_col: Column in metadata for grouping
        output_dir: Directory for output files
        output_format: Output file format (svg, png, pdf)
        dpi: DPI for raster formats
        feature_type: Type of features ("pathway" or "gene")
        top_n: Number of top features to include
        
    Returns:
        Path to output file
    """
    logger.info(f"Generating barplot of top {top_n} features...")
    
    # Get shared samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df[sample_id_col]))
    
    # Filter data to shared samples
    abundance_filtered = abundance_df[shared_samples]
    
    # Get group information for each sample
    sample_groups = metadata_df.set_index(sample_id_col).loc[shared_samples, group_col]
    
    # Calculate mean abundance per group
    group_means = {}
    for group in sample_groups.unique():
        group_samples = sample_groups[sample_groups == group].index
        group_means[group] = abundance_filtered[group_samples].mean(axis=1)
    
    # Combine group means into a DataFrame
    mean_df = pd.DataFrame(group_means)
    
    # Calculate overall mean abundance for each feature
    mean_df['overall_mean'] = mean_df.mean(axis=1)
    
    # Sort by overall mean abundance and get top N features
    top_features = mean_df.sort_values('overall_mean', ascending=False).head(top_n).drop('overall_mean', axis=1)
    
    # Transpose for easier plotting
    plot_df = top_features.transpose()
    
    # Create bar plot
    plt.figure(figsize=(12, 8))
    plot_df.plot(kind='bar', ax=plt.gca())
    
    # Customize plot
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('Group')
    plt.ylabel('Mean Abundance')
    
    feature_type_plural = "Pathways" if feature_type == "pathway" else "Gene Families"
    plt.title(f'Top {top_n} {feature_type_plural} by Mean Abundance')
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    # Save plot
    output_file = os.path.join(output_dir, f"barplot_top{top_n}_{feature_type}.{output_format}")
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Barplot saved to {output_file}")
    return output_file

def generate_feature_boxplot(
    abundance_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    feature: str,
    sample_id_col: str,
    group_col: str,
    output_dir: str = "./visualizations",
    output_format: str = "svg",
    dpi: int = 300,
    log_transform: bool = True,
    feature_type: str = "pathway"
) -> str:
    """
    Generate boxplot for a specific feature.
    
    Args:
        abundance_df: Original abundance data
        metadata_df: Metadata DataFrame
        feature: Feature to plot
        sample_id_col: Column in metadata for sample IDs
        group_col: Column in metadata for grouping
        output_dir: Directory for output files
        output_format: Output file format (svg, png, pdf)
        dpi: DPI for raster formats
        log_transform: Whether to apply log10(x+1) transformation
        feature_type: Type of features ("pathway" or "gene")
        
    Returns:
        Path to output file
    """
    logger.info(f"Generating boxplot for feature: {feature}")
    
    # Check if feature exists in abundance data
    if feature not in abundance_df.index:
        logger.error(f"Feature '{feature}' not found in abundance data")
        return ""
    
    # Get shared samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df[sample_id_col]))
    
    # Filter data to shared samples
    abundance_filtered = abundance_df[shared_samples]
    
    # Get feature abundance
    feature_abundance = abundance_filtered.loc[feature]
    
    # Apply log transformation if requested
    if log_transform:
        feature_abundance = np.log10(feature_abundance + 1)
        y_label = f"log10(Abundance + 1)"
    else:
        y_label = "Abundance"
    
    # Create DataFrame for plotting
    plot_df = pd.DataFrame({
        'Abundance': feature_abundance,
        group_col: [metadata_df.set_index(sample_id_col).loc[sample, group_col] for sample in shared_samples]
    })
    
    # Create plot
    plt.figure(figsize=(10, 6))
    
    # Draw boxplot
    sns.boxplot(x=group_col, y='Abundance', data=plot_df)
    
    # Add individual points
    sns.stripplot(x=group_col, y='Abundance', data=plot_df, 
                 color='black', alpha=0.5, jitter=True)
    
    # Set labels and title
    plt.xlabel(group_col)
    plt.ylabel(y_label)
    
    feature_type_capitalized = "Pathway" if feature_type == "pathway" else "Gene"
    plt.title(f"{feature_type_capitalized} Abundance: {feature}")
    
    # Use tight layout to ensure everything fits
    plt.tight_layout()
    
    # Save plot
    feature_name_safe = feature.replace('/', '_').replace('\\', '_').replace(' ', '_')
    output_file = os.path.join(output_dir, f"boxplot_{feature_name_safe}.{output_format}")
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Boxplot saved to {output_file}")
    return output_file

def generate_abundance_histogram(
    abundance_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    sample_id_col: str,
    group_col: str,
    output_dir: str = "./visualizations",
    output_format: str = "svg",
    dpi: int = 300,
    log_transform: bool = True,
    feature_type: str = "pathway"
) -> str:
    """
    Generate histograms of abundance distributions.
    
    Args:
        abundance_df: Original abundance data
        metadata_df: Metadata DataFrame
        sample_id_col: Column in metadata for sample IDs
        group_col: Column in metadata for grouping
        output_dir: Directory for output files
        output_format: Output file format (svg, png, pdf)
        dpi: DPI for raster formats
        log_transform: Whether to apply log10(x+1) transformation
        feature_type: Type of features ("pathway" or "gene")
        
    Returns:
        Path to output file
    """
    logger.info("Generating abundance histograms...")
    
    # Get shared samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df[sample_id_col]))
    
    # Filter data to shared samples
    abundance_filtered = abundance_df[shared_samples]
    
    # Apply log transformation if requested
    if log_transform:
        abundance_transformed = np.log10(abundance_filtered + 1)
        transform_label = "log10(Abundance + 1)"
    else:
        abundance_transformed = abundance_filtered
        transform_label = "Abundance"
    
    # Get group information
    sample_groups = metadata_df.set_index(sample_id_col).loc[shared_samples, group_col]
    
    # Create plot
    plt.figure(figsize=(12, 8))
    
    # Plot histograms for each group
    for group in sample_groups.unique():
        group_samples = sample_groups[sample_groups == group].index
        group_data = abundance_transformed[group_samples].values.flatten()
        
        # Filter out zeros and NaNs
        group_data = group_data[~np.isnan(group_data)]
        group_data = group_data[group_data > 0]
        
        sns.histplot(group_data, kde=True, label=group, alpha=0.5)
    
    # Set plot attributes
    plt.xlabel(transform_label)
    plt.ylabel('Count')
    
    feature_type_plural = "Pathways" if feature_type == "pathway" else "Gene Families"
    plt.title(f'Distribution of {feature_type_plural} Abundance by Group')
    
    plt.legend(title=group_col)
    plt.tight_layout()
    
    # Save plot
    output_file = os.path.join(output_dir, f"histogram_{feature_type}.{output_format}")
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Abundance histogram saved to {output_file}")
    return output_file

def parse_args():
    """Parse command line arguments for the Visualization module."""
    parser = argparse.ArgumentParser(
        description="Create visualizations for HUMAnN3 output files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Description:
  Visualization module creates publication-quality figures from HUMAnN3 output data to explore
  patterns, relationships, and differences between samples and groups. It supports multiple
  visualization types to help interpret functional metagenomic data.

Key Features:
  • PCA plots to visualize overall sample relationships and clustering
  • Heatmaps to display abundance patterns across samples
  • Barplots to compare top features between groups
  • Feature-specific boxplots to examine distributions of individual pathways/genes
  • Abundance histograms to visualize distribution characteristics
  • Customizable plot styling and output formats

Plot Types:
  • PCA: Principal Component Analysis for sample clustering and dimensionality reduction
  • Heatmap: Hierarchical clustering heatmap of top features across samples
  • Barplot: Grouped bar charts of top features by group
  • Boxplot: Distribution boxplots for specific features
  • Histogram: Distribution of abundance values across samples and groups

Common Usage:
  # Generate all default visualizations:
  humann4-tools viz --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv --metadata-file metadata.csv

  # Generate specific visualization types:
  humann4-tools viz --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv --metadata-file metadata.csv --pca --heatmap

  # Create boxplot for a specific pathway:
  humann4-tools viz --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv --metadata-file metadata.csv --feature "PWY-7219: adenosine ribonucleotides de novo biosynthesis"

  # Customize plot format:
  humann4-tools viz --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv --metadata-file metadata.csv --format pdf --dpi 600
"""
    )
    
    # Help info option
    parser.add_argument("--help-info", action="store_true",
                      help="Show detailed information about options and outputs")
    
    # Required arguments
    parser.add_argument("--abundance-file", required=True, 
                      help="Path to the unstratified abundance file (pathway or gene family)")
    parser.add_argument("--metadata-file", required=True,
                      help="Path to sample metadata CSV file")
    
    # Visualization options
    parser.add_argument("--output-dir", default="./Visualizations",
                      help="Directory for output files")
    parser.add_argument("--feature-type", choices=["pathway", "gene"], default="pathway",
                      help="Type of features in the abundance file (pathway or gene)")
    parser.add_argument("--group-col", default="Group",
                      help="Column name in metadata for coloring points")
    parser.add_argument("--shape-col", 
                      help="Column name in metadata for point shapes")
    parser.add_argument("--sample-id-col", 
                      help="Column name in metadata for sample IDs (autodetected if not specified)")
    parser.add_argument("--top-n", type=int, default=25,
                      help="Number of top features to include in bar plots (default: 25)")
    parser.add_argument("--format", default="svg", choices=["svg", "png", "pdf"],
                      help="Output format for plots (default: svg)")
    parser.add_argument("--dpi", type=int, default=300,
                      help="DPI for raster formats like PNG (default: 300)")
    
    # Plot selection
    parser.add_argument("--pca", action="store_true", default=True,
                      help="Generate PCA plot")
    parser.add_argument("--heatmap", action="store_true",
                      help="Generate heatmap of top features")
    parser.add_argument("--barplot", action="store_true", default=True,
                      help="Generate barplot of top features by group")
    parser.add_argument("--abundance-hist", action="store_true",
                      help="Generate histograms of abundance distributions")
    
    # Boxplot specific options
    parser.add_argument("--feature", 
                      help="Generate boxplot for a specific feature (pathway or gene)")
    
    # Additional options
    parser.add_argument("--log-transform", action="store_true", default=True,
                      help="Apply log10(x+1) transformation to abundance data")
    parser.add_argument("--log-file", 
                      help="Path to log file")
    parser.add_argument("--log-level", default="INFO", 
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                      help="Logging level")
    
    return parser.parse_args()

@track_peak_memory
def main():
    """Main function to create visualizations."""
    # Parse arguments
    args = parse_args()
    
    # Setup logging
    log_level = getattr(logging, args.log_level.upper())
    setup_logger(args.log_file, log_level)
    
    logger.info("Starting HUMAnN3 Tools Visualization Module")
    start_time = time.time()
    
    # Check if files exist
    for file_path, desc in [(args.abundance_file, "Abundance file"), (args.metadata_file, "Metadata file")]:
        if not os.path.exists(file_path):
            logger.error(f"ERROR: {desc} not found: {file_path}")
            return 1
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read and process data
    abundance_df, abundance_transformed, merged_df, groups, feature_col, sample_id_col = read_and_process_data(
        abundance_file=args.abundance_file,
        metadata_file=args.metadata_file,
        sample_id_col=args.sample_id_col,
        group_col=args.group_col,
        feature_type=args.feature_type,
        log_transform=args.log_transform
    )
    
    if abundance_df.empty:
        logger.error("Failed to process input data")
        return 1
    
    # Track which visualizations were generated
    generated_visualizations = []
    
    # Generate requested plots
    if args.pca:
        output_file = generate_pca_plot(
            abundance_transformed=abundance_transformed,
            metadata_df=merged_df,
            sample_id_col=sample_id_col,
            group_col=args.group_col,
            shape_col=args.shape_col,
            output_dir=args.output_dir,
            output_format=args.format,
            dpi=args.dpi,
            feature_type=args.feature_type
        )
        if output_file:
            generated_visualizations.append(("PCA Plot", output_file))
    
    if args.barplot:
        output_file = generate_barplot(
            abundance_df=abundance_df,
            metadata_df=merged_df,
            sample_id_col=sample_id_col,
            group_col=args.group_col,
            output_dir=args.output_dir,
            output_format=args.format,
            dpi=args.dpi,
            feature_type=args.feature_type,
            top_n=args.top_n
        )
        if output_file:
            generated_visualizations.append(("Barplot", output_file))
    
    if args.heatmap:
        output_file = generate_heatmap(
            abundance_transformed=abundance_transformed,
            metadata_df=merged_df,
            sample_id_col=sample_id_col,
            group_col=args.group_col,
            output_dir=args.output_dir,
            output_format=args.format,
            dpi=args.dpi,
            feature_type=args.feature_type,
            top_n=args.top_n
        )
        if output_file:
            generated_visualizations.append(("Heatmap", output_file))
    
    if args.abundance_hist:
        output_file = generate_abundance_histogram(
            abundance_df=abundance_df,
            metadata_df=merged_df,
            sample_id_col=sample_id_col,
            group_col=args.group_col,
            output_dir=args.output_dir,
            output_format=args.format,
            dpi=args.dpi,
            log_transform=args.log_transform,
            feature_type=args.feature_type
        )
        if output_file:
            generated_visualizations.append(("Abundance Histogram", output_file))
    
    # Handle feature boxplot if specified
    if args.feature:
        output_file = generate_feature_boxplot(
            abundance_df=abundance_df,
            metadata_df=merged_df,
            feature=args.feature,
            sample_id_col=sample_id_col,
            group_col=args.group_col,
            output_dir=args.output_dir,
            output_format=args.format,
            dpi=args.dpi,
            log_transform=args.log_transform,
            feature_type=args.feature_type
        )
        if output_file:
            generated_visualizations.append(("Feature Boxplot", output_file))
    
    # Print summary of generated visualizations
    if generated_visualizations:
        logger.info("\nGenerated Visualizations:")
        for viz_type, viz_path in generated_visualizations:
            logger.info(f"  {viz_type}: {viz_path}")
    else:
        logger.info("\nNo visualizations were generated. Use options like --pca, --heatmap, --feature, etc.")
    
    # Print elapsed time
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    
    logger.info(f"Total processing time: {int(minutes)}m {int(seconds)}s")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

def setup_logger(log_file=None, log_level=logging.INFO):
    """Set up the logger with console and optional file output."""
    # Remove any existing handlers to avoid duplication
    logger.handlers = []
    logger.setLevel(log_level)

    # Format for logs