# humann4-tools

A comprehensive Python package for assigning raw metagenomic sequence reads to microbial gene and pathway databases using HUMAnN4, followed by downstream processing and analysis.

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
  - [Step 1: Set up humann4 environment](#step-1-set-up-humann4-environment)
  - [Step 2: Install humann4-tools ](#step-2-install-humann4-tools)
  - [Step 3: Verify Installation](#step-3-verify-installation)
- [Workflow Overview](#workflow-overview)
- [Command Line Interface](#command-line-interface)
  - [1. KneadData](#1-kneaddata)
  - [2. HUMAnN4](#2-humann4)
  - [3. Join and Normalize](#3-join-and-normalize)
  - [4. Statistical Testing](#4-statistical-testing)
  - [5. Differential Abundance](#5-differential-abundance)
  - [6. Visualization](#6-visualization)
- [Input Methods](#input-methods)
- [Example Workflow](#example-workflow)
- [Troubleshooting](#troubleshooting)
- [Getting Help](#getting-help)
- [License](#license)

## Introduction

humann4-tools provides a complete workflow for metagenomic analysis, from quality control to visualization. The package is designed to simplify and standardize the analysis of metagenomic data using the HUMAnN4 pipeline, with modular components that can be run individually or as a complete workflow.

## Installation

### Step 1: Set up humann4 environment

**IMPORTANT**: humann4-tools, humann4-kneaddata, and humann4-humann4 tools require an installation of the humann4 workflow and is most conveniently installed using conda. After installation and configuration of humann dependencies and databases, the conda environement will include properly configured HUMAnN4, KneadData, and other dependencies. You must activate this environment before using these tools.

```bash
# Install humann4 environment (if not already installed)
conda create -n humann4 python=3.12
conda activate humann4
conda install -c biobakery -c conda-forge -c bioconda humann>=4.0 kneaddata=0.10.0

# Download and install necessary databases
# MetaPhlAn database (required for HUMAnN4)
metaphlan --install

# Install HUMAnN databases (if not already installed)
# humann_databases --download chocophlan full /path/to/databases
# humann_databases --download uniref uniref90_diamond /path/to/databases
# humann_databases --download utility_mapping full /path/to/databases
```

### Step 2: Install HUMAnN4 Tools

**ALWAYS ENSURE humann4 ENVIRONMENT IS ACTIVATED BEFORE RUNNING ANY COMMANDS**

```bash
# Make sure humann4 is activated
conda activate humann4

# Install humann4-tools from the repository
pip install git+https://github.com/haslamdb/humann4_tools.git
```

For development installation:

```bash
# Clone the repository
git clone https://github.com/haslamdb/humann4_tools.git
cd humann4_tools

# Make sure humann4 is activated
conda activate humann4

# Install in development mode
pip install -e .
```

After installation, you can use either:
- Individual commands: `humann4-kneaddata`, `humann4-humann4`, etc.
- The unified interface: `humann4-tools kneaddata`, `humann4-tools humann4`, etc.

### Step 3: Verify Installation

To verify that humann4-tools is installed correctly:

```bash
# Make sure humann4 is activated
conda activate humann4

# Check version
humann4-tools --version

# List available commands
humann4-tools --help
```

## Workflow Overview

humann4-tools provides a modular workflow for metagenomic analysis:

1. **KneadData**: Quality control and host depletion
2. **HUMAnN4**: Process cleaned sequences through HUMAnN4
3. **Join,  Normalize, Unstratify**: Combine, normalize, and split HUMAnN4 output files
4. **Statistical Testing**: Perform statistical tests across groups
5. **Differential Abundance**: Apply methods like ALDEx2, ANCOM, and ANCOM-BC
6. **Visualization**: Generate plots and figures for publication


## Command Line Interface

humann4-tools provides a consistent command-line interface for each step of the workflow. Each command can be accessed either through the main `humann4-tools` command or as individual commands.

**IMPORTANT: Always make sure the humann4 environment is activated before running these commands.**

```bash
conda activate humann4
```

### 1. KneadData

Quality control and host depletion:

```bash
# Using the main command
humann4-tools kneaddata --input-files sample_R1.fastq.gz sample_R2.fastq.gz --paired --reference-dbs human_db --output-dir kneaddata_output

# Using the standalone command
humann4-kneaddata --input-files sample_R1.fastq.gz sample_R2.fastq.gz --paired --reference-dbs human_db --output-dir kneaddata_output
```

Key options:
- `--input-files`: Input FASTQ file(s)
- `--paired`: Flag if input files are paired-end reads
- `--reference-dbs`: Path to reference database(s) for decontamination
- `--output-dir`: Directory for output files
- `--threads`: Number of threads to use

### 2. HUMAnN4

Run HUMAnN4 on cleaned sequence files:

```bash
# Using the main command
humann4-tools humann4 --input-dir kneaddata_output --output-dir humann4_output --nucleotide-db chocophlan --protein-db uniref

# Using the standalone command
humann4-humann4 --input-dir kneaddata_output --output-dir humann4_output --nucleotide-db chocophlan --protein-db uniref
```

Key options:
- `--input-dir`: Directory containing KneadData output files
- `--input-files`: Alternatively, specify input files directly. If --paired is selected, files will first be concatenated.
- `--nucleotide-db`: Path to nucleotide database
- `--protein-db`: Path to protein database
- `--output-dir`: Directory for output files
- `--threads`: Number of threads to use
- `--use-parallel`: Process multiple samples in parallel
- `--bypass-prescreen`: Skip MetaPhlAn taxonomic prescreen (useful if MetaPhlAn database isn't installed)

### 3. Join,  Normalize, Unstratify

Join, normalize, and unstratify HUMAnN4 output files:

```bash
# Using the main command
humann4-tools join --input-dir humann4_output/PathwayAbundance --pathabundance --output-dir joined_output --units cpm

# Using the standalone command
humann4-join --input-dir humann4_output/PathwayAbundance --pathabundance --output-dir joined_output --units cpm
```

Key options:
- `--input-dir`: Directory with HUMAnN4 output files
- One of: `--pathabundance`, `--pathcoverage`, or `--genefamilies` to specify file type
- `--output-dir`: Directory for output files
- `--units`: Units for normalization (cpm or relab)
- `--update-snames`: Update sample names during normalization

### 4. Statistical Testing

Run statistical tests on processed data:

```bash
# Using the main command
humann4-tools stats --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir statistical_results --group-col Treatment

# Using the standalone command
humann4-stats --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir statistical_results --group-col Treatment
```

Key options:
- `--abundance-file`: Path to unstratified abundance file
- `--metadata-file`: Path to metadata CSV file
- `--output-dir`: Directory for output files
- `--group-col`: Column name for grouping samples
- `--feature-type`: Type of features (pathway or gene)
- `--alpha`: Significance threshold (default: 0.05)

### 5. Differential Abundance

Run differential abundance analysis:

```bash
# Using the main command
humann4-tools diff --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir diff_abundance --group-col Treatment --methods aldex2,ancom,ancom-bc

# Using the standalone command
humann4-diff --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir diff_abundance --group-col Treatment --methods aldex2,ancom,ancom-bc
```

Key options:
- `--abundance-file`: Path to unstratified abundance file
- `--metadata-file`: Path to metadata CSV file
- `--output-dir`: Directory for output files
- `--group-col`: Column name for grouping samples
- `--methods`: Methods to use (aldex2, ancom, ancom-bc)
- `--filter-groups`: Filter groups for comparison (required for ALDEx2)
- `--exclude-unmapped`: Exclude unmapped features from analysis

### 6. Visualization

Create visualizations from processed data:

```bash
# Using the main command
humann4-tools viz --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir visualizations --pca --heatmap --barplot

# Using the standalone command
humann4-viz --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir visualizations --pca --heatmap --barplot
```

Key options:
- `--abundance-file`: Path to unstratified abundance file
- `--metadata-file`: Path to metadata CSV file
- `--output-dir`: Directory for output files
- `--group-col`: Column name for coloring points
- Plot selection: `--pca`, `--heatmap`, `--barplot`, `--abundance-hist`
- `--feature`: Generate boxplot for a specific feature
- `--format`: Output format (svg, png, pdf)

## Input Methods

humann4-tools supports three different input methods across all commands:

1. **Direct File Input**: Specify the input files directly with `--input-files`
2. **Sample List File**: Provide a tab-delimited file with sample IDs and file paths using `--samples-file`
3. **Metadata-Driven**: Use a metadata CSV file to locate sequence files with `--metadata-file` and `--seq-dir`

Example for metadata-driven workflow:

```bash
humann4-tools kneaddata --metadata-file metadata.csv --seq-dir /path/to/sequences --r1-suffix _R1.fastq.gz --r2-suffix _R2.fastq.gz --paired --reference-dbs human_db
```

## Example Workflow

Here's a complete example workflow:

```bash
# Activate humann4 environment (REQUIRED)
conda activate humann4

# 1. Quality control with KneadData
humann4-tools kneaddata --input-files sample1_R1.fastq.gz sample1_R2.fastq.gz --paired --reference-dbs human_db --output-dir kneaddata_output --threads 8

# 2. HUMAnN4 functional profiling
humann4-tools humann4 --input-dir kneaddata_output --output-dir humann4_output --nucleotide-db chocophlan --protein-db uniref --threads 8

# 3. Join and normalize pathway abundance files
humann4-tools join --input-dir humann4_output/PathwayAbundance --pathabundance --output-dir joined_output --units cpm

# 4. Join and normalize gene family files
humann4-tools join --input-dir humann4_output/GeneFamilies --genefamilies --output-dir joined_output --units cpm

# 5. Statistical testing
humann4-tools stats --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir statistical_results --group-col Treatment

# 6. Differential abundance analysis
humann4-tools diff --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir diff_abundance --group-col Treatment

# 7. Visualization
humann4-tools viz --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir visualizations --pca --heatmap --barplot
```

## Troubleshooting

### Common Issues

1. **Environment Issues**:
   - Make sure you've activated the humann4 environment: `conda activate humann4`
   - If you see errors about missing databases or tools, this is often because the humann4 environment is not activated

2. **Missing Databases**:
   - MetaPhlAn error: Run `metaphlan --install` to install the required database
   - HUMAnN4 database errors: Make sure to install the ChocoPhlAn and UniRef databases

3. **Missing Input Files**: 
   - Check file paths and naming patterns
   - Verify that you're using the correct directory structure

4. **Sample Key Issues**:
   - Ensure sample identifiers in the CSV match the file names
   - Check for duplicate sample IDs
   - Verify CSV encoding (use UTF-8)

5. **KneadData Errors**:
   - Ensure reference databases are properly built and indexed
   - For paired-end issues, try different `--decontaminate-pairs` options

6. **HUMAnN4 Errors**:
   - If you get a "No MetaPhlAn BowTie2 database found" error, run `metaphlan --install`
   - As a workaround, you can use `--bypass-prescreen` to skip the MetaPhlAn step
   - Ensure nucleotide and protein databases are correctly installed
   - Check that input files are in the correct format

7. **Statistical Test Errors**:
   - ALDEx2 requires exactly two groups
   - Check for missing values in abundance data

8. **Memory Issues with Large Datasets**:
   - Run steps separately to manage memory usage
   - Use `--threads` to control CPU usage

## Getting Help

If you encounter issues not covered in this documentation, please:

- Check the log file for detailed error messages
- Set `--log-level DEBUG` for more verbose output
- Verify you're using the humann4 environment
- Open an issue on the GitHub repository with a description of the problem and relevant log entries

## License

This project is licensed under the MIT License - see the LICENSE file for details.
