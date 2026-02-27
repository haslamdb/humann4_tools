# HUMAnN4 Tools - Installation Guide

This guide provides detailed instructions for installing HUMAnN4 Tools and all its dependencies.

## Prerequisites

Before installing HUMAnN4 Tools, you need to have the following dependencies installed:

1. **Python 3.7+**: The package is designed to work with Python 3.7 or newer.
2. **HUMAnN4**: The core functional profiling tool.
3. **KneadData**: For quality control and host sequence removal.

## Installation Methods

### Option 1: Using Conda (Recommended)

The recommended way to install HUMAnN4 Tools is through Conda, which will also handle all the dependencies.

```bash
# Create a new conda environment
conda create -n humann4-env python=3.9

# Activate the environment
conda activate humann4-env

# Install bioBakery tools from the biobakery channel
conda install -c biobakery humann=3.6 kneaddata=0.10.0

# Install HUMAnN4 Tools
pip install git+https://github.com/dhaslam/humann4_tools.git
```

### Option 2: Using Pip

If you already have HUMAnN4 and KneadData installed, you can install HUMAnN4 Tools directly with pip:

```bash
pip install git+https://github.com/dhaslam/humann4_tools.git
```

### Option 3: Development Installation

For development or to modify the code, clone the repository and install in development mode:

```bash
# Clone the repository
git clone https://github.com/dhaslam/humann4_tools.git

# Navigate to the directory
cd humann4_tools

# Install in development mode
pip install -e .
```

## Verifying Installation

After installation, you can verify that HUMAnN4 Tools is installed correctly:

```bash
# Check the version
humann4-tools --version

# Show the help message
humann4-tools --help
```

## Installing Additional Dependencies

Some features of HUMAnN4 Tools require additional packages:

```bash
# For differential abundance visualizations
pip install matplotlib-venn

# For statistical testing
pip install scikit-posthocs
```

## Database Setup

Before running analyses, you'll need to set up the databases for KneadData and HUMAnN4:

```bash
# Download KneadData human genome database
kneaddata_database --download human_genome bowtie2 /path/to/kneaddata_db

# Download HUMAnN4 databases
humann_databases --download chocophlan full /path/to/chocophlan
humann_databases --download uniref uniref90_diamond /path/to/uniref
```

## Troubleshooting

If you encounter issues with the installation:

1. Ensure all dependencies are correctly installed
2. Check that the correct environment is activated
3. Verify database paths are correct when running commands
4. Check log files for detailed error messages

For further assistance, please open an issue on the GitHub repository.

## Next Steps

Once installation is complete, see the [README.md](README.md) for usage instructions.
