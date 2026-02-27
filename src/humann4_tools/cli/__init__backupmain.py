# humann4_tools/__init__.py
"""
HUMAnN3 Tools - A comprehensive toolkit for metagenomic analysis with HUMAnN3.

This package provides a modular workflow for processing and analyzing metagenomic data:
1. Quality control and host depletion with KneadData
2. Functional profiling with HUMAnN3
3. Joining, normalizing, and unstratifying HUMAnN3 output files
4. Statistical testing
5. Differential abundance analysis
6. Visualization
"""

__version__ = "0.1.0"

# Import key classes and functions for easy access
from src.humann4_tools.utils.input_handler import (
    find_sample_files,
    collect_files_from_metadata,
    read_samples_file
)

from src.humann4_tools.utils.file_utils import (
    check_file_exists,
    strip_suffixes_from_file_headers,
    sanitize_filename
)
