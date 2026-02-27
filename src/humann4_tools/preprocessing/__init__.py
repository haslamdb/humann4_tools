# humann4_tools/preprocessing/__init__.py
"""
Preprocessing module for HUMAnN3 Tools.

This module provides functions to preprocess raw sequence files using KneadData and run HUMAnN3 on the results.
"""

# Import pipeline functions if available
try:
    from src.humann4_tools.preprocessing.pipeline import (
        run_preprocessing_pipeline, 
        run_preprocessing_pipeline_parallel
    )
except ImportError:
    pass
