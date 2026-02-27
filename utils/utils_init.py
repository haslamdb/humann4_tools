# src/humann4_tools/utils/__init__.py
"""
Utility modules for humann4_tools.

This package contains various utility modules:
- cmd_utils.py: Utilities for running shell commands
- file_utils.py: Utilities for file handling
- input_handler.py: Utilities for handling different input methods
- resource_utils.py: Utilities for monitoring and managing system resources
"""

# Import key functions
from src.humann4_tools.utils.file_utils import (
    check_file_exists,
    strip_suffixes_from_file_headers,
    sanitize_filename
)

from src.humann4_tools.utils.cmd_utils import run_cmd

from src.humann4_tools.utils.input_handler import (
    get_input_files,
    find_sample_files,
    collect_files_from_metadata,
    read_samples_file
)

from src.humann4_tools.utils.resource_utils import (
    track_peak_memory,
    limit_memory_usage
)
