# src/humann4_tools/cli/__init__.py
"""
Command-line interface modules for humann4_tools.

This package contains the command-line interface modules for each step of the workflow:
- kneaddata_cli.py: Quality control and host depletion
- humann4_cli.py: Functional profiling with HUMAnN3
- join_cli.py: Joining, normalizing, and unstratifying HUMAnN3 output files
- stats_cli.py: Statistical testing
- diff_cli.py: Differential abundance analysis
- viz_cli.py: Visualization
- main_cli.py: Main CLI interface that dispatches to the other modules
"""

# Import main CLI entry point
from src.humann4_tools.cli.main_cli import main
