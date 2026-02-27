#!/usr/bin/env python3
"""
Test script to verify imports for key modules in humann4_tools package.

This script attempts to import key modules in the humann4_tools package
and reports which imports succeed and which fail.
"""

import os
import sys
import importlib
import traceback

# Colors for terminal output
GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
RESET = "\033[0m"
BOLD = "\033[1m"

def attempt_import(module_path):
    """
    Attempt to import a module and return the result.
    
    Args:
        module_path: The full module path to import
        
    Returns:
        (success, error_message): Tuple with success flag and error message
    """
    try:
        importlib.import_module(module_path)
        return True, None
    except Exception as e:
        return False, f"{str(e)}\n{traceback.format_exc()}"

def main():
    """Main entry point for the script."""
    print(f"{BOLD}Testing key imports for humann4_tools package{RESET}")
    print("This will attempt to import key modules to verify the import structure is correct.\n")
    
    # Add the parent directory to path so we can import from src
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    
    # Define key modules to test
    key_modules = [
        # Core modules
        "src.humann4_tools.logger",
        
        # Utils modules
        "src.humann4_tools.utils.file_utils",
        "src.humann4_tools.utils.cmd_utils",
        "src.humann4_tools.utils.sample_utils",
        "src.humann4_tools.utils.resource_utils",
        "src.humann4_tools.utils.input_handler",
        
        # HUMAnN4 modules
        "src.humann4_tools.humann4.gene_processing",
        "src.humann4_tools.humann4.pathway_processing",
        "src.humann4_tools.humann4.join_unstratify",
        
        # Analysis modules
        "src.humann4_tools.analysis.metadata",
        "src.humann4_tools.analysis.statistical",
        "src.humann4_tools.analysis.differential_abundance",
        "src.humann4_tools.analysis.visualizations",
        
        # CLI modules
        "src.humann4_tools.cli.main_cli",
        "src.humann4_tools.cli.humann4_cli",
        "src.humann4_tools.cli.join_cli",
        "src.humann4_tools.cli.stats_cli",
        "src.humann4_tools.cli.diff_cli",
        "src.humann4_tools.cli.viz_cli",
        "src.humann4_tools.cli.kneaddata_cli"
    ]
    
    # Results tracking
    success_count = 0
    fail_count = 0
    
    # Test each module
    print(f"{BOLD}Testing individual modules...{RESET}")
    for module_path in key_modules:
        success, error = attempt_import(module_path)
        if success:
            success_count += 1
            print(f"{GREEN}✓ {module_path}{RESET}")
        else:
            fail_count += 1
            print(f"{RED}✗ {module_path}{RESET}")
            print(f"  Error: {error.split(chr(10))[0]}")  # Show only the first line of the error
    
    # Summary
    total = success_count + fail_count
    success_rate = success_count / total * 100 if total > 0 else 0
    
    print(f"\n{BOLD}Summary:{RESET}")
    print(f"Total modules tested: {total}")
    print(f"Successful imports: {success_count} ({success_rate:.1f}%)")
    print(f"Failed imports: {fail_count} ({100-success_rate:.1f}%)")
    
    # Result code
    return 0 if fail_count == 0 else 1

if __name__ == "__main__":
    sys.exit(main())