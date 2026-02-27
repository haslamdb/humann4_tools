#!/usr/bin/env python3
"""
Update Python imports in src directory.

This script finds and replaces Python import statements in the src directory:
- Changes 'from humann4_tools.' to 'from src.humann4_tools.'
- Changes 'import humann4_tools.' to 'import src.humann4_tools.'
"""

import os
import re
import sys
from typing import List, Tuple

def find_python_files(directory: str) -> List[str]:
    """Find all Python files in a directory recursively."""
    python_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.py'):
                python_files.append(os.path.join(root, file))
    return python_files

def update_imports(file_path: str) -> Tuple[int, int]:
    """
    Update imports in a Python file.
    
    Returns:
        Tuple of (from_import_count, import_count)
    """
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Replace 'from humann4_tools.' with 'from src.humann4_tools.'
    from_pattern = r'from\s+humann4_tools\.'
    updated_content, from_count = re.subn(from_pattern, 'from src.humann4_tools.', content)
    
    # Replace 'import humann4_tools.' with 'import src.humann4_tools.'
    import_pattern = r'import\s+humann4_tools\.'
    updated_content, import_count = re.subn(import_pattern, 'import src.humann4_tools.', updated_content)
    
    # Only write the file if changes were made
    if from_count > 0 or import_count > 0:
        with open(file_path, 'w') as f:
            f.write(updated_content)
    
    return from_count, import_count

def main():
    """Main function to update imports in Python files."""
    src_dir = "/home/david/Documents/Code/humann4_tools/src"
    
    if not os.path.isdir(src_dir):
        print(f"Error: Directory not found: {src_dir}")
        return 1
    
    print(f"Scanning for Python files in {src_dir}...")
    python_files = find_python_files(src_dir)
    print(f"Found {len(python_files)} Python files")
    
    total_from_count = 0
    total_import_count = 0
    modified_files = 0
    
    for file_path in python_files:
        from_count, import_count = update_imports(file_path)
        
        if from_count > 0 or import_count > 0:
            modified_files += 1
            rel_path = os.path.relpath(file_path, "/home/david/Documents/Code/humann4_tools")
            print(f"Updated {rel_path}: {from_count} 'from' imports, {import_count} 'import' statements")
        
        total_from_count += from_count
        total_import_count += import_count
    
    print("\nSummary:")
    print(f"Modified {modified_files} out of {len(python_files)} files")
    print(f"Updated {total_from_count} 'from humann4_tools.' imports")
    print(f"Updated {total_import_count} 'import humann4_tools.' statements")
    print(f"Total updates: {total_from_count + total_import_count}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())