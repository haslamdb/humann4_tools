#!/usr/bin/env python3
"""
Update Python imports and entry points in the project.

This script:
1. Updates imports in Python files in the src directory:
   - Changes 'from humann4_tools.' to 'from src.humann4_tools.'
   - Changes 'import humann4_tools.' to 'import src.humann4_tools.'
   - Handles various import patterns properly
2. Updates entry points in pyproject.toml
"""

import os
import re
import sys
from typing import List, Tuple, Dict

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
        lines = f.readlines()
    
    from_count = 0
    import_count = 0
    modified_lines = []
    
    # Process line by line to avoid inappropriate replacements
    for line in lines:
        # Skip lines that already have 'src.humann4_tools' to avoid double prefixing
        if 'src.humann4_tools' in line:
            modified_lines.append(line)
            continue
            
        # Match 'from humann4_tools.' pattern precisely
        from_pattern = r'^(\s*from\s+)(humann4_tools\.)(.+)$'
        from_match = re.match(from_pattern, line)
        if from_match:
            new_line = f"{from_match.group(1)}src.{from_match.group(2)}{from_match.group(3)}"
            modified_lines.append(new_line)
            from_count += 1
            continue
            
        # Match 'import humann4_tools.' pattern precisely
        import_pattern = r'^(\s*import\s+)(humann4_tools\.)(.+)$'
        import_match = re.match(import_pattern, line)
        if import_match:
            new_line = f"{import_match.group(1)}src.{import_match.group(2)}{import_match.group(3)}"
            modified_lines.append(new_line)
            import_count += 1
            continue
            
        # Match 'from humann4_tools import' pattern
        direct_import_pattern = r'^(\s*from\s+)(humann4_tools)(\s+import.+)$'
        direct_match = re.match(direct_import_pattern, line)
        if direct_match:
            new_line = f"{direct_match.group(1)}src.{direct_match.group(2)}{direct_match.group(3)}"
            modified_lines.append(new_line)
            from_count += 1
            continue
            
        # Match 'import humann4_tools' pattern
        simple_import_pattern = r'^(\s*import\s+)(humann4_tools)(\s|$)'
        simple_match = re.match(simple_import_pattern, line)
        if simple_match:
            new_line = f"{simple_match.group(1)}src.{simple_match.group(2)}{simple_match.group(3)}"
            modified_lines.append(new_line)
            import_count += 1
            continue
            
        # If no match, keep the line unchanged
        modified_lines.append(line)
    
    # Only write the file if changes were made
    if from_count > 0 or import_count > 0:
        with open(file_path, 'w') as f:
            f.writelines(modified_lines)
    
    return from_count, import_count

def update_entry_points(file_path: str) -> int:
    """
    Update entry points in pyproject.toml.
    
    Returns:
        Number of entry points updated
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    entry_point_section = False
    updates = 0
    modified_lines = []
    
    for line in lines:
        # Check if we're in the [project.scripts] section
        if '[project.scripts]' in line:
            entry_point_section = True
            modified_lines.append(line)
            continue
        
        # Check if we've left the section
        if entry_point_section and line.strip().startswith('['):
            entry_point_section = False
            modified_lines.append(line)
            continue
        
        # Update entry points 
        if entry_point_section and '=' in line and 'humann4_tools.' in line and 'src.humann4_tools.' not in line:
            # Replace 'humann4_tools.' with 'src.humann4_tools.'
            new_line = line.replace('humann4_tools.', 'src.humann4_tools.')
            modified_lines.append(new_line)
            updates += 1
            continue
        
        # Keep other lines unchanged
        modified_lines.append(line)
    
    # Only write back if changes were made
    if updates > 0:
        with open(file_path, 'w') as f:
            f.writelines(modified_lines)
    
    return updates

def update_script_files(script_dir: str) -> int:
    """Update shebang imports in script files."""
    if not os.path.exists(script_dir):
        print(f"Scripts directory not found: {script_dir}")
        return 0
        
    script_files = os.listdir(script_dir)
    updated_count = 0
    
    for file_name in script_files:
        file_path = os.path.join(script_dir, file_name)
        if not os.path.isfile(file_path) or not os.access(file_path, os.X_OK):
            continue
            
        with open(file_path, 'r') as f:
            content = f.read()
            
        # Look for imports of humann4_tools modules in executable scripts
        updated_content, count = re.subn(
            r'from humann4_tools\.', 
            'from src.humann4_tools.', 
            content
        )
        
        # Also update 'import humann4_tools' pattern
        further_updated, count2 = re.subn(
            r'import humann4_tools',
            'import src.humann4_tools',
            updated_content
        )
        
        total_count = count + count2
        
        if total_count > 0:
            with open(file_path, 'w') as f:
                f.write(further_updated)
            print(f"Updated script: {file_path} ({total_count} imports)")
            updated_count += total_count
            
    return updated_count

def test_imports() -> Dict[str, bool]:
    """
    Test if imports work correctly.
    
    Returns:
        Dict of test name to success status
    """
    results = {}
    
    # Create a temporary file for testing
    temp_file = "/tmp/import_test.py"
    
    # Test 1: Import from src.humann4_tools package
    test_code = "import src.humann4_tools\nprint('Success')"
    with open(temp_file, 'w') as f:
        f.write(test_code)
    
    try:
        result = os.system(f"python {temp_file}")
        results["Import package"] = result == 0
    except:
        results["Import package"] = False
    
    # Test 2: Import from a module
    test_code = "from src.humann4_tools.logger import setup_logger\nprint('Success')"
    with open(temp_file, 'w') as f:
        f.write(test_code)
    
    try:
        result = os.system(f"python {temp_file}")
        results["Import module"] = result == 0
    except:
        results["Import module"] = False
    
    # Cleanup
    if os.path.exists(temp_file):
        os.remove(temp_file)
    
    return results

def main():
    """Main function to update imports in Python files."""
    base_dir = "/home/david/Documents/Code/humann4_tools"
    src_dir = os.path.join(base_dir, "src")
    pyproject_toml = os.path.join(base_dir, "pyproject.toml")
    script_dir = os.path.join(base_dir, "scripts") if os.path.exists(os.path.join(base_dir, "scripts")) else None
    bin_dir = os.path.join(base_dir, "bin") if os.path.exists(os.path.join(base_dir, "bin")) else None
    
    if not os.path.isdir(src_dir):
        print(f"Error: Source directory not found: {src_dir}")
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
            rel_path = os.path.relpath(file_path, base_dir)
            print(f"Updated {rel_path}: {from_count} 'from' imports, {import_count} 'import' statements")
        
        total_from_count += from_count
        total_import_count += import_count
    
    # Update entry points in pyproject.toml
    if os.path.exists(pyproject_toml):
        print(f"\nUpdating entry points in {pyproject_toml}...")
        entry_points_updated = update_entry_points(pyproject_toml)
        print(f"Updated {entry_points_updated} entry points in pyproject.toml")
    else:
        print(f"Error: pyproject.toml not found at {pyproject_toml}")
        entry_points_updated = 0
    
    # Update executable scripts if they exist
    script_updates = 0
    if script_dir:
        print(f"\nChecking scripts in {script_dir}...")
        script_updates += update_script_files(script_dir)
        
    if bin_dir:
        print(f"\nChecking binaries in {bin_dir}...")
        script_updates += update_script_files(bin_dir)
    
    # Check for other Python files at the root level
    root_py_files = [f for f in os.listdir(base_dir) if f.endswith('.py') and os.path.isfile(os.path.join(base_dir, f))]
    if root_py_files:
        print(f"\nChecking Python files in the root directory...")
        for py_file in root_py_files:
            file_path = os.path.join(base_dir, py_file)
            from_count, import_count = update_imports(file_path)
            if from_count > 0 or import_count > 0:
                print(f"Updated {py_file}: {from_count} 'from' imports, {import_count} 'import' statements")
                total_from_count += from_count
                total_import_count += import_count
                modified_files += 1
    
    print("\nSummary:")
    print(f"Modified {modified_files} out of {len(python_files)} Python files")
    print(f"Updated {total_from_count} 'from humann4_tools.' imports")
    print(f"Updated {total_import_count} 'import humann4_tools.' statements")
    print(f"Updated {entry_points_updated} entry points in pyproject.toml")
    
    if script_updates:
        print(f"Updated {script_updates} imports in script files")
        
    print(f"Total updates: {total_from_count + total_import_count + entry_points_updated + script_updates}")
    
    print("\nTesting imports...")
    test_results = test_imports()
    for test_name, success in test_results.items():
        status = "✓ Success" if success else "✗ Failed"
        print(f"  {test_name}: {status}")
    
    print("\nNext steps:")
    print("1. Run the updated imports script with: python new_imports_fix.py")
    print("2. Reinstall the package with: pip install -e .")
    print("3. Test your entry points to ensure they work correctly")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())