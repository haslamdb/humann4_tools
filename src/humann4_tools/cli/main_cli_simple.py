#!/usr/bin/env python3
"""
HUMAnN3 Tools - Simplified Main CLI
"""

import argparse
import sys
import os

def setup_subparsers(parser):
    """
    Set up the subparsers for each command.
    
    Args:
        parser: The main argument parser
        
    Returns:
        The parser with subparsers added
    """
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Add humann4 subcommand
    humann4_parser = subparsers.add_parser('humann4', help='Run HUMAnN3')
    humann4_parser.add_argument('--input-dir', help='Input directory for HUMAnN3')
    humann4_parser.add_argument('--output-dir', help='Output directory for HUMAnN3')
    humann4_parser.add_argument('--threads', type=int, default=1, help='Number of threads')
    
    # Add join subcommand
    join_parser = subparsers.add_parser('join', help='Join HUMAnN3 output files')
    join_parser.add_argument('--input-dir', help='Input directory for Join')
    join_parser.add_argument('--output-dir', help='Output directory for Join')
    join_parser.add_argument('--pathabundance', action='store_true', help='Process pathway abundance files')
    
    return parser

def main():
    """
    Main entry point for HUMAnN3 Tools CLI.
    """
    # Create the main parser
    parser = argparse.ArgumentParser(
        description="HUMAnN3 Tools - A comprehensive toolkit for metagenomic analysis"
    )
    
    # Add version argument
    parser.add_argument('--version', action='version', version='HUMAnN3 Tools v0.1.0')
    
    # Add subparsers
    parser = setup_subparsers(parser)
    
    # Parse arguments
    args = parser.parse_args()
    
    # If no command was provided, show help
    if not hasattr(args, 'command') or not args.command:
        parser.print_help()
        return 0
    
    # Print the arguments
    print("Command:", args.command)
    for arg, value in vars(args).items():
        if arg != 'command':
            print(f"  {arg}: {value}")
    
    # Run the appropriate command
    if args.command == 'humann4':
        print("Running HUMAnN3...")
    elif args.command == 'join':
        print("Joining files...")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())