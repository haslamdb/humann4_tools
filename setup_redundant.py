#!/usr/bin/env python3
"""
Setup script for humann4_tools package.
"""

from setuptools import setup, find_packages

# Read the long description from README.md
try:
    with open("README.md", "r") as fh:
        long_description = fh.read()
except FileNotFoundError:
    long_description = "HUMAnN4 Tools - A comprehensive toolkit for metagenomic analysis"

# Define package metadata
setup(
    name="humann4_tools",
    version="0.2.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A comprehensive toolkit for metagenomic analysis with HUMAnN4",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/humann4_tools",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.12",
    install_requires=[
        "numpy>=1.16.0",
        "pandas>=1.0.0",
        "matplotlib>=3.1.0",
        "seaborn>=0.10.0",
        "scikit-learn>=0.22.0",
        "scipy>=1.4.0",
        "statsmodels>=0.11.0",
        "scikit-posthocs>=0.6.0",
    ],
    entry_points={
        "console_scripts": [
            "humann4-tools=humann4_tools.cli.main_cli:main",
            "humann4-kneaddata=humann4_tools.cli.kneaddata_cli:main",
            "humann4-humann4=humann4_tools.cli.humann4_cli:main",
            "humann4-join=humann4_tools.cli.join_cli:main",
            "humann4-stats=humann4_tools.cli.stats_cli:main",
            "humann4-diff=humann4_tools.cli.diff_cli:main",
            "humann4-viz=humann4_tools.cli.viz_cli:main",
        ],
    },
)