from setuptools import setup, find_packages

setup(
    name="humann4_tools",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        # Core data processing
        "pandas>=1.3.0",
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        
        # Statistical and scientific libraries
        "scikit-bio>=0.5.7",
        "scikit-learn>=1.0.0", 
        "scikit-posthocs>=0.7.0",
        "statsmodels>=0.13.0",
        
        # Visualization
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "matplotlib-venn>=0.11.7",
        
        # System monitoring and utilities
        "psutil>=5.9.0",
        "tqdm>=4.62.0",
        
        # Optional but recommended
        "humann>=4.0",
        "kneaddata"
    ],
    entry_points={
        'console_scripts': [
            'humann4-tools=humann4_tools.humann_tools.cli:main',
        ],
    },
    author="David Haslam",
    author_email="dbhaslam@gmail.com",
    description="A comprehensive Python package for processing and analyzing HUMAnN4 output data from metagenomic sequencing",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/dhaslam/humann4_tools",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
    python_requires=">=3.7",
)
