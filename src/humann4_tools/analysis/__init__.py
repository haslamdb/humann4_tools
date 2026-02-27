# humann4_tools/analysis/__init__.py
"""Analysis functions for processed HUMAnN3 data."""

from src.humann4_tools.analysis.differential_abundance import (
    aldex2_like,
    ancom,
    ancom_bc,
    run_differential_abundance_analysis
)

from src.humann4_tools.analysis.statistical import (
    kruskal_wallis_dunn,
    run_statistical_tests
)

from src.humann4_tools.analysis.visualizations import (
    read_and_process_gene_families,
    read_and_process_pathways
)
