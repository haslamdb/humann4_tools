# humann4_tools/humann4/__init__.py
"""HUMAnN3 processing functions."""

# Import placeholder for pathway and gene processing
# These need to be implemented
try:
    from src.humann4_tools.humann4.pathway_processing import process_pathway_abundance
    from src.humann4_tools.humann4.gene_processing import process_gene_families
    from src.humann4_tools.humann4.join_unstratify import process_join_unstratify
except ImportError:
    pass
