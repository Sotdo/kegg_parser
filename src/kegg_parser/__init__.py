"""
kegg_parser
===========

Tools to fetch, parse and flatten KEGG BRITE hierarchies and KEGG PATHWAY
annotations for a given organism into gene-level TSV or Parquet tables.

The package exposes two functional areas:

- ``kegg_parser.brite``: download and flatten BRITE JSON trees.
- ``kegg_parser.pathway``: fetch pathway/gene links and build mapping tables.

Author: Yusheng Yang (guidance) + Agent (implementation)
Date:   2026-09-22
Version: 1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from kegg_parser.brite import (
    default_brite_ids,
    download_brite_json,
    fetch_and_flatten_brite,
    flatten_brite_json,
    flatten_brite_tree,
    list_organism_brite_ids,
    process_brite_trees,
    select_brite_ids,
)
from kegg_parser.config import (
    BriteColumn,
    GeneSummaryColumn,
    OutputFormat,
    PathwayMappingColumn,
)
from kegg_parser.pathway import (
    build_gene_pathway_summary,
    build_pathway_gene_mapping,
    process_pathways,
)
from kegg_parser.utils import setup_logger

__version__ = "1.0.0"

__all__ = [
    "BriteColumn",
    "GeneSummaryColumn",
    "OutputFormat",
    "PathwayMappingColumn",
    "build_gene_pathway_summary",
    "build_pathway_gene_mapping",
    "default_brite_ids",
    "download_brite_json",
    "fetch_and_flatten_brite",
    "flatten_brite_json",
    "flatten_brite_tree",
    "list_organism_brite_ids",
    "process_brite_trees",
    "process_pathways",
    "select_brite_ids",
    "setup_logger",
]
