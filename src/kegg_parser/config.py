"""
KEGG Parser Configuration
=========================

Central place for API endpoints, request pacing, retry policy, output layout,
and the column-name enums shared across the ``kegg_parser`` package. Keeping
these values in one module makes the fetch and parse logic easy to tune without
touching business code.

This is a library module: it exposes constants and enums only, with no CLI.

Output
------
- Constant strings / floats / tuples consumed by ``kegg_parser`` modules.
- ``StrEnum`` classes describing output table schemas and formats.

Author: Yusheng Yang (guidance) + Agent (implementation)
Date:   2026-09-22
Version: 1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from enum import StrEnum

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
# --- API endpoints -----------------------------------------------------------
KEGG_HTEXT_URL = "https://www.kegg.jp/kegg-bin/download_htext"
KEGG_REST_BASE = "https://rest.kegg.jp"

# --- Default BRITE identifiers ----------------------------------------------
DEFAULT_BRITE_SUFFIX = "00001"
MULTI_VALUE_SEPARATOR = ";"
LEVEL_SEPARATOR = " > "
ORGANISM_NAME_SEPARATOR = " - "

# --- Request pacing & retry policy ------------------------------------------
REQUEST_CONNECT_TIMEOUT = 30.0
REQUEST_READ_TIMEOUT = 600.0
REQUEST_WRITE_TIMEOUT = 30.0
REQUEST_POOL_TIMEOUT = 30.0

MIN_REQUEST_DELAY = 0.3
MAX_REQUEST_DELAY = 0.5

# Retry policy consumed by the ``stamina``-backed decorator in ``utils``.
MAX_RETRIES = 3
RETRY_TIMEOUT = 60.0
RETRY_INITIAL_DELAY = 0.5
RETRY_MAX_DELAY = 10.0
RETRY_JITTER = 1.0

USER_AGENT = "kegg-parser/1.0 (bioinformatics annotation tool)"

# --- Output layout -----------------------------------------------------------
RAW_DIRNAME = "raw"
DERIVED_DIRNAME = "derived"
RAW_BRITE_SUBDIR = "brite"
RAW_PATHWAY_SUBDIR = "pathway"
RAW_MODULE_SUBDIR = "module"
RAW_GENE_SUBDIR = "gene"

BRITE_TABLE_STEM = "brite_flat"
PATHWAY_MAPPING_STEM = "pathway_gene_mapping"
GENE_SUMMARY_STEM = "gene_pathway_summary"
MODULE_MAPPING_STEM = "module_gene_mapping"
GENE_MODULE_SUMMARY_STEM = "gene_module_summary"

# --- Missing value / encoding conventions ------------------------------------
NA_VALUE = ""
TEXT_ENCODING = "utf-8"

# Prefix KEGG/PomBase adds to systematic gene names (e.g. SPOM_SPAC186.08C).
SYSTEMATIC_NAME_PREFIX = "SPOM_"


class OutputFormat(StrEnum):
    """Supported on-disk table formats."""

    TSV = "tsv"
    PARQUET = "parquet"


class BriteColumn(StrEnum):
    """Column names of the flattened BRITE gene table."""

    BRITE_ID = "BRITE_ID"
    LEVEL_A = "Level_A"
    LEVEL_B = "Level_B"
    LEVEL_C = "Level_C"
    LEVEL_D = "Level_D"
    LEVEL_E = "Level_E"
    KO_ID = "KO_ID"
    KO_NAME = "KO_Name"
    GENE_ID = "Gene_ID"
    GENE_SYMBOL = "Gene_Symbol"
    GENE_DESCRIPTION = "Gene_Description"
    EC_NUMBER = "EC_Number"


class PathwayMappingColumn(StrEnum):
    """Column names of the gene-to-pathway mapping table."""

    GENE_ID = "Gene_ID"
    GENE_SYMBOL = "Gene_Symbol"
    GENE_DESCRIPTION = "Gene_Description"
    PATHWAY_ID = "Pathway_ID"
    PATHWAY_NAME = "Pathway_Name"
    PATHWAY_CLASS = "Pathway_Class"


class GeneSummaryColumn(StrEnum):
    """Column names of the per-gene pathway aggregation table."""

    GENE_ID = "Gene_ID"
    GENE_SYMBOL = "Gene_Symbol"
    GENE_DESCRIPTION = "Gene_Description"
    PATHWAY_IDS = "Pathway_IDs"
    PATHWAY_NAMES = "Pathway_Names"
    PATHWAY_CLASSES = "Pathway_Classes"
    PATHWAY_COUNT = "Pathway_Count"


class ModuleMappingColumn(StrEnum):
    """Column names of the gene-to-module mapping table."""

    GENE_ID = "Gene_ID"
    GENE_SYMBOL = "Gene_Symbol"
    GENE_DESCRIPTION = "Gene_Description"
    MODULE_ID = "Module_ID"
    MODULE_NAME = "Module_Name"
    MODULE_CLASS = "Module_Class"


class GeneModuleSummaryColumn(StrEnum):
    """Column names of the per-gene module aggregation table."""

    GENE_ID = "Gene_ID"
    GENE_SYMBOL = "Gene_Symbol"
    GENE_DESCRIPTION = "Gene_Description"
    MODULE_IDS = "Module_IDs"
    MODULE_NAMES = "Module_Names"
    MODULE_CLASSES = "Module_Classes"
    MODULE_COUNT = "Module_Count"


BRITE_COLUMNS: list[str] = [member.value for member in BriteColumn]
PATHWAY_MAPPING_COLUMNS: list[str] = [member.value for member in PathwayMappingColumn]
GENE_SUMMARY_COLUMNS: list[str] = [member.value for member in GeneSummaryColumn]
MODULE_MAPPING_COLUMNS: list[str] = [member.value for member in ModuleMappingColumn]
GENE_MODULE_SUMMARY_COLUMNS: list[str] = [member.value for member in GeneModuleSummaryColumn]
