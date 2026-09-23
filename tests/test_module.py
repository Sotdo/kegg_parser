"""
Unit tests for ``kegg_parser.module``
=====================================

Validate the parsing of the global ``list/module``, ``link/module/<org>`` and
the ``br:ko00002`` hierarchy, plus the construction of the mapping and per-gene
summary tables on small synthetic fixtures.

Author: Yusheng Yang (guidance) + Agent (implementation)
Date:   2026-09-23
Version: 1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from kegg_parser.config import GENE_MODULE_SUMMARY_COLUMNS, MODULE_MAPPING_COLUMNS
from kegg_parser.kegg_rest import parse_gene_list
from kegg_parser.module import (
    build_gene_module_summary,
    build_module_gene_mapping,
    parse_module_class,
    parse_module_links,
    parse_module_list,
)

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
MODULE_LIST_TEXT = (
    "M00001\tGlycolysis (Embden-Meyerhof pathway), glucose => pyruvate\n"
    "M00002\tGlycolysis, core module involving three-carbon compounds\n"
)
MODULE_LINKS_TEXT = "spo:2538729\tmd:spo_M00001\nspo:111\tmd:spo_M00001\nspo:2538729\tmd:spo_M00002\n"
MODULE_CLASS_TEXT = (
    "+D\tModule\n"
    "!\n"
    "APathway modules\n"
    "B  Carbohydrate metabolism\n"
    "C    Central carbohydrate metabolism\n"
    "D      M00001  Glycolysis (Embden-Meyerhof pathway) [PATH:map00010]\n"
    "D      M00002  Glycolysis, core module [PATH:map00010]\n"
)
GENE_LIST_TEXT = (
    "spo:2538729\tCDS\tI:complement(<1..5662)\ttpi1, SPOM_SPCC24B10.21; triosephosphate isomerase\n"
    "spo:111\tCDS\tI:12158..12994\tSPOM_SPAC212.08C; GPI anchored protein\n"
)


# =============================================================================
# CORE LOGIC
# =============================================================================
def test_parse_module_list() -> None:
    """The global module list maps base ids to names."""
    names = parse_module_list(MODULE_LIST_TEXT)
    assert names["M00001"] == "Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate"
    assert names["M00002"] == "Glycolysis, core module involving three-carbon compounds"


def test_parse_module_links_strips_prefixes() -> None:
    """Gene and module ids lose their KEGG database prefixes."""
    pairs = parse_module_links(MODULE_LINKS_TEXT)
    assert ("2538729", "spo_M00001") in pairs
    assert ("111", "spo_M00001") in pairs


def test_parse_module_class_builds_hierarchy() -> None:
    """The br:ko00002 hierarchy yields 'A > B > C' class strings per module."""
    classes = parse_module_class(MODULE_CLASS_TEXT)
    expected = "Pathway modules > Carbohydrate metabolism > Central carbohydrate metabolism"
    assert classes["M00001"] == expected
    assert classes["M00002"] == expected


def test_build_module_mapping_and_summary() -> None:
    """Mapping joins names/classes/gene names and the summary aggregates per gene."""
    names = parse_module_list(MODULE_LIST_TEXT)
    links = parse_module_links(MODULE_LINKS_TEXT)
    classes = parse_module_class(MODULE_CLASS_TEXT)
    genes = parse_gene_list(GENE_LIST_TEXT)

    mapping = build_module_gene_mapping(names, links, classes, genes)
    assert list(mapping.columns) == MODULE_MAPPING_COLUMNS
    assert len(mapping) == 3

    row = mapping[(mapping["Gene_ID"] == "2538729") & (mapping["Module_ID"] == "spo_M00001")].iloc[0]
    assert row["Gene_Symbol"] == "tpi1"
    assert row["Gene_Description"] == "triosephosphate isomerase"
    assert row["Module_Name"] == "Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate"
    assert row["Module_Class"] == "Pathway modules > Carbohydrate metabolism > Central carbohydrate metabolism"

    summary = build_gene_module_summary(mapping)
    assert list(summary.columns) == GENE_MODULE_SUMMARY_COLUMNS
    gene_row = summary[summary["Gene_ID"] == "2538729"].iloc[0]
    assert gene_row["Gene_Symbol"] == "tpi1"
    assert gene_row["Module_Count"] == 2
    assert gene_row["Module_IDs"] == "spo_M00001;spo_M00002"


def test_build_gene_module_summary_empty() -> None:
    """An empty mapping yields an empty summary with the right schema."""
    empty_mapping = build_module_gene_mapping({}, [], {})
    summary = build_gene_module_summary(empty_mapping)
    assert list(summary.columns) == GENE_MODULE_SUMMARY_COLUMNS
    assert summary.empty
