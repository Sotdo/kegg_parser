"""
Unit tests for ``kegg_parser.pathway``
=====================================

Validate the parsing of KEGG ``list/pathway``, ``link/pathway`` and the
``br08901`` hierarchy, plus the construction of the mapping and per-gene
summary tables on small synthetic fixtures.

Author: Yusheng Yang (guidance) + Agent (implementation)
Date:   2026-09-22
Version: 1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from kegg_parser.config import GENE_SUMMARY_COLUMNS, PATHWAY_MAPPING_COLUMNS
from kegg_parser.kegg_rest import parse_gene_list
from kegg_parser.pathway import (
    build_gene_pathway_summary,
    build_pathway_gene_mapping,
    parse_pathway_class,
    parse_pathway_links,
    parse_pathway_list,
)

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
ORG_SUFFIX = " - Schizosaccharomyces pombe (fission yeast)"
LIST_TEXT = (
    f"spo01100\tMetabolic pathways{ORG_SUFFIX}\n"
    f"spo00533\tGlycosaminoglycan biosynthesis - keratan sulfate{ORG_SUFFIX}\n"
    f"spo00010\tGlycolysis / Gluconeogenesis{ORG_SUFFIX}\n"
)
LINKS_TEXT = "spo:2538729\tpath:spo00010\nspo:2538729\tpath:spo01100\nspo:111\tpath:spo00533\n"
CLASS_TEXT = (
    "AMetabolism\n"
    "B  Global and overview maps\n"
    "C    01100  Metabolic pathways\n"
    "C    00533  Glycosaminoglycan biosynthesis - keratan sulfate\n"
    "B  Carbohydrate metabolism\n"
    "C    00010  Glycolysis / Gluconeogenesis\n"
)
GENE_LIST_TEXT = (
    "spo:2538729\tCDS\tI:complement(<1..5662)\ttpi1, SPOM_SPCC24B10.21; triosephosphate isomerase\n"
    "spo:111\tCDS\tI:12158..12994\tSPOM_SPAC212.08C; GPI anchored protein\n"
)


# =============================================================================
# CORE LOGIC
# =============================================================================
def test_parse_pathway_list_strips_organism_suffix() -> None:
    """Only the repeated organism suffix is removed, not hyphens in real names."""
    names = parse_pathway_list(LIST_TEXT)
    assert names["spo00010"] == "Glycolysis / Gluconeogenesis"
    assert names["spo01100"] == "Metabolic pathways"
    assert names["spo00533"] == "Glycosaminoglycan biosynthesis - keratan sulfate"


def test_parse_pathway_links_strips_prefixes() -> None:
    """Gene and pathway ids lose their KEGG database prefixes."""
    pairs = parse_pathway_links(LINKS_TEXT)
    assert ("2538729", "spo00010") in pairs
    assert ("111", "spo00533") in pairs


def test_parse_pathway_class_builds_hierarchy() -> None:
    """The br08901 hierarchy yields 'A > B' class strings per map number."""
    classes = parse_pathway_class(CLASS_TEXT)
    assert classes["01100"] == "Metabolism > Global and overview maps"
    assert classes["00010"] == "Metabolism > Carbohydrate metabolism"


def test_parse_gene_list_extracts_symbol_and_description() -> None:
    """The gene list yields the primary alias and definition per gene id."""
    genes = parse_gene_list(GENE_LIST_TEXT)
    assert genes["2538729"] == ("tpi1", "triosephosphate isomerase")
    assert genes["111"] == ("SPAC212.08C", "GPI anchored protein")


def test_build_mapping_and_summary() -> None:
    """Mapping joins names/classes/gene names and the summary aggregates per gene."""
    names = parse_pathway_list(LIST_TEXT)
    links = parse_pathway_links(LINKS_TEXT)
    classes = parse_pathway_class(CLASS_TEXT)
    genes = parse_gene_list(GENE_LIST_TEXT)

    mapping = build_pathway_gene_mapping(names, links, classes, genes)
    assert list(mapping.columns) == PATHWAY_MAPPING_COLUMNS
    assert len(mapping) == 3

    row = mapping[(mapping["Gene_ID"] == "2538729") & (mapping["Pathway_ID"] == "spo00010")].iloc[0]
    assert row["Gene_Symbol"] == "tpi1"
    assert row["Gene_Description"] == "triosephosphate isomerase"
    assert row["Pathway_Name"] == "Glycolysis / Gluconeogenesis"
    assert row["Pathway_Class"] == "Metabolism > Carbohydrate metabolism"

    summary = build_gene_pathway_summary(mapping)
    assert list(summary.columns) == GENE_SUMMARY_COLUMNS
    gene_row = summary[summary["Gene_ID"] == "2538729"].iloc[0]
    assert gene_row["Gene_Symbol"] == "tpi1"
    assert gene_row["Gene_Description"] == "triosephosphate isomerase"
    assert gene_row["Pathway_Count"] == 2
    assert gene_row["Pathway_IDs"] == "spo00010;spo01100"


def test_build_gene_pathway_summary_empty() -> None:
    """An empty mapping yields an empty summary with the right schema."""
    empty_mapping = build_pathway_gene_mapping({}, [], {})
    summary = build_gene_pathway_summary(empty_mapping)
    assert list(summary.columns) == GENE_SUMMARY_COLUMNS
    assert summary.empty
