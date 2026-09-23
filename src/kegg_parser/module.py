"""
KEGG MODULE Fetching and Mapping
================================

Retrieve organism-specific KEGG MODULE annotations through
``bioservices.kegg.KEGG`` and turn them into two gene-centric tables:

- ``module_gene_mapping``: one row per ``(Gene_ID, Module_ID)`` association,
  enriched with the gene symbol/description and the module name and class.
- ``gene_module_summary``: one row per ``Gene_ID`` with all module ids, names
  and classes aggregated into semicolon-separated strings.

KEGG does not serve an organism-specific module list, so module names come from
the global ``list/module`` endpoint and module classes from the ``br:ko00002``
hierarchy; associations come from ``link/module/<org>``.

This is a library module: it exposes fetch/parse functions only, with no CLI.

Output
------
- ``pandas.DataFrame`` from ``build_module_gene_mapping`` and
  ``build_gene_module_summary``.
- Two table files written under ``<outdir>/derived/`` by ``process_modules``.

Author: Yusheng Yang (guidance) + Agent (implementation)
Date:   2026-09-23
Version: 1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# Standard library
import re
from pathlib import Path

# Data processing
import pandas as pd

# Third-party
from bioservices import KEGG
from loguru import logger

# Project imports
from kegg_parser.config import (
    DERIVED_DIRNAME,
    GENE_MODULE_SUMMARY_COLUMNS,
    GENE_MODULE_SUMMARY_STEM,
    MODULE_MAPPING_COLUMNS,
    MODULE_MAPPING_STEM,
    RAW_DIRNAME,
    RAW_MODULE_SUBDIR,
    OutputFormat,
)
from kegg_parser.kegg_rest import (
    create_kegg_client,
    download_cached_text,
    fetch_gene_list,
    first_value,
    gene_list_cache_path,
    join_unique,
    parse_gene_list,
    read_kegg_tsv,
    strip_kegg_prefix,
)
from kegg_parser.utils import ensure_dir, polite_delay, retryable, write_table

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
MODULE_BASE_PATTERN = re.compile(r"(M\d{5})$")
MODULE_CLASS_LINE_PATTERN = re.compile(r"^(M\d{5})\s+(.*)$")
MODULE_CLASS_FILE_ID = "br:ko00002"
MODULE_LIST_FILE_NAME = "module_list.txt"


# =============================================================================
# CORE LOGIC
# =============================================================================
@retryable
def _fetch_module_list(client: KEGG) -> str:
    """Fetch the global ``list/module`` (module ids and names) from KEGG."""
    polite_delay()
    logger.info("Requesting KEGG module list")
    return client.list("module") or ""


@retryable
def _fetch_module_links(client: KEGG, org: str) -> str:
    """Fetch ``link/module/<org>`` (gene to module associations) from KEGG."""
    polite_delay()
    logger.info(f"Requesting KEGG module-gene links for organism '{org}'")
    return client.link("module", org) or ""


@retryable
def _fetch_module_class(client: KEGG) -> str:
    """Fetch the ``br:ko00002`` module hierarchy used to classify modules."""
    polite_delay()
    logger.info(f"Requesting KEGG module hierarchy '{MODULE_CLASS_FILE_ID}'")
    return client.get(MODULE_CLASS_FILE_ID) or ""


def parse_module_list(text: str) -> dict[str, str]:
    """Parse ``list/module`` text into a ``Module_Base -> Module_Name`` mapping."""
    frame = read_kegg_tsv(text, ["Module_Base", "Module_Name"])
    frame["Module_Base"] = frame["Module_Base"].str.strip()
    frame["Module_Name"] = frame["Module_Name"].str.strip()
    frame = frame[(frame["Module_Base"] != "") & (frame["Module_Name"] != "")]
    names = dict(zip(frame["Module_Base"], frame["Module_Name"], strict=True))
    logger.info(f"Parsed {len(names):,} module names")
    return names


def parse_module_links(text: str) -> list[tuple[str, str]]:
    """Parse ``link/module`` text into ``(Gene_ID, Module_ID)`` pairs."""
    frame = read_kegg_tsv(text, ["Gene_ID", "Module_ID"])
    frame["Gene_ID"] = frame["Gene_ID"].str.strip().map(strip_kegg_prefix)
    frame["Module_ID"] = frame["Module_ID"].str.strip().map(strip_kegg_prefix)
    frame = frame[(frame["Gene_ID"] != "") & (frame["Module_ID"] != "")]
    pairs = list(frame.itertuples(index=False, name=None))
    logger.info(f"Parsed {len(pairs):,} gene-module associations")
    return pairs


def parse_module_class(text: str) -> dict[str, str]:
    """Parse ``br:ko00002`` into a ``Module_Base -> 'Class A > Class B > Class C'`` mapping."""
    classes: dict[str, str] = {}
    current_a = ""
    current_b = ""
    current_c = ""
    for line in text.splitlines():
        if not line.strip():
            continue
        marker = line[0]
        match marker:
            case "A":
                current_a = line[1:].strip()
                current_b = ""
                current_c = ""
            case "B":
                current_b = line[1:].strip()
                current_c = ""
            case "C":
                current_c = line[1:].strip()
            case "D":
                class_match = MODULE_CLASS_LINE_PATTERN.match(line[1:].strip())
                if class_match:
                    parts = [part for part in (current_a, current_b, current_c) if part]
                    classes[class_match.group(1)] = " > ".join(parts)
    logger.info(f"Parsed module class for {len(classes):,} modules")
    return classes


def fetch_and_cache_module_data(
    org: str,
    outdir: Path,
    force: bool = False,
) -> tuple[dict[str, str], list[tuple[str, str]], dict[str, str], dict[str, tuple[str, str]]]:
    """Fetch (with caching) and parse the module list, links, class hierarchy and gene names."""
    client = create_kegg_client(org)
    raw_dir = ensure_dir(outdir / RAW_DIRNAME / RAW_MODULE_SUBDIR)

    list_text = download_cached_text(
        raw_dir / MODULE_LIST_FILE_NAME, force, lambda: _fetch_module_list(client), "module list"
    )
    links_text = download_cached_text(
        raw_dir / f"{org}_module_links.txt", force, lambda: _fetch_module_links(client, org), "module-gene links"
    )
    class_text = download_cached_text(
        raw_dir / "ko00002_module_class.txt", force, lambda: _fetch_module_class(client), "module class hierarchy"
    )
    gene_text = download_cached_text(
        gene_list_cache_path(outdir, org), force, lambda: fetch_gene_list(client, org), "gene list"
    )

    names = parse_module_list(list_text)
    links = parse_module_links(links_text)
    classes = parse_module_class(class_text)
    genes = parse_gene_list(gene_text)
    return names, links, classes, genes


def build_module_gene_mapping(
    module_names: dict[str, str],
    links: list[tuple[str, str]],
    module_classes: dict[str, str],
    gene_names: dict[str, tuple[str, str]] | None = None,
) -> pd.DataFrame:
    """Build the deduplicated gene-to-module mapping table enriched with gene names."""
    gene_names = gene_names or {}
    records: list[dict[str, str]] = []
    for gene_id, module_id in links:
        base_match = MODULE_BASE_PATTERN.search(module_id)
        module_base = base_match.group(1) if base_match else ""
        symbol, description = gene_names.get(gene_id, ("", ""))
        records.append(
            {
                "Gene_ID": gene_id,
                "Gene_Symbol": symbol,
                "Gene_Description": description,
                "Module_ID": module_id,
                "Module_Name": module_names.get(module_base, ""),
                "Module_Class": module_classes.get(module_base, ""),
            }
        )
    frame = pd.DataFrame(records, columns=MODULE_MAPPING_COLUMNS).drop_duplicates().reset_index(drop=True)
    logger.info(f"Built module mapping table with {len(frame):,} unique gene-module rows")
    return frame


def build_gene_module_summary(mapping: pd.DataFrame) -> pd.DataFrame:
    """Aggregate the module mapping table into one row per gene."""
    if mapping.empty:
        logger.warning("Module mapping table is empty; returning an empty gene summary")
        return pd.DataFrame(columns=GENE_MODULE_SUMMARY_COLUMNS)

    grouped = (
        mapping.groupby("Gene_ID", sort=True)
        .agg(
            Gene_Symbol=("Gene_Symbol", first_value),
            Gene_Description=("Gene_Description", first_value),
            Module_IDs=("Module_ID", join_unique),
            Module_Names=("Module_Name", join_unique),
            Module_Classes=("Module_Class", join_unique),
            Module_Count=("Module_ID", "nunique"),
        )
        .reset_index()
    )
    grouped = grouped[GENE_MODULE_SUMMARY_COLUMNS]
    logger.info(f"Built gene module summary with {len(grouped):,} unique genes")
    return grouped


@logger.catch(reraise=True)
def process_modules(
    org: str,
    outdir: Path,
    fmt: OutputFormat | str = OutputFormat.TSV,
    force: bool = False,
) -> tuple[Path, Path]:
    """Fetch, parse and persist the module mapping and gene summary tables."""
    logger.info(f"Processing KEGG MODULE data for organism '{org}'")
    names, links, classes, genes = fetch_and_cache_module_data(org, outdir, force=force)
    mapping = build_module_gene_mapping(names, links, classes, genes)
    summary = build_gene_module_summary(mapping)

    derived_dir = outdir / DERIVED_DIRNAME
    mapping_path = write_table(mapping, derived_dir / MODULE_MAPPING_STEM, fmt)
    summary_path = write_table(summary, derived_dir / GENE_MODULE_SUMMARY_STEM, fmt)
    return mapping_path, summary_path
