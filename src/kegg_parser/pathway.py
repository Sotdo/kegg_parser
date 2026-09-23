"""
KEGG PATHWAY Fetching and Mapping
=================================

Retrieve organism-specific pathway annotations through
``bioservices.kegg.KEGG`` and turn them into two gene-centric tables:

- ``pathway_gene_mapping``: one row per ``(Gene_ID, Pathway_ID)`` association,
  enriched with the gene symbol/description and the pathway name and class.
- ``gene_pathway_summary``: one row per ``Gene_ID`` with all pathway ids,
  names and classes aggregated into semicolon-separated strings.

Raw responses are cached under ``<outdir>/raw/pathway/`` so repeated runs avoid
re-querying the KEGG REST API.

This is a library module: it exposes fetch/parse functions only, with no CLI.

Output
------
- ``pandas.DataFrame`` from ``build_pathway_gene_mapping`` and
  ``build_gene_pathway_summary``.
- Two table files written under ``<outdir>/derived/`` by ``process_pathways``.

Author: Yusheng Yang (guidance) + Agent (implementation)
Date:   2026-09-22
Version: 1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# Standard library
import re
from collections import Counter
from pathlib import Path

# Data processing
import pandas as pd

# Third-party
from bioservices import KEGG
from loguru import logger

# Project imports
from kegg_parser.config import (
    DERIVED_DIRNAME,
    GENE_SUMMARY_COLUMNS,
    GENE_SUMMARY_STEM,
    ORGANISM_NAME_SEPARATOR,
    PATHWAY_MAPPING_COLUMNS,
    PATHWAY_MAPPING_STEM,
    RAW_DIRNAME,
    RAW_PATHWAY_SUBDIR,
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
PATHWAY_NUMBER_PATTERN = re.compile(r"(\d{5})$")
PATHWAY_CLASS_LINE_PATTERN = re.compile(r"^(\d{5})\s+(.*)$")
PATHWAY_CLASS_FILE_ID = "br:br08901"


# =============================================================================
# CORE LOGIC
# =============================================================================
@retryable
def _fetch_pathway_list(client: KEGG, org: str) -> str:
    """Fetch ``list/pathway/<org>`` from KEGG."""
    polite_delay()
    logger.info(f"Requesting KEGG pathway list for organism '{org}'")
    return client.list("pathway", org) or ""


@retryable
def _fetch_pathway_links(client: KEGG, org: str) -> str:
    """Fetch ``link/pathway/<org>`` (gene to pathway associations) from KEGG."""
    polite_delay()
    logger.info(f"Requesting KEGG pathway-gene links for organism '{org}'")
    return client.link("pathway", org) or ""


@retryable
def _fetch_pathway_class(client: KEGG) -> str:
    """Fetch the global ``br08901`` pathway hierarchy used to classify pathways."""
    polite_delay()
    logger.info(f"Requesting KEGG pathway hierarchy '{PATHWAY_CLASS_FILE_ID}'")
    return client.get(PATHWAY_CLASS_FILE_ID) or ""


def _strip_common_organism_suffix(names: dict[str, str]) -> None:
    """Strip the repeated `` - <organism>`` tail from organism-specific pathway names in place."""
    if not names:
        return
    suffixes: Counter[str] = Counter()
    for name in names.values():
        if ORGANISM_NAME_SEPARATOR in name:
            suffixes[name.rsplit(ORGANISM_NAME_SEPARATOR, 1)[1].strip()] += 1
    if not suffixes:
        return
    common_suffix, occurrences = suffixes.most_common(1)[0]
    if occurrences < max(2, len(names) // 2):
        return
    marker = f"{ORGANISM_NAME_SEPARATOR}{common_suffix}"
    for key, name in names.items():
        if name.endswith(marker):
            names[key] = name[: -len(marker)].strip()
    logger.debug(f"Stripped organism suffix '{common_suffix}' from {occurrences} pathway name(s)")


def parse_pathway_list(text: str) -> dict[str, str]:
    """Parse ``list/pathway`` text into a ``Pathway_ID -> Pathway_Name`` mapping."""
    frame = read_kegg_tsv(text, ["Pathway_ID", "Pathway_Name"])
    frame["Pathway_ID"] = frame["Pathway_ID"].str.strip().map(strip_kegg_prefix)
    frame["Pathway_Name"] = frame["Pathway_Name"].str.strip()
    frame = frame[(frame["Pathway_ID"] != "") & (frame["Pathway_Name"] != "")]
    names = dict(zip(frame["Pathway_ID"], frame["Pathway_Name"], strict=True))
    _strip_common_organism_suffix(names)
    logger.info(f"Parsed {len(names):,} pathway names")
    return names


def parse_pathway_links(text: str) -> list[tuple[str, str]]:
    """Parse ``link/pathway`` text into ``(Gene_ID, Pathway_ID)`` pairs."""
    frame = read_kegg_tsv(text, ["Gene_ID", "Pathway_ID"])
    frame["Gene_ID"] = frame["Gene_ID"].str.strip().map(strip_kegg_prefix)
    frame["Pathway_ID"] = frame["Pathway_ID"].str.strip().map(strip_kegg_prefix)
    frame = frame[(frame["Gene_ID"] != "") & (frame["Pathway_ID"] != "")]
    pairs = list(frame.itertuples(index=False, name=None))
    logger.info(f"Parsed {len(pairs):,} gene-pathway associations")
    return pairs


def parse_pathway_class(text: str) -> dict[str, str]:
    """Parse the ``br08901`` hierarchy into a ``map number -> 'Class A > Class B'`` mapping."""
    classes: dict[str, str] = {}
    current_a = ""
    current_b = ""
    for line in text.splitlines():
        if not line.strip():
            continue
        marker = line[0]
        match marker:
            case "A":
                current_a = line[1:].strip()
                current_b = ""
            case "B":
                current_b = line[1:].strip()
            case "C":
                class_match = PATHWAY_CLASS_LINE_PATTERN.match(line[1:].strip())
                if class_match:
                    parts = [part for part in (current_a, current_b) if part]
                    classes[class_match.group(1)] = " > ".join(parts)
    logger.info(f"Parsed pathway class for {len(classes):,} map numbers")
    return classes


def fetch_and_cache_pathway_data(
    org: str,
    outdir: Path,
    force: bool = False,
) -> tuple[dict[str, str], list[tuple[str, str]], dict[str, str], dict[str, tuple[str, str]]]:
    """Fetch (with caching) and parse the pathway list, links, class hierarchy and gene names."""
    client = create_kegg_client(org)
    raw_dir = ensure_dir(outdir / RAW_DIRNAME / RAW_PATHWAY_SUBDIR)

    list_text = download_cached_text(
        raw_dir / f"{org}_pathway_list.txt", force, lambda: _fetch_pathway_list(client, org), "pathway list"
    )
    links_text = download_cached_text(
        raw_dir / f"{org}_pathway_links.txt", force, lambda: _fetch_pathway_links(client, org), "pathway-gene links"
    )
    class_text = download_cached_text(
        raw_dir / "br08901_pathway_class.txt", force, lambda: _fetch_pathway_class(client), "pathway class hierarchy"
    )
    gene_text = download_cached_text(
        gene_list_cache_path(outdir, org), force, lambda: fetch_gene_list(client, org), "gene list"
    )

    names = parse_pathway_list(list_text)
    links = parse_pathway_links(links_text)
    classes = parse_pathway_class(class_text)
    genes = parse_gene_list(gene_text)
    return names, links, classes, genes


def build_pathway_gene_mapping(
    pathway_names: dict[str, str],
    links: list[tuple[str, str]],
    pathway_classes: dict[str, str],
    gene_names: dict[str, tuple[str, str]] | None = None,
) -> pd.DataFrame:
    """Build the deduplicated gene-to-pathway mapping table enriched with gene names."""
    gene_names = gene_names or {}
    records: list[dict[str, str]] = []
    for gene_id, pathway_id in links:
        number_match = PATHWAY_NUMBER_PATTERN.search(pathway_id)
        map_number = number_match.group(1) if number_match else ""
        symbol, description = gene_names.get(gene_id, ("", ""))
        records.append(
            {
                "Gene_ID": gene_id,
                "Gene_Symbol": symbol,
                "Gene_Description": description,
                "Pathway_ID": pathway_id,
                "Pathway_Name": pathway_names.get(pathway_id, ""),
                "Pathway_Class": pathway_classes.get(map_number, ""),
            }
        )
    frame = pd.DataFrame(records, columns=PATHWAY_MAPPING_COLUMNS).drop_duplicates().reset_index(drop=True)
    logger.info(f"Built mapping table with {len(frame):,} unique gene-pathway rows")
    return frame


def build_gene_pathway_summary(mapping: pd.DataFrame) -> pd.DataFrame:
    """Aggregate the mapping table into one row per gene."""
    if mapping.empty:
        logger.warning("Mapping table is empty; returning an empty gene summary")
        return pd.DataFrame(columns=GENE_SUMMARY_COLUMNS)

    grouped = (
        mapping.groupby("Gene_ID", sort=True)
        .agg(
            Gene_Symbol=("Gene_Symbol", first_value),
            Gene_Description=("Gene_Description", first_value),
            Pathway_IDs=("Pathway_ID", join_unique),
            Pathway_Names=("Pathway_Name", join_unique),
            Pathway_Classes=("Pathway_Class", join_unique),
            Pathway_Count=("Pathway_ID", "nunique"),
        )
        .reset_index()
    )
    grouped = grouped[GENE_SUMMARY_COLUMNS]
    logger.info(f"Built gene summary with {len(grouped):,} unique genes")
    return grouped


@logger.catch(reraise=True)
def process_pathways(
    org: str,
    outdir: Path,
    fmt: OutputFormat | str = OutputFormat.TSV,
    force: bool = False,
) -> tuple[Path, Path]:
    """Fetch, parse and persist the pathway mapping and gene summary tables."""
    logger.info(f"Processing KEGG PATHWAY data for organism '{org}'")
    names, links, classes, genes = fetch_and_cache_pathway_data(org, outdir, force=force)
    mapping = build_pathway_gene_mapping(names, links, classes, genes)
    summary = build_gene_pathway_summary(mapping)

    derived_dir = outdir / DERIVED_DIRNAME
    mapping_path = write_table(mapping, derived_dir / PATHWAY_MAPPING_STEM, fmt)
    summary_path = write_table(summary, derived_dir / GENE_SUMMARY_STEM, fmt)
    return mapping_path, summary_path
