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
import io
import logging
import re
from collections import Counter
from collections.abc import Callable, Sequence
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
    MULTI_VALUE_SEPARATOR,
    ORGANISM_NAME_SEPARATOR,
    PATHWAY_MAPPING_COLUMNS,
    PATHWAY_MAPPING_STEM,
    RAW_DIRNAME,
    RAW_PATHWAY_SUBDIR,
    TEXT_ENCODING,
    OutputFormat,
)
from kegg_parser.utils import atomic_write_text, ensure_dir, polite_delay, retryable, write_table

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
PATHWAY_NUMBER_PATTERN = re.compile(r"(\d{5})$")
PATHWAY_CLASS_LINE_PATTERN = re.compile(r"^(\d{5})\s+(.*)$")
PATHWAY_CLASS_FILE_ID = "br:br08901"


# =============================================================================
# CORE LOGIC
# =============================================================================
def _create_kegg_client(org: str) -> KEGG:
    """Instantiate a bioservices KEGG client with the organism cache pre-seeded."""
    logging.getLogger("bioservices").setLevel(logging.ERROR)
    client = KEGG()
    # KEGG retired the `list/organism` endpoint that bioservices queries to validate
    # organism codes. Seeding its private caches lets `list("pathway", org)` work
    # without that lookup while keeping `link()`/`get()` fully functional.
    client._organisms = [org]
    client._organisms_tnumbers = [org]
    return client


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


@retryable
def _fetch_gene_list(client: KEGG, org: str) -> str:
    """Fetch ``list/<org>`` (all genes with symbol and definition) from KEGG."""
    polite_delay()
    logger.info(f"Requesting KEGG gene list for organism '{org}'")
    return client.list(org) or ""


def _strip_kegg_prefix(value: str) -> str:
    """Remove a leading KEGG database prefix such as ``path:`` or ``spo:``."""
    return value.split(":", 1)[1] if ":" in value else value


def _read_kegg_tsv(text: str, columns: Sequence[str]) -> pd.DataFrame:
    """Parse tab-separated KEGG REST text into a string DataFrame with the given columns."""
    if not text.strip():
        return pd.DataFrame(columns=list(columns))
    frame = pd.read_csv(
        io.StringIO(text),
        sep="\t",
        header=None,
        dtype=str,
        keep_default_na=False,
        na_filter=False,
    ).fillna("")
    frame = frame.reindex(columns=range(len(columns)), fill_value="")
    frame.columns = list(columns)
    return frame


def _download_cached_text(cache_path: Path, force: bool, fetcher: Callable[[], str], description: str) -> str:
    """Return cached text or fetch and cache it when missing/forced."""
    if cache_path.exists() and not force:
        logger.info(f"Cache hit for {description}: {cache_path}")
        return cache_path.read_text(encoding=TEXT_ENCODING)

    text = fetcher()
    if not text.strip():
        raise RuntimeError(f"Empty response while fetching {description}")
    atomic_write_text(cache_path, text)
    logger.info(f"Saved raw {description} to {cache_path}")
    return text


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
    frame = _read_kegg_tsv(text, ["Pathway_ID", "Pathway_Name"])
    frame["Pathway_ID"] = frame["Pathway_ID"].str.strip().map(_strip_kegg_prefix)
    frame["Pathway_Name"] = frame["Pathway_Name"].str.strip()
    frame = frame[(frame["Pathway_ID"] != "") & (frame["Pathway_Name"] != "")]
    names = dict(zip(frame["Pathway_ID"], frame["Pathway_Name"], strict=True))
    _strip_common_organism_suffix(names)
    logger.info(f"Parsed {len(names):,} pathway names")
    return names


def parse_pathway_links(text: str) -> list[tuple[str, str]]:
    """Parse ``link/pathway`` text into ``(Gene_ID, Pathway_ID)`` pairs."""
    frame = _read_kegg_tsv(text, ["Gene_ID", "Pathway_ID"])
    frame["Gene_ID"] = frame["Gene_ID"].str.strip().map(_strip_kegg_prefix)
    frame["Pathway_ID"] = frame["Pathway_ID"].str.strip().map(_strip_kegg_prefix)
    frame = frame[(frame["Gene_ID"] != "") & (frame["Pathway_ID"] != "")]
    pairs = list(frame.itertuples(index=False, name=None))
    logger.info(f"Parsed {len(pairs):,} gene-pathway associations")
    return pairs


def parse_gene_list(text: str) -> dict[str, tuple[str, str]]:
    """Parse ``list/<org>`` text into a ``Gene_ID -> (Gene_Symbol, Gene_Description)`` mapping."""
    frame = _read_kegg_tsv(text, ["Gene_ID", "Definition", "Position", "Name"])
    frame = frame[frame["Gene_ID"] != ""]
    genes: dict[str, tuple[str, str]] = {}
    for gene_raw, name_field in zip(frame["Gene_ID"], frame["Name"], strict=True):
        symbol_field, _, description = name_field.partition(";")
        aliases = [alias.strip() for alias in symbol_field.split(",") if alias.strip()]
        genes[_strip_kegg_prefix(gene_raw.strip())] = (aliases[0] if aliases else "", description.strip())
    logger.info(f"Parsed gene names for {len(genes):,} genes")
    return genes


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
    client = _create_kegg_client(org)
    raw_dir = ensure_dir(outdir / RAW_DIRNAME / RAW_PATHWAY_SUBDIR)

    list_text = _download_cached_text(
        raw_dir / f"{org}_pathway_list.txt", force, lambda: _fetch_pathway_list(client, org), "pathway list"
    )
    links_text = _download_cached_text(
        raw_dir / f"{org}_pathway_links.txt", force, lambda: _fetch_pathway_links(client, org), "pathway-gene links"
    )
    class_text = _download_cached_text(
        raw_dir / "br08901_pathway_class.txt", force, lambda: _fetch_pathway_class(client), "pathway class hierarchy"
    )
    gene_text = _download_cached_text(
        raw_dir / f"{org}_gene_list.txt", force, lambda: _fetch_gene_list(client, org), "gene list"
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


def _join_unique(series: pd.Series) -> str:
    """Join unique non-empty values of a Series with the shared separator."""
    return MULTI_VALUE_SEPARATOR.join(sorted({str(value) for value in series if value}))


def _first_value(series: pd.Series) -> str:
    """Return the first non-empty value of a Series."""
    return next((str(value) for value in series if value), "")


def build_gene_pathway_summary(mapping: pd.DataFrame) -> pd.DataFrame:
    """Aggregate the mapping table into one row per gene."""
    if mapping.empty:
        logger.warning("Mapping table is empty; returning an empty gene summary")
        return pd.DataFrame(columns=GENE_SUMMARY_COLUMNS)

    grouped = (
        mapping.groupby("Gene_ID", sort=True)
        .agg(
            Gene_Symbol=("Gene_Symbol", _first_value),
            Gene_Description=("Gene_Description", _first_value),
            Pathway_IDs=("Pathway_ID", _join_unique),
            Pathway_Names=("Pathway_Name", _join_unique),
            Pathway_Classes=("Pathway_Class", _join_unique),
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
