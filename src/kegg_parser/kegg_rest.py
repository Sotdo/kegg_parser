"""
KEGG REST Plumbing
==================

Shared KEGG REST plumbing used by the PATHWAY and MODULE extractors: cached
HTTP text fetching, KEGG prefix stripping, tabular text parsing, the organism
gene-name list, and small per-gene aggregation helpers.

This is a library module: it exposes helpers only, with no CLI.

Output
------
- Decoded KEGG REST response text, cached under ``<outdir>/raw/``.
- Parsed dictionaries / lists / ``pandas.DataFrame`` objects.

Author: Yusheng Yang (guidance) + Agent (implementation)
Date:   2026-09-23
Version: 1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# Standard library
import io
from collections.abc import Callable, Sequence
from pathlib import Path

# Data processing
import pandas as pd

# Third-party
from loguru import logger

# Project imports
from kegg_parser.config import (
    KEGG_REST_BASE,
    MULTI_VALUE_SEPARATOR,
    RAW_DIRNAME,
    RAW_GENE_SUBDIR,
    TEXT_ENCODING,
)
from kegg_parser.utils import atomic_write_text, http_get_bytes, strip_systematic_prefix

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
GENE_LIST_FILENAME_TEMPLATE = "{org}_gene_list.txt"


# =============================================================================
# CORE LOGIC
# =============================================================================
def fetch_kegg_text(endpoint: str) -> str:
    """Fetch a KEGG REST endpoint (path after the base URL) and decode the response text."""
    logger.debug(f"Requesting KEGG REST endpoint '{endpoint}'")
    payload = http_get_bytes(f"{KEGG_REST_BASE}/{endpoint}")
    return payload.decode(TEXT_ENCODING, errors="replace")


def fetch_gene_list(org: str) -> str:
    """Fetch ``list/<org>`` (all genes with symbol and definition) from KEGG."""
    logger.info(f"Requesting KEGG gene list for organism '{org}'")
    return fetch_kegg_text(f"list/{org}")


def gene_list_cache_path(outdir: Path, org: str) -> Path:
    """Return the shared on-disk cache path for an organism's gene list."""
    return outdir / RAW_DIRNAME / RAW_GENE_SUBDIR / GENE_LIST_FILENAME_TEMPLATE.format(org=org)


def download_cached_text(cache_path: Path, force: bool, fetcher: Callable[[], str], description: str) -> str:
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


def strip_kegg_prefix(value: str) -> str:
    """Remove a leading KEGG database prefix such as ``path:``, ``md:`` or ``spo:``."""
    return value.split(":", 1)[1] if ":" in value else value


def read_kegg_tsv(text: str, columns: Sequence[str]) -> pd.DataFrame:
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


def parse_gene_list(text: str) -> dict[str, tuple[str, str]]:
    """Parse ``list/<org>`` text into a ``Gene_ID -> (Gene_Symbol, Gene_Description)`` mapping."""
    frame = read_kegg_tsv(text, ["Gene_ID", "Definition", "Position", "Name"])
    frame = frame[frame["Gene_ID"] != ""]
    genes: dict[str, tuple[str, str]] = {}
    for gene_raw, name_field in zip(frame["Gene_ID"], frame["Name"], strict=True):
        symbol_field, _, description = name_field.partition(";")
        aliases = [alias.strip() for alias in symbol_field.split(",") if alias.strip()]
        symbol = strip_systematic_prefix(aliases[0]) if aliases else ""
        genes[strip_kegg_prefix(gene_raw.strip())] = (symbol, description.strip())
    logger.info(f"Parsed gene names for {len(genes):,} genes")
    return genes


def join_unique(series: pd.Series) -> str:
    """Join unique non-empty values of a Series with the shared separator."""
    return MULTI_VALUE_SEPARATOR.join(sorted({str(value) for value in series if value}))


def first_value(series: pd.Series) -> str:
    """Return the first non-empty value of a Series."""
    return next((str(value) for value in series if value), "")
