"""
KEGG BRITE Fetching and Flattening
==================================

Download organism-specific KEGG BRITE hierarchy trees (``<org>00001`` and
optionally specialized trees such as ``<org>01000``) as JSON from the KEGG
``download_htext`` endpoint, cache them under ``<outdir>/raw/brite/``, then
recursively flatten each multi-branch tree into a gene-level table.

Flattening rules
----------------
- ``Level_A`` / ``Level_B`` / ``Level_C`` / ``Level_D`` hold the ancestor labels
  from the root downwards; any classification deeper than four levels is joined
  into ``Level_D`` so no information is lost. Empty classification levels are
  forward-filled with the previous level so A-D are gap-free. The matching
  ``Level_A_ID`` / ``Level_B_ID`` / ``Level_C_ID`` / ``Level_D_ID`` columns hold
  the KEGG id parsed from each label (bracket ids such as ``[PATH:spo00010]``
  or ``[TC:1.A.8]``, or the leading 5-digit code as in ``09100 Metabolism``);
  the id prefix (``PATH:``, ``BR:``...) is stripped so ids are joinable with
  ``Pathway_ID`` values in the pathway tables.
- ``Level_E`` holds the KO entry text and ``Level_F`` the organism gene entry.
  When KEGG merges both into a tab-separated leaf, they are split apart.
- ``KO_ID`` / ``KO_Name`` are parsed from the KO entry; ``Gene_ID`` /
  ``Gene_Symbol`` / ``Gene_Description`` and ``EC_Number`` from the gene entry.
- Terminal nodes that are pure classification labels (no KO and no gene) are
  skipped, and shallow branches never raise index errors.

This is a library module: it exposes fetch/parse functions only, with no CLI.

Output
------
- ``pandas.DataFrame`` with columns defined by ``config.BriteColumn``.
- Raw JSON cached at ``<outdir>/raw/brite/<brite_id>.json``.

Author: Yusheng Yang (guidance) + Agent (implementation)
Date:   2026-09-22
Version: 1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# Standard library
import json
import re
from collections.abc import Sequence
from pathlib import Path

# Data processing
import pandas as pd

# Third-party
from loguru import logger

# Project imports
from kegg_parser.config import (
    BRITE_COLUMNS,
    BRITE_TABLE_STEM,
    DEFAULT_BRITE_SUFFIX,
    DERIVED_DIRNAME,
    KEGG_HTEXT_URL,
    KEGG_REST_BASE,
    LEVEL_SEPARATOR,
    RAW_BRITE_SUBDIR,
    RAW_DIRNAME,
    TEXT_ENCODING,
    OutputFormat,
)
from kegg_parser.utils import (
    atomic_write_bytes,
    ensure_dir,
    http_get_bytes,
    strip_systematic_prefix,
    write_table,
)

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
KO_ID_PATTERN = re.compile(r"^(K\d{5})\b")
SPECIES_PREFIX_PATTERN = re.compile(r"^[A-Za-z]{2,5}:")
EC_PATTERN = re.compile(r"\[EC:([0-9][0-9.\- ]*)\]")
LEVEL_BRACKET_ID_PATTERN = re.compile(r"\[([A-Z]{2,6}):([A-Za-z0-9._\/-]+)\]")
LEVEL_CODE_PATTERN = re.compile(r"^(\d{5})\s+(.+)$")


# =============================================================================
# CORE LOGIC
# =============================================================================
def list_organism_brite_ids(org: str) -> dict[str, str]:
    """Return the mapping of BRITE id to title for all trees available for an organism."""
    payload = http_get_bytes(f"{KEGG_REST_BASE}/list/brite/{org}")
    text = payload.decode(TEXT_ENCODING, errors="replace")
    result: dict[str, str] = {}
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) >= 2 and fields[0].strip():
            result[fields[0].strip()] = fields[1].strip()
    if not result:
        raise RuntimeError(f"KEGG returned no BRITE trees for organism '{org}'")
    return result


def default_brite_ids(org: str) -> list[str]:
    """Return every BRITE tree id KEGG advertises for an organism, falling back to the main tree."""
    fallback = [f"{org}{DEFAULT_BRITE_SUFFIX}"]
    try:
        discovered = list(list_organism_brite_ids(org).keys())
    except Exception as error:
        logger.warning(f"Could not discover BRITE trees for '{org}' ({error}); using {fallback}")
        return fallback
    if not discovered:
        return fallback
    logger.info(f"Discovered {len(discovered)} BRITE tree(s) for organism '{org}'")
    return discovered


def select_brite_ids(org: str, brite_ids: Sequence[str] | None = None) -> list[str]:
    """Resolve the BRITE ids to process from explicit ids, or every advertised tree by default."""
    if brite_ids:
        return [brite_id for brite_id in brite_ids if brite_id]
    return default_brite_ids(org)


def _validate_json_payload(payload: bytes, brite_id: str) -> None:
    """Raise a clear error when a downloaded BRITE payload is not valid JSON."""
    if not payload.strip():
        raise RuntimeError(f"Empty response while downloading BRITE '{brite_id}'")
    try:
        json.loads(payload.decode(TEXT_ENCODING))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RuntimeError(f"Invalid JSON received for BRITE '{brite_id}': {error}") from error


def download_brite_json(brite_id: str, outdir: Path, force: bool = False) -> Path:
    """Download (or reuse a cached copy of) the raw BRITE JSON tree."""
    target = ensure_dir(outdir / RAW_DIRNAME / RAW_BRITE_SUBDIR) / f"{brite_id}.json"
    if target.exists() and not force:
        try:
            json.loads(target.read_text(encoding=TEXT_ENCODING))
        except (UnicodeDecodeError, json.JSONDecodeError):
            logger.warning(f"Cached BRITE file is invalid and will be re-downloaded: {target}")
        else:
            logger.info(f"Cache hit for BRITE '{brite_id}': {target}")
            return target

    logger.info(f"Downloading BRITE tree '{brite_id}' from {KEGG_HTEXT_URL}")
    payload = http_get_bytes(KEGG_HTEXT_URL, params={"htext": brite_id, "format": "json"})
    _validate_json_payload(payload, brite_id)
    atomic_write_bytes(target, payload)
    logger.info(f"Saved raw BRITE JSON to {target} ({len(payload):,} bytes)")
    return target


def load_brite_tree(path: Path) -> dict:
    """Load a cached BRITE JSON tree from disk."""
    with Path(path).open("r", encoding=TEXT_ENCODING) as handle:
        return json.load(handle)


def _parse_ko(text: str) -> tuple[str, str]:
    """Split KO text into ``(KO_ID, KO_Name)``."""
    cleaned = (text or "").strip()
    if not cleaned:
        return "", ""
    match = KO_ID_PATTERN.match(cleaned)
    if not match:
        return "", cleaned
    return match.group(1), cleaned[match.end():].strip()


def _parse_gene(text: str) -> tuple[str, str, str]:
    """Split a gene entry into ``(Gene_ID, Gene_Symbol, Gene_Description)``."""
    cleaned = (text or "").strip()
    if not cleaned:
        return "", "", ""
    head, _, tail = cleaned.partition(";")
    tokens = head.split()
    if not tokens:
        return "", "", cleaned
    gene_id = SPECIES_PREFIX_PATTERN.sub("", tokens[0])
    symbol = strip_systematic_prefix(" ".join(tokens[1:]).strip())
    description = tail.strip()
    return gene_id, symbol, description


def _parse_ec(*texts: str) -> str:
    """Collect every EC number found across the provided text fragments."""
    numbers: list[str] = []
    for text in texts:
        for match in EC_PATTERN.finditer(text or ""):
            for number in match.group(1).split():
                if number not in numbers:
                    numbers.append(number)
    return ";".join(numbers)


def _classification_levels(path: Sequence[str]) -> list[str]:
    """Return classification levels A-D, joining deeper levels into D and forward-filling gaps."""
    levels = [
        path[0] if len(path) > 0 else "",
        path[1] if len(path) > 1 else "",
        path[2] if len(path) > 2 else "",
        LEVEL_SEPARATOR.join(path[3:]) if len(path) > 3 else "",
    ]
    filled: list[str] = []
    previous = ""
    for level in levels:
        previous = level or previous
        filled.append(previous)
    return filled


def _strip_leading_code(name: str) -> str:
    """Drop a redundant leading 5-digit KEGG code that is already captured by a bracket id."""
    code = LEVEL_CODE_PATTERN.match(name)
    return code.group(2) if code else name


def _parse_level_label(text: str) -> tuple[str, str]:
    """Split a classification label into ``(id, name)`` from a KEGG id bracket or leading code."""
    cleaned = (text or "").strip()
    if not cleaned:
        return "", ""
    bracket = LEVEL_BRACKET_ID_PATTERN.search(cleaned)
    if bracket:
        name = f"{cleaned[:bracket.start()]} {cleaned[bracket.end():]}".strip()
        return bracket.group(2), _strip_leading_code(name)
    code = LEVEL_CODE_PATTERN.match(cleaned)
    if code:
        return code.group(1), code.group(2)
    return "", cleaned


def _build_leaf_record(name: str, path: Sequence[str], brite_id: str) -> dict[str, str] | None:
    """Build one gene-level record from a terminal node, or ``None`` when it is a label."""
    cleaned = (name or "").strip()
    if not cleaned:
        return None

    if "\t" in cleaned:
        parts = [part.strip() for part in cleaned.split("\t") if part.strip()]
        ko_text = next((part for part in parts if KO_ID_PATTERN.match(part)), "")
        gene_text = next((part for part in parts if not KO_ID_PATTERN.match(part)), "")
    elif KO_ID_PATTERN.match(cleaned):
        ko_text, gene_text = cleaned, ""
    else:
        return None

    if not ko_text and not gene_text:
        return None

    level_a, level_b, level_c, level_d = _classification_levels(path)
    level_a_id, level_a = _parse_level_label(level_a)
    level_b_id, level_b = _parse_level_label(level_b)
    level_c_id, level_c = _parse_level_label(level_c)
    level_d_id, level_d = _parse_level_label(level_d)

    ko_id, ko_name = _parse_ko(ko_text)
    gene_id, gene_symbol, gene_description = _parse_gene(gene_text)

    return {
        "BRITE_ID": brite_id,
        "Level_A": level_a,
        "Level_A_ID": level_a_id,
        "Level_B": level_b,
        "Level_B_ID": level_b_id,
        "Level_C": level_c,
        "Level_C_ID": level_c_id,
        "Level_D": level_d,
        "Level_D_ID": level_d_id,
        "Level_E": ko_text,
        "Level_F": gene_text,
        "KO_ID": ko_id,
        "KO_Name": ko_name,
        "Gene_ID": gene_id,
        "Gene_Symbol": gene_symbol,
        "Gene_Description": gene_description,
        "EC_Number": _parse_ec(gene_text, ko_text),
    }


def _walk_node(node: dict, path: tuple[str, ...], brite_id: str, records: list[dict[str, str]]) -> None:
    """Recursively walk a BRITE node, appending gene-level records at the leaves."""
    name = (node.get("name") or "").strip()
    children = node.get("children") or []
    if children:
        for child in children:
            _walk_node(child, (*path, name), brite_id, records)
        return
    record = _build_leaf_record(name, path, brite_id)
    if record is not None:
        records.append(record)


def flatten_brite_tree(tree: dict, brite_id: str) -> pd.DataFrame:
    """Flatten a BRITE JSON tree into a gene-level DataFrame."""
    records: list[dict[str, str]] = []
    for child in tree.get("children") or []:
        _walk_node(child, (), brite_id, records)
    frame = pd.DataFrame(records, columns=BRITE_COLUMNS)
    gene_count = int((frame["Gene_ID"] != "").sum()) if not frame.empty else 0
    logger.info(
        f"Flattened BRITE '{brite_id}': {len(frame):,} rows "
        f"({gene_count:,} with a gene id, {len(frame) - gene_count:,} KO-only)"
    )
    return frame


def flatten_brite_json(path: Path, brite_id: str) -> pd.DataFrame:
    """Load a BRITE JSON file and flatten it into a gene-level DataFrame."""
    return flatten_brite_tree(load_brite_tree(path), brite_id)


def fetch_and_flatten_brite(brite_id: str, outdir: Path, force: bool = False) -> pd.DataFrame:
    """Download (with cache) and flatten a single BRITE tree."""
    json_path = download_brite_json(brite_id, outdir, force=force)
    return flatten_brite_json(json_path, brite_id)


@logger.catch(reraise=True)
def process_brite_trees(
    org: str,
    outdir: Path,
    fmt: OutputFormat | str = OutputFormat.TSV,
    brite_ids: Sequence[str] | None = None,
    force: bool = False,
) -> Path:
    """Download, flatten and persist every requested BRITE tree for an organism."""
    resolved_ids = select_brite_ids(org, brite_ids)
    logger.info(f"Processing {len(resolved_ids)} BRITE tree(s) for organism '{org}': {resolved_ids}")

    frames: list[pd.DataFrame] = []
    for brite_id in resolved_ids:
        frame = fetch_and_flatten_brite(brite_id, outdir, force=force)
        if frame.empty:
            logger.warning(f"BRITE '{brite_id}' produced no gene records and will be skipped")
            continue
        frames.append(frame)

    if frames:
        combined = pd.concat(frames, ignore_index=True)
    else:
        logger.warning(f"No BRITE records parsed for organism '{org}'")
        combined = pd.DataFrame(columns=BRITE_COLUMNS)

    stem_path = outdir / DERIVED_DIRNAME / BRITE_TABLE_STEM
    output_path = write_table(combined, stem_path, fmt)
    logger.info(
        f"BRITE flattening complete: {len(combined):,} rows from {len(frames)} tree(s); saved to {output_path}"
    )
    return output_path
