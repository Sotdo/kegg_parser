#!/usr/bin/env python3

"""
Fetch and Flatten KEGG BRITE Trees
==================================

Download organism-specific KEGG BRITE hierarchy trees as JSON from the KEGG
``download_htext`` endpoint, cache them under ``<outdir>/raw/brite/``, and
flatten each tree into a gene-level table (one row per organism gene entry).

By default the main classification tree ``<org>00001`` is processed. Use
``--brite-ids`` to target specific trees (for example ``spo01000,spo03000``) or
``--all-brite`` to process every BRITE tree KEGG lists for the organism.

Output
------
- ``<outdir>/derived/brite_flat.tsv`` (or ``.parquet``) with the columns
  ``BRITE_ID, Level_A, Level_B, Level_C, Level_D, Level_E, KO_ID, KO_Name,
  Gene_ID, Gene_Symbol, Gene_Description, EC_Number``.
- Raw JSON cached at ``<outdir>/raw/brite/<brite_id>.json``.

Usage
-----
    mamba run -n kegg_parser python scripts/fetch_brite.py --org spo
    mamba run -n kegg_parser python scripts/fetch_brite.py --org spo --format parquet
    mamba run -n kegg_parser python scripts/fetch_brite.py --org spo --brite-ids spo00001,spo01000
    mamba run -n kegg_parser python scripts/fetch_brite.py --org spo --all-brite --force --verbose

Author: Yusheng Yang (guidance) + Agent (implementation)
Date:   2026-09-22
Version: 1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# Standard library
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# Project path setup (must precede project imports)
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../src").resolve()))

# Project imports
from kegg_parser.brite import process_brite_trees, select_brite_ids  # noqa: E402
from kegg_parser.config import OutputFormat  # noqa: E402
from kegg_parser.utils import setup_logger  # noqa: E402

# Third-party
from loguru import logger  # noqa: E402

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
DEFAULT_OUTDIR_TEMPLATE = "data/{org}"

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================


@dataclass(kw_only=True, slots=True, frozen=True)
class BriteRunConfig:
    """Runtime configuration for a BRITE fetching run."""

    org: str
    outdir: Path
    fmt: OutputFormat
    brite_ids: list[str]
    force: bool


# =============================================================================
# CORE LOGIC
# =============================================================================
def _split_brite_ids(raw: str | None) -> list[str] | None:
    """Split a comma-separated BRITE id string into a clean list."""
    if not raw:
        return None
    return [item.strip() for item in raw.split(",") if item.strip()]


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Fetch and flatten KEGG BRITE trees for an organism.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--org", required=True, help="KEGG organism code (e.g. spo, hsa, eco).")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help=f"Output root directory (default: {DEFAULT_OUTDIR_TEMPLATE.format(org='<org>')}).",
    )
    parser.add_argument(
        "--format",
        type=OutputFormat,
        choices=list(OutputFormat),
        default=OutputFormat.TSV,
        help="Output table format (default: tsv).",
    )
    parser.add_argument(
        "--brite-ids",
        default=None,
        help="Comma-separated BRITE ids to process (default: <org>00001).",
    )
    parser.add_argument(
        "--all-brite",
        action="store_true",
        help="Process every BRITE tree KEGG lists for the organism.",
    )
    parser.add_argument("--force", action="store_true", help="Re-download raw data even when cached.")
    parser.add_argument("--verbose", action="store_true", help="Enable debug-level logging.")
    return parser.parse_args()


def main() -> int:
    """Main execution flow for the BRITE CLI."""
    args = parse_args()
    setup_logger(args.verbose)

    outdir = (args.outdir if args.outdir is not None else Path(DEFAULT_OUTDIR_TEMPLATE.format(org=args.org))).resolve()
    try:
        brite_ids = select_brite_ids(args.org, brite_ids=_split_brite_ids(args.brite_ids), all_brite=args.all_brite)
        config = BriteRunConfig(
            org=args.org,
            outdir=outdir,
            fmt=args.format,
            brite_ids=brite_ids,
            force=args.force,
        )
        logger.info(f"BRITE run started: org={config.org}, trees={config.brite_ids}, format={config.fmt}")

        output_path = process_brite_trees(
            config.org,
            config.outdir,
            fmt=config.fmt,
            brite_ids=config.brite_ids,
            force=config.force,
        )
        logger.success(f"BRITE processing finished: {output_path}")
    except Exception as error:
        logger.exception(f"BRITE processing failed: {error}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
