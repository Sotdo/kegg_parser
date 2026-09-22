#!/usr/bin/env python3

"""
Fetch and Summarise KEGG PATHWAY Annotations
============================================

Retrieve organism-specific pathway annotations via ``bioservices.kegg.KEGG``
and build two gene-centric tables:

- A gene-to-pathway mapping table (one row per ``Gene_ID``/``Pathway_ID`` pair)
  enriched with the pathway name and its KEGG class.
- A per-gene summary table aggregating all pathway ids, names and classes into
  semicolon-separated strings.

Raw KEGG responses are cached under ``<outdir>/raw/pathway/``.

Output
------
- ``<outdir>/derived/pathway_gene_mapping.tsv`` (or ``.parquet``) with columns
  ``Gene_ID, Pathway_ID, Pathway_Name, Pathway_Class``.
- ``<outdir>/derived/gene_pathway_summary.tsv`` (or ``.parquet``) with columns
  ``Gene_ID, Pathway_IDs, Pathway_Names, Pathway_Classes, Pathway_Count``.

Usage
-----
    mamba run -n kegg_parser python scripts/fetch_pathway.py --org spo
    mamba run -n kegg_parser python scripts/fetch_pathway.py --org spo --format parquet
    mamba run -n kegg_parser python scripts/fetch_pathway.py --org spo --force --verbose

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
from kegg_parser.config import OutputFormat  # noqa: E402
from kegg_parser.pathway import process_pathways  # noqa: E402
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
class PathwayRunConfig:
    """Runtime configuration for a pathway fetching run."""

    org: str
    outdir: Path
    fmt: OutputFormat
    force: bool


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Fetch and summarise KEGG PATHWAY annotations for an organism.",
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
    parser.add_argument("--force", action="store_true", help="Re-download raw data even when cached.")
    parser.add_argument("--verbose", action="store_true", help="Enable debug-level logging.")
    return parser.parse_args()


def main() -> int:
    """Main execution flow for the pathway CLI."""
    args = parse_args()
    setup_logger(args.verbose)

    outdir = (args.outdir if args.outdir is not None else Path(DEFAULT_OUTDIR_TEMPLATE.format(org=args.org))).resolve()
    config = PathwayRunConfig(org=args.org, outdir=outdir, fmt=args.format, force=args.force)

    try:
        logger.info(f"PATHWAY run started: org={config.org}, format={config.fmt}")
        mapping_path, summary_path = process_pathways(
            config.org,
            config.outdir,
            fmt=config.fmt,
            force=config.force,
        )
        logger.success(f"PATHWAY processing finished: {mapping_path} | {summary_path}")
    except Exception as error:
        logger.exception(f"PATHWAY processing failed: {error}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
