#!/usr/bin/env python3

"""
Run the Full KEGG Parsing Pipeline
==================================

Convenience entry point that runs the complete ``kegg_parser`` workflow for a
single organism: it first flattens the requested KEGG BRITE trees (all of them
by default), then fetches and summarises the KEGG PATHWAY annotations. Raw
responses and derived tables are written under a single output root (default
``./data/<org>``).

Output
------
- ``<outdir>/derived/brite_flat.tsv`` (or ``.parquet``).
- ``<outdir>/derived/pathway_gene_mapping.tsv`` (or ``.parquet``).
- ``<outdir>/derived/gene_pathway_summary.tsv`` (or ``.parquet``).
- Raw caches under ``<outdir>/raw/{brite,pathway,gene}/``.

Usage
-----
    mamba run -n kegg_parser python scripts/run_kegg_pipeline.py --org spo
    mamba run -n kegg_parser python scripts/run_kegg_pipeline.py --org spo --format parquet
    mamba run -n kegg_parser python scripts/run_kegg_pipeline.py --org spo --brite-ids spo00001,spo03000
    mamba run -n kegg_parser python scripts/run_kegg_pipeline.py --org spo --force --verbose

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
class PipelineConfig:
    """Runtime configuration for the combined KEGG pipeline."""

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
        description="Run the full KEGG BRITE + PATHWAY parsing pipeline for an organism.",
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
        help="Comma-separated BRITE ids to process (default: every BRITE tree for the organism).",
    )
    parser.add_argument("--force", action="store_true", help="Re-download raw data even when cached.")
    parser.add_argument("--verbose", action="store_true", help="Enable debug-level logging.")
    return parser.parse_args()


def main() -> int:
    """Main execution flow for the combined pipeline CLI."""
    args = parse_args()
    setup_logger(args.verbose)

    outdir = (args.outdir if args.outdir is not None else Path(DEFAULT_OUTDIR_TEMPLATE.format(org=args.org))).resolve()
    try:
        brite_ids = select_brite_ids(args.org, brite_ids=_split_brite_ids(args.brite_ids))
        config = PipelineConfig(
            org=args.org,
            outdir=outdir,
            fmt=args.format,
            brite_ids=brite_ids,
            force=args.force,
        )
        logger.info(f"KEGG pipeline started: org={config.org}, trees={config.brite_ids}, format={config.fmt}")

        brite_path = process_brite_trees(
            config.org,
            config.outdir,
            fmt=config.fmt,
            brite_ids=config.brite_ids,
            force=config.force,
        )
        mapping_path, summary_path = process_pathways(
            config.org,
            config.outdir,
            fmt=config.fmt,
            force=config.force,
        )
        logger.success(
            f"KEGG pipeline finished: BRITE -> {brite_path}; "
            f"PATHWAY -> {mapping_path} | {summary_path}"
        )
    except Exception as error:
        logger.exception(f"KEGG pipeline failed: {error}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
