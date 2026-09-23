"""
KEGG Parser Shared Utilities
============================

Reusable helpers shared by the ``kegg_parser`` modules:

- Loguru logging setup for the CLI scripts.
- A ``stamina``-backed retry decorator (exponential backoff with jitter) for
  flaky network calls.
- Polite request pacing to respect KEGG's rate limits.
- HTTP fetching via ``httpx``.
- Atomic writers for raw payloads and a table writer/reader that supports both
  TSV and Parquet outputs.

This is a library module: it exposes utilities only, with no CLI.

Output
------
- ``Path`` objects pointing at raw cache files and written tables.
- ``pandas.DataFrame`` objects read back from disk.

Author: Yusheng Yang (guidance) + Agent (implementation)
Date:   2026-09-22
Version: 1.1.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# Standard library
import random
import sys
import time
from collections.abc import Callable
from pathlib import Path

# Data processing
import pandas as pd

# Third-party
import httpx
import stamina
from loguru import logger

# Project imports
from kegg_parser.config import (
    MAX_REQUEST_DELAY,
    MAX_RETRIES,
    MIN_REQUEST_DELAY,
    NA_VALUE,
    OutputFormat,
    REQUEST_CONNECT_TIMEOUT,
    REQUEST_POOL_TIMEOUT,
    REQUEST_READ_TIMEOUT,
    REQUEST_WRITE_TIMEOUT,
    RETRY_INITIAL_DELAY,
    RETRY_JITTER,
    RETRY_MAX_DELAY,
    RETRY_TIMEOUT,
    SYSTEMATIC_NAME_PREFIX,
    TEXT_ENCODING,
    USER_AGENT,
)

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
# ``requests`` (used internally by bioservices) subclasses IOError/OSError, so
# this tuple also covers bioservices network failures.
RETRYABLE_EXCEPTIONS: tuple[type[Exception], ...] = (httpx.HTTPError, OSError)


# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(verbose: bool = False) -> None:
    """Configure the global Loguru logger to write to stderr."""
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    logger.add(
        sys.stderr,
        level=level,
        format=(
            "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
            "<level>{level:<8}</level> | <cyan>{function}</cyan> - {message}"
        ),
    )


# =============================================================================
# CORE LOGIC
# =============================================================================
def retryable[T](func: Callable[..., T]) -> Callable[..., T]:
    """Decorate a function with stamina retries using the shared backoff policy."""
    return stamina.retry(
        on=RETRYABLE_EXCEPTIONS,
        attempts=MAX_RETRIES,
        timeout=RETRY_TIMEOUT,
        wait_initial=RETRY_INITIAL_DELAY,
        wait_max=RETRY_MAX_DELAY,
        wait_jitter=RETRY_JITTER,
    )(func)


def polite_delay(min_delay: float = MIN_REQUEST_DELAY, max_delay: float = MAX_REQUEST_DELAY) -> None:
    """Sleep for a random short interval to respect KEGG rate limits."""
    time.sleep(random.uniform(min_delay, max_delay))


def strip_systematic_prefix(symbol: str, prefix: str = SYSTEMATIC_NAME_PREFIX) -> str:
    """Strip the organism systematic-name prefix (e.g. ``SPOM_``) from a gene symbol."""
    return symbol.removeprefix(prefix)


@retryable
def http_get_bytes(
    url: str,
    params: dict[str, str] | None = None,
    timeout: httpx.Timeout | None = None,
) -> bytes:
    """Fetch a URL and return the full body as bytes, with retries and pacing."""
    polite_delay()
    logger.debug(f"GET {url} params={params}")
    request_timeout = timeout or httpx.Timeout(
        connect=REQUEST_CONNECT_TIMEOUT,
        read=REQUEST_READ_TIMEOUT,
        write=REQUEST_WRITE_TIMEOUT,
        pool=REQUEST_POOL_TIMEOUT,
    )
    response = httpx.get(
        url,
        params=params,
        timeout=request_timeout,
        follow_redirects=True,
        headers={"User-Agent": USER_AGENT},
    )
    response.raise_for_status()
    return response.content


def ensure_dir(path: Path) -> Path:
    """Create a directory (and parents) if needed and return it."""
    path.mkdir(parents=True, exist_ok=True)
    return path


def atomic_write_bytes(path: Path, data: bytes) -> Path:
    """Write bytes to ``path`` atomically via a temporary sibling file."""
    ensure_dir(path.parent)
    temp_path = path.with_name(f"{path.name}.tmp")
    temp_path.write_bytes(data)
    temp_path.replace(path)
    return path


def atomic_write_text(path: Path, text: str) -> Path:
    """Write text to ``path`` atomically via a temporary sibling file."""
    ensure_dir(path.parent)
    temp_path = path.with_name(f"{path.name}.tmp")
    temp_path.write_text(text, encoding=TEXT_ENCODING)
    temp_path.replace(path)
    return path


def write_table(frame: pd.DataFrame, stem_path: Path, fmt: OutputFormat | str = OutputFormat.TSV) -> Path:
    """Write a DataFrame to ``<stem_path>.<fmt>`` and return the written path."""
    resolved = OutputFormat(fmt)
    ensure_dir(stem_path.parent)
    if resolved is OutputFormat.PARQUET:
        output_path = stem_path.with_name(f"{stem_path.name}.parquet")
        frame.to_parquet(output_path, index=False)
    else:
        output_path = stem_path.with_name(f"{stem_path.name}.tsv")
        frame.to_csv(output_path, sep="\t", index=False, na_rep=NA_VALUE)
    logger.info(f"Wrote {len(frame):,} rows x {len(frame.columns)} columns to {output_path}")
    return output_path


def read_table(path: Path, fmt: OutputFormat | str | None = None) -> pd.DataFrame:
    """Read a TSV or Parquet table, inferring the format from the suffix when omitted."""
    path = Path(path)
    resolved = OutputFormat(fmt) if fmt is not None else OutputFormat(path.suffix.lstrip("."))
    if resolved is OutputFormat.PARQUET:
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[NA_VALUE])
