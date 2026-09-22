"""
KEGG Parser Shared Utilities
============================

Reusable helpers shared by the ``kegg_parser`` modules:

- Loguru logging setup for the CLI scripts.
- A retry-with-exponential-backoff decorator for flaky network calls.
- Polite request pacing to respect KEGG's rate limits.
- Atomic writers for raw payloads and a table writer/reader that supports both
  TSV and Parquet outputs.

This is a library module: it exposes utilities only, with no CLI.

Output
------
- ``Path`` objects pointing at raw cache files and written tables.
- ``pandas.DataFrame`` objects read back from disk.

Author: Yusheng Yang (guidance) + Agent (implementation)
Date:   2026-09-22
Version: 1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# Standard library
import functools
import random
import sys
import time
from collections.abc import Callable
from pathlib import Path
from typing import TypeVar

# Data processing
import pandas as pd

# Third-party
import requests
from loguru import logger

# Project imports
from kegg_parser.config import (
    DOWNLOAD_CHUNK_SIZE,
    MAX_REQUEST_DELAY,
    MAX_RETRIES,
    MIN_REQUEST_DELAY,
    NA_VALUE,
    OutputFormat,
    REQUEST_TIMEOUT,
    RETRY_BACKOFF_FACTOR,
    RETRY_BASE_DELAY,
    RETRY_JITTER,
    TEXT_ENCODING,
    USER_AGENT,
)

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
T = TypeVar("T")
RETRYABLE_EXCEPTIONS: tuple[type[Exception], ...] = (
    requests.RequestException,
    ConnectionError,
    TimeoutError,
    OSError,
)


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
def retry_with_backoff(
    max_retries: int = MAX_RETRIES,
    base_delay: float = RETRY_BASE_DELAY,
    backoff_factor: float = RETRY_BACKOFF_FACTOR,
    jitter: float = RETRY_JITTER,
    exceptions: tuple[type[Exception], ...] = RETRYABLE_EXCEPTIONS,
) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """Decorate a function with retry logic using exponential backoff."""

    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @functools.wraps(func)
        def wrapper(*args: object, **kwargs: object) -> T:
            for attempt in range(1, max_retries + 1):
                try:
                    return func(*args, **kwargs)
                except exceptions as error:
                    if attempt >= max_retries:
                        logger.error(f"{func.__name__} failed after {max_retries} attempt(s): {error}")
                        raise
                    delay = base_delay * (backoff_factor ** (attempt - 1)) + random.uniform(0.0, jitter)
                    logger.warning(
                        f"{func.__name__} attempt {attempt}/{max_retries} failed ({error}); retrying in {delay:.1f}s"
                    )
                    time.sleep(delay)
            raise RuntimeError(f"{func.__name__} exhausted retries unexpectedly")

        return wrapper

    return decorator


def polite_delay(min_delay: float = MIN_REQUEST_DELAY, max_delay: float = MAX_REQUEST_DELAY) -> None:
    """Sleep for a random short interval to respect KEGG rate limits."""
    time.sleep(random.uniform(min_delay, max_delay))


@retry_with_backoff()
def http_get_bytes(
    url: str,
    params: dict[str, str] | None = None,
    timeout: tuple[float, float] = REQUEST_TIMEOUT,
) -> bytes:
    """Fetch a URL and return the full body as bytes, with retries and pacing."""
    polite_delay()
    logger.debug(f"GET {url} params={params}")
    with requests.get(
        url,
        params=params,
        timeout=timeout,
        headers={"User-Agent": USER_AGENT},
        stream=True,
    ) as response:
        response.raise_for_status()
        buffer = bytearray()
        for chunk in response.iter_content(chunk_size=DOWNLOAD_CHUNK_SIZE):
            if chunk:
                buffer.extend(chunk)
    return bytes(buffer)


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
