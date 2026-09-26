# AGENTS.md

KEGG BRITE / PATHWAY fetcher that flattens KEGG data into gene-centric
TSV or Parquet tables. See `README.md` for output schemas.

## Environment & commands

All tooling lives in the `kegg_parser` mamba/conda env (Python 3.12). The system
`python` does **not** have the deps — always prefix commands:

```bash
mamba run -n kegg_parser python -m pytest              # full suite (22 offline tests)
mamba run -n kegg_parser python -m pytest tests/test_brite.py::test_flatten_brite_tree_parses_levels_and_fields
mamba run -n kegg_parser python -m ruff check .        # lint (must pass)
mamba run -n kegg_parser python scripts/run_kegg_pipeline.py --org spo   # full pipeline
```

There is no CI, pre-commit, typecheck, or codegen. `ruff` and `pytest` are the
only gates. Run lint + tests after changes.

## Layout

- `src/kegg_parser/` — library only, no CLI. `config.py` (endpoints, retry/pacing,
  output column enums), `utils.py` (logging, retry, atomic/table I/O), `kegg_rest.py`
  (KEGG REST client + text parsing), `brite.py`, `pathway.py`.
- `scripts/` — thin argparse CLIs. Each inserts `../src` into `sys.path`, so scripts
  and tests run without `pip install`; imports after that setup carry `# noqa: E402`.
- `tests/` — `conftest.py` adds `src` to `sys.path`. Tests use synthetic fixtures and
  never hit the network; keep new tests offline.
- `data/` and `tmp/` are gitignored. Default output root is `data/<org>`; raw responses
  cache under `<outdir>/raw/`, derived tables under `<outdir>/derived/`. `--force`
  bypasses the cache.

## Conventions

- `pyproject.toml` sets ruff `line-length = 120`, select `E,F,UP,B,T201`; **isort (`I`)
  is intentionally disabled** (import order is hand-maintained) — do not "fix" imports.
  `T201` forbids `print`; use `loguru` `logger` (set up via `setup_logger`).
- Files follow a banner-section layout: IMPORTS, GLOBAL CONSTANTS & ENUMS,
  CONFIGURATION & DATACLASSES (scripts), CORE LOGIC, MAIN EXECUTION. Preserve it.
- Public entrypoints: `process_brite_trees` and `process_pathways`, re-exported from
  `kegg_parser`. Network calls go through `utils.retryable` / `polite_delay`; top-level
  `process_*` are wrapped in `@logger.catch(reraise=True)`.

## Gotchas

- All KEGG REST access goes through `utils.http_get_bytes` (retryable + paced);
  `kegg_rest.fetch_kegg_text` wraps it for text endpoints. Add new endpoints there
  rather than calling `httpx` directly. `bioservices` is no longer a dependency.
- BRITE discovers all trees via `list/brite/<org>` with a fallback to `<org>00001`.
- `stamina` retries wrap network errors (`httpx.HTTPError`, `OSError`); `OSError` is
  kept for low-level socket failures surfaced by the underlying network stack.
