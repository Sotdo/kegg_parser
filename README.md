# kegg_parser

Fetch and flatten **KEGG BRITE** hierarchy trees and **KEGG PATHWAY**
annotations for a given organism (e.g. `spo`, `hsa`, `eco`) into gene-centric **TSV** or **Parquet**
tables.

The package handles two data products:

| Product | Source | Output |
|---|---|---|
| BRITE | `download_htext?format=json` | `derived/brite_flat.<fmt>` |
| PATHWAY | KEGG REST API | `derived/pathway_gene_mapping.<fmt>`, `derived/gene_pathway_summary.<fmt>` |

Raw KEGG responses are cached under `<outdir>/raw/` so repeated runs avoid re-querying the API.

## Installation

```bash
mamba env create -f environment.yml
mamba activate kegg_parser
```

`environment.yml` installs the package in editable mode (`pip install -e .`). To install into an
existing environment instead:

```bash
pip install -e .
```

## Quick start

```bash
# Full pipeline (BRITE + PATHWAY) for Schizosaccharomyces pombe.
# By default every BRITE tree KEGG advertises for the organism is integrated.
mamba run -n kegg_parser python scripts/run_kegg_pipeline.py --org spo

# Only BRITE, as Parquet, restricted to specific trees
mamba run -n kegg_parser python scripts/fetch_brite.py \
    --org spo --format parquet --brite-ids spo00001,spo01000

# Only PATHWAY
mamba run -n kegg_parser python scripts/fetch_pathway.py --org spo
```

### CLI arguments

| Argument | Applies to | Description |
|---|---|---|
| `--org` | all | KEGG organism code (required). |
| `--outdir` | all | Output root (default `./data/<org>`). |
| `--format` | all | `tsv` (default) or `parquet`. |
| `--brite-ids` | BRITE, pipeline | Comma-separated BRITE ids (default: every BRITE tree KEGG advertises for the organism). |
| `--force` | all | Ignore the raw cache and re-download. |
| `--verbose` | all | Debug-level logging. |

## Output schema

### `brite_flat`

One row per terminal node that resolves to a gene or KO entry.

| Column | Meaning |
|---|---|
| `BRITE_ID` | Source tree id (e.g. `spo00001`). |
| `Level_A` / `Level_B` / `Level_C` / `Level_D` | Ancestor labels from the root downward. Any classification deeper than four levels is joined into `Level_D` with `" > "`; empty classification levels are forward-filled with the previous level so A–D are gap-free. |
| `Level_E` | KO entry text (`K00844 HK; hexokinase [EC:2.7.1.1]`). |
| `Level_F` | Organism gene entry text (`2542634 hxk1; hexokinase 1`). |
| `KO_ID` / `KO_Name` | Parsed from `Level_E`. |
| `Gene_ID` / `Gene_Symbol` / `Gene_Description` | Parsed from `Level_F`; the KEGG/species prefix is stripped from `Gene_ID` and the `SPOM_` systematic-name prefix from `Gene_Symbol`. |
| `EC_Number` | All EC numbers found in the leaf, `;`-joined. |

Terminal classification labels without a KO or gene are skipped; shallow branches never raise.

### `pathway_gene_mapping`

`Gene_ID`, `Gene_Symbol`, `Gene_Description`, `Pathway_ID`, `Pathway_Name`, `Level_A`, `Level_B` — one row
per unique gene/pathway pair. `Level_A` / `Level_B` are the two levels of the KEGG pathway class
(from `br08901`). Gene symbol and description are resolved from the shared KEGG
`list/<org>` gene list (cached under `raw/gene/<org>_gene_list.txt`); the `SPOM_` systematic-name
prefix is stripped from `Gene_Symbol`.

### `gene_pathway_summary`

`Gene_ID`, `Gene_Symbol`, `Gene_Description`, `Pathway_IDs`, `Pathway_Names`, `Level_As`, `Level_Bs`
(`;`-joined) plus `Pathway_Count`.

## Project structure

```
kegg_parser/
├── environment.yml
├── pyproject.toml
├── README.md
├── src/kegg_parser/
│   ├── config.py     # endpoints, delays, retry policy, column enums
│   ├── utils.py      # logging, retry decorator, pacing, table I/O
│   ├── kegg_rest.py  # shared KEGG REST client, caching, TSV/gene-list parsing
│   ├── brite.py      # BRITE download, tree traversal, flattening
│   └── pathway.py    # PATHWAY fetching and aggregation
├── scripts/
│   ├── fetch_brite.py
│   ├── fetch_pathway.py
│   └── run_kegg_pipeline.py
└── tests/
```

## Python API

```python
from pathlib import Path
from kegg_parser import process_brite_trees, process_pathways, OutputFormat

outdir = Path("data/spo")
process_brite_trees("spo", outdir, fmt=OutputFormat.PARQUET)
process_pathways("spo", outdir, fmt=OutputFormat.PARQUET)
```

## Rate limiting & retries

KEGG REST requests are paced with a random `0.3–0.5 s` delay and wrapped by a `stamina` retry
decorator (exponential backoff with jitter, 3 attempts, 60 s cap). HTTP fetching uses `httpx` with a
long read timeout for the large BRITE JSON downloads.

## Tests

```bash
pytest
```
