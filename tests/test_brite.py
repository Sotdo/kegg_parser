"""
Unit tests for ``kegg_parser.brite``
====================================

Exercise the recursive BRITE flattening logic on a synthetic tree that covers
the real KEGG variations: tab-merged gene/KO leaves, KO-only leaves, reversed
field order, shallow branches, forward-filled levels and pure classification
labels.

Author: Yusheng Yang (guidance) + Agent (implementation)
Date:   2026-09-22
Version: 1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import kegg_parser.brite as brite_module
from kegg_parser.brite import (
    _parse_ec,
    _parse_gene,
    _parse_ko,
    default_brite_ids,
    flatten_brite_tree,
    select_brite_ids,
)
from kegg_parser.config import BRITE_COLUMNS

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
SYNTHETIC_TREE: dict = {
    "name": "spo00001",
    "children": [
        {
            "name": "09100 Metabolism",
            "children": [
                {
                    "name": "09101 Carbohydrate metabolism",
                    "children": [
                        {
                            "name": "00010 Glycolysis / Gluconeogenesis [PATH:spo00010]",
                            "children": [
                                {"name": "2542634 hxk1; hexokinase 1\tK00844 HK; hexokinase [EC:2.7.1.1]"},
                                {
                                    "name": (
                                        "2539625 pgi1; glucose-6-phosphate isomerase\t"
                                        "K01810 GPI; glucose-6-phosphate isomerase [EC:5.3.1.9]"
                                    )
                                },
                                {"name": "K99999 SOME; enzyme absent from this organism"},
                            ],
                        },
                        {"name": "Empty category leaf"},
                    ],
                },
                {
                    "name": "09102 Energy metabolism",
                    "children": [
                        {"name": "2539000 atp1; ATP synthase\tK02132 ATPF1A; ATP synthase [EC:3.6.3.14]"},
                        {"name": "K00844 HK; hexokinase\t2542634 hxk1; hexokinase 1"},
                    ],
                },
            ],
        }
    ],
}


# =============================================================================
# CORE LOGIC
# =============================================================================
def test_flatten_brite_tree_has_expected_schema() -> None:
    """The flattened frame must match the declared column schema."""
    frame = flatten_brite_tree(SYNTHETIC_TREE, "spo00001")
    assert list(frame.columns) == BRITE_COLUMNS


def test_flatten_brite_tree_parses_levels_and_fields() -> None:
    """A tab-merged gene/KO leaf expands into the expected levels and fields."""
    frame = flatten_brite_tree(SYNTHETIC_TREE, "spo00001")
    row = frame[frame["Gene_ID"] == "2542634"].iloc[0]

    assert row["BRITE_ID"] == "spo00001"
    assert row["Level_A"] == "09100 Metabolism"
    assert row["Level_B"] == "09101 Carbohydrate metabolism"
    assert row["Level_C"] == "00010 Glycolysis / Gluconeogenesis [PATH:spo00010]"
    assert row["Level_D"] == "00010 Glycolysis / Gluconeogenesis [PATH:spo00010]"
    assert row["Level_E"] == "K00844 HK; hexokinase [EC:2.7.1.1]"
    assert row["Level_F"] == "2542634 hxk1; hexokinase 1"
    assert row["KO_ID"] == "K00844"
    assert row["KO_Name"] == "HK; hexokinase [EC:2.7.1.1]"
    assert row["Gene_Symbol"] == "hxk1"
    assert row["Gene_Description"] == "hexokinase 1"
    assert row["EC_Number"] == "2.7.1.1"


def test_flatten_brite_tree_forward_fills_shallow_branch() -> None:
    """A shallow gene leaf fills empty levels C/D with the previous level."""
    frame = flatten_brite_tree(SYNTHETIC_TREE, "spo00001")
    row = frame[frame["Gene_ID"] == "2539000"].iloc[0]

    assert row["Level_B"] == "09102 Energy metabolism"
    assert row["Level_C"] == "09102 Energy metabolism"
    assert row["Level_D"] == "09102 Energy metabolism"


def test_flatten_brite_tree_forward_fills_empty_nodes() -> None:
    """Empty intermediate classification nodes inherit the previous level."""
    tree = {
        "name": "spo00001",
        "children": [
            {
                "name": "Top",
                "children": [
                    {
                        "name": "",
                        "children": [
                            {
                                "name": "Leaf",
                                "children": [
                                    {"name": "2542634 hxk1; hexokinase 1\tK00844 HK; hexokinase [EC:2.7.1.1]"}
                                ],
                            }
                        ],
                    }
                ],
            }
        ],
    }
    frame = flatten_brite_tree(tree, "spo00001")
    row = frame.iloc[0]
    assert row["Level_A"] == "Top"
    assert row["Level_B"] == "Top"
    assert row["Level_C"] == "Leaf"
    assert row["Level_D"] == "Leaf"


def test_flatten_brite_tree_skips_pure_labels() -> None:
    """Classification labels without a KO or gene must be skipped."""
    frame = flatten_brite_tree(SYNTHETIC_TREE, "spo00001")
    assert "Empty category leaf" not in set(frame["Level_C"]) | set(frame["Level_D"])
    assert len(frame) == 5


def test_flatten_brite_tree_keeps_ko_only_leaf() -> None:
    """A leaf that only defines a KO is kept with an empty gene id and Level_F."""
    frame = flatten_brite_tree(SYNTHETIC_TREE, "spo00001")
    row = frame[frame["KO_ID"] == "K99999"].iloc[0]
    assert row["Gene_ID"] == ""
    assert row["KO_Name"] == "SOME; enzyme absent from this organism"
    assert row["Level_E"] == "K99999 SOME; enzyme absent from this organism"
    assert row["Level_F"] == ""


def test_flatten_brite_tree_joins_deep_levels_into_level_d() -> None:
    """Trees deeper than A/B/C/D keep the extra classification inside Level_D."""
    tree = {
        "name": "spo01000",
        "children": [
            {
                "name": "1. Oxidoreductases",
                "children": [
                    {
                        "name": "1.1 Acting on the CH-OH group of donors",
                        "children": [
                            {
                                "name": "1.1.1 With NAD+ or NADP+ as acceptor",
                                "children": [
                                    {
                                        "name": "1.1.1.1 alcohol dehydrogenase",
                                        "children": [
                                            {
                                                "name": (
                                                    "2538902 adh1; alcohol dehydrogenase Adh1\t"
                                                    "K13953 adhP; alcohol dehydrogenase [EC:1.1.1.1]"
                                                )
                                            }
                                        ],
                                    }
                                ],
                            }
                        ],
                    }
                ],
            }
        ],
    }
    frame = flatten_brite_tree(tree, "spo01000")
    row = frame.iloc[0]
    assert row["Level_A"] == "1. Oxidoreductases"
    assert row["Level_B"] == "1.1 Acting on the CH-OH group of donors"
    assert row["Level_C"] == "1.1.1 With NAD+ or NADP+ as acceptor"
    assert row["Level_D"] == "1.1.1.1 alcohol dehydrogenase"
    assert row["KO_ID"] == "K13953"
    assert row["Gene_ID"] == "2538902"


def test_flatten_brite_tree_handles_reversed_field_order() -> None:
    """KO-first leaves are parsed correctly regardless of tab order."""
    frame = flatten_brite_tree(SYNTHETIC_TREE, "spo00001")
    reversed_rows = frame[(frame["KO_ID"] == "K00844") & (frame["Gene_ID"] == "2542634")]
    assert len(reversed_rows) == 2


def test_select_brite_ids_prefers_explicit_ids() -> None:
    """Explicit ids win and blanks are dropped, without any network access."""
    assert select_brite_ids("spo", ["spo00001", "", "spo03000"]) == ["spo00001", "spo03000"]


def test_default_brite_ids_discovers_all(monkeypatch) -> None:
    """Without explicit ids, every advertised tree is selected."""
    monkeypatch.setattr(brite_module, "list_organism_brite_ids", lambda org: {"spo00001": "a", "spo03000": "b"})
    assert default_brite_ids("spo") == ["spo00001", "spo03000"]
    assert select_brite_ids("spo") == ["spo00001", "spo03000"]


def test_default_brite_ids_falls_back_on_discovery_error(monkeypatch) -> None:
    """A discovery failure falls back to the main tree instead of crashing."""

    def _boom(org: str) -> dict[str, str]:
        raise RuntimeError("offline")

    monkeypatch.setattr(brite_module, "list_organism_brite_ids", _boom)
    assert default_brite_ids("spo") == ["spo00001"]


def test_parse_ko_and_gene_helpers() -> None:
    """The low-level field parsers split text as documented."""
    assert _parse_ko("K00844 HK; hexokinase [EC:2.7.1.1]") == ("K00844", "HK; hexokinase [EC:2.7.1.1]")
    assert _parse_ko("not-a-ko") == ("", "not-a-ko")
    assert _parse_gene("hsa:3098 HK1; hexokinase-1") == ("3098", "HK1", "hexokinase-1")
    assert _parse_gene("2542569 SPOM_SPAC186.08C; L-lactate dehydrogenase") == (
        "2542569",
        "SPAC186.08C",
        "L-lactate dehydrogenase",
    )
    assert _parse_gene("") == ("", "", "")
    assert _parse_ec("K00001 X; enzyme [EC:1.2.3.4]", "gene [EC:5.6.7.-]") == "1.2.3.4;5.6.7.-"
