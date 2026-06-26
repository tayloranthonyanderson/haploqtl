"""Equivalence of ``haploqtl.contrast`` to the original R contrast functions.

``tests/fixtures/eb9_contrast_golden.csv`` was produced by ``tests/r/contrast_reference.R``
(a base-R reproduction of the contrast logic in the original ``visualize_haplotypes.Rmd``,
Anderson et al. 2024) run on the committed cluster table ``eb9_clusters_subset.csv``. This
test asserts the Python contrast reproduces it window-for-window — so CI proves faithfulness
to the published method with no R dependency.
"""

from pathlib import Path

import pandas as pd

from haploqtl.contrast import contrast_oneway, contrast_twoway

FIX = Path(__file__).parent / "fixtures"
RESISTANT = [
    "191163",
    "191164",
    "191167",
    "191172",
    "191175",
    "191174H",
    "CU3AllData",
    "191357",
    "201041",
]
SUSCEPTIBLE = ["191165", "ERR418112", "SRR1572598"]


def _clusters() -> pd.DataFrame:
    return pd.read_csv(FIX / "eb9_clusters_subset.csv", dtype={"sample": str})


def _golden() -> pd.DataFrame:
    return pd.read_csv(FIX / "eb9_contrast_golden.csv")


def _merge(py: pd.DataFrame, golden: pd.DataFrame, col: str) -> pd.DataFrame:
    py = py.assign(position=py["position"].round().astype(int))
    py["py"] = py["diagnostic"].fillna(False).astype(int)
    g = golden.assign(Positions=golden["Positions"].astype(int))
    merged = py.merge(g, left_on="position", right_on="Positions", how="outer")
    assert merged["py"].notna().all() and merged[col].notna().all(), "window sets differ"
    return merged


def test_oneway_matches_r_golden():
    merged = _merge(contrast_oneway(_clusters(), RESISTANT), _golden(), "oneway")
    assert (merged["py"] == merged["oneway"]).all()


def test_twoway_matches_r_golden():
    merged = _merge(contrast_twoway(_clusters(), RESISTANT, SUSCEPTIBLE), _golden(), "twoway")
    assert (merged["py"] == merged["twoway"]).all()
    assert int(merged["py"].sum()) == 4  # the EB-9 diagnostic block, 4 windows
