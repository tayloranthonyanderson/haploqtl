"""Locks the EB-5 generalization case study (examples/finemap_eb5.py).

EB-5 is a different QTL, donor (Hawaii 7998), and trait (foliar, not collar rot) than EB-9.
Reproducing its published ~56% narrowing from the same engine — with a disjoint resistant/
susceptible set — is the evidence that the method is not specialized to EB-9.
"""

from pathlib import Path

import pandas as pd

from haploqtl.contrast import contrast_twoway
from haploqtl.introgression import call_interval, interval_reduction

FIX = Path(__file__).parent / "fixtures"
PHENO = Path(__file__).parent.parent / "data" / "eb9_phenotypes.csv"

DONOR = "HA7998"  # Hawaii 7998
RESISTANT = ["HA7998", "OH7536", "OH7663"]  # Hawaii 7998 -> OH7536 -> OH08-7663
SUSCEPTIBLE = ["OH88119", "SRR1572598", "ERR418112"]  # OH88119, NC 84173, Brandywine
PRIOR = (62_700_265, 63_842_577)  # Anderson et al. (2024), Table S2
WINDOW = 250_000


def _clusters() -> pd.DataFrame:
    return pd.read_csv(FIX / "eb5_clusters_finemap.csv", dtype={"sample": str})


def _call(resistant: list[str], susceptible: list[str]):
    return call_interval(
        contrast_twoway(_clusters(), resistant, susceptible), window=WINDOW, max_gap=1
    )


def test_eb5_reproduces_paper_56pct():
    call = _call(RESISTANT, SUSCEPTIBLE)
    assert call is not None
    assert (call.start, call.stop) == (62_850_818, 63_350_818)
    red = interval_reduction((call.start, call.stop), PRIOR)
    assert red.refined_span == 500_000  # paper Table S2: 500,126 bp
    assert 0.55 <= red.reduction <= 0.57  # paper: 56.2%


def test_oh7536_constrains_the_left_boundary():
    """Dropping OH7536 from the pathway widens the interval — it is an informative recombinant."""
    full = _call(RESISTANT, SUSCEPTIBLE)
    without = _call(["HA7998", "OH7663"], SUSCEPTIBLE)
    assert without is not None and full is not None
    assert without.start < full.start  # OH7536 trims the left edge


def test_oh7663_is_foliar_resistant_but_collar_rot_susceptible():
    """The EB-5 donor OH08-7663 carries foliar (EB-5) resistance yet is stem (EB-9) susceptible —
    EB-5 and EB-9 are independent traits and donors, which is why their donor sets barely overlap.
    """
    call = _call(RESISTANT, SUSCEPTIBLE)
    assert "OH7663" in RESISTANT  # in the EB-5 (foliar) resistant pathway
    assert call is not None
    pheno = pd.read_csv(PHENO).set_index("vcf_id")
    assert pheno.loc["OH7663", "class"] == "susceptible"  # EB-9 collar rot
    assert pheno.loc["OH7663", "mean_stem_lesion_pct"] > 50
