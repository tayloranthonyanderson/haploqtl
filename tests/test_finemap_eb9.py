"""Locks the EB-9 fine-mapping case study (examples/finemap_eb9.py).

Reproduces the paper's ~70% interval and the tighter ~81% interval obtained by adding the
phenotyped-susceptible recombinant Ailsa Craig, and checks the published phenotype classes
that make the argument non-circular.
"""

from pathlib import Path

import pandas as pd

from haploqtl.contrast import contrast_twoway
from haploqtl.introgression import call_interval, interval_reduction

FIX = Path(__file__).parent / "fixtures"
PHENO = Path(__file__).parent.parent / "data" / "eb9_phenotypes.csv"
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
AILSA = "SRR1572460"  # Ailsa Craig
BGV = "SRR7279711"  # BGV007871
PRIOR = (61_819_509, 63_679_761)
WINDOW = 250_000


def _clusters() -> pd.DataFrame:
    return pd.read_csv(FIX / "eb9_clusters_finemap.csv", dtype={"sample": str})


def _call(susceptible: list[str]):
    call = call_interval(
        contrast_twoway(_clusters(), RESISTANT, susceptible), window=WINDOW, max_gap=1
    )
    return call, interval_reduction((call.start, call.stop), PRIOR).reduction


def test_paper_baseline_is_70pct():
    call, red = _call(SUSCEPTIBLE)
    assert (call.start, call.stop) == (62_400_075, 62_950_075)
    assert 0.69 <= red <= 0.72


def test_ailsa_craig_narrows_to_81pct():
    call, red = _call(SUSCEPTIBLE + [AILSA])
    assert (call.start, call.stop) == (62_600_075, 62_950_075)  # left edge trimmed 62.40 -> 62.60
    assert 0.80 <= red <= 0.82


def test_published_phenotypes_make_the_argument():
    pheno = pd.read_csv(PHENO).set_index("vcf_id")
    assert pheno.loc[AILSA, "class"] == "susceptible"  # the recombinant that narrows (left edge)
    assert (
        pheno.loc[BGV, "class"] == "resistant"
    )  # same haplotype, opposite phenotype (the confound)
