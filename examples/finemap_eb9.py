"""Fine-mapping EB-9 further with a phenotyped recombinant.

Reproduces, from the bundled 780-genome chr09 fixture:

  1. the published EB-9 interval and ~70% reduction (resistant pedigree vs. susceptible
     controls);
  2. a tighter interval (~81%) by adding *Ailsa Craig* — a phenotyped-susceptible line that
     carries the Devon Surprise haplotype across the left flank, so that flank cannot be
     sufficient for resistance (a robust susceptible-recombinant argument);
  3. the limit: a "resistant recombinant" (BGV007871) carries Ailsa Craig's haplotype yet is
     resistant — its resistance comes from another source — so it cannot trim the right edge.

Stem-lesion phenotypes (mean % diseased stem, in planta) are from Anderson et al. (2024),
*The Plant Journal*; the sequence data is the bundled fixture.

Run:  uv run python examples/finemap_eb9.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from haploqtl.cluster import cluster_haplotypes
from haploqtl.contrast import compare_to_benchmark, contrast_twoway
from haploqtl.introgression import IntrogressionCall, call_interval, interval_reduction
from haploqtl.io import load_genotypes

ROOT = Path(__file__).resolve().parent.parent
VCF = ROOT / "data" / "SL4.0ch09_subset.vcf.gz"
PHENO = ROOT / "data" / "eb9_phenotypes.csv"

# The EB-9 introgression pathway (Devon Surprise -> Campbell 1943 -> NC EBR 2 -> NC 1 CELBR ...)
# and the susceptible controls, per Anderson et al. (2024).
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
DONOR = "191163"  # Devon Surprise
AILSA = "SRR1572460"  # Ailsa Craig (LA2838A)
BGV = "SRR7279711"  # BGV007871
PRIOR = (61_819_509, 63_679_761)  # prior EB-9 interval (the reduction denominator)
WINDOW = 250_000


def narrow(clusters: pd.DataFrame, susceptible: list[str], label: str) -> IntrogressionCall:
    """Call the EB-9 interval for a given susceptible set and print the % reduction."""
    track = contrast_twoway(clusters, RESISTANT, susceptible)
    call = call_interval(track, window=WINDOW, max_gap=1)
    assert call is not None, f"no diagnostic interval for {label!r}"
    red = interval_reduction((call.start, call.stop), PRIOR)
    print(
        f"  {label:<34} {call.start:>12,.0f} - {call.stop:<12,.0f} "
        f"({call.span:,.0f} bp, {red.reduction * 100:.1f}% reduction)"
    )
    return call


def main() -> None:
    pheno = pd.read_csv(PHENO).set_index("vcf_id")
    data = load_genotypes(VCF)
    print(f"Clustering {data.n_samples} genomes over chr09 (paper params: 250 kb / 100 kb)...\n")
    clusters = cluster_haplotypes(
        data, "ch09", window=WINDOW, step=100_000, min_snps=10, d_grid=np.arange(2.0, 80.0, 10.0)
    )

    print("EB-9 interval (resistant pedigree vs. susceptible set):")
    narrow(clusters, SUSCEPTIBLE, "paper's 3 susceptible controls")
    narrow(clusters, SUSCEPTIBLE + [AILSA], "+ Ailsa Craig (recombinant)")

    # Ailsa Craig is a *valid, non-circular* recombinant because it is phenotyped susceptible.
    # Ailsa Craig's recombination breakpoint within the EB-9 region (it carries a separate
    # donor block far to the right, so restrict to the EB-9 neighborhood).
    shared = compare_to_benchmark(clusters, DONOR).query(
        "sample == @AILSA and shares_benchmark == True and 62_000_000 <= position <= 63_200_000"
    )["position"]
    a, b = pheno.loc[AILSA], pheno.loc[BGV]
    print(
        f"\nAilsa Craig  — {a['mean_stem_lesion_pct']}% stem lesion ({a['class']}); "
        f"carries the Devon Surprise haplotype across the EB-9 region to ~{shared.max() / 1e6:.2f} Mb, then recombines away."
    )
    print(
        "               A susceptible carrying the left flank proves it isn't sufficient -> left edge 62.60 Mb."
    )
    print(
        f"\nBGV007871    — {b['mean_stem_lesion_pct']}% stem lesion ({b['class']}); "
        "carries Ailsa Craig's haplotype yet is resistant."
    )
    print(
        "               Same haplotype, opposite phenotype => resistance from another source; the right edge can't be trimmed."
    )


if __name__ == "__main__":
    main()
