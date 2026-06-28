"""Generalization check: fine-map EB-5 — a *different* QTL, donor, and trait than EB-9.

EB-9 (the other worked example) is collar-rot/stem resistance from *Devon Surprise*. EB-5 is
**foliar** early blight (*Alternaria linariae*) resistance traced to the mid-century breeding
line **Hawaii 7998** (the hypothesized *S. pimpinellifolium* donor was rejected in the paper).
Nothing about the engine is specialized to EB-9: the same clustering + two-way contrast, pointed
at chromosome 5 with Hawaii 7998's pathway and a different susceptible set, independently
recovers the published EB-5 introgression and its ~56% narrowing of the prior QTL.

From the bundled 780-genome chr05 fixture this reproduces:

  1. the EB-5 introgression common to the Hawaii 7998 pathway (HA7998 -> OH7536, OH08-7663) and
     absent from the susceptible controls (OH88119, NC 84173, Brandywine);
  2. a 56.2% reduction of the prior QTL interval (Anderson et al. 2024, Table S2): the refined
     span matches the paper to within ~0.1 kb, though the recovered interval sits ~150 kb right of
     the paper's exact boundaries (threshold sensitivity at the right edge) -- the same ~500 kb
     core and the same reduction.

Two points worth noting:

  * **In the paper the introgression is larger than the QTL support interval.** Its detected
    haplotype (62.35-63.20 Mb) runs ~350 kb past the prior QTL's lower bound, so the *refined QTL*
    is the haplotype intersected with the prior -- which is why Table S2 reports a 500 kb refined
    span from an 850 kb haplotype. ``interval_reduction`` measures the reduction on that overlap;
    EB-9, by contrast, sat nested inside its prior. (Our recovered haplotype happens to land nested
    in the prior, so the intersection is a no-op here -- but the logic is what reconciles the
    paper's own numbers.)
  * **Resistance is tissue-specific.** OH08-7663 carries Hawaii 7998's *foliar* (EB-5) resistance
    yet is *collar-rot* (EB-9) susceptible -- 64.7% stem lesion in the EB-9 phenotype table. The
    EB-5 mapping cross (CU151095-146 x OH08-7663) deliberately paired the two complementary
    resistances, which is why the EB-5 and EB-9 donor sets barely overlap.

Stated candidate gene Solyc05g053980 (63,362,638-63,363,006) sits just past the right shoulder of
the narrowed core -- inside the broader introgression / prior QTL, as it is for the paper's own
250 kb interval.

Run:  uv run python examples/finemap_eb5.py   (~70 s: clusters 780 genomes over chr05)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from haploqtl.cluster import cluster_haplotypes
from haploqtl.contrast import compare_to_benchmark, contrast_twoway
from haploqtl.introgression import call_interval, interval_reduction
from haploqtl.io import load_genotypes

ROOT = Path(__file__).resolve().parent.parent
VCF = ROOT / "data" / "SL4.0ch05_subset.vcf.gz"
PHENO = ROOT / "data" / "eb9_phenotypes.csv"  # EB-9 stem-lesion scores, for the tissue contrast

# The EB-5 / Hawaii 7998 foliar-resistance pathway and the documented susceptible controls,
# per Anderson et al. (2024). Resistant lines are human-named in the VCF; susceptibles are SRA.
DONOR = "HA7998"  # Hawaii 7998
RESISTANT = ["HA7998", "OH7536", "OH7663"]  # Hawaii 7998 -> OH7536 -> OH08-7663
SUSCEPTIBLE = ["OH88119", "SRR1572598", "ERR418112"]  # OH88119, NC 84173, Brandywine
PRIOR = (62_700_265, 63_842_577)  # prior EB-5 QTL interval (Anderson et al. 2024, Table S2)
WINDOW, STEP = 250_000, 50_000


def main() -> None:
    data = load_genotypes(VCF)
    print(f"Clustering {data.n_samples} genomes over chr05 (250 kb / 50 kb window)...\n")
    clusters = cluster_haplotypes(
        data, "ch05", window=WINDOW, step=STEP, min_snps=10, d_grid=np.arange(2.0, 100.0, 8.0)
    )

    track = contrast_twoway(clusters, RESISTANT, SUSCEPTIBLE)
    call = call_interval(track, window=WINDOW, max_gap=1)
    assert call is not None, "no diagnostic EB-5 interval"
    red = interval_reduction((call.start, call.stop), PRIOR)

    print("EB-5 foliar resistance (Hawaii 7998 pathway vs. susceptible controls):")
    print(
        f"  introgression haplotype : {call.start:>12,.0f} - {call.stop:<12,.0f} "
        f"({call.span:,.0f} bp, {call.n_windows} windows)"
    )
    print(
        f"  prior QTL interval      : {PRIOR[0]:>12,.0f} - {PRIOR[1]:<12,.0f} "
        f"({red.prior_span:,.0f} bp)"
    )
    print(
        f"  refined QTL (∩ prior)   : {red.refined[0]:>12,.0f} - {red.refined[1]:<12,.0f} "
        f"({red.refined_span:,.0f} bp, {red.reduction * 100:.1f}% reduction)"
    )
    print("  paper (Anderson et al. 2024, Table S2): 500,126 bp refined, 56.2% reduction\n")

    # Tissue-specific resistance: the EB-5 donor OH08-7663 is collar-rot (EB-9) susceptible.
    pheno = pd.read_csv(PHENO).set_index("vcf_id")
    oh = pheno.loc["OH7663"]
    print(
        f"OH08-7663 — EB-5 foliar-resistant (above), but EB-9 collar-rot SUSCEPTIBLE: "
        f"{oh['mean_stem_lesion_pct']}% stem lesion ({oh['class']})."
    )
    print(
        "           Foliar (EB-5) and stem (EB-9) resistance are independent traits and donors.\n"
    )

    # Other public lines carrying the Hawaii 7998 haplotype across the core (predicted EB-5).
    core = compare_to_benchmark(clusters, DONOR).query(
        "shares_benchmark == True and @call.start <= position <= @call.stop"
    )
    carriers = core.groupby("sample").size()
    extra = sorted(s for s in carriers.index if s not in RESISTANT + SUSCEPTIBLE)
    print(
        f"Other lines carrying the Hawaii 7998 haplotype across the core: {len(extra)} "
        f"(e.g. {', '.join(extra[:4])}) — EB-5 remains rare in public germplasm."
    )


if __name__ == "__main__":
    main()
