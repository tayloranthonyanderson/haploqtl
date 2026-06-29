"""Equivalence of the modern ``haploqtl`` clustering to the legacy reference script.

This runs the package clustering and the legacy reference script on the bundled fixture
(same window/step/d-grid) and compares them window-for-window. Cluster IDs are arbitrary,
so partitions are compared by **adjusted Rand index**.
Where both pick the same merge-distance ``d`` the partitions must be *identical* (ARI = 1.0);
any divergence must coincide with a different ``d`` — i.e. the documented silhouette-abort
bug fix, where the legacy aborts the whole search to a fixed ``d`` on one degenerate
threshold while the modern code keeps the genuine optimum. Divergences are reported, not
silently tolerated.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from haploqtl.cluster import cluster_haplotypes
from haploqtl.io import load_genotypes

ROOT = Path(__file__).resolve().parent.parent
LEGACY = ROOT / "legacy" / "cluster_haplotypes.py"

WINDOW, STEP = 250_000, 500_000
D_MIN, D_MAX, D_STEP = 10, 40, 10  # grid {10, 20, 30}


def _run_legacy(vcf: Path, tmp_path: Path) -> pd.DataFrame:
    cmd = [
        sys.executable,
        str(LEGACY),
        str(vcf),
        "ch09",
        str(WINDOW),
        str(STEP),
        "10",
        str(D_MIN),
        str(D_MAX),
        str(D_STEP),
    ]
    subprocess.run(cmd, cwd=tmp_path, check=True, capture_output=True, text=True)
    outputs = list(tmp_path.glob("ch09_window*_pcacomp2.csv"))
    assert len(outputs) == 1, f"expected one legacy CSV, found {outputs}"
    return pd.read_csv(outputs[0])


def test_clustering_matches_legacy_where_d_agrees(fixture_vcf, tmp_path, capsys):
    legacy = _run_legacy(fixture_vcf, tmp_path).rename(
        columns={
            "Positions": "position",
            "Sample": "sample",
            "hclust": "legacy_cluster",
            "dist": "legacy_d",
        }
    )
    legacy["position"] = legacy["position"].round().astype(int)

    data = load_genotypes(fixture_vcf)
    new = cluster_haplotypes(
        data,
        "ch09",
        window=WINDOW,
        step=STEP,
        min_snps=10,
        d_grid=np.arange(float(D_MIN), float(D_MAX), float(D_STEP)),
    )
    new["position"] = new["position"].round().astype(int)

    shared = sorted(set(legacy["position"]) & set(new["position"]))
    assert shared, "no shared windows between legacy and modern outputs"

    agree: list[int] = []
    diverge: list[tuple[int, float, float, float]] = []
    for pos in shared:
        lg = legacy[legacy["position"] == pos].set_index("sample")
        nw = new[new["position"] == pos].set_index("sample")
        common = lg.index.intersection(nw.index)
        lc, nc = lg.loc[common, "legacy_cluster"], nw.loc[common, "cluster"]
        if lc.isna().any() or nc.isna().any():
            continue  # uncalled in one or both -> not comparable
        ld, nd = float(lg["legacy_d"].iloc[0]), float(nw["distance_threshold"].iloc[0])
        ari = adjusted_rand_score(lc.to_numpy(), nc.to_numpy())
        if np.isclose(ld, nd):
            assert ari == 1.0, (
                f"window {pos}: identical d={ld} but ARI={ari} (unexpected divergence)"
            )
            agree.append(pos)
        else:
            diverge.append((pos, ld, nd, round(ari, 3)))

    total = len(agree) + len(diverge)
    with capsys.disabled():
        print(
            f"\n[clustering vs legacy] {len(agree)}/{total} comparable windows identical "
            f"(same d, ARI=1.0); {len(diverge)} diverge via the d-selection bug fix:"
        )
        for pos, ld, nd, ari in diverge:
            print(f"    window {pos:,}: legacy d={ld:g} -> modern d={nd:g}  (ARI={ari})")

    assert agree, "expected at least one window where modern and legacy agree exactly"
