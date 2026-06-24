"""Baseline test for the vendored reference pipeline under ``legacy/``.

Guarantees the published reference script still runs end-to-end on the bundled fixture and
emits its original schema, so the modern ``haploqtl`` package can be compared against it.
Fast parameters keep CI quick; this is a contract test, not a scientific reproduction.
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
LEGACY = ROOT / "legacy" / "cluster_haplotypes.py"
LEGACY_COLS = {"PC1", "PC2", "dist", "hclust", "Positions", "Sample", "Chromosome"}


def test_legacy_pipeline_runs_and_schema_is_stable(fixture_vcf, tmp_path):
    # window=250kb, step=1Mb, d-grid {10,20} -> a few windows, runs in seconds.
    cmd = [
        sys.executable,
        str(LEGACY),
        str(fixture_vcf),
        "ch09",
        "250000",
        "1000000",
        "10",
        "10",
        "30",
        "10",
    ]
    subprocess.run(cmd, cwd=tmp_path, check=True, capture_output=True, text=True)

    outputs = list(tmp_path.glob("ch09_window*_pcacomp2.csv"))
    assert len(outputs) == 1, f"expected exactly one output csv, found {outputs}"

    df = pd.read_csv(outputs[0])
    assert LEGACY_COLS.issubset(df.columns)
    assert df["Sample"].nunique() == 780
    assert df["hclust"].notna().any()
