"""Smoke test for the bundled reference pipeline.

This proves the local-ancestry pipeline executes end-to-end on the committed chr09
fixture and emits a haplotype-cluster table with a stable schema. It uses fast
parameters (large step + tiny distance grid) so CI stays quick; it is a contract test,
not a scientific reproduction of the paper's fine-mapped intervals.
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
FIXTURE = ROOT / "data" / "SL4.0ch09_subset.vcf.gz"
LEGACY = ROOT / "legacy" / "cluster_haplotypes.py"
EXPECTED_COLS = {"PC1", "PC2", "dist", "hclust", "Positions", "Sample", "Chromosome"}
N_SAMPLES = 780


def test_fixture_present():
    assert FIXTURE.exists(), "bundled chr09 fixture is missing"


def test_pipeline_runs_and_schema_is_stable(tmp_path):
    # window=250kb, step=1Mb, d-grid {10,20} -> 3 windows, runs in a few seconds.
    cmd = [
        sys.executable,
        str(LEGACY),
        str(FIXTURE),
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
    missing = EXPECTED_COLS - set(df.columns)
    assert not missing, f"output is missing expected columns: {missing}"
    assert df["Sample"].nunique() == N_SAMPLES
    assert df["hclust"].notna().any(), "no haplotype clusters were assigned"
