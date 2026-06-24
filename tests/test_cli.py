"""Test the CLI entry point end-to-end via ``python -m haploqtl``."""

import subprocess
import sys

import pandas as pd


def test_cli_cluster_writes_output(fixture_vcf, tmp_path):
    out = tmp_path / "out.csv"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "haploqtl",
            "cluster",
            str(fixture_vcf),
            "--chrom",
            "ch09",
            "--window",
            "250000",
            "--step",
            "1000000",
            "--min-snps",
            "10",
            "--d-min",
            "10",
            "--d-max",
            "40",
            "--d-step",
            "10",
            "-o",
            str(out),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert out.exists(), result.stderr
    df = pd.read_csv(out)
    assert df["sample"].nunique() == 780
    assert {"chromosome", "position", "sample", "cluster", "PC1", "PC2"}.issubset(df.columns)
