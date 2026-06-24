"""Tests for the windowed clustering core, including the silhouette-search bug fix."""

import numpy as np

from haploqtl.cluster import cluster_haplotypes, select_distance_threshold
from haploqtl.io import load_genotypes


def test_select_skips_degenerate_thresholds():
    # Two well-separated blobs. A tiny threshold puts every sample in its own cluster, and
    # a huge threshold merges everything into one; both leave the silhouette undefined. The
    # original code aborted the whole search on the first such threshold and fell back to a
    # fixed d=10. The fixed version skips degenerate thresholds and finds the real optimum.
    rng = np.random.default_rng(0)
    X = np.vstack([rng.normal(0.0, 0.1, (20, 5)), rng.normal(10.0, 0.1, (20, 5))])
    d_grid = np.array([0.001, 5.0, 500.0])  # only 5.0 yields a valid 2-cluster partition

    result = select_distance_threshold(X, d_grid)

    assert result is not None
    assert np.unique(result.labels).size == 2
    assert result.distance_threshold == 5.0


def test_select_returns_none_when_all_degenerate():
    X = np.zeros((10, 4))  # identical samples -> always a single cluster
    assert select_distance_threshold(X, np.array([1.0, 5.0])) is None


def test_cluster_haplotypes_on_fixture(fixture_vcf):
    data = load_genotypes(fixture_vcf)
    df = cluster_haplotypes(
        data,
        "ch09",
        window=250_000,
        step=1_000_000,
        min_snps=10,
        d_grid=np.arange(10.0, 40.0, 10.0),
    )
    expected = {"chromosome", "position", "sample", "cluster", "distance_threshold", "PC1", "PC2"}
    assert expected.issubset(df.columns)
    assert df["sample"].nunique() == 780
    assert df["cluster"].dropna().nunique() >= 2, "expected multiple haplotypes in some window"
