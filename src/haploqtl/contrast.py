"""Contrast haplotype clusters across resistant / susceptible sample sets.

A typed, vectorized port of the three R functions in the original
``visualize_haplotypes.Rmd`` (Anderson et al. 2024) that turned the per-window cluster
table into introgression calls:

* :func:`compare_to_benchmark` <- ``compare_clusters`` (donor "chromosome painting")
* :func:`contrast_oneway`      <- ``contrast_oneway``  (one cluster shared across a set)
* :func:`contrast_twoway`      <- ``contrast_twoway``  (shared by resistant *and* absent
  from susceptible — the diagnostic-introgression pattern)

All operate on the long cluster table emitted by
:func:`haploqtl.cluster.cluster_haplotypes` (columns ``chromosome, position, sample,
cluster, distance_threshold, ...``). Cluster IDs are arbitrary per window, so contrasts
compare *co-membership*, never the ID values themselves.
"""

from __future__ import annotations

import pandas as pd


def _cluster_matrix(clusters: pd.DataFrame) -> pd.DataFrame:
    """Pivot the long cluster table to a (window position x sample) cluster-ID matrix."""
    return clusters.pivot(index="position", columns="sample", values="cluster")


def window_metadata(clusters: pd.DataFrame) -> pd.DataFrame:
    """One row per window: ``chromosome`` and ``distance_threshold``, indexed by position."""
    meta = clusters[["chromosome", "position", "distance_threshold"]].drop_duplicates(
        subset="position"
    )
    return meta.set_index("position").sort_index()


def _all_share(block: pd.DataFrame) -> pd.Series:
    """True per window when every column is non-NA and identical (one shared cluster)."""
    first = block.iloc[:, 0]
    return block.notna().all(axis=1) & block.eq(first, axis=0).all(axis=1)


def compare_to_benchmark(clusters: pd.DataFrame, benchmark: str) -> pd.DataFrame:
    """Per window, flag every sample that shares the benchmark's haplotype cluster.

    Long output ``[position, sample, shares_benchmark]``; ``shares_benchmark`` is NA in
    windows where either the sample or the benchmark was left uncalled (NaN cluster).
    """
    matrix = _cluster_matrix(clusters)
    if benchmark not in matrix.columns:
        raise KeyError(f"benchmark sample {benchmark!r} not in cluster table")
    bench = matrix[benchmark]
    valid = matrix.notna() & bench.notna().to_numpy()[:, None]
    shares = matrix.eq(bench, axis=0).where(valid).copy()  # defragment before melt
    return shares.reset_index().melt(
        id_vars="position", var_name="sample", value_name="shares_benchmark"
    )


def contrast_oneway(clusters: pd.DataFrame, resistant: list[str]) -> pd.DataFrame:
    """Per window, flag where all ``resistant`` samples fall in one shared cluster.

    Output ``[position, diagnostic]`` (bool), in ascending window order.
    """
    block = _cluster_matrix(clusters)[resistant]
    return _all_share(block).rename("diagnostic").reset_index()


def contrast_twoway(
    clusters: pd.DataFrame, resistant: list[str], susceptible: list[str]
) -> pd.DataFrame:
    """Per window, flag the diagnostic-introgression pattern.

    ``diagnostic`` is True where (a) all ``resistant`` samples share one cluster and
    (b) no ``susceptible`` sample is in that cluster — the two-way contrast the paper used
    to call donor introgressions. A susceptible left uncalled (NaN) does not block a call,
    matching the original ``%in%`` semantics. Output ``[position, diagnostic]``.
    """
    matrix = _cluster_matrix(clusters)
    res = matrix[resistant]
    sus = matrix[susceptible]
    res_ok = _all_share(res)
    res_cluster = res.iloc[:, 0]
    sus_clear = sus.ne(res_cluster, axis=0).all(axis=1)  # NA susceptible -> clear (no match)
    return (res_ok & sus_clear).rename("diagnostic").reset_index()
