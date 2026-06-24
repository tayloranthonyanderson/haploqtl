"""Windowed local-ancestry haplotype clustering.

For each sliding window the genomes are grouped into local haplotypes with Ward
hierarchical agglomerative clustering. The merge-distance threshold is auto-tuned per
window by maximizing the mean silhouette coefficient, so the number of haplotype clusters
emerges from the genetic variance present in that window.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import scale

from .io import GenotypeData
from .windows import iter_windows

log = logging.getLogger(__name__)

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int_]


@dataclass(frozen=True)
class WindowClustering:
    """Result of selecting and applying the best distance threshold in one window."""

    distance_threshold: float
    labels: IntArray
    silhouette: float


def select_distance_threshold(X: FloatArray, d_grid: FloatArray) -> WindowClustering | None:
    """Choose the Ward merge-distance threshold that maximizes the silhouette coefficient.

    Thresholds that collapse the samples into a single cluster, or that split them into one
    cluster per sample, leave the silhouette undefined and are skipped.

    Bug fix vs. the original implementation: the reference script wrapped the *entire*
    threshold search in one ``try/except ValueError``, so a single degenerate threshold
    (typically the small-distance end, which over-splits) aborted the whole search and fell
    back to a fixed ``d = 10``. Here each threshold is evaluated independently, so the search
    returns the genuine silhouette optimum among valid partitions.

    Returns:
        The best :class:`WindowClustering`, or ``None`` if no threshold in ``d_grid``
        yields a valid (2..n-1 cluster) partition.
    """
    n_samples = X.shape[0]
    best: WindowClustering | None = None
    for d in d_grid:
        labels = AgglomerativeClustering(distance_threshold=float(d), n_clusters=None).fit_predict(
            X
        )
        n_clusters = int(np.unique(labels).size)
        if n_clusters < 2 or n_clusters >= n_samples:
            continue
        score = float(silhouette_score(X, labels))
        if best is None or score > best.silhouette:
            best = WindowClustering(distance_threshold=float(d), labels=labels, silhouette=score)
    return best


def _window_features(dosage_window: npt.NDArray[np.int8]) -> FloatArray | None:
    """Build a scaled (samples x SNPs) feature matrix, dropping monomorphic SNPs.

    Monomorphic (zero-variance) sites would standardize to NaN and break clustering; they
    carry no haplotype information and are removed for numerical stability. Returns ``None``
    if no informative SNPs remain.
    """
    features = dosage_window.T.astype(np.float64)  # (samples, snps)
    informative = features.var(axis=0) > 0
    features = features[:, informative]
    if features.shape[1] == 0:
        return None
    return scale(features)


def cluster_haplotypes(
    data: GenotypeData,
    chromosome: str,
    *,
    window: int = 250_000,
    step: int = 100_000,
    min_snps: int = 10,
    d_grid: FloatArray,
    n_components: int = 2,
) -> pd.DataFrame:
    """Cluster genomes into local haplotypes along a sliding window.

    Returns a tidy long-format frame with one row per (window, sample) and columns:
    ``chromosome, position, sample, cluster, distance_threshold, PC1..PCk``. Windows with
    too few informative SNPs (or no valid clustering) emit rows with NaN cluster/PC values
    so the position axis stays complete for downstream plotting.
    """
    if min_snps < n_components:
        raise ValueError("min_snps must be >= n_components")
    pc_cols = [f"PC{i + 1}" for i in range(n_components)]
    columns = ["chromosome", "position", "sample", "cluster", "distance_threshold", *pc_cols]

    windows = list(iter_windows(data.positions, window, step))
    log.info(
        "Clustering %d windows on %s (%d samples, %d variants)",
        len(windows),
        chromosome,
        data.n_samples,
        data.n_variants,
    )

    records: list[dict[str, object]] = []
    n_called = 0
    for w in windows:
        mask = (data.positions >= w.start) & (data.positions <= w.stop)
        if int(mask.sum()) <= min_snps:
            _append_nan_rows(records, data.samples, chromosome, w.center, pc_cols)
            continue
        X = _window_features(data.dosage[mask])
        result = (
            None if X is None or X.shape[1] < n_components else select_distance_threshold(X, d_grid)
        )
        if X is None or result is None:
            _append_nan_rows(records, data.samples, chromosome, w.center, pc_cols)
            continue
        pcs = PCA(n_components=n_components, random_state=0).fit_transform(X)
        for s in range(data.n_samples):
            row: dict[str, object] = {
                "chromosome": chromosome,
                "position": w.center,
                "sample": data.samples[s],
                "cluster": int(result.labels[s]),
                "distance_threshold": result.distance_threshold,
            }
            for j, col in enumerate(pc_cols):
                row[col] = float(pcs[s, j])
            records.append(row)
        n_called += 1
        log.debug(
            "window %d-%d: k=%d, d=%.2f, silhouette=%.3f",
            w.start,
            w.stop,
            int(np.unique(result.labels).size),
            result.distance_threshold,
            result.silhouette,
        )

    log.info("Assigned haplotypes in %d/%d windows", n_called, len(windows))
    return pd.DataFrame.from_records(records, columns=columns)


def _append_nan_rows(
    records: list[dict[str, object]],
    samples: npt.NDArray[np.str_],
    chromosome: str,
    position: float,
    pc_cols: list[str],
) -> None:
    for sample in samples:
        row: dict[str, object] = {
            "chromosome": chromosome,
            "position": position,
            "sample": sample,
            "cluster": np.nan,
            "distance_threshold": np.nan,
        }
        for col in pc_cols:
            row[col] = np.nan
        records.append(row)
