"""haploqtl - reproducible, AI-augmented local-ancestry inference for QTL discovery.

The method groups genomes into local haplotypes along a sliding genomic window using
Ward hierarchical clustering, with the merge-distance threshold auto-tuned per window by
maximizing the silhouette coefficient. This traces cryptic ancestral introgressions,
fine-maps QTL, and predicts trait donors across large public sequence libraries.

Reference: Anderson et al. (2024), *The Plant Journal* 117(2):404-415.
https://doi.org/10.1111/tpj.16495
"""

from __future__ import annotations

__version__ = "0.2.0"

from .cluster import WindowClustering, cluster_haplotypes, select_distance_threshold
from .io import GenotypeData, load_genotypes
from .windows import Window, iter_windows

__all__ = [
    "__version__",
    "GenotypeData",
    "load_genotypes",
    "Window",
    "iter_windows",
    "WindowClustering",
    "cluster_haplotypes",
    "select_distance_threshold",
]
