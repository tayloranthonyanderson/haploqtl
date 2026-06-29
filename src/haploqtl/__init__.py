"""haploqtl - reproducible, typed local-ancestry inference for QTL discovery.

The method groups genomes into local haplotypes along a sliding genomic window using
Ward hierarchical clustering, with the merge-distance threshold auto-tuned per window by
maximizing the silhouette coefficient. This traces cryptic ancestral introgressions,
fine-maps QTL, and predicts trait donors across large public sequence libraries.

Reference: Anderson et al. (2024), *The Plant Journal* 117(2):404-415.
https://doi.org/10.1111/tpj.16495
"""

from __future__ import annotations

__version__ = "0.9.0"

from .accessions import load_name_map, resolve_samples
from .cluster import WindowClustering, cluster_haplotypes, select_distance_threshold
from .contrast import compare_to_benchmark, contrast_oneway, contrast_twoway
from .introgression import (
    BlockExtent,
    DonorBlockSummary,
    IntervalReduction,
    IntrogressionCall,
    call_interval,
    donor_block_summary,
    interval_reduction,
    refine_with_markers,
)
from .io import GenotypeData, load_genotypes
from .markers import DiagnosticMarker, find_diagnostic_markers
from .painting import Painting, build_painting, render_ascii, render_svg
from .windows import Window, iter_windows

__all__ = [
    "__version__",
    # io / windows / clustering
    "GenotypeData",
    "load_genotypes",
    "Window",
    "iter_windows",
    "WindowClustering",
    "cluster_haplotypes",
    "select_distance_threshold",
    # accessions
    "load_name_map",
    "resolve_samples",
    # contrast
    "compare_to_benchmark",
    "contrast_oneway",
    "contrast_twoway",
    # markers
    "DiagnosticMarker",
    "find_diagnostic_markers",
    # introgression
    "IntrogressionCall",
    "IntervalReduction",
    "BlockExtent",
    "DonorBlockSummary",
    "call_interval",
    "interval_reduction",
    "donor_block_summary",
    "refine_with_markers",
    # painting
    "Painting",
    "build_painting",
    "render_ascii",
    "render_svg",
]
