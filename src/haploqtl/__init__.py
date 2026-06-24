"""haploqtl — reproducible, AI-augmented local-ancestry inference for QTL discovery.

The method groups genomes into local haplotypes along a sliding genomic window using
Ward hierarchical clustering, with the merge-distance threshold auto-tuned per window by
maximizing the silhouette coefficient. This traces cryptic ancestral introgressions,
fine-maps QTL, and predicts trait donors across large public sequence libraries.

Status: **Phase 0 (foundations).** The clean, typed, tested package API lands in Phase 1.
Until then the published reference implementation is vendored under ``legacy/`` and
exercised by the demo (``scripts/run_demo.sh``) and the test suite. See the README roadmap.

Reference: Anderson et al. (2024), *The Plant Journal* 117(2):404-415.
https://doi.org/10.1111/tpj.16495
"""

__version__ = "0.1.0"
