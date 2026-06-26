"""Call introgression intervals and donor-block retention from the cluster contrast.

This is the algorithmic layer the original analysis lacked — it narrowed intervals *by
eye*. Here:

* :func:`call_interval` — the narrowed interval is the longest (gap-tolerant) run of
  windows passing the two-way diagnostic contrast.
* :func:`donor_block_summary` — per-line donor-block extents and their **intersection**
  (the fine-mapped core), plus which line defines each boundary. No generational/decay
  model is implied: these are observed extents, not a rate over time.
* :func:`refine_with_markers` — sharpen a window-level interval to SNP resolution using
  diagnostic markers (the layered cluster -> SNP step).
* :func:`interval_reduction` — % reduction of a called interval vs. a prior interval.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pandas as pd

from .accessions import BUNDLED_RENAME
from .contrast import compare_to_benchmark
from .markers import find_diagnostic_markers

FloatArray = npt.NDArray[np.float64]
BoolArray = npt.NDArray[np.bool_]


@dataclass(frozen=True)
class IntrogressionCall:
    """A narrowed interval from a contiguous (gap-tolerant) run of diagnostic windows.

    ``start``/``stop`` are **window-coverage** bounds: the first/last diagnostic window
    centers expanded by half the window width on each side, since a diagnostic window of
    width ``W`` supports the introgression across its whole ``[center +/- W/2]`` span.
    ``center_start``/``center_stop`` keep the raw window-center run for transparency.
    """

    start: float
    stop: float
    center_start: float
    center_stop: float
    n_windows: int
    max_gap: int

    @property
    def span(self) -> float:
        return self.stop - self.start


@dataclass(frozen=True)
class IntervalReduction:
    called: tuple[float, float]
    prior: tuple[float, float]
    called_span: float
    prior_span: float
    reduction: float  # fraction in [0, 1]


@dataclass(frozen=True)
class BlockExtent:
    """One line's retained donor-block: the diagnostic run overlapping the anchor."""

    sample: str
    start: float
    stop: float
    extent: float
    n_windows: int


@dataclass(frozen=True)
class DonorBlockSummary:
    benchmark: str
    per_line: list[BlockExtent]
    core_start: float
    core_stop: float
    left_boundary_line: str
    right_boundary_line: str

    @property
    def core_span(self) -> float:
        return self.core_stop - self.core_start


def _gap_tolerant_runs(flags: BoolArray, max_gap: int) -> list[tuple[int, int]]:
    """Group True positions into runs, bridging up to ``max_gap`` False windows.

    Returns a list of ``(first_true_index, last_true_index)`` pairs (array-index space).
    """
    idx = np.where(flags)[0]
    if idx.size == 0:
        return []
    runs: list[tuple[int, int]] = []
    start = prev = int(idx[0])
    for k in idx[1:]:
        k = int(k)
        if k - prev <= max_gap + 1:
            prev = k
        else:
            runs.append((start, prev))
            start = prev = k
    runs.append((start, prev))
    return runs


def call_interval(
    track: pd.DataFrame, *, window: float, max_gap: int = 0
) -> IntrogressionCall | None:
    """Narrow to the longest gap-tolerant run of diagnostic windows.

    ``track`` is the output of :func:`haploqtl.contrast.contrast_twoway` (or oneway):
    columns ``[position, diagnostic]``. ``window`` is the clustering window width (bp),
    used to expand the center-run to genomic-coverage bounds. The chosen run maximizes
    genomic span, then the number of diagnostic windows. Returns ``None`` if no window is
    diagnostic.
    """
    track = track.sort_values("position")
    positions = track["position"].to_numpy(dtype=np.float64)
    flags = track["diagnostic"].fillna(False).to_numpy(dtype=bool)
    runs = _gap_tolerant_runs(flags, max_gap)
    if not runs:
        return None
    i0, i1 = max(
        runs, key=lambda r: (positions[r[1]] - positions[r[0]], int(flags[r[0] : r[1] + 1].sum()))
    )
    half = window / 2.0
    return IntrogressionCall(
        start=float(positions[i0] - half),
        stop=float(positions[i1] + half),
        center_start=float(positions[i0]),
        center_stop=float(positions[i1]),
        n_windows=int(flags[i0 : i1 + 1].sum()),
        max_gap=int(max_gap),
    )


def interval_reduction(
    called: tuple[float, float], prior: tuple[float, float]
) -> IntervalReduction:
    """Fractional reduction in interval span of ``called`` relative to ``prior``."""
    called_span = float(called[1] - called[0])
    prior_span = float(prior[1] - prior[0])
    reduction = float(1.0 - called_span / prior_span) if prior_span else 0.0
    return IntervalReduction(
        called=(float(called[0]), float(called[1])),
        prior=(float(prior[0]), float(prior[1])),
        called_span=called_span,
        prior_span=prior_span,
        reduction=reduction,
    )


def _run_containing(
    positions: FloatArray, flags: BoolArray, anchor: float | None, max_gap: int
) -> tuple[int, int] | None:
    """The gap-tolerant run overlapping ``anchor`` (bp); else the longest-span run."""
    runs = _gap_tolerant_runs(flags, max_gap)
    if not runs:
        return None
    if anchor is not None:
        for i0, i1 in runs:
            if positions[i0] <= anchor <= positions[i1]:
                return i0, i1
    return max(runs, key=lambda r: positions[r[1]] - positions[r[0]])


def donor_block_summary(
    clusters: pd.DataFrame,
    benchmark: str,
    resistant: list[str],
    *,
    window: float,
    anchor: float | None = None,
    max_gap: int = 0,
) -> DonorBlockSummary:
    """Per-line retained donor-block extents and their intersection (fine-mapped core).

    For each line in ``resistant``, the retained block is the diagnostic run (sharing the
    ``benchmark`` haplotype cluster) overlapping ``anchor``, reported in window-coverage
    bounds (centers expanded by ``window``/2). The **core** is the intersection
    ``[max(starts), min(stops)]`` over lines whose block overlaps the anchor — the region
    every resistant line still carries from the donor, and (by construction) the same
    object as the two-way :func:`call_interval`.
    """
    wide = (
        compare_to_benchmark(clusters, benchmark)
        .pivot(index="position", columns="sample", values="shares_benchmark")
        .sort_index()
    )
    positions = wide.index.to_numpy(dtype=np.float64)
    half = window / 2.0
    extents: list[BlockExtent] = []
    for sample in resistant:
        if sample not in wide.columns:
            extents.append(BlockExtent(sample, float("nan"), float("nan"), 0.0, 0))
            continue
        flags = wide[sample].fillna(False).to_numpy(dtype=bool)
        run = _run_containing(positions, flags, anchor, max_gap)
        if run is None:
            extents.append(BlockExtent(sample, float("nan"), float("nan"), 0.0, 0))
            continue
        i0, i1 = run
        start, stop = float(positions[i0] - half), float(positions[i1] + half)
        extents.append(
            BlockExtent(
                sample=sample,
                start=start,
                stop=stop,
                extent=stop - start,
                n_windows=int(flags[i0 : i1 + 1].sum()),
            )
        )

    overlapping = [e for e in extents if e.n_windows > 0]
    if overlapping:
        core_start = max(e.start for e in overlapping)
        core_stop = min(e.stop for e in overlapping)
        left = max(overlapping, key=lambda e: e.start).sample
        right = min(overlapping, key=lambda e: e.stop).sample
    else:
        core_start = core_stop = float("nan")
        left = right = ""

    return DonorBlockSummary(
        benchmark=benchmark,
        per_line=extents,
        core_start=core_start,
        core_stop=core_stop,
        left_boundary_line=left,
        right_boundary_line=right,
    )


def refine_with_markers(
    call: IntrogressionCall,
    vcf_path: str | Path,
    seqid: str,
    resistant: list[str],
    susceptible: list[str],
    *,
    pad: int = 0,
    rename_map_path: str | Path = BUNDLED_RENAME,
    max_exceptions: int = 0,
) -> dict | None:
    """Sharpen a window-level call to SNP resolution via diagnostic markers.

    Searches ``[call.start - pad, call.stop + pad]`` for diagnostic markers and returns the
    span of their positions. Returns ``None`` if no diagnostic marker is found.
    """
    start = int(call.start - pad)
    end = int(call.stop + pad)
    markers, meta = find_diagnostic_markers(
        vcf_path,
        seqid,
        start,
        end,
        resistant,
        susceptible,
        rename_map_path=rename_map_path,
        max_exceptions=max_exceptions,
    )
    if not markers:
        return None
    positions = [m.position for m in markers]
    return {
        "start": min(positions),
        "stop": max(positions),
        "n_markers": len(markers),
        "search_region": (start, end),
        "meta": meta,
    }
