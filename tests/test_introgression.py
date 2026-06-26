"""Introgression interval-calling, reduction, donor-block retention, and SNP refinement.

Pins the paper's EB-9 numbers (recovered from code) and the internal coherence of the
layer: the two-way interval equals the donor-block intersection (the fine-mapped core), and
the SNP refinement falls inside the window-level call.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from haploqtl import contrast
from haploqtl.introgression import (
    _gap_tolerant_runs,
    call_interval,
    donor_block_summary,
    interval_reduction,
    refine_with_markers,
)

FIX = Path(__file__).parent / "fixtures"
VCF = Path(__file__).parent.parent / "data" / "SL4.0ch09_subset.vcf.gz"
RESISTANT = [
    "191163",
    "191164",
    "191167",
    "191172",
    "191175",
    "191174H",
    "CU3AllData",
    "191357",
    "201041",
]
SUSCEPTIBLE = ["191165", "ERR418112", "SRR1572598"]
PRIOR = (61_819_509, 63_679_761)  # prior EB-9 interval (the reduction denominator)
WINDOW = 250_000
DONOR = "191163"  # Devon Surprise


def _clusters() -> pd.DataFrame:
    return pd.read_csv(FIX / "eb9_clusters_subset.csv", dtype={"sample": str})


def _call(clusters: pd.DataFrame):
    two = contrast.contrast_twoway(clusters, RESISTANT, SUSCEPTIBLE)
    return call_interval(two, window=WINDOW, max_gap=1)


def test_eb9_interval_and_70pct_reduction():
    call = _call(_clusters())
    assert call is not None
    assert call.n_windows == 4
    assert (call.start, call.stop) == (62_400_075, 62_950_075)  # window-coverage bounds
    red = interval_reduction((call.start, call.stop), PRIOR)
    assert 0.69 <= red.reduction <= 0.72  # the paper's ~70% narrowing, recovered from code


def test_core_equals_two_way_interval():
    """The intersection of resistant donor-blocks is the same object as the two-way call."""
    clusters = _clusters()
    call = _call(clusters)
    anchor = (call.center_start + call.center_stop) / 2
    summary = donor_block_summary(
        clusters, DONOR, RESISTANT, window=WINDOW, anchor=anchor, max_gap=1
    )
    assert (summary.core_start, summary.core_stop) == (call.start, call.stop)


def test_donor_block_extents_ordering():
    clusters = _clusters()
    call = _call(clusters)
    anchor = (call.center_start + call.center_stop) / 2
    summary = donor_block_summary(
        clusters, DONOR, RESISTANT, window=WINDOW, anchor=anchor, max_gap=1
    )
    by_sample = {e.sample: e for e in summary.per_line}
    # The donor itself retains the largest block; the boundary-defining lines the smallest.
    assert by_sample[DONOR].extent == max(e.extent for e in summary.per_line)
    assert summary.left_boundary_line in RESISTANT and summary.right_boundary_line in RESISTANT
    assert summary.core_span > 0


def test_snp_refine_inside_call():
    call = _call(_clusters())
    refined = refine_with_markers(call, VCF, "ch09", RESISTANT, SUSCEPTIBLE, pad=0)
    assert refined is not None
    assert refined["n_markers"] > 100  # 133 on this set
    assert call.start <= refined["start"] <= refined["stop"] <= call.stop


def test_gap_tolerant_runs_bridge():
    flags = np.array([True, True, False, True, False, False, True])
    assert _gap_tolerant_runs(flags, 0) == [(0, 1), (3, 3), (6, 6)]
    assert _gap_tolerant_runs(flags, 1) == [(0, 3), (6, 6)]  # one-window gap bridged
    assert _gap_tolerant_runs(flags, 2) == [(0, 6)]  # two-window gap bridged


def test_call_interval_none_when_empty():
    track = pd.DataFrame({"position": [1.0, 2.0], "diagnostic": [False, False]})
    assert call_interval(track, window=WINDOW) is None
