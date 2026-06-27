"""Chromosome painting from haploqtl.painting (driven by compare_to_benchmark)."""

from pathlib import Path

import pandas as pd

from haploqtl.contrast import contrast_twoway
from haploqtl.painting import build_painting, render_ascii, render_svg

FIX = Path(__file__).parent / "fixtures"
DONOR = "191163"  # Devon Surprise
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
EB9 = (62_400_075, 62_950_075)


def _clusters() -> pd.DataFrame:
    return pd.read_csv(FIX / "eb9_clusters_subset.csv", dtype={"sample": str})


def test_donor_shares_itself_everywhere_called():
    painting = build_painting(_clusters(), DONOR)
    called = [v for v in painting.shares[DONOR] if v is not None]
    assert called and all(called)  # the benchmark always shares its own cluster
    assert painting.shared_fraction[DONOR] >= 0.9


def test_ordered_by_extent_donor_first():
    painting = build_painting(_clusters(), DONOR)
    assert painting.samples[0] == DONOR
    fracs = [painting.shared_fraction[s] for s in painting.samples]
    assert fracs == sorted(fracs, reverse=True)


def test_painting_consistent_with_two_way_contrast():
    """At every diagnostic window, all resistant lines are painted and no susceptible is."""
    clusters = _clusters()
    painting = build_painting(clusters, DONOR)
    two = contrast_twoway(clusters, RESISTANT, SUSCEPTIBLE)
    diagnostic = {round(p) for p in two[two["diagnostic"]]["position"]}
    assert diagnostic
    idx = {round(p): i for i, p in enumerate(painting.positions)}
    for dp in diagnostic:
        i = idx[dp]
        assert all(painting.shares[r][i] for r in RESISTANT)
        assert not any(painting.shares[s][i] for s in SUSCEPTIBLE)


def test_susceptibles_retain_less_than_donor():
    painting = build_painting(_clusters(), DONOR)
    for s in SUSCEPTIBLE:
        assert painting.shared_fraction[s] < painting.shared_fraction[DONOR]


def test_render_ascii_structure():
    painting = build_painting(_clusters(), DONOR)
    text = render_ascii(painting, tags={DONOR: "R"}, eb9=EB9)
    assert "█" in text and "Mb" in text and "EB-9 core" in text


def test_render_svg_wellformed():
    painting = build_painting(_clusters(), DONOR)
    svg = render_svg(painting, eb9=EB9, title="painting test")
    assert svg.startswith("<svg") and svg.rstrip().endswith("</svg>")
    assert "EB-9" in svg and "painting test" in svg
    assert svg.count("<rect") >= len(painting.samples)  # at least a background bar per row
