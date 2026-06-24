"""Tests for the diagnostic_variants skill script (offline, fixture VCF)."""

MARKER_COLS = {
    "chromosome",
    "position",
    "ref",
    "alt",
    "resistant_allele",
    "resistant_carrier_fraction",
    "susceptible_carrier_fraction",
}


def test_eb9_diagnostic_markers(load_script, fixture_vcf):
    mod = load_script("diagnostic_variants")
    markers, meta = mod.find_diagnostic_markers(
        fixture_vcf,
        "SL4.0ch09",
        62452852,
        63002852,
        mod.DEFAULT_RESISTANT,
        mod.DEFAULT_SUSCEPTIBLE,
    )
    assert meta["resistant_n"] == 6 and meta["susceptible_n"] == 3
    assert not meta["resistant_missing"] and not meta["susceptible_missing"]
    assert len(markers) > 50, "expected paper-scale diagnostic markers in EB-9 core"
    assert MARKER_COLS.issubset(markers[0])
    # every reported marker is resistant-fixed and susceptible-absent (max_exceptions=0)
    assert all(m["resistant_carrier_fraction"] == 1.0 for m in markers)
    assert all(m["susceptible_carrier_fraction"] == 0.0 for m in markers)


def test_resolve_samples_strips_x_prefix(load_script):
    mod = load_script("diagnostic_variants")
    resolved, missing = mod.resolve_samples(
        ["Devon Surprise"], {"Devon Surprise": "X191163"}, {"191163"}
    )
    assert resolved == ["191163"]
    assert not missing
