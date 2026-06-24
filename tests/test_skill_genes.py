"""Tests for the genes_in_interval skill script (offline, bundled ITAG4.1 slice)."""


def test_eb9_recovers_paper_candidate_families(load_script):
    mod = load_script("genes_in_interval")
    genes = mod.genes_in_interval(mod.DEFAULT_GFF, "SL4.0ch09", 62452852, 63002852)
    assert genes, "no genes found in EB-9 interval"
    descriptions = " ".join(g.description.lower() for g in genes)
    for term in ("potassium", "f-box", "cation efflux", "metal tolerance", "oxygenase"):
        assert term in descriptions, f"expected a {term!r} gene in the EB-9 interval"
    first = genes[0]
    assert first.gene_id.startswith("Solyc09g")
    assert first.seqid == "SL4.0ch09"
    assert first.start <= first.end


def test_eb5_includes_resistance_protein(load_script):
    mod = load_script("genes_in_interval")
    genes = mod.genes_in_interval(mod.DEFAULT_GFF, "SL4.0ch05", 62000000, 63500000)
    hits = [g for g in genes if g.gene_id.startswith("Solyc05g053980")]
    assert hits, "Solyc05g053980 not found in EB-5 region"
    assert "resistance" in hits[0].description.lower()


def test_normalize_chrom(load_script):
    mod = load_script("genes_in_interval")
    assert mod.normalize_chrom("9") == "SL4.0ch09"
    assert mod.normalize_chrom("ch9") == "SL4.0ch09"
    assert mod.normalize_chrom("ch09") == "SL4.0ch09"
    assert mod.normalize_chrom("SL4.0ch09") == "SL4.0ch09"
