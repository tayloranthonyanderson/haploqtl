"""Tests for the gene_function skill script (UniProt parsing, hermetic via injected fetch)."""

HEADER = "Entry\tProtein names\tGene Names\tKeywords"


def test_parses_uniprot_tsv(load_script):
    mod = load_script("gene_function")
    tsv = f"{HEADER}\nA0A3Q7J1B7\tPotassium transporter\tLOC1\tPotassium;Transport\n"
    result = mod.gene_function("Solyc09g074790.3", fetch=lambda _q: tsv)
    assert result["uniprot"] == "A0A3Q7J1B7"
    assert result["protein_name"] == "Potassium transporter"
    assert "Potassium" in result["keywords"]
    assert result["gene_id"] == "Solyc09g074790.3"  # original id preserved in output


def test_handles_no_hit(load_script):
    mod = load_script("gene_function")
    result = mod.gene_function("Solyc05g053980", fetch=lambda _q: f"{HEADER}\n")
    assert result["uniprot"] is None
    assert result["protein_name"] is None


def test_strips_version_in_query(load_script):
    mod = load_script("gene_function")
    captured = {}

    def fake(query: str) -> str:
        captured["query"] = query
        return f"{HEADER}\n"

    mod.gene_function("Solyc09g074790.3", fetch=fake)
    queried_id = captured["query"].split(" AND")[0]
    assert queried_id == "Solyc09g074790"  # ITAG version suffix stripped before querying
