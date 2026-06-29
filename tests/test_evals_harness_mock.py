"""Hermetic end-to-end test: the oracle provider scores high, the garbage provider low.

Exercises the full pipeline (gzip GFF parse -> prompt -> provider -> verifiers) with no
network and no API key.
"""

from evals.harness import load_items, run_item
from evals.providers import garbage_provider, oracle_provider


def _eb9():
    return next(i for i in load_items() if i["id"] == "eb9")


def _accept_oracle_pmid(pmid):
    return pmid == "10000001"


def test_oracle_scores_high():
    r = run_item(_eb9(), oracle_provider, "mock", pmid_resolver=_accept_oracle_pmid)
    assert r["n_interval_genes"] > 10  # confirms the gzip GFF was parsed
    assert r["gene_existence"]["existence_rate"] == 1.0
    assert r["gene_existence"]["n_named"] >= 1
    assert r["citations"]["fabrication_rate"] == 0.0
    assert r["calibration"]["calibration"] >= 0.75
    assert r["parse_error"] is False


def test_garbage_scores_low():
    r = run_item(_eb9(), garbage_provider, "mock", pmid_resolver=_accept_oracle_pmid)
    assert r["status"] == "ok"  # a real-but-bad answer, not a refusal
    assert r["gene_existence"]["existence_rate"] == 0.0
    assert r["gene_existence"]["n_hallucinated"] == 1
    assert r["citations"]["fabrication_rate"] == 1.0
    assert r["calibration"]["calibration"] == 0.0


def test_empty_response_flagged_as_no_response():
    def refuse(model, prompt):
        return ""  # a refusal returns empty content

    r = run_item(_eb9(), refuse, "mock", pmid_resolver=_accept_oracle_pmid)
    assert r["status"] == "no_response"  # distinct from a faithfulness failure
