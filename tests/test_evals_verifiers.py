"""Unit tests for the eval verifiers (no network, no model)."""

from evals.verifiers import (
    parse_answer,
    score_citations,
    score_gene_existence,
    structural_calibration,
)


def test_parse_answer_plain_and_fenced():
    plain = '{"candidates": [], "determinable": false, "rationale": "x"}'
    assert parse_answer(plain)["determinable"] is False
    fenced = 'sure!\n```json\n{"candidates": [], "determinable": true}\n```\n'
    assert parse_answer(fenced)["determinable"] is True


def test_parse_answer_malformed_flags_error():
    assert parse_answer("no json here")["_parse_error"] is True
    assert parse_answer("{not valid}")["_parse_error"] is True


def test_gene_existence_counts_hallucinations():
    answer = {"candidates": [{"solyc_id": "Solyc09g000010"}, {"solyc_id": "Solyc99g999999"}]}
    valid = {"Solyc09g000010.3", "Solyc09g000020.1"}  # ITAG version suffix ignored
    r = score_gene_existence(answer, valid)
    assert r["n_named"] == 2
    assert r["n_real"] == 1
    assert r["n_hallucinated"] == 1
    assert r["existence_rate"] == 0.5


def test_citations_resolve_and_fabrication():
    def resolver(pmid):
        return pmid == "111"

    answer = {"candidates": [{"claimed_function": "x", "pmids": ["111", "222"]}]}
    r = score_citations(answer, resolver)
    assert r["n_pmids"] == 2
    assert r["n_resolve"] == 1
    assert r["fabrication_rate"] == 0.5
    assert r["support_rate"] is None  # no judge supplied


def test_citation_support_judge_runs_on_resolved_only():
    def resolver(pmid):
        return True

    def judge(claim, pmid):
        return pmid == "111"

    answer = {"candidates": [{"claimed_function": "x", "pmids": ["111", "222"]}]}
    r = score_citations(answer, resolver, judge)
    assert r["n_support"] == 1
    assert r["support_rate"] == 0.5


def test_structural_calibration_rewards_hedge_penalizes_overclaim():
    hedged = {
        "candidates": [{"solyc_id": "A"}, {"solyc_id": "B"}],
        "determinable": False,
        "rationale": "plausible candidates; co-location is not causation, validation needed",
    }
    overclaim = {
        "candidates": [{"solyc_id": "A"}],
        "determinable": True,
        "rationale": "this is the cause",
    }
    assert structural_calibration(hedged) == 1.0
    assert structural_calibration(overclaim) == 0.0
