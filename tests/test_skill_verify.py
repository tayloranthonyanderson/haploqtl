"""Hermetic tests for the skill's Step 6 verify gate (no network, no API key).

Loads the standalone ``verify_report.py`` by path and exercises the pure ``verify``
function with an injected PMID resolver — the same injection pattern the eval uses.
"""

import importlib.util
import sys
from pathlib import Path

_VR = (
    Path(__file__).resolve().parent.parent
    / "skills"
    / "qtl-candidate-gene"
    / "scripts"
    / "verify_report.py"
)


def _load_verify():
    spec = importlib.util.spec_from_file_location("skill_verify_report", _VR)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["skill_verify_report"] = mod
    spec.loader.exec_module(mod)
    return mod


verify_report = _load_verify()


def _resolver(real):
    return lambda pmid: pmid in real


def test_drops_hallucinated_gene_and_dead_pmid():
    report = {
        "candidates": [
            {"solyc_id": "Solyc09g074790", "claimed_function": "x", "pmids": ["111", "222"]},
            {"solyc_id": "Solyc99g999999", "claimed_function": "fake", "pmids": ["333"]},
        ],
        "determinable": False,
        "rationale": "plausible candidates; co-location is not causation, validation needed",
    }
    valid = {"Solyc09g074790.1", "Solyc09g074510.2"}  # ITAG version suffix ignored
    result = verify_report.verify(report, valid, _resolver({"111"}))

    assert result["dropped_genes"] == ["Solyc99g999999"]
    kept = result["report"]["candidates"]
    assert [c["solyc_id"] for c in kept] == ["Solyc09g074790"]
    assert kept[0]["pmids"] == ["111"]  # dead PMID stripped
    assert kept[0]["unresolved_pmids_removed"] == ["222"]
    assert "222" in result["dead_pmids"]
    assert result["stamp"]["genes_in_interval"] == "1/2"
    assert result["stamp"]["pmids_resolved"] == "1/2"  # counted over the kept candidate only
    assert result["calibration_ok"] is True


def test_flags_overconfident_call():
    report = {
        "candidates": [{"solyc_id": "Solyc09g074790", "pmids": []}],
        "determinable": True,
        "rationale": "this is the cause",
    }
    result = verify_report.verify(report, {"Solyc09g074790"}, _resolver(set()))
    assert result["calibration_ok"] is False
    assert "flagged" in result["stamp"]["calibration"]


def test_render_stamp_reports_removals():
    report = {
        "candidates": [{"solyc_id": "Solyc09g074790", "pmids": ["999"]}],
        "determinable": False,
        "rationale": "plausible",
    }
    result = verify_report.verify(report, {"Solyc09g074790"}, _resolver(set()))
    stamp = verify_report.render_stamp(result)
    assert "PMIDs resolved 0/1" in stamp
    assert "999" in stamp  # the stripped PMID is named for re-retrieval


def test_strips_resolving_but_unsupported_pmid():
    report = {
        "candidates": [
            {
                "solyc_id": "Solyc09g074790",
                "claimed_function": "F-box defense",
                "pmids": ["111", "222"],
            },
        ],
        "determinable": False,
        "rationale": "plausible candidate",
    }

    def support(claim, pmid):  # both PMIDs resolve; only 111 actually supports the claim
        return {"supported": pmid == "111", "evidence": "a quote" if pmid == "111" else ""}

    result = verify_report.verify(
        report, {"Solyc09g074790"}, _resolver({"111", "222"}), support_judge=support
    )
    kept = result["report"]["candidates"][0]
    assert kept["pmids"] == ["111"]  # off-topic-but-real PMID stripped
    assert kept["unsupported_pmids_removed"] == ["222"]
    assert kept["support_evidence"] == {"111": "a quote"}
    assert "222" in result["unsupported_pmids"]
    assert result["stamp"]["pmids_resolved"] == "2/2"
    assert result["stamp"]["pmids_support_claim"] == "1/2"
    stamp = verify_report.render_stamp(result)
    assert "PMIDs supporting the claim 1/2" in stamp
    assert "222" in stamp


def test_no_support_judge_leaves_stamp_without_support_line():
    report = {
        "candidates": [{"solyc_id": "Solyc09g074790", "pmids": ["111"]}],
        "determinable": False,
        "rationale": "plausible",
    }
    result = verify_report.verify(report, {"Solyc09g074790"}, _resolver({"111"}))
    assert "pmids_support_claim" not in result["stamp"]  # deterministic path unchanged
    assert "supporting the claim" not in verify_report.render_stamp(result)
