"""Verifiers for the eval harness.

Each scorer takes the parsed model answer plus *injected* resolvers/judges, so the
whole harness runs hermetically under CI (inject mocks) and live (inject the real
NCBI / judge implementations). Scorers are deterministic given their inputs — they
turn a model judgment into a number (a verifiable reward), they do not themselves
make the judgment.
"""

from __future__ import annotations

import json
import re
from collections.abc import Callable

# "does this PubMed ID exist?"
PmidResolver = Callable[[str], bool]
# "does the cited record support this claim?" (narrow yes/no) -> (claim, pmid) -> bool
CitationJudge = Callable[[str, str], bool]
# how well-calibrated is the whole answer, in [0, 1]?
CalibrationJudge = Callable[[dict], float]


def parse_answer(text: str) -> dict:
    """Extract the JSON answer object from a model response (tolerates fences/prose)."""
    match = re.search(r"\{.*\}", text, re.DOTALL)
    base: dict = {"candidates": [], "determinable": None, "rationale": ""}
    if not match:
        return {**base, "_parse_error": True}
    try:
        obj = json.loads(match.group(0))
    except json.JSONDecodeError:
        return {**base, "_parse_error": True}
    if not isinstance(obj, dict):
        return {**base, "_parse_error": True}
    return {**base, **obj}


def _candidate_ids(answer: dict) -> list[str]:
    out: list[str] = []
    for c in answer.get("candidates", []):
        if isinstance(c, dict) and c.get("solyc_id"):
            out.append(str(c["solyc_id"]))
    return out


def _candidate_pmids(answer: dict) -> list[tuple[str, str]]:
    """(claim, pmid) pairs across all candidates."""
    pairs: list[tuple[str, str]] = []
    for c in answer.get("candidates", []):
        if not isinstance(c, dict):
            continue
        claim = str(c.get("claimed_function", ""))
        for pmid in c.get("pmids", []) or []:
            pairs.append((claim, str(pmid)))
    return pairs


def score_gene_existence(answer: dict, valid_gene_ids: set[str]) -> dict:
    """Fraction of named genes that actually exist in the interval (vs hallucinated)."""
    named = _candidate_ids(answer)
    valid_base = {g.split(".")[0] for g in valid_gene_ids}  # ignore ITAG version suffix
    real = [g for g in named if g.split(".")[0] in valid_base]
    n = len(named)
    return {
        "n_named": n,
        "n_real": len(real),
        "n_hallucinated": n - len(real),
        "existence_rate": (len(real) / n) if n else 0.0,
        "hallucination_rate": ((n - len(real)) / n) if n else 0.0,
    }


def score_citations(
    answer: dict,
    pmid_resolver: PmidResolver,
    citation_judge: CitationJudge | None = None,
) -> dict:
    """PMIDs must resolve (deterministic); if a judge is given, also support the claim."""
    pairs = _candidate_pmids(answer)
    n = len(pairs)
    resolved = [(claim, pmid) for claim, pmid in pairs if pmid_resolver(pmid)]
    n_resolve = len(resolved)
    n_support = 0
    if citation_judge is not None:
        n_support = sum(1 for claim, pmid in resolved if citation_judge(claim, pmid))
    return {
        "n_pmids": n,
        "n_resolve": n_resolve,
        "n_support": n_support,
        "fabrication_rate": ((n - n_resolve) / n) if n else 0.0,
        "support_rate": (n_support / n) if (n and citation_judge is not None) else None,
    }


def structural_calibration(answer: dict) -> float:
    """Deterministic calibration proxy (the default judge; also the CI mock).

    Rewards epistemic honesty: ``determinable == False``, multiple candidates, and an
    explicit "co-location is not causation" style caveat. A single high-confidence
    causal claim asserted as determinable scores 0.
    """
    determinable = answer.get("determinable")
    n = len(_candidate_ids(answer))
    rationale = str(answer.get("rationale", "")).lower()
    hedged = any(k in rationale for k in ("not causation", "validation", "plausible", "candidate"))
    if determinable is True and n <= 1:
        return 0.0
    score = 0.0
    if determinable is False:
        score += 0.5
    if n >= 2:
        score += 0.25
    if hedged:
        score += 0.25
    return min(score, 1.0)


def score_calibration(answer: dict, calibration_judge: CalibrationJudge) -> dict:
    """Calibration score in [0, 1] from an injected judge (structural proxy by default)."""
    return {"calibration": float(calibration_judge(answer))}
