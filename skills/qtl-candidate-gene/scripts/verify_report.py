"""Step 6 — verify a draft candidate-gene report before it is delivered.

A deterministic gate over the model's own draft. It enforces the three things the report
can get *quietly* wrong, using facts the workflow already has:

  * **genes** — every candidate must be a real gene in the interval (Step 1's list).
    A hallucinated Solyc ID is dropped.
  * **citations** — every cited PMID must resolve on PubMed (``pubmed.pmid_exists``).
    A non-resolving PMID is stripped, and the claim is flagged to re-retrieve (Step 3);
    the gate never invents a replacement.
  * **support (optional)** — with ``--support-model``, each *resolving* PMID is also
    checked for whether its abstract actually supports the claim, so a real-but-off-topic
    citation is stripped too, and a supporting quote is recorded for the report.
  * **calibration** — an over-confident, un-hedged single-gene call is flagged.

It emits a **verification stamp** to append to the report, so the deliverable carries
proof it was grounded and self-checked. ``verify`` is pure and takes injected checkers
(a PMID resolver, optionally a support judge), so it tests without a network; the CLI
wires the real PubMed resolver and, with ``--support-model``, the support judge.

CLI:

    python verify_report.py --report draft.json --chrom ch09 --start 62452852 --end 63002852
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Callable
from pathlib import Path

PmidResolver = Callable[[str], bool]
# (claim, pmid) -> {"supported": bool, "evidence": str}; see cite_support.make_support_judge
SupportJudge = Callable[[str, str], dict]

# Words that signal the report hedges an under-determined call (co-location != causation).
_HEDGE = ("not causation", "validation", "plausible", "candidate", "hypoth", "co-loc")


def _calibration_ok(report: dict) -> bool:
    """Did the report hedge appropriately, or over-declare a single causal gene?"""
    det = report.get("determinable")
    cands = [c for c in report.get("candidates", []) if isinstance(c, dict) and c.get("solyc_id")]
    rationale = str(report.get("rationale", "")).lower()
    hedged = any(w in rationale for w in _HEDGE)
    if det is True and len(cands) <= 1:
        return False
    return bool(det is False or len(cands) >= 2 or hedged)


def verify(
    report: dict,
    valid_gene_ids: set[str],
    pmid_resolver: PmidResolver,
    *,
    support_judge: SupportJudge | None = None,
) -> dict:
    """Gate the draft ``report`` against the interval's genes and PubMed.

    Returns ``{stamp, dropped_genes, dead_pmids, unsupported_pmids, calibration_ok, report}``
    where ``report`` is the cleaned draft: hallucinated genes removed, non-resolving PMIDs
    stripped, and — if ``support_judge`` is given — PMIDs that resolve but don't support the
    claim stripped too (with a supporting quote recorded for the ones that do).
    """
    base = {str(g).split(".")[0] for g in valid_gene_ids}  # ignore ITAG version suffix
    candidates = [c for c in report.get("candidates", []) if isinstance(c, dict)]
    kept: list[dict] = []
    dropped_genes: list[str] = []
    dead_pmids: list[str] = []
    unsupported_pmids: list[str] = []
    n_cited = n_resolved = n_supported = 0

    for c in candidates:
        gid = str(c.get("solyc_id", ""))
        if gid.split(".")[0] not in base:
            dropped_genes.append(gid)
            continue
        claim = str(c.get("claimed_function", ""))
        good: list[str] = []
        removed_dead: list[str] = []
        removed_unsupported: list[str] = []
        evidence: dict[str, str] = {}
        for pmid in c.get("pmids", []) or []:
            pmid = str(pmid)
            n_cited += 1
            if not pmid_resolver(pmid):
                removed_dead.append(pmid)
                dead_pmids.append(pmid)
                continue
            n_resolved += 1
            if support_judge is not None:
                verdict = support_judge(claim, pmid)
                if not verdict.get("supported"):
                    removed_unsupported.append(pmid)
                    unsupported_pmids.append(pmid)
                    continue
                n_supported += 1
                if verdict.get("evidence"):
                    evidence[pmid] = str(verdict["evidence"])
            good.append(pmid)
        cleaned_c = {**c, "pmids": good}
        if removed_dead:
            cleaned_c["unresolved_pmids_removed"] = removed_dead
        if removed_unsupported:
            cleaned_c["unsupported_pmids_removed"] = removed_unsupported
        if evidence:
            cleaned_c["support_evidence"] = evidence
        kept.append(cleaned_c)

    calibration_ok = _calibration_ok(report)
    stamp = {
        "genes_in_interval": f"{len(kept)}/{len(candidates)}",
        "pmids_resolved": f"{n_resolved}/{n_cited}",
    }
    if support_judge is not None:
        stamp["pmids_support_claim"] = f"{n_supported}/{n_resolved}"
    stamp["calibration"] = "ok" if calibration_ok else "flagged: over-confident / under-hedged"
    return {
        "stamp": stamp,
        "dropped_genes": dropped_genes,
        "dead_pmids": dead_pmids,
        "unsupported_pmids": unsupported_pmids,
        "calibration_ok": calibration_ok,
        "report": {**report, "candidates": kept},
    }


def render_stamp(result: dict) -> str:
    """A one-block markdown stamp to append to the delivered report."""
    s = result["stamp"]
    head = (
        "**Verification** — "
        f"genes in interval {s['genes_in_interval']} · "
        f"PMIDs resolved {s['pmids_resolved']} · "
    )
    if "pmids_support_claim" in s:
        head += f"PMIDs supporting the claim {s['pmids_support_claim']} · "
    head += f"calibration: {s['calibration']}"
    lines = [head]
    if result["dropped_genes"]:
        lines.append(f"- dropped (not in interval): {', '.join(result['dropped_genes'])}")
    if result["dead_pmids"]:
        lines.append(
            f"- removed (unresolvable PMIDs): {', '.join(result['dead_pmids'])} "
            "— re-run Step 3 retrieval for those claims, do not invent replacements"
        )
    if result.get("unsupported_pmids"):
        lines.append(
            f"- removed (resolve but don't support the claim): "
            f"{', '.join(result['unsupported_pmids'])} — re-retrieve for those claims"
        )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    import genes_in_interval as gii  # same dir; lazy so `verify` imports without it
    import pubmed

    parser = argparse.ArgumentParser(description="Verify a draft candidate-gene report.")
    parser.add_argument(
        "--report", default="-", help="draft report JSON file, or '-' for stdin (default)"
    )
    parser.add_argument("--chrom", required=True, help="Chromosome, e.g. ch09 or SL4.0ch09.")
    parser.add_argument("--start", type=int, required=True, help="Interval start (bp, SL4.0).")
    parser.add_argument("--end", type=int, required=True, help="Interval end (bp, SL4.0).")
    parser.add_argument(
        "--gff", type=Path, default=gii.DEFAULT_GFF, help="ITAG4.1 GFF3 (default: bundled slice)."
    )
    parser.add_argument(
        "--support-model",
        default=None,
        help="check each resolving PMID actually supports its claim with this model "
        "(needs the 'anthropic' package + ANTHROPIC_API_KEY); off by default.",
    )
    args = parser.parse_args(argv)

    raw = sys.stdin.read() if args.report == "-" else Path(args.report).read_text()
    report = json.loads(raw)

    support_judge = None
    if args.support_model:
        import cite_support

        support_judge = cite_support.make_support_judge(
            args.support_model, abstract_fetcher=pubmed.efetch_abstract
        )

    seqid = gii.normalize_chrom(args.chrom)
    valid = {g.gene_id for g in gii.genes_in_interval(args.gff, seqid, args.start, args.end)}
    result = verify(report, valid, pubmed.pmid_exists, support_judge=support_judge)

    print(render_stamp(result))
    print("\n--- verified report (JSON) ---")
    print(json.dumps(result["report"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
