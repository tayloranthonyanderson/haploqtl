"""Step 6 — verify a draft candidate-gene report before it is delivered.

A deterministic gate over the model's own draft. It enforces the three things the report
can get *quietly* wrong, using facts the workflow already has:

  * **genes** — every candidate must be a real gene in the interval (Step 1's list).
    A hallucinated Solyc ID is dropped.
  * **citations** — every cited PMID must resolve on PubMed (``pubmed.pmid_exists``).
    A non-resolving PMID is stripped, and the claim is flagged to re-retrieve (Step 3);
    the gate never invents a replacement.
  * **calibration** — an over-confident, un-hedged single-gene call is flagged.

It emits a **verification stamp** to append to the report, so the deliverable carries
proof it was grounded and self-checked. ``verify`` is pure and takes an injected
resolver, so it tests without a network; the CLI wires the real PubMed resolver.

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


def verify(report: dict, valid_gene_ids: set[str], pmid_resolver: PmidResolver) -> dict:
    """Gate the draft ``report`` against the interval's genes and PubMed.

    Returns ``{stamp, dropped_genes, dead_pmids, calibration_ok, report}`` where ``report``
    is the cleaned draft (hallucinated genes removed, non-resolving PMIDs stripped).
    """
    base = {str(g).split(".")[0] for g in valid_gene_ids}  # ignore ITAG version suffix
    candidates = [c for c in report.get("candidates", []) if isinstance(c, dict)]
    kept: list[dict] = []
    dropped_genes: list[str] = []
    dead_pmids: list[str] = []
    n_cited = n_resolved = 0

    for c in candidates:
        gid = str(c.get("solyc_id", ""))
        if gid.split(".")[0] not in base:
            dropped_genes.append(gid)
            continue
        good: list[str] = []
        removed: list[str] = []
        for pmid in c.get("pmids", []) or []:
            pmid = str(pmid)
            n_cited += 1
            if pmid_resolver(pmid):
                n_resolved += 1
                good.append(pmid)
            else:
                removed.append(pmid)
                dead_pmids.append(pmid)
        cleaned_c = {**c, "pmids": good}
        if removed:
            cleaned_c["unresolved_pmids_removed"] = removed
        kept.append(cleaned_c)

    calibration_ok = _calibration_ok(report)
    stamp = {
        "genes_in_interval": f"{len(kept)}/{len(candidates)}",
        "pmids_resolved": f"{n_resolved}/{n_cited}",
        "calibration": "ok" if calibration_ok else "flagged: over-confident / under-hedged",
    }
    return {
        "stamp": stamp,
        "dropped_genes": dropped_genes,
        "dead_pmids": dead_pmids,
        "calibration_ok": calibration_ok,
        "report": {**report, "candidates": kept},
    }


def render_stamp(result: dict) -> str:
    """A one-block markdown stamp to append to the delivered report."""
    s = result["stamp"]
    lines = [
        "**Verification** — "
        f"genes in interval {s['genes_in_interval']} · "
        f"PMIDs resolved {s['pmids_resolved']} · "
        f"calibration: {s['calibration']}"
    ]
    if result["dropped_genes"]:
        lines.append(f"- dropped (not in interval): {', '.join(result['dropped_genes'])}")
    if result["dead_pmids"]:
        lines.append(
            f"- removed (unresolvable PMIDs): {', '.join(result['dead_pmids'])} "
            "— re-run Step 3 retrieval for those claims, do not invent replacements"
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
    args = parser.parse_args(argv)

    raw = sys.stdin.read() if args.report == "-" else Path(args.report).read_text()
    report = json.loads(raw)

    seqid = gii.normalize_chrom(args.chrom)
    valid = {g.gene_id for g in gii.genes_in_interval(args.gff, seqid, args.start, args.end)}
    result = verify(report, valid, pubmed.pmid_exists)

    print(render_stamp(result))
    print("\n--- verified report (JSON) ---")
    print(json.dumps(result["report"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
