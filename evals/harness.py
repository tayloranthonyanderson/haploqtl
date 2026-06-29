"""Run the candidate-gene faithfulness + calibration eval over the bundled items.

Tail-interpreter design: precompute the authoritative gene list for each interval
(ITAG4.1 GFF), ask the model to interpret it, and score the answer against those
facts. See ``evals/README.md`` for scope and how to run.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

from evals.providers import CallFn
from evals.verifiers import (
    CalibrationJudge,
    CitationJudge,
    PmidResolver,
    parse_answer,
    score_calibration,
    score_citations,
    score_gene_existence,
    structural_calibration,
)

ROOT = Path(__file__).resolve().parent.parent
GFF = ROOT / "data" / "ITAG4.1_genes.gff3.gz"
ITEMS = Path(__file__).resolve().parent / "items.jsonl"
_SKILL = ROOT / "skills" / "qtl-candidate-gene" / "scripts" / "genes_in_interval.py"


def _load_genes_in_interval() -> Callable[..., list[Any]]:
    """Import ``genes_in_interval`` from the (non-package) skill script by path."""
    name = "evals_genes_in_interval"
    spec = importlib.util.spec_from_file_location(name, _SKILL)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # register so the module's @dataclass can resolve its module
    spec.loader.exec_module(mod)
    return mod.genes_in_interval


def load_items(path: Path = ITEMS) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def interval_genes(item: dict, gff: Path = GFF) -> list[Any]:
    fn = _load_genes_in_interval()
    return fn(gff, item["chrom"], int(item["start"]), int(item["end"]))


def build_prompt(item: dict, genes: list[Any]) -> str:
    table = "\n".join(f"{g.gene_id}\t{g.description}" for g in genes)
    return (
        "You are interpreting a fine-mapped QTL interval in tomato (SL4.0 assembly).\n"
        f"Trait: {item['trait']}\n"
        f"Interval: {item['chrom']}:{item['start']}-{item['end']}\n\n"
        f"Genes in the interval (ITAG4.1; gene_id <tab> description):\n{table}\n\n"
        "Identify the most plausible candidate gene(s) for the trait and support each.\n"
        "Return ONLY a JSON object of the form:\n"
        '{"candidates": [{"solyc_id": "...", "claimed_function": "...", '
        '"pmids": ["..."], "confidence": "high|medium|low"}], '
        '"determinable": <true if the causal gene can be determined from co-location + '
        'annotation alone, else false>, "rationale": "1-3 sentences"}\n'
        "Use only gene IDs from the list above. Cite real PubMed PMIDs. Be honest about "
        "determinability — interval co-location is not proof of causation."
    )


def run_item(
    item: dict,
    call_fn: CallFn,
    model: str,
    *,
    pmid_resolver: PmidResolver,
    calibration_judge: CalibrationJudge = structural_calibration,
    citation_judge: CitationJudge | None = None,
    gff: Path = GFF,
) -> dict:
    genes = interval_genes(item, gff)
    valid_ids = {g.gene_id for g in genes}
    prompt = build_prompt(item, genes)
    raw = call_fn(model, prompt)
    answer = parse_answer(raw)
    if not raw.strip():
        status = "no_response"  # e.g. a model refusal (stop_reason="refusal") -> empty content
    elif answer.get("_parse_error"):
        status = "parse_error"
    else:
        status = "ok"
    return {
        "id": item["id"],
        "trait": item["trait"],
        "n_interval_genes": len(genes),
        "status": status,
        "gene_existence": score_gene_existence(answer, valid_ids),
        "citations": score_citations(answer, pmid_resolver, citation_judge),
        "calibration": score_calibration(answer, calibration_judge),
        "parse_error": bool(answer.get("_parse_error", False)),
    }


def run_all(
    call_fn: CallFn,
    model: str,
    *,
    pmid_resolver: PmidResolver,
    calibration_judge: CalibrationJudge = structural_calibration,
    citation_judge: CitationJudge | None = None,
    items: list[dict] | None = None,
    gff: Path = GFF,
) -> list[dict]:
    rows = items if items is not None else load_items()
    return [
        run_item(
            item,
            call_fn,
            model,
            pmid_resolver=pmid_resolver,
            calibration_judge=calibration_judge,
            citation_judge=citation_judge,
            gff=gff,
        )
        for item in rows
    ]
