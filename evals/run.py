"""CLI to run the eval live and write a results report.

Hermetic CI does not use this — it drives the harness with mock providers (see
``tests/test_evals_harness_mock.py``). This entry point is for live runs against the
Anthropic API. Example:

    uv sync --extra eval
    export ANTHROPIC_API_KEY=...            # or a gitignored .env
    uv run python -m evals.run \\
        --model claude-opus-4-8 claude-sonnet-4-6 claude-haiku-4-5-20251001 \\
        -o evals/RESULTS.md
"""

from __future__ import annotations

import argparse
import statistics
from pathlib import Path

from evals.harness import load_items, run_all
from evals.live import ncbi_pmid_resolver
from evals.providers import anthropic_provider


def _summary(rows: list[dict]) -> dict:
    ok = [r for r in rows if r["status"] == "ok"]
    nan = float("nan")
    return {
        "n": len(rows),
        "n_ok": len(ok),
        "n_refused": sum(1 for r in rows if r["status"] == "no_response"),
        "n_parse_error": sum(1 for r in rows if r["status"] == "parse_error"),
        "existence": statistics.mean(r["gene_existence"]["existence_rate"] for r in ok)
        if ok
        else nan,
        "fabrication": statistics.mean(r["citations"]["fabrication_rate"] for r in ok)
        if ok
        else nan,
        "calibration": statistics.mean(r["calibration"]["calibration"] for r in ok) if ok else nan,
    }


def _cell(r: dict, group: str, key: str) -> str:
    """A score cell — dash for items the model didn't answer (so refusals aren't read as 0)."""
    return f"{r[group][key]:.2f}" if r["status"] == "ok" else "—"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the candidate-gene eval (live).")
    parser.add_argument(
        "--model", nargs="+", required=True, dest="models", help="one or more Anthropic model ids"
    )
    parser.add_argument(
        "--limit", type=int, default=None, help="run only the first N items (smoke)"
    )
    parser.add_argument("-o", "--output", type=Path, default=None, help="write a markdown report")
    args = parser.parse_args(argv)

    items = load_items()
    if args.limit:
        items = items[: args.limit]

    results: dict[str, list[dict]] = {}
    for model in args.models:
        rows = run_all(anthropic_provider, model, pmid_resolver=ncbi_pmid_resolver, items=items)
        results[model] = rows
        s = _summary(rows)
        print(
            f"\n# {model}  (responded {s['n_ok']}/{s['n']}, "
            f"refused {s['n_refused']}, parse_error {s['n_parse_error']})"
        )
        for r in rows:
            print(
                f"  {r['id']:5} [{r['status']:>12}] genes={r['n_interval_genes']:3} "
                f"existence={_cell(r, 'gene_existence', 'existence_rate')} "
                f"fabrication={_cell(r, 'citations', 'fabrication_rate')} "
                f"calibration={_cell(r, 'calibration', 'calibration')}"
            )
        if s["n_ok"]:
            print(
                f"  MEAN(answered)  existence={s['existence']:.2f} "
                f"fabrication={s['fabrication']:.2f} calibration={s['calibration']:.2f}"
            )

    if args.output:
        lines = [
            "# Eval results — candidate-gene interpretation",
            "",
            "Faithfulness + calibration of the candidate-gene interpretation step (**not** causal",
            "correctness — see [`README.md`](README.md)). Higher *gene existence* = fewer",
            "hallucinated genes; lower *citation fabrication* = fewer non-resolving PMIDs; higher",
            "*calibration* = better hedging of an under-determined call. Means are over items the",
            "model actually answered; refusals and unparseable replies are reported separately, not",
            "scored as faithfulness failures.",
            "",
            "## Summary",
            "",
            "| model | answered | refused | gene existence | citation fabrication | calibration |",
            "|---|---|---|---|---|---|",
        ]
        for model, rows in results.items():
            s = _summary(rows)
            cols = (
                (f"{s['existence']:.2f}", f"{s['fabrication']:.2f}", f"{s['calibration']:.2f}")
                if s["n_ok"]
                else ("—", "—", "—")
            )
            lines.append(
                f"| `{model}` | {s['n_ok']}/{s['n']} | {s['n_refused']} | "
                f"{cols[0]} | {cols[1]} | {cols[2]} |"
            )
        lines += ["", "## Per item", ""]
        for model, rows in results.items():
            lines += [
                f"### `{model}`",
                "",
                "| item | genes | status | existence | fabrication | calibration |",
                "|---|---|---|---|---|---|",
            ]
            lines += [
                f"| {r['id']} | {r['n_interval_genes']} | {r['status']} | "
                f"{_cell(r, 'gene_existence', 'existence_rate')} | "
                f"{_cell(r, 'citations', 'fabrication_rate')} | "
                f"{_cell(r, 'calibration', 'calibration')} |"
                for r in rows
            ]
            lines.append("")
        args.output.write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
