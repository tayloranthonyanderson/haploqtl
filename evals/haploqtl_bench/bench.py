#!/usr/bin/env python3
"""Run haploqtl-bench against one or more models, score deterministically, and report.

The model call is injected (`call_fn`), so the harness is hermetically testable with a fake
client and has no hard dependency on a provider SDK. The default runner uses the Anthropic
SDK (imported lazily) and reads ``ANTHROPIC_API_KEY`` from the environment.

Run:  uv run python evals/haploqtl_bench/bench.py            # all default models
      uv run python evals/haploqtl_bench/bench.py --limit 5  # quick smoke
"""

from __future__ import annotations

import argparse
import json
import re
from collections.abc import Callable
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATASET = HERE / "dataset.jsonl"
RESULTS = HERE / "results"
DEFAULT_MODELS = ["claude-opus-4-8", "claude-sonnet-4-6", "claude-haiku-4-5-20251001"]

CallFn = Callable[[str, str], str]  # (model, prompt) -> response text


def load_dataset(path: Path = DATASET) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def build_prompt(item: dict) -> str:
    task = item["task"]
    if task == "marker_diagnosticity":
        res = item["input"]["resistant_alt_dosage"]
        sus = item["input"]["susceptible_alt_dosage"]
        return (
            "You are scoring a SNP as a marker-assisted-selection marker.\n"
            "ALT-allele dosage (0, 1, or 2 copies) for each line:\n"
            f"  resistant lines:   {res}\n"
            f"  susceptible lines: {sus}\n"
            "A SNP is 'diagnostic' if the ALT allele is present (dosage >= 1) in EVERY "
            "resistant line AND absent (dosage == 0) in EVERY susceptible line.\n"
            'Reply with JSON only: {"answer": "yes"} or {"answer": "no"}.'
        )
    if task == "min_interval":
        lines = "\n".join(
            f"  {name}: {lo}-{hi}" for name, (lo, hi) in item["input"]["intervals"].items()
        )
        return (
            "Each line below carries an introgression spanning a base-pair interval on one "
            "chromosome:\n"
            f"{lines}\n"
            "Report the interval shared by ALL lines (the intersection).\n"
            'Reply with JSON only: {"start": <int>, "end": <int>}.'
        )
    raise ValueError(f"unknown task: {task}")


def parse_answer(task: str, text: str) -> object:
    obj = None
    match = re.search(r"\{.*\}", text, re.S)
    if match:
        try:
            obj = json.loads(match.group(0))
        except json.JSONDecodeError:
            obj = None
    if task == "marker_diagnosticity":
        if isinstance(obj, dict) and "answer" in obj:
            return str(obj["answer"]).strip().lower()
        low = text.strip().lower()
        if "yes" in low and "no" not in low:
            return "yes"
        if "no" in low and "yes" not in low:
            return "no"
        return low[:3]
    if task == "min_interval":
        if isinstance(obj, dict) and "start" in obj and "end" in obj:
            try:
                return [int(obj["start"]), int(obj["end"])]
            except (TypeError, ValueError):
                pass
        nums = re.findall(r"\d+", text.replace(",", ""))
        if len(nums) >= 2:
            return [int(nums[0]), int(nums[1])]
        return None
    raise ValueError(f"unknown task: {task}")


def score(item: dict, parsed: object) -> bool:
    return parsed == item["answer"]


def run(models: list[str], call_fn: CallFn, dataset: list[dict] | None = None) -> dict:
    items = dataset if dataset is not None else load_dataset()
    results: dict[str, list[dict]] = {}
    for model in models:
        records = []
        for item in items:
            text = call_fn(model, build_prompt(item))
            parsed = parse_answer(item["task"], text)
            records.append({"id": item["id"], "task": item["task"], "correct": score(item, parsed)})
        results[model] = records
    return results


def summarize(results: dict) -> dict:
    summary: dict[str, dict] = {}
    for model, records in results.items():
        tasks: dict[str, dict] = {}
        for rec in records:
            t = tasks.setdefault(rec["task"], {"correct": 0, "total": 0})
            t["total"] += 1
            t["correct"] += int(rec["correct"])
        for t in tasks.values():
            t["accuracy"] = round(t["correct"] / t["total"], 3) if t["total"] else 0.0
        correct = sum(r["correct"] for r in records)
        total = len(records)
        summary[model] = {
            "by_task": tasks,
            "overall": {
                "correct": correct,
                "total": total,
                "accuracy": round(correct / total, 3) if total else 0.0,
            },
        }
    return summary


def markdown_table(summary: dict) -> str:
    tasks = sorted({t for s in summary.values() for t in s["by_task"]})
    header = "| model | " + " | ".join(tasks) + " | overall |"
    sep = "|" + "---|" * (len(tasks) + 2)
    rows = [header, sep]
    for model, s in summary.items():
        cells = [f"{s['by_task'].get(t, {}).get('accuracy', float('nan')):.0%}" for t in tasks]
        rows.append(f"| `{model}` | " + " | ".join(cells) + f" | {s['overall']['accuracy']:.0%} |")
    return "\n".join(rows)


_CLIENT: object | None = None


def _anthropic_call(model: str, prompt: str) -> str:
    from anthropic import Anthropic  # lazy: only needed for a live run

    global _CLIENT
    if _CLIENT is None:
        _CLIENT = Anthropic()
    resp = _CLIENT.messages.create(
        model=model,
        max_tokens=300,
        temperature=0,
        messages=[{"role": "user", "content": prompt}],
    )
    return "".join(b.text for b in resp.content if getattr(b, "type", "") == "text")


def _plot(summary: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tasks = sorted({t for s in summary.values() for t in s["by_task"]}) + ["overall"]
    models = list(summary)
    width = 0.8 / max(len(models), 1)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for i, model in enumerate(models):
        vals = [summary[model]["by_task"].get(t, {}).get("accuracy", 0.0) for t in tasks[:-1]]
        vals.append(summary[model]["overall"]["accuracy"])
        ax.bar([x + i * width for x in range(len(tasks))], vals, width, label=model)
    ax.set_xticks([x + width * (len(models) - 1) / 2 for x in range(len(tasks))])
    ax.set_xticklabels(tasks, rotation=15, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("accuracy")
    ax.set_title("haploqtl-bench")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--models", nargs="+", default=DEFAULT_MODELS)
    parser.add_argument("--limit", type=int, default=None, help="Only run the first N items.")
    args = parser.parse_args(argv)

    items = load_dataset()
    if args.limit:
        items = items[: args.limit]
    results = run(args.models, _anthropic_call, dataset=items)
    summary = summarize(results)
    table = markdown_table(summary)
    print(table)

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "results.json").write_text(json.dumps(results, indent=2))
    (RESULTS / "summary.md").write_text(
        f"# haploqtl-bench results\n\n{len(items)} items.\n\n{table}\n"
    )
    _plot(summary, RESULTS / "accuracy.png")
    print(f"\nwrote results to {RESULTS}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
