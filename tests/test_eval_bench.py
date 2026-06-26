"""Hermetic tests for haploqtl-bench — no network or model calls.

A fake "oracle" model reconstructs the correct answer from the prompt, exercising the full
build_prompt -> parse_answer -> score -> summarize pipeline; a fake "wrong" model checks that
scoring actually discriminates.
"""

import importlib.util
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT / "evals" / "haploqtl_bench" / "bench.py"


def _load_bench():
    spec = importlib.util.spec_from_file_location("haploqtl_bench", BENCH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _oracle(_model: str, prompt: str) -> str:
    if "ALT-allele dosage" in prompt:
        lists = re.findall(r"\[([0-9, ]+)\]", prompt)
        res = [int(x) for x in lists[0].split(",")]
        sus = [int(x) for x in lists[1].split(",")]
        diagnostic = all(v >= 1 for v in res) and all(v == 0 for v in sus)
        return json.dumps({"answer": "yes" if diagnostic else "no"})
    pairs = re.findall(r":\s*(\d+)-(\d+)", prompt)
    lo = max(int(a) for a, _ in pairs)
    hi = min(int(b) for _, b in pairs)
    return json.dumps({"start": lo, "end": hi})


def test_dataset_well_formed():
    bench = _load_bench()
    items = bench.load_dataset()
    assert len(items) >= 40
    for item in items:
        assert {"id", "task", "input", "answer"} <= set(item)
        assert item["task"] in {"marker_diagnosticity", "min_interval"}


def test_parse_answer_handles_json_and_loose_text():
    bench = _load_bench()
    assert bench.parse_answer("marker_diagnosticity", '{"answer": "yes"}') == "yes"
    assert bench.parse_answer("marker_diagnosticity", "The answer is No.") == "no"
    assert bench.parse_answer("min_interval", '{"start": 62530000, "end": 62830000}') == [
        62530000,
        62830000,
    ]
    assert bench.parse_answer("min_interval", "shared: 62,530,000-62,830,000") == [
        62530000,
        62830000,
    ]


def test_oracle_scores_perfectly_and_wrong_model_does_not():
    bench = _load_bench()
    oracle = bench.summarize(bench.run(["oracle"], _oracle))
    assert oracle["oracle"]["overall"]["accuracy"] == 1.0

    wrong = bench.summarize(
        bench.run(["wrong"], lambda _m, _p: '{"answer": "yes", "start": 0, "end": 0}')
    )
    assert wrong["wrong"]["overall"]["accuracy"] < 1.0
