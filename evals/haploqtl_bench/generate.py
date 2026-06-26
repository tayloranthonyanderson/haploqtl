#!/usr/bin/env python3
"""Generate haploqtl-bench — a small, verifiable benchmark of biology-data reasoning tasks
derived from the EB-9 case study. Deterministic (seeded); writes dataset.jsonl.

Tasks (both have verifiable ground truth computed from the data, no human judgement):

  marker_diagnosticity
      Given anonymized resistant- and susceptible-line ALT-allele dosages (0/1/2) at a SNP,
      decide whether the ALT allele is *diagnostic* of the resistant haplotype — present in
      every resistant line and absent from every susceptible line. Tests genotype reasoning.
      Ground truth is computed directly from the genotypes.

  min_interval
      Given several lines' shared-haplotype intervals, return the interval common to all of
      them (the intersection). Mirrors how overlapping introgressions fine-map a QTL. Ground
      truth is the interval arithmetic.

Run:  uv run python evals/haploqtl_bench/generate.py
"""

from __future__ import annotations

import json
import random
from pathlib import Path

import allel
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
VCF = ROOT / "data" / "SL4.0ch09_subset.vcf.gz"
OUT = Path(__file__).resolve().parent / "dataset.jsonl"

EB9 = (62452852, 63003000)
RESISTANT = [
    "191163",
    "191164",
    "191167",
    "191172",
    "191175",
    "191174H",
]  # Devon Surprise + pathway
SUSCEPTIBLE = ["191165", "SRR1572598", "ERR418112"]  # NC EBR 1, NC 84173, Brandywine


def marker_items(rnd: random.Random, n_pos: int = 20, n_neg: int = 20) -> list[dict]:
    callset = allel.read_vcf(
        str(VCF), fields=["samples", "variants/CHROM", "variants/POS", "calldata/GT"]
    )
    samples = list(callset["samples"])
    index = {s: i for i, s in enumerate(samples)}
    pos = callset["variants/POS"]
    chrom = np.char.upper(callset["variants/CHROM"].astype(str))
    dosage = allel.GenotypeArray(callset["calldata/GT"]).to_n_alt(fill=-1)
    res_idx = [index[s] for s in RESISTANT]
    sus_idx = [index[s] for s in SUSCEPTIBLE]

    region = np.where((chrom == "SL4.0CH09") & (pos >= EB9[0]) & (pos <= EB9[1]))[0]
    positives, negatives = [], []
    for vi in region:
        rd, sd = dosage[vi, res_idx], dosage[vi, sus_idx]
        if (rd < 0).any() or (sd < 0).any():
            continue
        diagnostic = bool((rd >= 1).all() and (sd == 0).all())
        item = (
            [int(x) for x in rd],
            [int(x) for x in sd],
            diagnostic,
            int(pos[vi]),
        )
        (positives if diagnostic else negatives).append(item)

    rnd.shuffle(positives)
    rnd.shuffle(negatives)
    chosen = positives[:n_pos] + negatives[:n_neg]
    rnd.shuffle(chosen)

    items = []
    for k, (rd, sd, label, position) in enumerate(chosen):
        items.append(
            {
                "id": f"marker_{k:03d}",
                "task": "marker_diagnosticity",
                "input": {"resistant_alt_dosage": rd, "susceptible_alt_dosage": sd},
                "answer": "yes" if label else "no",
                "meta": {"position": position},
            }
        )
    return items


def interval_items(rnd: random.Random, n: int = 12) -> list[dict]:
    # item 0 uses the real EB-9 lines' shared-haplotype extents (bp), as computed by haploqtl
    # clustering against Devon Surprise on the bundled panel.
    real = {
        "line A": [61130000, 65030000],
        "line B": [61130000, 63530000],
        "line C": [62530000, 63030000],
        "line D": [62330000, 62830000],
    }
    items = [
        {
            "id": "interval_000",
            "task": "min_interval",
            "input": {"intervals": real},
            "answer": [max(v[0] for v in real.values()), min(v[1] for v in real.values())],
        }
    ]
    for k in range(1, n):
        m = rnd.randint(3, 5)
        core_s = rnd.randint(62_000_000, 62_600_000)
        core_e = rnd.randint(62_800_000, 63_400_000)
        ivs = {}
        for j in range(m):
            ivs[f"line {j + 1}"] = [
                core_s - rnd.randint(0, 800_000),
                core_e + rnd.randint(0, 800_000),
            ]
        items.append(
            {
                "id": f"interval_{k:03d}",
                "task": "min_interval",
                "input": {"intervals": ivs},
                "answer": [max(v[0] for v in ivs.values()), min(v[1] for v in ivs.values())],
            }
        )
    return items


def main() -> None:
    rnd = random.Random(0)
    items = marker_items(rnd) + interval_items(rnd)
    with open(OUT, "w") as handle:
        for item in items:
            handle.write(json.dumps(item) + "\n")
    by_task: dict[str, int] = {}
    for item in items:
        by_task[item["task"]] = by_task.get(item["task"], 0) + 1
    print(f"wrote {len(items)} items to {OUT}")
    for task, count in by_task.items():
        print(f"  {task}: {count}")


if __name__ == "__main__":
    main()
