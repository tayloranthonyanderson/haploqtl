# `evals/` — candidate-gene interpretation eval

A small, hermetically-tested harness that measures how a model interprets a fine-mapped
QTL interval — the job the [`qtl-candidate-gene`](../skills/qtl-candidate-gene/) skill
specifies. It is the first component in this repo that actually invokes and scores a
model.

## What it measures — and what it does not

It scores **faithfulness** and **calibration**, *not* causal correctness:

- **Gene existence** — does the model name only genes that are actually in the interval,
  vs. hallucinating IDs? Deterministic, checked against the ITAG4.1 GFF.
- **Citation grounding** — do the PMIDs it cites resolve to real PubMed records (and,
  with the optional judge, support the claim)? Fabrication rate.
- **Calibration** — does it hedge appropriately (multiple candidates, "co-location is
  not causation", validation needed), or over-declare a single causal gene?

It deliberately does **not** score "did it find the true causal gene." Candidate genes
are hypotheses, and many functional alleles are wild-species introgressions absent from
the SL4.0 (Heinz) reference — so there is no trustworthy answer key, here or genome-wide.
This is an eval-engineering artifact, not a claim about tomato genetics.

## Design

Tail-interpreter: the harness precomputes the authoritative gene list for each interval
(from the bundled [`data/ITAG4.1_genes.gff3.gz`](../data/)) and asks the model to
interpret it; verifiers score the structured answer against those facts. The scorers are
deterministic given their inputs — they turn a model judgment into a number (a verifiable
reward). The agentic variant (give the model the tools, score the trajectory) is future
work.

## Items

[`items.jsonl`](items.jsonl) — 4 `(trait, interval)` items across 3 chromosomes and 3
trait types: the two early-blight loci from the source paper (EB-9, EB-5) plus two
independent loci (Bs4 bacterial spot, fw2.2 fruit weight). The literature supplies the
*questions* only; no paper's candidate gene is used as an answer key.

## Running it

**Hermetic (no key, no cost)** — the CI path, via mock providers:

```bash
uv run pytest tests/test_evals_verifiers.py tests/test_evals_harness_mock.py
```

**Live** — needs the optional dependency and a key:

```bash
uv sync --extra eval
export ANTHROPIC_API_KEY=...                          # or a gitignored .env
uv run python -m evals.run --model claude-opus-4-8 -o evals/RESULTS.md
```

Live runs use `evals.live.ncbi_pmid_resolver` (NCBI E-utilities) for citation resolution
and the deterministic `structural_calibration` scorer; the optional
`evals.live.llm_citation_judge` adds citation-support checking.

## Worked example

**Input** — the harness feeds the trait + the gene list for one item; the model returns a
ranked candidate set as JSON. EB-9 (`SL4.0ch09:62,400,000–62,950,000`, trait = early-blight
collar rot; 50 genes in the interval).

**What Opus 4.8 returned** (an actual run): `determinable: false` — *"The interval contains no
canonical resistance gene … co-location plus annotation alone cannot pinpoint a causal gene, so
functional validation is required."*

| candidate | conf | claimed function | cited PMIDs |
|---|---|---|---|
| Solyc09g074890 | low | RALF immune peptide | 28096186, 32554493 |
| Solyc09g074620 | low | secreted zinc metalloproteinase | 18931680 |
| Solyc09g074850 | low | glutathione S-transferase | 19010102 |
| Solyc09g074430 | low | flavin-monooxygenase (SAR) | 16778020 |
| Solyc09g074390 | low | F-box / SCF component | 24043848 |

**What the harness checks, and the result:**

- **Gene existence → 1.00** ✅ — all 5 named genes are really in the interval (checked against the
  ITAG4.1 GFF); none hallucinated.
- **Calibration → 1.00** ✅ — it hedged: `determinable=false`, several *low*-confidence candidates,
  and an explicit "not causation, needs validation" caveat (a single over-confident causal call
  would score 0).
- **Citation grounding → fabrication 0.50** ❌ — of the 6 cited PMIDs, only 3 resolve on PubMed;
  the other 3 do not exist.

So an answer that reasons *well* — sensible defense-related gene families, correctly hedged — still
**fabricated half its citations**. That is the quiet failure the harness exists to surface and
quantify, and it's why faithfulness and calibration are scored separately. (Illustrative single
run; numbers vary — see Caveats.)

## Caveats

- [`RESULTS.md`](RESULTS.md) is a **single run** per model. These models are not fully
  deterministic and `temperature` is not settable on them, so expect run-to-run variance
  (≈ ±0.1 on fabrication). Multi-seed confidence intervals are future work.
- Calibration uses a deterministic rubric (`structural_calibration`), not an LLM judge — by
  design, for reproducibility. The LLM citation-support judge is opt-in.
- Four items is a small set: read the numbers as a demonstration of the harness, not a
  definitive model ranking.
