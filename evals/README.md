# `evals/` — candidate-gene interpretation eval

A small, hermetically-tested harness that measures how a model interprets a fine-mapped
QTL interval — the job the [`qtl-candidate-gene`](../skills/qtl-candidate-gene/) skill
specifies. It is the one component in this repo that actually invokes and scores a model.

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
This is an eval-engineering artifact, not a claim about tomato genetics, and it ships
**without a published model leaderboard** — run it yourself.

## Design

Tail-interpreter: the harness precomputes the authoritative gene list for each interval
(from the bundled [`data/ITAG4.1_genes.gff3.gz`](../data/)) and asks the model to
interpret it; verifiers score the structured answer against those facts. The scorers are
deterministic given their inputs — they turn a model judgment into a number (a verifiable
reward). The agentic variant (give the model the tools, score the trajectory) is future
work.

Two arms share the same scoring:

- **closed-book** — the model cites PMIDs from memory.
- **retrieval-augmented** — the model is given a live `search_pubmed` tool (NCBI
  E-utilities) and instructed to cite only PMIDs the tool returns, so citations are
  grounded by construction rather than recalled.

## Items

[`items.jsonl`](items.jsonl) — 4 `(trait, interval)` items across 3 chromosomes and 3
trait types: the two early-blight loci from the source paper (EB-9, EB-5) plus two
independent loci (Bs4 bacterial spot, fw2.2 fruit weight). The literature supplies the
*questions* only; no paper's candidate gene is used as an answer key.

## Running it

**Hermetic (no key, no cost)** — the CI path, via mock providers:

```bash
uv run pytest tests/test_evals_verifiers.py tests/test_evals_harness_mock.py \
    tests/test_evals_retrieval.py
```

**Live** — needs the optional dependency and a key:

```bash
uv sync --extra eval
export ANTHROPIC_API_KEY=...                          # or a gitignored .env
uv run python -m evals.run --model claude-opus-4-8 --arm closed retrieval -o report.md
```

`--arm closed retrieval` runs both arms so you can compare grounded vs. recalled
citations; `-o` writes a local markdown report (not committed). Live runs use
`evals.live.ncbi_pmid_resolver` (NCBI E-utilities, throttled to NCBI's anonymous rate
limit so real PMIDs aren't miscounted as fabricated under burst load) and the
deterministic `structural_calibration` scorer; the optional `evals.live.llm_citation_judge`
adds citation-support checking.

## Notes

- Four items is a small **demonstration** set — the harness is the artifact, not a model
  ranking.
- Calibration uses a deterministic rubric (`structural_calibration`), not an LLM judge —
  by design, for reproducibility. The LLM citation-support judge is opt-in.
- These models are not fully deterministic and `temperature` is not settable on them, so
  per-run scores vary; that variance is one reason no results table is checked in here.
