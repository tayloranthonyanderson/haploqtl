# haploqtl-bench

A small, **verifiable** benchmark of biology-data reasoning, derived from the EB-9 case
study. Every item has ground truth computed directly from the data — no human judgement, no
LLM-as-judge — so scores are exact and reproducible. It measures the kind of reasoning a
breeder/geneticist needs from a model working over genotype and mapping data.

## Tasks

| task | items | what it tests | ground truth |
|---|---|---|---|
| `marker_diagnosticity` | 40 (balanced) | reasoning over a genotype table: is a SNP's ALT allele present in **every** resistant line and absent from **every** susceptible line (i.e. diagnostic for marker-assisted selection)? | computed from the bundled VCF over the EB-9 interval |
| `min_interval` | 12 | fine-mapping logic: given several lines' introgression intervals, return the interval common to all (the intersection) | interval arithmetic; item 0 uses the real EB-9 line extents |

Lines are anonymized (`resistant`/`susceptible`, `line A…`) so the task is genotype/interval
reasoning, not recall of which named cultivar is resistant.

## Run

```bash
uv sync --extra eval
export ANTHROPIC_API_KEY=...        # or put it in a gitignored .env and `source` it
uv run python evals/haploqtl_bench/bench.py            # all default models
uv run python evals/haploqtl_bench/bench.py --limit 5  # quick smoke test
```

Defaults to current Claude models (`claude-opus-4-8`, `claude-sonnet-4-6`,
`claude-haiku-4-5-20251001`); override with `--models`. Calls run at `temperature=0` for
reproducibility. Outputs land in `results/`: `results.json`, `summary.md`, and `accuracy.png`.

Regenerate the dataset (deterministic) with `uv run python evals/haploqtl_bench/generate.py`.

## Why these tasks

Both tasks isolate reasoning a breeder actually needs from a model working over genomic
data — marker diagnosticity and fine-mapping interval logic — while staying fully
verifiable, so a score is exact rather than a judgement call. The harness is
provider-agnostic (the model call is injected), so adding other providers or tasks is a
small change.

## Results

_Pending a live run — populated by `bench.py` into `results/summary.md`._
