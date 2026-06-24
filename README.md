# haploqtl

[![CI](https://github.com/tayloranthonyanderson/haploqtl/actions/workflows/ci.yml/badge.svg)](https://github.com/tayloranthonyanderson/haploqtl/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)

**Reproducible, AI-augmented local-ancestry inference for fine-mapping QTL and predicting trait donors from large genomic sequence libraries.**

`haploqtl` turns a library of whole-genome sequences into actionable breeding intelligence. It detects cryptic ancestral introgressions, narrows the genomic intervals of quantitative trait loci (QTL), traces those loci to their historical donor cultivars, and predicts the trait in untested gene-bank accessions — without the need for purpose-built mapping populations.

> **Provenance.** This project modernizes and extends the method published in
> Anderson *et al.* (2024), *The Plant Journal* — [doi:10.1111/tpj.16495](https://doi.org/10.1111/tpj.16495),
> on which I am first author. The original research scripts live at
> [masudermann/HaplotypeAnalysis_Visualization](https://github.com/masudermann/HaplotypeAnalysis_Visualization);
> the reference clustering script (which I authored) is vendored verbatim under [`legacy/`](legacy/)
> and is the baseline this repository is being rebuilt around.

## The method

Genomes are grouped into **local haplotypes** along a stepped, sliding genomic window using **Ward hierarchical agglomerative clustering**. The merge-distance threshold is not fixed: for each window it is **auto-tuned by maximizing the mean silhouette coefficient**, so the number of haplotype clusters emerges from the genetic variance present in that window. This needs no genetic map, no reference panel, and no pre-specified number of ancestral groups — and it scales to hundreds of genomes in hours. In the source paper it fine-mapped two early-blight resistance QTL (a 70% and 56% reduction in interval size), traced them to the cultivars *Devon Surprise* and *Hawaii 7998*, and predicted resistance that was then **experimentally confirmed** in gene-bank accessions.

## Status & roadmap

This repository is under active development. **Phase 0 is complete**: the pipeline is reproducible from a clean `git clone` against a bundled fixture, with packaging, tests, and CI in place.

- [x] **Phase 0 — Foundations & provenance.** Packaged project, pinned environment, CI, bundled chr09 fixture, reproducible demo, vendored reference implementation.
- [ ] **Phase 1 — Modernized core.** Refactor the reference script into a typed, tested, documented `haploqtl` package with a real CLI.
- [ ] **Phase 2 — Agent Skill.** A candidate-gene interpretation skill: interval → genes/annotations → literature-grounded hypotheses → marker-assisted-selection markers → breeder report.
- [ ] **Phase 3 — Database connector.** Wire the candidate-gene workflow to standard genomics databases (Ensembl / NCBI / Sol Genomics Network).
- [ ] **Phase 4 — Evaluation benchmark.** A verifiable-reward benchmark scoring agents on interval recovery, resistance prediction, and candidate-gene identification against the paper's validated ground truth.

## Quickstart

Requires [uv](https://docs.astral.sh/uv/) (`curl -LsSf https://astral.sh/uv/install.sh | sh`).

```bash
git clone https://github.com/tayloranthonyanderson/haploqtl
cd haploqtl
uv sync --extra dev                 # creates .venv and installs everything
uv run bash scripts/run_demo.sh     # reproduce a minimal EB-9 result on the fixture
```

The demo runs windowed haplotype clustering over the bundled **780-genome chr09 fixture** (~4 Mb spanning the EB-9 QTL) and writes per-window haplotype-cluster tables to `output/`. Run the test suite with `uv run pytest`.

## Repository layout

```
haploqtl/
├── src/haploqtl/      # package home (Phase 1 clean implementation lands here)
├── legacy/            # vendored, attributed reference script from the published paper
├── data/              # bundled chr09 fixture (780 genomes) + accession name map
├── scripts/           # run_demo.sh — reproduce a minimal EB-9 result
├── tests/             # pipeline smoke / schema-contract tests
└── .github/workflows/ # CI: lint + format + test on Python 3.11 & 3.12
```

## Citation

If you use this software or method, please cite the paper:

```bibtex
@article{anderson2024haploqtl,
  title   = {Detection of trait donors and {QTL} boundaries for early blight resistance
             using local ancestry inference in a library of genomic sequences for tomato},
  author  = {Anderson, Taylor A. and Sudermann, Martha A. and DeJong, Darlene M. and
             Francis, David M. and Smart, Christine D. and Mutschler, Martha A.},
  journal = {The Plant Journal},
  volume  = {117},
  number  = {2},
  pages   = {404--415},
  year    = {2024},
  doi     = {10.1111/tpj.16495}
}
```

## License

[MIT](LICENSE) © 2026 Taylor A. Anderson. The bundled fixture derives from publicly
available sequence data (NCBI SRA BioProject PRJNA790656 and other public accessions; see
[`data/README.md`](data/README.md)).
