<p align="center">
  <img src="assets/banner.svg" alt="haploqtl — local-ancestry inference for fine-mapping QTL and finding trait donors" width="100%">
</p>

# haploqtl

[![CI](https://github.com/tayloranthonyanderson/haploqtl/actions/workflows/ci.yml/badge.svg)](https://github.com/tayloranthonyanderson/haploqtl/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)

Local-ancestry inference that turns large whole-genome sequence libraries into fine-mapped QTL and predicted trait donors. Reproducible, typed, and tested.

`haploqtl` detects cryptic ancestral introgressions, narrows the genomic intervals of quantitative trait loci (QTL), traces those loci to their historical donor cultivars, and predicts the trait in untested gene-bank accessions — without the need for purpose-built mapping populations.

> *Personal project, built on my own time and equipment using publicly available or self-collected data. Not affiliated with, funded by, or derived from any employer's work, data, or systems.*

> **Provenance.** This project modernizes and extends the method published in
> Anderson *et al.* (2024), *The Plant Journal* — [doi:10.1111/tpj.16495](https://doi.org/10.1111/tpj.16495),
> on which I am first author. The original research scripts live at
> [masudermann/HaplotypeAnalysis_Visualization](https://github.com/masudermann/HaplotypeAnalysis_Visualization);
> the reference clustering script (which I authored) is vendored verbatim under [`legacy/`](legacy/)
> and is the baseline this repository is being rebuilt around.

## The method

Genomes are grouped into **local haplotypes** along a stepped, sliding genomic window using **Ward hierarchical agglomerative clustering**. The merge-distance threshold is not fixed: for each window it is **auto-tuned by maximizing the mean silhouette coefficient**, so the number of haplotype clusters emerges from the genetic variance present in that window. This needs no genetic map, no reference panel, and no pre-specified number of ancestral groups — and it scales to hundreds of genomes in hours. In the source paper it fine-mapped two early-blight resistance QTL (a 70% and 56% reduction in interval size), traced them to the cultivars *Devon Surprise* and *Hawaii 7998*, and predicted resistance that was then **experimentally confirmed** in gene-bank accessions.

## Pipeline

```mermaid
flowchart LR
  V["780-genome VCF"] --> C["Windowed Ward clustering (silhouette-tuned)"]
  C --> H["Local haplotypes"]
  H --> Q["Fine-mapped QTL + introgression decay"]
  Q --> S["qtl-candidate-gene skill"]
  S --> G["Genes - ITAG4.1"]
  S --> U["Function - UniProt"]
  S --> M["Diagnostic MAS markers"]
  G --> R["Breeder-facing report"]
  U --> R
  M --> R
```

## Results

Reproduced by this repository on the bundled 780-genome panels (chromosomes 9 and 5), now end-to-end from `haploqtl` — clustering, introgression calling, and interval-narrowing all in code:

- **EB-9 fine-mapped to a 70% reduction, from code.** The two-way diagnostic contrast narrows EB-9 to `SL4.0ch09:62.40–62.95 Mb` (550 kb) — a 70.4% reduction of the prior ~1.86 Mb interval, matching the paper's ~70%. The call is validated two ways: window-for-window against the original R contrast functions, and the underlying clustering against the legacy script (identical partitions wherever the merge-distance agrees).
- **Introgression decay traces EB-9 to a 1920s heirloom.** The locus traces to *Devon Surprise*; its haplotype block is shared across the resistant pedigree and absent from susceptibles. The donor block decays from ~4.15 Mb down to ~0.75 Mb in the modern *NC 1 CELBR* recombinants — the recombination that erodes the introgression over breeding history is exactly what fine-maps the locus: those shorter blocks define the narrowed boundaries, and the intersection of all the blocks is the fine-mapped core. (Per-line extents are measured directly; we don't fit a per-generation decay rate, as the pedigree depth isn't reliably known.)
- **Diagnostic MAS markers** populate the interval — SNPs present in every resistant line and absent from every susceptible control: 133 across the full nine-line two-way set, 185 across the smaller six-line set used by the candidate-gene skill; the SNP-level span refines the boundary to `62.60–62.94 Mb`.
- **Candidate genes** recovered from ITAG4.1 in the interval: potassium transporters, F-box proteins, a cation-efflux protein, a metal-tolerance protein, and an Fe(II)-dependent oxygenase.

**Experimental validation** (from the source paper, Anderson *et al.* 2024 — not recomputed here):

- Haplotype-predicted resistance was **confirmed *in planta***: predicted-resistant accessions had far lower *A. linariae* stem disease than predicted-susceptible (mean 3.2 vs 28.2; Welch *t*(31.3) = −5.8, *P* < 0.001).

## Sharpening EB-9 with a phenotyped recombinant

<p align="center">
  <img src="assets/eb9_finemapping_summary.svg" alt="EB-9 fine-mapping: narrowing the interval from 70.4% to 81.2% reduction using a phenotyped recombinant" width="100%">
</p>

The interval-narrowing isn't only a reproduction — it can go further than the published interval, entirely from the bundled data. [`examples/finemap_eb9.py`](examples/finemap_eb9.py) walks through it:

- **Reproduce the paper.** Resistant pedigree vs. the susceptible controls gives `SL4.0ch09:62.40–62.95 Mb` — a 70.4% reduction of the prior interval.
- **Then narrow it.** *Ailsa Craig* is a phenotyped-susceptible heirloom (32.7% mean stem lesion *in planta*; [`data/eb9_phenotypes.csv`](data/eb9_phenotypes.csv)) that nonetheless carries the *Devon Surprise* haplotype across the left flank. A susceptible carrying it proves that flank can't be sufficient — so adding it as a susceptible trims the left edge to 62.60 Mb: 81.2% reduction (350 kb).
- **Know when to stop.** Nine *resistant* lines look like they'd trim the right edge, but `BGV007871` carries Ailsa Craig's exact haplotype while being resistant (0.1% stem lesion) — its resistance comes from another source. Resistant recombinants are confounded, so the right edge holds. The susceptible-recombinant argument is robust; the resistant one isn't.

Phenotypes are from Anderson *et al.* (2024); the sequence data is the bundled fixture. Run it with `uv run python examples/finemap_eb9.py`.

## Generalizing to a second locus: EB-5

The engine isn't specialized to EB-9. [`examples/finemap_eb5.py`](examples/finemap_eb5.py) points the *same* clustering and two-way contrast at chromosome 5, with a different donor, trait, and panel — *Hawaii 7998*'s foliar (not collar-rot) resistance pathway (`HA7998 → OH7536 → OH08-7663`) against a disjoint susceptible set (`OH88119`, `NC 84173`, `Brandywine`) — and recovers the published EB-5 introgression: a 56.2% reduction of the prior QTL (refined span 500 kb, matching Anderson *et al.* 2024 Table S2 to within ~0.1 kb — though the recovered interval sits ~150 kb right of the paper's exact boundaries, a threshold sensitivity at the right edge), on a bundled 780-genome chr05 fixture.

Two details make it a genuine test rather than a copy of EB-9:

- **The introgression is larger than the statistical QTL interval.** In the paper the EB-5 haplotype runs ~350 kb past the prior QTL's lower bound, so the *refined QTL* is the haplotype intersected with the prior — `interval_reduction` measures the reduction on that overlap (a 500 kb refined span from an 850 kb haplotype). EB-9, by contrast, sat nested inside its prior.
- **Resistance is tissue-specific.** `OH08-7663` carries Hawaii 7998's foliar (EB-5) resistance yet is collar-rot (EB-9) susceptible (64.7% stem lesion). The EB-5 mapping cross (`CU151095-146 × OH08-7663`) deliberately paired the two complementary resistances, which is why the EB-5 and EB-9 donor sets barely overlap.

Run it with `uv run python examples/finemap_eb5.py` (~70 s).

## Status & roadmap

This repository is under active development. Phases 0–2 are complete: a typed, tested `haploqtl` package with a real CLI, plus an Agent Skill that interprets a fine-mapped interval into candidate genes and MAS markers — all reproducible from a clean `git clone`.

- [x] **Phase 0 — Foundations & provenance.** Packaged project, pinned environment, CI, bundled chr09 fixture, reproducible demo, vendored reference implementation.
- [x] **Phase 1 — Modernized core.** Reference script refactored into a typed, tested, documented `haploqtl` package with a real CLI. Two latent bugs in the reference fixed: the silhouette search no longer aborts to a fixed fallback threshold on a single degenerate distance, and the final genomic window is no longer dropped.
- [x] **Phase 1.5 — Downstream introgression layer (R → Python).** The published downstream analysis lived entirely in R (`visualize_haplotypes.Rmd`). Ported into the typed, tested package: `contrast` (the one/two-way diagnostic contrast), `introgression` (algorithmic interval-narrowing — the original did it by eye — plus per-line donor-block retention and the fine-mapped core), and `markers` (diagnostic SNPs). Exposed as `haploqtl introgression`. Validated by **window-for-window equivalence to the original R** and a new **clustering ↔ legacy equivalence test** (identical partitions where the merge-distance agrees).
- [x] **Phase 2 — Agent Skill.** [`qtl-candidate-gene`](skills/qtl-candidate-gene/) — interval → candidate genes (ITAG4.1) → protein function (live UniProt) → diagnostic MAS markers → breeder report. Authored in Anthropic's Agent Skill (`SKILL.md`) format.
- [ ] **Next — Evaluation benchmark.** A small, verifiable benchmark scoring model *judgment* on the candidate-gene workflow — hallucinated genes/functions, mis-ranked candidates, misread genotypes, ungrounded citations (in design).

## Quickstart

Requires [uv](https://docs.astral.sh/uv/) (`curl -LsSf https://astral.sh/uv/install.sh | sh`).

```bash
git clone https://github.com/tayloranthonyanderson/haploqtl
cd haploqtl
uv sync --extra dev                 # creates .venv and installs everything
uv run bash scripts/run_demo.sh     # reproduce a minimal EB-9 result on the fixture
```

The demo runs windowed haplotype clustering over the bundled **780-genome chr09 fixture** (~4 Mb spanning the EB-9 QTL) and writes per-window haplotype-cluster tables to `output/`. Run the test suite with `uv run pytest`.

### CLI

```bash
uv run haploqtl cluster data/SL4.0ch09_subset.vcf.gz \
    --chrom ch09 --window 250000 --step 100000 \
    --d-min 2 --d-max 80 --d-step 10 \
    --output output/ch09_haplotypes.csv
```

The merge-distance threshold is auto-tuned per window (`--d-min/--d-max/--d-step` set the search grid). Output is a tidy long table — one row per (window, sample) with columns `chromosome, position, sample, cluster, distance_threshold, PC1..PCk`. See `uv run haploqtl cluster --help` for all options.

Then call the introgression interval, % reduction, and per-line donor-block retention straight from the VCF:

```bash
uv run haploqtl introgression data/SL4.0ch09_subset.vcf.gz --chrom ch09 \
    --resistant "Devon Surprise,Cambell 1943,NC EBR 2,NC 1 CELBR A,NC 1 CELBR B,CU151095-146,CU151011-146,CU191357,CU201041" \
    --susceptible "NC EBR 1,Brandywine,NC 84173" \
    --prior 61819509-63679761
```

This clusters, runs the two-way diagnostic contrast, narrows EB-9 to the longest gap-tolerant diagnostic run (550 kb, **70.4% reduction**), summarizes per-line donor-block retention and the fine-mapped core, and refines the boundary to SNP resolution. Accession names are resolved to VCF IDs automatically. See `uv run haploqtl introgression --help`.

Or **paint** each accession where it shares a donor's haplotype along the chromosome — the modern equivalent of the original R chromosome-painting figure (prints a terminal painting; `--svg` writes a to-scale figure):

```bash
uv run haploqtl paint data/SL4.0ch09_subset.vcf.gz --chrom ch09 \
    --benchmark "Devon Surprise" --highlight 62400075-62950075 --svg eb9_painting.svg
```

<p align="center">
  <img src="assets/eb9_painting.svg" alt="Chromosome painting of the EB-9 region versus Devon Surprise, ordered by donor-haplotype extent" width="100%">
</p>

## Agent Skill: `qtl-candidate-gene`

An [Agent Skill](skills/qtl-candidate-gene/) (in Anthropic's `SKILL.md` format) that interprets a fine-mapped interval the way a geneticist would — turning coordinates into biology:

1. **Genes in the interval** — from a bundled SGN **ITAG4.1** slice (SL4.0), paper-exact `Solyc` IDs + functional descriptions.
2. **Protein function** — live **UniProt** REST query (the genes-to-function database step).
3. **Diagnostic MAS markers** — variants present in all resistant donors and absent from all susceptible controls, computed from the VCF (reproducing the paper's contrast).
4. **Synthesis** — a ranked, mechanistically-reasoned candidate-gene report + marker table.

On the EB-9 interval it recovers exactly the gene families the paper highlighted (potassium transporters, F-box, cation efflux, metal-tolerance, Fe(II)-oxygenase) and 185 diagnostic markers whose first position coincides with the paper's chromosome-painting boundary. See the [worked example](skills/qtl-candidate-gene/EXAMPLE.md).

## Repository layout

```
haploqtl/
├── src/haploqtl/      # io (VCF→dosage), windows, cluster (silhouette-tuned Ward),
│                      # contrast + introgression + markers (interval-narrowing &
│                      # donor-block retention), painting, accessions, cli
├── skills/            # qtl-candidate-gene Agent Skill (SKILL.md + scripts + references)
├── legacy/            # vendored, attributed reference script from the published paper
├── data/              # bundled chr09 + chr05 fixtures (780 genomes) + names + phenotypes
├── examples/          # finemap_eb9.py (sharpen EB-9) + finemap_eb5.py (generalize to EB-5)
├── scripts/           # run_demo.sh — reproduce a minimal EB-9 result
├── tests/             # unit, CLI, skill, legacy-baseline, R-equivalence (tests/r +
│                      # tests/fixtures golden) and clustering↔legacy tests
└── .github/workflows/ # CI: lint + format + type-check + test on Python 3.11 & 3.12
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
