<p align="center">
  <img src="assets/banner.svg" alt="haploqtl — local-ancestry inference for fine-mapping QTL and finding trait donors" width="100%">
</p>

# haploqtl

[![CI](https://github.com/tayloranthonyanderson/haploqtl/actions/workflows/ci.yml/badge.svg)](https://github.com/tayloranthonyanderson/haploqtl/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)

Local-ancestry inference that turns large whole-genome sequence libraries into fine-mapped QTL and predicted trait donors. Reproducible, typed, and tested.

`haploqtl` detects cryptic ancestral introgressions, narrows the genomic intervals of quantitative trait loci (QTL), traces those loci to their historical donor cultivars, and predicts the trait in untested gene-bank accessions — without the need for purpose-built mapping populations.

It doesn't stop at mapping. Each fine-mapped interval feeds an Agent Skill that turns it into a **grounded, self-verifying candidate-gene report** — driving real genome and literature databases, then checking its own genes and citations against them — with an evaluation harness that scores that step.

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
  S --> G["Genes · ITAG4.1"]
  S --> U["Function · UniProt"]
  S --> L["Literature · PubMed"]
  S --> M["Diagnostic MAS markers"]
  G --> VG{"Verify gate"}
  U --> VG
  L --> VG
  M --> VG
  VG --> R["Breeder-facing report + stamp"]
```

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

## Results

On the bundled 780-genome chr09 panel, `haploqtl` reproduces the paper's EB-9 result end to end — clustering, introgression calling, and interval-narrowing all in code:

- **EB-9 fine-mapped to a 70.4% reduction** (to `SL4.0ch09:62.40–62.95 Mb`, 550 kb from the prior ~1.86 Mb interval), validated two ways: window-for-window against the original R contrast, and the clustering against the legacy script (identical partitions wherever the merge-distance agrees).
- **Traced to the 1920s heirloom *Devon Surprise*** — its donor haplotype is shared across the resistant pedigree and absent from susceptibles; 185 diagnostic MAS markers refine the boundary to `62.60–62.94 Mb`.

In the source paper (Anderson *et al.* 2024, not recomputed here), haplotype-predicted resistance was **confirmed *in planta*** (predicted-resistant vs -susceptible stem disease, mean 3.2 vs 28.2; Welch *t*(31.3) = −5.8, *P* < 0.001).

Two runnable examples go further, entirely from the bundled data (each script's docstring walks through the reasoning):

- [`examples/finemap_eb9.py`](examples/finemap_eb9.py) — sharpens EB-9 from a 70% to an **81% reduction** by adding a phenotyped-susceptible recombinant (*Ailsa Craig*), and shows why the resistant recombinants can't tighten the other edge.
- [`examples/finemap_eb5.py`](examples/finemap_eb5.py) — points the same clustering and two-way contrast at a **second locus** (a different donor, trait, and chromosome) and recovers the published EB-5 introgression (~56% reduction).

Run with `uv run python examples/finemap_eb9.py` (or `finemap_eb5.py`).

## Candidate-gene interpretation (Agent Skill + eval)

An [Agent Skill](skills/qtl-candidate-gene/) (in Anthropic's `SKILL.md` format) that takes a fine-mapped interval **plus the trait it was mapped for** and returns two deliverables:

- **Candidate genes** — a shortlist of the interval's genes that could plausibly *underlie* the trait, with functions and mechanistic rationale. Narrows ~50 genes to a few hypotheses worth chasing, automating the gene-by-gene function and literature lookup a geneticist would do by hand.
- **Selection markers** — diagnostic SNPs that *track* the trait (present in the trait-positive lines, absent from the negatives), so a breeder can select by genotype alone — **even if the causal gene is never identified.**

The trait is an input, so this generalizes to any tomato QTL. What makes it usable rather than a plausible-sounding guess: the agent **drives real databases as tools** — ITAG4.1 for genes, UniProt for protein function, PubMed for the literature — and then **verifies its own draft** before returning it.

```mermaid
flowchart TB
  IN["Input: interval + trait (+ VCF)"]
  subgraph TOOLS["Agent drives real databases as tools"]
    direction TB
    S1["1 · Genes in interval — ITAG4.1"]
    S2["2 · Protein function — UniProt (live)"]
    S3["3 · Retrieve literature — PubMed (live)"]
    S4["4 · Diagnostic MAS markers — VCF"]
    S1 --> S2 --> S3 --> S4
  end
  IN --> S1
  S4 --> S5["5 · Synthesize: rank candidates, cite retrieved PMIDs"]
  S5 --> VG{"6 · Verify gate"}
  VG -->|"drop hallucinated genes · strip dead<br/>and unsupported PMIDs · flag overconfidence"| S5
  VG -->|pass| OUT["Report + stamp<br/>(grounded, cited, self-checked)"]
  VG -.->|"same verifiers, run as a benchmark"| EV["evals/"]
```

**Grounding by construction.** Step 3 retrieves real PMIDs from PubMed instead of recalling them from memory (recalled citations are frequently fabricated). The Step 6 **verify gate** then drops any gene not in the interval, strips any cited PMID that doesn't resolve on PubMed, optionally checks with an LLM that each abstract actually *supports* its claim — keeping a verbatim quote for the ones that do — flags an over-confident call, and stamps the report. It never invents a replacement: a stripped citation is sent back to Step 3.

**The eval is those same checks, measured.** [`evals/`](evals/) runs the verifiers as an offline benchmark rather than a per-run guardrail — **closed-book** (cite from memory) vs **retrieval-augmented** (live PubMed tool), scoring gene-existence, citation-grounding, and calibration. It is deliberately **not** a correctness check (candidates stay hypotheses; co-location isn't causation) and ships without a leaderboard, since per-run scores vary. One set of verifiers, two roles.

On the EB-9 interval it recovers exactly the gene families the paper highlighted (potassium transporters, F-box, cation efflux, metal-tolerance, Fe(II)-oxygenase) and 185 diagnostic markers whose first position coincides with the paper's chromosome-painting boundary. See the [worked example](skills/qtl-candidate-gene/EXAMPLE.md).

## Repository layout

```
haploqtl/
├── src/haploqtl/      # io (VCF→dosage), windows, cluster (silhouette-tuned Ward),
│                      # contrast + introgression + markers (interval-narrowing &
│                      # donor-block retention), painting, accessions, cli
├── skills/            # qtl-candidate-gene Agent Skill (SKILL.md + scripts + references)
├── evals/             # candidate-gene interpretation eval (faithfulness + calibration)
├── legacy/            # vendored, attributed reference script from the published paper
├── data/              # bundled chr09 + chr05 fixtures + ITAG4.1 GFF + names + phenotypes
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
