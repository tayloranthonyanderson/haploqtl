---
name: qtl-candidate-gene
description: Interpret a fine-mapped QTL / introgression interval for a stated trait in tomato. Given an interval (SL4.0 coordinates) and the trait it was mapped for, lists the genes in the interval from SGN ITAG4.1, enriches them with protein function from UniProt, ranks candidate genes by plausibility for that trait, grounds every literature citation with a live PubMed search, optionally flags marker-assisted-selection (MAS) markers from a two-group phenotype contrast, then self-verifies the draft — dropping any gene not in the interval and any citation that doesn't resolve — and delivers a breeder-facing candidate-gene report with a verification stamp. Works for any tomato trait — disease resistance, fruit size/quality, color, plant architecture, abiotic-stress tolerance. Use when a user has a QTL interval plus a trait and wants candidate or causal genes, mechanistic hypotheses, or MAS/KASP markers, or wants to interpret haploqtl clustering output.
---

# Candidate-Gene Interpretation for Tomato QTL

Turn a fine-mapped QTL interval into a ranked, mechanistically-reasoned set of candidate genes plus marker-assisted-selection (MAS) markers — the manual triage a breeder/geneticist does with a genome browser, an annotation database, and the literature, done as a repeatable workflow.

**Target users:** plant breeders and geneticists who have localized a trait (e.g. via `haploqtl` local-ancestry clustering or traditional QTL mapping) and need to know *which genes* are in the interval, *what they do*, and *which markers* to use for selection.

## When to use this skill

Use when the user:
- Has a genomic interval on the SL4.0 tomato assembly (e.g. `ch09:62,452,852-63,002,852`) **plus the trait it was mapped for**, and wants the candidate genes.
- Asks for candidate or causal genes, candidate-gene hypotheses, or "what's in this QTL" for a given trait.
- Wants MAS / KASP markers diagnostic of a trait-associated haplotype.
- Has `haploqtl` cluster output and wants to interpret a locus biologically.

The trait can be anything mapped in tomato — disease resistance, fruit size/quality, color, architecture, abiotic-stress tolerance. It determines which gene families are weighted as candidates, so it is required, not assumed.

## Inputs

- **Required:** an interval — chromosome + start + end (SL4.0 bp).
- **Required:** the **trait** the QTL was mapped for (free text, e.g. "early-blight resistance", "fruit weight", "soluble solids content"). This drives candidate ranking.
- **Optional (for MAS markers):** a VCF of the relevant accessions plus two contrasting phenotype groups — e.g. resistant donors vs susceptible controls, or high- vs low-trait lines.

## Workflow

```
[ ] Step 1  Genes in the interval        (offline, ITAG4.1)
[ ] Step 2  Protein function             (live UniProt)        -- needs network
[ ] Step 3  Retrieve literature          (live PubMed)         -- needs network
[ ] Step 4  Diagnostic MAS markers       (offline, from VCF)   -- if accessions given
[ ] Step 5  Synthesize candidate report  (reasoning)
[ ] Step 6  Verify + stamp               (PubMed + ITAG4.1)    -- gate before delivery
```

### Step 1 — Genes in the interval
List the protein-coding genes overlapping the interval, with ITAG4.1 functional descriptions:

```bash
python scripts/genes_in_interval.py --chrom ch09 --start 62452852 --end 63002852
```

Returns JSON (`gene_id, seqid, start, end, strand, description`). The bundled source is an ITAG4.1 slice (SL4.0), so IDs and coordinates match the published annotation. For an interval outside the bundled EB-5/EB-9 regions, pass a full ITAG4.1 GFF3 with `--gff` (see `references/data-sources.md`).

### Step 2 — Protein function (UniProt)
Enrich the most promising genes with protein-level function from UniProtKB (live):

```bash
python scripts/gene_function.py Solyc09g074790 Solyc09g074510
```

**DECISION POINT:** this requires network access. If offline, skip and rely on the ITAG4.1 descriptions from Step 1. Not every Solyc gene has a UniProt entry; `null` fields are expected and fine.

### Step 3 — Retrieve literature (PubMed)
Search PubMed for evidence linking the candidate gene families (Step 1–2) to **the trait** — e.g. for disease resistance *"tomato F-box protein defense Alternaria"*; for fruit size *"tomato cell number regulator fruit weight"*:

```bash
python scripts/pubmed.py "tomato F-box protein defense Alternaria" --retmax 5
```

Returns real `{pmid, title}` records from NCBI E-utilities. **Cite only PMIDs this returns — never one recalled from memory** (recalled PMIDs are frequently fabricated). The Step 6 gate strips any cited PMID that doesn't resolve, so memory-citations get caught — but retrieve up front so the claim keeps its evidence. Needs network; if truly offline, cite no PMIDs rather than guessing — the report still stands on the ITAG4.1 + UniProt evidence.

### Step 4 — Diagnostic MAS markers (if accessions provided)
Find variants whose allele is present in one phenotype group and absent from the contrasting group — diagnostic markers for marker-assisted selection. The two groups are defined by the trait (resistant vs susceptible for a disease; high vs low for a quantitative trait). The script's flags are named `--resistant`/`--susceptible` for the canonical disease case, but they simply mean "carries the trait allele" vs "does not":

```bash
python scripts/diagnostic_variants.py --vcf <accessions.vcf.gz> \
    --chrom ch09 --start 62452852 --end 63002852 \
    --resistant "<trait-positive lines>" --susceptible "<trait-negative lines>"
```

With no group flags the defaults reproduce the EB-9 contrast from Anderson et al. (2024) (Devon Surprise pathway vs NC EBR 1 / NC 84173 / Brandywine). Override with `--resistant`/`--susceptible` (comma-separated accession names) and relax with `--max-exceptions`.

### Step 5 — Synthesize the candidate-gene report
Using `references/interpretation-guide.md`, reason over the assembled evidence and produce a report:
1. **Ranked candidate genes** — name only genes from Step 1 and cite only PMIDs retrieved in Step 3; weight gene families plausibly linked to **the stated trait** above housekeeping genes; justify each with its function + literature. Match the families to the trait (don't default to disease): e.g. disease resistance → R-genes/NBS-LRR, receptor-like kinases, F-box/ubiquitin signaling, defense transporters/TFs; fruit size/quality → cell-number/expansion regulators, invertases and sugar/acid transport, hormone signaling; abiotic stress → transporters, LEA/dehydrins, stress TFs. See `references/interpretation-guide.md`.
2. **Mechanistic hypothesis** per top candidate (how it could plausibly affect the trait).
3. **MAS markers** — the diagnostic marker table from Step 4, ready for assay design.
4. **Caveats** — overlapping loci, annotation version, that interval co-location is not causation.

### Step 6 — Verify the report (gate before delivery)
Run the deterministic gate on your draft before returning it. Save the draft as a JSON object (the `candidates` / `determinable` / `rationale` structure) and check it against the interval and PubMed:

```bash
python scripts/verify_report.py --report draft.json --chrom ch09 --start 62452852 --end 62950075
```

It **drops** any candidate whose Solyc ID is not in the interval, **strips** any cited PMID that doesn't resolve on PubMed (flagging that claim to re-retrieve — it never invents a replacement), **flags** an over-confident un-hedged call, and prints a **verification stamp**. Append the stamp to the report. If PMIDs were stripped, return to Step 3, retrieve real ones for those claims, and re-verify.

For the strongest check, add `--support-model <model>` (e.g. `claude-haiku-4-5-20251001`): for every PMID that resolves, the gate fetches its abstract and asks the model whether it *actually supports* the specific claim — stripping citations that are real but off-topic, and recording the supporting quote for the ones that hold. Without the flag the gate only confirms a PMID *resolves*, not that it *backs* the claim. This needs network, the `anthropic` package, and an `ANTHROPIC_API_KEY`; the rest of the gate stays deterministic and offline-capable.

```bash
python scripts/verify_report.py --report draft.json --chrom ch09 --start 62452852 --end 62950075 \
    --support-model claude-haiku-4-5-20251001
```

The stamp then also reports *PMIDs supporting the claim*, and each surviving citation carries the abstract quote that backs it — put that quote in the report so a reader can verify it at a glance. (To read an abstract yourself: `python scripts/pubmed.py --abstract <PMID>`.)

## Output

A markdown report: ranked candidates with evidence, mechanistic hypotheses, a MAS marker table, explicit caveats, and the Step 6 **verification stamp** (genes in interval · PMIDs resolved · *PMIDs supporting the claim* · calibration) confirming the report was grounded and self-checked. See `EXAMPLE.md` for a worked EB-9 run.

## References
- `references/data-sources.md` — data provenance (ITAG4.1, UniProt), coordinates, regenerating the gene slice, why Ensembl Plants REST is not used.
- `references/interpretation-guide.md` — how to rank candidate genes for the trait at hand (disease resistance worked through in detail, other traits sketched), and MAS-marker conventions.
