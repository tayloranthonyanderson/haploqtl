---
name: qtl-candidate-gene
description: Interpret a fine-mapped QTL / introgression interval in tomato. Lists candidate genes in the interval from SGN ITAG4.1 (SL4.0), enriches them with protein function from UniProt, flags marker-assisted-selection (MAS) markers whose allele tracks the resistant haplotype, and drafts a breeder-facing candidate-gene report. Use when a user has a QTL interval or introgression boundaries (SL4.0 coordinates), wants candidate or causal genes, mechanistic hypotheses, MAS/KASP markers, or wants to interpret haplotype-clustering output from haploqtl.
---

# Candidate-Gene Interpretation for Tomato QTL

Turn a fine-mapped QTL interval into a ranked, mechanistically-reasoned set of candidate genes plus marker-assisted-selection (MAS) markers — the manual triage a breeder/geneticist does with a genome browser, an annotation database, and the literature, done as a repeatable workflow.

**Target users:** plant breeders and geneticists who have localized a trait (e.g. via `haploqtl` local-ancestry clustering or traditional QTL mapping) and need to know *which genes* are in the interval, *what they do*, and *which markers* to use for selection.

## When to use this skill

Use when the user:
- Has a genomic interval on the SL4.0 tomato assembly (e.g. `ch09:62,452,852-63,002,852`) and wants the genes in it.
- Asks for candidate or causal genes, candidate-gene hypotheses, or "what's in this QTL".
- Wants MAS / KASP markers diagnostic of a resistant haplotype.
- Has `haploqtl` cluster output and wants to interpret a locus biologically.

## Inputs

- **Required:** an interval — chromosome + start + end (SL4.0 bp).
- **Optional (for MAS markers):** a VCF of the relevant accessions plus the names of resistant donors and susceptible controls.

## Workflow

```
[ ] Step 1  Genes in the interval        (offline, ITAG4.1)
[ ] Step 2  Protein function             (live UniProt)        -- needs network
[ ] Step 3  Literature context           (PubMed)              -- optional
[ ] Step 4  Diagnostic MAS markers       (offline, from VCF)   -- if accessions given
[ ] Step 5  Synthesize candidate report  (reasoning)
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

### Step 3 — Literature context (optional)
If you have PubMed access, search for evidence linking the candidate gene families (Step 1–2) to disease resistance — e.g. *"tomato F-box protein defense Alternaria"*, *"potassium transporter stomatal immunity"*. Cite PMIDs in the report.

### Step 4 — Diagnostic MAS markers (if accessions provided)
Find variants whose allele is present in all resistant donors and absent from all susceptible controls — diagnostic markers for marker-assisted selection:

```bash
python scripts/diagnostic_variants.py --vcf <accessions.vcf.gz> \
    --chrom ch09 --start 62452852 --end 63002852
```

Defaults reproduce the EB-9 contrast from Anderson et al. (2024) (resistant = Devon Surprise pathway; susceptible = NC EBR 1, NC 84173, Brandywine). Override with `--resistant`/`--susceptible` (comma-separated accession names) and relax with `--max-exceptions`.

### Step 5 — Synthesize the candidate-gene report
Using `references/interpretation-guide.md`, reason over the assembled evidence and produce a report:
1. **Ranked candidate genes** — weight gene families implicated in (quantitative) disease resistance (R-genes/NBS-LRR, receptor-like kinases, F-box/ubiquitin signaling, transporters tied to stomatal/ion defense, defense TFs) above housekeeping genes; justify each with its function + literature.
2. **Mechanistic hypothesis** per top candidate (how it could plausibly affect the trait).
3. **MAS markers** — the diagnostic marker table from Step 4, ready for assay design.
4. **Caveats** — overlapping loci, annotation version, that interval co-location is not causation.

## Output

A markdown report: ranked candidates with evidence, mechanistic hypotheses, a MAS marker table, and explicit caveats. See `EXAMPLE.md` for a worked EB-9 run.

## References
- `references/data-sources.md` — data provenance (ITAG4.1, UniProt), coordinates, regenerating the gene slice, why Ensembl Plants REST is not used.
- `references/interpretation-guide.md` — how to rank candidate genes for quantitative disease resistance, and MAS-marker conventions.
