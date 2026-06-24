# Worked example — EB-9 collar-rot resistance (chr09)

A complete run of the skill on the **EB-9** interval `ch09:62,452,852-63,002,852` (the refined
boundary from Anderson et al. 2024, donor *Devon Surprise*), using the repo's bundled fixture.

## Step 1 — Genes in the interval
```bash
python scripts/genes_in_interval.py --chrom ch09 --start 62452852 --end 63002852 --format tsv
```
Candidate-relevant genes (excerpt):

| gene_id | start | description |
|---|---|---|
| Solyc09g074510.3 | 62,519,836 | Tubby-like F-box protein |
| Solyc09g160100.1 | 62,625,587 | F-box protein |
| Solyc09g074740.2 | 62,770,894 | Cation efflux family protein |
| Solyc09g074750.3 | 62,776,208 | Metal tolerance protein C1 |
| Solyc09g074790.3 | 62,817,971 | Potassium transporter |
| Solyc09g074800.2 | 62,835,326 | Potassium transporter |
| Solyc09g074820.4 | 62,866,311 | Potassium transporter |
| Solyc09g074920.3 | 62,933,032 | 2-oxoglutarate / Fe(II)-dependent oxygenase |

These recover exactly the gene families the paper highlighted for EB-9.

## Step 2 — Protein function (live UniProt)
```bash
python scripts/gene_function.py Solyc09g074790
```
```json
[
  {
    "gene_id": "Solyc09g074790",
    "uniprot": "A0A3Q7J1B7",
    "protein_name": "Potassium transporter",
    "keywords": "Ion transport;Membrane;Potassium;Potassium transport;Reference proteome;Transmembrane;Transport"
  }
]
```

## Step 4 — Diagnostic MAS markers
```bash
python scripts/diagnostic_variants.py --vcf ../../data/SL4.0ch09_subset.vcf.gz \
    --chrom ch09 --start 62452852 --end 63002852
```
→ **185 markers** (resistant n=6, susceptible n=3). First rows:

| position | ref | alt | resistant_allele | resistant_carrier_fraction | susceptible_carrier_fraction |
|---|---|---|---|---|---|
| 62,599,611 | C | A | A | 1.0 | 0.0 |
| 62,599,671 | C | T | T | 1.0 | 0.0 |
| 62,600,497 | C | T | T | 1.0 | 0.0 |

The first marker, `62,599,611`, coincides with the paper's pairwise chromosome-painting lower
bound for EB-9 (`62,599,611–62,945,798`) — an independent check that the markers track the real
introgression.

## Step 5 — Synthesized candidate-gene report (example)

**EB-9 (chr09:62.45–63.00 Mb), donor Devon Surprise — candidate genes for collar-rot resistance**

1. **Potassium transporters (`Solyc09g074790/074800/074820`)** — three tandem K⁺ transporters
   (UniProt: ion transport, membrane). K⁺ flux drives stomatal closure, a front-line defense
   against foliar/stem pathogens; a strong, mechanistically coherent candidate for quantitative
   resistance.
2. **F-box proteins (`Solyc09g074510`, `Solyc09g160100`)** — ubiquitin-mediated turnover in
   hormone/defense signaling; plausible regulators of a defense response.
3. **Metal tolerance / cation efflux (`Solyc09g074750`, `Solyc09g074740`)** — metal homeostasis
   at the infection site; secondary candidates.
4. **2OG/Fe(II) oxygenase (`Solyc09g074920`)** — defense-related secondary metabolism; secondary.

**MAS markers:** 185 diagnostic SNPs (resistant-fixed, susceptible-absent); see table above.
Markers nearest the K⁺-transporter cluster are the most useful for selection.

**Caveats:** co-location is not causation; candidates are ITAG4.1/SL4.0; quantitative resistance
is unlikely to be a canonical R gene. Expression/functional work is needed to confirm causality.
