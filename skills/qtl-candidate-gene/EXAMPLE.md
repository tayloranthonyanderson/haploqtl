# Worked example — EB-9 collar-rot resistance (chr09)

A complete run of the skill on the **EB-9** interval `ch09:62,452,852-63,002,852` (the refined
boundary from Anderson et al. 2024, donor *Devon Surprise*), trait = *early-blight collar rot*,
using the repo's bundled fixture. The point of the example is the last step: the report is
**grounded** (real genes, real citations) and **self-verified** before it ships.

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

## Step 3 — Retrieve literature (live PubMed)
```bash
python scripts/pubmed.py "F-box protein plant immunity defense ubiquitin" --retmax 3
```
```json
[
  {"pmid": "42010305", "title": "Salicylic acid modulates its catabolic enzymes via proteasomal degradation linked to SCF..."},
  {"pmid": "41288434", "title": "Viral Silencing Suppressor Activity in Plants Modifies Aphid Antiviral Immunity..."},
  {"pmid": "40810633", "title": "Post-Translational Modifications of TOE3 Regulate Antiviral Defense in Tobacco."}
]
```
Only PMIDs the search returns may be cited — never one recalled from memory. (Run the same
query per candidate family; the K⁺-transporter draft below cites `12481099` from its own search.)

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

## Step 5 — Synthesize the draft report
Rank candidates for the trait, citing only the PMIDs retrieved in Step 3. As a JSON draft:

```json
{
  "candidates": [
    {"solyc_id": "Solyc09g074790",
     "claimed_function": "Potassium transporter; K+ flux drives stomatal closure, a front-line defense against foliar/stem pathogens",
     "pmids": ["12481099"], "confidence": "low"},
    {"solyc_id": "Solyc09g074510",
     "claimed_function": "Tubby-like F-box protein; ubiquitin-mediated turnover in hormone/defense signaling",
     "pmids": ["42010305", "41288434"], "confidence": "low"}
  ],
  "determinable": false,
  "rationale": "Several plausible candidates; co-location is not causation, functional validation required."
}
```

## Step 6 — Verify the report (gate before delivery)
```bash
python scripts/verify_report.py --report draft.json --chrom ch09 --start 62452852 --end 63002852 \
    --support-model claude-haiku-4-5-20251001
```
Actual output of this run:

```
**Verification** — genes in interval 2/2 · PMIDs resolved 3/3 · PMIDs supporting the claim 2/3 · calibration: ok
- removed (resolve but don't support the claim): 12481099 — re-retrieve for those claims
```

What happened, and why it matters:

- **Both genes are real** in the interval (2/2) and **all 3 PMIDs resolve** (3/3) — no hallucinated
  IDs, no fabricated citations.
- The K⁺-transporter citation `12481099` resolves but its abstract is about calcium/ABA elicitor
  signaling, **not** a K⁺-transporter role in defense — so the support check **stripped it** and
  flagged the claim to re-retrieve. (Resolving ≠ supporting; this is the check that catches a
  real-but-off-topic citation.)
- The two F-box citations **survived with verbatim supporting quotes** attached to the report, e.g.
  `42010305`: *"…DMR6-ASSOCIATED F-BOX 1 (DAF1)…SCF-type E3 ligase-mediated proteasomal turnover of
  DMR6…modulating SA-mediated cell death."*

The delivered report carries that stamp, so a reader sees at a glance that every gene is in the
interval and every surviving citation actually backs its claim.

**Caveats** (always stated): co-location is not causation; candidates are ITAG4.1/SL4.0; a
quantitative QTL is unlikely to be a canonical R gene. Expression/functional work confirms causality
— the skill produces a vetted starting set, not an answer.
