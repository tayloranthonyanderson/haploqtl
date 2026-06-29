# Interpreting candidate genes for a QTL trait

Guidance for Step 5 (synthesis). The goal is a *ranked, justified* shortlist — not a gene dump.

**The principle is trait-specific:** weight gene families with a plausible mechanistic link to *the trait you were given*, and down-weight housekeeping genes. Disease resistance is worked through below in detail (it is the source paper's case); the same logic applies to any trait — several others are sketched at the end.

## Ranking heuristics — disease resistance (worked example)
Weight gene families with plausible roles in (quantitative) disease resistance above housekeeping genes:

- **NBS-LRR / "R" / disease-resistance proteins** — classic resistance genes. Strong candidates, but note: canonical R genes usually give *qualitative* (hypersensitive) resistance; a *quantitative* QTL may instead be a weak/atypical allele or a different mechanism. Flag the tension explicitly.
- **Receptor-like kinases (RLKs) / LRR-kinases** — pattern recognition and defense signaling.
- **F-box / ubiquitin-ligase components** — regulate hormone (JA/SA) signaling and defense protein turnover.
- **Transporters tied to defense physiology** — e.g. potassium transporters (stomatal closure, the stomatal-immunity route exploited by foliar pathogens), metal/cation transporters (metal homeostasis at the infection site).
- **Oxidoreductases / 2OG-Fe(II) oxygenases, PR proteins, defense transcription factors** — secondary metabolism and defense gene regulation.
- **Down-weight** housekeeping/structural genes unless variant or expression evidence points to them.

Use the UniProt keywords and the ITAG description together; corroborate with literature (PubMed) where possible.

## Other traits (same logic)
The heuristic is always "families mechanistically linked to *this* trait." For non-disease traits, weight accordingly and down-weight housekeeping genes:
- **Fruit size / weight** — cell-number/expansion and cell-cycle regulators (e.g. CNR/`fw2.2`, KLUH/`fw3.2`, CSR/`fw11.3`), CLAVATA–WUSCHEL meristem genes (locule number/`fas`, `lc`).
- **Sugar / soluble solids, acidity** — invertases (e.g. LIN5), sugar and organic-acid transporters, central carbon metabolism.
- **Fruit shape** — OFP/`ovate`, IQD/`sun`, fasciation genes.
- **Color / ripening** — carotenoid-pathway enzymes (e.g. CYC-B/`Beta`), MADS-box ripening regulators (RIN), GOLDEN2-LIKE (`u`).
- **Abiotic stress** — ion/water transporters, LEA/dehydrins, stress-responsive TFs (DREB/NAC/WRKY).

Always corroborate with UniProt keywords + literature, and match the proposed mechanism to the trait's genetic architecture.

## Caveats to state in every report
- **Co-location is not causation.** Genes are candidates because they sit in the interval, not because they are proven causal.
- **Annotation version.** Candidates are ITAG4.1 on SL4.0; other versions differ.
- **Overlapping loci.** A resistance QTL can co-localize with an unrelated R gene from the same introgression. (In the source paper, EB-5 overlaps the *Rx-3* bacterial-spot locus from the same Hawaii 7998 introgression — so a "resistance protein" in the EB-5 interval may belong to *Rx-3*, not early-blight resistance.)
- **Quantitative vs qualitative.** Match the proposed mechanism to the trait's genetic architecture.

## MAS-marker conventions
- A diagnostic marker (`diagnostic_variants.py`) carries an allele present in all of one phenotype group (e.g. resistant donors, or high-trait lines) and absent from the contrasting group — i.e. it tracks the trait-associated haplotype.
- Report `chromosome, position, ref, alt, resistant_allele` and the carrier fractions. Markers nearest the strongest candidate genes (and spread across the interval) are most useful.
- For assay design (KASP/CAPS) the user will need ~50 bp of flanking sequence per marker from the SL4.0 reference FASTA — a downstream step, not produced here.

## Worked context from the source paper (Anderson et al. 2024)
- **EB-9** (stem/collar-rot resistance, donor *Devon Surprise*, `ch09`): candidate genes included an **F-box**, **potassium transporters**, a **metal tolerance protein**, a **cation efflux protein**, and an **Fe(II)-dependent oxygenase** — all recoverable from Step 1 in the EB-9 interval.
- **EB-5** (foliar resistance, donor *Hawaii 7998*, `ch05`): a **plant resistance protein** (`Solyc05g053980`) — but possibly the overlapping *Rx-3* gene rather than the early-blight QTL.
