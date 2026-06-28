# Data sources & provenance

## Gene models — SGN ITAG4.1 (offline, source of truth)
- Genes and coordinates come from the **ITAG4.1** annotation on the **SL4.0** assembly — the exact reference used in Anderson et al. (2024), so `Solyc` IDs and positions match the paper.
- The repo bundles a **slice** of the ITAG4.1 gene models over the EB-5 and EB-9 regions at `data/ITAG4.1_EB_regions.gff3` (gene + mRNA lines only). Functional descriptions are the AHRD annotations carried on the mRNA `Note=` attribute.
- Coordinates are 1-based bp on SL4.0. Seqids are `SL4.0ch00`…`SL4.0ch12`. Note that a matching VCF may spell the chromosome `SL4.0CH09` (uppercase) — the scripts match case-insensitively.

### Regenerating / extending the slice
Download the full annotation and slice the region(s) you need:
```bash
curl -O https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.1_release/ITAG4.1_gene_models.gff
awk -F'\t' '($3=="gene"||$3=="mRNA") && $1=="SL4.0ch09" && $4<=63200000 && $5>=62000000' \
    ITAG4.1_gene_models.gff
```
For intervals outside the bundled regions, pass a full ITAG4.1 GFF3 to `genes_in_interval.py --gff`.

## Protein function — UniProt (live REST)
- `gene_function.py` queries **UniProtKB** (`https://rest.uniprot.org`) for the tomato proteome (`organism_id:4081`).
- UniProt entries cross-reference the ITAG/EnsemblPlants `Solyc` transcript IDs, so a free-text search on the base Solyc ID (version stripped) resolves the entry. Not every ITAG gene is in UniProt; missing entries return `null` and the workflow falls back to the ITAG4.1 description.

## Literature — PubMed
- Look up supporting literature in PubMed for Step 3.

## Why not Ensembl Plants REST?
A live Ensembl Plants region query would be a natural fit, but as of this writing the EnsemblGenomes REST host (`rest.ensemblgenomes.org`) is not resolving, and the main `rest.ensembl.org` does not serve plant species (`solanum_lycopersicum` is unknown there). The authoritative, version-exact path is therefore the bundled SGN ITAG4.1 slice; UniProt provides the live-database integration.
