# `data/` — bundled fixture

### `SL4.0ch09_subset.vcf.gz`
A phased, imputed VCF of **780 tomato and wild-relative genomes**, restricted to a ~4 Mb
region of chromosome 9 (SL4.0 reference) spanning the **EB-9** early-blight resistance QTL
(≈ chr09:62.4–63.0 Mb). 18,224 biallelic variants. This is the minimal slice needed to
reproduce a fine-mapped EB-9 result and to keep the demo and CI fast.

- Reference genome: SL4.0 (Sol Genomics Network).
- Sequence sources: NCBI SRA BioProject **PRJNA790656** (newly sequenced accessions) plus
  publicly available accessions used in the study (SRP094624, SRP045767, ERP004618,
  SRP150040).
- Redistributed from the public analysis repository
  [masudermann/HaplotypeAnalysis_Visualization](https://github.com/masudermann/HaplotypeAnalysis_Visualization).

### `AccessionRename.txt`
Maps the numeric sample IDs in the VCF header to human-readable accession names.

See the source paper for full methods: Anderson *et al.* (2024), *The Plant Journal*
117(2):404–415, [doi:10.1111/tpj.16495](https://doi.org/10.1111/tpj.16495).
