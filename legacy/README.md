# `legacy/` — vendored reference implementation

`cluster_haplotypes.py` is the original research script from the paper, vendored here
**verbatim and unmodified** for provenance and reproducibility:

- Source: [masudermann/HaplotypeAnalysis_Visualization](https://github.com/masudermann/HaplotypeAnalysis_Visualization)
- Authored by **Taylor A. Anderson**.
- Published in Anderson *et al.* (2024), *The Plant Journal* 117(2):404–415, [doi:10.1111/tpj.16495](https://doi.org/10.1111/tpj.16495).

It is intentionally left as-is. The clean, typed, tested reimplementation now lives in
[`../src/haploqtl/`](../src/haploqtl/) (Phase 1 of the roadmap, complete); the test suite
checks that modern implementation against this reference script — identical clustering
partitions wherever the merge-distance agrees.

**Usage** (8 positional arguments):

```bash
python cluster_haplotypes.py [vcf] [chrom_basename] [window] [step] [min_snps] [min_d] [max_d] [step_d]
```
