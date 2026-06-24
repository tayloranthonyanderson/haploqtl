#!/usr/bin/env bash
# Reproduce a minimal EB-9 local-ancestry result from the bundled chr09 fixture
# (780 tomato genomes, ~4 Mb around the EB-9 QTL on chromosome 9) using the haploqtl CLI.
#
# Usage:  uv run bash scripts/run_demo.sh
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "${repo_root}/output"

echo ">> haploqtl: windowed Ward haplotype clustering on the chr09 fixture (paper-style 250kb/100kb window)..."
haploqtl cluster \
    "${repo_root}/data/SL4.0ch09_subset.vcf.gz" \
    --chrom ch09 \
    --window 250000 --step 100000 --min-snps 10 \
    --d-min 2 --d-max 80 --d-step 10 \
    --output "${repo_root}/output/ch09_haplotypes.csv"

echo ""
echo ">> Done. Per-window haplotype clusters written to:"
ls -la "${repo_root}/output/ch09_haplotypes.csv"
