#!/usr/bin/env bash
# Reproduce a minimal EB-9 local-ancestry result from the bundled chr09 fixture
# (780 tomato genomes, ~4 Mb around the EB-9 QTL on chromosome 9).
#
# Usage:  uv run bash scripts/run_demo.sh
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "${repo_root}/output"
cd "${repo_root}/output"

echo ">> Windowed Ward haplotype clustering on the chr09 fixture (paper-style 250kb/100kb window)..."
python "${repo_root}/legacy/cluster_haplotypes.py" \
    "${repo_root}/data/SL4.0ch09_subset.vcf.gz" ch09 \
    250000 100000 10 2 80 10

echo ""
echo ">> Done. Per-window haplotype clusters written to output/:"
ls -la "${repo_root}/output"/ch09_*.csv
