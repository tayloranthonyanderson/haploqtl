#!/usr/bin/env python3
"""Enrich tomato genes with protein-level function from UniProt (live REST query).

Demonstrates wiring to a standard biology database (UniProtKB) by gene identifier. UniProt
entries cross-reference the ITAG/EnsemblPlants Solyc IDs, so a free-text Solyc query against
the tomato proteome (organism 4081) resolves the entry and returns its accession, protein
name, and functional keywords. Use this to add protein-level context to the gene list from
genes_in_interval.py.

Network access required.

Usage:
    python gene_function.py Solyc09g074590 Solyc05g053980
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.parse
import urllib.request
from collections.abc import Callable

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
TOMATO_TAXON = 4081
FIELDS = "accession,protein_name,gene_names,keyword"
USER_AGENT = "haploqtl (github.com/tayloranthonyanderson/haploqtl)"

Fetcher = Callable[[str], str]


def fetch_uniprot(query: str, timeout: float = 30.0) -> str:
    params = urllib.parse.urlencode(
        {"query": query, "fields": FIELDS, "format": "tsv", "size": "1"}
    )
    request = urllib.request.Request(
        f"{UNIPROT_SEARCH}?{params}", headers={"User-Agent": USER_AGENT}
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode()


def gene_function(gene_id: str, fetch: Fetcher = fetch_uniprot) -> dict[str, str | None]:
    """Return UniProt function for one Solyc gene ID (None fields if no entry is found)."""
    query_id = gene_id.split(".")[0]  # strip the ITAG version suffix; UniProt indexes the base ID
    tsv = fetch(f"{query_id} AND organism_id:{TOMATO_TAXON}")
    rows = [r for r in tsv.splitlines() if r.strip()]
    if len(rows) < 2:
        return {"gene_id": gene_id, "uniprot": None, "protein_name": None, "keywords": None}
    record = dict(zip(rows[0].split("\t"), rows[1].split("\t"), strict=False))
    return {
        "gene_id": gene_id,
        "uniprot": record.get("Entry"),
        "protein_name": record.get("Protein names"),
        "keywords": record.get("Keywords"),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("gene_ids", nargs="+", help="Solyc gene IDs, e.g. Solyc09g074590.")
    parser.add_argument("-o", "--output", type=str, default=None, help="Output JSON path.")
    args = parser.parse_args(argv)

    results = []
    for gene_id in args.gene_ids:
        try:
            results.append(gene_function(gene_id))
        except Exception as exc:  # network/parse failure: report and continue
            print(f"warning: UniProt lookup failed for {gene_id}: {exc}", file=sys.stderr)
            results.append(
                {"gene_id": gene_id, "uniprot": None, "protein_name": None, "keywords": None}
            )

    text = json.dumps(results, indent=2)
    if args.output:
        with open(args.output, "w") as handle:
            handle.write(text + "\n")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
