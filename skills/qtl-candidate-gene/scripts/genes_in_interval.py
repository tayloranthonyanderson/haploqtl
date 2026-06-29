#!/usr/bin/env python3
"""List protein-coding genes overlapping a genomic interval, with functional annotation.

Source of truth: a bundled SGN ITAG4.1 gene-model slice (SL4.0 assembly / ITAG4.1
annotation), so gene IDs and coordinates match Anderson et al. (2024) exactly. Functional
descriptions are taken from the AHRD annotations on the mRNA records.

Usage:
    python genes_in_interval.py --chrom ch09 --start 62452852 --end 63002852
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

DEFAULT_GFF = Path(__file__).resolve().parent.parent / "data" / "ITAG4.1_EB_regions.gff3"


@dataclass
class Gene:
    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    description: str


def normalize_chrom(chrom: str) -> str:
    """Accept '9', 'ch9', 'ch09', or 'SL4.0ch09' and return the GFF seqid 'SL4.0ch09'."""
    if chrom.startswith("SL4.0ch"):
        return chrom
    digits = re.sub(r"\D", "", chrom)
    if not digits:
        raise ValueError(f"could not parse a chromosome number from {chrom!r}")
    return f"SL4.0ch{int(digits):02d}"


def _attr(attrs: str, key: str) -> str | None:
    match = re.search(rf"{key}=([^;]+)", attrs)
    return match.group(1) if match else None


def _open_text(path: Path):
    """Open a GFF for reading text, transparently handling gzip (``.gz``)."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def genes_in_interval(gff_path: Path, seqid: str, start: int, end: int) -> list[Gene]:
    """Return genes overlapping ``seqid:start-end``, annotated from their mRNA descriptions."""
    genes: dict[str, Gene] = {}
    descriptions: dict[str, str] = {}
    with _open_text(gff_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] != seqid:
                continue
            fstart, fend = int(fields[3]), int(fields[4])
            if fstart > end or fend < start:  # no overlap with the requested interval
                continue
            attrs = fields[8]
            if fields[2] == "gene":
                gene_id = (_attr(attrs, "ID") or "").replace("gene:", "")
                genes[gene_id] = Gene(gene_id, fields[0], fstart, fend, fields[6], "")
            elif fields[2] == "mRNA":
                parent = (_attr(attrs, "Parent") or "").replace("gene:", "")
                note = re.sub(r"\s*\(AHRD.*?\)", "", _attr(attrs, "Note") or "").strip()
                if parent and note and parent not in descriptions:
                    descriptions[parent] = note
    for gene_id, description in descriptions.items():
        if gene_id in genes:
            genes[gene_id].description = description
    return sorted(genes.values(), key=lambda g: g.start)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--chrom", required=True, help="Chromosome, e.g. ch09 or SL4.0ch09.")
    parser.add_argument("--start", type=int, required=True, help="Interval start (bp, SL4.0).")
    parser.add_argument("--end", type=int, required=True, help="Interval end (bp, SL4.0).")
    parser.add_argument(
        "--gff", type=Path, default=DEFAULT_GFF, help="ITAG4.1 GFF3 (default: bundled EB slice)."
    )
    parser.add_argument("--format", choices=["json", "tsv"], default="json")
    parser.add_argument("-o", "--output", type=Path, default=None, help="Output (default: stdout).")
    args = parser.parse_args(argv)

    if args.start > args.end:
        parser.error("--start must be <= --end")
    seqid = normalize_chrom(args.chrom)
    genes = genes_in_interval(args.gff, seqid, args.start, args.end)

    if args.format == "json":
        text = json.dumps([asdict(g) for g in genes], indent=2)
    else:
        header = ["gene_id", "seqid", "start", "end", "strand", "description"]
        lines = ["\t".join(header)]
        lines += [
            "\t".join([g.gene_id, g.seqid, str(g.start), str(g.end), g.strand, g.description])
            for g in genes
        ]
        text = "\n".join(lines)

    if args.output:
        args.output.write_text(text + "\n")
    else:
        print(text)
    print(f"{len(genes)} genes in {seqid}:{args.start}-{args.end}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
