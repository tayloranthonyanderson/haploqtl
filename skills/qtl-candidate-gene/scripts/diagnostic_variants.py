#!/usr/bin/env python3
"""Find marker-assisted-selection (MAS) markers whose allele tracks the resistant haplotype.

Replicates the contrast from Anderson et al. (2024, Fig. 2): a diagnostic variant carries an
allele present in every resistant donor but absent from every susceptible control. By default
the resistant set is the EB-9 introgression pathway (Devon Surprise and its descendants) and
the susceptible controls are NC EBR 1, NC 84173, and Brandywine. ``--max-exceptions`` allows a
few lines to violate the rule (recombinants / genotyping error).

Usage:
    python diagnostic_variants.py --vcf SL4.0ch09_subset.vcf.gz \
        --chrom ch09 --start 62452852 --end 63002852
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import allel
import numpy as np

DEFAULT_RENAME = Path(__file__).resolve().parent.parent / "data" / "AccessionRename.txt"
DEFAULT_RESISTANT = [
    "Devon Surprise",
    "Cambell 1943",
    "NC EBR 2",
    "NC 1 CELBR A",
    "NC 1 CELBR B",
    "CU151095-146",
]
DEFAULT_SUSCEPTIBLE = ["NC EBR 1", "NC 84173", "Brandywine"]
FIELDNAMES = [
    "chromosome",
    "position",
    "ref",
    "alt",
    "resistant_allele",
    "resistant_carrier_fraction",
    "susceptible_carrier_fraction",
]


def normalize_chrom(chrom: str) -> str:
    if chrom.startswith("SL4.0ch"):
        return chrom
    digits = re.sub(r"\D", "", chrom)
    return f"SL4.0ch{int(digits):02d}"


def load_name_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                mapping[parts[1]] = parts[0]
    return mapping


def resolve_samples(
    names: list[str], name_map: dict[str, str], vcf_samples: set[str]
) -> tuple[list[str], list[str]]:
    """Map accession names to VCF sample IDs (the rename map uses an X-prefix the VCF strips)."""
    resolved, missing = [], []
    for name in names:
        code = name_map.get(name, name)
        if code in vcf_samples:
            resolved.append(code)
        elif code.startswith("X") and code[1:] in vcf_samples:
            resolved.append(code[1:])
        else:
            missing.append(name)
    return resolved, missing


def find_diagnostic_markers(
    vcf_path: Path,
    seqid: str,
    start: int,
    end: int,
    resistant_names: list[str],
    susceptible_names: list[str],
    rename_map_path: Path = DEFAULT_RENAME,
    max_exceptions: int = 0,
) -> tuple[list[dict], dict]:
    """Return (markers, meta).

    A marker carries an allele present in all resistant samples and absent from all susceptible
    samples, allowing up to ``max_exceptions`` violators on each side.
    """
    callset = allel.read_vcf(
        str(vcf_path),
        fields=[
            "samples",
            "variants/CHROM",
            "variants/POS",
            "variants/REF",
            "variants/ALT",
            "calldata/GT",
        ],
    )
    samples = list(callset["samples"])
    name_map = load_name_map(rename_map_path)
    res_codes, res_missing = resolve_samples(resistant_names, name_map, set(samples))
    sus_codes, sus_missing = resolve_samples(susceptible_names, name_map, set(samples))
    if not res_codes or not sus_codes:
        raise ValueError(
            f"no samples resolved (resistant missing={res_missing}, susceptible missing={sus_missing})"
        )

    index = {s: i for i, s in enumerate(samples)}
    res_idx = [index[c] for c in res_codes]
    sus_idx = [index[c] for c in sus_codes]
    n_res, n_sus = len(res_idx), len(sus_idx)

    pos = callset["variants/POS"]
    chrom = callset["variants/CHROM"]
    ref = callset["variants/REF"]
    alt = callset["variants/ALT"]
    dosage = allel.GenotypeArray(callset["calldata/GT"]).to_n_alt(fill=-1)

    # VCF CHROM casing (SL4.0CH09) can differ from the GFF seqid (SL4.0ch09); match case-insensitively
    chrom_upper = np.char.upper(chrom.astype(str))
    region_idx = np.where((chrom_upper == seqid.upper()) & (pos >= start) & (pos <= end))[0]
    markers = []
    for vi in region_idx:
        res_dose, sus_dose = dosage[vi, res_idx], dosage[vi, sus_idx]
        if (res_dose < 0).any() or (sus_dose < 0).any():  # skip missing in target sets
            continue
        # carrier counts for the ALT allele (dosage >= 1) and the REF allele (dosage <= 1)
        res_alt, sus_alt = int((res_dose >= 1).sum()), int((sus_dose >= 1).sum())
        res_ref, sus_ref = int((res_dose <= 1).sum()), int((sus_dose <= 1).sum())
        if res_alt >= n_res - max_exceptions and sus_alt <= max_exceptions:
            allele, res_carr, sus_carr = str(alt[vi][0]), res_alt, sus_alt
        elif res_ref >= n_res - max_exceptions and sus_ref <= max_exceptions:
            allele, res_carr, sus_carr = str(ref[vi]), res_ref, sus_ref
        else:
            continue
        markers.append(
            {
                "chromosome": seqid,
                "position": int(pos[vi]),
                "ref": str(ref[vi]),
                "alt": str(alt[vi][0]),
                "resistant_allele": allele,
                "resistant_carrier_fraction": round(res_carr / n_res, 3),
                "susceptible_carrier_fraction": round(sus_carr / n_sus, 3),
            }
        )
    meta = {
        "resistant_n": n_res,
        "susceptible_n": n_sus,
        "resistant_missing": res_missing,
        "susceptible_missing": sus_missing,
        "variants_scanned": int(region_idx.size),
    }
    return markers, meta


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--resistant", default=",".join(DEFAULT_RESISTANT))
    parser.add_argument("--susceptible", default=",".join(DEFAULT_SUSCEPTIBLE))
    parser.add_argument("--rename-map", type=Path, default=DEFAULT_RENAME)
    parser.add_argument(
        "--max-exceptions", type=int, default=0, help="Lines allowed to violate the rule per side."
    )
    parser.add_argument("-o", "--output", type=Path, default=None)
    args = parser.parse_args(argv)

    markers, meta = find_diagnostic_markers(
        args.vcf,
        normalize_chrom(args.chrom),
        args.start,
        args.end,
        [s.strip() for s in args.resistant.split(",") if s.strip()],
        [s.strip() for s in args.susceptible.split(",") if s.strip()],
        args.rename_map,
        args.max_exceptions,
    )
    if meta["resistant_missing"] or meta["susceptible_missing"]:
        missing = meta["resistant_missing"] + meta["susceptible_missing"]
        print(f"warning: not found in VCF: {missing}", file=sys.stderr)

    handle = open(args.output, "w", newline="") if args.output else sys.stdout
    try:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(markers)
    finally:
        if args.output:
            handle.close()
    print(
        f"{len(markers)} diagnostic MAS markers in {args.chrom}:{args.start}-{args.end} "
        f"(resistant n={meta['resistant_n']}, susceptible n={meta['susceptible_n']})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
