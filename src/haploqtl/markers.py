"""Diagnostic (marker-assisted-selection) markers from a VCF over a genomic interval.

Lifted from the ``qtl-candidate-gene`` skill's ``diagnostic_variants.py`` so the package
owns the science. A diagnostic variant carries an allele present in every resistant sample
and absent from every susceptible sample (allowing a configurable number of exceptions for
recombinants / genotyping error) — the SNP-level form of
:func:`haploqtl.contrast.contrast_twoway`, and the high-resolution layer used to refine
introgression boundaries.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import allel
import numpy as np

from .accessions import BUNDLED_RENAME, load_name_map, resolve_samples


@dataclass(frozen=True)
class DiagnosticMarker:
    chromosome: str
    position: int
    ref: str
    alt: str
    resistant_allele: str
    resistant_carrier_fraction: float
    susceptible_carrier_fraction: float


def normalize_seqid(chrom: str) -> str:
    """Map a short chromosome label (``ch09``) to the SL4.0 VCF seqid (``SL4.0ch09``)."""
    if chrom.upper().startswith("SL4.0CH"):
        return chrom
    digits = re.sub(r"\D", "", chrom)
    return f"SL4.0ch{int(digits):02d}"


def find_diagnostic_markers(
    vcf_path: str | Path,
    seqid: str,
    start: int,
    end: int,
    resistant: list[str],
    susceptible: list[str],
    *,
    rename_map_path: str | Path = BUNDLED_RENAME,
    max_exceptions: int = 0,
) -> tuple[list[DiagnosticMarker], dict]:
    """Return ``(markers, meta)`` for diagnostic variants in ``seqid:start-end``.

    A marker carries an allele present in all resistant samples and absent from all
    susceptible samples, allowing up to ``max_exceptions`` violators per side.
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
    if callset is None or "calldata/GT" not in callset:
        raise ValueError(f"no genotype data could be read from {vcf_path}")

    samples = list(callset["samples"])
    available = set(samples)
    name_map = load_name_map(rename_map_path)
    res_ids, res_missing = resolve_samples(resistant, available, name_map)
    sus_ids, sus_missing = resolve_samples(susceptible, available, name_map)
    if not res_ids or not sus_ids:
        raise ValueError(
            f"no samples resolved (resistant missing={res_missing}, susceptible missing={sus_missing})"
        )

    index = {s: i for i, s in enumerate(samples)}
    res_idx = [index[c] for c in res_ids]
    sus_idx = [index[c] for c in sus_ids]
    n_res, n_sus = len(res_idx), len(sus_idx)

    seqid = normalize_seqid(seqid)
    pos = callset["variants/POS"]
    chrom = callset["variants/CHROM"]
    ref = callset["variants/REF"]
    alt = callset["variants/ALT"]
    dosage = allel.GenotypeArray(callset["calldata/GT"]).to_n_alt(fill=-1)

    # VCF CHROM casing (SL4.0CH09) can differ from the seqid (SL4.0ch09); match upper-cased.
    chrom_upper = np.char.upper(chrom.astype(str))
    region_idx = np.where((chrom_upper == seqid.upper()) & (pos >= start) & (pos <= end))[0]

    markers: list[DiagnosticMarker] = []
    for vi in region_idx:
        res_dose, sus_dose = dosage[vi, res_idx], dosage[vi, sus_idx]
        if (res_dose < 0).any() or (sus_dose < 0).any():  # missing in a target set
            continue
        res_alt, sus_alt = int((res_dose >= 1).sum()), int((sus_dose >= 1).sum())
        res_ref, sus_ref = int((res_dose <= 1).sum()), int((sus_dose <= 1).sum())
        if res_alt >= n_res - max_exceptions and sus_alt <= max_exceptions:
            allele, res_carr, sus_carr = str(alt[vi][0]), res_alt, sus_alt
        elif res_ref >= n_res - max_exceptions and sus_ref <= max_exceptions:
            allele, res_carr, sus_carr = str(ref[vi]), res_ref, sus_ref
        else:
            continue
        markers.append(
            DiagnosticMarker(
                chromosome=seqid,
                position=int(pos[vi]),
                ref=str(ref[vi]),
                alt=str(alt[vi][0]),
                resistant_allele=allele,
                resistant_carrier_fraction=round(res_carr / n_res, 3),
                susceptible_carrier_fraction=round(sus_carr / n_sus, 3),
            )
        )

    meta = {
        "resistant_n": n_res,
        "susceptible_n": n_sus,
        "resistant_missing": res_missing,
        "susceptible_missing": sus_missing,
        "variants_scanned": int(region_idx.size),
    }
    return markers, meta
