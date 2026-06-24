"""Loading genotypes from VCF into an alt-allele dosage matrix."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import allel
import numpy as np
import numpy.typing as npt

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class GenotypeData:
    """Genome-wide (within-chromosome) genotypes for a set of samples.

    Attributes:
        samples: sample identifiers, shape ``(n_samples,)``.
        positions: variant positions in bp, shape ``(n_variants,)``, ascending.
        dosage: alt-allele counts (0/1/2), shape ``(n_variants, n_samples)``.
    """

    samples: npt.NDArray[np.str_]
    positions: npt.NDArray[np.int64]
    dosage: npt.NDArray[np.int8]

    @property
    def n_samples(self) -> int:
        return int(self.samples.size)

    @property
    def n_variants(self) -> int:
        return int(self.positions.size)


def load_genotypes(vcf_path: str | Path) -> GenotypeData:
    """Read a (optionally gzipped) single-chromosome VCF into a dosage matrix.

    The method assumes phased/imputed, biallelic input with no missing calls, as produced
    by the upstream pipeline. Any missing calls are treated as homozygous reference.
    """
    path = str(vcf_path)
    log.info("Reading VCF: %s", path)
    callset = allel.read_vcf(path)
    if callset is None or "calldata/GT" not in callset:
        raise ValueError(f"no genotype data could be read from {path}")
    genotypes = allel.GenotypeArray(callset["calldata/GT"])
    dosage = genotypes.to_n_alt(fill=0).astype(np.int8)
    samples = np.asarray(callset["samples"], dtype=np.str_)
    positions = np.asarray(callset["variants/POS"], dtype=np.int64)
    log.info("Loaded %d variants x %d samples", positions.size, samples.size)
    return GenotypeData(samples=samples, positions=positions, dosage=dosage)
