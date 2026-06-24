"""Command-line interface for haploqtl."""

from __future__ import annotations

import logging
from pathlib import Path

import click
import numpy as np

from . import __version__
from .cluster import cluster_haplotypes
from .io import load_genotypes


@click.group()
@click.version_option(__version__, prog_name="haploqtl")
@click.option("-v", "--verbose", is_flag=True, help="Enable debug logging.")
def main(verbose: bool) -> None:
    """haploqtl - local-ancestry haplotype clustering for QTL discovery."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(message)s",
        datefmt="%H:%M:%S",
    )


@main.command()
@click.argument("vcf", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--chrom", required=True, help="Chromosome label for the output (e.g. ch09).")
@click.option("--window", default=250_000, show_default=True, help="Sliding window size (bp).")
@click.option("--step", default=100_000, show_default=True, help="Window step size (bp).")
@click.option("--min-snps", default=10, show_default=True, help="Minimum SNPs per window.")
@click.option("--d-min", default=2.0, show_default=True, help="Minimum merge-distance threshold.")
@click.option(
    "--d-max", default=80.0, show_default=True, help="Maximum merge-distance threshold (exclusive)."
)
@click.option("--d-step", default=10.0, show_default=True, help="Merge-distance grid step.")
@click.option(
    "--n-components", default=2, show_default=True, help="Number of principal components to record."
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=None,
    help="Output CSV path [default: <chrom>_haploqtl_w<window>_s<step>.csv].",
)
def cluster(
    vcf: Path,
    chrom: str,
    window: int,
    step: int,
    min_snps: int,
    d_min: float,
    d_max: float,
    d_step: float,
    n_components: int,
    output: Path | None,
) -> None:
    """Cluster genomes into local haplotypes along a sliding window over a VCF."""
    if window <= 0 or step <= 0:
        raise click.BadParameter("window and step must be positive")
    if d_min >= d_max:
        raise click.BadParameter("--d-min must be less than --d-max")
    if min_snps < n_components:
        raise click.BadParameter("--min-snps must be >= --n-components")
    d_grid = np.arange(d_min, d_max, d_step, dtype=np.float64)
    if d_grid.size == 0:
        raise click.BadParameter("empty distance grid; check --d-min / --d-max / --d-step")

    data = load_genotypes(vcf)
    frame = cluster_haplotypes(
        data,
        chrom,
        window=window,
        step=step,
        min_snps=min_snps,
        d_grid=d_grid,
        n_components=n_components,
    )
    out_path = output or Path(f"{chrom}_haploqtl_w{window}_s{step}.csv")
    frame.to_csv(out_path, index=False)
    n_windows = int(frame["position"].nunique())
    click.echo(
        f"Wrote {len(frame):,} rows ({n_windows} windows x {data.n_samples} samples) to {out_path}"
    )
