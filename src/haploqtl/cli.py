"""Command-line interface for haploqtl."""

from __future__ import annotations

import logging
from dataclasses import asdict
from pathlib import Path

import click
import numpy as np
import pandas as pd

from . import __version__
from .accessions import load_name_map, resolve_samples
from .cluster import cluster_haplotypes
from .contrast import contrast_twoway
from .introgression import (
    call_interval,
    donor_block_summary,
    interval_reduction,
    refine_with_markers,
)
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


@main.command()
@click.argument("vcf", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--chrom", required=True, help="Chromosome label (e.g. ch09).")
@click.option(
    "--resistant", required=True, help="Comma-separated resistant accession names or IDs."
)
@click.option("--susceptible", required=True, help="Comma-separated susceptible names or IDs.")
@click.option(
    "--benchmark", default=None, help="Donor sample for block extents [default: first resistant]."
)
@click.option(
    "--clusters",
    "clusters_csv",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Precomputed cluster CSV; skips clustering the VCF.",
)
@click.option("--window", default=250_000, show_default=True, help="Sliding window size (bp).")
@click.option("--step", default=100_000, show_default=True, help="Window step size (bp).")
@click.option("--min-snps", default=10, show_default=True, help="Minimum SNPs per window.")
@click.option("--d-min", default=2.0, show_default=True, help="Minimum merge-distance threshold.")
@click.option(
    "--d-max", default=80.0, show_default=True, help="Maximum merge-distance (exclusive)."
)
@click.option("--d-step", default=10.0, show_default=True, help="Merge-distance grid step.")
@click.option("--max-gap", default=1, show_default=True, help="Window gaps tolerated within a run.")
@click.option(
    "--max-exceptions", default=0, show_default=True, help="Marker-refine exceptions per side."
)
@click.option("--prior", default=None, help="Prior interval 'start-end' for % reduction.")
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Write diagnostic-track + block-extents CSVs here.",
)
def introgression(
    vcf: Path,
    chrom: str,
    resistant: str,
    susceptible: str,
    benchmark: str | None,
    clusters_csv: Path | None,
    window: int,
    step: int,
    min_snps: int,
    d_min: float,
    d_max: float,
    d_step: float,
    max_gap: int,
    max_exceptions: int,
    prior: str | None,
    output_dir: Path | None,
) -> None:
    """Call a narrowed introgression interval + donor-block retention from a VCF.

    Clusters the VCF (or loads --clusters), runs the two-way diagnostic contrast, narrows the
    interval to the longest gap-tolerant diagnostic run, reports the % reduction vs a prior
    interval, summarizes per-line donor-block retention (and the fine-mapped core), and
    refines the boundary to SNP resolution with diagnostic markers.
    """
    if clusters_csv is not None:
        clusters = pd.read_csv(clusters_csv, dtype={"sample": str})
    else:
        if d_min >= d_max:
            raise click.BadParameter("--d-min must be less than --d-max")
        d_grid = np.arange(d_min, d_max, d_step, dtype=np.float64)
        if d_grid.size == 0:
            raise click.BadParameter("empty distance grid; check --d-min / --d-max / --d-step")
        data = load_genotypes(vcf)
        clusters = cluster_haplotypes(
            data, chrom, window=window, step=step, min_snps=min_snps, d_grid=d_grid
        )

    available = set(clusters["sample"].astype(str))
    name_map = load_name_map()
    res_ids, res_missing = resolve_samples(
        [s.strip() for s in resistant.split(",") if s.strip()], available, name_map
    )
    sus_ids, sus_missing = resolve_samples(
        [s.strip() for s in susceptible.split(",") if s.strip()], available, name_map
    )
    if res_missing or sus_missing:
        click.echo(f"warning: not in cluster table: {res_missing + sus_missing}", err=True)
    if not res_ids or not sus_ids:
        raise click.ClickException("no resistant or susceptible samples resolved")
    bench_id = res_ids[0]
    if benchmark:
        resolved, _ = resolve_samples([benchmark], available, name_map)
        if resolved:
            bench_id = resolved[0]

    track = contrast_twoway(clusters, res_ids, sus_ids)
    n_diag = int(track["diagnostic"].fillna(False).sum())
    call = call_interval(track, window=window, max_gap=max_gap)
    if call is None:
        click.echo(f"No diagnostic window found ({n_diag} diagnostic / {len(track)} windows).")
        return

    click.echo(f"\nTwo-way diagnostic windows: {n_diag} / {len(track)}")
    click.echo(
        f"Introgression interval ({chrom}): {call.start:,.0f}-{call.stop:,.0f} "
        f"({call.span:,.0f} bp, {call.n_windows} windows, max_gap={call.max_gap})"
    )

    if prior:
        parts = prior.split("-")
        if len(parts) != 2:
            raise click.BadParameter("--prior must be 'start-end'")
        p0, p1 = int(parts[0]), int(parts[1])
        red = interval_reduction((call.start, call.stop), (p0, p1))
        click.echo(
            f"Prior interval: {p0:,}-{p1:,} ({red.prior_span:,.0f} bp) -> "
            f"reduction {red.reduction * 100:.1f}%"
        )

    anchor = (call.center_start + call.center_stop) / 2
    summary = donor_block_summary(
        clusters, bench_id, res_ids, window=window, anchor=anchor, max_gap=max_gap
    )
    click.echo(
        f"\nDonor-block core (intersection, benchmark={bench_id}): "
        f"{summary.core_start:,.0f}-{summary.core_stop:,.0f} ({summary.core_span:,.0f} bp); "
        f"boundaries L={summary.left_boundary_line} R={summary.right_boundary_line}"
    )
    click.echo("Per-line retained donor-block extent:")
    for e in sorted(summary.per_line, key=lambda x: -x.extent):
        click.echo(
            f"  {e.sample:>12}  {e.start:>13,.0f}-{e.stop:<13,.0f} extent {e.extent:>10,.0f} ({e.n_windows} win)"
        )

    refined = refine_with_markers(call, vcf, chrom, res_ids, sus_ids, max_exceptions=max_exceptions)
    if refined:
        click.echo(
            f"\nSNP-refined interval: {refined['start']:,}-{refined['stop']:,} "
            f"({refined['n_markers']} diagnostic markers)"
        )

    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        track.to_csv(output_dir / "diagnostic_track.csv", index=False)
        pd.DataFrame([asdict(e) for e in summary.per_line]).to_csv(
            output_dir / "block_extents.csv", index=False
        )
        click.echo(f"\nWrote diagnostic_track.csv + block_extents.csv to {output_dir}/")
