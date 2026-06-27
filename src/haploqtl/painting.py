"""Chromosome painting: per-line sharing of a donor (benchmark) haplotype across windows.

Driven by :func:`haploqtl.contrast.compare_to_benchmark` (the port of the original Rmd's
``compare_clusters``). Builds a tidy per-line / per-window sharing matrix and renders it as
either a terminal painting or a self-contained, to-scale SVG — the modern equivalent of the
ggplot facets in the original ``visualize_haplotypes.Rmd``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .contrast import compare_to_benchmark


@dataclass(frozen=True)
class Painting:
    """Per-line sharing of the benchmark haplotype across ordered windows.

    ``shares[sample]`` is a list aligned to ``positions`` with ``True`` (shares the
    benchmark cluster), ``False`` (different cluster), or ``None`` (uncalled window).
    """

    benchmark: str
    positions: list[float]
    samples: list[str]
    shares: dict[str, list[bool | None]]
    shared_fraction: dict[str, float]


def build_painting(
    clusters: pd.DataFrame,
    benchmark: str,
    *,
    samples: list[str] | None = None,
    order_by_extent: bool = True,
) -> Painting:
    """Build a :class:`Painting` of every sample (or ``samples``) vs the ``benchmark``."""
    wide = (
        compare_to_benchmark(clusters, benchmark)
        .pivot(index="sample", columns="position", values="shares_benchmark")
        .sort_index(axis=1)
    )
    positions = [float(p) for p in wide.columns]
    rows = list(wide.index) if samples is None else [s for s in samples if s in wide.index]

    shares: dict[str, list[bool | None]] = {}
    fraction: dict[str, float] = {}
    for s in rows:
        vals = [None if pd.isna(v) else bool(v) for v in wide.loc[s].tolist()]
        shares[s] = vals
        n_true = sum(1 for v in vals if v)
        fraction[s] = n_true / len(vals) if vals else 0.0

    if order_by_extent:
        rows = sorted(rows, key=lambda s: fraction[s], reverse=True)
    return Painting(benchmark, positions, rows, shares, fraction)


def _spacing(positions: list[float]) -> float:
    if len(positions) < 2:
        return 1.0
    return float(np.median(np.diff(positions)))


def render_ascii(
    painting: Painting,
    *,
    labels: dict[str, str] | None = None,
    tags: dict[str, str] | None = None,
    eb9: tuple[int, int] | None = None,
    label_width: int = 20,
) -> str:
    """Render the painting as a fixed-width terminal figure (``█`` = shares the benchmark)."""
    labels = labels or {}
    tags = tags or {}
    pos = painting.positions
    lines: list[str] = []

    # Mb ruler.
    marks = [" "] * len(pos)
    nums = [" "] * len(pos)
    for mb in range(int(pos[0] // 1e6) + 1, int(pos[-1] // 1e6) + 1):
        idx = min(range(len(pos)), key=lambda i: abs(pos[i] - mb * 1e6))
        marks[idx] = "|"
        for k, ch in enumerate(str(mb)):
            if idx + k < len(nums):
                nums[idx + k] = ch
    indent = " " * (label_width + 4)
    lines.append(indent + "".join(marks))
    lines.append(indent + "".join(nums) + "  Mb")

    def cell(v: bool | None) -> str:
        return " " if v is None else ("█" if v else "·")

    for s in painting.samples:
        bar = "".join(cell(v) for v in painting.shares[s])
        tag = tags.get(s, " ")
        name = labels.get(s, s)
        pct = round(painting.shared_fraction[s] * 100)
        lines.append(f"{tag} {name:>{label_width}}  │{bar}│ {pct:>3d}%")

    if eb9 is not None:
        band = "".join("▔" if eb9[0] <= p <= eb9[1] else " " for p in pos)
        lines.append(indent + band + f"  EB-9 core ({eb9[0] / 1e6:.2f}–{eb9[1] / 1e6:.2f} Mb)")
    return "\n".join(lines)


def render_svg(
    painting: Painting,
    *,
    labels: dict[str, str] | None = None,
    tags: dict[str, str] | None = None,
    eb9: tuple[int, int] | None = None,
    title: str | None = None,
) -> str:
    """Render the painting as a self-contained, to-scale SVG (teal = shares the benchmark)."""
    labels = labels or {}
    tags = tags or {}
    pos = painting.positions
    pmin, pmax = pos[0], pos[-1]
    half = _spacing(pos) / 2.0
    span = (pmax + half) - (pmin - half) or 1.0

    left, right, top, bottom, row_h = 196, 96, 70, 44, 19
    plot_w = 820
    width = left + plot_w + right
    n = len(painting.samples)
    height = top + n * row_h + bottom

    def x(p: float) -> float:
        return left + (p - (pmin - half)) / span * plot_w

    teal, faint, coral = "#2dd4bf", "#1d2c44", "#ef6a55"
    out: list[str] = [
        f'<svg width="100%" viewBox="0 0 {width} {height}" role="img" '
        'xmlns="http://www.w3.org/2000/svg" '
        'font-family="ui-sans-serif, system-ui, -apple-system, Segoe UI, Helvetica, Arial, sans-serif">',
        f'<rect width="{width}" height="{height}" fill="#0c1426"/>',
    ]
    ttl = title or f"Chromosome painting vs {labels.get(painting.benchmark, painting.benchmark)}"
    out.append(
        f'<text x="{left}" y="30" font-size="18" font-weight="700" fill="#eaf1f8">{ttl}</text>'
    )
    out.append(
        f'<text x="{left}" y="48" font-size="12" fill="#9bb0c4">teal = shares the donor haplotype'
        f" &#183; ordered by extent &#183; drawn to SL4.0 scale</text>"
    )

    # EB-9 highlight band (drawn under the rows).
    if eb9 is not None:
        bx0, bx1 = x(eb9[0]), x(eb9[1])
        out.append(
            f'<rect x="{bx0:.1f}" y="{top - 6}" width="{bx1 - bx0:.1f}" height="{n * row_h + 6}" '
            f'fill="{coral}" opacity="0.12"/>'
        )
        out.append(
            f'<text x="{(bx0 + bx1) / 2:.1f}" y="{top - 10}" text-anchor="middle" font-size="11" '
            f'fill="#f0a08f">EB-9</text>'
        )

    # Mb axis ticks.
    for mb in range(int(pmin // 1e6) + 1, int(pmax // 1e6) + 1):
        tx = x(mb * 1e6)
        out.append(
            f'<line x1="{tx:.1f}" y1="{top - 4}" x2="{tx:.1f}" y2="{top + n * row_h}" stroke="#243348" stroke-width="1"/>'
        )
        out.append(
            f'<text x="{tx:.1f}" y="{top + n * row_h + 16}" text-anchor="middle" font-size="11" fill="#7f93a8">{mb} Mb</text>'
        )

    # Rows: merge contiguous shared windows into single bars (compact and clean).
    bar_h = row_h - 5
    for i, s in enumerate(painting.samples):
        ry = top + i * row_h
        cy = ry + row_h / 2
        out.append(
            f'<rect x="{left}" y="{ry + 2.5:.1f}" width="{plot_w}" height="{bar_h:.1f}" fill="{faint}" opacity="0.5"/>'
        )
        shares = painting.shares[s]
        j = 0
        while j < len(pos):
            if shares[j]:
                k = j
                while k + 1 < len(pos) and shares[k + 1]:
                    k += 1
                bx0, bx1 = x(pos[j] - half), x(pos[k] + half)
                out.append(
                    f'<rect x="{bx0:.1f}" y="{ry + 2.5:.1f}" width="{bx1 - bx0:.1f}" height="{bar_h:.1f}" fill="{teal}"/>'
                )
                j = k + 1
            else:
                j += 1
        tag = tags.get(s)
        if tag in ("R", "S"):
            out.append(
                f'<circle cx="10" cy="{cy:.1f}" r="3.5" fill="{teal if tag == "R" else coral}"/>'
            )
        name = labels.get(s, s)
        out.append(
            f'<text x="{left - 10}" y="{cy + 4:.1f}" text-anchor="end" font-size="12" fill="#cdd9e6">{name}</text>'
        )
        out.append(
            f'<text x="{left + plot_w + 10}" y="{cy + 4:.1f}" font-size="11" fill="#7f93a8">{round(painting.shared_fraction[s] * 100)}%</text>'
        )

    # Legend.
    ly = height - 14
    out.append(
        f'<circle cx="{left + 6}" cy="{ly - 4}" r="3.5" fill="{teal}"/><text x="{left + 16}" y="{ly}" font-size="11" fill="#9bb0c4">resistant pedigree</text>'
    )
    out.append(
        f'<circle cx="{left + 150}" cy="{ly - 4}" r="3.5" fill="{coral}"/><text x="{left + 160}" y="{ly}" font-size="11" fill="#9bb0c4">susceptible control</text>'
    )
    out.append("</svg>")
    return "\n".join(out)
