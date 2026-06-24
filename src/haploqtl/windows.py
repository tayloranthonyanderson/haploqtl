"""Sliding-window generation over genomic positions."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt


@dataclass(frozen=True)
class Window:
    """A half-open-ish genomic window ``[start, stop]`` (both bounds inclusive of SNPs)."""

    start: int
    stop: int

    @property
    def center(self) -> float:
        """Window midpoint, used as the window's representative position."""
        return (self.start + self.stop) / 2


def iter_windows(positions: npt.NDArray[np.integer], window: int, step: int) -> Iterator[Window]:
    """Yield sliding windows spanning the observed positions.

    A window starts at every ``step`` from the first observed position up to and
    including the last observed position.

    Bug fix vs. the original implementation: the reference script computed the window
    count as ``int(((last - first) - window) / step)``, which silently dropped the final
    partial window and left the chromosome tail uncovered. Here windows are emitted while
    ``start <= last``, so the entire arm is covered (the last window may extend past the
    last SNP, which is harmless — it simply captures fewer SNPs).
    """
    if window <= 0 or step <= 0:
        raise ValueError("window and step must be positive integers")
    if positions.size == 0:
        return
    first = int(positions[0])
    last = int(positions[-1])
    start = first
    while start <= last:
        yield Window(start=start, stop=start + window)
        start += step
