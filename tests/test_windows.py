"""Tests for sliding-window generation, including the dropped-last-window bug fix."""

import numpy as np
import pytest

from haploqtl.windows import Window, iter_windows


def test_window_center():
    assert Window(100, 300).center == 200


def test_iter_windows_covers_the_last_position():
    # The original formula int(((last - first) - window) / step) yields 2 windows here
    # (starts 0, 300) and never covers the last position (1000). The fixed iterator emits
    # windows while start <= last, so the tail is covered.
    positions = np.array([0, 100, 1000], dtype=np.int64)
    windows = list(iter_windows(positions, window=200, step=300))
    assert [w.start for w in windows] == [0, 300, 600, 900]
    assert any(w.start <= 1000 <= w.stop for w in windows), "last position left uncovered"


def test_iter_windows_empty_positions():
    assert list(iter_windows(np.array([], dtype=np.int64), 100, 50)) == []


@pytest.mark.parametrize(("window", "step"), [(0, 10), (10, 0), (-1, 5)])
def test_iter_windows_rejects_nonpositive(window, step):
    with pytest.raises(ValueError):
        list(iter_windows(np.array([0, 10], dtype=np.int64), window, step))
