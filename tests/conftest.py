"""Shared test fixtures."""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="session")
def fixture_vcf() -> Path:
    """Path to the bundled 780-genome chr09 (EB-9 region) fixture."""
    path = ROOT / "data" / "SL4.0ch09_subset.vcf.gz"
    assert path.exists(), "bundled chr09 fixture is missing"
    return path
