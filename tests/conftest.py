"""Shared test fixtures."""

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parent.parent
SKILL_SCRIPTS = ROOT / "skills" / "qtl-candidate-gene" / "scripts"


@pytest.fixture(scope="session")
def fixture_vcf() -> Path:
    """Path to the bundled 780-genome chr09 (EB-9 region) fixture."""
    path = ROOT / "data" / "SL4.0ch09_subset.vcf.gz"
    assert path.exists(), "bundled chr09 fixture is missing"
    return path


@pytest.fixture(scope="session")
def load_script():
    """Import a skill script (which is not part of the installable package) by file name."""

    def _load(name: str) -> ModuleType:
        path = SKILL_SCRIPTS / f"{name}.py"
        spec = importlib.util.spec_from_file_location(f"skill_{name}", path)
        assert spec and spec.loader, f"could not load {path}"
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module  # dataclasses/typing resolve annotations via sys.modules
        spec.loader.exec_module(module)
        return module

    return _load
