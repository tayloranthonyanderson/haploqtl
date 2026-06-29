"""Live resolvers/judges for real eval runs (never used in CI).

The NCBI E-utilities code (PubMed search + PMID existence + abstract fetch) and the
citation-support judge are owned by the **skill** — it is the portable artifact that must
run standalone — so the eval imports them from there rather than keeping its own copies,
the same way ``harness`` imports ``genes_in_interval``. This module re-exports those
helpers under the names the eval scorers expect.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any

_SCRIPTS = Path(__file__).resolve().parent.parent / "skills" / "qtl-candidate-gene" / "scripts"


def _load_skill(filename: str, modname: str) -> Any:
    """Import a standalone skill script by path (its canonical home; no package)."""
    if modname in sys.modules:
        return sys.modules[modname]
    spec = importlib.util.spec_from_file_location(modname, _SCRIPTS / filename)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


_pubmed = _load_skill("pubmed.py", "qtl_skill_pubmed")

# Re-exports: same NCBI code as the skill's Step 3 / verify gate, under eval-facing names.
pubmed_search = _pubmed.pubmed_search
pmid_exists = _pubmed.pmid_exists
ncbi_pmid_resolver = _pubmed.pmid_exists  # "does this PubMed ID exist?" — the citation scorer


def llm_citation_judge(model: str = "claude-haiku-4-5-20251001"):
    """Citation-support judge ``(claim, pmid) -> bool``, backed by the skill's ``cite_support``.

    Exactly the judge the skill's verify gate uses with ``--support-model`` (which returns a
    richer ``{supported, evidence}``), narrowed here to the bool the eval's citation scorer
    wants. Optional — the default eval run is resolve-only (no judge).
    """
    support = _load_skill("cite_support.py", "qtl_skill_cite_support").make_support_judge(
        model, abstract_fetcher=_pubmed.efetch_abstract
    )
    return lambda claim, pmid: bool(support(claim, pmid)["supported"])
