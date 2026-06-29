"""Live resolvers/judges for real eval runs (never used in CI).

The NCBI E-utilities code (PubMed search + PMID existence + abstract fetch) is owned by
the skill — ``skills/qtl-candidate-gene/scripts/pubmed.py`` — because the skill is the
portable artifact that must run standalone. The eval imports it from there rather than
keeping its own copy, the same way ``harness`` imports ``genes_in_interval``. This module
re-exports those helpers under the names the eval scorers expect, and adds the optional
Anthropic-backed citation-support judge (the only eval-only, API-key-dependent piece).
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path
from typing import Any

_PUBMED = (
    Path(__file__).resolve().parent.parent
    / "skills"
    / "qtl-candidate-gene"
    / "scripts"
    / "pubmed.py"
)


def _load_pubmed() -> Any:
    """Import the skill's standalone ``pubmed.py`` by path (canonical NCBI helpers)."""
    name = "qtl_skill_pubmed"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _PUBMED)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_pubmed = _load_pubmed()

# Re-exports: same NCBI code as the skill's Step 3 / verify gate, under eval-facing names.
pubmed_search = _pubmed.pubmed_search
pmid_exists = _pubmed.pmid_exists
ncbi_pmid_resolver = _pubmed.pmid_exists  # "does this PubMed ID exist?" — the citation scorer


def llm_citation_judge(model: str = "claude-haiku-4-5-20251001"):
    """Build a citation-support judge: ``(claim, pmid) -> bool`` backed by Anthropic.

    Fetches the PubMed title+abstract and asks the model whether it supports the claim.
    Optional — the default eval run is resolve-only (no judge).
    """
    import anthropic  # lazy; only for live runs

    client = anthropic.Anthropic(api_key=os.environ["ANTHROPIC_API_KEY"])

    def judge(claim: str, pmid: str) -> bool:
        abstract = _pubmed.efetch_abstract(pmid)
        if abstract is None:
            return False
        msg = client.messages.create(
            model=model,
            max_tokens=8,
            messages=[
                {
                    "role": "user",
                    "content": (
                        f"Claim: {claim}\n\nPubMed record:\n{abstract}\n\n"
                        "Does the record support the claim? Answer only YES or NO."
                    ),
                }
            ],
        )
        text = "".join(b.text for b in msg.content if b.type == "text").strip().upper()
        return text.startswith("YES")

    return judge
