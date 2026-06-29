"""Live resolvers/judges for real eval runs (never used in CI).

Network + API-key dependent: NCBI E-utilities for PMID existence, and an optional
Anthropic-backed judge for whether a cited record supports a claim. Calibration uses
the deterministic ``structural_calibration`` from ``verifiers`` (no LLM judge needed).
"""

from __future__ import annotations

import json
import os
import urllib.parse
import urllib.request

_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def ncbi_pmid_resolver(pmid: str, timeout: float = 20.0) -> bool:
    """True if a PubMed ID resolves to a real record (NCBI E-utilities esummary)."""
    pmid = str(pmid).strip()
    if not pmid.isdigit():
        return False
    params = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json"})
    try:
        with urllib.request.urlopen(f"{_ESUMMARY}?{params}", timeout=timeout) as resp:
            data = json.loads(resp.read().decode())
    except Exception:
        return False
    rec = data.get("result", {}).get(pmid)
    return bool(rec) and "error" not in rec and bool(rec.get("title"))


def llm_citation_judge(model: str = "claude-haiku-4-5-20251001"):
    """Build a citation-support judge: ``(claim, pmid) -> bool`` backed by Anthropic.

    Fetches the PubMed title+abstract and asks the model whether it supports the claim.
    Optional — the v1 smoke test runs resolve-only (no judge).
    """
    import anthropic  # lazy; only for live runs

    client = anthropic.Anthropic(api_key=os.environ["ANTHROPIC_API_KEY"])

    def _abstract(pmid: str) -> str:
        params = urllib.parse.urlencode(
            {"db": "pubmed", "id": pmid, "rettype": "abstract", "retmode": "text"}
        )
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{params}"
        with urllib.request.urlopen(url, timeout=20.0) as resp:
            return resp.read().decode()[:4000]

    def judge(claim: str, pmid: str) -> bool:
        try:
            abstract = _abstract(pmid)
        except Exception:
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
