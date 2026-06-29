"""Live resolvers/judges for real eval runs (never used in CI).

Network + API-key dependent: NCBI E-utilities for PMID existence, and an optional
Anthropic-backed judge for whether a cited record supports a claim. Calibration uses
the deterministic ``structural_calibration`` from ``verifiers`` (no LLM judge needed).
"""

from __future__ import annotations

import json
import os
import threading
import time
import urllib.parse
import urllib.request

_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# NCBI caps anonymous E-utilities at 3 requests/sec. The eval bursts many calls
# (a search arm fires esearch+esummary per query, then the resolver re-checks every
# cited PMID), so without spacing a real PMID can return HTTP 429 and be miscounted as
# fabricated. A process-wide throttle + retry keeps the fabrication metric honest.
_NCBI_MIN_INTERVAL = 0.4
_ncbi_lock = threading.Lock()
_ncbi_last = 0.0


def _ncbi_throttle() -> None:
    global _ncbi_last
    with _ncbi_lock:
        wait = _NCBI_MIN_INTERVAL - (time.monotonic() - _ncbi_last)
        if wait > 0:
            time.sleep(wait)
        _ncbi_last = time.monotonic()


def _eutils(url: str, timeout: float, attempts: int = 3) -> str | None:
    """Throttled GET against E-utilities with backoff; returns the body or None.

    None means transient/network failure after all retries — distinct from a valid
    response that simply reports no record. Only the former should be retried; a real
    "PMID doesn't exist" comes back as a 200 the caller inspects.
    """
    for i in range(attempts):
        _ncbi_throttle()
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return resp.read().decode()
        except Exception:
            if i + 1 == attempts:
                return None
            time.sleep(0.5 * (i + 1))
    return None


def pubmed_search(query: str, retmax: int = 5, timeout: float = 20.0) -> list[dict]:
    """Search PubMed (NCBI E-utilities esearch -> esummary) for *real* records.

    Returns ``[{"pmid", "title"}]`` for the top hits, or ``[]`` on no match / error.
    This is the tool the retrieval-augmented eval arm gives the model so it cites
    PMIDs it has actually seen instead of recalling them from memory.
    """
    query = str(query).strip()
    if not query:
        return []
    sp = urllib.parse.urlencode(
        {"db": "pubmed", "term": query, "retmax": max(1, int(retmax)), "retmode": "json"}
    )
    body = _eutils(f"{_ESEARCH}?{sp}", timeout)
    if body is None:
        return []
    ids = json.loads(body).get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    pp = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids), "retmode": "json"})
    body = _eutils(f"{_ESUMMARY}?{pp}", timeout)
    if body is None:
        return [{"pmid": i, "title": ""} for i in ids]
    result = json.loads(body).get("result", {})
    return [{"pmid": i, "title": result.get(i, {}).get("title", "")} for i in ids]


def ncbi_pmid_resolver(pmid: str, timeout: float = 20.0) -> bool:
    """True if a PubMed ID resolves to a real record (NCBI E-utilities esummary)."""
    pmid = str(pmid).strip()
    if not pmid.isdigit():
        return False
    params = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json"})
    body = _eutils(f"{_ESUMMARY}?{params}", timeout)
    if body is None:
        return False
    rec = json.loads(body).get("result", {}).get(pmid)
    return bool(rec) and "error" not in rec and bool(rec.get("title"))


def llm_citation_judge(model: str = "claude-haiku-4-5-20251001"):
    """Build a citation-support judge: ``(claim, pmid) -> bool`` backed by Anthropic.

    Fetches the PubMed title+abstract and asks the model whether it supports the claim.
    Optional — the v1 smoke test runs resolve-only (no judge).
    """
    import anthropic  # lazy; only for live runs

    client = anthropic.Anthropic(api_key=os.environ["ANTHROPIC_API_KEY"])

    def _abstract(pmid: str) -> str | None:
        params = urllib.parse.urlencode(
            {"db": "pubmed", "id": pmid, "rettype": "abstract", "retmode": "text"}
        )
        body = _eutils(f"{_EFETCH}?{params}", 20.0)
        return body[:4000] if body is not None else None

    def judge(claim: str, pmid: str) -> bool:
        abstract = _abstract(pmid)
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
