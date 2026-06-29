"""PubMed retrieval via NCBI E-utilities — the skill's Step 3 grounding tool.

Standalone and stdlib-only, so it runs from the skill bundle with no extra deps.
The candidate-gene workflow uses it to cite PMIDs it has actually retrieved instead of
recalling them from memory (recalled PMIDs are frequently fabricated), and the Step 6
verify gate uses ``pmid_exists`` to strip any citation that doesn't resolve.

CLI (Step 3):

    python pubmed.py "tomato F-box protein defense Alternaria" --retmax 5

prints a JSON list of ``{"pmid", "title"}``. Importable helpers: ``pubmed_search``,
``pmid_exists``, ``efetch_abstract``.
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
import time
import urllib.parse
import urllib.request

_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_TOOL = "haploqtl-qtl-candidate-gene"

# NCBI caps anonymous E-utilities at 3 requests/sec. The workflow bursts calls (a search
# fires esearch+esummary, then the verify gate re-checks every cited PMID), so without
# spacing a real PMID can return HTTP 429 and look like it doesn't exist. A process-wide
# throttle + retry keeps existence checks honest.
_MIN_INTERVAL = 0.4
_lock = threading.Lock()
_last = 0.0


def _throttle() -> None:
    global _last
    with _lock:
        wait = _MIN_INTERVAL - (time.monotonic() - _last)
        if wait > 0:
            time.sleep(wait)
        _last = time.monotonic()


def _get(url: str, timeout: float, attempts: int = 3) -> str | None:
    """Throttled GET with backoff; returns the body, or None on transient failure.

    None means network/HTTP failure after all retries — distinct from a valid 200 that
    simply reports no record (which the caller inspects). Only the former is retried.
    """
    sep = "&" if "?" in url else "?"
    url = f"{url}{sep}tool={_TOOL}"
    for i in range(attempts):
        _throttle()
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return resp.read().decode()
        except Exception:
            if i + 1 == attempts:
                return None
            time.sleep(0.5 * (i + 1))
    return None


def pubmed_search(query: str, retmax: int = 5, timeout: float = 20.0) -> list[dict]:
    """Search PubMed (esearch -> esummary) and return ``[{"pmid", "title"}]`` for the
    top hits, or ``[]`` on no match / error."""
    query = str(query).strip()
    if not query:
        return []
    sp = urllib.parse.urlencode(
        {"db": "pubmed", "term": query, "retmax": max(1, int(retmax)), "retmode": "json"}
    )
    body = _get(f"{_ESEARCH}?{sp}", timeout)
    if body is None:
        return []
    ids = json.loads(body).get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    pp = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids), "retmode": "json"})
    body = _get(f"{_ESUMMARY}?{pp}", timeout)
    if body is None:
        return [{"pmid": i, "title": ""} for i in ids]
    result = json.loads(body).get("result", {})
    return [{"pmid": i, "title": result.get(i, {}).get("title", "")} for i in ids]


def pmid_exists(pmid: str, timeout: float = 20.0) -> bool:
    """True if a PubMed ID resolves to a real record (esummary)."""
    pmid = str(pmid).strip()
    if not pmid.isdigit():
        return False
    params = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json"})
    body = _get(f"{_ESUMMARY}?{params}", timeout)
    if body is None:
        return False
    rec = json.loads(body).get("result", {}).get(pmid)
    return bool(rec) and "error" not in rec and bool(rec.get("title"))


def efetch_abstract(pmid: str, max_chars: int = 4000, timeout: float = 20.0) -> str | None:
    """Fetch a record's title+abstract as plain text (efetch), or None on failure."""
    params = urllib.parse.urlencode(
        {"db": "pubmed", "id": str(pmid), "rettype": "abstract", "retmode": "text"}
    )
    body = _get(f"{_EFETCH}?{params}", timeout)
    return body[:max_chars] if body is not None else None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Search PubMed for real records (NCBI).")
    parser.add_argument("query", help="PubMed query, e.g. 'tomato F-box protein defense'")
    parser.add_argument("--retmax", type=int, default=5, help="max records (default 5)")
    args = parser.parse_args(argv)
    json.dump(pubmed_search(args.query, args.retmax), sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
