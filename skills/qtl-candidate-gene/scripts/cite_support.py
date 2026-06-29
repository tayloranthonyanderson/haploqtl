"""Citation-support judge — does a cited paper actually *support* the claim?

The deterministic gate (`verify_report.py`) can only check that a PMID *resolves*; a real
but off-topic paper still passes. This adds the missing check: fetch the record's abstract
and ask a model whether it supports the specific claim, returning a verdict plus a short
supporting quote to put in the report.

This is the one piece that needs an LLM, so it is **opt-in** — `anthropic` is imported
lazily and only when a judge is actually built, and an `ANTHROPIC_API_KEY` is required.
The abstract fetcher is injected (so the eval can reuse this without a hard import of the
skill's `pubmed` module), and so is the client, which keeps the judge unit-testable.
"""

from __future__ import annotations

import json
import os
import re
from collections.abc import Callable
from typing import Any

AbstractFetcher = Callable[[str], "str | None"]
# (claim, pmid) -> {"supported": bool, "evidence": str}
SupportJudge = Callable[[str, str], dict]

_PROMPT = (
    "Claim: {claim}\n\nPubMed abstract (PMID {pmid}):\n{abstract}\n\n"
    "Does this abstract specifically support the claim? Reply with ONLY a JSON object:\n"
    '{{"supported": true|false, "evidence": "<short verbatim quote from the abstract that '
    'supports it, or empty if not supported>"}}'
)


def _parse(text: str) -> dict:
    match = re.search(r"\{.*\}", text, re.DOTALL)
    if match:
        try:
            obj = json.loads(match.group(0))
            return {
                "supported": bool(obj.get("supported")),
                "evidence": str(obj.get("evidence", ""))[:300],
            }
        except json.JSONDecodeError:
            pass
    return {"supported": text.strip().upper().startswith("YES"), "evidence": ""}


def make_support_judge(
    model: str,
    *,
    abstract_fetcher: AbstractFetcher,
    client: Any = None,
) -> SupportJudge:
    """Build a ``(claim, pmid) -> {"supported", "evidence"}`` judge backed by Anthropic.

    ``abstract_fetcher`` turns a PMID into abstract text (or None). ``client`` is any object
    with ``.messages.create(...)`` — pass a fake to test, or leave None to lazily construct
    a real Anthropic client (needs ``anthropic`` + ``ANTHROPIC_API_KEY``).
    """

    def judge(claim: str, pmid: str) -> dict:
        abstract = abstract_fetcher(pmid)
        if not abstract:
            return {"supported": False, "evidence": ""}
        nonlocal client
        if client is None:
            import anthropic  # lazy; only when a real judge actually runs

            client = anthropic.Anthropic(api_key=os.environ["ANTHROPIC_API_KEY"])
        msg = client.messages.create(
            model=model,
            max_tokens=200,
            messages=[
                {
                    "role": "user",
                    "content": _PROMPT.format(claim=claim, pmid=pmid, abstract=abstract),
                }
            ],
        )
        text = "".join(b.text for b in msg.content if b.type == "text")
        return _parse(text)

    return judge
