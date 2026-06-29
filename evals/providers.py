"""Model providers for the eval harness.

A provider is a ``CallFn``: ``(model, prompt) -> response_text``. The harness is
provider-agnostic, so CI runs a hermetic mock (no key, no cost, deterministic) while
live runs hit the real API.
"""

from __future__ import annotations

import json
import os
import re
from collections.abc import Callable

CallFn = Callable[[str, str], str]

_SOLYC = re.compile(r"Solyc\d{2}g\d{6}")


def _gene_ids_in_prompt(prompt: str) -> list[str]:
    """The Solyc gene IDs the harness listed in the prompt (deduplicated, in order)."""
    seen: dict[str, None] = {}
    for gid in _SOLYC.findall(prompt):
        seen.setdefault(gid, None)
    return list(seen)


def oracle_provider(model: str, prompt: str) -> str:
    """A faithful, well-calibrated mock answer reconstructed from the prompt.

    Names only genes actually present in the prompt, cites a placeholder PMID (the
    hermetic resolver accepts it), and hedges on determinability. Confirms the harness
    scores a *good* answer highly.
    """
    top = _gene_ids_in_prompt(prompt)[:3]
    answer = {
        "candidates": [
            {
                "solyc_id": gid,
                "claimed_function": "plausible defense/trait-related function",
                "pmids": ["10000001"],
                "confidence": "medium",
            }
            for gid in top
        ],
        "determinable": False,
        "rationale": (
            "Several genes in the interval are plausible candidates; co-location is "
            "not causation and functional validation is required."
        ),
    }
    return json.dumps(answer)


def garbage_provider(model: str, prompt: str) -> str:
    """A bad answer: a hallucinated gene, a fabricated PMID, and an overclaim.

    Confirms the harness penalizes hallucination, fabrication, and over-confidence.
    """
    answer = {
        "candidates": [
            {
                "solyc_id": "Solyc99g999999",  # not in any interval
                "claimed_function": "definitely THE causal gene",
                "pmids": ["00000000"],  # does not resolve
                "confidence": "high",
            }
        ],
        "determinable": True,
        "rationale": "This single gene is certainly the cause.",
    }
    return json.dumps(answer)


def anthropic_provider(model: str, prompt: str) -> str:
    """Live Anthropic provider. Needs the optional ``anthropic`` dep + ``ANTHROPIC_API_KEY``."""
    try:
        import anthropic
    except ImportError as exc:  # pragma: no cover - exercised only in live runs
        raise RuntimeError(
            "the 'anthropic' package is required for live runs: uv sync --extra eval"
        ) from exc
    client = anthropic.Anthropic(api_key=os.environ["ANTHROPIC_API_KEY"])
    message = client.messages.create(
        model=model,
        max_tokens=1024,
        messages=[{"role": "user", "content": prompt}],
    )
    return "".join(block.text for block in message.content if block.type == "text")
