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
from typing import Any

CallFn = Callable[[str, str], str]
# Execute a model tool call: (tool_name, tool_input) -> result text fed back to the model.
ToolHandler = Callable[[str, dict], str]

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


def run_tool_loop(
    client: Any,
    model: str,
    prompt: str,
    tools: list[dict],
    handler: ToolHandler,
    *,
    max_iters: int = 6,
    max_tokens: int = 2048,
) -> str:
    """Drive an Anthropic tool-use loop to completion and return the final text.

    ``client`` is anything exposing ``.messages.create(...)`` — the real SDK client in
    live runs, a scripted fake in tests — so the loop mechanics are verified without a
    key or network. ``handler(name, input) -> str`` runs a tool call; its result is fed
    back as a ``tool_result``. Returns ``""`` on a refusal (``stop_reason == "refusal"``),
    matching the harness's ``no_response`` handling.
    """
    messages: list[dict] = [{"role": "user", "content": prompt}]
    for _ in range(max_iters):
        resp = client.messages.create(
            model=model, max_tokens=max_tokens, tools=tools, messages=messages
        )
        if resp.stop_reason == "refusal":
            return ""
        if resp.stop_reason != "tool_use":
            return "".join(b.text for b in resp.content if b.type == "text")
        messages.append({"role": "assistant", "content": resp.content})
        results = [
            {"type": "tool_result", "tool_use_id": b.id, "content": handler(b.name, b.input)}
            for b in resp.content
            if b.type == "tool_use"
        ]
        messages.append({"role": "user", "content": results})
    # Iteration budget exhausted — make one final call with no tools to force an answer.
    resp = client.messages.create(model=model, max_tokens=max_tokens, messages=messages)
    return "".join(b.text for b in resp.content if b.type == "text")


_PUBMED_TOOL: dict = {
    "name": "search_pubmed",
    "description": (
        "Search PubMed for real publications. Returns a JSON list of {pmid, title}. "
        "Call this to find PMIDs that actually exist before citing them."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "query": {
                "type": "string",
                "description": "PubMed query, e.g. 'tomato F-box protein Alternaria defense'",
            },
            "retmax": {"type": "integer", "description": "max results to return (default 5)"},
        },
        "required": ["query"],
    },
}

_RETRIEVAL_SUFFIX = (
    "\n\nYou have a `search_pubmed` tool. Before citing ANY PMID, call it and cite ONLY "
    "PMIDs it returns — do not cite PMIDs from memory, as unverified PMIDs are frequently "
    "wrong. If a search finds nothing relevant, cite no PMID for that claim rather than "
    "guessing. End with the JSON object requested above."
)


def anthropic_retrieval_provider(model: str, prompt: str) -> str:
    """Retrieval-augmented arm: the model searches PubMed (live) for real PMIDs.

    Same ``(model, prompt) -> text`` contract as :func:`anthropic_provider`, but gives
    the model a live ``search_pubmed`` tool so it grounds citations instead of recalling
    them — the v2 fix for the fabrication the closed-book arm measures.
    """
    try:
        import anthropic
    except ImportError as exc:  # pragma: no cover - exercised only in live runs
        raise RuntimeError(
            "the 'anthropic' package is required for live runs: uv sync --extra eval"
        ) from exc
    from evals.live import pubmed_search

    client = anthropic.Anthropic(api_key=os.environ["ANTHROPIC_API_KEY"])

    def handler(name: str, tool_input: dict) -> str:
        if name != "search_pubmed":
            return json.dumps({"error": f"unknown tool {name!r}"})
        hits = pubmed_search(str(tool_input.get("query", "")), int(tool_input.get("retmax") or 5))
        return json.dumps(hits)

    return run_tool_loop(client, model, prompt + _RETRIEVAL_SUFFIX, [_PUBMED_TOOL], handler)
