"""Hermetic tests for the retrieval-arm tool loop (no network, no API key).

``run_tool_loop`` is the only new live-path mechanic, so it is exercised here against a
scripted fake client: a tool call is dispatched to the handler, its result is threaded
back as a ``tool_result``, and the final text is returned. Refusal and the iteration-cap
fallback are checked too.
"""

import json

from evals.providers import run_tool_loop


class _Block:
    def __init__(self, type, **kw):
        self.type = type
        for k, v in kw.items():
            setattr(self, k, v)


class _Resp:
    def __init__(self, stop_reason, content):
        self.stop_reason = stop_reason
        self.content = content


class _FakeClient:
    """Returns scripted responses in order; records every create() call."""

    def __init__(self, scripted):
        self._scripted = list(scripted)
        self.calls = []

    @property
    def messages(self):
        return self

    def create(self, **kw):
        self.calls.append(kw)
        return self._scripted.pop(0)


def test_tool_loop_searches_then_answers():
    tool_call = _Block("tool_use", name="search_pubmed", input={"query": "tomato F-box"}, id="t1")
    final = _Block(
        "text",
        text='{"candidates": [{"solyc_id": "Solyc09g000010", "pmids": ["111"]}], '
        '"determinable": false, "rationale": "plausible"}',
    )
    client = _FakeClient([_Resp("tool_use", [tool_call]), _Resp("end_turn", [final])])
    seen = {}

    def handler(name, tool_input):
        seen["name"] = name
        seen["query"] = tool_input["query"]
        return json.dumps([{"pmid": "111", "title": "a real paper"}])

    out = run_tool_loop(client, "m", "PROMPT", [{"name": "search_pubmed"}], handler)

    assert seen == {"name": "search_pubmed", "query": "tomato F-box"}
    assert "111" in out  # final answer came back
    assert len(client.calls) == 2
    # the tool_result was threaded into the second request's message list
    second = client.calls[1]["messages"]
    assert any(
        m["role"] == "user"
        and isinstance(m["content"], list)
        and m["content"][0].get("type") == "tool_result"
        and m["content"][0]["tool_use_id"] == "t1"
        for m in second
    )


def test_tool_loop_refusal_returns_empty():
    client = _FakeClient([_Resp("refusal", [])])
    out = run_tool_loop(client, "m", "P", [{"name": "search_pubmed"}], lambda n, i: "")
    assert out == ""  # maps to status="no_response" in the harness, not a faithfulness 0


def test_tool_loop_caps_iterations_then_forces_answer():
    # Always asks to call a tool; the loop must stop and make a final tool-free call.
    looping = _Resp("tool_use", [_Block("tool_use", name="search_pubmed", input={}, id="x")])
    forced = _Resp("end_turn", [_Block("text", text="done")])
    client = _FakeClient([looping] * 6 + [forced])
    out = run_tool_loop(
        client, "m", "P", [{"name": "search_pubmed"}], lambda n, i: "[]", max_iters=6
    )
    assert out == "done"
    assert len(client.calls) == 7  # 6 capped iterations + 1 forced final call
    assert "tools" not in client.calls[-1]  # final call drops the tools to force an answer
