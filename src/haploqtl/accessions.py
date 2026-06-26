"""Resolve human-readable accession names to VCF / cluster-table sample IDs.

The bundled rename map (``data/AccessionRename.txt``) is a two-column
``<sample_id>\\t<accession name>`` table. Sample IDs in the VCF header (and hence in the
cluster table emitted by :mod:`haploqtl.cluster`) sometimes drop a leading ``X`` present
in the map (map ``X191163`` -> sample ``191163``); other IDs (``SRR``/``ERR``/``HA``)
appear verbatim. :func:`resolve_samples` reconciles all three cases.
"""

from __future__ import annotations

from pathlib import Path

# Repo-bundled default. Core functions accept an explicit path so the package does not
# hard-depend on the repository layout when installed elsewhere.
BUNDLED_RENAME = Path(__file__).resolve().parents[2] / "data" / "AccessionRename.txt"


def load_name_map(path: str | Path = BUNDLED_RENAME) -> dict[str, str]:
    """Return an ``{accession_name: sample_id}`` mapping from the rename table."""
    mapping: dict[str, str] = {}
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                sample_id, name = parts
                mapping[name] = sample_id
    return mapping


def resolve_samples(
    names: list[str], available: set[str], name_map: dict[str, str] | None = None
) -> tuple[list[str], list[str]]:
    """Map accession names (or raw IDs) onto the IDs actually present in ``available``.

    Tries, in order: the rename map, an ``X``-stripped form of the mapped code, then the
    name verbatim (already an ID). Returns ``(resolved_ids, missing_names)``.
    """
    name_map = name_map or {}
    resolved: list[str] = []
    missing: list[str] = []
    for name in names:
        code = name_map.get(name, name)
        if code in available:
            resolved.append(code)
        elif code.startswith("X") and code[1:] in available:
            resolved.append(code[1:])
        elif name in available:
            resolved.append(name)
        else:
            missing.append(name)
    return resolved, missing
