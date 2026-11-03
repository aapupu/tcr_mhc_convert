"""Conservative MHC name normalization helpers.

The public registry indexes both canonical names and aliases.  These helpers
normalize typography only; they never infer a nearest allele.
"""

from __future__ import annotations

import re
from typing import Iterable, Tuple

from .errors import MHCTransformError


_HLA_LOCUS = r"(?:A|B|C|E|F|G|DRA|DRB[1-5]|DQA1|DQB1|DPA1|DPB1)"
_HLA_COLON_RE = re.compile(
    rf"^(?:HLA[-_]?)?({_HLA_LOCUS})\*?([0-9]+(?::[0-9]+){{1,3}}[A-Z]?)$",
    re.IGNORECASE,
)
_HLA_LEGACY_RE = re.compile(
    rf"^(?:HLA[-_]?)?({_HLA_LOCUS})\*?([0-9]{{4}}|[0-9]{{6}}|[0-9]{{8}})([A-Z]?)$",
    re.IGNORECASE,
)


def normalize_mhc_class(value: str) -> str:
    normalized = str(value).strip().upper().replace("MHC", "")
    normalized = normalized.replace("CLASS", "").replace("-", "")
    if normalized in {"1", "I"}:
        return "I"
    if normalized in {"2", "II"}:
        return "II"
    raise MHCTransformError("mhc_class must be I or II")


def normalize_sequence(sequence: str, *, allow_x: bool = False) -> str:
    normalized = re.sub(r"[\s-]+", "", str(sequence)).upper()
    alphabet = set("ACDEFGHIKLMNPQRSTVWY" + ("X" if allow_x else ""))
    invalid = sorted(set(normalized) - alphabet)
    if not normalized:
        raise MHCTransformError("sequence is empty")
    if invalid:
        raise MHCTransformError(
            "sequence contains unsupported residue(s): " + ", ".join(invalid)
        )
    return normalized


def canonical_hla(name: str) -> str | None:
    """Return a canonical HLA spelling when ``name`` is unambiguous."""

    compact = re.sub(r"\s+", "", str(name)).upper()
    match = _HLA_COLON_RE.fullmatch(compact)
    if match:
        return f"HLA-{match.group(1)}*{match.group(2)}"
    match = _HLA_LEGACY_RE.fullmatch(compact)
    if not match:
        return None
    digits, suffix = match.group(2), match.group(3)
    fields = [digits[index : index + 2] for index in range(0, len(digits), 2)]
    fields[-1] += suffix
    return f"HLA-{match.group(1)}*{':'.join(fields)}"


def two_field_hla(name: str) -> str | None:
    canonical = canonical_hla(name)
    if canonical is None:
        return None
    locus, fields = canonical[4:].split("*", 1)
    parts = fields.split(":")
    if len(parts) < 2:
        return None
    second = re.match(r"[0-9]+", parts[1])
    if second is None:
        return None
    suffix = re.search(r"([A-Z])$", fields)
    return (
        f"HLA-{locus}*{parts[0]}:{second.group(0)}"
        f"{suffix.group(1) if suffix else ''}"
    )


def _component_key(value: str) -> str:
    value = value.strip().upper().replace("HLA-", "").replace("HLA_", "")
    return re.sub(r"[^A-Z0-9]", "", value)


def exact_name_key(name: str) -> str:
    """Normalize spelling while retaining the number of HLA allele fields."""

    text = str(name).strip()
    components = [part for part in re.split(r"\s*(?:/|\+)\s*", text) if part]
    if len(components) > 1:
        return "/".join(exact_name_key(part) for part in components)
    if ":" in text:
        hla = canonical_hla(text)
        if hla:
            return hla.removeprefix("HLA-")
    compact_hla = re.fullmatch(
        rf"(?:HLA[-_]?)?({_HLA_LOCUS})\*?[0-9]+[A-Z]?",
        re.sub(r"\s+", "", text),
        flags=re.IGNORECASE,
    )
    if compact_hla:
        # Separator-free HLA names live in their own namespace. They become
        # ambiguous once fields exceed two digits; keeping them separate
        # prevents legacy B9530 from contaminating explicit B*95:30.
        return "COMPACT:" + _component_key(text)
    return _component_key(text)


def name_keys(name: str) -> Tuple[str, ...]:
    """Return typography-insensitive keys without changing biological content."""

    text = str(name).strip()
    if not text:
        return ()
    components = [part for part in re.split(r"\s*(?:/|\+)\s*", text) if part]
    keys = {exact_name_key(text)}
    if len(components) > 1:
        keys.add("/".join(_component_key(part) for part in components))
    hla = canonical_hla(text)
    if hla:
        keys.add(_component_key(hla))
        two_field = two_field_hla(hla)
        if two_field:
            keys.add(_component_key(two_field))
    return tuple(sorted(key for key in keys if key))


def all_name_keys(names: Iterable[str]) -> Tuple[str, ...]:
    keys = set()
    for name in names:
        keys.update(name_keys(name))
    return tuple(sorted(keys))

