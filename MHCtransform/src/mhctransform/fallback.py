"""Shared helpers for auditable partial and pseudo-scaffold fallbacks.

The strict registry never imports fallback rows into its exact indexes. This
module contains deterministic sequence/alignment utilities and policy
constants used by the offline fallback builder and runtime filtering.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, Iterable, Mapping, Sequence, Tuple

from .errors import MHCTransformError


FALLBACK_POLICIES = ("none", "partial", "imputed", "pseudo-scaffold")
PARTIAL_REPRESENTATION = "source-partial-v1"
IMPUTED_PRECURSOR_DEFINITION = "precursor-imputed-v1"
IMPUTED_AF3_DEFINITION = "af3-soluble-extracellular-imputed-v1"
EXACT_ECTODOMAIN_EVIDENCE = "source_exact_ectodomain"
PARTIAL_EVIDENCE = "source_partial"
IMPUTED_EVIDENCE = "partial_guided_imputed"
PSEUDO_SCAFFOLD_EVIDENCE = "pseudo_only_scaffold"

POLICY_EVIDENCE = {
    "none": frozenset(),
    "partial": frozenset({PARTIAL_EVIDENCE, EXACT_ECTODOMAIN_EVIDENCE}),
    "imputed": frozenset(
        {PARTIAL_EVIDENCE, EXACT_ECTODOMAIN_EVIDENCE, IMPUTED_EVIDENCE}
    ),
    "pseudo-scaffold": frozenset(
        {
            PARTIAL_EVIDENCE,
            EXACT_ECTODOMAIN_EVIDENCE,
            IMPUTED_EVIDENCE,
            PSEUDO_SCAFFOLD_EVIDENCE,
        }
    ),
}


def normalize_fallback_policy(value: str | None) -> str:
    policy = str(value or "none").strip().lower().replace("_", "-")
    if policy not in FALLBACK_POLICIES:
        raise MHCTransformError(
            "fallback must be one of " + ", ".join(FALLBACK_POLICIES)
        )
    return policy


def evidence_allowed(policy: str, evidence_level: str) -> bool:
    return evidence_level in POLICY_EVIDENCE[normalize_fallback_policy(policy)]


def evidence_priority(policy: str, evidence_level: str) -> int:
    """Return a lower-is-better evidence rank for a requested policy."""

    policy = normalize_fallback_policy(policy)
    if policy == "partial":
        order = (EXACT_ECTODOMAIN_EVIDENCE, PARTIAL_EVIDENCE)
    elif policy in {"imputed", "pseudo-scaffold"}:
        order = (
            EXACT_ECTODOMAIN_EVIDENCE,
            IMPUTED_EVIDENCE,
            PSEUDO_SCAFFOLD_EVIDENCE,
            PARTIAL_EVIDENCE,
        )
    else:
        return 999
    try:
        return order.index(evidence_level)
    except ValueError:
        return 999


def hamming_distance(left: str, right: str) -> int:
    """Hamming distance while treating ``X`` as unspecified."""

    if len(left) != len(right):
        raise ValueError("pseudosequences must have equal length")
    return sum(a != b for a, b in zip(left, right) if a != "X" and b != "X")


def compress_positions(positions: Iterable[int]) -> str:
    ordered = sorted(set(int(position) for position in positions))
    if not ordered:
        return ""
    spans = []
    start = previous = ordered[0]
    for position in ordered[1:]:
        if position == previous + 1:
            previous = position
            continue
        spans.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = position
    spans.append(str(start) if start == previous else f"{start}-{previous}")
    return ";".join(spans)


def expand_positions(spans: str) -> Tuple[int, ...]:
    positions = []
    for item in str(spans).split(";"):
        item = item.strip()
        if not item:
            continue
        if "-" in item:
            start, end = (int(value) for value in item.split("-", 1))
            positions.extend(range(start, end + 1))
        else:
            positions.append(int(item))
    return tuple(sorted(set(positions)))


@lru_cache(maxsize=4096)
def global_alignment_map(reference: str, target: str) -> Mapping[int, int]:
    """Map 1-based reference residues onto a homologous target sequence."""

    n, m = len(reference), len(target)
    # MHC alleles are substitution-rich but short internal indels are rare.
    # A weak gap penalty can invent a compensating delete/insert pair merely
    # to avoid a few polymorphisms (for example around HLA-B positions 80-81),
    # shifting otherwise homologous contact coordinates.  Keep this in step
    # with ``local_alignment_map``: a real one-residue indel still wins over
    # hundreds of downstream mismatches, while substitution-only alleles stay
    # on the same coordinate frame.
    gap = -8
    scores = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        scores[i][0] = i * gap
        trace[i][0] = 1
    for j in range(1, m + 1):
        scores[0][j] = j * gap
        trace[0][j] = 2
    for i in range(1, n + 1):
        ref_residue = reference[i - 1]
        previous = scores[i - 1]
        current = scores[i]
        for j in range(1, m + 1):
            values = (
                previous[j - 1] + (2 if ref_residue == target[j - 1] else -1),
                previous[j] + gap,
                current[j - 1] + gap,
            )
            move = max(range(3), key=values.__getitem__)
            current[j] = values[move]
            trace[i][j] = move
    mapping: Dict[int, int] = {}
    i, j = n, m
    while i or j:
        move = trace[i][j]
        if move == 0:
            mapping[i] = j
            i -= 1
            j -= 1
        elif move == 1:
            i -= 1
        else:
            j -= 1
    return mapping


@dataclass(frozen=True)
class ContactProjection:
    """Auditable projection of pseudo contact sites onto a protein sequence."""

    positions: Tuple[int, ...]
    anchors: Tuple[int, ...]
    alignment_gap_pseudo_positions: Tuple[int, ...]
    shifted_pseudo_positions: Tuple[int, ...]
    coordinate_offsets: Tuple[int, ...]
    coordinate_method: str
    max_abs_shift: int


def _interpolated_alignment_position(
    reference_position: int,
    mapping: Mapping[int, int],
    target_length: int,
) -> int:
    """Interpolate a target coordinate when a reference contact maps to a gap."""

    left = max(
        (position for position in mapping if position < reference_position),
        default=None,
    )
    right = min(
        (position for position in mapping if position > reference_position),
        default=None,
    )
    if left is None and right is None:
        raise ValueError(
            f"cannot interpolate reference contact position {reference_position}"
        )
    estimates = []
    if left is not None:
        estimates.append(mapping[left] + reference_position - left)
    if right is not None:
        estimates.append(mapping[right] - (right - reference_position))
    # All coordinates are positive.  Integer division after adding half the
    # denominator implements deterministic round-half-up for two estimates.
    estimate = (sum(estimates) + len(estimates) // 2) // len(estimates)
    return max(1, min(estimate, target_length))


def project_contact_positions(
    reference_sequence: str,
    target_sequence: str,
    reference_positions: Sequence[int],
    target_pseudo_segment: str,
) -> ContactProjection:
    """Project contacts with bounded, explicitly audited residue recovery.

    Alignment is the coordinate authority.  A reference contact aligned to a
    gap is linearly interpolated from its nearest mapped neighbours.  A non-X
    pseudo residue may recover at most one position away, and only when the
    bounded candidate is unique.  ``X`` is a wildcard used to certify the
    coordinate; it is never inserted into a protein sequence.
    """

    if len(reference_positions) != len(target_pseudo_segment):
        raise ValueError("contact position and pseudosequence lengths differ")
    if len(target_sequence) < len(reference_positions):
        raise ValueError("target sequence is shorter than the contact set")

    mapping = global_alignment_map(reference_sequence, target_sequence)
    anchors = []
    gap_positions = []
    for pseudo_index, reference_position in enumerate(reference_positions, start=1):
        if reference_position in mapping:
            anchor = mapping[reference_position]
        else:
            anchor = _interpolated_alignment_position(
                reference_position, mapping, len(target_sequence)
            )
            gap_positions.append(pseudo_index)
        anchors.append(anchor)

    candidate_sets = []
    for index, (anchor, residue) in enumerate(
        zip(anchors, target_pseudo_segment), start=1
    ):
        lower = max(1, anchor - 1)
        upper = min(len(target_sequence), anchor + 1)
        if residue == "X":
            candidates = tuple(range(lower, upper + 1))
        elif target_sequence[anchor - 1] == residue:
            # The alignment anchor is authoritative when it already matches;
            # nearby repeated residues must not create artificial ambiguity.
            candidates = (anchor,)
        else:
            candidates = tuple(
                position
                for position in range(lower, upper + 1)
                if target_sequence[position - 1] == residue
            )
        if not candidates:
            raise ValueError(
                f"cannot locate pseudo residue {index}={residue} within 1 aa "
                f"of alignment anchor {anchor}"
            )
        candidate_sets.append(candidates)

    # Resolve repeated residues and wildcard coordinates as one global,
    # strictly increasing assignment.  This avoids a greedy choice at one
    # contact stealing the only valid coordinate from the next.  An assignment
    # is accepted only when its minimum-total-shift solution is unique.
    # state value: (cost, coordinate path, number of equal-cost paths capped at 2)
    states: Dict[int, Tuple[int, Tuple[int, ...], int]] = {0: (0, (), 1)}
    for index, (anchor, candidates) in enumerate(
        zip(anchors, candidate_sets), start=1
    ):
        next_states: Dict[int, Tuple[int, Tuple[int, ...], int]] = {}
        for previous, (cost, path, path_count) in states.items():
            for chosen in candidates:
                if chosen <= previous:
                    continue
                new_cost = cost + abs(chosen - anchor)
                new_path = (*path, chosen)
                existing = next_states.get(chosen)
                if existing is None or new_cost < existing[0]:
                    next_states[chosen] = (new_cost, new_path, path_count)
                elif new_cost == existing[0]:
                    next_states[chosen] = (
                        existing[0],
                        existing[1],
                        min(2, existing[2] + path_count),
                    )
        if not next_states:
            raise ValueError(
                f"cannot preserve contact order at pseudo position {index}"
            )
        states = next_states

    minimum_cost = min(value[0] for value in states.values())
    best = [value for value in states.values() if value[0] == minimum_cost]
    solution_count = sum(value[2] for value in best)
    if solution_count != 1:
        raise ValueError(
            "ambiguous bounded contact recovery: "
            f"{solution_count} minimum-shift assignments"
        )
    located = best[0][1]
    offsets = tuple(
        chosen - anchor for chosen, anchor in zip(located, anchors)
    )
    shifted_positions = tuple(
        index
        for index, offset in enumerate(offsets, start=1)
        if offset
    )

    observed = "".join(target_sequence[position - 1] for position in located)
    for expected_residue, observed_residue in zip(target_pseudo_segment, observed):
        if expected_residue != "X" and expected_residue != observed_residue:
            raise ValueError("projected positions do not reproduce the donor pseudo")

    if gap_positions and shifted_positions:
        method = "global_alignment_gap_interpolation_and_bounded_recovery"
    elif gap_positions:
        method = "global_alignment_gap_interpolation"
    elif shifted_positions:
        method = "global_alignment_bounded_recovery"
    else:
        method = "global_alignment"
    return ContactProjection(
        positions=tuple(located),
        anchors=tuple(anchors),
        alignment_gap_pseudo_positions=tuple(gap_positions),
        shifted_pseudo_positions=shifted_positions,
        coordinate_offsets=offsets,
        coordinate_method=method,
        max_abs_shift=max((abs(value) for value in offsets), default=0),
    )


def locate_contact_positions(
    reference_sequence: str,
    target_sequence: str,
    reference_positions: Sequence[int],
    target_pseudo_segment: str,
) -> Tuple[int, ...]:
    """Backward-compatible coordinate-only contact projection helper."""

    return project_contact_positions(
        reference_sequence,
        target_sequence,
        reference_positions,
        target_pseudo_segment,
    ).positions


@lru_cache(maxsize=4096)
def local_alignment_map(query: str, target: str) -> Mapping[int, int]:
    """Return a 1-based Smith-Waterman query-to-target residue map."""

    n, m = len(query), len(target)
    # MHC allele sequences are substitution-rich but rarely contain short
    # internal indels. A strong gap penalty prevents local alignment from
    # inventing a one-residue shift merely to avoid two polymorphisms.
    gap = -8
    scores = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[3] * (m + 1) for _ in range(n + 1)]
    best_score = best_i = best_j = 0
    for i in range(1, n + 1):
        query_residue = query[i - 1]
        previous = scores[i - 1]
        current = scores[i]
        for j in range(1, m + 1):
            diagonal = previous[j - 1] + (
                3 if query_residue == target[j - 1] else -1
            )
            values = (diagonal, previous[j] + gap, current[j - 1] + gap, 0)
            move = max(range(4), key=values.__getitem__)
            current[j] = values[move]
            trace[i][j] = move
            if current[j] > best_score:
                best_score, best_i, best_j = current[j], i, j

    mapping: Dict[int, int] = {}
    i, j = best_i, best_j
    while i and j and scores[i][j] > 0:
        move = trace[i][j]
        if move == 0:
            mapping[i] = j
            i -= 1
            j -= 1
        elif move == 1:
            i -= 1
        elif move == 2:
            j -= 1
        else:
            break
    return mapping


@lru_cache(maxsize=4096)
def apply_observed_sequence(
    scaffold: str,
    observed_sequence: str,
) -> Tuple[str, Tuple[int, ...], Tuple[str, ...]]:
    """Overlay a homologous deposited sequence onto a complete scaffold."""

    mapping = local_alignment_map(observed_sequence, scaffold)
    residues = list(scaffold)
    observed_positions = []
    for query_position, scaffold_position in mapping.items():
        residue = observed_sequence[query_position - 1]
        if residue == "X":
            continue
        residues[scaffold_position - 1] = residue
        observed_positions.append(scaffold_position)

    flags = []
    if scaffold:
        coverage = len(set(observed_positions)) / len(scaffold)
        if coverage < 0.5:
            flags.append("observed_target_coverage_below_50_percent")
        elif coverage < 0.9:
            flags.append("observed_target_coverage_below_90_percent")
    return "".join(residues), tuple(sorted(set(observed_positions))), tuple(flags)


def apply_pseudosequence(
    scaffold: str,
    contact_positions: Sequence[int],
    pseudo_segment: str,
    observed_positions: Iterable[int] = (),
) -> Tuple[str, Tuple[int, ...], Tuple[int, ...]]:
    """Apply non-X pseudo residues without overwriting conflicting observations."""

    if len(contact_positions) != len(pseudo_segment):
        raise ValueError("contact position and pseudosequence lengths differ")
    residues = list(scaffold)
    observed = set(observed_positions)
    unknown = []
    for pseudo_index, (position, residue) in enumerate(
        zip(contact_positions, pseudo_segment), start=1
    ):
        if residue == "X":
            unknown.append(pseudo_index)
            continue
        if position in observed and residues[position - 1] != residue:
            raise ValueError(
                f"deposited residue conflicts with pseudo at pseudo position {pseudo_index}"
            )
        residues[position - 1] = residue
        observed.add(position)
    return "".join(residues), tuple(sorted(observed)), tuple(unknown)


def first_field(name: str) -> str:
    match = re.search(r"\*([0-9]+)", str(name))
    return match.group(1) if match else ""
