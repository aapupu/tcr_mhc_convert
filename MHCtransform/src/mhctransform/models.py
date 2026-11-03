"""Immutable public result objects used by :mod:`mhctransform`."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Tuple


@dataclass(frozen=True)
class SequenceCandidate:
    record_id: str
    allele: str
    aliases: Tuple[str, ...]
    species: str
    mhc_class: str
    chain_role: str
    locus: str
    sequence_kind: str
    sequence: str
    source: str
    source_version: str
    source_accession: str
    is_complete: bool


@dataclass(frozen=True)
class PseudoCandidate:
    record_id: str
    canonical_name: str
    aliases: Tuple[str, ...]
    species: str
    mhc_class: str
    definition_id: str
    pseudo_sequence: str
    alpha_allele: str
    beta_allele: str
    source: str
    source_version: str


@dataclass(frozen=True)
class AF3SequenceCandidate:
    record_id: str
    allele: str
    aliases: Tuple[str, ...]
    species: str
    mhc_class: str
    chain_role: str
    locus: str
    definition_id: str
    sequence: str
    source_full_record_id: str
    precursor_start_1based: int
    precursor_end_1based: int
    boundary_method: str
    boundary_source: str
    source_accession: str


@dataclass(frozen=True)
class MoleculeCandidate:
    canonical_name: str
    species: str
    mhc_class: str
    chains: Tuple[SequenceCandidate, ...]


@dataclass(frozen=True)
class AF3MoleculeCandidate:
    canonical_name: str
    species: str
    mhc_class: str
    chains: Tuple[AF3SequenceCandidate, ...]


@dataclass(frozen=True)
class QueryResult:
    """A non-lossy lookup result.

    ``status`` is one of ``exact_unique``, ``exact_ambiguous``,
    ``projected_unique``, ``projected_ambiguous`` or ``not_found``.
    Ambiguity is deliberate: neither a two-field allele name nor a 34-residue
    pseudosequence is guaranteed to identify exactly one full protein/name.
    Reverse full-sequence results may compact many deposited allele names into
    one sequence candidate; those names remain visible in ``aliases``.
    """

    query: str
    representation: str
    status: str
    candidates: Tuple[Any, ...] = ()
    notes: Tuple[str, ...] = ()

    @property
    def candidate_count(self) -> int:
        return len(self.candidates)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "query": self.query,
            "representation": self.representation,
            "status": self.status,
            "candidate_count": self.candidate_count,
            "candidates": [asdict(candidate) for candidate in self.candidates],
            "notes": list(self.notes),
        }


def result_status(count: int, *, projected: bool = False) -> str:
    if count == 0:
        return "not_found"
    if count == 1:
        return "projected_unique" if projected else "exact_unique"
    return "projected_ambiguous" if projected else "exact_ambiguous"

