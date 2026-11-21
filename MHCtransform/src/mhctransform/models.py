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
    evidence_level: str
    official_support: bool
    validation_status: str
    source_record_ids: Tuple[str, ...]
    quality_flags: Tuple[str, ...]


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
class PBRSequenceCandidate:
    """A contiguous peptide-binding structural domain selection.

    MHC-I contains the heavy-chain alpha1+alpha2 domains. MHC-II uses one
    candidate for alpha1 and another for beta1; these chains are never
    concatenated.
    """

    record_id: str
    allele: str
    aliases: Tuple[str, ...]
    species: str
    mhc_class: str
    chain_role: str
    locus: str
    definition_id: str
    sequence: str
    source_af3_record_id: str
    af3_start_1based: int
    af3_end_1based: int
    source_full_record_id: str
    precursor_start_1based: int
    precursor_end_1based: int
    boundary_method: str
    boundary_source: str
    source_accession: str


@dataclass(frozen=True)
class FallbackSequenceCandidate:
    """A deposited fragment or explicitly synthetic precursor construct.

    These candidates intentionally live outside :class:`SequenceCandidate` so
    a partial or imputed construct can never enter the strict full-sequence
    index by accident.
    """

    record_id: str
    allele: str
    aliases: Tuple[str, ...]
    species: str
    mhc_class: str
    chain_role: str
    locus: str
    representation: str
    sequence: str
    evidence_level: str
    observed_spans_1based: str
    imputed_spans_1based: str
    observed_residue_count: int
    imputed_residue_count: int
    source_record_ids: Tuple[str, ...]
    template_record_id: str
    template_allele: str
    target_pseudo_record_id: str
    pseudo_distance: int
    unknown_pseudo_positions: Tuple[int, ...]
    quality_flags: Tuple[str, ...]


@dataclass(frozen=True)
class FallbackAF3SequenceCandidate:
    """An exact ectodomain or auditable AF3 fallback construct."""

    record_id: str
    allele: str
    aliases: Tuple[str, ...]
    species: str
    mhc_class: str
    chain_role: str
    locus: str
    definition_id: str
    sequence: str
    evidence_level: str
    observed_spans_1based: str
    imputed_spans_1based: str
    observed_residue_count: int
    imputed_residue_count: int
    source_record_ids: Tuple[str, ...]
    template_record_id: str
    template_allele: str
    target_pseudo_record_id: str
    pseudo_distance: int
    unknown_pseudo_positions: Tuple[int, ...]
    quality_flags: Tuple[str, ...]


@dataclass(frozen=True)
class FallbackPBRSequenceCandidate:
    """An opt-in non-strict PBR representation.

    ``sequence_representation`` separates a template-completed construct from
    a partial chain whose genuinely unknown residues are encoded as ``X``.
    ``unknown_*`` is deliberately distinct from ``imputed_*``: an unknown
    residue has not been asserted to equal a template residue. For a masked
    candidate, ``source_af3_record_id`` identifies the chain used to establish
    the PBR coordinate frame only; the masked sequence is not asserted to be a
    literal AF3-chain substring. ``masked_from_record_id`` links the separate
    template-completed PBR record without reusing its identity.
    """

    record_id: str
    allele: str
    aliases: Tuple[str, ...]
    species: str
    mhc_class: str
    chain_role: str
    locus: str
    definition_id: str
    sequence_representation: str
    sequence: str
    source_af3_record_id: str
    af3_start_1based: int
    af3_end_1based: int
    evidence_level: str
    observed_spans_1based: str
    unknown_spans_1based: str
    imputed_spans_1based: str
    observed_residue_count: int
    unknown_residue_count: int
    imputed_residue_count: int
    source_record_ids: Tuple[str, ...]
    masked_from_record_id: str
    template_record_id: str
    template_allele: str
    target_pseudo_record_id: str
    pseudo_distance: int
    unknown_pseudo_positions: Tuple[int, ...]
    quality_flags: Tuple[str, ...]


@dataclass(frozen=True)
class MoleculeCandidate:
    canonical_name: str
    species: str
    mhc_class: str
    chains: Tuple[SequenceCandidate | FallbackSequenceCandidate, ...]


@dataclass(frozen=True)
class AF3MoleculeCandidate:
    canonical_name: str
    species: str
    mhc_class: str
    chains: Tuple[AF3SequenceCandidate | FallbackAF3SequenceCandidate, ...]


@dataclass(frozen=True)
class PBRMoleculeCandidate:
    canonical_name: str
    species: str
    mhc_class: str
    chains: Tuple[PBRSequenceCandidate | FallbackPBRSequenceCandidate, ...]


@dataclass(frozen=True)
class QueryResult:
    """A non-lossy lookup result.

    ``status`` is one of ``exact_unique``, ``exact_ambiguous``,
    ``projected_unique``, ``projected_ambiguous``, ``compatible_unique``,
    ``compatible_ambiguous`` or ``not_found``. ``compatible_*`` is reserved
    for reverse matches to partial or engineered fallback representations;
    it must not be interpreted as proof of a biological allele sequence.
    For forward fallback queries, ``exact_*`` describes exact name resolution;
    ``evidence`` and the candidate's ``evidence_level`` describe sequence
    provenance.
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
    fallback_used: str = "none"
    evidence: str = "source_exact"

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
            "fallback_used": self.fallback_used,
            "evidence": self.evidence,
        }


def result_status(count: int, *, projected: bool = False) -> str:
    if count == 0:
        return "not_found"
    if count == 1:
        return "projected_unique" if projected else "exact_unique"
    return "projected_ambiguous" if projected else "exact_ambiguous"

