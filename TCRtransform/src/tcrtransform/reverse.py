"""Infer V/J genes and CDR regions from a full TCR protein sequence."""

from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass

from .errors import TCRParseError
from .junctions import JunctionMetadata, get_junction_metadata, load_junction_reference
from .references import (
    VGeneLoopReference,
    choose_trbc_gene,
    clean_sequence,
    load_gene_sequences,
    normalize_species,
)
from .regions import CDR3Boundary, define_cdr3_boundary


@dataclass(frozen=True)
class VCandidate:
    gene: str
    identity: float
    evidence_length: int
    chain_start: int
    anchor_chain_index: int


@dataclass(frozen=True)
class JCandidate:
    gene: str
    identity: float
    evidence_length: int
    anchor_chain_index: int
    constant_gene: str


@dataclass(frozen=True)
class SpeciesReference:
    species: str
    gene_sequences: dict[str, str]
    junctions: dict[str, JunctionMetadata]


@dataclass(frozen=True)
class ChainAnnotation:
    species: str
    chain: str
    v_genes: str
    j_genes: str
    constant_genes: str
    cdr1: str
    cdr2: str
    cdr2_5: str
    cdr3: str
    v_identity: float
    j_identity: float
    cdr3_boundary: CDR3Boundary
    status: str

    def to_dict(self, prefix: str = "") -> dict[str, object]:
        values = asdict(self)
        values.pop("cdr3_boundary")
        result = {f"{prefix}{key}": value for key, value in values.items()}
        result.update(self.cdr3_boundary.to_dict(f"{prefix}cdr3_"))
        return result


def load_species_reference(species: object) -> SpeciesReference:
    species = normalize_species(species)
    return SpeciesReference(
        species=species,
        gene_sequences=load_gene_sequences(species),
        junctions=load_junction_reference(species),
    )


def _identity(left: str, right: str) -> float:
    if len(left) != len(right) or not left:
        return 0.0
    return sum(a == b for a, b in zip(left, right)) / len(left)


def _constant_start(chain: str, constant_sequence: str) -> int | None:
    if chain.endswith(constant_sequence):
        return len(chain) - len(constant_sequence)
    position = chain.find(constant_sequence)
    if position >= 0:
        return position
    probe = constant_sequence[: min(30, len(constant_sequence))]
    positions = [
        match.start() for match in re.finditer(f"(?={re.escape(probe)})", chain)
    ]
    return positions[0] if len(positions) == 1 else None


def _best_v_candidates(
    chain: str,
    chain_letter: str,
    reference: SpeciesReference,
    max_leader_length: int,
    min_identity: float,
) -> list[VCandidate]:
    prefix = "TRAV" if chain_letter == "A" else "TRBV"
    candidates: list[VCandidate] = []
    for gene, sequence in reference.gene_sequences.items():
        if not gene.startswith(prefix):
            continue
        junction = get_junction_metadata(reference.junctions, gene)
        if junction is None or junction.region != "V" or junction.cdr3_anchor != "C":
            continue
        anchor_index = junction.cdr3_anchor_index
        if not 0 <= anchor_index < len(sequence) or sequence[anchor_index] != "C":
            continue
        germline_prefix = sequence[: anchor_index + 1]
        if len(germline_prefix) < 50 or len(chain) < len(germline_prefix):
            continue
        max_offset = min(max_leader_length, len(chain) - len(germline_prefix))
        best: VCandidate | None = None
        for offset in range(max_offset + 1):
            observed = chain[offset : offset + len(germline_prefix)]
            identity = _identity(germline_prefix, observed)
            if identity < min_identity or observed[-1] != "C":
                continue
            candidate = VCandidate(
                gene=gene,
                identity=identity,
                evidence_length=len(germline_prefix),
                chain_start=offset,
                anchor_chain_index=offset + anchor_index,
            )
            if best is None or (
                candidate.identity,
                candidate.evidence_length,
                -candidate.chain_start,
            ) > (best.identity, best.evidence_length, -best.chain_start):
                best = candidate
        if best is not None:
            candidates.append(best)
    if not candidates:
        raise TCRParseError(
            f"no {prefix} reference reaches its CDR3 C with identity >= {min_identity:.3f}"
        )
    best_score = max((item.identity, item.evidence_length) for item in candidates)
    return sorted(
        [
            item
            for item in candidates
            if math.isclose(item.identity, best_score[0], abs_tol=1e-12)
            and item.evidence_length == best_score[1]
        ],
        key=lambda item: item.gene,
    )


def _best_j_candidates(
    chain: str,
    chain_letter: str,
    reference: SpeciesReference,
    min_identity: float,
) -> list[JCandidate]:
    prefix = "TRAJ" if chain_letter == "A" else "TRBJ"
    candidates: list[JCandidate] = []
    for gene, sequence in reference.gene_sequences.items():
        if not gene.startswith(prefix):
            continue
        junction = get_junction_metadata(reference.junctions, gene)
        if junction is None or junction.region != "J":
            continue
        anchor_index = junction.cdr3_anchor_index
        if not 0 <= anchor_index < len(sequence):
            continue
        constant_gene = "TRAC" if chain_letter == "A" else choose_trbc_gene(gene)
        constant_sequence = reference.gene_sequences.get(constant_gene, "")
        if not constant_sequence:
            continue
        bridge = junction.jc_bridge if chain_letter == "A" else ""
        if chain_letter == "A" and not bridge:
            continue
        constant_start = _constant_start(chain, constant_sequence)
        if constant_start is None:
            continue
        germline_tail = sequence[anchor_index:] + bridge
        tail_start = constant_start - len(germline_tail)
        if tail_start < 0:
            continue
        identity = _identity(germline_tail, chain[tail_start:constant_start])
        if identity < min_identity:
            continue
        upstream_overlap = 1
        max_upstream = min(anchor_index + 1, tail_start + 1, 30)
        for length in range(max_upstream, 0, -1):
            observed = chain[tail_start - length + 1 : tail_start + 1]
            germline = sequence[anchor_index - length + 1 : anchor_index + 1]
            if observed == germline:
                upstream_overlap = length
                break
        candidates.append(
            JCandidate(
                gene=gene,
                identity=identity,
                evidence_length=len(germline_tail) + upstream_overlap - 1,
                anchor_chain_index=tail_start,
                constant_gene=constant_gene,
            )
        )
    if not candidates:
        raise TCRParseError(
            f"no {prefix} anchor-to-constant tail matches with identity >= {min_identity:.3f}"
        )
    best_score = max((item.identity, item.evidence_length) for item in candidates)
    return sorted(
        [
            item
            for item in candidates
            if math.isclose(item.identity, best_score[0], abs_tol=1e-12)
            and item.evidence_length == best_score[1]
        ],
        key=lambda item: item.gene,
    )


def _observed_v_loops(
    chain: str,
    chain_letter: str,
    candidates: list[VCandidate],
    reference: SpeciesReference,
    loop_reference: VGeneLoopReference,
) -> tuple[str, str, str]:
    observed_sets = [set(), set(), set()]
    for candidate in candidates:
        try:
            loops = loop_reference.resolve(
                reference.species, chain_letter, candidate.gene
            )
        except ValueError:
            continue
        gene_sequence = reference.gene_sequences[candidate.gene]
        for target, loop in zip(
            observed_sets, (loops.cdr1, loops.cdr2, loops.cdr2_5)
        ):
            position = gene_sequence.find(loop)
            if position < 0:
                continue
            start = candidate.chain_start + position
            observed = chain[start : start + len(loop)]
            if len(observed) == len(loop):
                target.add(observed)
    if any(len(values) != 1 for values in observed_sets):
        raise TCRParseError(
            "best V candidate(s) do not yield one consensus CDR1/CDR2/CDR2.5"
        )
    return tuple(next(iter(values)) for values in observed_sets)  # type: ignore[return-value]


def _consensus_boundary(
    cdr3: str,
    chain_letter: str,
    v_candidates: list[VCandidate],
    j_candidates: list[JCandidate],
    reference: SpeciesReference,
) -> CDR3Boundary:
    boundaries = [
        define_cdr3_boundary(
            cdr3=cdr3,
            v_gene=v.gene,
            j_gene=j.gene,
            chain=chain_letter,
            species=reference.species,
            gene_sequences=reference.gene_sequences,
            junctions=reference.junctions,
        )
        for v in v_candidates
        for j in j_candidates
    ]
    signatures = {
        (
            boundary.v_supported_prefix,
            boundary.j_supported_suffix,
            boundary.non_vj_region,
            boundary.fixed_mask,
            boundary.status,
        )
        for boundary in boundaries
    }
    if len(signatures) != 1:
        raise TCRParseError(
            "ambiguous V/J candidates do not yield one consensus CDR3 partition"
        )
    return boundaries[0]


def parse_chain(
    chain_sequence: object,
    chain: str,
    reference: SpeciesReference | None = None,
    loop_reference: VGeneLoopReference | None = None,
    *,
    species: object = "human",
    max_leader_length: int = 60,
    min_v_identity: float = 0.95,
    min_j_identity: float = 0.90,
) -> ChainAnnotation:
    chain = str(chain).strip().upper()
    if chain not in {"A", "B"}:
        raise TCRParseError(f"unsupported chain {chain!r}; expected A or B")
    if max_leader_length < 0:
        raise TCRParseError("max_leader_length must be non-negative")
    if not 0 < min_v_identity <= 1 or not 0 < min_j_identity <= 1:
        raise TCRParseError("identity thresholds must be in (0, 1]")
    reference = reference or load_species_reference(species)
    loops = loop_reference or VGeneLoopReference.from_tsv()
    sequence = clean_sequence(chain_sequence, f"TCR {chain} full sequence")
    v_candidates = _best_v_candidates(
        sequence, chain, reference, max_leader_length, min_v_identity
    )
    j_candidates = _best_j_candidates(sequence, chain, reference, min_j_identity)

    cdr3s: set[str] = set()
    for v_candidate in v_candidates:
        for j_candidate in j_candidates:
            start = v_candidate.anchor_chain_index
            end = j_candidate.anchor_chain_index
            if start > end:
                continue
            cdr3 = sequence[start : end + 1]
            junction = get_junction_metadata(reference.junctions, j_candidate.gene)
            if (
                junction is not None
                and cdr3.startswith("C")
                and cdr3.endswith(junction.cdr3_anchor)
            ):
                cdr3s.add(cdr3)
    if len(cdr3s) != 1:
        raise TCRParseError(
            "best V/J candidates do not define one inclusive-anchor CDR3"
        )
    cdr3 = next(iter(cdr3s))
    cdr1, cdr2, cdr2_5 = _observed_v_loops(
        sequence, chain, v_candidates, reference, loops
    )
    boundary = _consensus_boundary(
        cdr3, chain, v_candidates, j_candidates, reference
    )
    status_parts = []
    # Score all available alleles, but report gene-level ambiguity: synonymous
    # proteins and the canonical alias must not create spurious extra genes.
    v_genes = sorted({item.gene.split("*", 1)[0] for item in v_candidates})
    j_genes = sorted({item.gene.split("*", 1)[0] for item in j_candidates})
    if len(v_genes) > 1:
        status_parts.append("ambiguous_v")
    if len(j_genes) > 1:
        status_parts.append("ambiguous_j")
    status = "+".join(status_parts) if status_parts else "ok"
    return ChainAnnotation(
        species=reference.species,
        chain=chain,
        v_genes="|".join(v_genes),
        j_genes="|".join(j_genes),
        constant_genes="|".join(
            sorted({item.constant_gene for item in j_candidates})
        ),
        cdr1=cdr1,
        cdr2=cdr2,
        cdr2_5=cdr2_5,
        cdr3=cdr3,
        v_identity=max(item.identity for item in v_candidates),
        j_identity=max(item.identity for item in j_candidates),
        cdr3_boundary=boundary,
        status=status,
    )

