"""CDR loop annotation and evidence-based CDR3 V/J boundary definition."""

from __future__ import annotations

from dataclasses import asdict, dataclass

from .errors import TCRTransformError
from .junctions import (
    JunctionMetadata,
    TerminalMatch,
    get_junction_metadata,
    load_junction_reference,
    match_j_terminal,
    match_v_terminal,
)
from .references import (
    VGeneLoopReference,
    VGeneLoops,
    clean_sequence,
    get_gene_sequence,
    load_gene_sequences,
    normalize_gene,
    normalize_species,
)


@dataclass(frozen=True)
class CDR3Boundary:
    """Protein-level evidence for the V/J-supported parts of one CDR3.

    ``non_vj_region`` is the contiguous CDR3 interval not covered by the
    accepted V-prefix or J-suffix match. It is not proof of nucleotide-level
    N/P addition because convergent amino acids and codon phase are unresolved.
    """

    cdr3: str
    species: str
    chain: str
    v_gene: str
    j_gene: str
    v_supported_prefix: str
    j_supported_suffix: str
    non_vj_region: str
    keep_prefix: int
    keep_suffix: int
    non_vj_start: int
    non_vj_end: int
    non_vj_length: int
    fixed_mask: str
    v_gene_start: int
    v_gene_end: int
    v_offset_from_end: int
    j_gene_start: int
    j_gene_end: int
    j_offset_from_start: int
    v_short_match_length: int
    j_short_match_length: int
    v_j_support_overlap: int
    safe_for_optimization: bool
    status: str

    def to_dict(self, prefix: str = "") -> dict[str, object]:
        return {f"{prefix}{key}": value for key, value in asdict(self).items()}


def annotate_v_gene(
    species: object,
    chain: str,
    v_gene: object,
    *,
    loop_reference: VGeneLoopReference | None = None,
) -> VGeneLoops:
    reference = loop_reference or VGeneLoopReference.from_tsv()
    return reference.resolve(species, chain, v_gene)


def _short_v_match(
    sequence: str,
    cdr3: str,
    window: int,
    metadata: JunctionMetadata | None,
) -> TerminalMatch | None:
    return match_v_terminal(
        sequence,
        cdr3,
        min_overlap=1,
        max_overlap=2,
        terminal_window=window,
        junction=metadata,
    )


def _short_j_match(
    cdr3: str,
    sequence: str,
    window: int,
    metadata: JunctionMetadata | None,
) -> TerminalMatch | None:
    return match_j_terminal(
        cdr3,
        sequence,
        min_overlap=1,
        max_overlap=2,
        start_window=window,
        junction=metadata,
    )


def define_cdr3_boundary(
    *,
    cdr3: object,
    v_gene: object,
    j_gene: object,
    chain: str,
    species: object = "human",
    gene_sequences: dict[str, str] | None = None,
    junctions: dict[str, JunctionMetadata] | None = None,
    min_overlap: int = 3,
    max_overlap: int = 30,
    v_terminal_window: int = 8,
    j_start_window: int = 10,
) -> CDR3Boundary:
    """Locate V-supported, J-supported, and intervening CDR3 residues.

    The default ``min_overlap=3`` is deliberately strict. If it rejects an
    otherwise guarded 1-2 aa match, that shorter match is reported in the
    diagnostic fields but is not used as an accepted boundary. Set
    ``min_overlap`` to 1 or 2 to accept only those short matches that pass the
    anchor, expected-position, and uniqueness safeguards.
    """

    species = normalize_species(species)
    chain = str(chain).strip().upper()
    if chain not in {"A", "B"}:
        raise TCRTransformError(f"unsupported chain {chain!r}; expected A or B")
    if min_overlap < 1 or max_overlap < min_overlap:
        raise TCRTransformError("overlap limits must satisfy 1 <= min <= max")
    cdr3 = clean_sequence(cdr3, f"TCR {chain} CDR3")
    v_gene = normalize_gene(v_gene, f"TCR {chain} V gene")
    j_gene = normalize_gene(j_gene, f"TCR {chain} J gene")
    expected_v = "TRAV" if chain == "A" else "TRBV"
    expected_j = "TRAJ" if chain == "A" else "TRBJ"
    if not v_gene.split("*", 1)[0].startswith(expected_v):
        raise TCRTransformError(f"{v_gene!r} is not a {expected_v} gene")
    if not j_gene.split("*", 1)[0].startswith(expected_j):
        raise TCRTransformError(f"{j_gene!r} is not a {expected_j} gene")

    genes = gene_sequences or load_gene_sequences(species)
    junction_records = junctions or load_junction_reference(species)
    v_sequence = get_gene_sequence(genes, v_gene)
    j_sequence = get_gene_sequence(genes, j_gene)
    if not v_sequence or not j_sequence:
        missing = v_gene if not v_sequence else j_gene
        raise TCRTransformError(f"reference sequence missing for {missing}")
    v_metadata = get_junction_metadata(junction_records, v_gene)
    j_metadata = get_junction_metadata(junction_records, j_gene)

    v_match = match_v_terminal(
        v_sequence,
        cdr3,
        min_overlap=min_overlap,
        max_overlap=max_overlap,
        terminal_window=v_terminal_window,
        junction=v_metadata,
    )
    j_match = match_j_terminal(
        cdr3,
        j_sequence,
        min_overlap=min_overlap,
        max_overlap=max_overlap,
        start_window=j_start_window,
        junction=j_metadata,
    )
    v_short = None if v_match is not None or min_overlap <= 2 else _short_v_match(
        v_sequence, cdr3, v_terminal_window, v_metadata
    )
    j_short = None if j_match is not None or min_overlap <= 2 else _short_j_match(
        cdr3, j_sequence, j_start_window, j_metadata
    )

    keep_prefix = v_match.length if v_match else 0
    keep_suffix = j_match.length if j_match else 0
    mask = [False] * len(cdr3)
    for index in range(min(keep_prefix, len(cdr3))):
        mask[index] = True
    for index in range(max(0, len(cdr3) - keep_suffix), len(cdr3)):
        mask[index] = True
    editable = [index for index, fixed in enumerate(mask) if not fixed]
    if editable:
        non_start = editable[0]
        non_end = editable[-1] + 1
        non_vj = cdr3[non_start:non_end]
    else:
        non_start = min(keep_prefix, len(cdr3))
        non_end = non_start
        non_vj = ""
    overlap = max(0, keep_prefix + keep_suffix - len(cdr3))

    status_parts: list[str] = []
    if v_match is None:
        status_parts.append(
            "v_below_min_overlap" if v_short is not None else "v_unresolved"
        )
    if j_match is None:
        status_parts.append(
            "j_below_min_overlap" if j_short is not None else "j_unresolved"
        )
    if not status_parts:
        if overlap:
            status_parts.append("overlapping_v_j_support")
        elif not non_vj:
            status_parts.append("fully_vj_supported")
        else:
            status_parts.append("ok")
    safe = v_match is not None and j_match is not None

    def coordinate(match: TerminalMatch | None, field: str) -> int:
        return int(getattr(match, field)) if match is not None else -1

    return CDR3Boundary(
        cdr3=cdr3,
        species=species,
        chain=chain,
        v_gene=v_gene,
        j_gene=j_gene,
        v_supported_prefix=v_match.sequence if v_match else "",
        j_supported_suffix=j_match.sequence if j_match else "",
        non_vj_region=non_vj,
        keep_prefix=keep_prefix,
        keep_suffix=keep_suffix,
        non_vj_start=non_start,
        non_vj_end=non_end,
        non_vj_length=len(non_vj),
        fixed_mask="".join("1" if fixed else "0" for fixed in mask),
        v_gene_start=coordinate(v_match, "gene_start"),
        v_gene_end=coordinate(v_match, "gene_end"),
        v_offset_from_end=coordinate(v_match, "offset"),
        j_gene_start=coordinate(j_match, "gene_start"),
        j_gene_end=coordinate(j_match, "gene_end"),
        j_offset_from_start=coordinate(j_match, "offset"),
        v_short_match_length=v_short.length if v_short else 0,
        j_short_match_length=j_short.length if j_short else 0,
        v_j_support_overlap=overlap,
        safe_for_optimization=safe,
        status="+".join(status_parts),
    )

