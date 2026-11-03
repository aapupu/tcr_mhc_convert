"""Forward reconstruction of full alpha and beta TCR protein chains."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .errors import TCRTransformError
from .junctions import (
    JunctionMetadata,
    TerminalMatch,
    alpha_jc_bridge,
    get_junction_metadata,
    load_junction_reference,
    match_j_terminal,
    match_v_terminal,
)
from .references import (
    DEFAULT_JUNCTION_REFERENCE,
    choose_trbc_gene,
    clean_sequence,
    load_gene_sequences,
    normalize_gene,
    normalize_species,
)


@dataclass(frozen=True)
class ReconstructedChain:
    sequence: str
    species: str
    chain: str
    cdr3: str
    v_gene: str
    j_gene: str
    constant_gene: str
    v_match: TerminalMatch
    j_match: TerminalMatch
    jc_bridge: str
    v_reference_gene: str
    j_reference_gene: str
    constant_reference_gene: str
    reference_resolution: str


@dataclass(frozen=True)
class ReconstructedPair:
    alpha: ReconstructedChain
    beta: ReconstructedChain


def find_overlap(
    left: str,
    right: str,
    *,
    min_overlap: int = 1,
    max_overlap: int = 30,
) -> tuple[int, str]:
    """Return the longest literal suffix-prefix overlap."""

    maximum = min(len(left), len(right), max_overlap)
    for length in range(maximum, min_overlap - 1, -1):
        if left[-length:] == right[:length]:
            return length, left[-length:]
    return 0, ""


def append_constant_region(
    variable_chain: str,
    constant_sequence: str,
    bridge_residue: str = "",
) -> str:
    variable_chain = variable_chain + bridge_residue
    overlap_length, _ = find_overlap(variable_chain, constant_sequence)
    return variable_chain + constant_sequence[overlap_length:]


def _validate_chain_genes(chain: str, v_gene: str, j_gene: str) -> None:
    expected_v = "TRAV" if chain == "A" else "TRBV"
    expected_j = "TRAJ" if chain == "A" else "TRBJ"
    if not v_gene.split("*", 1)[0].startswith(expected_v):
        raise TCRTransformError(f"{v_gene!r} is not a {expected_v} gene")
    if not j_gene.split("*", 1)[0].startswith(expected_j):
        raise TCRTransformError(f"{j_gene!r} is not a {expected_j} gene")


def _resolve_sequence_with_key(
    gene_sequences: dict[str, str],
    gene_name: str,
) -> tuple[str, str, str]:
    if gene_name in gene_sequences:
        return gene_name, gene_sequences[gene_name], "exact"
    base = gene_name.split("*", 1)[0]
    if base in gene_sequences:
        mode = "gene_level_fallback" if base != gene_name else "exact"
        return base, gene_sequences[base], mode
    return "", "", "missing"


def _missing_boundary_message(
    *,
    side: str,
    gene: str,
    min_overlap: int,
    short_match: TerminalMatch | None,
) -> str:
    if short_match is not None:
        return (
            f"only a boundary-valid {short_match.length}-aa {side} match was found "
            f"for {gene}; it is below the configured min_overlap={min_overlap} QC "
            "threshold and is not classified as a biologically invalid TCR"
        )
    return (
        f"no boundary-valid {side} match was found for {gene}; guarded 1-2 aa "
        "matches require the annotated anchor, expected boundary, and a unique position"
    )


def reconstruct_tcr_chain(
    *,
    cdr3: object,
    v_gene: object,
    j_gene: object,
    chain: str,
    species: object = "human",
    constant_gene: object | None = None,
    gene_sequences: dict[str, str] | None = None,
    junctions: dict[str, JunctionMetadata] | None = None,
    gene_reference_path: str | Path | None = None,
    junction_reference_path: str | Path = DEFAULT_JUNCTION_REFERENCE,
    sequence_column: str = "amino_acid",
    min_overlap: int = 3,
    max_overlap: int = 30,
    v_terminal_window: int = 8,
    j_start_window: int = 10,
) -> ReconstructedChain:
    """Build one trimming-aware, constant-completed TCR chain.

    CDR3 includes the V-side conserved C and the J-side anchor. The default
    three-residue overlap is a reconstruction QC threshold. Setting it to one
    or two invokes the stricter anchor, boundary, and uniqueness safeguards in
    :mod:`tcrtransform.junctions`.
    """

    species = normalize_species(species)
    chain = str(chain).strip().upper()
    if chain not in {"A", "B"}:
        raise TCRTransformError(f"unsupported chain {chain!r}; expected A or B")
    if min_overlap < 1 or max_overlap < min_overlap:
        raise TCRTransformError("overlap limits must satisfy 1 <= min <= max")
    if v_terminal_window < 0 or j_start_window < 0:
        raise TCRTransformError("boundary windows must be non-negative")

    cdr3 = clean_sequence(cdr3, f"TCR {chain} CDR3")
    v_gene = normalize_gene(v_gene, f"TCR {chain} V gene")
    j_gene = normalize_gene(j_gene, f"TCR {chain} J gene")
    _validate_chain_genes(chain, v_gene, j_gene)
    if constant_gene is None:
        constant_gene = "TRAC" if chain == "A" else choose_trbc_gene(j_gene)
    constant_gene = normalize_gene(constant_gene, f"TCR {chain} constant gene")
    expected_constant = "TRAC" if chain == "A" else "TRBC"
    if not constant_gene.split("*", 1)[0].startswith(expected_constant):
        raise TCRTransformError(
            f"{constant_gene!r} is not a {expected_constant} constant gene"
        )

    genes = gene_sequences or load_gene_sequences(
        species, path=gene_reference_path, sequence_column=sequence_column
    )
    junction_records = junctions or load_junction_reference(
        species, path=junction_reference_path
    )
    v_reference_gene, v_sequence, v_resolution = _resolve_sequence_with_key(
        genes, v_gene
    )
    j_reference_gene, j_sequence, j_resolution = _resolve_sequence_with_key(
        genes, j_gene
    )
    constant_reference_gene, constant_sequence, c_resolution = _resolve_sequence_with_key(
        genes, constant_gene
    )
    missing = [
        name
        for name, sequence in (
            (v_gene, v_sequence),
            (j_gene, j_sequence),
            (constant_gene, constant_sequence),
        )
        if not sequence
    ]
    if missing:
        raise TCRTransformError(
            f"reference sequence missing for: {', '.join(missing)}"
        )

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
    if v_match is None:
        short_match = match_v_terminal(
            v_sequence,
            cdr3,
            min_overlap=1,
            max_overlap=min(2, max_overlap),
            terminal_window=v_terminal_window,
            junction=v_metadata,
        )
        raise TCRTransformError(
            _missing_boundary_message(
                side="V-CDR3", gene=v_gene, min_overlap=min_overlap,
                short_match=short_match,
            )
        )

    j_match = match_j_terminal(
        cdr3,
        j_sequence,
        min_overlap=min_overlap,
        max_overlap=max_overlap,
        start_window=j_start_window,
        junction=j_metadata,
    )
    if j_match is None:
        short_match = match_j_terminal(
            cdr3,
            j_sequence,
            min_overlap=1,
            max_overlap=min(2, max_overlap),
            start_window=j_start_window,
            junction=j_metadata,
        )
        raise TCRTransformError(
            _missing_boundary_message(
                side="CDR3-J", gene=j_gene, min_overlap=min_overlap,
                short_match=short_match,
            )
        )

    variable_chain = (
        v_sequence[: v_match.gene_end]
        + cdr3[v_match.length :]
        + j_sequence[j_match.gene_end :]
    )
    bridge = ""
    if chain == "A":
        bridge = alpha_jc_bridge(junction_records, j_gene)
        if not bridge:
            raise TCRTransformError(
                f"no reliable productive TRAJ-TRAC split-codon residue for {j_gene}"
            )
    sequence = append_constant_region(variable_chain, constant_sequence, bridge)
    if cdr3 not in sequence:
        raise TCRTransformError("internal reconstruction error: CDR3 was not retained")
    if not sequence.endswith(constant_sequence):
        raise TCRTransformError("internal reconstruction error: constant region mismatch")
    return ReconstructedChain(
        sequence=sequence,
        species=species,
        chain=chain,
        cdr3=cdr3,
        v_gene=v_gene,
        j_gene=j_gene,
        constant_gene=constant_gene,
        v_match=v_match,
        j_match=j_match,
        jc_bridge=bridge,
        v_reference_gene=v_reference_gene,
        j_reference_gene=j_reference_gene,
        constant_reference_gene=constant_reference_gene,
        reference_resolution=(
            f"V:{v_resolution}({v_reference_gene});"
            f"J:{j_resolution}({j_reference_gene});"
            f"C:{c_resolution}({constant_reference_gene})"
        ),
    )


def build_tcr_chain(**kwargs: object) -> str:
    """Convenience wrapper returning only the reconstructed sequence."""

    return reconstruct_tcr_chain(**kwargs).sequence


def build_paired_tcr(
    *,
    cdr3a: object,
    trav: object,
    traj: object,
    cdr3b: object,
    trbv: object,
    trbj: object,
    species: object = "human",
    **kwargs: object,
) -> ReconstructedPair:
    species = normalize_species(species)
    genes = kwargs.pop("gene_sequences", None)
    if genes is None:
        genes = load_gene_sequences(
            species,
            path=kwargs.get("gene_reference_path"),
            sequence_column=str(kwargs.get("sequence_column", "amino_acid")),
        )
    junctions = kwargs.pop("junctions", None)
    if junctions is None:
        junctions = load_junction_reference(
            species,
            path=kwargs.get("junction_reference_path", DEFAULT_JUNCTION_REFERENCE),
        )
    shared = dict(
        species=species,
        gene_sequences=genes,
        junctions=junctions,
        **kwargs,
    )
    alpha = reconstruct_tcr_chain(
        cdr3=cdr3a, v_gene=trav, j_gene=traj, chain="A", **shared
    )
    beta = reconstruct_tcr_chain(
        cdr3=cdr3b, v_gene=trbv, j_gene=trbj, chain="B", **shared
    )
    return ReconstructedPair(alpha=alpha, beta=beta)

