"""Trimming-aware V/CDR3 and CDR3/J protein boundary matching."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

from .errors import TCRReferenceError
from .references import (AA20, DEFAULT_JUNCTION_REFERENCE, normalize_species,
                         normalize_gene, imgt_records, imgt_representatives, reconstruction_candidate)


@dataclass(frozen=True)
class JunctionMetadata:
    species: str
    gene_name: str
    reference_allele: str
    region: str
    cdr3_anchor: str
    cdr3_anchor_index: int
    jc_bridge: str = ""
    jc_bridge_codon: str = ""


@dataclass(frozen=True)
class TerminalMatch:
    length: int
    sequence: str
    gene_start: int
    gene_end: int
    offset: int


def normalize_gene_name(value: object) -> str:
    return normalize_gene(value)


def load_junction_reference(
    species: object,
    path: str | Path = DEFAULT_JUNCTION_REFERENCE,
) -> dict[str, JunctionMetadata]:
    species = normalize_species(species)
    path = Path(path)
    records: dict[str, JunctionMetadata] = {}
    if path.resolve() == DEFAULT_JUNCTION_REFERENCE.resolve():
        rows = {normalize_gene(row["allele"]): row for row in imgt_records(species)
                if reconstruction_candidate(row)}
        rows.update(imgt_representatives(species))
        for name, row in rows.items():
            records[name] = JunctionMetadata(species, name, normalize_gene(row["allele"]),
                row["region"], row["cdr3_anchor"], int(row["cdr3_anchor_index"]),
                row["jc_bridge"], row["jc_bridge_codon"])
        return records
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "species", "gene_name", "reference_allele", "region",
            "cdr3_anchor", "cdr3_anchor_index", "jc_bridge", "jc_bridge_codon",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise TCRReferenceError(f"{path} is missing columns: {sorted(missing)}")
        for row_number, row in enumerate(reader, start=2):
            if normalize_species(row["species"]) != species:
                continue
            gene = normalize_gene_name(row["gene_name"])
            region = str(row["region"]).strip().upper()
            anchor = str(row["cdr3_anchor"]).strip().upper()
            bridge = str(row["jc_bridge"]).strip().upper()
            bridge_codon = str(row["jc_bridge_codon"]).strip().upper()
            try:
                anchor_index = int(row["cdr3_anchor_index"])
            except ValueError as exc:
                raise TCRReferenceError(
                    f"{path} row {row_number} has invalid cdr3_anchor_index"
                ) from exc
            if region not in {"V", "J"}:
                raise TCRReferenceError(f"{path} row {row_number} has invalid region")
            if len(anchor) != 1 or anchor not in AA20:
                raise TCRReferenceError(f"{path} row {row_number} has invalid anchor")
            if bridge and (len(bridge) != 1 or bridge not in AA20):
                raise TCRReferenceError(f"{path} row {row_number} has invalid bridge")
            if bridge_codon and not re.fullmatch(r"[ACGT]{3}", bridge_codon):
                raise TCRReferenceError(
                    f"{path} row {row_number} has invalid bridge codon"
                )
            record = JunctionMetadata(
                species=species,
                gene_name=gene,
                reference_allele=normalize_gene_name(row["reference_allele"]),
                region=region,
                cdr3_anchor=anchor,
                cdr3_anchor_index=anchor_index,
                jc_bridge=bridge,
                jc_bridge_codon=bridge_codon,
            )
            previous = records.setdefault(gene, record)
            if previous != record:
                raise TCRReferenceError(
                    f"conflicting junction records for {species} {gene}"
                )
    if not records:
        raise TCRReferenceError(f"{path} contains no junction records for {species}")
    return records


def get_junction_metadata(
    records: dict[str, JunctionMetadata],
    gene_name: object,
) -> JunctionMetadata | None:
    name = normalize_gene_name(gene_name)
    return records.get(name) or records.get(name.split("*", 1)[0])


def _find_v_candidates(
    v_sequence: str,
    query: str,
    terminal_window: int,
) -> list[tuple[int, int]]:
    candidates: list[tuple[int, int]] = []
    start = 0
    while True:
        position = v_sequence.find(query, start)
        if position < 0:
            break
        end = position + len(query)
        offset = len(v_sequence) - end
        if 0 <= offset <= terminal_window:
            candidates.append((position, offset))
        start = position + 1
    return candidates


def _find_j_candidates(j_sequence: str, query: str, start_window: int) -> list[int]:
    candidates: list[int] = []
    start = 0
    while True:
        position = j_sequence.find(query, start)
        if position < 0:
            break
        if 0 <= position <= start_window:
            candidates.append(position)
        start = position + 1
    return candidates


def _valid_short_v_match(
    v_sequence: str,
    cdr3: str,
    query: str,
    candidates: list[tuple[int, int]],
    junction: JunctionMetadata | None,
) -> bool:
    if not cdr3.startswith("C") or not query.startswith("C"):
        return False
    if junction is None or junction.region != "V" or junction.cdr3_anchor != "C":
        return False
    index = junction.cdr3_anchor_index
    return (
        0 <= index < len(v_sequence)
        and v_sequence[index] == "C"
        and len(candidates) == 1
        and candidates[0][0] == index
    )


def _valid_short_j_match(
    j_sequence: str,
    cdr3: str,
    query: str,
    candidates: list[int],
    junction: JunctionMetadata | None,
) -> bool:
    if junction is None or junction.region != "J":
        return False
    anchor = junction.cdr3_anchor
    index = junction.cdr3_anchor_index
    if not cdr3.endswith(anchor) or not query.endswith(anchor):
        return False
    return (
        0 <= index < len(j_sequence)
        and j_sequence[index] == anchor
        and len(candidates) == 1
        and candidates[0] + len(query) - 1 == index
    )


def match_v_terminal(
    v_sequence: str,
    cdr3: str,
    *,
    min_overlap: int = 3,
    max_overlap: int = 30,
    terminal_window: int = 8,
    junction: JunctionMetadata | None = None,
) -> TerminalMatch | None:
    """Return the longest trimming-aware V-supported CDR3 prefix.

    One- and two-residue matches are accepted only when the conserved C,
    annotated boundary, expected terminal window, and uniqueness all agree.
    """

    max_length = min(len(v_sequence), len(cdr3), max_overlap)
    for length in range(max_length, min_overlap - 1, -1):
        query = cdr3[:length]
        candidates = _find_v_candidates(v_sequence, query, terminal_window)
        if not candidates:
            continue
        if length <= 2 and not _valid_short_v_match(
            v_sequence, cdr3, query, candidates, junction
        ):
            continue
        position, offset = min(candidates, key=lambda item: (item[1], -item[0]))
        return TerminalMatch(length, query, position, position + length, offset)
    return None


def match_j_terminal(
    cdr3: str,
    j_sequence: str,
    *,
    min_overlap: int = 3,
    max_overlap: int = 30,
    start_window: int = 10,
    junction: JunctionMetadata | None = None,
) -> TerminalMatch | None:
    """Return the longest trimming-aware J-supported CDR3 suffix."""

    max_length = min(len(cdr3), len(j_sequence), max_overlap)
    for length in range(max_length, min_overlap - 1, -1):
        query = cdr3[-length:]
        candidates = _find_j_candidates(j_sequence, query, start_window)
        if not candidates:
            continue
        if length <= 2 and not _valid_short_j_match(
            j_sequence, cdr3, query, candidates, junction
        ):
            continue
        position = min(candidates)
        return TerminalMatch(length, query, position, position + length, position)
    return None


def alpha_jc_bridge(
    records: dict[str, JunctionMetadata],
    j_gene: object,
) -> str:
    record = get_junction_metadata(records, j_gene)
    if record is None or record.region != "J":
        return ""
    return record.jc_bridge

