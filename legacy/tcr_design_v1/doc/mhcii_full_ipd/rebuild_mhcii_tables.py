#!/usr/bin/env python3
"""Rebuild structure-ready MHC-II chain tables from raw IPD FASTA files.

The model sequence is the mature, membrane-free alpha1/alpha2 or beta1/beta2
core. Boundaries are transferred from a complete human reference by local
protein alignment. Partial records are rejected rather than padded or guessed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import pandas as pd

try:
    from Bio.Align import PairwiseAligner
except ImportError as exc:  # pragma: no cover - environment-specific message
    raise SystemExit(
        "Biopython is required to rebuild the MHC-II tables: "
        "python -m pip install biopython pandas"
    ) from exc


SCRIPT_VERSION = "1.0"
AA20 = frozenset("ACDEFGHIKLMNPQRSTVWY")
HYDROPHOBIC = frozenset("AILMFWVY")
EXPRESSION_SUFFIXES = frozenset("NLSCAQ")

MIN_REFERENCE_COVERAGE = 0.97
MIN_ALIGNMENT_IDENTITY = 0.72
MAX_UNALIGNED_REFERENCE_N = 2
MAX_UNALIGNED_REFERENCE_C = 3
MAX_MATURE_LENGTH_DEVIATION = 5
HYDROPHOBIC_WINDOW = 19
MAX_HYDROPHOBIC_IN_WINDOW = 13

SPECIES_BY_PREFIX = {
    "HLA": "human",
    "Mamu": "rhesus_macaque",
    "Patr": "chimpanzee",
    "Gogo": "gorilla",
    "BoLA": "cattle",
    "SLA": "pig",
    "DLA": "dog",
    "Eqca": "horse",
}

# Keep classical DP/DQ loci and the DR loci represented in the current project.
# DQA2/DQB2/DPB2 are intentionally outside this modeling reference set.
MHCII_GENE_RE = re.compile(
    r"^(?:DPA1|DPB1|DQA1|DQB1|DQA|DQB|DRA|DRB\d*)$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class ReferenceSlice:
    allele: str
    start_0based: int
    end_exclusive: int


# The end points stop at the conserved membrane-proximal domain terminus and
# omit signal peptide, connecting stalk, transmembrane helix, and cytoplasmic
# tail. Coordinates are zero-based, half-open in the raw precursor sequence.
REFERENCE_SLICES = {
    ("DR", "alpha"): ReferenceSlice("HLA-DRA*01:01:01:01", 25, 207),
    ("DR", "beta"): ReferenceSlice("HLA-DRB1*01:01:01:01", 29, 220),
    ("DQ", "alpha"): ReferenceSlice("HLA-DQA1*01:01:01:01", 23, 207),
    ("DQ", "beta"): ReferenceSlice("HLA-DQB1*05:01:01:01", 32, 223),
    ("DP", "alpha"): ReferenceSlice("HLA-DPA1*01:03:01:01", 31, 211),
    ("DP", "beta"): ReferenceSlice("HLA-DPB1*01:01:01:01", 29, 217),
}

SOURCE_FILES = (
    ("IPD-IMGT-HLA_hla_prot.fasta", "IPD-IMGT/HLA"),
    ("IPD-MHC_MHC_prot.fasta", "IPD-MHC"),
)

INTEGER_COLUMNS = (
    "source_index",
    "length",
    "declared_length",
    "reference_length",
    "aligned_residues",
    "unaligned_reference_n",
    "unaligned_reference_c",
    "alignment_target_start_0based",
    "alignment_target_end_exclusive",
    "mature_start_0based",
    "mature_end_exclusive",
    "mature_length",
    "n_terminal_trim",
    "c_terminal_trim",
    "max_hydrophobic_residues_19aa",
    "two_field_sequence_count",
    "source_allele_count",
    "source_record_count",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory containing raw/ (default: directory of this script).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: --data-dir).",
    )
    return parser.parse_args()


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    sequence_parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts).upper()
                header = line[1:].strip()
                sequence_parts = []
            else:
                if header is None:
                    raise ValueError(f"{path}:{line_number}: sequence before FASTA header")
                sequence_parts.append(line)
    if header is not None:
        yield header, "".join(sequence_parts).upper()


def parse_header_length(header: str) -> int | None:
    fields = header.split()
    if len(fields) >= 3 and fields[2].isdigit():
        return int(fields[2])
    return None


def classify_gene(gene: str) -> tuple[str, str]:
    gene_upper = gene.upper()
    family = gene_upper[:2]
    chain = "alpha" if gene_upper.startswith(("DPA", "DQA", "DRA")) else "beta"
    return family, chain


def read_source_records(data_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    source_index = 0
    for filename, database in SOURCE_FILES:
        path = data_dir / "raw" / filename
        if not path.is_file():
            raise FileNotFoundError(path)
        for header, sequence in iter_fasta(path):
            fields = header.split()
            if len(fields) < 2:
                raise ValueError(f"Malformed FASTA header in {path}: {header!r}")
            source_accession, allele_token = fields[:2]
            if database == "IPD-IMGT/HLA":
                species_prefix = "HLA"
                name = f"HLA-{allele_token}"
            else:
                if "-" not in allele_token:
                    continue
                species_prefix = allele_token.split("-", 1)[0]
                name = allele_token

            if species_prefix not in SPECIES_BY_PREFIX or "-" not in name:
                continue
            gene = name.split("-", 1)[1].split("*", 1)[0]
            if not MHCII_GENE_RE.fullmatch(gene):
                continue

            family, chain = classify_gene(gene)
            declared_length = parse_header_length(header)
            invalid_chars = "".join(sorted(set(sequence) - AA20))
            rows.append(
                {
                    "source_index": source_index,
                    "name": name,
                    "gene": gene,
                    "family": family,
                    "chain": chain,
                    "species": SPECIES_BY_PREFIX[species_prefix],
                    "species_prefix": species_prefix,
                    "database": database,
                    "source_accession": source_accession,
                    "length": len(sequence),
                    "declared_length": declared_length,
                    "length_matches_header": (
                        declared_length is None or declared_length == len(sequence)
                    ),
                    "invalid_chars": invalid_chars,
                    "sequence": sequence,
                    "header": header,
                }
            )
            source_index += 1
    if not rows:
        raise ValueError("No supported MHC-II records were found in the raw FASTA files")
    return pd.DataFrame(rows)


def expression_suffix(name: str) -> str:
    match = re.search(r"\d([A-Z])$", name.upper())
    if match and match.group(1) in EXPRESSION_SUFFIXES:
        return match.group(1)
    return ""


def two_field_name(name: str) -> str:
    if "*" not in name:
        return name
    prefix, allele = name.split("*", 1)
    fields = allele.split(":")
    if len(fields) < 2:
        return name
    return f"{prefix}*{fields[0]}:{fields[1]}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def max_hydrophobic_count(sequence: str, window: int = HYDROPHOBIC_WINDOW) -> int:
    if not sequence:
        return 0
    if len(sequence) <= window:
        return sum(residue in HYDROPHOBIC for residue in sequence)
    return max(
        sum(residue in HYDROPHOBIC for residue in sequence[start : start + window])
        for start in range(len(sequence) - window + 1)
    )


def build_reference_sequences(
    source: pd.DataFrame,
) -> dict[tuple[str, str], tuple[ReferenceSlice, str]]:
    references: dict[tuple[str, str], tuple[ReferenceSlice, str]] = {}
    for key, reference_slice in REFERENCE_SLICES.items():
        matches = source.loc[source["name"].eq(reference_slice.allele), "sequence"]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one raw reference record for {reference_slice.allele}, "
                f"found {len(matches)}"
            )
        precursor = str(matches.iloc[0])
        if not (0 <= reference_slice.start_0based < reference_slice.end_exclusive <= len(precursor)):
            raise ValueError(f"Invalid reference slice for {reference_slice.allele}")
        references[key] = (
            reference_slice,
            precursor[reference_slice.start_0based : reference_slice.end_exclusive],
        )
    return references


def new_aligner() -> PairwiseAligner:
    return PairwiseAligner(
        mode="local",
        match_score=2.0,
        mismatch_score=-1.0,
        open_gap_score=-5.0,
        extend_gap_score=-0.5,
    )


def assess_record(
    row: object,
    references: dict[tuple[str, str], tuple[ReferenceSlice, str]],
    aligner: PairwiseAligner,
) -> dict[str, object]:
    sequence = str(row.sequence)
    suffix = expression_suffix(str(row.name))
    result: dict[str, object] = {
        "expression_suffix": suffix,
        "reference_allele": "",
        "reference_length": None,
        "alignment_identity": None,
        "reference_coverage": None,
        "aligned_residues": None,
        "unaligned_reference_n": None,
        "unaligned_reference_c": None,
        "alignment_target_start_0based": None,
        "alignment_target_end_exclusive": None,
        "mature_start_0based": None,
        "mature_end_exclusive": None,
        "mature_length": None,
        "n_terminal_trim": None,
        "c_terminal_trim": None,
        "max_hydrophobic_residues_19aa": None,
        "full_aa": "",
        "model_ready": False,
        "qc_reason": "",
    }
    reasons: list[str] = []
    if str(row.invalid_chars):
        reasons.append("invalid_amino_acid")
    if not bool(row.length_matches_header):
        reasons.append("header_length_mismatch")
    if suffix:
        reasons.append(f"expression_suffix_{suffix}")

    key = (str(row.family), str(row.chain))
    reference_slice, reference = references[key]
    result["reference_allele"] = reference_slice.allele
    result["reference_length"] = len(reference)

    if str(row.invalid_chars) or not sequence:
        result["qc_reason"] = ";".join(reasons or ["empty_sequence"])
        return result

    alignments = aligner.align(reference, sequence)
    if len(alignments) == 0:
        reasons.append("no_reference_alignment")
        result["qc_reason"] = ";".join(reasons)
        return result

    alignment = alignments[0]
    reference_blocks, target_blocks = alignment.aligned
    reference_start = int(reference_blocks[0][0])
    reference_end = int(reference_blocks[-1][1])
    target_start = int(target_blocks[0][0])
    target_end = int(target_blocks[-1][1])

    aligned_residues = 0
    matches = 0
    for (ref_start, ref_end), (seq_start, seq_end) in zip(
        reference_blocks, target_blocks
    ):
        ref_start = int(ref_start)
        ref_end = int(ref_end)
        seq_start = int(seq_start)
        seq_end = int(seq_end)
        if ref_end - ref_start != seq_end - seq_start:
            raise AssertionError("Unexpected unequal gap-free alignment blocks")
        aligned_residues += ref_end - ref_start
        matches += sum(
            left == right
            for left, right in zip(
                reference[ref_start:ref_end], sequence[seq_start:seq_end]
            )
        )

    missing_n = reference_start
    missing_c = len(reference) - reference_end
    mature_start = target_start - missing_n
    mature_end = target_end + missing_c
    coverage = aligned_residues / len(reference)
    identity = matches / aligned_residues if aligned_residues else 0.0

    mature_sequence = ""
    if 0 <= mature_start < mature_end <= len(sequence):
        mature_sequence = sequence[mature_start:mature_end]
    mature_length = len(mature_sequence)
    hydrophobic_count = max_hydrophobic_count(mature_sequence)

    if missing_n > MAX_UNALIGNED_REFERENCE_N:
        reasons.append("reference_n_terminus_not_aligned")
    if missing_c > MAX_UNALIGNED_REFERENCE_C:
        reasons.append("reference_c_terminus_not_aligned")
    if mature_start < 0:
        reasons.append("incomplete_mature_n_terminus")
    if mature_end > len(sequence):
        reasons.append("incomplete_mature_c_terminus")
    if coverage < MIN_REFERENCE_COVERAGE:
        reasons.append("low_reference_coverage")
    if identity < MIN_ALIGNMENT_IDENTITY:
        reasons.append("low_alignment_identity")
    if not (
        len(reference) - MAX_MATURE_LENGTH_DEVIATION
        <= mature_length
        <= len(reference) + MAX_MATURE_LENGTH_DEVIATION
    ):
        reasons.append("unexpected_mature_length")
    if hydrophobic_count > MAX_HYDROPHOBIC_IN_WINDOW:
        reasons.append("transmembrane_like_segment_in_mature_sequence")

    result.update(
        {
            "alignment_identity": round(identity, 6),
            "reference_coverage": round(coverage, 6),
            "aligned_residues": aligned_residues,
            "unaligned_reference_n": missing_n,
            "unaligned_reference_c": missing_c,
            "alignment_target_start_0based": target_start,
            "alignment_target_end_exclusive": target_end,
            "mature_start_0based": mature_start,
            "mature_end_exclusive": mature_end,
            "mature_length": mature_length,
            "n_terminal_trim": mature_start if mature_start >= 0 else None,
            "c_terminal_trim": (
                len(sequence) - mature_end if mature_end <= len(sequence) else None
            ),
            "max_hydrophobic_residues_19aa": hydrophobic_count,
            "full_aa": mature_sequence if not reasons else "",
            "model_ready": not reasons,
            "qc_reason": ";".join(reasons) if reasons else "accepted",
        }
    )
    return result


def add_lookup_safety_columns(curated: pd.DataFrame) -> pd.DataFrame:
    result = curated.copy()
    result["two_field_name"] = result["name"].map(two_field_name)
    result["two_field_key"] = result["two_field_name"].str.upper()
    sequence_counts = result.groupby("two_field_key")["full_aa"].nunique()
    result["two_field_sequence_count"] = result["two_field_key"].map(sequence_counts)
    result["two_field_fallback_safe"] = result["two_field_sequence_count"].eq(1)
    return result


def select_model_input(curated: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    safe = curated[curated["two_field_fallback_safe"]].copy()
    conflicts = curated[~curated["two_field_fallback_safe"]].copy()

    selected_rows: list[pd.Series] = []
    for _, group in safe.groupby("two_field_key", sort=True):
        ordered = group.sort_values(
            ["name", "database", "source_accession", "source_index"],
            kind="stable",
        )
        selected = ordered.iloc[0].copy()
        selected["source_name"] = selected["name"]
        selected["source_allele_count"] = group["name"].nunique()
        selected["source_record_count"] = len(group)
        selected["name"] = selected["two_field_name"]
        selected["lookup_ready"] = True
        selected_rows.append(selected)

    model_input = pd.DataFrame(selected_rows)
    if not model_input.empty:
        model_input = model_input.sort_values(
            ["two_field_key", "source_name"], kind="stable"
        ).reset_index(drop=True)

    if not conflicts.empty:
        conflicts = conflicts.copy()
        conflicts["conflict_reason"] = "multiple_mature_protein_sequences_for_lookup_key"
        conflicts = conflicts.sort_values(
            ["two_field_key", "name", "source_accession"], kind="stable"
        ).reset_index(drop=True)
    return model_input, conflicts


def normalize_integer_columns(*tables: pd.DataFrame) -> None:
    for table in tables:
        for column in INTEGER_COLUMNS:
            if column in table.columns:
                table[column] = pd.to_numeric(table[column], errors="raise").astype("Int64")


def build_tables(
    data_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    source = read_source_records(data_dir)
    references = build_reference_sequences(source)
    aligner = new_aligner()
    assessments = pd.DataFrame(
        [assess_record(row, references, aligner) for row in source.itertuples(index=False)]
    )
    qc = pd.concat([source.reset_index(drop=True), assessments], axis=1)
    qc["two_field_name"] = qc["name"].map(two_field_name)

    curated = qc[qc["model_ready"]].copy()
    curated["full_aa_original"] = curated["sequence"]
    curated["af_ready_note"] = (
        "Alignment-derived mature MHC-II core; signal peptide, connecting stalk, "
        "transmembrane helix, and cytoplasmic tail omitted"
    )
    curated = add_lookup_safety_columns(curated)
    curated = curated.sort_values(
        ["two_field_key", "name", "source_accession"], kind="stable"
    ).reset_index(drop=True)

    model_input, conflicts = select_model_input(curated)
    normalize_integer_columns(qc, curated, model_input, conflicts)
    validate_tables(qc, curated, model_input, conflicts, references)

    reason_counts: Counter[str] = Counter()
    for reason_string in qc.loc[~qc["model_ready"], "qc_reason"]:
        reason_counts.update(str(reason_string).split(";"))

    raw_hashes = {
        filename: sha256(data_dir / "raw" / filename)
        for filename, _ in SOURCE_FILES
    }
    summary: dict[str, object] = {
        "script_version": SCRIPT_VERSION,
        "raw_sha256": raw_hashes,
        "thresholds": {
            "minimum_reference_coverage": MIN_REFERENCE_COVERAGE,
            "minimum_alignment_identity": MIN_ALIGNMENT_IDENTITY,
            "maximum_unaligned_reference_n": MAX_UNALIGNED_REFERENCE_N,
            "maximum_unaligned_reference_c": MAX_UNALIGNED_REFERENCE_C,
            "maximum_mature_length_deviation": MAX_MATURE_LENGTH_DEVIATION,
            "hydrophobic_window": HYDROPHOBIC_WINDOW,
            "maximum_hydrophobic_residues_per_window": MAX_HYDROPHOBIC_IN_WINDOW,
            "excluded_expression_suffixes": sorted(EXPRESSION_SUFFIXES),
        },
        "counts": {
            "raw_supported_records": len(qc),
            "model_ready_source_records": len(curated),
            "safe_lookup_records": len(model_input),
            "ambiguous_lookup_source_records": len(conflicts),
            "ambiguous_lookup_keys": int(conflicts["two_field_key"].nunique())
            if not conflicts.empty
            else 0,
        },
        "model_ready_by_species": {
            str(key): int(value)
            for key, value in curated.groupby("species").size().sort_index().items()
        },
        "model_ready_by_family_chain": {
            f"{family}_{chain}": int(value)
            for (family, chain), value in curated.groupby(["family", "chain"])
            .size()
            .sort_index()
            .items()
        },
        "rejection_reason_counts": dict(sorted(reason_counts.items())),
        "ambiguous_lookup_keys": sorted(conflicts["two_field_name"].unique().tolist())
        if not conflicts.empty
        else [],
        "references": {
            f"{family}_{chain}": {
                "allele": reference_slice.allele,
                "start_0based": reference_slice.start_0based,
                "end_exclusive": reference_slice.end_exclusive,
                "mature_length": len(reference_sequence),
                "mature_start": reference_sequence[:12],
                "mature_end": reference_sequence[-12:],
            }
            for (family, chain), (reference_slice, reference_sequence) in references.items()
        },
    }
    return qc, curated, model_input, conflicts, summary


def validate_tables(
    qc: pd.DataFrame,
    curated: pd.DataFrame,
    model_input: pd.DataFrame,
    conflicts: pd.DataFrame,
    references: dict[tuple[str, str], tuple[ReferenceSlice, str]],
) -> None:
    if qc.empty or curated.empty or model_input.empty:
        raise AssertionError("Rebuilt MHC-II tables must not be empty")
    if qc["header"].duplicated().any():
        raise AssertionError("Raw FASTA headers are not unique")
    if int(qc["model_ready"].sum()) != len(curated):
        raise AssertionError("Curated row count does not match QC model_ready count")
    if not curated["model_ready"].all():
        raise AssertionError("Curated table contains a rejected row")
    if curated["full_aa"].eq("").any():
        raise AssertionError("Curated table contains an empty mature sequence")
    if curated["expression_suffix"].ne("").any():
        raise AssertionError("Expression-suffix alleles entered the curated table")
    if model_input["name"].str.upper().duplicated().any():
        raise AssertionError("Model-input lookup names are not unique")
    if not model_input["two_field_sequence_count"].eq(1).all():
        raise AssertionError("Ambiguous protein sequence entered model-input lookup")
    if set(model_input["two_field_key"]) & set(conflicts.get("two_field_key", [])):
        raise AssertionError("Ambiguous lookup key entered model-input table")

    for row in curated.itertuples(index=False):
        mature = str(row.full_aa)
        original = str(row.sequence)
        start = int(row.mature_start_0based)
        end = int(row.mature_end_exclusive)
        if mature != original[start:end]:
            raise AssertionError(f"Boundary/sequence mismatch for {row.name}")
        if set(mature) - AA20:
            raise AssertionError(f"Invalid mature amino acid for {row.name}")
        if max_hydrophobic_count(mature) > MAX_HYDROPHOBIC_IN_WINDOW:
            raise AssertionError(f"Hydrophobic segment remains in {row.name}")

    model_by_name = model_input.set_index(model_input["name"].str.upper())
    for _, (reference_slice, reference_sequence) in references.items():
        key = two_field_name(reference_slice.allele).upper()
        if key not in model_by_name.index:
            raise AssertionError(f"Reference lookup key missing: {key}")
        if str(model_by_name.loc[key, "full_aa"]) != reference_sequence:
            raise AssertionError(f"Reference mature sequence changed: {key}")


def ordered_columns(table: pd.DataFrame, kind: str) -> list[str]:
    common = [
        "name",
        "gene",
        "family",
        "chain",
        "species",
        "species_prefix",
        "database",
        "source_accession",
        "length",
        "declared_length",
        "length_matches_header",
        "invalid_chars",
        "expression_suffix",
        "model_ready",
        "qc_reason",
        "reference_allele",
        "reference_length",
        "alignment_identity",
        "reference_coverage",
        "aligned_residues",
        "unaligned_reference_n",
        "unaligned_reference_c",
        "alignment_target_start_0based",
        "alignment_target_end_exclusive",
        "mature_start_0based",
        "mature_end_exclusive",
        "mature_length",
        "n_terminal_trim",
        "c_terminal_trim",
        "max_hydrophobic_residues_19aa",
        "two_field_name",
    ]
    if kind == "qc":
        desired = common + ["full_aa", "sequence", "header", "source_index"]
    else:
        desired = common + [
            "two_field_sequence_count",
            "two_field_fallback_safe",
            "full_aa",
            "full_aa_original",
            "af_ready_note",
            "sequence",
            "header",
            "source_index",
        ]
        if kind == "model":
            desired[1:1] = [
                "source_name",
                "source_allele_count",
                "source_record_count",
                "lookup_ready",
            ]
        elif kind == "conflict":
            desired.insert(1, "conflict_reason")
    return [column for column in desired if column in table.columns]


def atomic_write_csv(table: pd.DataFrame, path: Path, columns: Iterable[str]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    table.to_csv(temporary, index=False, columns=list(columns))
    temporary.replace(path)


def atomic_write_json(data: dict[str, object], path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def write_outputs(
    output_dir: Path,
    qc: pd.DataFrame,
    curated: pd.DataFrame,
    model_input: pd.DataFrame,
    conflicts: pd.DataFrame,
    summary: dict[str, object],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    atomic_write_csv(
        qc,
        output_dir / "mhcii_chain_qc_all.csv",
        ordered_columns(qc, "qc"),
    )
    atomic_write_csv(
        curated,
        output_dir / "mhcii_chain_curated.csv",
        ordered_columns(curated, "curated"),
    )
    atomic_write_csv(
        model_input,
        output_dir / "mhcii_model_input.csv",
        ordered_columns(model_input, "model"),
    )
    atomic_write_csv(
        conflicts,
        output_dir / "mhcii_two_field_conflicts.csv",
        ordered_columns(conflicts, "conflict"),
    )
    atomic_write_json(summary, output_dir / "mhcii_cleaning_summary.json")


def main() -> None:
    args = parse_args()
    data_dir = args.data_dir.expanduser().resolve()
    output_dir = (
        args.output_dir.expanduser().resolve() if args.output_dir else data_dir
    )
    qc, curated, model_input, conflicts, summary = build_tables(data_dir)
    write_outputs(output_dir, qc, curated, model_input, conflicts, summary)

    counts = summary["counts"]
    print(f"Raw supported records: {counts['raw_supported_records']}")
    print(f"Model-ready source records: {counts['model_ready_source_records']}")
    print(f"Safe lookup records: {counts['safe_lookup_records']}")
    print(f"Ambiguous lookup keys excluded: {counts['ambiguous_lookup_keys']}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
