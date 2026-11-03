#!/usr/bin/env python3
"""Convert TCR gene/CDR3 and MHC allele calls to five protein chains.

The output columns are fixed as:
  A = TCR alpha
  B = TCR beta
  C = MHC alpha (MHC-I heavy chain or MHC-II alpha)
  D = beta-2-microglobulin for MHC-I, or MHC-II beta
  E = peptide

TCR chains are reconstructed with trimming-aware V-CDR3 and CDR3-J
matches. They are not built by directly concatenating V + CDR3 + J + C.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

from dataprocess import choose_trbc_gene, find_overlap


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_GENE_AA = PROJECT_ROOT / "doc/tcr_gene_aa_human_TRA_TRB_gene.csv"
DEFAULT_MHCI = (
    PROJECT_ROOT
    / "doc/mhci_full_ipd/mhci_full_ipd_curated_model_input_species_matched_b2m.csv"
)
DEFAULT_MHCII = PROJECT_ROOT / "doc/mhcii_full_ipd/mhcii_model_input.csv"

OUTPUT_COLUMNS = ["A", "B", "C", "D", "E"]
AA_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYX")

COLUMN_ALIASES = {
    "mhc": ("MHC", "mhc", "HLA", "hla", "mhc_name"),
    "peptide": ("peptide", "Peptide", "PEPTIDE", "epitope", "Epitope"),
    "cdr3a": ("CDR3A", "cdr3a", "CDR3a", "cdr3a_aa"),
    "cdr3b": ("CDR3B", "cdr3b", "CDR3b", "cdr3b_aa"),
    "trav": ("TRAV", "trav", "av_gene"),
    "traj": ("TRAJ", "traj", "aj_gene"),
    "trbv": ("TRBV", "trbv", "bv_gene", "v_gene"),
    "trbj": ("TRBJ", "trbj", "bj_gene", "j_gene"),
}

HUMAN_MHC_LOCI = (
    "DPA1",
    "DPB1",
    "DQA1",
    "DQB1",
    "DRA",
    "DRB1",
    "DRB2",
    "DRB3",
    "DRB4",
    "DRB5",
    "DRB6",
    "DRB7",
    "DRB8",
    "DRB9",
    "A",
    "B",
    "C",
    "E",
    "F",
    "G",
)


class ConversionError(ValueError):
    """A row cannot be converted without making an unsafe assumption."""


@dataclass(frozen=True)
class TerminalMatch:
    length: int
    sequence: str
    gene_start: int
    gene_end: int
    offset: int


@dataclass(frozen=True)
class MHCIRecord:
    name: str
    alpha_sequence: str
    b2m_sequence: str


@dataclass(frozen=True)
class MHCIIRecord:
    name: str
    sequence: str
    chain: str
    family: str
    species: str


class AlleleLookup:
    """Prefer an exact allele and fall back deterministically to two fields."""

    def __init__(self) -> None:
        self.exact: dict[str, object] = {}
        self.two_field: dict[str, object] = {}

    def add(self, allele_name: str, record: object) -> None:
        exact_key = normalize_allele_name(allele_name)
        if not exact_key:
            return
        self.exact.setdefault(exact_key, record)
        self.two_field.setdefault(two_field_allele(exact_key), record)

    def get(self, allele_name: str) -> object | None:
        key = normalize_allele_name(allele_name)
        if not key:
            return None
        return self.exact.get(key) or self.two_field.get(two_field_allele(key))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert MHC/peptide/CDR3A/CDR3B/TRAV/TRAJ/TRBV/TRBJ CSV rows "
            "to A=TCRa, B=TCRb, C=MHCa, D=B2M-or-MHCb, E=peptide."
        )
    )
    parser.add_argument("input_csv", help="Input CSV path.")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output CSV path. Default: <input_stem>_full_sequences.csv",
    )
    parser.add_argument("--gene-aa-path", default=str(DEFAULT_GENE_AA))
    parser.add_argument(
        "--gene-seq-column",
        default="amino_acid",
        help=(
            "TCR reference sequence column. The default is the current "
            "structure-model-ready sequence; use amino_acid_original if needed."
        ),
    )
    parser.add_argument("--mhci-path", default=str(DEFAULT_MHCI))
    parser.add_argument(
        "--mhci-seq-column",
        default="full_aa",
        help="MHC-I alpha sequence column (default: model-ready full_aa).",
    )
    parser.add_argument("--mhcii-path", default=str(DEFAULT_MHCII))
    parser.add_argument(
        "--mhcii-seq-column",
        default="full_aa",
        help="MHC-II chain sequence column (default: model-ready full_aa).",
    )
    parser.add_argument(
        "--default-dr-alpha",
        default="HLA-DRA*01:01",
        help=(
            "Alpha allele used only when a single human HLA-DRB allele is given. "
            "Set to an empty string to require an explicit DRA/DRB pair."
        ),
    )
    parser.add_argument("--min-overlap", type=int, default=3)
    parser.add_argument("--max-overlap", type=int, default=30)
    parser.add_argument("--v-terminal-window", type=int, default=8)
    parser.add_argument("--j-start-window", type=int, default=10)
    parser.add_argument(
        "--on-error",
        choices=("fail", "skip"),
        default="fail",
        help="fail: write nothing if any row is invalid; skip: write valid rows and a rejection CSV.",
    )
    parser.add_argument(
        "--rejected-output",
        default=None,
        help="Rejection CSV used with --on-error skip. Default: <output_stem>_rejected.csv",
    )
    return parser.parse_args()


def clean_sequence(value: object, label: str, allow_empty: bool = False) -> str:
    if value is None or pd.isna(value):
        if allow_empty:
            return ""
        raise ConversionError(f"{label} is empty")

    raw = str(value).strip().upper()
    if not raw:
        if allow_empty:
            return ""
        raise ConversionError(f"{label} is empty")

    sequence = re.sub(r"[\s.\-]", "", raw)
    invalid = sorted(set(sequence) - AA_ALPHABET)
    if invalid:
        raise ConversionError(f"{label} contains invalid amino-acid characters: {''.join(invalid)}")
    if not sequence and not allow_empty:
        raise ConversionError(f"{label} is empty after cleaning")
    return sequence


def clean_gene_name(value: object, label: str) -> str:
    if value is None or pd.isna(value):
        raise ConversionError(f"{label} is empty")
    name = re.sub(r"\s+", "", str(value).upper())
    if not name:
        raise ConversionError(f"{label} is empty")
    return name


def require_gene_family(gene_name: str, expected_prefix: str, label: str) -> None:
    base_name = gene_name.split("*", 1)[0]
    if not base_name.startswith(expected_prefix):
        raise ConversionError(
            f"{label} value {gene_name!r} is not an {expected_prefix} gene"
        )


def resolve_input_columns(columns: Iterable[object]) -> dict[str, str]:
    stripped = [str(column).strip() for column in columns]
    lower_to_columns: dict[str, list[str]] = {}
    for column in stripped:
        lower_to_columns.setdefault(column.lower(), []).append(column)

    resolved: dict[str, str] = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        matches: list[str] = []
        for alias in aliases:
            matches.extend(lower_to_columns.get(alias.lower(), []))
        matches = list(dict.fromkeys(matches))
        if not matches:
            raise ConversionError(
                f"missing input column for {canonical}; accepted names: {', '.join(aliases)}"
            )
        if len(matches) > 1:
            raise ConversionError(
                f"ambiguous columns for {canonical}: {matches}; keep only one accepted alias"
            )
        resolved[canonical] = matches[0]
    return resolved


def load_gene_sequences(path: Path, sequence_column: str) -> dict[str, str]:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"gene_name", sequence_column}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ConversionError(f"{path} is missing columns: {missing}")

    sequences: dict[str, str] = {}
    for row_number, row in table.iterrows():
        name = clean_gene_name(row["gene_name"], f"gene table row {row_number + 2} gene_name")
        sequence = clean_sequence(
            row[sequence_column],
            f"gene table row {row_number + 2} {sequence_column}",
        )
        sequences.setdefault(name, sequence)
    return sequences


def get_gene_sequence(gene_sequences: dict[str, str], gene_name: str) -> str:
    name = clean_gene_name(gene_name, "TCR gene")
    sequence = gene_sequences.get(name)
    if sequence:
        return sequence
    return gene_sequences.get(name.split("*", 1)[0], "")


def match_v_terminal(
    v_sequence: str,
    cdr3: str,
    min_overlap: int,
    max_overlap: int,
    terminal_window: int,
) -> TerminalMatch | None:
    max_k = min(len(v_sequence), len(cdr3), max_overlap)
    for length in range(max_k, min_overlap - 1, -1):
        query = cdr3[:length]
        candidates: list[tuple[int, int, int]] = []
        start = 0
        while True:
            position = v_sequence.find(query, start)
            if position < 0:
                break
            end = position + length
            offset = len(v_sequence) - end
            if 0 <= offset <= terminal_window:
                candidates.append((offset, -position, position))
            start = position + 1
        if candidates:
            offset, _, position = min(candidates)
            return TerminalMatch(length, query, position, position + length, offset)
    return None


def match_j_terminal(
    cdr3: str,
    j_sequence: str,
    min_overlap: int,
    max_overlap: int,
    start_window: int,
) -> TerminalMatch | None:
    max_k = min(len(cdr3), len(j_sequence), max_overlap)
    for length in range(max_k, min_overlap - 1, -1):
        query = cdr3[-length:]
        candidates: list[tuple[int, int]] = []
        start = 0
        while True:
            position = j_sequence.find(query, start)
            if position < 0:
                break
            if 0 <= position <= start_window:
                candidates.append((position, position + length))
            start = position + 1
        if candidates:
            position, end = min(candidates)
            return TerminalMatch(length, query, position, end, position)
    return None


def append_constant_region(variable_chain: str, constant_sequence: str) -> str:
    overlap_length, _ = find_overlap(
        variable_chain,
        constant_sequence,
        min_overlap=1,
        max_overlap=30,
    )
    return variable_chain + constant_sequence[overlap_length:]


def build_tcr_chain(
    *,
    cdr3: str,
    v_gene: str,
    j_gene: str,
    constant_gene: str,
    gene_sequences: dict[str, str],
    min_overlap: int,
    max_overlap: int,
    v_terminal_window: int,
    j_start_window: int,
    chain_label: str,
) -> str:
    v_sequence = get_gene_sequence(gene_sequences, v_gene)
    j_sequence = get_gene_sequence(gene_sequences, j_gene)
    constant_sequence = get_gene_sequence(gene_sequences, constant_gene)

    missing = []
    if not v_sequence:
        missing.append(v_gene)
    if not j_sequence:
        missing.append(j_gene)
    if not constant_sequence:
        missing.append(constant_gene)
    if missing:
        raise ConversionError(f"{chain_label} reference gene sequence missing: {', '.join(missing)}")

    v_match = match_v_terminal(
        v_sequence,
        cdr3,
        min_overlap=min_overlap,
        max_overlap=max_overlap,
        terminal_window=v_terminal_window,
    )
    if v_match is None:
        raise ConversionError(
            f"{chain_label} CDR3 prefix has no trimming-aware overlap with {v_gene}"
        )

    j_match = match_j_terminal(
        cdr3,
        j_sequence,
        min_overlap=min_overlap,
        max_overlap=max_overlap,
        start_window=j_start_window,
    )
    if j_match is None:
        raise ConversionError(
            f"{chain_label} CDR3 suffix has no trimming-aware overlap with {j_gene}"
        )

    # Keep the matched V-side CDR3 prefix once, insert the remaining CDR3,
    # then continue J after its matched CDR3 suffix.
    variable_chain = (
        v_sequence[: v_match.gene_end]
        + cdr3[v_match.length :]
        + j_sequence[j_match.gene_end :]
    )
    full_chain = append_constant_region(variable_chain, constant_sequence)
    if cdr3 not in full_chain:
        raise ConversionError(f"internal error: {chain_label} output does not contain its CDR3")
    return full_chain


def normalize_allele_name(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    name = str(value).strip().upper()
    if not name or name == "NAN":
        return ""
    name = name.replace(" ", "").replace("_", "-").replace("：", ":")
    name = re.sub(r"^HLA(?=[A-Z])", "HLA-", name)

    if not name.startswith("HLA-"):
        for locus in HUMAN_MHC_LOCI:
            if re.match(rf"^{re.escape(locus)}(?:\*|\d)", name):
                name = "HLA-" + name
                break

    if "*" not in name:
        loci_pattern = "|".join(re.escape(locus) for locus in HUMAN_MHC_LOCI)
        match = re.match(rf"^(HLA-(?:{loci_pattern}))(\d+(?::\d+)+[A-Z]?)$", name)
        if match:
            name = f"{match.group(1)}*{match.group(2)}"
    return name


def two_field_allele(value: object) -> str:
    name = normalize_allele_name(value)
    if "*" not in name:
        return name
    prefix, allele = name.split("*", 1)
    fields = allele.split(":")
    if len(fields) < 2:
        return name
    return f"{prefix}*{fields[0]}:{fields[1]}"


def load_mhci_lookup(path: Path, sequence_column: str) -> AlleleLookup:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"MHC", sequence_column, "b2m_aa"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ConversionError(f"{path} is missing columns: {missing}")

    lookup = AlleleLookup()
    for row_number, row in table.iterrows():
        name = normalize_allele_name(row["MHC"])
        if not name:
            continue
        alpha = clean_sequence(
            row[sequence_column],
            f"MHC-I table row {row_number + 2} {sequence_column}",
            allow_empty=True,
        )
        b2m = clean_sequence(
            row["b2m_aa"],
            f"MHC-I table row {row_number + 2} b2m_aa",
            allow_empty=True,
        )
        if alpha and b2m:
            lookup.add(name, MHCIRecord(name, alpha, b2m))
    return lookup


def parse_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def load_mhcii_lookup(path: Path, sequence_column: str) -> AlleleLookup:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"name", "chain", "family", "species", sequence_column}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ConversionError(f"{path} is missing columns: {missing}")

    lookup = AlleleLookup()
    for row_number, row in table.iterrows():
        if "model_ready" in table.columns and not parse_bool(row["model_ready"]):
            continue
        name = normalize_allele_name(row["name"])
        if not name:
            continue
        sequence = clean_sequence(
            row[sequence_column],
            f"MHC-II table row {row_number + 2} {sequence_column}",
            allow_empty=True,
        )
        chain = str(row["chain"]).strip().lower()
        family = str(row["family"]).strip().upper()
        species = str(row["species"]).strip().lower()
        if sequence and chain in {"alpha", "beta"}:
            lookup.add(name, MHCIIRecord(name, sequence, chain, family, species))
    return lookup


def split_mhc_alleles(mhc_value: object) -> list[str]:
    if mhc_value is None or pd.isna(mhc_value):
        raise ConversionError("MHC is empty")
    raw = str(mhc_value).strip()
    if not raw:
        raise ConversionError("MHC is empty")
    parts = [part.strip() for part in re.split(r"[/+;|,&]", raw) if part.strip()]
    if not parts:
        raise ConversionError("MHC is empty")
    return parts


def resolve_mhc_chains(
    mhc_value: object,
    mhci_lookup: AlleleLookup,
    mhcii_lookup: AlleleLookup,
    default_dr_alpha: str,
) -> tuple[str, str]:
    alleles = split_mhc_alleles(mhc_value)

    if len(alleles) == 1:
        mhci = mhci_lookup.get(alleles[0])
        if isinstance(mhci, MHCIRecord):
            return mhci.alpha_sequence, mhci.b2m_sequence

        mhcii = mhcii_lookup.get(alleles[0])
        if not isinstance(mhcii, MHCIIRecord):
            raise ConversionError(f"MHC allele not found in MHC-I or MHC-II tables: {alleles[0]}")
        if mhcii.chain == "beta" and mhcii.family == "DR" and mhcii.species == "human":
            if not str(default_dr_alpha).strip():
                raise ConversionError(
                    "a single HLA-DRB allele requires an explicit DRA allele or --default-dr-alpha"
                )
            alpha = mhcii_lookup.get(default_dr_alpha)
            if not isinstance(alpha, MHCIIRecord) or alpha.chain != "alpha":
                raise ConversionError(f"default DR alpha allele not found: {default_dr_alpha}")
            if alpha.family != "DR" or alpha.species != mhcii.species:
                raise ConversionError(
                    f"default DR alpha {alpha.name} is incompatible with beta {mhcii.name}"
                )
            return alpha.sequence, mhcii.sequence
        if mhcii.chain == "alpha":
            raise ConversionError(
                f"single MHC-II alpha allele {mhcii.name} does not determine its beta partner"
            )
        raise ConversionError(
            f"single {mhcii.name} does not determine a unique MHC-II alpha partner; "
            "put both alleles in MHC separated by '/' or '+'"
        )

    if len(alleles) != 2:
        raise ConversionError(
            f"MHC-II requires exactly two allele names, but found {len(alleles)}: {alleles}"
        )

    records = [mhcii_lookup.get(allele) for allele in alleles]
    missing = [allele for allele, record in zip(alleles, records) if not isinstance(record, MHCIIRecord)]
    if missing:
        raise ConversionError(f"MHC-II allele not found: {', '.join(missing)}")

    typed_records = [record for record in records if isinstance(record, MHCIIRecord)]
    alpha_records = [record for record in typed_records if record.chain == "alpha"]
    beta_records = [record for record in typed_records if record.chain == "beta"]
    if len(alpha_records) != 1 or len(beta_records) != 1:
        raise ConversionError("MHC-II pair must contain one alpha allele and one beta allele")

    alpha = alpha_records[0]
    beta = beta_records[0]
    if alpha.family != beta.family:
        raise ConversionError(
            f"MHC-II family mismatch: {alpha.name} ({alpha.family}) vs {beta.name} ({beta.family})"
        )
    if alpha.species != beta.species:
        raise ConversionError(
            f"MHC-II species mismatch: {alpha.name} ({alpha.species}) vs "
            f"{beta.name} ({beta.species})"
        )
    return alpha.sequence, beta.sequence


def convert_row(
    row: pd.Series,
    columns: dict[str, str],
    gene_sequences: dict[str, str],
    mhci_lookup: AlleleLookup,
    mhcii_lookup: AlleleLookup,
    args: argparse.Namespace,
) -> dict[str, str]:
    cdr3a = clean_sequence(row[columns["cdr3a"]], "CDR3A")
    cdr3b = clean_sequence(row[columns["cdr3b"]], "CDR3B")
    if not cdr3a.startswith("C"):
        raise ConversionError("CDR3A must start with the conserved C residue")
    if not cdr3b.startswith("C"):
        raise ConversionError("CDR3B must start with the conserved C residue")
    peptide = clean_sequence(row[columns["peptide"]], "peptide")
    trav = clean_gene_name(row[columns["trav"]], "TRAV")
    traj = clean_gene_name(row[columns["traj"]], "TRAJ")
    trbv = clean_gene_name(row[columns["trbv"]], "TRBV")
    trbj = clean_gene_name(row[columns["trbj"]], "TRBJ")
    require_gene_family(trav, "TRAV", "TRAV")
    require_gene_family(traj, "TRAJ", "TRAJ")
    require_gene_family(trbv, "TRBV", "TRBV")
    require_gene_family(trbj, "TRBJ", "TRBJ")

    tcra = build_tcr_chain(
        cdr3=cdr3a,
        v_gene=trav,
        j_gene=traj,
        constant_gene="TRAC",
        gene_sequences=gene_sequences,
        min_overlap=args.min_overlap,
        max_overlap=args.max_overlap,
        v_terminal_window=args.v_terminal_window,
        j_start_window=args.j_start_window,
        chain_label="TCR alpha",
    )
    tcrb = build_tcr_chain(
        cdr3=cdr3b,
        v_gene=trbv,
        j_gene=trbj,
        constant_gene=choose_trbc_gene(trbj),
        gene_sequences=gene_sequences,
        min_overlap=args.min_overlap,
        max_overlap=args.max_overlap,
        v_terminal_window=args.v_terminal_window,
        j_start_window=args.j_start_window,
        chain_label="TCR beta",
    )
    mhca, mhc_partner = resolve_mhc_chains(
        row[columns["mhc"]],
        mhci_lookup=mhci_lookup,
        mhcii_lookup=mhcii_lookup,
        default_dr_alpha=args.default_dr_alpha,
    )
    return {"A": tcra, "B": tcrb, "C": mhca, "D": mhc_partner, "E": peptide}


def inferred_output_path(input_path: Path) -> Path:
    return input_path.with_name(f"{input_path.stem}_full_sequences.csv")


def inferred_rejected_path(output_path: Path) -> Path:
    return output_path.with_name(f"{output_path.stem}_rejected.csv")


def format_error_summary(errors: list[dict[str, object]], max_examples: int = 12) -> str:
    lines = [f"Conversion failed for {len(errors)} row(s):"]
    for error in errors[:max_examples]:
        lines.append(f"  CSV row {error['csv_row']}: {error['conversion_error']}")
    if len(errors) > max_examples:
        lines.append(f"  ... and {len(errors) - max_examples} more")
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    if args.min_overlap < 1:
        raise SystemExit("--min-overlap must be at least 1")
    if args.max_overlap < args.min_overlap:
        raise SystemExit("--max-overlap must be greater than or equal to --min-overlap")
    if args.v_terminal_window < 0 or args.j_start_window < 0:
        raise SystemExit("terminal windows must be non-negative")

    input_path = Path(args.input_csv).expanduser().resolve()
    output_path = (
        Path(args.output).expanduser().resolve()
        if args.output
        else inferred_output_path(input_path)
    )

    try:
        input_table = pd.read_csv(input_path, dtype=str, keep_default_na=False)
        if input_table.empty:
            raise ConversionError(f"{input_path} contains no data rows")
        input_table.columns = [str(column).strip() for column in input_table.columns]
        if input_table.columns.duplicated().any():
            duplicates = input_table.columns[input_table.columns.duplicated()].tolist()
            raise ConversionError(f"duplicate input column names after stripping spaces: {duplicates}")
        columns = resolve_input_columns(input_table.columns)
        gene_sequences = load_gene_sequences(Path(args.gene_aa_path), args.gene_seq_column)
        mhci_lookup = load_mhci_lookup(Path(args.mhci_path), args.mhci_seq_column)
        mhcii_lookup = load_mhcii_lookup(Path(args.mhcii_path), args.mhcii_seq_column)
    except (
        ConversionError,
        FileNotFoundError,
        pd.errors.EmptyDataError,
        pd.errors.ParserError,
    ) as exc:
        raise SystemExit(f"Input/reference error: {exc}") from exc

    converted_rows: list[dict[str, str]] = []
    errors: list[dict[str, object]] = []
    rejected_indices: list[int] = []
    for index, row in input_table.iterrows():
        try:
            converted_rows.append(
                convert_row(
                    row,
                    columns=columns,
                    gene_sequences=gene_sequences,
                    mhci_lookup=mhci_lookup,
                    mhcii_lookup=mhcii_lookup,
                    args=args,
                )
            )
        except ConversionError as exc:
            rejected_indices.append(index)
            errors.append({"csv_row": index + 2, "conversion_error": str(exc)})

    if errors and args.on_error == "fail":
        raise SystemExit(format_error_summary(errors))
    if not converted_rows:
        raise SystemExit(format_error_summary(errors) + "\nNo valid rows remain; no output was written.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_table = pd.DataFrame(converted_rows, columns=OUTPUT_COLUMNS)
    output_table.to_csv(output_path, index=False)

    rejected_path: Path | None = None
    if errors and args.on_error == "skip":
        rejected_path = (
            Path(args.rejected_output).expanduser().resolve()
            if args.rejected_output
            else inferred_rejected_path(output_path)
        )
        rejected_path.parent.mkdir(parents=True, exist_ok=True)
        rejected = input_table.loc[rejected_indices].copy()
        rejected.insert(0, "csv_row", [error["csv_row"] for error in errors])
        rejected["conversion_error"] = [error["conversion_error"] for error in errors]
        rejected.to_csv(rejected_path, index=False)

    print(f"Input rows: {len(input_table)}")
    print(f"Converted rows: {len(output_table)}")
    print(f"Rejected rows: {len(errors)}")
    print(f"Output: {output_path}")
    if rejected_path is not None:
        print(f"Rejected output: {rejected_path}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Interrupted", file=sys.stderr)
        raise SystemExit(130)
