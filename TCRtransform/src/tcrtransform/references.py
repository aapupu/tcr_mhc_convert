"""Bundled reference paths and human/mouse TCR lookup helpers."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

from .errors import TCRReferenceError, TCRTransformError


DATA_DIR = Path(__file__).resolve().parent / "data"
DEFAULT_TCRDIST_REFERENCE = DATA_DIR / "combo_xcr_2024-03-05.tsv"
DEFAULT_JUNCTION_REFERENCE = DATA_DIR / "tcr_gene_junction_reference.csv"
TCR_GENE_AA_PATHS = {
    "human": DATA_DIR / "tcr_gene_aa_human_TRA_TRB_gene.csv",
    "mouse": DATA_DIR / "tcr_gene_aa_mouse_TRA_TRB_gene.csv",
}

AA20 = frozenset("ACDEFGHIKLMNPQRSTVWY")
AA_INPUT = AA20 | {"X"}
SPECIES_ALIASES = {
    "human": "human",
    "homo sapiens": "human",
    "homo_sapiens": "human",
    "hs": "human",
    "mouse": "mouse",
    "mus musculus": "mouse",
    "mus_musculus": "mouse",
    "mm": "mouse",
}


def normalize_species(value: object) -> str:
    key = re.sub(r"\s+", " ", str(value).strip().lower().replace("+", " "))
    species = SPECIES_ALIASES.get(key)
    if species is None:
        raise TCRTransformError(
            f"unsupported species {value!r}; expected human or mouse"
        )
    return species


def normalize_gene(value: object, label: str = "TCR gene") -> str:
    gene = re.sub(r"\s+", "", str(value).strip().upper())
    if not gene:
        raise TCRTransformError(f"{label} is empty")
    return gene


def clean_sequence(
    value: object,
    label: str = "amino-acid sequence",
    *,
    allow_empty: bool = False,
    allow_x: bool = True,
) -> str:
    raw = "" if value is None else str(value).strip().upper()
    if not raw:
        if allow_empty:
            return ""
        raise TCRTransformError(f"{label} is empty")
    sequence = re.sub(r"[\s.\-]", "", raw)
    alphabet = AA_INPUT if allow_x else AA20
    invalid = sorted(set(sequence) - alphabet)
    if invalid:
        raise TCRTransformError(
            f"{label} contains invalid amino-acid characters: {''.join(invalid)}"
        )
    if not sequence and not allow_empty:
        raise TCRTransformError(f"{label} is empty after cleaning")
    return sequence


def resolve_gene_reference_path(
    species: object,
    override: str | Path | None = None,
) -> Path:
    if override is not None and str(override).strip():
        return Path(override).expanduser().resolve()
    return TCR_GENE_AA_PATHS[normalize_species(species)]


def load_gene_sequences(
    species: object,
    *,
    path: str | Path | None = None,
    sequence_column: str = "amino_acid",
) -> dict[str, str]:
    reference_path = resolve_gene_reference_path(species, path)
    sequences: dict[str, str] = {}
    with reference_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"gene_name", sequence_column}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise TCRReferenceError(
                f"{reference_path} is missing columns: {sorted(missing)}"
            )
        for row_number, row in enumerate(reader, start=2):
            gene = normalize_gene(row["gene_name"], f"reference row {row_number} gene")
            sequence = clean_sequence(
                row[sequence_column],
                f"reference row {row_number} {sequence_column}",
                allow_x=False,
            )
            previous = sequences.setdefault(gene, sequence)
            if previous != sequence:
                raise TCRReferenceError(
                    f"conflicting sequences for {gene} in {reference_path}"
                )
    if not sequences:
        raise TCRReferenceError(f"{reference_path} contains no gene sequences")
    return sequences


def get_gene_sequence(gene_sequences: dict[str, str], gene_name: object) -> str:
    name = normalize_gene(gene_name)
    return gene_sequences.get(name) or gene_sequences.get(name.split("*", 1)[0], "")


def choose_trbc_gene(trbj: object, default: str = "TRBC1") -> str:
    """Map TRBJ1-* to TRBC1 and TRBJ2-* to TRBC2."""

    name = normalize_gene(trbj, "TRBJ")
    base = name.split("*", 1)[0]
    if base.startswith("TRBJ2"):
        return "TRBC2"
    if base.startswith("TRBJ1"):
        return "TRBC1"
    return default


def _ungap_loop(value: str) -> str:
    sequence = re.sub(r"[.\-\s]", "", str(value).upper())
    if not sequence or set(sequence) - AA20:
        raise TCRReferenceError(f"invalid CDR loop {value!r}")
    return sequence


def _has_standard_sequence(value: str) -> bool:
    sequence = re.sub(r"[.\-\s]", "", str(value).upper())
    return bool(sequence) and not (set(sequence) - AA20)


def _natural_allele_key(allele: str) -> tuple[object, ...]:
    suffix = allele.split("*", 1)[1] if "*" in allele else ""
    return tuple(
        int(part) if part.isdigit() else part
        for part in re.split(r"(\d+)", suffix)
    )


@dataclass(frozen=True)
class VGeneLoops:
    organism: str
    chain: str
    allele: str
    gene: str
    cdr1: str
    cdr2: str
    cdr2_5: str
    full_sequence_usable: bool


class VGeneLoopReference:
    """Allele-aware CDR1/CDR2/CDR2.5 lookup from the bundled tcrdist table."""

    def __init__(self, records: list[VGeneLoops], skipped_records: int = 0) -> None:
        self.skipped_records = skipped_records
        self.exact: dict[tuple[str, str, str], VGeneLoops] = {}
        grouped: dict[tuple[str, str, str], list[VGeneLoops]] = {}
        for record in records:
            key = (record.organism, record.chain, record.allele)
            previous = self.exact.setdefault(key, record)
            if previous != record:
                raise TCRReferenceError(
                    f"conflicting tcrdist rows for {record.organism} {record.allele}"
                )
            grouped.setdefault(
                (record.organism, record.chain, record.gene), []
            ).append(record)

        self.canonical: dict[tuple[str, str, str], VGeneLoops] = {}
        for key, candidates in grouped.items():
            usable = [record for record in candidates if record.full_sequence_usable]
            pool = usable or candidates
            pool.sort(
                key=lambda record: (
                    0 if record.allele.endswith("*01") else 1,
                    _natural_allele_key(record.allele),
                )
            )
            self.canonical[key] = pool[0]

    @classmethod
    def from_tsv(
        cls,
        path: str | Path = DEFAULT_TCRDIST_REFERENCE,
    ) -> "VGeneLoopReference":
        path = Path(path)
        records: list[VGeneLoops] = []
        skipped = 0
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "id", "organism", "chain", "region", "aligned_protseq", "cdrs"
            }
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise TCRReferenceError(f"{path} is missing columns: {sorted(missing)}")
            for row in reader:
                organism = str(row["organism"]).strip().lower()
                chain = str(row["chain"]).strip().upper()
                if organism not in {"human", "mouse"}:
                    continue
                if chain not in {"A", "B"} or str(row["region"]).strip().upper() != "V":
                    continue
                allele = normalize_gene(row["id"])
                expected = "TRAV" if chain == "A" else "TRBV"
                if not allele.startswith(expected):
                    continue
                loops = str(row["cdrs"]).split(";")
                if len(loops) < 3:
                    skipped += 1
                    continue
                try:
                    cdr1, cdr2, cdr2_5 = (_ungap_loop(loop) for loop in loops[:3])
                except TCRReferenceError:
                    skipped += 1
                    continue
                records.append(
                    VGeneLoops(
                        organism=organism,
                        chain=chain,
                        allele=allele,
                        gene=allele.split("*", 1)[0],
                        cdr1=cdr1,
                        cdr2=cdr2,
                        cdr2_5=cdr2_5,
                        full_sequence_usable=_has_standard_sequence(
                            row["aligned_protseq"]
                        ),
                    )
                )
        if not records:
            raise TCRReferenceError(f"{path} contains no usable TRAV/TRBV rows")
        return cls(records, skipped_records=skipped)

    def resolve(self, species: object, chain: str, gene_name: object) -> VGeneLoops:
        organism = normalize_species(species)
        chain = str(chain).strip().upper()
        if chain not in {"A", "B"}:
            raise TCRTransformError(f"unsupported chain {chain!r}; expected A or B")
        gene_name = normalize_gene(gene_name, "V gene")
        expected = "TRAV" if chain == "A" else "TRBV"
        if not gene_name.startswith(expected):
            raise TCRTransformError(
                f"{gene_name!r} is not a {expected} gene for chain {chain}"
            )
        if "*" in gene_name:
            record = self.exact.get((organism, chain, gene_name))
            if record is None:
                raise TCRReferenceError(
                    f"unrecognized {organism} allele-specific V gene {gene_name!r}"
                )
            return record
        record = self.canonical.get((organism, chain, gene_name))
        if record is None:
            raise TCRReferenceError(f"unrecognized {organism} V gene {gene_name!r}")
        return record

