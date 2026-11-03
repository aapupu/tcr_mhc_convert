#!/usr/bin/env python3
"""Build gene-level human/mouse CDR3-anchor and TRAJ-TRAC bridge metadata."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import tempfile
from pathlib import Path

from tcrtransform.references import (
    AA20,
    DATA_DIR,
    DEFAULT_TCRDIST_REFERENCE,
    DEFAULT_JUNCTION_REFERENCE,
    resolve_gene_reference_path,
)


DEFAULT_METADATA = DATA_DIR / "tcr_gene_junction_reference.metadata.json"
TRAC_SPLICE_ACCEPTOR_PREFIX = {"human": "AT", "mouse": "AC"}
CODON_TABLE = {
    "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
    "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
    "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
    "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tcrdist-reference", default=str(DEFAULT_TCRDIST_REFERENCE))
    parser.add_argument("--human-gene-reference", default=None)
    parser.add_argument("--mouse-gene-reference", default=None)
    parser.add_argument("--output", default=str(DEFAULT_JUNCTION_REFERENCE))
    parser.add_argument("--metadata-output", default=str(DEFAULT_METADATA))
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def natural_key(value: str) -> tuple[object, ...]:
    return tuple(int(part) if part.isdigit() else part for part in re.split(r"(\d+)", value))


def clean_aligned_sequence(value: object) -> str:
    return re.sub(r"[.\-*\s]", "", str(value).upper())


def load_gene_table(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"gene_name", "amino_acid_original"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        return {
            str(row["gene_name"]).strip().upper(): str(row["amino_acid_original"]).strip().upper()
            for row in reader
            if str(row["gene_name"]).strip().upper().startswith(("TRAV", "TRBV", "TRAJ", "TRBJ"))
        }


def load_tcrdist_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "id", "organism", "chain", "region", "nucseq", "frame",
            "aligned_protseq", "cdrs",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        return [dict(row) for row in reader]


def sequence_identity(left: str, right: str) -> tuple[int, int]:
    matches = sum(a == b for a, b in zip(left, right))
    return matches, -abs(len(left) - len(right))


def choose_reference_row(
    candidates: list[dict[str, str]],
    gene_sequence: str,
) -> dict[str, str]:
    usable = []
    for row in candidates:
        protein = clean_aligned_sequence(row["aligned_protseq"])
        if not protein or set(protein) - AA20:
            continue
        usable.append((sequence_identity(protein, gene_sequence), row))
    if not usable:
        raise ValueError(f"no standard tcrdist protein candidates for {candidates[0]['id']}")
    best_score = max(score for score, _ in usable)
    best = [row for score, row in usable if score == best_score]
    best.sort(
        key=lambda row: (
            0 if str(row["id"]).upper().endswith("*01") else 1,
            natural_key(str(row["id"]).upper()),
        )
    )
    return best[0]


def v_anchor(sequence: str, row: dict[str, str]) -> tuple[str, int]:
    parts = str(row["cdrs"]).split(";")
    fragment = clean_aligned_sequence(parts[3]) if len(parts) >= 4 else ""
    if fragment.startswith("C"):
        position = sequence.rfind(fragment)
        if position >= 0:
            return "C", position
    position = sequence.rfind("C", max(0, len(sequence) - 20))
    if position < 0:
        raise ValueError(f"cannot locate conserved V CDR3 C in {row['id']}")
    return "C", position


def j_anchor(sequence: str, row: dict[str, str]) -> tuple[str, int]:
    fragment = clean_aligned_sequence(str(row["cdrs"]).split(";", 1)[0])
    if not fragment:
        raise ValueError(f"missing J CDR3 annotation for {row['id']}")
    anchor = fragment[-1]
    position = sequence.find(fragment)
    if position >= 0:
        return anchor, position + len(fragment) - 1

    reference = clean_aligned_sequence(row["aligned_protseq"])
    reference_position = reference.find(fragment)
    expected = reference_position + len(fragment) - 1
    if 0 <= expected < len(sequence) and sequence[expected] == anchor:
        return anchor, expected

    candidates = [index for index, residue in enumerate(sequence[:15]) if residue == anchor]
    if len(candidates) == 1:
        return anchor, candidates[0]
    raise ValueError(f"cannot locate annotated J CDR3 anchor in {row['id']}")


def traj_bridge(species: str, row: dict[str, str]) -> tuple[str, str]:
    try:
        frame_offset = int(row["frame"]) - 1
    except ValueError:
        return "", ""
    nucleotide = re.sub(r"\s+", "", str(row["nucseq"]).upper())
    if frame_offset not in {0, 1, 2} or set(nucleotide) - set("ACGT"):
        return "", ""
    coding = nucleotide[frame_offset:]
    remainder = coding[(len(coding) // 3) * 3 :]
    # Productive TRAJ genes represented by the bundled reference retain one J
    # nucleotide; the first two TRAC exon nucleotides complete this codon.
    if len(remainder) != 1:
        return "", ""
    codon = remainder + TRAC_SPLICE_ACCEPTOR_PREFIX[species]
    residue = CODON_TABLE.get(codon, "")
    if not residue or residue == "*":
        return "", ""
    return residue, codon


def build_rows(reference_path: Path, gene_paths: dict[str, Path]) -> tuple[list[dict[str, object]], dict[str, object]]:
    tcrdist_rows = load_tcrdist_rows(reference_path)
    output: list[dict[str, object]] = []
    unsupported_traj_bridges: dict[str, list[str]] = {"human": [], "mouse": []}
    unsupported_cdr3_anchors: dict[str, list[str]] = {"human": [], "mouse": []}

    for species in ("human", "mouse"):
        genes = load_gene_table(gene_paths[species])
        by_gene: dict[str, list[dict[str, str]]] = {}
        for row in tcrdist_rows:
            gene = str(row["id"]).strip().upper().split("*", 1)[0]
            if (
                str(row["organism"]).strip().lower() == species
                and str(row["chain"]).strip().upper() in {"A", "B"}
                and str(row["region"]).strip().upper() in {"V", "J"}
                and gene in genes
            ):
                by_gene.setdefault(gene, []).append(row)

        for gene_name, gene_sequence in genes.items():
            candidates = by_gene.get(gene_name, [])
            if not candidates:
                raise ValueError(f"no tcrdist row for {species} {gene_name}")
            row = choose_reference_row(candidates, gene_sequence)
            region = str(row["region"]).strip().upper()
            try:
                if region == "V":
                    anchor, anchor_index = v_anchor(gene_sequence, row)
                else:
                    anchor, anchor_index = j_anchor(gene_sequence, row)
            except ValueError:
                unsupported_cdr3_anchors[species].append(gene_name)
                continue
            bridge = ""
            bridge_codon = ""
            if gene_name.startswith("TRAJ"):
                bridge, bridge_codon = traj_bridge(species, row)
                if not bridge:
                    unsupported_traj_bridges[species].append(gene_name)
            output.append(
                {
                    "species": species,
                    "gene_name": gene_name,
                    "reference_allele": str(row["id"]).strip().upper(),
                    "region": region,
                    "cdr3_anchor": anchor,
                    "cdr3_anchor_index": anchor_index,
                    "jc_bridge": bridge,
                    "jc_bridge_codon": bridge_codon,
                }
            )

    output.sort(key=lambda row: (str(row["species"]), natural_key(str(row["gene_name"]))))
    metadata = {
        "tcrdist_reference": reference_path.name,
        "tcrdist_reference_sha256": hashlib.sha256(reference_path.read_bytes()).hexdigest(),
        "gene_references": {
            species: {
                "file": path.name,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for species, path in gene_paths.items()
        },
        "coordinate_system": "cdr3_anchor_index is zero-based in amino_acid_original",
        "traj_trac_split_codon_rule": (
            "Translate the final in-frame TRAJ nucleotide together with the first two "
            "TRAC exon nucleotides (human AT; mouse AC)."
        ),
        "unsupported_traj_bridges": unsupported_traj_bridges,
        "unsupported_cdr3_anchors": unsupported_cdr3_anchors,
        "rows": len(output),
    }
    return output, metadata


def atomic_write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=path.parent,
        prefix=f".{path.name}.", suffix=".tmp", delete=False,
    ) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def atomic_write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent,
        prefix=f".{path.name}.", suffix=".tmp", delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def main() -> None:
    args = parse_args()
    output_path = Path(args.output).expanduser().resolve()
    metadata_path = Path(args.metadata_output).expanduser().resolve()
    for path in (output_path, metadata_path):
        if path.exists() and not args.overwrite:
            raise SystemExit(f"{path} exists; use --overwrite to replace it")
    gene_paths = {
        "human": Path(args.human_gene_reference).expanduser().resolve()
        if args.human_gene_reference else resolve_gene_reference_path("human"),
        "mouse": Path(args.mouse_gene_reference).expanduser().resolve()
        if args.mouse_gene_reference else resolve_gene_reference_path("mouse"),
    }
    rows, metadata = build_rows(
        Path(args.tcrdist_reference).expanduser().resolve(), gene_paths
    )
    atomic_write_csv(output_path, rows)
    atomic_write_json(metadata_path, metadata)
    print(f"Rows: {len(rows)}")
    print(f"Reference: {output_path}")
    print(f"Metadata: {metadata_path}")


if __name__ == "__main__":
    main()
