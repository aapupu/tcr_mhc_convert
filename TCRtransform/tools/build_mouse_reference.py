#!/usr/bin/env python3
"""Build a gene-level mouse TRA/TRB amino-acid reference from tcrdist data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import tempfile
from collections import defaultdict
from pathlib import Path

from tcrtransform.references import AA20, DATA_DIR, DEFAULT_TCRDIST_REFERENCE


DEFAULT_OUTPUT = DATA_DIR / "tcr_gene_aa_mouse_TRA_TRB_gene.csv"
DEFAULT_METADATA = DATA_DIR / "tcr_gene_aa_mouse_TRA_TRB_gene.metadata.json"
TCRDIST_COMMIT = "0e6efa39092e2c64c9c0b52a72c2f5961c161ebb"

# IMGT mouse C-domain displays. The model-ready column intentionally excludes
# connecting, transmembrane, and cytoplasmic regions, matching the human table.
MOUSE_TRAC = (
    "IQNPEPAVYQLKDPRSQDSTLCLFTDFDSQINVPKTMESGTFITDKTVLDMKAMDSKSNGAIA"
    "WSNQTSFTCQDIFKETNATYPSS"
)
MOUSE_TRBC_DOMAIN = (
    "EDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARGFFPDHVELSWWVNGKEVHSGVSTDPQAY"
    "KESNYSYCLSSRLRVSATFWHNPRNHFRCQVQFHGLSEEDKWPEGSPKPVTQNISAEAWGRA"
)
MOUSE_TRBC1_FULL = (
    MOUSE_TRBC_DOMAIN
    + "DCGITSASYQQGVLSATILYEILLGKATLYAVLVSTLVVMAMVKRKNS"
)
MOUSE_TRBC2_FULL = (
    MOUSE_TRBC_DOMAIN
    + "DCGITSASYHQGVLSATILYEILLGKATLYAVLVSGLVLMAMVKKKNS"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", default=str(DEFAULT_TCRDIST_REFERENCE))
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--metadata-output", default=str(DEFAULT_METADATA))
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def natural_key(value: str) -> tuple[object, ...]:
    return tuple(
        int(part) if part.isdigit() else part
        for part in re.split(r"(\d+)", value)
    )


def clean_aligned_sequence(value: str) -> str:
    sequence = re.sub(r"[.\-\s]", "", str(value).upper())
    if not sequence or set(sequence) - AA20:
        return ""
    return sequence


def atomic_write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temp_path = Path(handle.name)
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "gene_name",
                "amino_acid",
                "amino_acid_original",
                "af_trimmed_len",
                "af_trim_method",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp_path, path)


def atomic_write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temp_path = Path(handle.name)
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temp_path, path)


def build_rows(reference_path: Path) -> tuple[list[dict[str, object]], dict[str, object]]:
    groups: dict[str, list[tuple[str, str]]] = defaultdict(list)
    seen: dict[str, str] = {}
    excluded_nonstandard: list[str] = []
    specifications = {
        ("A", "V"): "TRAV",
        ("A", "J"): "TRAJ",
        ("B", "V"): "TRBV",
        ("B", "J"): "TRBJ",
    }

    with reference_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"id", "organism", "chain", "region", "aligned_protseq"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{reference_path} is missing columns: {sorted(missing)}"
            )

        for row in reader:
            if row["organism"].strip().lower() != "mouse":
                continue
            chain_region = (
                row["chain"].strip().upper(),
                row["region"].strip().upper(),
            )
            prefix = specifications.get(chain_region)
            allele = row["id"].strip().upper()
            if prefix is None or not allele.startswith(prefix):
                continue
            sequence = clean_aligned_sequence(row["aligned_protseq"])
            previous_sequence = seen.get(allele)
            if previous_sequence is not None:
                if previous_sequence != sequence:
                    raise ValueError(
                        f"conflicting duplicate tcrdist sequences for mouse {allele}"
                    )
                continue
            seen[allele] = sequence
            if not sequence:
                excluded_nonstandard.append(allele)
                continue
            groups[allele.split("*", 1)[0]].append((allele, sequence))

    selected_alleles: dict[str, str] = {}
    rows: list[dict[str, object]] = []
    for gene, candidates in groups.items():
        candidates.sort(
            key=lambda item: (
                0 if item[0].endswith("*01") else 1,
                natural_key(item[0]),
            )
        )
        allele, sequence = candidates[0]
        selected_alleles[gene] = allele
        rows.append(
            {
                "gene_name": gene,
                "amino_acid": sequence,
                "amino_acid_original": sequence,
                "af_trimmed_len": len(sequence),
                "af_trim_method": "unchanged",
            }
        )

    rows.extend(
        [
            {
                "gene_name": "TRAC",
                "amino_acid": MOUSE_TRAC,
                "amino_acid_original": MOUSE_TRAC,
                "af_trimmed_len": len(MOUSE_TRAC),
                "af_trim_method": "unchanged",
            },
            {
                "gene_name": "TRBC1",
                "amino_acid": MOUSE_TRBC_DOMAIN,
                "amino_acid_original": MOUSE_TRBC1_FULL,
                "af_trimmed_len": len(MOUSE_TRBC_DOMAIN),
                "af_trim_method": "TRBC_structural_domain_end_TQNISAEAWGRA",
            },
            {
                "gene_name": "TRBC2",
                "amino_acid": MOUSE_TRBC_DOMAIN,
                "amino_acid_original": MOUSE_TRBC2_FULL,
                "af_trimmed_len": len(MOUSE_TRBC_DOMAIN),
                "af_trim_method": "TRBC_structural_domain_end_TQNISAEAWGRA",
            },
        ]
    )
    rows.sort(key=lambda row: natural_key(str(row["gene_name"])))

    sha256 = hashlib.sha256(reference_path.read_bytes()).hexdigest()
    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        match = re.match(r"^(TRAC|TRAJ|TRAV|TRBC|TRBJ|TRBV)", str(row["gene_name"]))
        counts[match.group(1) if match else "other"] += 1
    metadata = {
        "species": "Mus musculus",
        "tcrdist_reference": reference_path.name,
        "tcrdist_reference_sha256": sha256,
        "tcrdist_git_commit": TCRDIST_COMMIT,
        "canonical_allele_rule": (
            "Prefer a usable *01 allele; otherwise use the naturally lowest "
            "allele with a standard amino-acid sequence."
        ),
        "constant_region_sources": {
            "TRAC_TRBC": "IMGT mouse TRA/TRB C-REGION protein display",
            "TRBC2_full_length_check": "UniProtKB P01851",
        },
        "gene_counts": dict(sorted(counts.items())),
        "selected_alleles": dict(sorted(selected_alleles.items())),
        "excluded_nonstandard_alleles": sorted(set(excluded_nonstandard)),
    }
    return rows, metadata


def main() -> None:
    args = parse_args()
    reference_path = Path(args.reference).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()
    metadata_path = Path(args.metadata_output).expanduser().resolve()
    for path in (output_path, metadata_path):
        if path.exists() and not args.overwrite:
            raise SystemExit(f"{path} exists; use --overwrite to replace it")

    try:
        rows, metadata = build_rows(reference_path)
        atomic_write_csv(output_path, rows)
        atomic_write_json(metadata_path, metadata)
    except (OSError, csv.Error, ValueError) as exc:
        raise SystemExit(f"Mouse reference build failed: {exc}") from exc

    print(f"Rows: {len(rows)}")
    print(f"Reference: {output_path}")
    print(f"Metadata: {metadata_path}")


if __name__ == "__main__":
    main()
