#!/usr/bin/env python3
"""Build auditable peptide-binding-region (PBR) sequence registries.

The PBR is deliberately distinct from the full AF3 ectodomain and from the
non-contiguous 34-residue pseudosequence.  MHC-I contains one contiguous
alpha1+alpha2 heavy-chain region.  MHC-II contains two independent regions,
alpha1 and beta1, which are never concatenated.

Each boundary is projected within a species/locus from a reviewed UniProt
reference. MHC-II uses the directly annotated Alpha-1/Beta-1 ``REGION``;
MHC-I uses the directly annotated adjacent Alpha-1 and Alpha-2 ``REGION``
coordinates. This avoids applying one raw residue count across unlike loci.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from dataclasses import dataclass
from difflib import SequenceMatcher
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = PROJECT_ROOT / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from mhctransform.fallback import (  # noqa: E402
    compress_positions,
    expand_positions,
    global_alignment_map,
)
from mhctransform.names import exact_name_key  # noqa: E402
from mhctransform.pbr import pbr_definition  # noqa: E402


@dataclass(frozen=True)
class BoundaryReference:
    species: str
    locus: str
    allele: str
    pbr_end_1based: int
    uniprot_accession: str
    source_note: str


def _class_i_reference(
    species: str,
    locus: str,
    allele: str,
    accession: str,
    region_start: int,
    region_end: int,
) -> BoundaryReference:
    pbr_end_1based = region_end - region_start + 1
    return BoundaryReference(
        species=species,
        locus=locus,
        allele=allele,
        pbr_end_1based=pbr_end_1based,
        uniprot_accession=accession,
        source_note=(
            f"UniProtKB {accession}: reviewed adjacent Alpha-1 and Alpha-2 "
            f"REGIONs span {region_start}-{region_end}"
        ),
    )


def _region_reference(
    species: str,
    locus: str,
    allele: str,
    accession: str,
    region_start: int,
    region_end: int,
    region_name: str,
) -> BoundaryReference:
    if region_start < 1 or region_end < region_start:
        raise ValueError(f"invalid UniProt REGION for {locus}")
    return BoundaryReference(
        species=species,
        locus=locus,
        allele=allele,
        pbr_end_1based=region_end - region_start + 1,
        uniprot_accession=accession,
        source_note=(
            f"UniProtKB {accession}: reviewed {region_name} REGION "
            f"{region_start}-{region_end}"
        ),
    )


# The coordinates below are the reviewed UniProt feature coordinates used by
# the strict AF3 builder.  Allele-specific indels are handled by alignment
# projection from the named same-locus AF3 reference, not by reusing the raw
# reference length.
BOUNDARY_REFERENCES: Tuple[BoundaryReference, ...] = (
    _class_i_reference("human", "A", "HLA-A*03:01", "P04439", 25, 206),
    _class_i_reference("human", "B", "HLA-B*07:02", "P01889", 25, 206),
    _class_i_reference("human", "C", "HLA-C*07:02", "P10321", 25, 206),
    _class_i_reference("human", "E", "HLA-E*01:01", "P13747", 22, 203),
    _class_i_reference("human", "F", "HLA-F*01:01", "P30511", 22, 203),
    _class_i_reference("human", "G", "HLA-G*01:01", "P17693", 25, 206),
    _region_reference("human", "DRA", "HLA-DRA*01:01", "P01903", 26, 109, "Alpha-1"),
    _region_reference("human", "DRB1", "HLA-DRB1*15:01", "P01911", 30, 124, "Beta-1"),
    _region_reference("human", "DRB3", "HLA-DRB3*01:01", "P79483", 30, 124, "Beta-1"),
    _region_reference("human", "DRB4", "HLA-DRB4*01:01", "P13762", 30, 124, "Beta-1"),
    _region_reference("human", "DRB5", "HLA-DRB5*01:01", "Q30154", 30, 124, "Beta-1"),
    _region_reference("human", "DQA1", "HLA-DQA1*01:01", "P01909", 24, 119, "Alpha-1"),
    _region_reference("human", "DQB1", "HLA-DQB1*05:01", "P01920", 33, 126, "Beta-1"),
    _region_reference("human", "DPA1", "HLA-DPA1*01:03", "P20036", 29, 115, "Alpha-1"),
    _region_reference("human", "DPB1", "HLA-DPB1*01:01", "P04440", 30, 121, "Beta-1"),
    _class_i_reference("mouse", "H2-D1", "H-2-Db", "P01899", 25, 206),
    _class_i_reference("mouse", "H2-K1", "H-2-Kb", "P01901", 22, 203),
    _class_i_reference("mouse", "H2-L", "H-2-Ld", "P01897", 25, 206),
    _region_reference("mouse", "H2-Aa", "H2-Aa*b", "P14434", 24, 111, "Alpha-1"),
    _region_reference("mouse", "H2-Ab1", "H2-Ab1*b", "P14483", 28, 122, "Beta-1"),
    _region_reference("mouse", "H2-Ea", "H2-Ea*d", "P01904", 26, 109, "Alpha-1"),
    _region_reference("mouse", "H2-Eb1", "H2-Eb1*d", "P01915", 32, 121, "Beta-1"),
)

REFERENCE_BY_KEY = {
    (item.species, item.locus): item for item in BOUNDARY_REFERENCES
}

PBR_FIELDS = (
    "record_id",
    "allele",
    "aliases",
    "species",
    "mhc_class",
    "chain_role",
    "locus",
    "definition_id",
    "sequence",
    "source_af3_record_id",
    "af3_span_1based",
    "source_full_record_id",
    "precursor_span_1based",
    "boundary_method",
    "boundary_source",
    "source_accession",
)

FALLBACK_PBR_FIELDS = (
    "record_id",
    "allele",
    "aliases",
    "species",
    "mhc_class",
    "chain_role",
    "locus",
    "representation",
    "sequence",
    "source_af3_record_id",
    "af3_span_1based",
    "evidence_level",
    "observed_spans_1based",
    "imputed_spans_1based",
    "observed_residue_count",
    "imputed_residue_count",
    "source_record_ids",
    "template_record_id",
    "template_allele",
    "target_pseudo_record_id",
    "pseudo_distance",
    "unknown_pseudo_positions",
    "quality_flags",
)


def _record_id(*parts: str) -> str:
    return hashlib.sha256("\x1f".join(parts).encode("utf-8")).hexdigest()[:20]


def _read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _write_csv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(
        rows,
        key=lambda row: (
            str(row["species"]),
            str(row["mhc_class"]),
            str(row["locus"]),
            str(row["allele"]),
            str(row["record_id"]),
        ),
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(ordered)


def _parse_span(value: str) -> Tuple[int, int]:
    start, end = (int(item) for item in value.split("-", 1))
    if start < 1 or end < start:
        raise ValueError(f"invalid 1-based span: {value}")
    return start, end


def _fast_project_cut(
    reference: str,
    target: str,
    cut: int,
    opcodes: Sequence[Tuple[str, int, int, int, int]] | None = None,
) -> int:
    if opcodes is None:
        opcodes = SequenceMatcher(None, reference, target, autojunk=False).get_opcodes()
    for tag, left_start, left_end, right_start, right_end in opcodes:
        if cut < left_start:
            return right_start
        if left_start <= cut <= left_end:
            if tag == "equal" or (
                tag == "replace"
                and left_end - left_start == right_end - right_start
            ):
                return right_start + cut - left_start
            if tag == "delete":
                return right_start
            if tag == "insert":
                continue
            # Unequal replacement blocks are rare short indels. Project the
            # cut proportionally within that local block, then validate below.
            left_width = left_end - left_start
            right_width = right_end - right_start
            return round(
                right_start + (cut - left_start) * right_width / left_width
            )
    return len(target)


def _project_cut(reference: str, target: str, cut: int) -> int:
    """Project a domain cut through a deterministic same-locus alignment.

    A fast exact-block alignment handles ordinary substitutions. Rare indel
    cases are detected by reciprocal projection and resolved with the full
    Needleman-Wunsch mapping used elsewhere in the reference builder.
    """

    if not 0 < cut < len(reference):
        raise ValueError(f"reference cut {cut} outside sequence length {len(reference)}")
    matcher = SequenceMatcher(None, reference, target, autojunk=False)
    matching = sum(block.size for block in matcher.get_matching_blocks())
    if matching / max(len(reference), len(target)) < 0.55:
        raise ValueError("target is too divergent from its locus PBR reference")
    projected = _fast_project_cut(reference, target, cut, matcher.get_opcodes())
    if _fast_project_cut(target, reference, projected) == cut:
        return projected
    mapping = global_alignment_map(reference, target)
    if cut not in mapping:
        raise ValueError(f"PBR boundary residue {cut} maps to a target gap")
    return mapping[cut]


def _reference_sequences(
    af3_rows: Sequence[Mapping[str, str]],
) -> Mapping[Tuple[str, str], Tuple[BoundaryReference, str]]:
    by_name: MutableMapping[str, List[Mapping[str, str]]] = {}
    for row in af3_rows:
        for name in (row["allele"], *row["aliases"].split(";")):
            if name:
                by_name.setdefault(exact_name_key(name), []).append(row)
    result = {}
    for spec in BOUNDARY_REFERENCES:
        candidates = [
            row
            for row in by_name.get(exact_name_key(spec.allele), ())
            if row["species"] == spec.species and row["locus"] == spec.locus
        ]
        sequences = {row["sequence"] for row in candidates}
        if len(sequences) != 1:
            raise ValueError(
                f"PBR reference {spec.allele} resolves to {len(sequences)} sequences"
            )
        sequence = next(iter(sequences))
        if len(sequence) <= spec.pbr_end_1based:
            raise ValueError(f"PBR reference sequence too short: {spec.allele}")
        result[(spec.species, spec.locus)] = (spec, sequence)
    return result


def _pbr_end(
    row: Mapping[str, str],
    references: Mapping[Tuple[str, str], Tuple[BoundaryReference, str]],
) -> Tuple[int, BoundaryReference]:
    key = (row["species"], row["locus"])
    if key not in references:
        raise ValueError(f"no PBR boundary reference for {key}")
    spec, reference_sequence = references[key]
    end = _project_cut(reference_sequence, row["sequence"], spec.pbr_end_1based)
    minimum = 171 if row["mhc_class"] == "I" else (73 if row["chain_role"] == "alpha" else 90)
    if not minimum <= end < len(row["sequence"]):
        raise ValueError(
            f"implausible PBR boundary for {row['allele']}: {end}/{len(row['sequence'])}"
        )
    return end, spec


def build_strict_pbr_file(data_dir: Path) -> int:
    af3_rows = _read_csv(data_dir / "af3_sequences.csv")
    references = _reference_sequences(af3_rows)
    rows = []
    for row in af3_rows:
        if row["chain_role"] == "companion":
            continue
        end, spec = _pbr_end(row, references)
        precursor_start, _ = _parse_span(row["precursor_span_1based"])
        precursor_end = precursor_start + end - 1
        definition = pbr_definition(row["mhc_class"], row["chain_role"])
        sequence = row["sequence"][:end]
        rows.append(
            {
                "record_id": _record_id(
                    "pbr", definition, row["record_id"], sequence
                ),
                "allele": row["allele"],
                "aliases": row["aliases"],
                "species": row["species"],
                "mhc_class": row["mhc_class"],
                "chain_role": row["chain_role"],
                "locus": row["locus"],
                "definition_id": definition,
                "sequence": sequence,
                "source_af3_record_id": row["record_id"],
                "af3_span_1based": f"1-{end}",
                "source_full_record_id": row["source_full_record_id"],
                "precursor_span_1based": f"{precursor_start}-{precursor_end}",
                "boundary_method": (
                    "same-species/locus sequence-alignment projection of reviewed "
                    "UniProt Alpha-domain REGION boundary"
                ),
                "boundary_source": (
                    f"reference {spec.allele}; {spec.source_note}"
                ),
                "source_accession": row["source_accession"],
            }
        )
    _write_csv(data_dir / "pbr_sequences.csv", PBR_FIELDS, rows)
    return len(rows)


def build_fallback_pbr_file(data_dir: Path) -> int:
    fallback_path = data_dir / "fallback_af3_sequences.csv"
    if not fallback_path.is_file():
        return 0
    af3_rows = _read_csv(data_dir / "af3_sequences.csv")
    references = _reference_sequences(af3_rows)
    rows = []
    for row in _read_csv(fallback_path):
        if row["chain_role"] == "companion":
            continue
        end, _spec = _pbr_end(row, references)
        definition = pbr_definition(row["mhc_class"], row["chain_role"])
        observed = tuple(
            position
            for position in expand_positions(row["observed_spans_1based"])
            if position <= end
        )
        imputed = tuple(
            position
            for position in expand_positions(row["imputed_spans_1based"])
            if position <= end
        )
        if set(observed) & set(imputed):
            raise ValueError(f"overlapping PBR provenance spans: {row['record_id']}")
        if set(observed) | set(imputed) != set(range(1, end + 1)):
            raise ValueError(f"incomplete PBR provenance spans: {row['record_id']}")
        sequence = row["sequence"][:end]
        flags = {
            value for value in row["quality_flags"].split(";") if value
        }
        flags.add("pbr_cropped_from_af3_representation")
        rows.append(
            {
                "record_id": _record_id(
                    "fallback-pbr", definition, row["record_id"], sequence
                ),
                "allele": row["allele"],
                "aliases": row["aliases"],
                "species": row["species"],
                "mhc_class": row["mhc_class"],
                "chain_role": row["chain_role"],
                "locus": row["locus"],
                "representation": definition,
                "sequence": sequence,
                "source_af3_record_id": row["record_id"],
                "af3_span_1based": f"1-{end}",
                "evidence_level": row["evidence_level"],
                "observed_spans_1based": compress_positions(observed),
                "imputed_spans_1based": compress_positions(imputed),
                "observed_residue_count": len(observed),
                "imputed_residue_count": len(imputed),
                "source_record_ids": row["source_record_ids"],
                "template_record_id": row["template_record_id"],
                "template_allele": row["template_allele"],
                "target_pseudo_record_id": row["target_pseudo_record_id"],
                "pseudo_distance": row["pseudo_distance"],
                "unknown_pseudo_positions": row["unknown_pseudo_positions"],
                "quality_flags": ";".join(sorted(flags)),
            }
        )
    _write_csv(
        data_dir / "fallback_pbr_sequences.csv", FALLBACK_PBR_FIELDS, rows
    )
    return len(rows)


def pbr_metadata_fields() -> Dict[str, object]:
    return {
        "pbr_definitions": {
            "mhci-alpha1-alpha2-v1": {
                "chains": ["MHC-I heavy"],
                "description": "Contiguous mature heavy-chain alpha1+alpha2 domains; B2M is not included.",
            },
            "mhcii-alpha1-v1": {
                "chains": ["MHC-II alpha"],
                "description": "Contiguous mature MHC-II alpha1 domain.",
            },
            "mhcii-beta1-v1": {
                "chains": ["MHC-II beta"],
                "description": "Contiguous mature MHC-II beta1 domain.",
            },
            "boundary_policy": (
                    "Same-species/locus alignment projection from reviewed UniProt "
                    "Alpha-1/Beta-1 REGION features for MHC-II and adjacent "
                    "Alpha-1 plus Alpha-2 REGION features for MHC-I."
            ),
            "reference_records": [
                {
                    "species": item.species,
                    "locus": item.locus,
                    "allele": item.allele,
                    "pbr_end_1based": item.pbr_end_1based,
                    "uniprot_accession": item.uniprot_accession,
                }
                for item in BOUNDARY_REFERENCES
            ],
        }
    }


def update_pbr_metadata(data_dir: Path, strict_count: int, fallback_count: int) -> None:
    path = data_dir / "metadata.json"
    if not path.is_file():
        return
    with path.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    metadata.update(pbr_metadata_fields())
    metadata.setdefault("counts", {}).update(
        {"pbr_total": strict_count, "fallback_pbr": fallback_count}
    )
    with path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def build_pbr_files(data_dir: Path) -> Mapping[str, int]:
    strict_count = build_strict_pbr_file(data_dir)
    fallback_count = build_fallback_pbr_file(data_dir)
    update_pbr_metadata(data_dir, strict_count, fallback_count)
    return {"pbr_total": strict_count, "fallback_pbr": fallback_count}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=PROJECT_ROOT / "src" / "mhctransform" / "data",
    )
    args = parser.parse_args()
    counts = build_pbr_files(args.data_dir)
    print(json.dumps(counts, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
