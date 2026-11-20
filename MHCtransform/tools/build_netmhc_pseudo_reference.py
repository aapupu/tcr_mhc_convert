#!/usr/bin/env python3
"""Build the NetMHC-supported and sequence-derived 34-aa registry.

The public NetMHCpan/NetMHCIIpan allele-name lists define *supported names*.
For names with a strict PBR, contact residues are independently extracted from
the pinned IPD/UniProt-derived sequence.  Public DTU training-data pseudo files
are supplied only as local validation oracles: neither the archives nor rows
available only from those files are copied into the distributable registry.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date
from difflib import SequenceMatcher
from itertools import product
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = PROJECT_ROOT / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from mhctransform.fallback import locate_contact_positions  # noqa: E402


PSEUDO_I_DEFINITION = "netmhcpan-34-v1"
PSEUDO_II_DEFINITION = "netmhciipan-34-v1"
OFFICIAL_I_VERSION = "NetMHCpan-4.2f"
OFFICIAL_II_VERSION = "NetMHCIIpan-4.3k"
OFFICIAL_I_LIST_URL = (
    "https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/allele.list"
)
OFFICIAL_II_LIST_URL = (
    "https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/alleles_name.txt"
)
OFFICIAL_I_LIST_SHA256 = (
    "09167ed2e1f372ae732b978e36e7000a15c85cef2d9eae6c484012f3467b8e18"
)
OFFICIAL_II_LIST_SHA256 = (
    "143930820fcc5d395c78280bb82a07a6c0ee943db2fc3f7ab0f7a78b1cacfcb3"
)
OFFICIAL_I_TRAINING_PSEUDO_SHA256 = (
    "bb8fe764cd187e55d7f83a465eef152419151ad9924825c411ba4e8ed84ef5a7"
)
OFFICIAL_II_TRAINING_PSEUDO_SHA256 = (
    "daf3b7d693e76d8e8102c374d130467ff1ef7f9311064a730d56a753a3d514b5"
)

PSEUDO_FIELDS = (
    "record_id",
    "canonical_name",
    "aliases",
    "species",
    "mhc_class",
    "definition_id",
    "pseudo_sequence",
    "alpha_allele",
    "beta_allele",
    "source",
    "source_version",
    "evidence_level",
    "official_support",
    "validation_status",
    "source_record_ids",
    "quality_flags",
)
CONFLICT_FIELDS = PSEUDO_FIELDS + ("reason",)
COMPONENT_FIELDS = (
    "record_id",
    "allele",
    "aliases",
    "species",
    "mhc_class",
    "chain_role",
    "locus",
    "definition_id",
    "pseudo_segment",
    "segment_start_1based",
    "source",
    "source_version",
    "evidence_level",
    "official_support",
    "validation_status",
    "source_record_ids",
    "quality_flags",
)

MHC_I_POSITIONS = (
    7, 9, 24, 45, 59, 62, 63, 66, 67, 69, 70, 73, 74, 76, 77,
    80, 81, 84, 95, 97, 99, 114, 116, 118, 143, 147, 150, 152,
    156, 158, 159, 163, 167, 171,
)

# Published MHC-I contact positions are alignment/reference coordinates, not
# safe raw string offsets for every allele.  Rare in-domain indels therefore
# have to be projected from an ordinary same-locus reference before residues
# are selected.
MHC_I_ANCHORS: Mapping[Tuple[str, str], str] = {
    ("human", "A"): "HLA-A*03:01",
    ("human", "B"): "HLA-B*07:02",
    ("human", "C"): "HLA-C*07:02",
    ("human", "E"): "HLA-E*01:01",
    ("human", "F"): "HLA-F*01:01",
    ("human", "G"): "HLA-G*01:01",
    ("mouse", "H2-D1"): "H-2-Db",
    ("mouse", "H2-K1"): "H-2-Kb",
    ("mouse", "H2-L"): "H-2-Ld",
}

# Positions are 1-based within the PBR.  DQ/DP values were projected from the
# published DR reference positions through locus-matched structural sequence
# alignments and validated against the 2026-07-06 public 4.3 training artifact.
MHC_II_CONTACTS: Mapping[str, Tuple[str, Tuple[int, ...]]] = {
    "DRA": (
        "HLA-DRA*01:01",
        (9, 11, 22, 24, 31, 52, 53, 58, 59, 61, 65, 66, 68, 72, 73),
    ),
    "DQA1": (
        "HLA-DQA1*01:01",
        (11, 14, 25, 27, 34, 55, 56, 61, 62, 64, 68, 69, 71, 75, 76),
    ),
    "DPA1": (
        "HLA-DPA1*01:03",
        (12, 14, 25, 27, 34, 55, 56, 61, 62, 64, 68, 69, 71, 75, 76),
    ),
    "DRB1": (
        "HLA-DRB1*15:01",
        (9, 11, 13, 26, 28, 30, 47, 57, 67, 70, 71, 74, 77, 78, 81, 85, 86, 89, 90),
    ),
    "DRB3": (
        "HLA-DRB3*01:01",
        (9, 11, 13, 26, 28, 30, 47, 57, 67, 70, 71, 74, 77, 78, 81, 85, 86, 89, 90),
    ),
    "DRB4": (
        "HLA-DRB4*01:01",
        (9, 11, 13, 26, 28, 30, 47, 57, 67, 70, 71, 74, 77, 78, 81, 85, 86, 89, 90),
    ),
    "DRB5": (
        "HLA-DRB5*01:01",
        (9, 11, 13, 26, 28, 30, 47, 57, 67, 70, 71, 74, 77, 78, 81, 85, 86, 89, 90),
    ),
    "DQB1": (
        "HLA-DQB1*02:01",
        (9, 11, 13, 26, 28, 30, 47, 57, 67, 70, 71, 74, 77, 78, 81, 85, 86, 89, 90),
    ),
    "DPB1": (
        "HLA-DPB1*01:01",
        (9, 11, 13, 24, 26, 28, 45, 55, 65, 68, 69, 72, 75, 76, 79, 83, 84, 87, 88),
    ),
}

MOUSE_PAIRS: Mapping[str, Tuple[str, str]] = {
    "H-2-IAb": ("H2-Aa*b", "H2-Ab1*b"),
    "H-2-IAd": ("H2-Aa*d", "H2-Ab1*d"),
    "H-2-IAk": ("H2-Aa*k", "H2-Ab1*k"),
    "H-2-IAq": ("H2-Aa*q", "H2-Ab1*q"),
    "H-2-IAs": ("H2-Aa*s", "H2-Ab1*s"),
    "H-2-IAu": ("H2-Aa*u", "H2-Ab1*u"),
    "H-2-IEd": ("H2-Ea*d", "H2-Eb1*d"),
    "H-2-IEk": ("H2-Ea*k", "H2-Eb1*k"),
}

# NetMHC's public names predate two current IPD expression suffixes.  Preserve
# the present IPD canonical name and expose the historical unsuffixed form as
# an alias; never erase the biologically meaningful Q suffix.
OFFICIAL_TO_CURRENT_COMPONENT: Mapping[str, str] = {
    "HLA-DQA1*01:07": "HLA-DQA1*01:07Q",
    "HLA-DRB1*15:13": "HLA-DRB1*15:13Q",
}
CURRENT_TO_OFFICIAL_COMPONENT = {
    current: official for official, current in OFFICIAL_TO_CURRENT_COMPONENT.items()
}


def current_component_name(name: str) -> str:
    return OFFICIAL_TO_CURRENT_COMPONENT.get(name, name)


def official_component_name(name: str) -> str:
    return CURRENT_TO_OFFICIAL_COMPONENT.get(name, name)


@dataclass(frozen=True)
class ComponentValue:
    sequence: str
    evidence: str
    source: str
    source_version: str
    source_ids: Tuple[str, ...]
    validation: str
    official_support: bool
    aliases: Tuple[str, ...] = ()
    quality_flags: Tuple[str, ...] = ()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record_id(*parts: object) -> str:
    return hashlib.sha256(
        "\x1f".join(str(part) for part in parts).encode("utf-8")
    ).hexdigest()[:20]


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def split_values(value: str) -> Tuple[str, ...]:
    return tuple(item for item in str(value).split(";") if item)


def contact_alignment_map(reference: str, target: str) -> Mapping[int, int]:
    """Conservatively map aligned PBR coordinates without inventing gaps.

    Equal-length IPD PBR chains use the published mature-chain coordinates
    directly.  For a genuine length difference, only exact/equal-width
    SequenceMatcher blocks are mapped; a contact inside a deletion or an
    unequal replacement remains unmapped and the derived pseudo is rejected.
    """

    if len(reference) == len(target):
        return {position: position for position in range(1, len(reference) + 1)}
    mapping: Dict[int, int] = {}
    matcher = SequenceMatcher(None, reference, target, autojunk=False)
    for tag, left_start, left_end, right_start, right_end in matcher.get_opcodes():
        if tag == "equal" or (
            tag == "replace"
            and left_end - left_start == right_end - right_start
        ):
            for offset in range(left_end - left_start):
                mapping[left_start + offset + 1] = right_start + offset + 1
    return mapping


def write_csv(
    path: Path,
    fields: Sequence[str],
    rows: Iterable[Mapping[str, object]],
    *,
    key_fields: Sequence[str],
) -> None:
    ordered = sorted(
        rows,
        key=lambda row: tuple(str(row.get(field, "")) for field in key_fields),
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in ordered:
            writer.writerow({field: row.get(field, "") for field in fields})


def canonical_hla_component(gene: str, digits: str) -> str:
    if ":" in digits:
        first, second = digits.split(":", 1)
    else:
        if len(digits) < 4:
            raise ValueError(f"invalid compact HLA field: {gene}{digits}")
        first, second = digits[:-2], digits[-2:]
    return f"HLA-{gene}*{first}:{second}"


def parse_official_i_name(raw: str) -> str | None:
    match = re.fullmatch(r"HLA-([ABCEFG])([0-9]+):([0-9]+)", raw)
    if match:
        return f"HLA-{match.group(1)}*{match.group(2)}:{match.group(3)}"
    if re.fullmatch(r"H-2-[DKL][A-Za-z0-9]+", raw):
        return raw
    return None


def parse_official_ii_name(raw: str) -> Tuple[str, str, str] | None:
    dr = re.fullmatch(r"(DRB[1-5])_([0-9]+)", raw)
    if dr:
        beta = current_component_name(
            canonical_hla_component(dr.group(1), dr.group(2))
        )
        alpha = "HLA-DRA*01:01"
        return f"{alpha}/{beta}", alpha, beta
    pair = re.fullmatch(
        r"HLA-(D[QP]A1)([0-9]+)-(D[QP]B1)([0-9]+)", raw
    )
    if pair:
        alpha = current_component_name(
            canonical_hla_component(pair.group(1), pair.group(2))
        )
        beta = current_component_name(
            canonical_hla_component(pair.group(3), pair.group(4))
        )
        return f"{alpha}/{beta}", alpha, beta
    if raw in MOUSE_PAIRS:
        alpha, beta = MOUSE_PAIRS[raw]
        return raw, alpha, beta
    return None


def parse_training_pseudo_i(path: Path) -> Mapping[str, Tuple[str, str]]:
    result: Dict[str, Tuple[str, str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) < 2:
            continue
        name = parse_official_i_name(fields[0])
        sequence = fields[1].upper()
        if name and len(sequence) == 34:
            if name in result and result[name][0] != sequence:
                raise ValueError(f"conflicting official MHC-I pseudo: {name}")
            result[name] = (sequence, fields[0])
    return result


def parse_training_pseudo_ii(
    path: Path,
) -> Mapping[str, Tuple[str, str, str, str]]:
    result: Dict[str, Tuple[str, str, str, str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) < 2:
            continue
        parsed = parse_official_ii_name(fields[0])
        sequence = fields[1].upper()
        if parsed and len(sequence) == 34:
            name, alpha, beta = parsed
            if name in result and result[name][0] != sequence:
                raise ValueError(f"conflicting official MHC-II pseudo: {name}")
            result[name] = (sequence, alpha, beta, fields[0])
    return result


def parse_supported_i(path: Path) -> Tuple[set[str], Mapping[str, set[str]]]:
    if sha256_path(path) != OFFICIAL_I_LIST_SHA256:
        raise ValueError("NetMHCpan 4.2 allele.list hash changed")
    names: set[str] = set()
    aliases: MutableMapping[str, set[str]] = defaultdict(set)
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if not fields:
            continue
        name = parse_official_i_name(fields[0].strip())
        if name:
            names.add(name)
            aliases[name].add(fields[0].strip())
    return names, aliases


def parse_supported_ii(path: Path) -> Mapping[str, set[str]]:
    if sha256_path(path) != OFFICIAL_II_LIST_SHA256:
        raise ValueError("NetMHCIIpan 4.3 alleles_name.txt hash changed")
    result: MutableMapping[str, set[str]] = defaultdict(set)
    starts = (0, 16, 32, 48, 64, 80, 94, 110)
    columns = ("dr_beta", "dq_alpha", "dq_beta", "dp_alpha", "dp_beta", "mouse", "bola")
    lines = path.read_text(encoding="utf-8").splitlines()
    for line in lines[1:]:
        padded = line.ljust(110)
        for label, start, end in zip(columns, starts, starts[1:]):
            value = padded[start:end].strip()
            if not value or label == "bola":
                continue
            if label == "mouse":
                result[label].add(value)
            else:
                result[label].add(current_component_name(f"HLA-{value}"))
    return result


def aliases_for_chain(name: str, official_aliases: Iterable[str] = ()) -> Tuple[str, ...]:
    aliases = set(official_aliases)
    if name.startswith("HLA-"):
        aliases.add(name.removeprefix("HLA-"))
        aliases.add(name.replace("*", ""))
        official_name = official_component_name(name)
        if official_name != name:
            aliases.update(
                {
                    official_name,
                    official_name.removeprefix("HLA-"),
                    official_name.replace("*", ""),
                }
            )
    aliases.discard(name)
    return tuple(sorted(aliases))


def secondary_component_values(
    rows: Sequence[Mapping[str, str]],
) -> Mapping[str, Tuple[str, Tuple[str, ...]]]:
    segments: MutableMapping[str, List[str]] = defaultdict(list)
    aliases: MutableMapping[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["mhc_class"] != "II":
            continue
        segments[row["alpha_allele"]].append(row["pseudo_sequence"][:15])
        segments[row["beta_allele"]].append(row["pseudo_sequence"][15:])
        aliases[row["alpha_allele"]].add(row["alpha_allele"])
        aliases[row["beta_allele"]].add(row["beta_allele"])
    result = {}
    for allele, values in segments.items():
        counts = Counter(values)
        best = counts.most_common()
        if len(best) > 1 and best[0][1] == best[1][1]:
            continue
        result[allele] = (best[0][0], tuple(sorted(aliases[allele])))
    return result


def derive_from_strict_pbr(
    pbr_rows: Sequence[Mapping[str, str]],
) -> Tuple[Mapping[str, ComponentValue], Mapping[str, ComponentValue]]:
    by_name = {row["allele"]: row for row in pbr_rows}
    derived_i: Dict[str, ComponentValue] = {}
    derived_ii: Dict[str, ComponentValue] = {}
    for row in pbr_rows:
        sequence = row["sequence"]
        if row["mhc_class"] == "I" and row["chain_role"] == "heavy":
            anchor_name = MHC_I_ANCHORS.get((row["species"], row["locus"]))
            anchor = by_name.get(anchor_name or "")
            if anchor is None:
                raise ValueError(
                    f"missing MHC-I PBR contact anchor for "
                    f"{row['species']} {row['locus']}"
                )
            mapping = contact_alignment_map(anchor["sequence"], sequence)
            if any(position not in mapping for position in MHC_I_POSITIONS):
                continue
            pseudo = "".join(
                sequence[mapping[position] - 1] for position in MHC_I_POSITIONS
            )
            if "X" in pseudo:
                continue
            derived_i[row["allele"]] = ComponentValue(
                pseudo,
                "derived_exact_pbr",
                "pinned IPD/UniProt strict PBR contact extraction",
                "IPD-IMGT/HLA-3.66.0",
                (row["record_id"],),
                "derived_from_exact_pbr",
                False,
                split_values(row["aliases"]),
            )
            continue
        if row["mhc_class"] != "II" or row["locus"] not in MHC_II_CONTACTS:
            continue
        anchor_name, positions = MHC_II_CONTACTS[row["locus"]]
        anchor = by_name.get(anchor_name)
        if anchor is None:
            raise ValueError(f"missing PBR contact anchor {anchor_name}")
        mapping = contact_alignment_map(anchor["sequence"], sequence)
        if any(position not in mapping for position in positions):
            continue
        pseudo = "".join(sequence[mapping[position] - 1] for position in positions)
        if "X" in pseudo:
            continue
        derived_ii[row["allele"]] = ComponentValue(
            pseudo,
            "derived_exact_pbr",
            "pinned IPD/UniProt strict PBR contact extraction",
            "IPD-IMGT/HLA-3.66.0",
            (row["record_id"],),
            "derived_from_exact_pbr",
            False,
            split_values(row["aliases"]),
        )
    return derived_i, derived_ii


def derive_from_deposited_sequences(
    full_rows: Sequence[Mapping[str, str]],
    pbr_rows: Sequence[Mapping[str, str]],
    *,
    target_i: set[str],
    target_ii: set[str],
) -> Tuple[
    Mapping[str, ComponentValue],
    Mapping[str, ComponentValue],
    Mapping[str, Tuple[str, ...]],
]:
    """Extract contacts from any deposited record that covers every position.

    Local alignment is used because IPD partial protein records commonly start
    inside the mature PBR (for example at PBR residue 2) and omit signal/TM
    context.  A row is accepted only when every published contact position maps
    to a concrete non-X deposited residue.  No residue is imputed.
    """

    pbr_by_name = {row["allele"]: row for row in pbr_rows}
    candidates_i: MutableMapping[str, List[ComponentValue]] = defaultdict(list)
    candidates_ii: MutableMapping[str, List[ComponentValue]] = defaultdict(list)
    for row in full_rows:
        if row["chain_role"] == "companion" or not row["sequence"]:
            continue
        if row["mhc_class"] == "I" and row["allele"] not in target_i:
            continue
        if row["mhc_class"] == "II" and row["allele"] not in target_ii:
            continue
        sequence = row["sequence"]
        if row["mhc_class"] == "I" and row["chain_role"] == "heavy":
            anchor_name = MHC_I_ANCHORS.get((row["species"], row["locus"]))
            anchor = pbr_by_name.get(anchor_name or "")
            if anchor is None:
                continue
            positions = MHC_I_POSITIONS
        elif row["mhc_class"] == "II" and row["locus"] in MHC_II_CONTACTS:
            anchor_name, positions = MHC_II_CONTACTS[row["locus"]]
            anchor = pbr_by_name.get(anchor_name)
            if anchor is None:
                continue
        else:
            continue
        matcher = SequenceMatcher(None, anchor["sequence"], sequence, autojunk=False)
        # Partial records may contain little beyond the groove (DPB1*113:01
        # is a validated boundary case at ~0.549).  The permissive homology
        # screen is followed by stricter requirements: same locus, all contact
        # columns mapped, no X, and cross-check against the training oracle.
        if matcher.ratio() < 0.50:
            continue
        mapping: Dict[int, int] = {}
        for tag, left_start, left_end, right_start, right_end in matcher.get_opcodes():
            if tag == "equal" or (
                tag == "replace"
                and left_end - left_start == right_end - right_start
            ):
                for offset in range(left_end - left_start):
                    mapping[left_start + offset + 1] = right_start + offset + 1
        if any(position not in mapping for position in positions):
            continue
        residues = tuple(sequence[mapping[position] - 1] for position in positions)
        if any(residue == "X" for residue in residues):
            continue
        pseudo = "".join(residues)
        exact_full = row.get("is_complete", "false").lower() == "true"
        evidence = (
            "derived_contact_complete_full"
            if exact_full
            else "derived_contact_complete_partial"
        )
        value = ComponentValue(
            pseudo,
            evidence,
            "pinned IPD/UniProt deposited contact-complete protein",
            "IPD-IMGT/HLA-3.66.0",
            (row["record_id"],),
            evidence,
            False,
            split_values(row["aliases"]),
        )
        (candidates_i if row["mhc_class"] == "I" else candidates_ii)[
            row["allele"]
        ].append(value)

    conflicts: Dict[str, Tuple[str, ...]] = {}

    def collapse(
        values: Mapping[str, Sequence[ComponentValue]], label: str
    ) -> Dict[str, ComponentValue]:
        result = {}
        for allele, items in values.items():
            sequences = {item.sequence for item in items}
            if len(sequences) != 1:
                conflicts[f"{label}:{allele}"] = tuple(sorted(sequences))
                continue
            sequence = next(iter(sequences))
            matching = [item for item in items if item.sequence == sequence]
            matching.sort(
                key=lambda item: (
                    0 if item.evidence == "derived_contact_complete_full" else 1,
                    item.source_ids,
                )
            )
            chosen = matching[0]
            result[allele] = ComponentValue(
                chosen.sequence,
                chosen.evidence,
                chosen.source,
                chosen.source_version,
                tuple(sorted({value for item in matching for value in item.source_ids})),
                chosen.validation,
                chosen.official_support,
                tuple(sorted({value for item in matching for value in item.aliases})),
            )
        return result

    return collapse(candidates_i, "I"), collapse(candidates_ii, "II"), conflicts


def recover_official_validated_mhci_indels(
    full_rows: Sequence[Mapping[str, str]],
    pbr_rows: Sequence[Mapping[str, str]],
    official_i: Mapping[str, Tuple[str, str]],
    unavailable_names: set[str],
) -> Mapping[str, ComponentValue]:
    """Resolve internal indels only when every oracle residue is observed.

    The oracle is used to select alignment columns, never to supply a residue.
    Terminally missing contact regions remain unavailable (for example H-2-Kq).
    """

    pbr_by_name = {row["allele"]: row for row in pbr_rows}
    recovered: Dict[str, ComponentValue] = {}
    for row in full_rows:
        name = row["allele"]
        if (
            name not in unavailable_names
            or name not in official_i
            or row["mhc_class"] != "I"
            or row["chain_role"] != "heavy"
        ):
            continue
        official_sequence = official_i[name][0]
        if "X" in official_sequence:
            continue
        anchor_name = MHC_I_ANCHORS.get((row["species"], row["locus"]))
        anchor = pbr_by_name.get(anchor_name or "")
        if anchor is None:
            continue
        conservative = contact_alignment_map(anchor["sequence"], row["sequence"])
        if not conservative:
            continue
        # Reject terminal truncations even if a coincidental nearby residue
        # could satisfy the oracle.  Only internal alignment ambiguity may be
        # resolved this way.
        if min(conservative) > min(MHC_I_POSITIONS) or max(conservative) < max(
            MHC_I_POSITIONS
        ):
            continue
        try:
            observed_positions = locate_contact_positions(
                anchor["sequence"],
                row["sequence"],
                MHC_I_POSITIONS,
                official_sequence,
            )
        except ValueError:
            continue
        observed = "".join(
            row["sequence"][position - 1] for position in observed_positions
        )
        if observed != official_sequence:
            continue
        evidence = (
            "derived_contact_complete_full"
            if row.get("is_complete", "false").lower() == "true"
            else "derived_contact_complete_partial"
        )
        recovered[name] = ComponentValue(
            observed,
            evidence,
            "pinned IPD/UniProt deposited protein with official alignment validation",
            "IPD-IMGT/HLA-3.66.0",
            (row["record_id"],),
            "official_training_exact",
            True,
            split_values(row["aliases"]),
            ("official_oracle_resolved_internal_indel_columns",),
        )
    return recovered


def locus_for_component(name: str) -> str:
    if name.startswith("HLA-"):
        return name.removeprefix("HLA-").split("*", 1)[0]
    return name.split("*", 1)[0]


def role_for_component(name: str) -> str:
    locus = locus_for_component(name).upper()
    return "alpha" if locus in {"DRA", "DQA1", "DPA1", "H2-AA", "H2-EA"} or name.startswith(("H2-Aa", "H2-Ea")) else "beta"


def component_row(
    allele: str,
    value: ComponentValue,
    supported_components: Mapping[str, set[str]],
) -> Dict[str, object]:
    role = role_for_component(allele)
    locus = locus_for_component(allele)
    official = value.official_support or any(
        allele in values for values in supported_components.values()
    )
    aliases = set(value.aliases) | set(aliases_for_chain(allele))
    aliases.discard(allele)
    quality_flags = set(value.quality_flags)
    if allele in CURRENT_TO_OFFICIAL_COMPONENT:
        quality_flags.add("official_name_maps_to_current_q_suffix")
    return {
        "record_id": record_id("pseudo-component", allele, value.sequence),
        "allele": allele,
        "aliases": ";".join(sorted(aliases)),
        "species": "human" if allele.startswith("HLA-") else "mouse",
        "mhc_class": "II",
        "chain_role": role,
        "locus": locus,
        "definition_id": PSEUDO_II_DEFINITION,
        "pseudo_segment": value.sequence,
        "segment_start_1based": 1 if role == "alpha" else 16,
        "source": value.source,
        "source_version": value.source_version,
        "evidence_level": value.evidence,
        "official_support": str(official).lower(),
        "validation_status": value.validation,
        "source_record_ids": ";".join(value.source_ids),
        "quality_flags": ";".join(sorted(quality_flags)),
    }


def normalized_base_row(row: Mapping[str, str]) -> Dict[str, object]:
    return {
        **{field: row.get(field, "") for field in PSEUDO_FIELDS},
        "evidence_level": row.get("evidence_level", "legacy_secondary"),
        "official_support": row.get("official_support", "false"),
        "validation_status": row.get("validation_status", "legacy_secondary_only"),
        "source_record_ids": row.get("source_record_ids", ""),
        "quality_flags": row.get("quality_flags", "legacy_secondary_reference"),
    }


def build(args: argparse.Namespace) -> Mapping[str, object]:
    data_dir = args.data_dir.resolve()
    secondary_rows = read_csv(data_dir / "secondary_pseudo_sequences.csv")
    pbr_rows = read_csv(data_dir / "pbr_sequences.csv")
    full_rows = read_csv(data_dir / "full_sequences.csv")
    supported_i, supported_i_aliases = parse_supported_i(args.supported_i)
    supported_ii = parse_supported_ii(args.supported_ii)
    official_i = parse_training_pseudo_i(args.official_pseudo_i)
    official_ii = parse_training_pseudo_ii(args.official_pseudo_ii)
    official_i_hash = sha256_path(args.official_pseudo_i)
    official_ii_hash = sha256_path(args.official_pseudo_ii)
    if official_i_hash != OFFICIAL_I_TRAINING_PSEUDO_SHA256:
        raise ValueError("NetMHCpan 4.2 training pseudosequence hash changed")
    if official_ii_hash != OFFICIAL_II_TRAINING_PSEUDO_SHA256:
        raise ValueError("NetMHCIIpan 4.3 training pseudosequence hash changed")
    derived_i, derived_ii = derive_from_strict_pbr(pbr_rows)
    secondary_i_names = {
        row["canonical_name"] for row in secondary_rows if row["mhc_class"] == "I"
    }
    secondary_ii_components = {
        value
        for row in secondary_rows
        if row["mhc_class"] == "II"
        for value in (row["alpha_allele"], row["beta_allele"])
    }
    official_ii_components = {
        value
        for _sequence, alpha, beta, _raw in official_ii.values()
        for value in (alpha, beta)
    }
    supported_ii_components = set().union(*supported_ii.values())
    supported_ii_components.add("HLA-DRA*01:01")
    deposited_i, deposited_ii, deposited_conflicts = derive_from_deposited_sequences(
        full_rows,
        pbr_rows,
        target_i=set(supported_i) | set(official_i) | secondary_i_names,
        target_ii=(
            supported_ii_components
            | official_ii_components
            | secondary_ii_components
        ),
    )
    official_guided_i = recover_official_validated_mhci_indels(
        full_rows,
        pbr_rows,
        official_i,
        set(official_i) - set(derived_i) - set(deposited_i),
    )
    deposited_i = {**deposited_i, **official_guided_i}

    secondary_i: MutableMapping[str, List[Mapping[str, str]]] = defaultdict(list)
    secondary_ii: MutableMapping[str, List[Mapping[str, str]]] = defaultdict(list)
    for row in secondary_rows:
        (secondary_i if row["mhc_class"] == "I" else secondary_ii)[
            row["canonical_name"]
        ].append(row)

    conflict_rows: List[Dict[str, object]] = []
    validation_conflicts_i = []
    final_i: List[Dict[str, object]] = []
    all_i_names = (
        set(supported_i)
        | set(official_i)
        | set(derived_i)
        | set(deposited_i)
        | set(secondary_i)
    )
    for name in sorted(all_i_names):
        chosen_sequence = ""
        source = ""
        source_version = ""
        evidence = ""
        validation = ""
        source_ids: Tuple[str, ...] = ()
        selected_quality_flags: Tuple[str, ...] = ()
        aliases = set(supported_i_aliases.get(name, ()))
        if name in derived_i:
            value = derived_i[name]
            chosen_sequence = value.sequence
            source = value.source
            source_version = value.source_version
            evidence = value.evidence
            validation = value.validation
            source_ids = value.source_ids
            selected_quality_flags = value.quality_flags
        elif name in deposited_i:
            value = deposited_i[name]
            chosen_sequence = value.sequence
            source = value.source
            source_version = value.source_version
            evidence = value.evidence
            validation = value.validation
            source_ids = value.source_ids
            selected_quality_flags = value.quality_flags
        elif secondary_i.get(name):
            rows = secondary_i[name]
            sequences = {row["pseudo_sequence"] for row in rows}
            if len(sequences) != 1:
                continue
            chosen_sequence = next(iter(sequences))
            source = rows[0]["source"]
            source_version = rows[0]["source_version"]
            evidence = "legacy_secondary"
            validation = "legacy_secondary_only"
        else:
            continue
        if name in official_i:
            official_sequence, raw_alias = official_i[name]
            aliases.add(raw_alias)
            if chosen_sequence != official_sequence:
                validation_conflicts_i.append(name)
                # The public training artifact is a validation oracle, not a
                # redistributable source database.  Do not publish a row that
                # fails that check and do not copy the oracle sequence in.
                continue
            validation = "official_training_exact"
        for old in secondary_i.get(name, ()):
            aliases.update(split_values(old["aliases"]))
            if old["pseudo_sequence"] != chosen_sequence:
                loser = normalized_base_row(old)
                loser["reason"] = (
                    "secondary pseudo conflicts with a higher-priority "
                    "independently derived PBR/contact-complete reference"
                )
                conflict_rows.append(loser)
        aliases.update(aliases_for_chain(name))
        aliases.discard(name)
        species = "human" if name.startswith("HLA-") else "mouse"
        row = {
            "record_id": record_id("pseudo", PSEUDO_I_DEFINITION, name, chosen_sequence),
            "canonical_name": name,
            "aliases": ";".join(sorted(aliases)),
            "species": species,
            "mhc_class": "I",
            "definition_id": PSEUDO_I_DEFINITION,
            "pseudo_sequence": chosen_sequence,
            "alpha_allele": name,
            "beta_allele": "",
            "source": source,
            "source_version": source_version,
            "evidence_level": evidence,
            "official_support": str(name in supported_i).lower(),
            "validation_status": validation,
            "source_record_ids": ";".join(source_ids),
            "quality_flags": ";".join(
                sorted(
                    {
                        *selected_quality_flags,
                        *(
                            ("outside_current_netmhcpan_support_list",)
                            if name not in supported_i
                            else ()
                        ),
                    }
                )
            ),
        }
        final_i.append(row)

    official_components: MutableMapping[str, set[str]] = defaultdict(set)
    for _name, (sequence, alpha, beta, _raw_alias) in official_ii.items():
        official_components[alpha].add(sequence[:15])
        official_components[beta].add(sequence[15:])
    unresolved_official_components = {
        allele: values for allele, values in official_components.items() if len(values) != 1
    }
    if unresolved_official_components:
        raise ValueError(
            f"official MHC-II component conflicts: {unresolved_official_components}"
        )

    secondary_components = secondary_component_values(secondary_rows)
    component_values: Dict[str, ComponentValue] = {}
    component_validation_conflicts = []
    component_secondary_conflicts = []
    all_component_names = (
        set(official_components)
        | set(derived_ii)
        | set(deposited_ii)
        | set(secondary_components)
    )
    supported_component_names = set().union(*supported_ii.values())
    supported_component_names.add("HLA-DRA*01:01")
    supported_component_names.update(
        component
        for molecule in supported_ii["mouse"]
        if molecule in MOUSE_PAIRS
        for component in MOUSE_PAIRS[molecule]
    )
    for allele in sorted(all_component_names):
        if allele in derived_ii:
            derived = derived_ii[allele]
            value = ComponentValue(
                derived.sequence,
                derived.evidence,
                derived.source,
                derived.source_version,
                derived.source_ids,
                derived.validation,
                allele in supported_component_names,
                derived.aliases,
            )
        elif allele in deposited_ii:
            deposited = deposited_ii[allele]
            value = ComponentValue(
                deposited.sequence,
                deposited.evidence,
                deposited.source,
                deposited.source_version,
                deposited.source_ids,
                deposited.validation,
                allele in supported_component_names,
                deposited.aliases,
            )
        elif allele in secondary_components:
            sequence, aliases = secondary_components[allele]
            value = ComponentValue(
                sequence,
                "legacy_secondary",
                "PUFFIN NetMHCIIpan-format secondary mirror",
                "pinned-secondary-snapshot",
                (),
                "legacy_secondary_only",
                allele in supported_component_names,
                aliases,
            )
        else:
            # An official-training-only component is intentionally not copied
            # into the distributable registry.  It remains part of the local
            # validation report and can be revisited when a primary sequence
            # covers every contact column.
            continue

        if allele in official_components:
            official_sequence = next(iter(official_components[allele]))
            if value.sequence != official_sequence:
                component_validation_conflicts.append(
                    {
                        "allele": allele,
                        "official": official_sequence,
                        "derived": value.sequence,
                    }
                )
                secondary = secondary_components.get(allele)
                if secondary is None or secondary[0] != official_sequence:
                    continue
                # Preserve an already-bundled legacy mirror only when it is an
                # exact match to the validation oracle; never copy the oracle.
                sequence, aliases = secondary
                value = ComponentValue(
                    sequence,
                    "legacy_secondary",
                    "PUFFIN NetMHCIIpan-format secondary mirror",
                    "pinned-secondary-snapshot",
                    (),
                    "official_training_exact",
                    allele in supported_component_names,
                    aliases,
                    ("conflicts_with_current_ipd_derived_sequence",),
                )
            else:
                value = ComponentValue(
                    value.sequence,
                    value.evidence,
                    value.source,
                    value.source_version,
                    value.source_ids,
                    "official_training_exact",
                    allele in supported_component_names,
                    value.aliases,
                    value.quality_flags,
                )
        if allele in secondary_components:
            secondary_sequence = secondary_components[allele][0]
            if secondary_sequence != value.sequence:
                component_secondary_conflicts.append(
                    {
                        "allele": allele,
                        "selected": value.sequence,
                        "secondary": secondary_sequence,
                    }
                )
        component_values[allele] = value

    component_rows = [
        component_row(allele, value, supported_ii)
        for allele, value in component_values.items()
    ]
    component_ids = {row["allele"]: row["record_id"] for row in component_rows}

    official_pairs: Dict[str, Tuple[str, str, Tuple[str, ...]]] = {}
    dra = "HLA-DRA*01:01"
    for beta in supported_ii["dr_beta"]:
        flags = {"dr_alpha_assumed"}
        if beta in CURRENT_TO_OFFICIAL_COMPONENT:
            flags.add("official_name_maps_to_current_q_suffix")
        official_pairs[f"{dra}/{beta}"] = (dra, beta, tuple(sorted(flags)))
    for alpha, beta in product(supported_ii["dq_alpha"], supported_ii["dq_beta"]):
        flags = {"computational_pairing_not_haplotype_evidence"}
        if alpha in CURRENT_TO_OFFICIAL_COMPONENT or beta in CURRENT_TO_OFFICIAL_COMPONENT:
            flags.add("official_name_maps_to_current_q_suffix")
        official_pairs[f"{alpha}/{beta}"] = (
            alpha,
            beta,
            tuple(sorted(flags)),
        )
    for alpha, beta in product(supported_ii["dp_alpha"], supported_ii["dp_beta"]):
        official_pairs[f"{alpha}/{beta}"] = (
            alpha,
            beta,
            ("computational_pairing_not_haplotype_evidence",),
        )
    for mouse in supported_ii["mouse"]:
        if mouse in MOUSE_PAIRS:
            official_pairs[mouse] = (*MOUSE_PAIRS[mouse], ())

    pair_specs: Dict[str, Tuple[str, str, Tuple[str, ...]]] = dict(official_pairs)
    for name, rows in secondary_ii.items():
        if rows:
            pair_specs.setdefault(
                name,
                (rows[0]["alpha_allele"], rows[0]["beta_allele"], ("legacy_pair",)),
            )
    for name, (_sequence, alpha, beta, _raw) in official_ii.items():
        pair_specs.setdefault(name, (alpha, beta, ("training_reference_pair",)))

    final_ii: List[Dict[str, object]] = []
    missing_official_pairs = []
    official_pair_conflicts = []
    for name, (alpha, beta, pair_flags) in sorted(pair_specs.items()):
        official_supported = name in official_pairs
        aliases = {f"{alpha}+{beta}"}
        official_alpha = official_component_name(alpha)
        official_beta = official_component_name(beta)
        if (official_alpha, official_beta) != (alpha, beta):
            aliases.update(
                {
                    f"{official_alpha}/{official_beta}",
                    f"{official_alpha}+{official_beta}",
                }
            )
        if beta.startswith("HLA-DRB"):
            aliases.update({beta, beta.removeprefix("HLA-")})
            aliases.update(
                {official_beta, official_beta.removeprefix("HLA-")}
            )
        training_reference = official_ii.get(name)
        if training_reference is not None:
            aliases.add(training_reference[3])

        if alpha in component_values and beta in component_values:
            alpha_value = component_values[alpha]
            beta_value = component_values[beta]
            sequence = alpha_value.sequence + beta_value.sequence
            component_evidence = {alpha_value.evidence, beta_value.evidence}
            if "legacy_secondary" in component_evidence:
                evidence = "legacy_secondary"
            elif "derived_contact_complete_partial" in component_evidence:
                evidence = "derived_contact_complete_partial"
            elif "derived_contact_complete_full" in component_evidence:
                evidence = "derived_contact_complete_full"
            elif "derived_exact_pbr" in component_evidence:
                evidence = "derived_exact_pbr"
            else:
                evidence = "legacy_secondary"
            validation = (
                "official_supported_component_composition"
                if official_supported
                else "legacy_component_composition"
            )
            source = "component-wise 15-aa alpha plus 19-aa beta pseudo composition"
            source_version = f"{OFFICIAL_II_VERSION}-support-space"
            source_ids = tuple(
                value
                for value in (component_ids.get(alpha, ""), component_ids.get(beta, ""))
                if value
            )

            if training_reference is not None:
                official_sequence = training_reference[0]
                if sequence == official_sequence:
                    validation = "official_training_exact"
                else:
                    official_pair_conflicts.append(name)
                    matching_secondary = [
                        row
                        for row in secondary_ii.get(name, ())
                        if row["pseudo_sequence"] == official_sequence
                    ]
                    if not matching_secondary:
                        continue
                    # This exact value was already present in the pinned
                    # secondary snapshot; the validation oracle is not copied.
                    old = matching_secondary[0]
                    sequence = old["pseudo_sequence"]
                    evidence = "legacy_secondary"
                    validation = "official_training_exact"
                    source = old["source"]
                    source_version = old["source_version"]
                    source_ids = ()
        elif secondary_ii.get(name):
            sequences = {row["pseudo_sequence"] for row in secondary_ii[name]}
            if len(sequences) != 1:
                continue
            sequence = next(iter(sequences))
            old = secondary_ii[name][0]
            evidence = "legacy_secondary"
            validation = (
                "official_training_exact"
                if training_reference is not None
                and sequence == training_reference[0]
                else "legacy_secondary_only"
            )
            source = old["source"]
            source_version = old["source_version"]
            source_ids = ()
        else:
            if official_supported:
                missing_official_pairs.append(name)
            continue
        for old in secondary_ii.get(name, ()):
            aliases.update(split_values(old["aliases"]))
            if old["pseudo_sequence"] != sequence:
                loser = normalized_base_row(old)
                loser["reason"] = (
                    "secondary pseudo conflicts with a higher-priority "
                    "independently reconstructed component reference"
                )
                conflict_rows.append(loser)
        aliases.discard(name)
        species = "human" if alpha.startswith("HLA-") else "mouse"
        flags = set(pair_flags)
        if alpha in component_values:
            flags.update(component_values[alpha].quality_flags)
        if beta in component_values:
            flags.update(component_values[beta].quality_flags)
        if "X" in sequence:
            flags.add(
                "contains_official_alignment_gap_token"
                if validation == "official_training_exact"
                else "contains_alignment_gap_or_unknown_token"
            )
        if not official_supported:
            flags.add("outside_current_netmhciipan_support_list")
        final_ii.append(
            {
                "record_id": record_id("pseudo", PSEUDO_II_DEFINITION, name, sequence),
                "canonical_name": name,
                "aliases": ";".join(sorted(aliases)),
                "species": species,
                "mhc_class": "II",
                "definition_id": PSEUDO_II_DEFINITION,
                "pseudo_sequence": sequence,
                "alpha_allele": alpha,
                "beta_allele": beta,
                "source": source,
                "source_version": source_version,
                "evidence_level": evidence,
                "official_support": str(official_supported).lower(),
                "validation_status": validation,
                "source_record_ids": ";".join(source_ids),
                "quality_flags": ";".join(sorted(flags)),
            }
        )

    final_rows = final_i + final_ii
    official_i_available = {row["canonical_name"] for row in final_i if row["official_support"] == "true"}
    official_ii_available = {row["canonical_name"] for row in final_ii if row["official_support"] == "true"}
    missing_supported_i = sorted(supported_i - official_i_available)
    missing_supported_ii = sorted(set(official_pairs) - official_ii_available)

    existing_quarantine_path = data_dir / "quarantined_pseudosequences.csv"
    existing_quarantine = (
        read_csv(existing_quarantine_path) if existing_quarantine_path.is_file() else []
    )
    normalized_quarantine = []
    for row in existing_quarantine:
        if row.get("reason", "").startswith("secondary pseudo conflicts with"):
            # These rows are generated by this overlay builder.  Recompute
            # them from the immutable secondary snapshot on every run so stale
            # conflicts cannot accumulate after a derivation rule is fixed.
            continue
        normalized = normalized_base_row(row)
        normalized["reason"] = row.get("reason", "legacy source conflict")
        normalized_quarantine.append(normalized)
    quarantine_unique = {
        (str(row.get("record_id", "")), str(row.get("reason", ""))): row
        for row in (*normalized_quarantine, *conflict_rows)
    }

    write_csv(
        data_dir / "pseudo_sequences.csv",
        PSEUDO_FIELDS,
        final_rows,
        key_fields=("mhc_class", "species", "canonical_name", "pseudo_sequence"),
    )
    write_csv(
        data_dir / "pseudo_components.csv",
        COMPONENT_FIELDS,
        component_rows,
        key_fields=("species", "chain_role", "locus", "allele"),
    )
    write_csv(
        existing_quarantine_path,
        CONFLICT_FIELDS,
        quarantine_unique.values(),
        key_fields=("mhc_class", "canonical_name", "reason", "pseudo_sequence"),
    )

    by_evidence = Counter(str(row["evidence_level"]) for row in final_rows)
    report = {
        "official_versions": {
            "netmhcpan": OFFICIAL_I_VERSION,
            "netmhciipan": OFFICIAL_II_VERSION,
        },
        "official_training_validation": {
            "mhc_i_scope": len(official_i),
            "mhc_i_final_exact_matches": sum(
                1
                for name, (sequence, _raw) in official_i.items()
                if any(
                    row["canonical_name"] == name and row["pseudo_sequence"] == sequence
                    for row in final_i
                )
            ),
            "mhc_i_derived_conflicts": validation_conflicts_i,
            "mhc_ii_scope": len(official_ii),
            "mhc_ii_final_exact_matches": sum(
                1
                for name, (sequence, _alpha, _beta, _raw) in official_ii.items()
                if any(
                    row["canonical_name"] == name and row["pseudo_sequence"] == sequence
                    for row in final_ii
                )
            ),
            "mhc_ii_component_derived_conflicts": component_validation_conflicts,
            "mhc_ii_composition_conflicts": official_pair_conflicts,
            "training_pseudo_files_are_not_bundled": True,
            "training_pseudo_rows_are_not_copied": True,
        },
        "support_coverage": {
            "mhc_i_supported_names": len(supported_i),
            "mhc_i_available": len(official_i_available),
            "mhc_i_missing": missing_supported_i,
            "mhc_ii_supported_molecules": len(official_pairs),
            "mhc_ii_available": len(official_ii_available),
            "mhc_ii_missing": missing_supported_ii,
        },
        "counts": {
            "pseudo_total": len(final_rows),
            "pseudo_i": len(final_i),
            "pseudo_ii": len(final_ii),
            "pseudo_components": len(component_rows),
            "quarantined": len(quarantine_unique),
            "by_evidence": dict(sorted(by_evidence.items())),
            "secondary_component_conflicts": len(component_secondary_conflicts),
            "deposited_derivation_conflicts": len(deposited_conflicts),
            "official_guided_mhci_indel_recoveries": len(official_guided_i),
        },
        "source_hashes": {
            "netmhcpan_support_list": sha256_path(args.supported_i),
            "netmhciipan_support_list": sha256_path(args.supported_ii),
            "netmhcpan_training_pseudo_validation": official_i_hash,
            "netmhciipan_training_pseudo_validation": official_ii_hash,
        },
    }
    validation_path = data_dir / "pseudo_validation.json"
    with validation_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    metadata_path = data_dir / "metadata.json"
    with metadata_path.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    metadata["schema_version"] = max(int(metadata.get("schema_version", 1)), 4)
    metadata["generated_on"] = date.today().isoformat()
    pseudo_by_species_class = Counter(
        f"{row['species']}/{row['mhc_class']}" for row in final_rows
    )
    metadata.setdefault("counts", {}).update(
        {
            "pseudo_total": len(final_rows),
            "pseudo_components": len(component_rows),
            "pseudo_quarantined": len(quarantine_unique),
            "pseudo_by_evidence": dict(sorted(by_evidence.items())),
            "pseudo_by_species_class": dict(
                sorted(pseudo_by_species_class.items())
            ),
            "netmhcpan_supported_available": len(official_i_available),
            "netmhcpan_supported_total": len(supported_i),
            "netmhciipan_supported_available": len(official_ii_available),
            "netmhciipan_supported_total": len(official_pairs),
        }
    )
    metadata["quarantine"] = {
        "count": len(quarantine_unique),
        "file": "quarantined_pseudosequences.csv",
        "rule": (
            "Preserve fixed legacy-source conflicts and recompute current "
            "secondary-vs-derived conflicts on every overlay build; never "
            "silently select two sequences for one canonical name."
        ),
    }
    metadata.setdefault("sources", {}).update(
        {
            "netmhcpan_4_2_supported_names": {
                "name": f"{OFFICIAL_I_VERSION} public allele.list",
                "url": OFFICIAL_I_LIST_URL,
                "sha256": OFFICIAL_I_LIST_SHA256,
                "retrieved_on": "2026-11-20",
            },
            "netmhciipan_4_3_supported_names": {
                "name": f"{OFFICIAL_II_VERSION} public alleles_name.txt",
                "url": OFFICIAL_II_LIST_URL,
                "sha256": OFFICIAL_II_LIST_SHA256,
                "retrieved_on": "2026-11-20",
            },
            "netmhc_training_pseudo_validation": {
                "netmhcpan_sha256": official_i_hash,
                "netmhciipan_sha256": official_ii_hash,
                "upstream_last_modified": "2026-07-06",
                "bundled": False,
                "rows_copied_into_registry": False,
                "purpose": "local validation gate only",
            },
        }
    )
    metadata["pseudo_derivation"] = {
        "priority": [
            "official training validation of independently reconstructed rows",
            "derived exact PBR",
            "contact-complete deposited sequence",
            "pinned secondary mirror compatibility row",
        ],
        "mhc_i_contact_positions_1based": list(MHC_I_POSITIONS),
        "mhc_ii_component_contact_positions_1based": {
            locus: list(positions)
            for locus, (_anchor, positions) in MHC_II_CONTACTS.items()
        },
        "mhc_ii_contact_anchors": {
            locus: anchor for locus, (anchor, _positions) in MHC_II_CONTACTS.items()
        },
        "dynamic_pairing_note": (
            "NetMHCIIpan-supported DQ/DP alpha-beta products describe computational "
            "support, not population haplotype evidence."
        ),
        "validation_file": "pseudo_validation.json",
    }
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return report


def parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=PROJECT_ROOT / "src" / "mhctransform" / "data",
    )
    parser.add_argument(
        "--supported-i",
        type=Path,
        default=PROJECT_ROOT / "tools" / "reference" / "NetMHCpan-4.2-allele.list",
    )
    parser.add_argument(
        "--supported-ii",
        type=Path,
        default=PROJECT_ROOT / "tools" / "reference" / "NetMHCIIpan-4.3-alleles_name.txt",
    )
    parser.add_argument("--official-pseudo-i", type=Path, required=True)
    parser.add_argument("--official-pseudo-ii", type=Path, required=True)
    return parser


def main() -> int:
    args = parser().parse_args()
    report = build(args)
    print(json.dumps(report["counts"], indent=2, sort_keys=True))
    compact_coverage = {
        key: len(value) if isinstance(value, list) else value
        for key, value in report["support_coverage"].items()
    }
    print(json.dumps(compact_coverage, indent=2, sort_keys=True))
    print("Full missing-name lists: pseudo_validation.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
