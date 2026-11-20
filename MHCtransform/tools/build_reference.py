#!/usr/bin/env python3
"""Build the bundled MHC registry from pinned/local primary data.

Inputs
------
* IPD-IMGT/HLA ``hla_prot.fasta`` for human precursor proteins.
* Pinned IPD-IMGT/HLA protein alignments plus UniProt topology annotations
  for mature, soluble extracellular chains suitable for AF3 input.
* UniProt entries for a curated set of common mouse molecules/B2M.  Reviewed
  records are preferred; explicitly labelled fragment/unreviewed entries are
  retained only with ``is_complete=false`` unless their sequence is complete.
* PUFFIN's public NetMHCpan/NetMHCIIpan-format pseudosequence mirrors.

The generated CSV files are runtime assets.  Runtime lookup never needs a
network connection and never derives a pseudosequence by slicing an unaligned
precursor protein.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Sequence, Tuple


PSEUDO_I_DEFINITION = "netmhcpan-34-v1"
PSEUDO_II_DEFINITION = "netmhciipan-34-v1"
AF3_DEFINITION = "af3-soluble-extracellular-v1"
PUFFIN_COMMIT = "eddf17d972285a89e80e2837adc1ee6308f5b96b"
PSEUDO_I_URL = (
    f"https://raw.githubusercontent.com/gifford-lab/PUFFIN/{PUFFIN_COMMIT}/"
    "data/MHC_pseudo.dat"
)
PSEUDO_II_URL = (
    f"https://raw.githubusercontent.com/gifford-lab/PUFFIN/{PUFFIN_COMMIT}/data/"
    "pseudosequence.2016.all.X.dat"
)
HLA_RELEASE = "3.66.0"
HLA_COMMIT = "b2055ea"
HLA_FASTA_URL = (
    f"https://raw.githubusercontent.com/ANHIG/IMGTHLA/{HLA_COMMIT}/hla_prot.fasta"
)
HLA_STATUS_URL = (
    f"https://raw.githubusercontent.com/ANHIG/IMGTHLA/{HLA_COMMIT}/Allele_status.txt"
)
HLA_2009_URL = (
    f"https://raw.githubusercontent.com/ANHIG/IMGTHLA/{HLA_COMMIT}/Nomenclature_2009.txt"
)
HLA_DELETED_URL = (
    f"https://raw.githubusercontent.com/ANHIG/IMGTHLA/{HLA_COMMIT}/Deleted_alleles.txt"
)

# UniProt REST accession batches are mutable even when the accession list is
# fixed.  Pin their content hashes so a later upstream update cannot silently
# alter full proteins or AF3 boundaries during a rebuild.
UNIPROT_FASTA_SHA256 = "ab801102d204e11362f236bac1553c3c850dad9abf26e84ad8898ee486d77fa1"
UNIPROT_TOPOLOGY_SHA256 = "e04ec670bb801ce850d941bc760c0092944d89cc3107accb89fd106e8989e3e0"

HLA_ALIGNMENT_FILES: Mapping[str, str] = {
    "A": "A_prot.txt",
    "B": "B_prot.txt",
    "C": "C_prot.txt",
    "E": "E_prot.txt",
    "F": "F_prot.txt",
    "G": "G_prot.txt",
    "DRA": "DRA_prot.txt",
    "DRB1": "DRB_prot.txt",
    "DRB3": "DRB_prot.txt",
    "DRB4": "DRB_prot.txt",
    "DRB5": "DRB_prot.txt",
    "DQA1": "DQA1_prot.txt",
    "DQB1": "DQB1_prot.txt",
    "DPA1": "DPA1_prot.txt",
    "DPB1": "DPB1_prot.txt",
}

# Reviewed UniProt records provide the signal-peptide and transmembrane
# topology used to define the reference extracellular span for each IPD
# alignment.  The actual allele sequence is always taken from the pinned IPD
# release, not from these reference records.
HUMAN_TOPOLOGY_ACCESSIONS: Mapping[str, str] = {
    "A": "P04439",
    "B": "P01889",
    "C": "P10321",
    "E": "P13747",
    "F": "P30511",
    "G": "P17693",
    "DRA": "P01903",
    "DRB1": "P01911",
    "DRB3": "P79483",
    "DRB4": "P13762",
    "DRB5": "Q30154",
    "DQA1": "P01909",
    "DQB1": "P01920",
    "DPA1": "P20036",
    "DPB1": "P04440",
}

AA = set("ACDEFGHIKLMNPQRSTVWY")
AA_X = AA | {"X"}

HUMAN_LOCI: Mapping[str, Tuple[str, str, int]] = {
    "A": ("I", "heavy", 350),
    "B": ("I", "heavy", 350),
    "C": ("I", "heavy", 350),
    "E": ("I", "heavy", 340),
    "F": ("I", "heavy", 330),
    "G": ("I", "heavy", 330),
    "DRA": ("II", "alpha", 245),
    "DRB1": ("II", "beta", 250),
    "DRB3": ("II", "beta", 250),
    "DRB4": ("II", "beta", 250),
    "DRB5": ("II", "beta", 250),
    "DQA1": ("II", "alpha", 245),
    "DQB1": ("II", "beta", 250),
    "DPA1": ("II", "alpha", 250),
    "DPB1": ("II", "beta", 250),
}


# UniProt accessions and explicit names are intentionally curated rather than
# inferred from a free-text FASTA header.  Partial entries are retained
# with is_complete=false, so a caller must opt in to use them.
MOUSE_PROTEINS: Tuple[Mapping[str, object], ...] = (
    # MHC-I heavy chains.
    {"accession": "P01899", "allele": "H-2-Db", "aliases": "H2-D1*b", "class": "I", "role": "heavy", "locus": "H2-D1"},
    {"accession": "P01900", "allele": "H-2-Dd", "aliases": "H2-D1*d", "class": "I", "role": "heavy", "locus": "H2-D1"},
    {"accession": "P14426", "allele": "H-2-Dk", "aliases": "H2-D1*k", "class": "I", "role": "heavy", "locus": "H2-D1"},
    {"accession": "P14427", "allele": "H-2-Dp", "aliases": "H2-D1*p", "class": "I", "role": "heavy", "locus": "H2-D1"},
    {"accession": "Q31145", "allele": "H-2-Dq", "aliases": "H2-D1*q", "class": "I", "role": "heavy", "locus": "H2-D1"},
    {"accession": "P01901", "allele": "H-2-Kb", "aliases": "H2-K1*b", "class": "I", "role": "heavy", "locus": "H2-K1"},
    {"accession": "P01902", "allele": "H-2-Kd", "aliases": "H2-K1*d", "class": "I", "role": "heavy", "locus": "H2-K1"},
    {"accession": "P04223", "allele": "H-2-Kk", "aliases": "H2-K1*k", "class": "I", "role": "heavy", "locus": "H2-K1"},
    {"accession": "P14428", "allele": "H-2-Kq", "aliases": "H2-K1*q", "class": "I", "role": "heavy", "locus": "H2-K1"},
    {"accession": "P01897", "allele": "H-2-Ld", "aliases": "H2-L*d", "class": "I", "role": "heavy", "locus": "H2-L"},
    {"accession": "Q31151", "allele": "H-2-Lq", "aliases": "H2-L*q", "class": "I", "role": "heavy", "locus": "H2-L"},
    # MHC-II I-A chains.
    {"accession": "P14434", "allele": "H2-Aa*b", "aliases": "H-2-IAb-alpha", "class": "II", "role": "alpha", "locus": "H2-Aa"},
    {"accession": "P14483", "allele": "H2-Ab1*b", "aliases": "H-2-IAb-beta", "class": "II", "role": "beta", "locus": "H2-Ab1"},
    {"accession": "P04228", "allele": "H2-Aa*d", "aliases": "H-2-IAd-alpha", "class": "II", "role": "alpha", "locus": "H2-Aa"},
    {"accession": "P01921", "allele": "H2-Ab1*d", "aliases": "H-2-IAd-beta", "class": "II", "role": "beta", "locus": "H2-Ab1"},
    {"accession": "P01910", "allele": "H2-Aa*k", "aliases": "H-2-IAk-alpha", "class": "II", "role": "alpha", "locus": "H2-Aa"},
    {"accession": "P06343", "allele": "H2-Ab1*k", "aliases": "H-2-IAk-beta", "class": "II", "role": "beta", "locus": "H2-Ab1"},
    {"accession": "P04227", "allele": "H2-Aa*q", "aliases": "H-2-IAq-alpha", "class": "II", "role": "alpha", "locus": "H2-Aa"},
    {"accession": "P06342", "allele": "H2-Ab1*q", "aliases": "H-2-IAq-beta", "class": "II", "role": "beta", "locus": "H2-Ab1"},
    {"accession": "P14437", "allele": "H2-Aa*s", "aliases": "H-2-IAs-alpha", "class": "II", "role": "alpha", "locus": "H2-Aa"},
    {"accession": "P06345", "allele": "H2-Ab1*s", "aliases": "H-2-IAs-beta", "class": "II", "role": "beta", "locus": "H2-Ab1"},
    {"accession": "P14438", "allele": "H2-Aa*u", "aliases": "H-2-IAu-alpha", "class": "II", "role": "alpha", "locus": "H2-Aa"},
    {"accession": "P06344", "allele": "H2-Ab1*u", "aliases": "H-2-IAu-beta", "class": "II", "role": "beta", "locus": "H2-Ab1"},
    # MHC-II I-E chains.  The accessions are kept explicit so provenance is auditable.
    {"accession": "P01904", "allele": "H2-Ea*d", "aliases": "H-2-IEd-alpha", "class": "II", "role": "alpha", "locus": "H2-Ea"},
    {"accession": "P01915", "allele": "H2-Eb1*d", "aliases": "H-2-IEd-beta", "class": "II", "role": "beta", "locus": "H2-Eb1"},
    {"accession": "P04224", "allele": "H2-Ea*k", "aliases": "H-2-IEk-alpha", "class": "II", "role": "alpha", "locus": "H2-Ea"},
    {"accession": "P04230", "allele": "H2-Eb1*b", "aliases": "H-2-IEb-beta", "class": "II", "role": "beta", "locus": "H2-Eb1"},
    {"accession": "Q31163", "allele": "H2-Eb1*k", "aliases": "H-2-IEk-beta", "class": "II", "role": "beta", "locus": "H2-Eb1"},
    # Species-matched class-I companion.
    {"accession": "P01887", "allele": "B2M-mouse", "aliases": "mouse-B2M;H2-B2m", "class": "I", "role": "companion", "locus": "B2m"},
    # Human companion is fetched from the same reviewed source.
    {"accession": "P61769", "allele": "B2M-human", "aliases": "human-B2M;HLA-B2M", "class": "I", "role": "companion", "locus": "B2M", "species": "human"},
)


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

# UniProt P14437 is a complete mature membrane protein without a deposited
# signal peptide.  It is not a precursor_full record, but its annotated
# extracellular region is exact and is therefore valid for the strict AF3
# representation.
MOUSE_EXACT_ECTODOMAIN_SPANS: Mapping[str, Tuple[int, int]] = {
    "P14437": (1, 195),
}


FULL_FIELDS = (
    "record_id",
    "allele",
    "aliases",
    "species",
    "mhc_class",
    "chain_role",
    "locus",
    "sequence_kind",
    "sequence",
    "source",
    "source_version",
    "source_accession",
    "is_complete",
)
AF3_FIELDS = (
    "record_id",
    "allele",
    "aliases",
    "species",
    "mhc_class",
    "chain_role",
    "locus",
    "definition_id",
    "sequence",
    "source_full_record_id",
    "precursor_span_1based",
    "boundary_method",
    "boundary_source",
    "source_accession",
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
)
PAIR_FIELDS = (
    "canonical_name",
    "aliases",
    "species",
    "mhc_class",
    "alpha_or_heavy",
    "beta_or_companion",
    "source",
)
CONFLICT_FIELDS = PSEUDO_FIELDS + ("reason",)
AF3_EXCLUSION_FIELDS = (
    "full_record_id",
    "allele",
    "aliases",
    "species",
    "mhc_class",
    "chain_role",
    "locus",
    "source_accession",
    "reason",
)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_sha256(path: Path, expected: str, label: str) -> Path:
    observed = sha256_path(path)
    if observed != expected:
        raise ValueError(
            f"{label} snapshot hash changed: {observed} != {expected}; "
            "review the upstream changes before updating the pinned hash"
        )
    return path


def record_id(*parts: str) -> str:
    payload = "\x1f".join(parts).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:20]


def parse_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header = None
    sequence: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence).upper()
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence).upper()


def two_field_name(full_name: str) -> str:
    name = full_name.upper().replace("HLA-", "")
    locus, fields_text = name.split("*", 1)
    fields = fields_text.split(":")
    second = re.match(r"[0-9]+", fields[1])
    if second is None:
        raise ValueError(f"cannot parse HLA allele: {full_name}")
    suffix = re.search(r"([A-Z])$", fields_text)
    return (
        f"HLA-{locus}*{fields[0]}:{second.group(0)}"
        f"{suffix.group(1) if suffix else ''}"
    )


def human_aliases(full_name: str, canonical: str) -> Tuple[str, ...]:
    unprefixed = full_name.removeprefix("HLA-")
    canonical_unprefixed = canonical.removeprefix("HLA-")
    # Do not create separator-free aliases.  Once an HLA field contains three
    # or more digits, compact spellings are intrinsically ambiguous (for
    # example A*02:11:02 versus A*02:1102).
    return tuple(sorted({full_name, unprefixed, canonical_unprefixed} - {canonical}))


def parse_hla_status(path: Path) -> Mapping[str, Mapping[str, str]]:
    lines = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    if not lines:
        raise ValueError(f"empty HLA status file: {path}")
    return {row["Allele"].upper(): row for row in csv.DictReader(lines)}


def parse_nomenclature_2009(path: Path) -> Mapping[str, str]:
    """Map pre-2010 compact names to current colon-delimited names."""

    mapping = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) < 2 or "*" not in fields[0]:
            continue
        mapping[fields[0].upper()] = fields[1].upper()
    if "A*9201" not in mapping or mapping["A*9201"] != "A*02:101":
        raise ValueError("unexpected Nomenclature_2009.txt format")
    return mapping


def parse_deleted_renames(path: Path) -> Mapping[str, str]:
    lines = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    mapping = {}
    for row in csv.DictReader(lines):
        match = re.search(
            r"(?:\brenamed|\bidentical\s+to(?:\s+corrected\s+sequence\s+of)?|"
            r"\bsuffix\s+added\s+to\s+allele\s+becoming)\s+"
            r"([A-Z0-9]+\*[0-9:]+[A-Z]?)\b",
            row["Description"],
            flags=re.IGNORECASE,
        )
        if match:
            mapping[row["Allele"].upper()] = match.group(1).upper()
    return mapping


def hla_two_field(name: str) -> str:
    normalized = name.upper().removeprefix("HLA-")
    locus, fields_text = normalized.split("*", 1)
    if locus == "CW":
        locus = "C"
    fields = fields_text.split(":")
    if len(fields) < 2:
        raise ValueError(f"HLA name has fewer than two fields: {name}")
    second = re.match(r"[0-9]+", fields[1])
    if second is None:
        raise ValueError(f"invalid HLA second field: {name}")
    suffix = re.search(r"([A-Z])$", fields_text)
    return (
        f"HLA-{locus}*{fields[0]}:{second.group(0)}"
        f"{suffix.group(1) if suffix else ''}"
    )


def official_two_field_names(status_path: Path) -> set[str]:
    names = set()
    for name in parse_hla_status(status_path):
        if "*" not in name or name.split("*", 1)[0] not in HUMAN_LOCI:
            continue
        if ":" not in name.split("*", 1)[1]:
            continue
        names.add(hla_two_field(name))
    return names


def build_human_full(path: Path, status_path: Path) -> List[Dict[str, str]]:
    source_hash = sha256_path(path)
    status_hash = sha256_path(status_path)
    status_by_allele = parse_hla_status(status_path)
    source_version = (
        f"IPD-IMGT/HLA {HLA_RELEASE} commit:{HLA_COMMIT} "
        f"protein-sha256:{source_hash[:12]} status-sha256:{status_hash[:12]}"
    )
    grouped: MutableMapping[Tuple[str, str], Dict[str, object]] = {}
    header_re = re.compile(r"^(HLA:[^\s]+)\s+([^\s]+)")
    for header, sequence in parse_fasta(path):
        match = header_re.match(header)
        if not match or "*" not in match.group(2):
            continue
        accession, deposited = match.groups()
        locus = deposited.split("*", 1)[0].upper()
        if locus not in HUMAN_LOCI:
            continue
        mhc_class, role, _ = HUMAN_LOCI[locus]
        full_name = f"HLA-{deposited.upper()}"
        # Expression suffixes apply to the exact official allele name and can
        # originate in non-coding fields (for example A*02:01:01:02L).  Do not
        # invent a synthetic two-field name such as A*02:01L or merge it with
        # normally expressed alleles.  Unsuffixed records are still grouped by
        # their two-field protein designation plus exact precursor sequence.
        canonical = (
            full_name
            if re.search(r"[NLSQCA]$", full_name)
            else two_field_name(full_name)
        )
        key = (canonical, sequence)
        item = grouped.setdefault(
            key,
            {
                "canonical": canonical,
                "sequence": sequence,
                "locus": locus,
                "class": mhc_class,
                "role": role,
                "aliases": set(),
                "accessions": set(),
                "complete": False,
            },
        )
        item["aliases"].update(human_aliases(full_name, canonical))  # type: ignore[union-attr]
        item["accessions"].add(accession)  # type: ignore[union-attr]
        status = status_by_allele.get(deposited.upper())
        officially_full = bool(status and status.get("Partial") == "Full")
        item["complete"] = bool(item["complete"] or (officially_full and set(sequence) <= AA))

    rows = []
    for item in grouped.values():
        aliases = sorted(item["aliases"])  # type: ignore[arg-type]
        sequence = str(item["sequence"])
        canonical = str(item["canonical"])
        rows.append(
            {
                "record_id": record_id("full", "human", canonical, sequence),
                "allele": canonical,
                "aliases": ";".join(aliases),
                "species": "human",
                "mhc_class": str(item["class"]),
                "chain_role": str(item["role"]),
                "locus": str(item["locus"]),
                "sequence_kind": (
                    "precursor_full" if bool(item["complete"]) else "partial_deposited_protein"
                ),
                "sequence": sequence,
                "source": "IPD-IMGT/HLA hla_prot.fasta",
                "source_version": source_version,
                "source_accession": ",".join(sorted(item["accessions"])),  # type: ignore[arg-type]
                "is_complete": str(bool(item["complete"])).lower(),
            }
        )
    return rows


def download(url: str, path: Path, *, offline: bool) -> Path:
    if path.exists():
        return path
    if offline:
        raise FileNotFoundError(f"offline mode: missing cached source {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "MHCtransform-reference-builder/0.1"})
    with urllib.request.urlopen(request, timeout=120) as response, path.open("wb") as handle:
        handle.write(response.read())
    return path


def fetch_uniprot(cache_dir: Path, *, offline: bool) -> Path:
    accessions = sorted({str(item["accession"]) for item in MOUSE_PROTEINS})
    query = urllib.parse.urlencode(
        {"accessions": ",".join(accessions), "format": "fasta"}, safe=","
    )
    url = f"https://rest.uniprot.org/uniprotkb/accessions?{query}"
    accession_hash = hashlib.sha256(",".join(accessions).encode("ascii")).hexdigest()[:12]
    return verify_sha256(
        download(
            url,
            cache_dir / f"uniprot_mhc_mouse_b2m_{accession_hash}.fasta",
            offline=offline,
        ),
        UNIPROT_FASTA_SHA256,
        "UniProt full-protein",
    )


def fetch_uniprot_topology(cache_dir: Path, *, offline: bool) -> Path:
    accessions = sorted(
        {str(item["accession"]) for item in MOUSE_PROTEINS}
        | set(HUMAN_TOPOLOGY_ACCESSIONS.values())
    )
    fields = ",".join(
        (
            "accession",
            "id",
            "organism_name",
            "gene_primary",
            "length",
            "sequence",
            "ft_signal",
            "ft_chain",
            "ft_transmem",
            "ft_topo_dom",
        )
    )
    query = urllib.parse.urlencode(
        {"accessions": ",".join(accessions), "format": "tsv", "fields": fields},
        safe=",",
    )
    url = f"https://rest.uniprot.org/uniprotkb/accessions?{query}"
    accession_hash = hashlib.sha256(",".join(accessions).encode("ascii")).hexdigest()[:12]
    return verify_sha256(
        download(
            url,
            cache_dir / f"uniprot_mhc_topology_{accession_hash}.tsv",
            offline=offline,
        ),
        UNIPROT_TOPOLOGY_SHA256,
        "UniProt topology",
    )


def fetch_hla_alignments(cache_dir: Path, *, offline: bool) -> Mapping[str, Path]:
    paths = {}
    for filename in sorted(set(HLA_ALIGNMENT_FILES.values())):
        url = (
            f"https://raw.githubusercontent.com/ANHIG/IMGTHLA/{HLA_COMMIT}/"
            f"alignments/{filename}"
        )
        paths[filename] = download(
            url,
            cache_dir / f"IMGTHLA-{HLA_RELEASE}-{filename}",
            offline=offline,
        )
    return paths


def parse_uniprot_topology(path: Path) -> Mapping[str, Mapping[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {row["Entry"]: row for row in rows}
    if not result:
        raise ValueError(f"empty UniProt topology table: {path}")
    return result


def _exact_feature_bounds(text: str, feature: str) -> Tuple[int, int] | None:
    match = re.search(rf"(?:^|; )?{feature} ([0-9]+)\.\.([0-9]+)(?:;|$)", text)
    if match is None:
        return None
    return int(match.group(1)), int(match.group(2))


def extracellular_span_from_topology(
    topology: Mapping[str, str], *, companion: bool = False
) -> Tuple[int, int] | None:
    topological_domain = topology.get("Topological domain", "")
    extracellular = re.search(
        r"TOPO_DOM ([0-9]+)\.\.([0-9]+); /note=\"Extracellular\"",
        topological_domain,
    )
    if extracellular:
        return int(extracellular.group(1)), int(extracellular.group(2))

    if companion:
        return _exact_feature_bounds(topology.get("Chain", ""), "CHAIN")

    signal = _exact_feature_bounds(topology.get("Signal peptide", ""), "SIGNAL")
    transmembrane = _exact_feature_bounds(topology.get("Transmembrane", ""), "TRANSMEM")
    if signal is None or transmembrane is None:
        return None
    return signal[1] + 1, transmembrane[0] - 1


def parse_hla_protein_alignment(path: Path) -> Tuple[str, Mapping[str, str]]:
    """Expand an IPD text alignment into full aligned allele sequences."""
    blocks: List[Dict[str, str]] = []
    block: Dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if raw_line.startswith(" Prot"):
            if block:
                blocks.append(block)
            block = {}
            continue
        if len(raw_line) < 19:
            continue
        name = raw_line[:19].strip().upper()
        if "*" not in name or " " in name:
            continue
        if name in block:
            raise ValueError(f"duplicate alignment row in one block: {name}")
        block[name] = raw_line[19:].upper()
    if block:
        blocks.append(block)
    if not blocks:
        raise ValueError(f"no protein alignment rows parsed from {path}")

    # The first allele is IPD's global reference.  Most displayed blocks list
    # it first, while a final insertion block can contain only alleles that
    # extend beyond the reference.  Display spaces separate residue groups;
    # only leading spaces in a block represent empty alignment columns.
    reference_name = next(iter(blocks[0]))
    block_segments: List[Tuple[str, Mapping[str, str]]] = []
    all_names: set[str] = set()
    for current in blocks:
        has_global_reference = reference_name in current
        if has_global_reference:
            local_reference_name = reference_name
        else:
            # A terminal insertion block may contain only the few alleles that
            # extend past the global reference, and their extensions need not
            # have equal lengths.  The longest displayed row defines the
            # block width; absent trailing positions are alignment gaps.
            local_reference_name = max(
                current,
                key=lambda name: len(re.sub(r"\s+", "", current[name])),
            )
        reference_display = current[local_reference_name].rstrip()
        leading_gaps = len(reference_display) - len(reference_display.lstrip(" "))
        local_reference = "." * leading_gaps + re.sub(
            r"\s+", "", reference_display[leading_gaps:]
        )
        if not local_reference:
            continue
        resolved_rows: Dict[str, str] = {}
        for name, display in current.items():
            prefix = display[:leading_gaps].ljust(leading_gaps)
            prefix_segment = "".join(
                "." if residue == " " else residue for residue in prefix
            )
            body = re.sub(r"\s+", "", display[leading_gaps:])
            segment = prefix_segment + body
            if len(segment) > len(local_reference):
                raise ValueError(
                    f"alignment width overflow in {path} for {name}: "
                    f"{len(segment)} > {len(local_reference)}"
                )
            if len(segment) < len(local_reference):
                if has_global_reference:
                    segment += "".join(
                        "." if residue == "." else "*"
                        for residue in local_reference[len(segment) :]
                    )
                else:
                    segment += "." * (len(local_reference) - len(segment))
            resolved_rows[name] = "".join(
                ref_residue if residue == "-" else residue
                for residue, ref_residue in zip(segment, local_reference)
            )
        canonical_segment = resolved_rows.get(reference_name, "." * len(local_reference))
        block_segments.append((canonical_segment, resolved_rows))
        all_names.update(resolved_rows)

    aligned = {}
    for name in all_names:
        parts = []
        for canonical_segment, rows in block_segments:
            parts.append(
                rows.get(
                    name,
                    "".join(
                        "." if residue == "." else "*"
                        for residue in canonical_segment
                    ),
                )
            )
        aligned[name] = "".join(parts)
    return reference_name, {
        name: sequence for name, sequence in aligned.items()
    }


def _aligned_plain_sequence(aligned: str) -> str:
    return "".join(
        "X" if residue == "*" else residue
        for residue in aligned
        if residue != "."
    )


def _reference_span_columns(
    aligned_reference: str, start_1based: int, end_1based: int
) -> Tuple[int, int]:
    residue_number = 0
    start_column = None
    end_column = None
    for column, residue in enumerate(aligned_reference):
        if residue not in AA_X:
            continue
        residue_number += 1
        if residue_number == start_1based:
            start_column = column
        if residue_number == end_1based:
            end_column = column
            break
    if start_column is None or end_column is None:
        raise ValueError(
            f"reference alignment does not cover topology span {start_1based}-{end_1based}"
        )
    return start_column, end_column


def _af3_segment_is_plausible(
    full_sequence: str,
    segment: str,
    start_1based: int,
    end_1based: int,
    *,
    mhc_class: str,
    chain_role: str,
) -> bool:
    if not segment or set(segment) - AA:
        return False
    minimum = 90 if chain_role == "companion" else (270 if mhc_class == "I" else 175)
    if len(segment) < minimum:
        return False
    if full_sequence[start_1based - 1 : end_1based] != segment:
        return False
    if chain_role == "companion":
        return start_1based > 1 and end_1based == len(full_sequence)
    transmembrane_probe = full_sequence[end_1based : end_1based + 20]
    hydrophobic = set("ACFGILMVWY")
    return (
        15 <= start_1based - 1 <= 35
        and len(transmembrane_probe) >= 15
        and sum(residue in hydrophobic for residue in transmembrane_probe)
        / len(transmembrane_probe)
        >= 0.55
    )


def uniprot_by_accession(path: Path) -> Mapping[str, Tuple[str, str]]:
    result = {}
    for header, sequence in parse_fasta(path):
        parts = header.split("|", 2)
        accession = parts[1] if len(parts) >= 3 else header.split()[0]
        result[accession] = (header, sequence)
    return result


def mouse_complete(item: Mapping[str, object], header: str, sequence: str) -> bool:
    if set(sequence) - AA:
        return False
    if "(Fragment)" in header:
        return False
    if item["role"] == "companion":
        return len(sequence) >= 110
    if item["class"] == "I":
        return len(sequence) >= 320
    return len(sequence) >= 240


def build_mouse_and_b2m(path: Path) -> List[Dict[str, str]]:
    source_hash = sha256_path(path)
    entries = uniprot_by_accession(path)
    rows = []
    missing = []
    for item in MOUSE_PROTEINS:
        accession = str(item["accession"])
        if accession not in entries:
            missing.append(accession)
            continue
        header, sequence = entries[accession]
        species = str(item.get("species", "mouse"))
        allele = str(item["allele"])
        rows.append(
            {
                "record_id": record_id("full", species, allele, sequence),
                "allele": allele,
                "aliases": str(item["aliases"]),
                "species": species,
                "mhc_class": str(item["class"]),
                "chain_role": str(item["role"]),
                "locus": str(item["locus"]),
                "sequence_kind": (
                    "precursor_full"
                    if mouse_complete(item, header, sequence)
                    else "partial_deposited_protein"
                ),
                "sequence": sequence,
                "source": "UniProtKB",
                "source_version": f"download-sha256:{source_hash[:16]}",
                "source_accession": accession,
                "is_complete": str(mouse_complete(item, header, sequence)).lower(),
            }
        )
    if missing:
        print("warning: UniProt returned no record for " + ", ".join(missing), file=sys.stderr)
    return rows


def _af3_exclusion(full_row: Mapping[str, str], reason: str) -> Dict[str, str]:
    return {
        "full_record_id": full_row["record_id"],
        "allele": full_row["allele"],
        "aliases": full_row["aliases"],
        "species": full_row["species"],
        "mhc_class": full_row["mhc_class"],
        "chain_role": full_row["chain_role"],
        "locus": full_row["locus"],
        "source_accession": full_row["source_accession"],
        "reason": reason,
    }


def _add_af3_candidate(
    grouped: MutableMapping[Tuple[str, ...], Dict[str, object]],
    full_row: Mapping[str, str],
    sequence: str,
    start_1based: int,
    end_1based: int,
    *,
    boundary_method: str,
    boundary_source: str,
) -> None:
    key = (
        full_row["species"],
        full_row["mhc_class"],
        full_row["chain_role"],
        full_row["locus"],
        full_row["allele"],
        sequence,
        full_row["record_id"],
        str(start_1based),
        str(end_1based),
    )
    item = grouped.setdefault(
        key,
        {
            "aliases": set(),
            "methods": set(),
            "boundary_sources": set(),
            "source_accession": full_row["source_accession"],
        },
    )
    item["aliases"].update(part for part in full_row["aliases"].split(";") if part)  # type: ignore[union-attr]
    item["methods"].add(boundary_method)  # type: ignore[union-attr]
    item["boundary_sources"].add(boundary_source)  # type: ignore[union-attr]


def _finalize_af3_candidates(
    grouped: Mapping[Tuple[str, ...], Mapping[str, object]]
) -> List[Dict[str, str]]:
    rows = []
    for key, item in grouped.items():
        (
            species,
            mhc_class,
            chain_role,
            locus,
            allele,
            sequence,
            source_full_record_id,
            start_1based,
            end_1based,
        ) = key
        aliases = set(item["aliases"])  # type: ignore[arg-type]
        aliases.discard(allele)
        rows.append(
            {
                "record_id": record_id(
                    "af3",
                    AF3_DEFINITION,
                    species,
                    allele,
                    sequence,
                    source_full_record_id,
                ),
                "allele": allele,
                "aliases": ";".join(sorted(aliases)),
                "species": species,
                "mhc_class": mhc_class,
                "chain_role": chain_role,
                "locus": locus,
                "definition_id": AF3_DEFINITION,
                "sequence": sequence,
                "source_full_record_id": source_full_record_id,
                "precursor_span_1based": f"{start_1based}-{end_1based}",
                "boundary_method": "; ".join(sorted(item["methods"])),  # type: ignore[arg-type]
                "boundary_source": "; ".join(sorted(item["boundary_sources"])),  # type: ignore[arg-type]
                "source_accession": str(item["source_accession"]),
            }
        )
    return rows


def build_af3_sequences(
    full_rows: Sequence[Mapping[str, str]],
    alignment_paths: Mapping[str, Path],
    topology_path: Path,
) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    """Build mature soluble extracellular chains with auditable boundaries."""

    topology = parse_uniprot_topology(topology_path)
    parsed_alignments: Dict[str, Tuple[str, Mapping[str, str]]] = {}
    locus_alignment: Dict[
        str, Tuple[Mapping[str, str], int, int, int, int, str, str]
    ] = {}

    for locus, filename in HLA_ALIGNMENT_FILES.items():
        if filename not in parsed_alignments:
            parsed_alignments[filename] = parse_hla_protein_alignment(
                alignment_paths[filename]
            )
        _, aligned = parsed_alignments[filename]
        topology_accession = HUMAN_TOPOLOGY_ACCESSIONS[locus]
        topology_row = topology.get(topology_accession)
        if topology_row is None:
            raise ValueError(f"missing UniProt topology record {topology_accession}")
        span = extracellular_span_from_topology(topology_row)
        if span is None:
            raise ValueError(f"no exact extracellular span for {topology_accession}")
        topology_length = int(topology_row["Length"])
        topology_sequence = topology_row["Sequence"]
        topology_matches = sorted(
            name
            for name, sequence in aligned.items()
            if name.startswith(f"{locus}*")
            and _aligned_plain_sequence(sequence) == topology_sequence
        )
        if not topology_matches:
            raise ValueError(
                f"no exact IPD alignment sequence matches UniProt topology "
                f"reference {topology_accession} for {locus}"
            )
        reference_name = topology_matches[0]
        reference_aligned = aligned[reference_name]
        if len(topology_sequence) != topology_length:
            raise ValueError(f"invalid UniProt length for {topology_accession}")
        start_column, end_column = _reference_span_columns(
            reference_aligned, span[0], span[1]
        )
        locus_alignment[locus] = (
            aligned,
            start_column,
            end_column,
            span[0],
            span[1],
            topology_accession,
            reference_name,
        )

    grouped: MutableMapping[Tuple[str, ...], Dict[str, object]] = {}
    exclusions: List[Dict[str, str]] = []
    excluded_suffixes = {"N", "L", "S", "Q", "C", "A"}

    for full_row in full_rows:
        curated_accessions = [
            accession
            for accession in full_row["source_accession"].split(",")
            if accession in MOUSE_EXACT_ECTODOMAIN_SPANS
        ]
        if curated_accessions:
            if len(curated_accessions) != 1:
                raise ValueError("multiple curated exact ectodomain accessions")
            accession = curated_accessions[0]
            topology_row = topology.get(accession)
            expected_span = MOUSE_EXACT_ECTODOMAIN_SPANS[accession]
            if topology_row is None or topology_row.get("Sequence") != full_row["sequence"]:
                raise ValueError(f"missing exact topology sequence for {accession}")
            observed_span = extracellular_span_from_topology(topology_row)
            if observed_span != expected_span:
                raise ValueError(
                    f"curated ectodomain span changed for {accession}: {observed_span}"
                )
            segment = full_row["sequence"][expected_span[0] - 1 : expected_span[1]]
            if not segment or set(segment) - AA:
                raise ValueError(f"invalid curated ectodomain for {accession}")
            _add_af3_candidate(
                grouped,
                full_row,
                segment,
                expected_span[0],
                expected_span[1],
                boundary_method=(
                    "reviewed UniProt exact extracellular topology on complete mature membrane sequence"
                ),
                boundary_source=f"UniProtKB {accession}",
            )
            continue
        if full_row["is_complete"] != "true":
            continue
        suffix_match = re.search(r"([A-Z])$", full_row["allele"])
        if full_row["species"] == "human" and suffix_match:
            if suffix_match.group(1) in excluded_suffixes:
                exclusions.append(
                    _af3_exclusion(
                        full_row,
                        f"HLA expression suffix {suffix_match.group(1)} is excluded from the default AF3 registry",
                    )
                )
                continue

        if full_row["species"] == "human" and full_row["locus"] in locus_alignment:
            (
                aligned,
                start_column,
                end_column,
                _,
                _,
                topology_accession,
                topology_reference_name,
            ) = locus_alignment[full_row["locus"]]
            names = [full_row["allele"], *full_row["aliases"].split(";")]
            segment_candidates: set[Tuple[str, int, int]] = set()
            for raw_name in names:
                name = raw_name.upper().removeprefix("HLA-")
                target_aligned = aligned.get(name)
                if target_aligned is None:
                    continue
                if _aligned_plain_sequence(target_aligned) != full_row["sequence"]:
                    continue
                segment = "".join(
                    residue
                    for residue in target_aligned[start_column : end_column + 1]
                    if residue in AA
                )
                start_index = full_row["sequence"].find(segment)
                if start_index < 0 or full_row["sequence"].count(segment) != 1:
                    continue
                start_1based = start_index + 1
                end_1based = start_index + len(segment)
                if _af3_segment_is_plausible(
                    full_row["sequence"],
                    segment,
                    start_1based,
                    end_1based,
                    mhc_class=full_row["mhc_class"],
                    chain_role=full_row["chain_role"],
                ):
                    segment_candidates.add((segment, start_1based, end_1based))
            if len(segment_candidates) != 1:
                reason = (
                    "no exact complete IPD alignment/topology projection"
                    if not segment_candidates
                    else "conflicting IPD alignment/topology projections"
                )
                exclusions.append(_af3_exclusion(full_row, reason))
                continue
            segment, start_1based, end_1based = next(iter(segment_candidates))
            filename = HLA_ALIGNMENT_FILES[full_row["locus"]]
            _add_af3_candidate(
                grouped,
                full_row,
                segment,
                start_1based,
                end_1based,
                boundary_method="IPD aligned-allele projection of reviewed UniProt topology",
                boundary_source=(
                    f"IPD-IMGT/HLA {HLA_RELEASE} {filename}; "
                    f"UniProtKB {topology_accession} matched to IPD {topology_reference_name}"
                ),
            )
            continue

        direct_candidates = []
        for accession in full_row["source_accession"].split(","):
            topology_row = topology.get(accession)
            if topology_row is None or topology_row.get("Sequence") != full_row["sequence"]:
                continue
            span = extracellular_span_from_topology(
                topology_row, companion=full_row["chain_role"] == "companion"
            )
            if span is None:
                continue
            segment = full_row["sequence"][span[0] - 1 : span[1]]
            if _af3_segment_is_plausible(
                full_row["sequence"],
                segment,
                span[0],
                span[1],
                mhc_class=full_row["mhc_class"],
                chain_role=full_row["chain_role"],
            ):
                direct_candidates.append((accession, segment, span[0], span[1]))
        if len(direct_candidates) != 1:
            exclusions.append(
                _af3_exclusion(full_row, "no unique exact UniProt extracellular topology span")
            )
            continue
        accession, segment, start_1based, end_1based = direct_candidates[0]
        _add_af3_candidate(
            grouped,
            full_row,
            segment,
            start_1based,
            end_1based,
            boundary_method="reviewed/annotated UniProt extracellular topology",
            boundary_source=f"UniProtKB {accession}",
        )

    return _finalize_af3_candidates(grouped), exclusions


def parse_pseudo_lines(path: Path) -> Iterator[Tuple[str, str]]:
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = re.split(r"\s+", line)
        sequence_indexes = [
            index
            for index, field in enumerate(fields)
            if len(field) == 34 and set(field.upper()) <= AA_X
        ]
        if not sequence_indexes:
            continue
        index = sequence_indexes[-1]
        sequence = fields[index].upper()
        names = [field for field_index, field in enumerate(fields) if field_index != index]
        if not names:
            continue
        yield names[0], sequence


def resolve_hla_component(
    gene: str,
    allele_text: str,
    legacy_names: Mapping[str, str],
    official_two_field: set[str],
) -> str | None:
    """Resolve compact/legacy model-artifact names to current two-field HLA names."""

    gene = gene.upper()
    text = allele_text.upper().replace("*", "").replace("_", "")
    suffix_match = re.fullmatch(r"([0-9:]+)([A-Z]?)", text)
    if suffix_match is None:
        return None
    numeric, suffix = suffix_match.groups()
    legacy_keys = [f"{gene}*{numeric}{suffix}"]
    if gene == "C":
        legacy_keys.extend(
            [
                f"CW*{numeric}{suffix}",
                f"CW*{numeric.replace(':', '')}{suffix}",
            ]
        )
    legacy = next((legacy_names[key] for key in legacy_keys if key in legacy_names), None)
    if legacy == "NONE":
        return None
    if legacy:
        seen = set()
        while legacy in legacy_names and legacy not in seen:
            seen.add(legacy)
            next_name = legacy_names[legacy]
            if next_name == "NONE":
                return None
            legacy = next_name
        return hla_two_field(legacy)

    if ":" in numeric:
        fields = numeric.split(":")
        if len(fields) < 2:
            return None
        return f"HLA-{gene}*{fields[0]}:{fields[1]}{suffix}"

    if len(numeric) < 4:
        return None
    if gene == "DPB1" and len(numeric) >= 5:
        # The 2016 NetMHCIIpan artifact removes the colon from names such as
        # DPB1*100:01.  Prefer the official 3+2 split when it exists; this is
        # distinct from the also-valid modern name DPB1*10:001.
        three_digit = f"HLA-{gene}*{numeric[:3]}:{numeric[3:]}{suffix}"
        two_digit = f"HLA-{gene}*{numeric[:2]}:{numeric[2:]}{suffix}"
        if three_digit in official_two_field:
            return three_digit
        if two_digit in official_two_field:
            return two_digit
        return three_digit
    return f"HLA-{gene}*{numeric[:2]}:{numeric[2:]}{suffix}"


def canonical_hla_compact(
    raw_name: str,
    legacy_names: Mapping[str, str],
    official_two_field: set[str],
) -> str | None:
    raw = raw_name.strip().upper().replace("HLA_", "HLA-")
    match = re.fullmatch(r"(?:HLA-)?(A|B|C|E|F|G)\*?([0-9:]+[A-Z]?)", raw)
    if not match:
        return None
    return resolve_hla_component(
        match.group(1), match.group(2), legacy_names, official_two_field
    )


def canonical_mouse_class_i(raw_name: str) -> str | None:
    match = re.fullmatch(r"H[-_]?2[-_]?([DKL])[-_]?([A-Za-z0-9]+)", raw_name.strip())
    if not match:
        return None
    suffix = match.group(2)
    if suffix:
        suffix = suffix[0].lower() + suffix[1:]
    return f"H-2-{match.group(1).upper()}{suffix}"


def build_pseudo_i(
    path: Path,
    legacy_names: Mapping[str, str],
    official_two_field: set[str],
) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    source_hash = sha256_path(path)
    grouped: MutableMapping[Tuple[str, str], set[str]] = defaultdict(set)
    species_by_name = {}
    for raw_name, sequence in parse_pseudo_lines(path):
        canonical = canonical_hla_compact(raw_name, legacy_names, official_two_field)
        species = "human"
        if canonical is None:
            canonical = canonical_mouse_class_i(raw_name)
            species = "mouse"
        if canonical is None:
            continue
        grouped[(canonical, sequence)].add(raw_name)
        species_by_name[canonical] = species

    rows = []
    for (canonical, sequence), aliases in grouped.items():
        aliases.discard(canonical)
        rows.append(
            {
                "record_id": record_id("pseudo", PSEUDO_I_DEFINITION, canonical, sequence),
                "canonical_name": canonical,
                "aliases": ";".join(sorted(aliases)),
                "species": species_by_name[canonical],
                "mhc_class": "I",
                "definition_id": PSEUDO_I_DEFINITION,
                "pseudo_sequence": sequence,
                "alpha_allele": canonical,
                "beta_allele": "",
                "source": "PUFFIN data/MHC_pseudo.dat (NetMHCpan-format mirror)",
                "source_version": f"sha256:{source_hash[:16]}",
            }
        )
    by_name: MutableMapping[str, List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_name[row["canonical_name"]].append(row)
    accepted = []
    quarantined = []
    for name, candidates in by_name.items():
        if len(candidates) == 1:
            accepted.append(candidates[0])
            continue
        # The mirror contains both historical compact and current colon names.
        # When they disagree after official renaming, prefer a row explicitly
        # carrying the current colon-delimited spelling and quarantine the
        # historical row.  A tie is not safe to resolve automatically.
        canonical_artifact_name = name.replace("*", "")
        scores = []
        for row in candidates:
            aliases = row["aliases"].split(";")
            scores.append(
                100 * int(canonical_artifact_name in aliases)
                + sum(":" in alias for alias in aliases)
            )
        best = max(scores)
        winners = [row for row, score in zip(candidates, scores) if score == best]
        if len(winners) != 1:
            detail = {row["pseudo_sequence"]: row["aliases"] for row in candidates}
            raise ValueError(f"unresolved MHC-I pseudosequence conflict for {name}: {detail}")
        winner = winners[0]
        accepted.append(winner)
        for row in candidates:
            if row is winner:
                continue
            conflict = dict(row)
            conflict["reason"] = (
                "historical compact-name row conflicts with the current colon-delimited row"
            )
            quarantined.append(conflict)
    return accepted, quarantined


def hla_class_ii_components(
    raw_name: str,
    legacy_names: Mapping[str, str],
    official_two_field: set[str],
) -> Tuple[str, str] | None:
    raw = raw_name.upper().replace("HLA-", "").replace("HLA_", "")
    matches = re.findall(
        r"(DRA|DRB[1-5]|DQA1|DQB1|DPA1|DPB1)[_*]?([0-9]+(?::[0-9]+)?[A-Z]?)",
        raw,
    )
    if not matches:
        return None
    names = [
        resolve_hla_component(gene, allele, legacy_names, official_two_field)
        for gene, allele in matches
    ]
    if any(name is None for name in names):
        return None
    resolved = [str(name) for name in names]
    if len(resolved) == 1 and matches[0][0].startswith("DRB"):
        return "HLA-DRA*01:01", resolved[0]
    alpha = [name for name in resolved if any(tag in name for tag in ("-DRA*", "-DQA1*", "-DPA1*"))]
    beta = [name for name in resolved if any(tag in name for tag in ("-DRB", "-DQB1*", "-DPB1*"))]
    if len(alpha) == 1 and len(beta) == 1:
        return alpha[0], beta[0]
    return None


def canonical_mouse_class_ii(raw_name: str) -> str | None:
    compact = re.sub(r"[^A-Za-z0-9]", "", raw_name)
    match = re.fullmatch(r"H2I([AE])([A-Za-z0-9]+)", compact, re.IGNORECASE)
    if not match:
        return None
    suffix = match.group(2)
    suffix = suffix[0].lower() + suffix[1:] if suffix else suffix
    return f"H-2-I{match.group(1).upper()}{suffix}"


def _unique_modal_segment(values: Iterable[str], label: str) -> str:
    counts = Counter(values)
    ranked = counts.most_common()
    if not ranked:
        raise ValueError(f"no pseudosequence segments for {label}")
    if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
        raise ValueError(f"unresolved tied pseudosequence conflict for {label}: {counts}")
    return ranked[0][0]


def build_pseudo_ii(
    path: Path,
    legacy_names: Mapping[str, str],
    official_two_field: set[str],
) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    source_hash = sha256_path(path)
    grouped: MutableMapping[Tuple[str, str, str, str], set[str]] = defaultdict(set)
    species_by_name = {}
    for raw_name, sequence in parse_pseudo_lines(path):
        mouse_name = canonical_mouse_class_ii(raw_name)
        if mouse_name is not None:
            if mouse_name not in MOUSE_PAIRS:
                continue
            alpha, beta = MOUSE_PAIRS[mouse_name]
            canonical = mouse_name
            species = "mouse"
        else:
            components = hla_class_ii_components(
                raw_name, legacy_names, official_two_field
            )
            if components is None:
                continue
            alpha, beta = components
            canonical = f"{alpha}/{beta}"
            species = "human"
        grouped[(canonical, sequence, alpha, beta)].add(raw_name)
        species_by_name[canonical] = species

    rows = []
    for (canonical, sequence, alpha, beta), aliases in grouped.items():
        aliases.update({f"{alpha}+{beta}"})
        if beta.startswith("HLA-DRB"):
            aliases.update({beta, beta.removeprefix("HLA-")})
        aliases.discard(canonical)
        rows.append(
            {
                "record_id": record_id("pseudo", PSEUDO_II_DEFINITION, canonical, sequence),
                "canonical_name": canonical,
                "aliases": ";".join(sorted(aliases)),
                "species": species_by_name[canonical],
                "mhc_class": "II",
                "definition_id": PSEUDO_II_DEFINITION,
                "pseudo_sequence": sequence,
                "alpha_allele": alpha,
                "beta_allele": beta,
                "source": "PUFFIN data/pseudosequence.2016.all.X.dat (NetMHCIIpan-format mirror)",
                "source_version": f"sha256:{source_hash[:16]}",
            }
        )
    alpha_segments: MutableMapping[str, List[str]] = defaultdict(list)
    beta_segments: MutableMapping[str, List[str]] = defaultdict(list)
    for row in rows:
        alpha_segments[row["alpha_allele"]].append(row["pseudo_sequence"][:15])
        beta_segments[row["beta_allele"]].append(row["pseudo_sequence"][15:])
    expected_alpha = {
        allele: _unique_modal_segment(values, f"alpha {allele}")
        for allele, values in alpha_segments.items()
    }
    expected_beta = {
        allele: _unique_modal_segment(values, f"beta {allele}")
        for allele, values in beta_segments.items()
    }

    accepted = []
    quarantined = []
    for row in rows:
        reasons = []
        if row["pseudo_sequence"][:15] != expected_alpha[row["alpha_allele"]]:
            reasons.append("alpha segment conflicts with the modal segment for this allele")
        if row["pseudo_sequence"][15:] != expected_beta[row["beta_allele"]]:
            reasons.append("beta segment conflicts with the modal segment for this allele")
        if reasons:
            conflict = dict(row)
            conflict["reason"] = "; ".join(reasons)
            quarantined.append(conflict)
        else:
            accepted.append(row)

    by_pair: MutableMapping[str, set[str]] = defaultdict(set)
    for row in accepted:
        by_pair[row["canonical_name"]].add(row["pseudo_sequence"])
    unresolved = {name: values for name, values in by_pair.items() if len(values) > 1}
    if unresolved:
        raise ValueError(f"unresolved MHC-II pseudosequence conflicts: {unresolved}")
    return accepted, quarantined


def build_pairs() -> List[Dict[str, str]]:
    rows = []
    for canonical, (alpha, beta) in MOUSE_PAIRS.items():
        superscript = canonical.replace("H-2-", "H2-")
        rows.append(
            {
                "canonical_name": canonical,
                "aliases": superscript,
                "species": "mouse",
                "mhc_class": "II",
                "alpha_or_heavy": alpha,
                "beta_or_companion": beta,
                "source": "IMGT mouse MHC nomenclature; chains mapped to explicit UniProt accessions",
            }
        )
    return rows


def write_csv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def count_by(rows: Iterable[Mapping[str, str]], *keys: str) -> Mapping[str, int]:
    counts: MutableMapping[str, int] = defaultdict(int)
    for row in rows:
        label = "/".join(row[key] for key in keys)
        counts[label] += 1
    return dict(sorted(counts.items()))


def build(args: argparse.Namespace) -> None:
    cache_dir = args.cache_dir.resolve()
    output_dir = args.output_dir.resolve()
    if bool(args.hla_fasta) != bool(args.hla_status):
        raise ValueError("--hla-fasta and --hla-status must be supplied together")
    if args.hla_fasta:
        hla_path = args.hla_fasta.resolve()
        hla_status_path = args.hla_status.resolve()
        if not hla_path.exists():
            raise FileNotFoundError(hla_path)
        if not hla_status_path.exists():
            raise FileNotFoundError(hla_status_path)
    else:
        hla_path = download(
            HLA_FASTA_URL,
            cache_dir / f"IMGTHLA-{HLA_RELEASE}-hla_prot.fasta",
            offline=args.offline,
        )
        hla_status_path = download(
            HLA_STATUS_URL,
            cache_dir / f"IMGTHLA-{HLA_RELEASE}-Allele_status.txt",
            offline=args.offline,
        )
    hla_2009_path = download(
        HLA_2009_URL,
        cache_dir / f"IMGTHLA-{HLA_RELEASE}-Nomenclature_2009.txt",
        offline=args.offline,
    )
    hla_deleted_path = download(
        HLA_DELETED_URL,
        cache_dir / f"IMGTHLA-{HLA_RELEASE}-Deleted_alleles.txt",
        offline=args.offline,
    )

    pseudo_i_path = download(PSEUDO_I_URL, cache_dir / "MHC_pseudo.dat", offline=args.offline)
    pseudo_ii_path = download(
        PSEUDO_II_URL,
        cache_dir / "pseudosequence.2016.all.X.dat",
        offline=args.offline,
    )
    uniprot_path = fetch_uniprot(cache_dir, offline=args.offline)
    uniprot_topology_path = fetch_uniprot_topology(cache_dir, offline=args.offline)
    hla_alignment_paths = fetch_hla_alignments(cache_dir, offline=args.offline)

    legacy_names = dict(parse_nomenclature_2009(hla_2009_path))
    legacy_names.update(parse_deleted_renames(hla_deleted_path))
    official_two_field = official_two_field_names(hla_status_path)
    full_rows = build_human_full(hla_path, hla_status_path) + build_mouse_and_b2m(uniprot_path)
    af3_rows, excluded_af3_rows = build_af3_sequences(
        full_rows, hla_alignment_paths, uniprot_topology_path
    )
    pseudo_i_rows, quarantined_pseudo_i = build_pseudo_i(
        pseudo_i_path, legacy_names, official_two_field
    )
    pseudo_ii_rows, quarantined_pseudo_ii = build_pseudo_ii(
        pseudo_ii_path, legacy_names, official_two_field
    )
    quarantined_pseudos = quarantined_pseudo_i + quarantined_pseudo_ii
    pseudo_rows = pseudo_i_rows + pseudo_ii_rows
    pair_rows = build_pairs()
    full_rows.sort(key=lambda row: (row["species"], row["mhc_class"], row["allele"], row["sequence"]))
    af3_rows.sort(
        key=lambda row: (
            row["species"],
            row["mhc_class"],
            row["allele"],
            row["sequence"],
        )
    )
    excluded_af3_rows.sort(
        key=lambda row: (row["species"], row["mhc_class"], row["allele"], row["reason"])
    )
    pseudo_rows.sort(
        key=lambda row: (
            row["mhc_class"],
            row["species"],
            row["canonical_name"],
            row["pseudo_sequence"],
        )
    )
    pair_rows.sort(key=lambda row: row["canonical_name"])

    write_csv(output_dir / "full_sequences.csv", FULL_FIELDS, full_rows)
    write_csv(output_dir / "af3_sequences.csv", AF3_FIELDS, af3_rows)
    write_csv(output_dir / "pseudo_sequences.csv", PSEUDO_FIELDS, pseudo_rows)
    write_csv(output_dir / "molecule_pairs.csv", PAIR_FIELDS, pair_rows)
    write_csv(
        output_dir / "quarantined_pseudosequences.csv",
        CONFLICT_FIELDS,
        quarantined_pseudos,
    )
    write_csv(
        output_dir / "excluded_af3_sequences.csv",
        AF3_EXCLUSION_FIELDS,
        excluded_af3_rows,
    )

    metadata = {
        "schema_version": 3,
        "generated_on": date.today().isoformat(),
        "sequence_kind": (
            "precursor_full (default full lookup); partial_deposited_protein (opt-in); "
            f"{AF3_DEFINITION} (mature soluble extracellular chain)"
        ),
        "counts": {
            "full_total": len(full_rows),
            "full_complete": sum(row["is_complete"] == "true" for row in full_rows),
            "full_by_species_class": count_by(full_rows, "species", "mhc_class"),
            "af3_total": len(af3_rows),
            "af3_excluded_complete": len(excluded_af3_rows),
            "af3_by_species_class": count_by(af3_rows, "species", "mhc_class"),
            "pseudo_total": len(pseudo_rows),
            "pseudo_quarantined": len(quarantined_pseudos),
            "pseudo_by_species_class": count_by(pseudo_rows, "species", "mhc_class"),
            "pairs": len(pair_rows),
        },
        "sources": {
            "human_full": {
                "name": f"IPD-IMGT/HLA {HLA_RELEASE} protein records with official sequence status",
                "url": "https://www.ebi.ac.uk/ipd/imgt/hla/",
                "commit": HLA_COMMIT,
                "sha256": sha256_path(hla_path),
                "status_url": HLA_STATUS_URL,
                "status_sha256": sha256_path(hla_status_path),
                "legacy_nomenclature_url": HLA_2009_URL,
                "legacy_nomenclature_sha256": sha256_path(hla_2009_path),
                "deleted_alleles_url": HLA_DELETED_URL,
                "deleted_alleles_sha256": sha256_path(hla_deleted_path),
            },
            "mouse_and_b2m_full": {
                "name": "UniProtKB accession batch",
                "url": "https://rest.uniprot.org/",
                "sha256": sha256_path(uniprot_path),
                "accessions": sorted({str(item["accession"]) for item in MOUSE_PROTEINS}),
            },
            "af3_topology": {
                "name": "UniProtKB topology feature batch",
                "url": "https://rest.uniprot.org/",
                "sha256": sha256_path(uniprot_topology_path),
                "human_reference_accessions": dict(sorted(HUMAN_TOPOLOGY_ACCESSIONS.items())),
            },
            "af3_human_alignments": {
                "name": f"IPD-IMGT/HLA {HLA_RELEASE} protein alignments",
                "commit": HLA_COMMIT,
                "files": {
                    filename: {
                        "url": (
                            f"https://raw.githubusercontent.com/ANHIG/IMGTHLA/{HLA_COMMIT}/"
                            f"alignments/{filename}"
                        ),
                        "sha256": sha256_path(path),
                    }
                    for filename, path in sorted(hla_alignment_paths.items())
                },
            },
            "pseudo_i": {
                "name": "PUFFIN MHC_pseudo.dat",
                "url": PSEUDO_I_URL,
                "commit": PUFFIN_COMMIT,
                "sha256": sha256_path(pseudo_i_path),
            },
            "pseudo_ii": {
                "name": "PUFFIN pseudosequence.2016.all.X.dat",
                "url": PSEUDO_II_URL,
                "commit": PUFFIN_COMMIT,
                "sha256": sha256_path(pseudo_ii_path),
            },
        },
        "pseudo_definitions": {
            PSEUDO_I_DEFINITION: {
                "length": 34,
                "chains": ["MHC-I heavy"],
                "scope_note": "Defined from classical peptide-presenting MHC-I; nonclassical loci/species are extrapolations.",
                "mature_reference_positions_1_based": [7, 9, 24, 45, 59, 62, 63, 66, 67, 69, 70, 73, 74, 76, 77, 80, 81, 84, 95, 97, 99, 114, 116, 118, 143, 147, 150, 152, 156, 158, 159, 163, 167, 171],
            },
            PSEUDO_II_DEFINITION: {
                "length": 34,
                "chains": ["MHC-II alpha (15)", "MHC-II beta (19)"],
                "alpha_mature_reference_positions_1_based": [9, 11, 22, 24, 31, 52, 53, 58, 59, 61, 65, 66, 68, 72, 73],
                "beta_mature_reference_positions_1_based": [9, 11, 13, 26, 28, 30, 47, 57, 67, 70, 71, 74, 77, 78, 81, 85, 86, 89, 90],
                "alphabet_note": "X is retained as the source artifact's reference-alignment gap/unknown token; it is not a canonical amino acid.",
                "pairing_note": "DR beta shorthand assumes DRA*01:01. DQ and DP require explicit alpha/beta sequence pairs; artifact presence does not establish biological haplotype evidence.",
            },
        },
        "af3_definition": {
            "definition_id": AF3_DEFINITION,
            "description": (
                "Mature soluble extracellular chain: signal peptide, transmembrane helix "
                "and cytoplasmic tail removed; all extracellular folding domains and the "
                "membrane-proximal connector retained."
            ),
            "official_vs_engineering_note": (
                "AlphaFold 3 specifies separate protein entities but no pMHC-TCR crop. "
                "This ectodomain profile is an auditable engineering default."
            ),
            "human_boundary_method": (
                "Reviewed UniProt locus topology projected through the pinned official "
                "IPD protein alignment to each exact deposited allele sequence."
            ),
            "mouse_and_b2m_boundary_method": (
                "Exact accession-specific UniProt topology features; P14437 is a "
                "reviewed signal-free mature membrane chain whose residues 1-195 are "
                "the exact extracellular domain."
            ),
            "excluded_expression_suffixes": ["N", "L", "S", "Q", "C", "A"],
            "excluded_file": "excluded_af3_sequences.csv",
            "pseudo34_is_not_a_structure_chain": True,
        },
        "quarantine": {
            "count": len(quarantined_pseudos),
            "file": "quarantined_pseudosequences.csv",
            "rule": "Prefer an explicit current colon-delimited artifact name over a conflicting historical alias; for MHC-II also enforce one modal alpha-15 and beta-19 segment per allele. Ties fail the build.",
        },
        "completeness_rules": {
            "human": {
                "method": "Allele_status.txt Partial == Full and protein alphabet contains no X",
                "release": HLA_RELEASE,
                "commit": HLA_COMMIT,
            },
            "mouse": {
                "method": "explicit UniProt accession, no Fragment flag, valid alphabet, plus conservative chain-length screen",
                "MHC-I_heavy_min_length": 320,
                "MHC-II_chain_min_length": 240,
                "B2M_min_length": 110,
            },
            "warning": "Complete means the bundled source represents the complete coding protein record. It is not an assertion of cell-surface expression; biologically truncating/null alleles can be complete records.",
        },
    }
    with (output_dir / "metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    # Fallback construction intentionally runs only after the strict registry
    # and its preliminary metadata are complete.  It reads those immutable
    # assets and writes separate indexes.
    from build_fallback_reference import build_fallback_files, fallback_metadata_fields

    fallback_counts = build_fallback_files(output_dir)
    metadata.update(fallback_metadata_fields())
    metadata["counts"].update(fallback_counts)
    with (output_dir / "metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    print(json.dumps(metadata["counts"], indent=2, sort_keys=True))


def parser() -> argparse.ArgumentParser:
    project_root = Path(__file__).resolve().parents[1]
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--hla-fasta",
        type=Path,
        help="override the pinned official hla_prot.fasta (requires --hla-status)",
    )
    result.add_argument(
        "--hla-status",
        type=Path,
        help="matching Allele_status.txt override (requires --hla-fasta)",
    )
    result.add_argument(
        "--cache-dir", type=Path, default=project_root / ".cache/reference"
    )
    result.add_argument(
        "--output-dir", type=Path, default=project_root / "src/mhctransform/data"
    )
    result.add_argument("--offline", action="store_true")
    return result


if __name__ == "__main__":
    build(parser().parse_args())

