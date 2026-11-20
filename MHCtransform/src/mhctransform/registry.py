"""Reference-backed, non-lossy MHC lookup API."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from importlib import resources
from itertools import product
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Mapping, Sequence, Tuple

from .errors import MHCTransformError
from .models import (
    AF3MoleculeCandidate,
    AF3SequenceCandidate,
    FallbackAF3SequenceCandidate,
    FallbackPBRSequenceCandidate,
    FallbackSequenceCandidate,
    MoleculeCandidate,
    PBRMoleculeCandidate,
    PBRSequenceCandidate,
    PseudoCandidate,
    QueryResult,
    SequenceCandidate,
    result_status,
)
from .fallback import (
    IMPUTED_AF3_DEFINITION,
    IMPUTED_EVIDENCE,
    PSEUDO_SCAFFOLD_EVIDENCE,
    evidence_allowed,
    evidence_priority,
    normalize_fallback_policy,
)
from .names import (
    canonical_hla,
    exact_name_key,
    normalize_mhc_class,
    normalize_sequence,
    two_field_hla,
)
from .pbr import PBR_DEFINITIONS


PSEUDO_I_DEFINITION = "netmhcpan-34-v1"
PSEUDO_II_DEFINITION = "netmhciipan-34-v1"
AF3_DEFINITION = "af3-soluble-extracellular-v1"
CURRENT_Q_COMPATIBILITY_ALIASES: Mapping[str, str] = {
    "HLA-DQA1*01:07Q": "HLA-DQA1*01:07",
    "HLA-DRB1*15:13Q": "HLA-DRB1*15:13",
}


def _split_aliases(value: str) -> Tuple[str, ...]:
    return tuple(part for part in value.split(";") if part)


def _index_names(primary: str, aliases: Iterable[str]) -> Tuple[str, ...]:
    """Return stored names plus two explicit NetMHC historical aliases."""

    names = {primary, *aliases}
    keys = {exact_name_key(name) for name in names}
    for current, historical in CURRENT_Q_COMPATIBILITY_ALIASES.items():
        if exact_name_key(current) in keys:
            names.add(historical)
    return tuple(sorted(names))


def _deduplicate(candidates: Iterable[object]) -> Tuple[object, ...]:
    by_id = {}
    for candidate in candidates:
        by_id[getattr(candidate, "record_id")] = candidate
    return tuple(by_id[key] for key in sorted(by_id))


def _candidate_evidence(candidate: object) -> str:
    if isinstance(candidate, SequenceCandidate) and not candidate.is_complete:
        return "source_partial"
    return str(getattr(candidate, "evidence_level", "source_exact"))


def _mhcii_family(locus: str, role: str) -> str | None:
    locus = locus.upper()
    if role == "alpha":
        if locus == "DRA":
            return "DR"
        if locus == "DQA1":
            return "DQ"
        if locus == "DPA1":
            return "DP"
        if locus == "H2-AA":
            return "IA"
        if locus == "H2-EA":
            return "IE"
    if role == "beta":
        if re.fullmatch(r"DRB[1-5]", locus):
            return "DR"
        if locus == "DQB1":
            return "DQ"
        if locus == "DPB1":
            return "DP"
        if locus == "H2-AB1":
            return "IA"
        if locus == "H2-EB1":
            return "IE"
    return None


class MHCRegistry:
    """Load the bundled human/mouse MHC registry.

    The registry never guesses a nearest allele.  Reverse queries return every
    exact reference match, which is essential because MHC pseudosequences and
    even protein sequences are frequently shared by several allele names.
    """

    def __init__(self, data_dir: str | Path | None = None) -> None:
        if data_dir is None:
            data_root = resources.files("mhctransform").joinpath("data")
        else:
            data_root = Path(data_dir)

        self._sequences = self._load_sequences(data_root.joinpath("full_sequences.csv"))
        self._af3_sequences = self._load_af3_sequences(
            data_root.joinpath("af3_sequences.csv")
        )
        self._af3_exclusions = self._load_pairs(
            data_root.joinpath("excluded_af3_sequences.csv")
        )
        self._pseudos = self._load_pseudos(data_root.joinpath("pseudo_sequences.csv"))
        pbr_path = data_root.joinpath("pbr_sequences.csv")
        self._pbr_sequences = (
            self._load_pbr_sequences(pbr_path) if pbr_path.is_file() else ()
        )
        self._pairs = self._load_pairs(data_root.joinpath("molecule_pairs.csv"))
        fallback_full_path = data_root.joinpath("fallback_full_sequences.csv")
        fallback_af3_path = data_root.joinpath("fallback_af3_sequences.csv")
        fallback_pbr_path = data_root.joinpath("fallback_pbr_sequences.csv")
        self._fallback_sequences = (
            self._load_fallback_sequences(fallback_full_path, af3=False)
            if fallback_full_path.is_file()
            else ()
        )
        self._fallback_af3_sequences = (
            self._load_fallback_sequences(fallback_af3_path, af3=True)
            if fallback_af3_path.is_file()
            else ()
        )
        self._fallback_pbr_sequences = (
            self._load_fallback_pbr_sequences(fallback_pbr_path)
            if fallback_pbr_path.is_file()
            else ()
        )
        pseudo_pbr_projection_path = data_root.joinpath("pseudo_pbr_projection.csv")
        self._pseudo_pbr_projections = (
            self._load_pairs(pseudo_pbr_projection_path)
            if pseudo_pbr_projection_path.is_file()
            else ()
        )
        with data_root.joinpath("metadata.json").open("r", encoding="utf-8") as handle:
            self.metadata: Mapping[str, object] = json.load(handle)

        self._sequence_name_index: DefaultDict[str, List[SequenceCandidate]] = defaultdict(list)
        self._sequence_value_index: DefaultDict[str, List[SequenceCandidate]] = defaultdict(list)
        for candidate in self._sequences:
            for name in _index_names(candidate.allele, candidate.aliases):
                self._sequence_name_index[exact_name_key(name)].append(candidate)
            self._sequence_value_index[candidate.sequence].append(candidate)

        self._af3_name_index: DefaultDict[str, List[AF3SequenceCandidate]] = defaultdict(list)
        self._af3_value_index: DefaultDict[str, List[AF3SequenceCandidate]] = defaultdict(list)
        for candidate in self._af3_sequences:
            for name in _index_names(candidate.allele, candidate.aliases):
                self._af3_name_index[exact_name_key(name)].append(candidate)
            self._af3_value_index[candidate.sequence].append(candidate)

        self._pbr_name_index: DefaultDict[str, List[PBRSequenceCandidate]] = defaultdict(list)
        self._pbr_value_index: DefaultDict[str, List[PBRSequenceCandidate]] = defaultdict(list)
        for candidate in self._pbr_sequences:
            for name in _index_names(candidate.allele, candidate.aliases):
                self._pbr_name_index[exact_name_key(name)].append(candidate)
            self._pbr_value_index[candidate.sequence].append(candidate)

        self._af3_exclusion_index: DefaultDict[str, List[Mapping[str, str]]] = defaultdict(list)
        for exclusion in self._af3_exclusions:
            for name in _index_names(
                exclusion["allele"], _split_aliases(exclusion["aliases"])
            ):
                self._af3_exclusion_index[exact_name_key(name)].append(exclusion)

        self._pseudo_name_index: DefaultDict[str, List[PseudoCandidate]] = defaultdict(list)
        self._pseudo_value_index: DefaultDict[str, List[PseudoCandidate]] = defaultdict(list)
        for candidate in self._pseudos:
            for name in _index_names(candidate.canonical_name, candidate.aliases):
                self._pseudo_name_index[exact_name_key(name)].append(candidate)
            self._pseudo_value_index[candidate.pseudo_sequence].append(candidate)

        self._pair_index: DefaultDict[str, List[Mapping[str, str]]] = defaultdict(list)
        for pair in self._pairs:
            for name in (pair["canonical_name"], *_split_aliases(pair["aliases"])):
                self._pair_index[exact_name_key(name)].append(pair)

        self._fallback_sequence_name_index: DefaultDict[
            str, List[FallbackSequenceCandidate]
        ] = defaultdict(list)
        self._fallback_sequence_value_index: DefaultDict[
            str, List[FallbackSequenceCandidate]
        ] = defaultdict(list)
        for candidate in self._fallback_sequences:
            for name in _index_names(candidate.allele, candidate.aliases):
                self._fallback_sequence_name_index[exact_name_key(name)].append(candidate)
            self._fallback_sequence_value_index[candidate.sequence].append(candidate)

        self._fallback_af3_name_index: DefaultDict[
            str, List[FallbackAF3SequenceCandidate]
        ] = defaultdict(list)
        self._fallback_af3_value_index: DefaultDict[
            str, List[FallbackAF3SequenceCandidate]
        ] = defaultdict(list)
        for candidate in self._fallback_af3_sequences:
            for name in _index_names(candidate.allele, candidate.aliases):
                self._fallback_af3_name_index[exact_name_key(name)].append(candidate)
            self._fallback_af3_value_index[candidate.sequence].append(candidate)

        self._fallback_pbr_name_index: DefaultDict[
            str, List[FallbackPBRSequenceCandidate]
        ] = defaultdict(list)
        self._fallback_pbr_value_index: DefaultDict[
            str, List[FallbackPBRSequenceCandidate]
        ] = defaultdict(list)
        for candidate in self._fallback_pbr_sequences:
            for name in _index_names(candidate.allele, candidate.aliases):
                self._fallback_pbr_name_index[exact_name_key(name)].append(candidate)
            self._fallback_pbr_value_index[candidate.sequence].append(candidate)

        if self._pseudo_pbr_projections:
            self._validate_pseudo_pbr_projections()

        # The projection table is a build-time certificate, not another
        # sequence source.  It records the one strict or fallback PBR record
        # whose projected contact coordinates reproduce each pseudo component.
        # Indexing by the component name prevents a strict record certified for
        # one allele from accidentally satisfying a conflicting allele that
        # happens to share the same deposited protein sequence.
        self._pseudo_pbr_projection_name_index: DefaultDict[
            str, List[Mapping[str, str]]
        ] = defaultdict(list)
        projection_by_allele: DefaultDict[str, List[Mapping[str, str]]] = defaultdict(list)
        for projection in self._pseudo_pbr_projections:
            projection_by_allele[exact_name_key(projection["allele"])].append(projection)
            for name in _index_names(projection["allele"], ()):
                self._pseudo_pbr_projection_name_index[exact_name_key(name)].append(
                    projection
                )

        # Preserve explicit historical/legacy aliases for MHC-I pseudo names.
        # MHC-II projection queries are component queries and are indexed by
        # alpha_allele/beta_allele directly above.
        for pseudo in self._pseudos:
            if pseudo.mhc_class != "I":
                continue
            rows = projection_by_allele.get(exact_name_key(pseudo.canonical_name), ())
            for alias in _index_names(pseudo.canonical_name, pseudo.aliases):
                self._pseudo_pbr_projection_name_index[exact_name_key(alias)].extend(rows)

    def _validate_pseudo_pbr_projections(self) -> None:
        """Reject stale or internally inconsistent pseudo/PBR certificates.

        The projection CSV spans three generated registries.  Validating it at
        load time prevents a partially rebuilt data directory from making
        ``require_pseudo_match=True`` trust a certificate for an older pseudo
        or PBR file.
        """

        expected: DefaultDict[
            Tuple[str, str, str, str], Dict[str, object]
        ] = defaultdict(lambda: {"segments": set(), "record_ids": set()})
        for pseudo in self._pseudos:
            if pseudo.mhc_class == "I":
                components = (
                    ("heavy", pseudo.canonical_name, pseudo.pseudo_sequence, 1),
                )
            else:
                components = (
                    ("alpha", pseudo.alpha_allele, pseudo.pseudo_sequence[:15], 1),
                    ("beta", pseudo.beta_allele, pseudo.pseudo_sequence[15:], 16),
                )
            for role, allele, segment, start in components:
                key = (pseudo.species, pseudo.mhc_class, role, exact_name_key(allele))
                item = expected[key]
                item["segments"].add((segment, start))  # type: ignore[union-attr]
                item["record_ids"].add(pseudo.record_id)  # type: ignore[union-attr]

        strict_by_id = {item.record_id: item for item in self._pbr_sequences}
        fallback_by_id = {
            item.record_id: item for item in self._fallback_pbr_sequences
        }
        seen: Dict[Tuple[str, str, str, str], str] = {}

        def fail(detail: str) -> None:
            raise MHCTransformError(
                "invalid or stale pseudo/PBR projection certificate: "
                f"{detail}. Rebuild the fallback/PBR reference data."
            )

        for row in self._pseudo_pbr_projections:
            key = (
                row["species"],
                row["mhc_class"],
                row["chain_role"],
                exact_name_key(row["allele"]),
            )
            if key in seen:
                fail(f"duplicate component {row['allele']}")
            seen[key] = row["projection_id"]
            target = expected.get(key)
            if target is None:
                fail(f"unknown component {row['allele']}")
            segments = target["segments"]
            if len(segments) != 1:
                fail(f"conflicting current pseudo segments for {row['allele']}")
            segment, start = next(iter(segments))  # type: ignore[arg-type]
            if row["pseudo_segment"] != segment:
                fail(f"pseudo segment changed for {row['allele']}")
            if int(row["pseudo_start_1based"]) != start:
                fail(f"pseudo offset changed for {row['allele']}")

            source = row["pbr_source"]
            source_index = strict_by_id if source == "strict" else fallback_by_id
            if source not in {"strict", "fallback"}:
                fail(f"unknown PBR source {source!r} for {row['allele']}")
            candidate = source_index.get(row["pbr_record_id"])
            if candidate is None:
                fail(f"missing PBR record for {row['allele']}")
            if (
                candidate.species != row["species"]
                or candidate.mhc_class != row["mhc_class"]
                or candidate.chain_role != row["chain_role"]
                or candidate.definition_id != row["definition_id"]
            ):
                fail(f"PBR metadata changed for {row['allele']}")

            raw_coordinates = [
                value
                for value in row["contact_positions_1based"].split(";")
                if value
            ]
            try:
                coordinates = tuple(int(value) for value in raw_coordinates)
            except ValueError:
                fail(f"invalid PBR coordinates for {row['allele']}")
                return
            if (
                len(coordinates) != len(segment)
                or tuple(sorted(set(coordinates))) != coordinates
                or not coordinates
                or coordinates[0] < 1
                or coordinates[-1] > len(candidate.sequence)
            ):
                fail(f"incomplete PBR coordinates for {row['allele']}")
            projected = "".join(
                candidate.sequence[position - 1] for position in coordinates
            )
            if projected != row["projected_segment"]:
                fail(f"projected PBR residues changed for {row['allele']}")
            if "X" in projected:
                fail(f"PBR projection contains X for {row['allele']}")
            if any(
                pseudo_residue != "X" and pseudo_residue != pbr_residue
                for pseudo_residue, pbr_residue in zip(segment, projected)
            ):
                fail(f"non-X PBR mismatch for {row['allele']}")
            wildcard_positions = ";".join(
                str(start + index)
                for index, residue in enumerate(segment)
                if residue == "X"
            )
            if row["wildcard_pseudo_positions"] != wildcard_positions:
                fail(f"X wildcard positions changed for {row['allele']}")
            expected_coverage = (
                "exact_non_x_with_explicit_wildcards"
                if "X" in segment
                else "exact"
            )
            if row["coverage_status"] != expected_coverage:
                fail(f"coverage status changed for {row['allele']}")

            try:
                offsets = tuple(
                    int(value) for value in row["coordinate_offsets"].split(";")
                )
                shifted = tuple(
                    int(value)
                    for value in row["shifted_pseudo_positions"].split(";")
                    if value
                )
                gaps = tuple(
                    int(value)
                    for value in row["alignment_gap_pseudo_positions"].split(";")
                    if value
                )
                maximum_shift = int(row["max_abs_shift"])
            except (KeyError, ValueError):
                fail(f"invalid coordinate audit for {row['allele']}")
                return
            if len(offsets) != len(segment) or any(abs(value) > 1 for value in offsets):
                fail(f"out-of-contract coordinate offsets for {row['allele']}")
            expected_shifted = tuple(
                start + index
                for index, offset in enumerate(offsets)
                if offset
            )
            if shifted != expected_shifted:
                fail(f"shifted coordinate audit changed for {row['allele']}")
            valid_pseudo_positions = set(range(start, start + len(segment)))
            if not set(gaps).issubset(valid_pseudo_positions):
                fail(f"alignment-gap audit is out of range for {row['allele']}")
            if maximum_shift != max((abs(value) for value in offsets), default=0):
                fail(f"maximum coordinate shift changed for {row['allele']}")
            expected_method = (
                "global_alignment_gap_interpolation_and_bounded_recovery"
                if gaps and shifted
                else (
                    "global_alignment_gap_interpolation"
                    if gaps
                    else (
                        "global_alignment_bounded_recovery"
                        if shifted
                        else "global_alignment"
                    )
                )
            )
            if row["coordinate_method"] != expected_method:
                fail(f"coordinate method changed for {row['allele']}")

            record_ids = target["record_ids"]
            if int(row["pseudo_record_count"]) != len(record_ids):  # type: ignore[arg-type]
                fail(f"pseudo usage count changed for {row['allele']}")
            if row["representative_pseudo_record_id"] not in record_ids:
                fail(f"representative pseudo changed for {row['allele']}")

        missing = sorted(set(expected) - set(seen))
        extra = sorted(set(seen) - set(expected))
        if missing or extra:
            detail = (
                f"component coverage changed (missing={len(missing)}, "
                f"extra={len(extra)})"
            )
            fail(detail)

        counts = self.metadata.get("counts", {})
        if isinstance(counts, Mapping):
            declared_components = counts.get("pseudo_pbr_projection_components")
            declared_pseudos = counts.get("pseudo_pbr_projection_pseudos")
            if declared_components is not None and int(declared_components) != len(seen):
                fail("metadata projection component count changed")
            if declared_pseudos is not None and int(declared_pseudos) != len(self._pseudos):
                fail("metadata projection pseudo count changed")

    @staticmethod
    def _load_sequences(path: object) -> Tuple[SequenceCandidate, ...]:
        candidates = []
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                candidates.append(
                    SequenceCandidate(
                        record_id=row["record_id"],
                        allele=row["allele"],
                        aliases=_split_aliases(row["aliases"]),
                        species=row["species"],
                        mhc_class=row["mhc_class"],
                        chain_role=row["chain_role"],
                        locus=row["locus"],
                        sequence_kind=row["sequence_kind"],
                        sequence=row["sequence"],
                        source=row["source"],
                        source_version=row["source_version"],
                        source_accession=row["source_accession"],
                        is_complete=row["is_complete"].lower() == "true",
                    )
                )
        return tuple(candidates)

    @staticmethod
    def _load_pseudos(path: object) -> Tuple[PseudoCandidate, ...]:
        candidates = []
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                candidates.append(
                    PseudoCandidate(
                        record_id=row["record_id"],
                        canonical_name=row["canonical_name"],
                        aliases=_split_aliases(row["aliases"]),
                        species=row["species"],
                        mhc_class=row["mhc_class"],
                        definition_id=row["definition_id"],
                        pseudo_sequence=row["pseudo_sequence"],
                        alpha_allele=row["alpha_allele"],
                        beta_allele=row["beta_allele"],
                        source=row["source"],
                        source_version=row["source_version"],
                        evidence_level=row.get("evidence_level", "legacy_secondary"),
                        official_support=row.get("official_support", "false").lower()
                        == "true",
                        validation_status=row.get(
                            "validation_status", "legacy_secondary_only"
                        ),
                        source_record_ids=_split_aliases(
                            row.get("source_record_ids", "")
                        ),
                        quality_flags=_split_aliases(row.get("quality_flags", "")),
                    )
                )
        return tuple(candidates)

    @staticmethod
    def _load_af3_sequences(path: object) -> Tuple[AF3SequenceCandidate, ...]:
        candidates = []
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                span_match = re.fullmatch(r"([0-9]+)-([0-9]+)", row["precursor_span_1based"])
                if span_match is None:
                    raise MHCTransformError(
                        f"invalid AF3 precursor span: {row['precursor_span_1based']}"
                    )
                candidates.append(
                    AF3SequenceCandidate(
                        record_id=row["record_id"],
                        allele=row["allele"],
                        aliases=_split_aliases(row["aliases"]),
                        species=row["species"],
                        mhc_class=row["mhc_class"],
                        chain_role=row["chain_role"],
                        locus=row["locus"],
                        definition_id=row["definition_id"],
                        sequence=row["sequence"],
                        source_full_record_id=row["source_full_record_id"],
                        precursor_start_1based=int(span_match.group(1)),
                        precursor_end_1based=int(span_match.group(2)),
                        boundary_method=row["boundary_method"],
                        boundary_source=row["boundary_source"],
                        source_accession=row["source_accession"],
                    )
                )
        return tuple(candidates)

    @staticmethod
    def _load_pbr_sequences(path: object) -> Tuple[PBRSequenceCandidate, ...]:
        candidates = []
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                af3_span = re.fullmatch(r"([0-9]+)-([0-9]+)", row["af3_span_1based"])
                precursor_span = re.fullmatch(
                    r"([0-9]+)-([0-9]+)", row["precursor_span_1based"]
                )
                if af3_span is None or precursor_span is None:
                    raise MHCTransformError(
                        f"invalid PBR span for {row['record_id']}"
                    )
                candidates.append(
                    PBRSequenceCandidate(
                        record_id=row["record_id"],
                        allele=row["allele"],
                        aliases=_split_aliases(row["aliases"]),
                        species=row["species"],
                        mhc_class=row["mhc_class"],
                        chain_role=row["chain_role"],
                        locus=row["locus"],
                        definition_id=row["definition_id"],
                        sequence=row["sequence"],
                        source_af3_record_id=row["source_af3_record_id"],
                        af3_start_1based=int(af3_span.group(1)),
                        af3_end_1based=int(af3_span.group(2)),
                        source_full_record_id=row["source_full_record_id"],
                        precursor_start_1based=int(precursor_span.group(1)),
                        precursor_end_1based=int(precursor_span.group(2)),
                        boundary_method=row["boundary_method"],
                        boundary_source=row["boundary_source"],
                        source_accession=row["source_accession"],
                    )
                )
        return tuple(candidates)

    @staticmethod
    def _load_pairs(path: object) -> Tuple[Mapping[str, str], ...]:
        with path.open("r", encoding="utf-8", newline="") as handle:
            return tuple(dict(row) for row in csv.DictReader(handle))

    @staticmethod
    def _load_fallback_sequences(
        path: object, *, af3: bool
    ) -> Tuple[FallbackSequenceCandidate | FallbackAF3SequenceCandidate, ...]:
        candidates = []
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                common = {
                    "record_id": row["record_id"],
                    "allele": row["allele"],
                    "aliases": _split_aliases(row["aliases"]),
                    "species": row["species"],
                    "mhc_class": row["mhc_class"],
                    "chain_role": row["chain_role"],
                    "locus": row["locus"],
                    "sequence": row["sequence"],
                    "evidence_level": row["evidence_level"],
                    "observed_spans_1based": row["observed_spans_1based"],
                    "imputed_spans_1based": row["imputed_spans_1based"],
                    "observed_residue_count": int(row["observed_residue_count"]),
                    "imputed_residue_count": int(row["imputed_residue_count"]),
                    "source_record_ids": _split_aliases(row["source_record_ids"]),
                    "template_record_id": row["template_record_id"],
                    "template_allele": row["template_allele"],
                    "target_pseudo_record_id": row["target_pseudo_record_id"],
                    "pseudo_distance": int(row["pseudo_distance"]),
                    "unknown_pseudo_positions": tuple(
                        int(value)
                        for value in row["unknown_pseudo_positions"].split(";")
                        if value
                    ),
                    "quality_flags": _split_aliases(row["quality_flags"]),
                }
                if af3:
                    candidates.append(
                        FallbackAF3SequenceCandidate(
                            **common, definition_id=row["representation"]
                        )
                    )
                else:
                    candidates.append(
                        FallbackSequenceCandidate(
                            **common, representation=row["representation"]
                        )
                    )
        return tuple(candidates)

    @staticmethod
    def _load_fallback_pbr_sequences(
        path: object,
    ) -> Tuple[FallbackPBRSequenceCandidate, ...]:
        candidates = []
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                af3_span = re.fullmatch(r"([0-9]+)-([0-9]+)", row["af3_span_1based"])
                if af3_span is None:
                    raise MHCTransformError(
                        f"invalid fallback PBR span: {row['af3_span_1based']}"
                    )
                candidates.append(
                    FallbackPBRSequenceCandidate(
                        record_id=row["record_id"],
                        allele=row["allele"],
                        aliases=_split_aliases(row["aliases"]),
                        species=row["species"],
                        mhc_class=row["mhc_class"],
                        chain_role=row["chain_role"],
                        locus=row["locus"],
                        definition_id=row["representation"],
                        sequence=row["sequence"],
                        source_af3_record_id=row["source_af3_record_id"],
                        af3_start_1based=int(af3_span.group(1)),
                        af3_end_1based=int(af3_span.group(2)),
                        evidence_level=row["evidence_level"],
                        observed_spans_1based=row["observed_spans_1based"],
                        imputed_spans_1based=row["imputed_spans_1based"],
                        observed_residue_count=int(row["observed_residue_count"]),
                        imputed_residue_count=int(row["imputed_residue_count"]),
                        source_record_ids=_split_aliases(row["source_record_ids"]),
                        template_record_id=row["template_record_id"],
                        template_allele=row["template_allele"],
                        target_pseudo_record_id=row["target_pseudo_record_id"],
                        pseudo_distance=int(row["pseudo_distance"]),
                        unknown_pseudo_positions=tuple(
                            int(value)
                            for value in row["unknown_pseudo_positions"].split(";")
                            if value
                        ),
                        quality_flags=_split_aliases(row["quality_flags"]),
                    )
                )
        return tuple(candidates)

    @staticmethod
    def _filter_sequence_candidates(
        candidates: Iterable[SequenceCandidate],
        *,
        mhc_class: str | None,
        species: str | None,
        include_partial: bool,
    ) -> Tuple[SequenceCandidate, ...]:
        expected_class = normalize_mhc_class(mhc_class) if mhc_class else None
        expected_species = species.lower().strip() if species else None
        if expected_species and expected_species not in {"human", "mouse"}:
            raise MHCTransformError("species must be human or mouse")
        filtered = (
            candidate
            for candidate in candidates
            if (expected_class is None or candidate.mhc_class == expected_class)
            and (expected_species is None or candidate.species == expected_species)
            and (include_partial or candidate.is_complete)
        )
        return _deduplicate(filtered)  # type: ignore[return-value]

    @staticmethod
    def _filter_af3_candidates(
        candidates: Iterable[AF3SequenceCandidate],
        *,
        mhc_class: str | None,
        species: str | None,
        definition_id: str,
    ) -> Tuple[AF3SequenceCandidate, ...]:
        if definition_id != AF3_DEFINITION:
            raise MHCTransformError(
                f"unsupported AF3 definition_id: {definition_id}; expected {AF3_DEFINITION}"
            )
        expected_class = normalize_mhc_class(mhc_class) if mhc_class else None
        expected_species = species.lower().strip() if species else None
        if expected_species and expected_species not in {"human", "mouse"}:
            raise MHCTransformError("species must be human or mouse")
        filtered = (
            candidate
            for candidate in candidates
            if candidate.definition_id == definition_id
            and (expected_class is None or candidate.mhc_class == expected_class)
            and (expected_species is None or candidate.species == expected_species)
        )
        return _deduplicate(filtered)  # type: ignore[return-value]

    @staticmethod
    def _filter_pbr_candidates(
        candidates: Iterable[PBRSequenceCandidate],
        *,
        mhc_class: str | None,
        species: str | None,
        definition_id: str | None,
    ) -> Tuple[PBRSequenceCandidate, ...]:
        if definition_id is not None and definition_id not in PBR_DEFINITIONS:
            raise MHCTransformError(
                "unsupported PBR definition_id: "
                f"{definition_id}; expected one of {', '.join(PBR_DEFINITIONS)}"
            )
        expected_class = normalize_mhc_class(mhc_class) if mhc_class else None
        expected_species = species.lower().strip() if species else None
        if expected_species and expected_species not in {"human", "mouse"}:
            raise MHCTransformError("species must be human or mouse")
        filtered = (
            candidate
            for candidate in candidates
            if (definition_id is None or candidate.definition_id == definition_id)
            and (expected_class is None or candidate.mhc_class == expected_class)
            and (expected_species is None or candidate.species == expected_species)
        )
        return _deduplicate(filtered)  # type: ignore[return-value]

    @staticmethod
    def _filter_fallback_candidates(
        candidates: Iterable[object],
        *,
        mhc_class: str | None,
        species: str | None,
        fallback_policy: str,
        allow_expression_qualified: bool = False,
        prefer_best: bool = True,
    ) -> Tuple[object, ...]:
        policy = normalize_fallback_policy(fallback_policy)
        expected_class = normalize_mhc_class(mhc_class) if mhc_class else None
        expected_species = species.lower().strip() if species else None
        if expected_species and expected_species not in {"human", "mouse"}:
            raise MHCTransformError("species must be human or mouse")
        eligible = [
            candidate
            for candidate in candidates
            if evidence_allowed(policy, getattr(candidate, "evidence_level"))
            and (expected_class is None or getattr(candidate, "mhc_class") == expected_class)
            and (expected_species is None or getattr(candidate, "species") == expected_species)
            and (
                allow_expression_qualified
                or "expression_qualified_target"
                not in getattr(candidate, "quality_flags")
            )
        ]
        if not eligible:
            return ()
        if not prefer_best:
            return _deduplicate(eligible)
        best_rank = min(
            evidence_priority(policy, getattr(candidate, "evidence_level"))
            for candidate in eligible
        )
        return _deduplicate(
            candidate
            for candidate in eligible
            if evidence_priority(policy, getattr(candidate, "evidence_level"))
            == best_rank
        )

    def _pseudo_pbr_certificates(
        self,
        name_key: str,
        *,
        mhc_class: str | None,
        species: str | None,
        definition_id: str | None,
    ) -> Tuple[Mapping[str, str], ...]:
        """Return build-certified pseudo/PBR projections for one component."""

        expected_class = normalize_mhc_class(mhc_class) if mhc_class else None
        expected_species = species.lower().strip() if species else None
        if expected_species and expected_species not in {"human", "mouse"}:
            raise MHCTransformError("species must be human or mouse")
        unique: Dict[str, Mapping[str, str]] = {}
        for row in self._pseudo_pbr_projection_name_index.get(name_key, ()):
            if expected_class is not None and row["mhc_class"] != expected_class:
                continue
            if expected_species is not None and row["species"] != expected_species:
                continue
            if definition_id is not None and row["definition_id"] != definition_id:
                continue
            unique[row["projection_id"]] = row
        return tuple(unique[key] for key in sorted(unique))

    def allele_to_sequence(
        self,
        allele: str,
        *,
        mhc_class: str | None = None,
        species: str | None = None,
        include_partial: bool = False,
        allow_two_field_projection: bool = False,
        fallback: str = "none",
    ) -> QueryResult:
        """Convert a chain name to strict precursor or explicit fallback sequence(s)."""

        fallback_policy = normalize_fallback_policy(fallback)
        if include_partial and fallback_policy != "none":
            raise MHCTransformError(
                "include_partial and fallback are separate modes; use only one"
            )
        exact = self._sequence_name_index.get(exact_name_key(allele), ())
        notes: List[str] = []
        projected = False
        fallback_key = exact_name_key(allele)
        candidates = self._filter_sequence_candidates(
            exact,
            mhc_class=mhc_class,
            species=species,
            include_partial=include_partial,
        )
        if not candidates and not exact and allow_two_field_projection:
            two_field = two_field_hla(allele)
            if two_field and exact_name_key(two_field) != exact_name_key(allele):
                candidates = self._filter_sequence_candidates(
                    self._sequence_name_index.get(exact_name_key(two_field), ()),
                    mhc_class=mhc_class,
                    species=species,
                    include_partial=include_partial,
                )
                if candidates:
                    projected = True
                    fallback_key = exact_name_key(two_field)
                    notes.append(
                        "No exact deposited name was found; explicit two-field projection was enabled, so all matching protein-group candidates are returned."
                    )
        if not candidates and not include_partial:
            partial = self._filter_sequence_candidates(
                exact,
                mhc_class=mhc_class,
                species=species,
                include_partial=True,
            )
            if partial:
                notes.append(
                    "Matching deposited sequence(s) exist but fail the precursor-completeness rule; retry with include_partial=True."
                )
        fallback_candidates: Tuple[object, ...] = ()
        if not candidates and fallback_policy != "none":
            fallback_candidates = self._filter_fallback_candidates(
                self._fallback_sequence_name_index.get(fallback_key, ()),
                mhc_class=mhc_class,
                species=species,
                fallback_policy=fallback_policy,
                allow_expression_qualified=True,
            )
            if fallback_candidates:
                notes.append(
                    "Strict complete sequence unavailable; an explicitly labelled fallback representation was returned."
                )
        returned_candidates = fallback_candidates or candidates
        evidence_levels = {
            _candidate_evidence(candidate)
            for candidate in returned_candidates
        }
        evidence = (
            next(iter(evidence_levels))
            if len(evidence_levels) == 1
            else ("mixed" if evidence_levels else "source_exact")
        )
        representations = {
            getattr(candidate, "sequence_kind", "precursor_full")
            for candidate in candidates
        }
        representation = (
            next(iter(representations))
            if len(representations) == 1
            else ("mixed_deposited_sequence" if representations else "precursor_full")
        )
        if fallback_candidates:
            representations = {
                getattr(candidate, "representation") for candidate in fallback_candidates
            }
            representation = (
                next(iter(representations))
                if len(representations) == 1
                else "mixed_fallback_sequence"
            )
        return QueryResult(
            query=allele,
            representation=representation,
            status=result_status(len(returned_candidates), projected=projected),
            candidates=returned_candidates,
            notes=tuple(notes),
            fallback_used=fallback_policy if fallback_candidates else "none",
            evidence=evidence,
        )

    def sequence_to_alleles(
        self,
        sequence: str,
        *,
        mhc_class: str | None = None,
        species: str | None = None,
        include_partial: bool = False,
        include_fallback: str = "none",
    ) -> QueryResult:
        """Return strict exact and opt-in compatible names for a protein sequence."""

        fallback_policy = normalize_fallback_policy(include_fallback)
        if include_partial and fallback_policy != "none":
            raise MHCTransformError(
                "include_partial and include_fallback are separate modes; use only one"
            )
        normalized = normalize_sequence(
            sequence,
            allow_x=include_partial or fallback_policy != "none",
        )
        candidates = self._filter_sequence_candidates(
            self._sequence_value_index.get(normalized, ()),
            mhc_class=mhc_class,
            species=species,
            include_partial=include_partial,
        )
        fallback_candidates: Tuple[object, ...] = ()
        if fallback_policy != "none":
            fallback_candidates = self._filter_fallback_candidates(
                self._fallback_sequence_value_index.get(normalized, ()),
                mhc_class=mhc_class,
                species=species,
                fallback_policy=fallback_policy,
                allow_expression_qualified=True,
                prefer_best=False,
            )
        returned_candidates = _deduplicate((*candidates, *fallback_candidates))
        deposited_names = {
            alias
            for candidate in returned_candidates
            for alias in candidate.aliases
            if alias.startswith("HLA-") and "*" in alias
        }
        reverse_status = result_status(len(returned_candidates))
        notes: Tuple[str, ...] = ()
        if len(deposited_names) > 1 and not fallback_candidates:
            reverse_status = "exact_ambiguous"
            notes = (
                "Multiple deposited HLA allele names share this exact protein; they are compactly listed in candidate aliases.",
            )
        if fallback_candidates:
            reverse_status = (
                "compatible_unique"
                if len(returned_candidates) == 1
                else "compatible_ambiguous"
            )
            notes = (*notes, "Sequence matches an engineered/partial fallback construct, not proof of an exact biological full allele.")
        evidence_levels = {
            _candidate_evidence(candidate)
            for candidate in returned_candidates
        }
        deposited_representations = {
            candidate.sequence_kind for candidate in candidates
        }
        strict_representation = (
            next(iter(deposited_representations))
            if len(deposited_representations) == 1
            else (
                "mixed_deposited_sequence"
                if deposited_representations
                else "precursor_full"
            )
        )
        return QueryResult(
            query=normalized,
            representation=(
                "fallback_sequence_match"
                if fallback_candidates
                else strict_representation
            ),
            status=reverse_status,
            candidates=returned_candidates,
            notes=notes,
            fallback_used=fallback_policy if fallback_candidates else "none",
            evidence=(
                next(iter(evidence_levels))
                if len(evidence_levels) == 1
                else ("mixed" if evidence_levels else "source_exact")
            ),
        )

    def allele_to_af3(
        self,
        allele: str,
        *,
        mhc_class: str | None = None,
        species: str | None = None,
        definition_id: str = AF3_DEFINITION,
        allow_two_field_projection: bool = False,
        fallback: str = "none",
        allow_expression_qualified: bool = False,
    ) -> QueryResult:
        """Convert a chain allele/name to a mature soluble AF3 input chain."""

        fallback_policy = normalize_fallback_policy(fallback)
        query_key = exact_name_key(allele)
        fallback_key = query_key
        exact = self._af3_name_index.get(query_key, ())
        candidates = self._filter_af3_candidates(
            exact,
            mhc_class=mhc_class,
            species=species,
            definition_id=definition_id,
        )
        notes: List[str] = []
        projected = False
        deposited_full = self._sequence_name_index.get(query_key, ())
        exclusions = self._af3_exclusion_index.get(query_key, ())
        if (
            not candidates
            and not exact
            and not deposited_full
            and not exclusions
            and allow_two_field_projection
        ):
            two_field = two_field_hla(allele)
            if two_field and exact_name_key(two_field) != query_key:
                candidates = self._filter_af3_candidates(
                    self._af3_name_index.get(exact_name_key(two_field), ()),
                    mhc_class=mhc_class,
                    species=species,
                    definition_id=definition_id,
                )
                if candidates:
                    projected = True
                    fallback_key = exact_name_key(two_field)
                    notes.append(
                        "No exact deposited name was found; explicit two-field projection was enabled, so all matching AF3 chain candidates are returned."
                    )
        if not candidates:
            reasons = sorted({item["reason"] for item in exclusions})
            if reasons:
                notes.extend(f"Not AF3-ready: {reason}" for reason in reasons)
            elif deposited_full:
                notes.append(
                    "An exact deposited protein name exists but is partial or otherwise not AF3-ready; two-field projection is suppressed."
                )
        fallback_candidates: Tuple[object, ...] = ()
        if not candidates and fallback_policy != "none":
            fallback_candidates = self._filter_fallback_candidates(
                self._fallback_af3_name_index.get(fallback_key, ()),
                mhc_class=mhc_class,
                species=species,
                fallback_policy=fallback_policy,
                allow_expression_qualified=allow_expression_qualified,
            )
            if fallback_candidates:
                notes.append(
                    "Strict AF3 chain unavailable; an explicitly labelled ectodomain fallback/construct was returned."
                )
            elif not allow_expression_qualified:
                unfiltered = self._filter_fallback_candidates(
                    self._fallback_af3_name_index.get(fallback_key, ()),
                    mhc_class=mhc_class,
                    species=species,
                    fallback_policy=fallback_policy,
                    allow_expression_qualified=True,
                )
                if unfiltered:
                    notes.append(
                        "Fallback exists but the allele has an expression suffix; set allow_expression_qualified=True to model it explicitly."
                    )
        returned_candidates = fallback_candidates or candidates
        evidence_levels = {
            getattr(candidate, "evidence_level", "source_exact")
            for candidate in returned_candidates
        }
        representations = {
            getattr(candidate, "definition_id", definition_id)
            for candidate in returned_candidates
        }
        return QueryResult(
            query=allele,
            representation=(
                next(iter(representations))
                if len(representations) == 1
                else definition_id
            ),
            status=result_status(len(returned_candidates), projected=projected),
            candidates=returned_candidates,
            notes=tuple(notes),
            fallback_used=fallback_policy if fallback_candidates else "none",
            evidence=(
                next(iter(evidence_levels))
                if len(evidence_levels) == 1
                else ("mixed" if evidence_levels else "source_exact")
            ),
        )

    def af3_to_alleles(
        self,
        sequence: str,
        *,
        mhc_class: str | None = None,
        species: str | None = None,
        definition_id: str = AF3_DEFINITION,
        include_fallback: str = "none",
    ) -> QueryResult:
        """Return strict exact and opt-in compatible names for an AF3 chain."""

        fallback_policy = normalize_fallback_policy(include_fallback)
        normalized = normalize_sequence(sequence)
        candidates = self._filter_af3_candidates(
            self._af3_value_index.get(normalized, ()),
            mhc_class=mhc_class,
            species=species,
            definition_id=definition_id,
        )
        fallback_candidates: Tuple[object, ...] = ()
        if fallback_policy != "none":
            fallback_candidates = self._filter_fallback_candidates(
                self._fallback_af3_value_index.get(normalized, ()),
                mhc_class=mhc_class,
                species=species,
                fallback_policy=fallback_policy,
                allow_expression_qualified=True,
                prefer_best=False,
            )
        returned_candidates = _deduplicate((*candidates, *fallback_candidates))
        deposited_names = {
            alias
            for candidate in returned_candidates
            for alias in candidate.aliases
            if alias.startswith("HLA-") and "*" in alias
        }
        reverse_status = result_status(len(returned_candidates))
        notes: Tuple[str, ...] = ()
        if len(deposited_names) > 1 and not fallback_candidates:
            reverse_status = "exact_ambiguous"
            notes = (
                "Multiple deposited HLA allele names share this exact AF3 chain; they are retained in candidate aliases.",
            )
        if fallback_candidates:
            reverse_status = (
                "compatible_unique"
                if len(returned_candidates) == 1
                else "compatible_ambiguous"
            )
            notes = (*notes, "Sequence matches an engineered AF3 fallback construct; this is compatibility, not proof of an exact biological allele.")
        evidence_levels = {
            getattr(candidate, "evidence_level", "source_exact")
            for candidate in returned_candidates
        }
        return QueryResult(
            query=normalized,
            representation=(
                "fallback_af3_construct_match" if fallback_candidates else definition_id
            ),
            status=reverse_status,
            candidates=returned_candidates,
            notes=notes,
            fallback_used=fallback_policy if fallback_candidates else "none",
            evidence=(
                next(iter(evidence_levels))
                if len(evidence_levels) == 1
                else ("mixed" if evidence_levels else "source_exact")
            ),
        )

    def allele_to_pbr(
        self,
        allele: str,
        *,
        mhc_class: str | None = None,
        species: str | None = None,
        definition_id: str | None = None,
        allow_two_field_projection: bool = False,
        fallback: str = "none",
        allow_expression_qualified: bool = False,
        require_pseudo_match: bool = False,
    ) -> QueryResult:
        """Convert a chain allele/name to a contiguous peptide-binding region.

        ``require_pseudo_match`` is an explicit compatibility mode.  It returns
        only the strict or fallback PBR record selected by the build-time
        pseudo/PBR projection contract.  The flag itself authorizes a certified
        synthetic fallback when the strict biological PBR conflicts with the
        DTU pseudo; normal calls keep the historical strict-first behaviour.
        """

        fallback_policy = normalize_fallback_policy(fallback)
        query_key = exact_name_key(allele)
        fallback_key = query_key
        exact = self._pbr_name_index.get(query_key, ())
        candidates = self._filter_pbr_candidates(
            exact,
            mhc_class=mhc_class,
            species=species,
            definition_id=definition_id,
        )
        notes: List[str] = []
        projected = False
        deposited_full = self._sequence_name_index.get(query_key, ())
        exclusions = self._af3_exclusion_index.get(query_key, ())
        if (
            not candidates
            and not exact
            and not deposited_full
            and not exclusions
            and allow_two_field_projection
        ):
            two_field = two_field_hla(allele)
            if two_field and exact_name_key(two_field) != query_key:
                candidates = self._filter_pbr_candidates(
                    self._pbr_name_index.get(exact_name_key(two_field), ()),
                    mhc_class=mhc_class,
                    species=species,
                    definition_id=definition_id,
                )
                if candidates:
                    projected = True
                    fallback_key = exact_name_key(two_field)
                    notes.append(
                        "No exact deposited name was found; explicit two-field projection was enabled, so all matching PBR candidates are returned."
                    )
        if not candidates:
            reasons = sorted({item["reason"] for item in exclusions})
            if reasons:
                notes.extend(f"Not strict-PBR-ready: {reason}" for reason in reasons)
            elif deposited_full:
                notes.append(
                    "An exact deposited protein name exists but is partial or otherwise lacks a strict PBR; two-field projection is suppressed."
                )

        fallback_candidates: Tuple[object, ...] = ()
        contract_policy = "pseudo-scaffold"
        if require_pseudo_match:
            certificates = self._pseudo_pbr_certificates(
                fallback_key,
                mhc_class=mhc_class,
                species=species,
                definition_id=definition_id,
            )
            strict_ids = {
                row["pbr_record_id"]
                for row in certificates
                if row["pbr_source"] == "strict"
            }
            fallback_ids = {
                row["pbr_record_id"]
                for row in certificates
                if row["pbr_source"] == "fallback"
            }
            candidates = tuple(
                candidate
                for candidate in candidates
                if candidate.record_id in strict_ids
            )
            certified_fallback_pool = tuple(
                candidate
                for candidate in self._fallback_pbr_name_index.get(fallback_key, ())
                if candidate.record_id in fallback_ids
            )
            fallback_candidates = self._filter_fallback_candidates(
                certified_fallback_pool,
                mhc_class=mhc_class,
                species=species,
                fallback_policy=contract_policy,
                allow_expression_qualified=allow_expression_qualified,
            )
            if definition_id is not None:
                fallback_candidates = tuple(
                    candidate
                    for candidate in fallback_candidates
                    if getattr(candidate, "definition_id") == definition_id
                )
            if candidates:
                notes.append(
                    "Pseudo-match contract: returned a certified strict PBR; every non-X pseudo residue matches its projected PBR coordinate exactly, while X remains an explicit wildcard."
                )
            elif fallback_candidates:
                synthetic = any(
                    getattr(candidate, "evidence_level", "")
                    in {IMPUTED_EVIDENCE, PSEUDO_SCAFFOLD_EVIDENCE}
                    for candidate in fallback_candidates
                )
                notes.append(
                    "Pseudo-match contract: returned the certified fallback PBR selected for this component; this is compatibility provenance, not proof of a strict deposited biological PBR."
                )
                if synthetic:
                    notes.append(
                        "The certified PBR is a synthetic engineered construct constrained by the pseudo; unknown X residues were treated as wildcards and were not inserted into the protein chain."
                    )
            else:
                expression_blocked = False
                if certified_fallback_pool and not allow_expression_qualified:
                    expression_blocked = bool(
                        self._filter_fallback_candidates(
                            certified_fallback_pool,
                            mhc_class=mhc_class,
                            species=species,
                            fallback_policy=contract_policy,
                            allow_expression_qualified=True,
                        )
                    )
                if expression_blocked:
                    notes.append(
                        "A pseudo-match-certified PBR fallback exists, but the allele has an expression suffix; set allow_expression_qualified=True to model it explicitly."
                    )
                elif certificates:
                    notes.append(
                        "A pseudo/PBR projection certificate exists, but no certified record satisfies the requested filters."
                    )
                else:
                    notes.append(
                        "No pseudo/PBR projection certificate exists for this exact component name."
                    )
        elif not candidates and fallback_policy != "none":
            fallback_candidates = self._filter_fallback_candidates(
                self._fallback_pbr_name_index.get(fallback_key, ()),
                mhc_class=mhc_class,
                species=species,
                fallback_policy=fallback_policy,
                allow_expression_qualified=allow_expression_qualified,
            )
            if definition_id is not None:
                fallback_candidates = tuple(
                    candidate
                    for candidate in fallback_candidates
                    if getattr(candidate, "definition_id") == definition_id
                )
            if fallback_candidates:
                notes.append(
                    "Strict PBR unavailable; an explicitly labelled PBR fallback/construct was returned."
                )
            elif not allow_expression_qualified:
                unfiltered = self._filter_fallback_candidates(
                    self._fallback_pbr_name_index.get(fallback_key, ()),
                    mhc_class=mhc_class,
                    species=species,
                    fallback_policy=fallback_policy,
                    allow_expression_qualified=True,
                )
                if definition_id is not None:
                    unfiltered = tuple(
                        candidate
                        for candidate in unfiltered
                        if getattr(candidate, "definition_id") == definition_id
                    )
                if unfiltered:
                    notes.append(
                        "PBR fallback exists but the allele has an expression suffix; set allow_expression_qualified=True to model it explicitly."
                    )
        returned_candidates = fallback_candidates or candidates
        evidence_levels = {
            getattr(candidate, "evidence_level", "source_exact")
            for candidate in returned_candidates
        }
        representations = {
            getattr(candidate, "definition_id") for candidate in returned_candidates
        }
        if require_pseudo_match:
            synthetic_contract = bool(
                fallback_candidates
                and evidence_levels & {IMPUTED_EVIDENCE, PSEUDO_SCAFFOLD_EVIDENCE}
            )
            representation = (
                "pseudo_compatible_pbr"
                if not returned_candidates
                else (
                    "pseudo_compatible_synthetic_pbr"
                    if synthetic_contract
                    else (
                        "pseudo_compatible_fallback_pbr"
                        if fallback_candidates
                        else "pseudo_compatible_strict_pbr"
                    )
                )
            )
        else:
            representation = (
                next(iter(representations))
                if len(representations) == 1
                else "mixed_peptide_binding_region"
            )
        return QueryResult(
            query=allele,
            representation=representation,
            status=result_status(len(returned_candidates), projected=projected),
            candidates=returned_candidates,
            notes=tuple(notes),
            fallback_used=(
                contract_policy
                if require_pseudo_match and fallback_candidates
                else (fallback_policy if fallback_candidates else "none")
            ),
            evidence=(
                next(iter(evidence_levels))
                if len(evidence_levels) == 1
                else ("mixed" if evidence_levels else "source_exact")
            ),
        )

    def pbr_to_alleles(
        self,
        sequence: str,
        *,
        mhc_class: str | None = None,
        species: str | None = None,
        definition_id: str | None = None,
        include_fallback: str = "none",
    ) -> QueryResult:
        """Return strict exact and opt-in compatible names for a PBR chain."""

        fallback_policy = normalize_fallback_policy(include_fallback)
        normalized = normalize_sequence(sequence)
        candidates = self._filter_pbr_candidates(
            self._pbr_value_index.get(normalized, ()),
            mhc_class=mhc_class,
            species=species,
            definition_id=definition_id,
        )
        fallback_candidates: Tuple[object, ...] = ()
        if fallback_policy != "none":
            fallback_candidates = self._filter_fallback_candidates(
                self._fallback_pbr_value_index.get(normalized, ()),
                mhc_class=mhc_class,
                species=species,
                fallback_policy=fallback_policy,
                allow_expression_qualified=True,
                prefer_best=False,
            )
            if definition_id is not None:
                fallback_candidates = tuple(
                    candidate
                    for candidate in fallback_candidates
                    if getattr(candidate, "definition_id") == definition_id
                )
        returned_candidates = _deduplicate((*candidates, *fallback_candidates))
        reverse_status = result_status(len(returned_candidates))
        notes: Tuple[str, ...] = ()
        if fallback_candidates:
            reverse_status = (
                "compatible_unique"
                if len(returned_candidates) == 1
                else "compatible_ambiguous"
            )
            notes = (
                "Sequence matches an engineered/partial PBR fallback construct; this is compatibility, not proof of an exact biological allele.",
            )
        evidence_levels = {
            getattr(candidate, "evidence_level", "source_exact")
            for candidate in returned_candidates
        }
        representations = {
            getattr(candidate, "definition_id") for candidate in returned_candidates
        }
        return QueryResult(
            query=normalized,
            representation=(
                "fallback_pbr_construct_match"
                if fallback_candidates
                else (
                    next(iter(representations))
                    if len(representations) == 1
                    else "mixed_peptide_binding_region"
                )
            ),
            status=reverse_status,
            candidates=returned_candidates,
            notes=notes,
            fallback_used=fallback_policy if fallback_candidates else "none",
            evidence=(
                next(iter(evidence_levels))
                if len(evidence_levels) == 1
                else ("mixed" if evidence_levels else "source_exact")
            ),
        )

    def allele_to_pseudo(
        self,
        allele_or_pair: str,
        *,
        mhc_class: str,
        definition_id: str,
        species: str | None = None,
    ) -> QueryResult:
        """Convert an allele/molecule name to exact 34-aa pseudosequence(s)."""

        expected_class = normalize_mhc_class(mhc_class)
        if not definition_id:
            raise MHCTransformError("definition_id is required for pseudosequence lookup")
        expected_species = species.lower().strip() if species else None
        if expected_species and expected_species not in {"human", "mouse"}:
            raise MHCTransformError("species must be human or mouse")

        query_key = exact_name_key(allele_or_pair)
        indexed = self._pseudo_name_index.get(query_key, ())
        notes: List[str] = []
        projected = False
        if not indexed:
            components = [
                part for part in re.split(r"\s*(?:/|\+)\s*", allele_or_pair) if part
            ]
            if len(components) == 2:
                fallback_name = "/".join(two_field_hla(part) or part for part in components)
            else:
                fallback_name = two_field_hla(allele_or_pair) or allele_or_pair
            if exact_name_key(fallback_name) != query_key:
                deposited = all(
                    exact_name_key(component)
                    in self._sequence_name_index
                    for component in components
                )
                if deposited:
                    indexed = self._pseudo_name_index.get(
                        exact_name_key(fallback_name), ()
                    )
                    if indexed:
                        projected = True
                        notes.append(
                            "The verified deposited higher-field name was projected to its two-field protein group before pseudosequence lookup."
                        )

        candidates = []
        for candidate in indexed:
            if candidate.mhc_class != expected_class:
                continue
            if candidate.definition_id != definition_id:
                continue
            if expected_species and candidate.species != expected_species:
                continue
            candidates.append(candidate)
        candidates = _deduplicate(candidates)
        evidence_levels = {
            candidate.evidence_level for candidate in candidates
        }
        return QueryResult(
            query=allele_or_pair,
            representation=definition_id,
            status=result_status(len(candidates), projected=projected),
            candidates=candidates,
            notes=tuple(notes),
            evidence=(
                next(iter(evidence_levels))
                if len(evidence_levels) == 1
                else ("mixed" if evidence_levels else "source_exact")
            ),
        )

    def pseudo_to_alleles(
        self,
        pseudo_sequence: str,
        *,
        mhc_class: str,
        definition_id: str,
        species: str | None = None,
    ) -> QueryResult:
        """Return every exact allele/molecule name sharing a 34-aa descriptor."""

        expected_class = normalize_mhc_class(mhc_class)
        if not definition_id:
            raise MHCTransformError("definition_id is required for pseudosequence lookup")
        normalized = normalize_sequence(pseudo_sequence, allow_x=True)
        if len(normalized) != 34:
            raise MHCTransformError("a 34-aa pseudosequence must contain exactly 34 residues")
        expected_species = species.lower().strip() if species else None
        if expected_species and expected_species not in {"human", "mouse"}:
            raise MHCTransformError("species must be human or mouse")
        candidates = []
        for candidate in self._pseudo_value_index.get(normalized, ()):
            if candidate.mhc_class != expected_class:
                continue
            if candidate.definition_id != definition_id:
                continue
            if expected_species and candidate.species != expected_species:
                continue
            candidates.append(candidate)
        candidates = _deduplicate(candidates)
        evidence_levels = {
            candidate.evidence_level for candidate in candidates
        }
        return QueryResult(
            query=normalized,
            representation=definition_id,
            status=result_status(len(candidates)),
            candidates=candidates,
            evidence=(
                next(iter(evidence_levels))
                if len(evidence_levels) == 1
                else ("mixed" if evidence_levels else "source_exact")
            ),
        )

    def _pair_specs(self, name: str, mhc_class: str) -> Tuple[Tuple[str, str, str], ...]:
        key = exact_name_key(name)
        specs = []
        for pair in self._pair_index.get(key, ()):
            if pair["mhc_class"] == mhc_class:
                specs.append(
                    (
                        pair["canonical_name"],
                        pair["alpha_or_heavy"],
                        pair["beta_or_companion"],
                    )
                )
        if specs:
            return tuple(specs)

        if mhc_class == "II":
            matching_pseudos = [
                candidate
                for candidate in self._pseudo_name_index.get(key, ())
                if candidate.mhc_class == "II"
                and candidate.definition_id == PSEUDO_II_DEFINITION
                and candidate.alpha_allele
                and candidate.beta_allele
            ]
            if matching_pseudos:
                return tuple(
                    (candidate.canonical_name, candidate.alpha_allele, candidate.beta_allele)
                    for candidate in _deduplicate(matching_pseudos)
                )

            components = [part for part in re.split(r"\s*(?:/|\+)\s*", name) if part]
            if len(components) == 2:
                return (("/".join(components), components[0], components[1]),)
            deposited_name = canonical_hla(name)
            protein_group = two_field_hla(name)
            if deposited_name and protein_group and "-DRB" in protein_group:
                return (
                    (
                        f"HLA-DRA*01:01/{deposited_name}",
                        "HLA-DRA*01:01",
                        deposited_name,
                    ),
                )
        return ()

    def molecule_to_chains(
        self,
        name: str,
        *,
        mhc_class: str,
        species: str | None = None,
        include_partial: bool = False,
        allow_two_field_projection: bool = False,
        fallback: str = "none",
    ) -> QueryResult:
        """Resolve a complete MHC molecule into its two precursor chains.

        MHC-I returns heavy chain plus species-matched B2M.  MHC-II accepts an
        explicit alpha/beta pair, a bundled mouse molecule alias, or a DRB
        shorthand (which is explicitly paired with HLA-DRA*01:01).
        """

        expected_class = normalize_mhc_class(mhc_class)
        fallback_policy = normalize_fallback_policy(fallback)
        notes: List[str] = []
        molecules: List[MoleculeCandidate] = []
        projected = False
        if expected_class == "I":
            heavy = self.allele_to_sequence(
                name,
                mhc_class="I",
                species=species,
                include_partial=include_partial,
                allow_two_field_projection=allow_two_field_projection,
                fallback=fallback_policy,
            )
            projected = projected or heavy.status.startswith("projected_")
            notes.extend(f"heavy chain: {note}" for note in heavy.notes)
            for heavy_candidate in heavy.candidates:
                if heavy_candidate.chain_role != "heavy":
                    continue
                companion = self.allele_to_sequence(
                    f"B2M-{heavy_candidate.species}",
                    mhc_class="I",
                    species=heavy_candidate.species,
                    include_partial=include_partial,
                )
                for b2m_candidate in companion.candidates:
                    molecules.append(
                        MoleculeCandidate(
                            canonical_name=heavy_candidate.allele,
                            species=heavy_candidate.species,
                            mhc_class="I",
                            chains=(heavy_candidate, b2m_candidate),
                        )
                    )
        else:
            specs = self._pair_specs(name, expected_class)
            if len(re.split(r"\s*(?:/|\+)\s*", name)) == 1 and "DRB" in name.upper():
                notes.append("DRB shorthand was paired with the conventional HLA-DRA*01:01 alpha chain.")
            for canonical_name, alpha_name, beta_name in specs:
                alpha = self.allele_to_sequence(
                    alpha_name,
                    mhc_class="II",
                    species=species,
                    include_partial=include_partial,
                    allow_two_field_projection=allow_two_field_projection,
                    fallback=fallback_policy,
                )
                beta = self.allele_to_sequence(
                    beta_name,
                    mhc_class="II",
                    species=species,
                    include_partial=include_partial,
                    allow_two_field_projection=allow_two_field_projection,
                    fallback=fallback_policy,
                )
                projected = projected or alpha.status.startswith("projected_")
                projected = projected or beta.status.startswith("projected_")
                notes.extend(f"alpha chain: {note}" for note in alpha.notes)
                notes.extend(f"beta chain: {note}" for note in beta.notes)
                for alpha_candidate, beta_candidate in product(alpha.candidates, beta.candidates):
                    if alpha_candidate.species != beta_candidate.species:
                        continue
                    if alpha_candidate.chain_role != "alpha" or beta_candidate.chain_role != "beta":
                        continue
                    alpha_family = _mhcii_family(alpha_candidate.locus, "alpha")
                    beta_family = _mhcii_family(beta_candidate.locus, "beta")
                    if not alpha_family or alpha_family != beta_family:
                        continue
                    molecules.append(
                        MoleculeCandidate(
                            canonical_name=canonical_name,
                            species=alpha_candidate.species,
                            mhc_class="II",
                            chains=(alpha_candidate, beta_candidate),
                        )
                    )

        unique = {}
        for molecule in molecules:
            key = (molecule.canonical_name, tuple(chain.record_id for chain in molecule.chains))
            unique[key] = molecule
        candidates = tuple(unique[key] for key in sorted(unique))
        deduplicated_notes = tuple(dict.fromkeys(notes))
        evidence_levels = {
            _candidate_evidence(chain)
            for molecule in candidates
            for chain in molecule.chains
        }
        fallback_returned = any(
            hasattr(chain, "evidence_level")
            for molecule in candidates
            for chain in molecule.chains
        )
        partial_deposited_returned = any(
            isinstance(chain, SequenceCandidate) and not chain.is_complete
            for molecule in candidates
            for chain in molecule.chains
        )
        return QueryResult(
            query=name,
            representation=(
                "fallback_precursor_molecule"
                if fallback_returned
                else (
                    "deposited_molecule_with_partial_chain"
                    if partial_deposited_returned
                    else "precursor_full_molecule"
                )
            ),
            status=result_status(len(candidates), projected=projected),
            candidates=candidates,
            notes=deduplicated_notes,
            fallback_used=fallback_policy if fallback_returned else "none",
            evidence=(
                next(iter(evidence_levels))
                if len(evidence_levels) == 1
                else ("mixed" if evidence_levels else "source_exact")
            ),
        )

    def molecule_to_af3_chains(
        self,
        name: str,
        *,
        mhc_class: str,
        species: str | None = None,
        definition_id: str = AF3_DEFINITION,
        allow_two_field_projection: bool = False,
        fallback: str = "none",
        allow_expression_qualified: bool = False,
    ) -> QueryResult:
        """Resolve an MHC molecule into the two separate AF3 protein chains."""

        expected_class = normalize_mhc_class(mhc_class)
        fallback_policy = normalize_fallback_policy(fallback)
        notes: List[str] = []
        molecules: List[AF3MoleculeCandidate] = []
        projected = False
        if expected_class == "I":
            heavy = self.allele_to_af3(
                name,
                mhc_class="I",
                species=species,
                definition_id=definition_id,
                allow_two_field_projection=allow_two_field_projection,
                fallback=fallback_policy,
                allow_expression_qualified=allow_expression_qualified,
            )
            projected = heavy.status.startswith("projected_")
            notes.extend(f"heavy chain: {note}" for note in heavy.notes)
            for heavy_candidate in heavy.candidates:
                if heavy_candidate.chain_role != "heavy":
                    continue
                companion = self.allele_to_af3(
                    f"B2M-{heavy_candidate.species}",
                    mhc_class="I",
                    species=heavy_candidate.species,
                    definition_id=definition_id,
                )
                for b2m_candidate in companion.candidates:
                    molecules.append(
                        AF3MoleculeCandidate(
                            canonical_name=heavy_candidate.allele,
                            species=heavy_candidate.species,
                            mhc_class="I",
                            chains=(heavy_candidate, b2m_candidate),
                        )
                    )
        else:
            specs = self._pair_specs(name, expected_class)
            if len(re.split(r"\s*(?:/|\+)\s*", name)) == 1 and "DRB" in name.upper():
                notes.append(
                    "DRB shorthand was paired with the conventional HLA-DRA*01:01 alpha chain."
                )
            for canonical_name, alpha_name, beta_name in specs:
                alpha = self.allele_to_af3(
                    alpha_name,
                    mhc_class="II",
                    species=species,
                    definition_id=definition_id,
                    allow_two_field_projection=allow_two_field_projection,
                    fallback=fallback_policy,
                    allow_expression_qualified=allow_expression_qualified,
                )
                beta = self.allele_to_af3(
                    beta_name,
                    mhc_class="II",
                    species=species,
                    definition_id=definition_id,
                    allow_two_field_projection=allow_two_field_projection,
                    fallback=fallback_policy,
                    allow_expression_qualified=allow_expression_qualified,
                )
                projected = projected or alpha.status.startswith("projected_")
                projected = projected or beta.status.startswith("projected_")
                notes.extend(f"alpha chain: {note}" for note in alpha.notes)
                notes.extend(f"beta chain: {note}" for note in beta.notes)
                for alpha_candidate, beta_candidate in product(
                    alpha.candidates, beta.candidates
                ):
                    if alpha_candidate.species != beta_candidate.species:
                        continue
                    if (
                        alpha_candidate.chain_role != "alpha"
                        or beta_candidate.chain_role != "beta"
                    ):
                        continue
                    alpha_family = _mhcii_family(alpha_candidate.locus, "alpha")
                    beta_family = _mhcii_family(beta_candidate.locus, "beta")
                    if not alpha_family or alpha_family != beta_family:
                        continue
                    molecules.append(
                        AF3MoleculeCandidate(
                            canonical_name=canonical_name,
                            species=alpha_candidate.species,
                            mhc_class="II",
                            chains=(alpha_candidate, beta_candidate),
                        )
                    )

        unique = {}
        for molecule in molecules:
            key = (
                molecule.canonical_name,
                tuple(chain.record_id for chain in molecule.chains),
            )
            unique[key] = molecule
        candidates = tuple(unique[key] for key in sorted(unique))
        evidence_levels = {
            getattr(chain, "evidence_level", "source_exact")
            for molecule in candidates
            for chain in molecule.chains
        }
        fallback_returned = any(
            hasattr(chain, "evidence_level")
            for molecule in candidates
            for chain in molecule.chains
        )
        synthetic_returned = bool(
            evidence_levels & {IMPUTED_EVIDENCE, PSEUDO_SCAFFOLD_EVIDENCE}
        )
        return QueryResult(
            query=name,
            representation=(
                f"{IMPUTED_AF3_DEFINITION}_molecule"
                if synthetic_returned
                else f"{definition_id}_molecule"
            ),
            status=result_status(len(candidates), projected=projected),
            candidates=candidates,
            notes=tuple(dict.fromkeys(notes)),
            fallback_used=fallback_policy if fallback_returned else "none",
            evidence=(
                next(iter(evidence_levels))
                if len(evidence_levels) == 1
                else ("mixed" if evidence_levels else "source_exact")
            ),
        )

    def molecule_to_pbr_chains(
        self,
        name: str,
        *,
        mhc_class: str,
        species: str | None = None,
        allow_two_field_projection: bool = False,
        fallback: str = "none",
        allow_expression_qualified: bool = False,
        require_pseudo_match: bool = False,
    ) -> QueryResult:
        """Resolve an MHC molecule into peptide-binding structural regions.

        MHC-I returns only the heavy-chain alpha1+alpha2 region (no B2M).
        MHC-II returns alpha1 and beta1 as two independent protein chains.
        """

        expected_class = normalize_mhc_class(mhc_class)
        fallback_policy = normalize_fallback_policy(fallback)
        notes: List[str] = []
        molecules: List[PBRMoleculeCandidate] = []
        projected = False
        if expected_class == "I":
            heavy = self.allele_to_pbr(
                name,
                mhc_class="I",
                species=species,
                allow_two_field_projection=allow_two_field_projection,
                fallback=fallback_policy,
                allow_expression_qualified=allow_expression_qualified,
                require_pseudo_match=require_pseudo_match,
            )
            projected = heavy.status.startswith("projected_")
            notes.extend(f"heavy chain: {note}" for note in heavy.notes)
            for heavy_candidate in heavy.candidates:
                if heavy_candidate.chain_role != "heavy":
                    continue
                molecules.append(
                    PBRMoleculeCandidate(
                        canonical_name=heavy_candidate.allele,
                        species=heavy_candidate.species,
                        mhc_class="I",
                        chains=(heavy_candidate,),
                    )
                )
        else:
            specs = self._pair_specs(name, expected_class)
            if len(re.split(r"\s*(?:/|\+)\s*", name)) == 1 and "DRB" in name.upper():
                notes.append(
                    "DRB shorthand was paired with the conventional HLA-DRA*01:01 alpha chain."
                )
            for canonical_name, alpha_name, beta_name in specs:
                alpha = self.allele_to_pbr(
                    alpha_name,
                    mhc_class="II",
                    species=species,
                    allow_two_field_projection=allow_two_field_projection,
                    fallback=fallback_policy,
                    allow_expression_qualified=allow_expression_qualified,
                    require_pseudo_match=require_pseudo_match,
                )
                beta = self.allele_to_pbr(
                    beta_name,
                    mhc_class="II",
                    species=species,
                    allow_two_field_projection=allow_two_field_projection,
                    fallback=fallback_policy,
                    allow_expression_qualified=allow_expression_qualified,
                    require_pseudo_match=require_pseudo_match,
                )
                projected = projected or alpha.status.startswith("projected_")
                projected = projected or beta.status.startswith("projected_")
                notes.extend(f"alpha chain: {note}" for note in alpha.notes)
                notes.extend(f"beta chain: {note}" for note in beta.notes)
                for alpha_candidate, beta_candidate in product(
                    alpha.candidates, beta.candidates
                ):
                    if alpha_candidate.species != beta_candidate.species:
                        continue
                    if (
                        alpha_candidate.chain_role != "alpha"
                        or beta_candidate.chain_role != "beta"
                    ):
                        continue
                    alpha_family = _mhcii_family(alpha_candidate.locus, "alpha")
                    beta_family = _mhcii_family(beta_candidate.locus, "beta")
                    if not alpha_family or alpha_family != beta_family:
                        continue
                    molecules.append(
                        PBRMoleculeCandidate(
                            canonical_name=canonical_name,
                            species=alpha_candidate.species,
                            mhc_class="II",
                            chains=(alpha_candidate, beta_candidate),
                        )
                    )

        unique = {}
        for molecule in molecules:
            key = (
                molecule.canonical_name,
                tuple(chain.record_id for chain in molecule.chains),
            )
            unique[key] = molecule
        candidates = tuple(unique[key] for key in sorted(unique))
        evidence_levels = {
            getattr(chain, "evidence_level", "source_exact")
            for molecule in candidates
            for chain in molecule.chains
        }
        fallback_returned = any(
            isinstance(chain, FallbackPBRSequenceCandidate)
            for molecule in candidates
            for chain in molecule.chains
        )
        synthetic_returned = any(
            getattr(chain, "evidence_level", "")
            in {IMPUTED_EVIDENCE, PSEUDO_SCAFFOLD_EVIDENCE}
            for molecule in candidates
            for chain in molecule.chains
        )
        if require_pseudo_match and expected_class == "II" and candidates:
            notes.append(
                "Pseudo-match contract was independently certified for both the MHC-II alpha1 and beta1 PBR chains."
            )
        if require_pseudo_match:
            representation = (
                "pseudo_compatible_synthetic_pbr_molecule"
                if synthetic_returned
                else "pseudo_compatible_pbr_molecule"
            )
        else:
            representation = (
                "peptide_binding_region_fallback_molecule"
                if fallback_returned
                else "peptide_binding_region_molecule"
            )
        return QueryResult(
            query=name,
            representation=representation,
            status=result_status(len(candidates), projected=projected),
            candidates=candidates,
            notes=tuple(dict.fromkeys(notes)),
            fallback_used=(
                "pseudo-scaffold"
                if require_pseudo_match and fallback_returned
                else (fallback_policy if fallback_returned else "none")
            ),
            evidence=(
                next(iter(evidence_levels))
                if len(evidence_levels) == 1
                else ("mixed" if evidence_levels else "source_exact")
            ),
        )

    def stats(self) -> Dict[str, object]:
        return {
            "sequence_records": len(self._sequences),
            "complete_sequence_records": sum(item.is_complete for item in self._sequences),
            "af3_sequence_records": len(self._af3_sequences),
            "pbr_sequence_records": len(self._pbr_sequences),
            "pseudo_pbr_projection_records": len(self._pseudo_pbr_projections),
            "pseudo_records": len(self._pseudos),
            "pseudo_official_supported_records": sum(
                item.official_support for item in self._pseudos
            ),
            "pseudo_by_evidence": dict(
                sorted(Counter(item.evidence_level for item in self._pseudos).items())
            ),
            "pair_records": len(self._pairs),
            "fallback_sequence_records": len(self._fallback_sequences),
            "fallback_af3_records": len(self._fallback_af3_sequences),
            "fallback_pbr_records": len(self._fallback_pbr_sequences),
            "fallback_by_evidence": dict(
                sorted(
                    Counter(
                        item.evidence_level
                        for item in (
                            *self._fallback_sequences,
                            *self._fallback_af3_sequences,
                            *self._fallback_pbr_sequences,
                        )
                    ).items()
                )
            ),
            "metadata": self.metadata,
        }

    def coverage(
        self,
        name: str,
        *,
        mhc_class: str,
        species: str | None = None,
        allow_expression_qualified: bool = False,
    ) -> Dict[str, object]:
        """Summarize strict, partial, imputed and pseudo-scaffold availability."""

        expected_class = normalize_mhc_class(mhc_class)
        pseudo_definition = (
            PSEUDO_I_DEFINITION if expected_class == "I" else PSEUDO_II_DEFINITION
        )
        pseudo = self.allele_to_pseudo(
            name,
            mhc_class=expected_class,
            definition_id=pseudo_definition,
            species=species,
        )
        if expected_class == "I":
            full_calls = {
                "strict_full": self.allele_to_sequence(
                    name, mhc_class="I", species=species
                ),
                "partial": self.allele_to_sequence(
                    name, mhc_class="I", species=species, fallback="partial"
                ),
                "imputed_full": self.allele_to_sequence(
                    name, mhc_class="I", species=species, fallback="imputed"
                ),
                "pseudo_scaffold_full": self.allele_to_sequence(
                    name,
                    mhc_class="I",
                    species=species,
                    fallback="pseudo-scaffold",
                ),
            }
            af3_calls = {
                "strict_af3": self.allele_to_af3(name, mhc_class="I", species=species),
                "partial_af3": self.allele_to_af3(
                    name,
                    mhc_class="I",
                    species=species,
                    fallback="partial",
                    allow_expression_qualified=allow_expression_qualified,
                ),
                "imputed_af3": self.allele_to_af3(
                    name,
                    mhc_class="I",
                    species=species,
                    fallback="imputed",
                    allow_expression_qualified=allow_expression_qualified,
                ),
                "pseudo_scaffold_af3": self.allele_to_af3(
                    name,
                    mhc_class="I",
                    species=species,
                    fallback="pseudo-scaffold",
                    allow_expression_qualified=allow_expression_qualified,
                ),
            }
            pbr_calls = {
                "strict_pbr": self.allele_to_pbr(
                    name, mhc_class="I", species=species
                ),
                "partial_pbr": self.allele_to_pbr(
                    name,
                    mhc_class="I",
                    species=species,
                    fallback="partial",
                    allow_expression_qualified=allow_expression_qualified,
                ),
                "imputed_pbr": self.allele_to_pbr(
                    name,
                    mhc_class="I",
                    species=species,
                    fallback="imputed",
                    allow_expression_qualified=allow_expression_qualified,
                ),
                "pseudo_scaffold_pbr": self.allele_to_pbr(
                    name,
                    mhc_class="I",
                    species=species,
                    fallback="pseudo-scaffold",
                    allow_expression_qualified=allow_expression_qualified,
                ),
            }
        else:
            full_calls = {
                "strict_full": self.molecule_to_chains(
                    name, mhc_class="II", species=species
                ),
                "partial": self.molecule_to_chains(
                    name, mhc_class="II", species=species, fallback="partial"
                ),
                "imputed_full": self.molecule_to_chains(
                    name, mhc_class="II", species=species, fallback="imputed"
                ),
                "pseudo_scaffold_full": self.molecule_to_chains(
                    name,
                    mhc_class="II",
                    species=species,
                    fallback="pseudo-scaffold",
                ),
            }
            af3_calls = {
                "strict_af3": self.molecule_to_af3_chains(
                    name, mhc_class="II", species=species
                ),
                "partial_af3": self.molecule_to_af3_chains(
                    name,
                    mhc_class="II",
                    species=species,
                    fallback="partial",
                    allow_expression_qualified=allow_expression_qualified,
                ),
                "imputed_af3": self.molecule_to_af3_chains(
                    name,
                    mhc_class="II",
                    species=species,
                    fallback="imputed",
                    allow_expression_qualified=allow_expression_qualified,
                ),
                "pseudo_scaffold_af3": self.molecule_to_af3_chains(
                    name,
                    mhc_class="II",
                    species=species,
                    fallback="pseudo-scaffold",
                    allow_expression_qualified=allow_expression_qualified,
                ),
            }
            pbr_calls = {
                "strict_pbr": self.molecule_to_pbr_chains(
                    name, mhc_class="II", species=species
                ),
                "partial_pbr": self.molecule_to_pbr_chains(
                    name,
                    mhc_class="II",
                    species=species,
                    fallback="partial",
                    allow_expression_qualified=allow_expression_qualified,
                ),
                "imputed_pbr": self.molecule_to_pbr_chains(
                    name,
                    mhc_class="II",
                    species=species,
                    fallback="imputed",
                    allow_expression_qualified=allow_expression_qualified,
                ),
                "pseudo_scaffold_pbr": self.molecule_to_pbr_chains(
                    name,
                    mhc_class="II",
                    species=species,
                    fallback="pseudo-scaffold",
                    allow_expression_qualified=allow_expression_qualified,
                ),
            }

        def summarize(result: QueryResult) -> Dict[str, object]:
            return {
                "status": result.status,
                "candidate_count": result.candidate_count,
                "representation": result.representation,
                "fallback_used": result.fallback_used,
                "evidence": result.evidence,
                "notes": list(result.notes),
            }

        return {
            "query": name,
            "mhc_class": expected_class,
            "pseudo": summarize(pseudo),
            **{key: summarize(value) for key, value in full_calls.items()},
            **{key: summarize(value) for key, value in af3_calls.items()},
            **{key: summarize(value) for key, value in pbr_calls.items()},
        }

