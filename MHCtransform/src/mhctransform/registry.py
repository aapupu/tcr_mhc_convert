"""Reference-backed, non-lossy MHC lookup API."""

from __future__ import annotations

import csv
import json
import re
from collections import defaultdict
from importlib import resources
from itertools import product
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Mapping, Sequence, Tuple

from .errors import MHCTransformError
from .models import (
    AF3MoleculeCandidate,
    AF3SequenceCandidate,
    MoleculeCandidate,
    PseudoCandidate,
    QueryResult,
    SequenceCandidate,
    result_status,
)
from .names import (
    canonical_hla,
    exact_name_key,
    normalize_mhc_class,
    normalize_sequence,
    two_field_hla,
)


PSEUDO_I_DEFINITION = "netmhcpan-34-v1"
PSEUDO_II_DEFINITION = "netmhciipan-34-v1"
AF3_DEFINITION = "af3-soluble-extracellular-v1"


def _split_aliases(value: str) -> Tuple[str, ...]:
    return tuple(part for part in value.split(";") if part)


def _deduplicate(candidates: Iterable[object]) -> Tuple[object, ...]:
    by_id = {}
    for candidate in candidates:
        by_id[getattr(candidate, "record_id")] = candidate
    return tuple(by_id[key] for key in sorted(by_id))


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
        self._pairs = self._load_pairs(data_root.joinpath("molecule_pairs.csv"))
        with data_root.joinpath("metadata.json").open("r", encoding="utf-8") as handle:
            self.metadata: Mapping[str, object] = json.load(handle)

        self._sequence_name_index: DefaultDict[str, List[SequenceCandidate]] = defaultdict(list)
        self._sequence_value_index: DefaultDict[str, List[SequenceCandidate]] = defaultdict(list)
        for candidate in self._sequences:
            for name in (candidate.allele, *candidate.aliases):
                self._sequence_name_index[exact_name_key(name)].append(candidate)
            self._sequence_value_index[candidate.sequence].append(candidate)

        self._af3_name_index: DefaultDict[str, List[AF3SequenceCandidate]] = defaultdict(list)
        self._af3_value_index: DefaultDict[str, List[AF3SequenceCandidate]] = defaultdict(list)
        for candidate in self._af3_sequences:
            for name in (candidate.allele, *candidate.aliases):
                self._af3_name_index[exact_name_key(name)].append(candidate)
            self._af3_value_index[candidate.sequence].append(candidate)

        self._af3_exclusion_index: DefaultDict[str, List[Mapping[str, str]]] = defaultdict(list)
        for exclusion in self._af3_exclusions:
            for name in (exclusion["allele"], *_split_aliases(exclusion["aliases"])):
                self._af3_exclusion_index[exact_name_key(name)].append(exclusion)

        self._pseudo_name_index: DefaultDict[str, List[PseudoCandidate]] = defaultdict(list)
        self._pseudo_value_index: DefaultDict[str, List[PseudoCandidate]] = defaultdict(list)
        for candidate in self._pseudos:
            for name in (candidate.canonical_name, *candidate.aliases):
                self._pseudo_name_index[exact_name_key(name)].append(candidate)
            self._pseudo_value_index[candidate.pseudo_sequence].append(candidate)

        self._pair_index: DefaultDict[str, List[Mapping[str, str]]] = defaultdict(list)
        for pair in self._pairs:
            for name in (pair["canonical_name"], *_split_aliases(pair["aliases"])):
                self._pair_index[exact_name_key(name)].append(pair)

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
    def _load_pairs(path: object) -> Tuple[Mapping[str, str], ...]:
        with path.open("r", encoding="utf-8", newline="") as handle:
            return tuple(dict(row) for row in csv.DictReader(handle))

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

    def allele_to_sequence(
        self,
        allele: str,
        *,
        mhc_class: str | None = None,
        species: str | None = None,
        include_partial: bool = False,
        allow_two_field_projection: bool = False,
    ) -> QueryResult:
        """Convert a chain allele/name to deposited precursor protein sequence(s)."""

        exact = self._sequence_name_index.get(exact_name_key(allele), ())
        notes: List[str] = []
        projected = False
        candidates = self._filter_sequence_candidates(
            exact,
            mhc_class=mhc_class,
            species=species,
            include_partial=include_partial,
        )
        if not candidates and not exact and allow_two_field_projection:
            fallback = two_field_hla(allele)
            if fallback and exact_name_key(fallback) != exact_name_key(allele):
                candidates = self._filter_sequence_candidates(
                    self._sequence_name_index.get(exact_name_key(fallback), ()),
                    mhc_class=mhc_class,
                    species=species,
                    include_partial=include_partial,
                )
                if candidates:
                    projected = True
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
        return QueryResult(
            query=allele,
            representation="precursor_full",
            status=result_status(len(candidates), projected=projected),
            candidates=candidates,
            notes=tuple(notes),
        )

    def sequence_to_alleles(
        self,
        sequence: str,
        *,
        mhc_class: str | None = None,
        species: str | None = None,
        include_partial: bool = False,
    ) -> QueryResult:
        """Return all exact allele/name candidates for a precursor sequence."""

        normalized = normalize_sequence(sequence, allow_x=include_partial)
        candidates = self._filter_sequence_candidates(
            self._sequence_value_index.get(normalized, ()),
            mhc_class=mhc_class,
            species=species,
            include_partial=include_partial,
        )
        deposited_names = {
            alias
            for candidate in candidates
            for alias in candidate.aliases
            if alias.startswith("HLA-") and "*" in alias
        }
        reverse_status = result_status(len(candidates))
        notes: Tuple[str, ...] = ()
        if len(deposited_names) > 1:
            reverse_status = "exact_ambiguous"
            notes = (
                "Multiple deposited HLA allele names share this exact protein; they are compactly listed in candidate aliases.",
            )
        return QueryResult(
            query=normalized,
            representation="precursor_full",
            status=reverse_status,
            candidates=candidates,
            notes=notes,
        )

    def allele_to_af3(
        self,
        allele: str,
        *,
        mhc_class: str | None = None,
        species: str | None = None,
        definition_id: str = AF3_DEFINITION,
        allow_two_field_projection: bool = False,
    ) -> QueryResult:
        """Convert a chain allele/name to a mature soluble AF3 input chain."""

        query_key = exact_name_key(allele)
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
            fallback = two_field_hla(allele)
            if fallback and exact_name_key(fallback) != query_key:
                candidates = self._filter_af3_candidates(
                    self._af3_name_index.get(exact_name_key(fallback), ()),
                    mhc_class=mhc_class,
                    species=species,
                    definition_id=definition_id,
                )
                if candidates:
                    projected = True
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
        return QueryResult(
            query=allele,
            representation=definition_id,
            status=result_status(len(candidates), projected=projected),
            candidates=candidates,
            notes=tuple(notes),
        )

    def af3_to_alleles(
        self,
        sequence: str,
        *,
        mhc_class: str | None = None,
        species: str | None = None,
        definition_id: str = AF3_DEFINITION,
    ) -> QueryResult:
        """Return every exact chain allele/name sharing an AF3 input sequence."""

        normalized = normalize_sequence(sequence)
        candidates = self._filter_af3_candidates(
            self._af3_value_index.get(normalized, ()),
            mhc_class=mhc_class,
            species=species,
            definition_id=definition_id,
        )
        deposited_names = {
            alias
            for candidate in candidates
            for alias in candidate.aliases
            if alias.startswith("HLA-") and "*" in alias
        }
        reverse_status = result_status(len(candidates))
        notes: Tuple[str, ...] = ()
        if len(deposited_names) > 1:
            reverse_status = "exact_ambiguous"
            notes = (
                "Multiple deposited HLA allele names share this exact AF3 chain; they are retained in candidate aliases.",
            )
        return QueryResult(
            query=normalized,
            representation=definition_id,
            status=reverse_status,
            candidates=candidates,
            notes=notes,
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
        return QueryResult(
            query=allele_or_pair,
            representation=definition_id,
            status=result_status(len(candidates), projected=projected),
            candidates=candidates,
            notes=tuple(notes),
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
        return QueryResult(
            query=normalized,
            representation=definition_id,
            status=result_status(len(candidates)),
            candidates=candidates,
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
        include_partial: bool = False,
        allow_two_field_projection: bool = False,
    ) -> QueryResult:
        """Resolve a complete MHC molecule into its two precursor chains.

        MHC-I returns heavy chain plus species-matched B2M.  MHC-II accepts an
        explicit alpha/beta pair, a bundled mouse molecule alias, or a DRB
        shorthand (which is explicitly paired with HLA-DRA*01:01).
        """

        expected_class = normalize_mhc_class(mhc_class)
        notes: List[str] = []
        molecules: List[MoleculeCandidate] = []
        projected = False
        if expected_class == "I":
            heavy = self.allele_to_sequence(
                name,
                mhc_class="I",
                include_partial=include_partial,
                allow_two_field_projection=allow_two_field_projection,
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
                    include_partial=include_partial,
                    allow_two_field_projection=allow_two_field_projection,
                )
                beta = self.allele_to_sequence(
                    beta_name,
                    mhc_class="II",
                    include_partial=include_partial,
                    allow_two_field_projection=allow_two_field_projection,
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
        return QueryResult(
            query=name,
            representation="precursor_full_molecule",
            status=result_status(len(candidates), projected=projected),
            candidates=candidates,
            notes=deduplicated_notes,
        )

    def molecule_to_af3_chains(
        self,
        name: str,
        *,
        mhc_class: str,
        definition_id: str = AF3_DEFINITION,
        allow_two_field_projection: bool = False,
    ) -> QueryResult:
        """Resolve an MHC molecule into the two separate AF3 protein chains."""

        expected_class = normalize_mhc_class(mhc_class)
        notes: List[str] = []
        molecules: List[AF3MoleculeCandidate] = []
        projected = False
        if expected_class == "I":
            heavy = self.allele_to_af3(
                name,
                mhc_class="I",
                definition_id=definition_id,
                allow_two_field_projection=allow_two_field_projection,
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
                    definition_id=definition_id,
                    allow_two_field_projection=allow_two_field_projection,
                )
                beta = self.allele_to_af3(
                    beta_name,
                    mhc_class="II",
                    definition_id=definition_id,
                    allow_two_field_projection=allow_two_field_projection,
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
        return QueryResult(
            query=name,
            representation=f"{definition_id}_molecule",
            status=result_status(len(candidates), projected=projected),
            candidates=candidates,
            notes=tuple(dict.fromkeys(notes)),
        )

    def stats(self) -> Dict[str, object]:
        return {
            "sequence_records": len(self._sequences),
            "complete_sequence_records": sum(item.is_complete for item in self._sequences),
            "af3_sequence_records": len(self._af3_sequences),
            "pseudo_records": len(self._pseudos),
            "pair_records": len(self._pairs),
            "metadata": self.metadata,
        }

