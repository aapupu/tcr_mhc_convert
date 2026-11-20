#!/usr/bin/env python3
"""Build explicit partial/imputed fallback registries from strict assets.

This builder never edits the strict ``full_sequences.csv`` or
``af3_sequences.csv`` registries.  The generated constructs are isolated in
two additional files and carry residue-level provenance.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = PROJECT_ROOT / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from mhctransform import MHCRegistry  # noqa: E402
from mhctransform.fallback import (  # noqa: E402
    EXACT_ECTODOMAIN_EVIDENCE,
    IMPUTED_AF3_DEFINITION,
    IMPUTED_EVIDENCE,
    IMPUTED_PRECURSOR_DEFINITION,
    PARTIAL_EVIDENCE,
    PARTIAL_REPRESENTATION,
    PSEUDO_SCAFFOLD_EVIDENCE,
    apply_observed_sequence,
    apply_pseudosequence,
    compress_positions,
    first_field,
    hamming_distance,
    locate_contact_positions,
    local_alignment_map,
)
from mhctransform.names import exact_name_key  # noqa: E402
from mhctransform.registry import (  # noqa: E402
    AF3_DEFINITION,
    PSEUDO_I_DEFINITION,
    PSEUDO_II_DEFINITION,
)


FALLBACK_FIELDS = (
    "record_id",
    "allele",
    "aliases",
    "species",
    "mhc_class",
    "chain_role",
    "locus",
    "representation",
    "sequence",
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


def fallback_metadata_fields() -> Dict[str, object]:
    """Return the shared schema fragment for strict and standalone builders."""

    return {
        "schema_version": 5,
        "sequence_kind": (
            "precursor_full (strict default); partial_deposited_protein (legacy opt-in); "
            f"{AF3_DEFINITION} (strict mature extracellular chain); source-partial-v1; "
            "precursor-imputed-v1; af3-soluble-extracellular-imputed-v1 "
            "(all fallback representations require explicit opt-in)"
        ),
        "fallback_definitions": {
            "source-partial-v1": (
                "Exact deposited fragment; returned only through an explicit fallback policy."
            ),
            "source_exact_ectodomain": (
                "A deposited fragment that covers every residue of the selected exact "
                "extracellular chain; no AF3 residue is imputed."
            ),
            "precursor-imputed-v1": (
                "Synthetic full-length scaffold retaining aligned deposited residues and "
                "all non-X target pseudosequence residues; remaining positions come from "
                "a same-species, same-locus complete template."
            ),
            "af3-soluble-extracellular-imputed-v1": (
                "Synthetic AF3 ectodomain scaffold with deposited-observed and "
                "non-observed/imputed spans."
            ),
            "policy_order": ["none", "partial", "imputed", "pseudo-scaffold"],
            "warning": (
                "A 34-aa pseudosequence is not invertible. Imputed constructs are structural "
                "scaffolds, not official allele protein sequences."
            ),
        },
        "pseudo_pbr_contract": {
            "projection_file": "pseudo_pbr_projection.csv",
            "rule": (
                "Every pseudo component must map all contact indices to one strict or "
                "explicit fallback PBR. Every non-X pseudo residue must match the "
                "projected PBR residue exactly. DTU X alignment/unknown tokens are "
                "recorded as explicit wildcards and are never inserted into a protein chain."
            ),
            "failure_policy": "Any uncovered component or non-X mismatch fails the build.",
        },
    }


def update_fallback_metadata(data_dir: Path, counts: Mapping[str, object]) -> None:
    """Synchronize metadata when this builder is run on its own."""

    path = data_dir / "metadata.json"
    if not path.is_file():
        return
    with path.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    metadata.update(fallback_metadata_fields())
    metadata.setdefault("counts", {}).update(counts)
    metadata["derived_artifacts_pending_rebuild"] = False
    with path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def _record_id(*parts: object) -> str:
    digest = hashlib.sha256("\x1f".join(str(part) for part in parts).encode()).hexdigest()
    return digest[:20]


def _write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FALLBACK_FIELDS, lineterminator="\n"
        )
        writer.writeheader()
        for row in sorted(
            rows,
            key=lambda item: (
                str(item["species"]),
                str(item["mhc_class"]),
                str(item["locus"]),
                str(item["allele"]),
                str(item["evidence_level"]),
                str(item["record_id"]),
            ),
        ):
            writer.writerow(row)


def _target_locus(
    allele: str,
    species: str,
    mhc_class: str,
    role: str,
    deposited: Sequence[object],
) -> str:
    loci = {str(getattr(candidate, "locus")) for candidate in deposited}
    if len(loci) == 1:
        return next(iter(loci))
    text = allele.upper().removeprefix("HLA-")
    if species == "human":
        return text.split("*", 1)[0]
    if mhc_class == "I":
        match = re.search(r"H[-_]?2[-_]?([DKL])", allele, flags=re.IGNORECASE)
        if not match:
            raise ValueError(f"cannot infer mouse MHC-I locus for {allele}")
        return {"D": "H2-D1", "K": "H2-K1", "L": "H2-L"}[match.group(1).upper()]
    return text.split("*", 1)[0]


def _quality_flags(
    allele: str,
    evidence: str,
    unknown: Sequence[int],
    alignment_flags: Iterable[str] = (),
) -> Tuple[str, ...]:
    flags = set(alignment_flags)
    if evidence in {IMPUTED_EVIDENCE, PSEUDO_SCAFFOLD_EVIDENCE}:
        flags.add("synthetic_construct_not_official_sequence")
        flags.add("template_derived_non_observed_residues")
    if evidence == PARTIAL_EVIDENCE:
        flags.add("deposited_fragment_not_full_length")
    if unknown:
        flags.add("pseudo_contains_unspecified_x")
    if re.search(r"[NLSQCA](?:/|$)", allele):
        flags.add("expression_qualified_target")
    return tuple(sorted(flags))


def _row(
    *,
    allele: str,
    aliases: Iterable[str],
    species: str,
    mhc_class: str,
    chain_role: str,
    locus: str,
    representation: str,
    sequence: str,
    evidence_level: str,
    observed_positions: Iterable[int],
    imputed_positions: Iterable[int] | None = None,
    source_record_ids: Iterable[str],
    template_record_id: str,
    template_allele: str,
    target_pseudo_record_id: str,
    pseudo_distance: int,
    unknown_pseudo_positions: Iterable[int],
    quality_flags: Iterable[str],
) -> Dict[str, object]:
    observed = tuple(sorted(set(int(value) for value in observed_positions)))
    if imputed_positions is None:
        observed_set = set(observed)
        imputed = tuple(
            position
            for position in range(1, len(sequence) + 1)
            if position not in observed_set
        )
    else:
        imputed = tuple(sorted(set(int(value) for value in imputed_positions)))
    source_ids = tuple(sorted(set(value for value in source_record_ids if value)))
    unknown = tuple(sorted(set(int(value) for value in unknown_pseudo_positions)))
    clean_aliases = tuple(sorted(set(value for value in aliases if value and value != allele)))
    flags = tuple(sorted(set(quality_flags)))
    return {
        "record_id": _record_id(
            "fallback",
            representation,
            species,
            allele,
            evidence_level,
            sequence,
            source_ids,
            template_record_id,
        ),
        "allele": allele,
        "aliases": ";".join(clean_aliases),
        "species": species,
        "mhc_class": mhc_class,
        "chain_role": chain_role,
        "locus": locus,
        "representation": representation,
        "sequence": sequence,
        "evidence_level": evidence_level,
        "observed_spans_1based": compress_positions(observed),
        "imputed_spans_1based": compress_positions(imputed),
        "observed_residue_count": len(observed),
        "imputed_residue_count": len(imputed),
        "source_record_ids": ";".join(source_ids),
        "template_record_id": template_record_id,
        "template_allele": template_allele,
        "target_pseudo_record_id": target_pseudo_record_id,
        "pseudo_distance": pseudo_distance,
        "unknown_pseudo_positions": ";".join(str(value) for value in unknown),
        "quality_flags": ";".join(flags),
    }


def _component_targets(registry: MHCRegistry) -> List[Dict[str, object]]:
    targets: MutableMapping[Tuple[str, str, str, str], Dict[str, object]] = {}
    for pseudo in registry._pseudos:
        if pseudo.mhc_class == "I":
            components = (("heavy", pseudo.canonical_name, pseudo.pseudo_sequence, 0),)
        else:
            components = (
                ("alpha", pseudo.alpha_allele, pseudo.pseudo_sequence[:15], 0),
                ("beta", pseudo.beta_allele, pseudo.pseudo_sequence[15:], 15),
            )
        for role, allele, segment, pseudo_offset in components:
            key = (pseudo.species, pseudo.mhc_class, role, exact_name_key(allele))
            deposited = registry._sequence_name_index.get(exact_name_key(allele), ())
            target = targets.setdefault(
                key,
                {
                    "allele": allele,
                    "aliases": set(),
                    "species": pseudo.species,
                    "mhc_class": pseudo.mhc_class,
                    "chain_role": role,
                    "locus": _target_locus(
                        allele, pseudo.species, pseudo.mhc_class, role, deposited
                    ),
                    "pseudo_segment": segment,
                    "pseudo_offset": pseudo_offset,
                    "pseudo_record_ids": set(),
                },
            )
            if target["pseudo_segment"] != segment:
                raise ValueError(f"conflicting component pseudo for {allele}")
            target["pseudo_record_ids"].add(pseudo.record_id)  # type: ignore[union-attr]
            if pseudo.mhc_class == "I":
                target["aliases"].update(pseudo.aliases)  # type: ignore[union-attr]
            for candidate in deposited:
                target["aliases"].update(candidate.aliases)  # type: ignore[union-attr]
    return list(targets.values())


def _template_catalog(
    registry: MHCRegistry,
    targets: Sequence[Mapping[str, object]],
) -> Mapping[Tuple[str, str, str, str], Tuple[Tuple[object, str], ...]]:
    full_by_id = {candidate.record_id: candidate for candidate in registry._sequences}
    target_by_key = {
        (
            str(item["species"]),
            str(item["mhc_class"]),
            str(item["chain_role"]),
            exact_name_key(str(item["allele"])),
        ): item
        for item in targets
    }
    grouped: MutableMapping[Tuple[str, str, str, str], List[Tuple[object, str]]] = defaultdict(list)
    for candidate in registry._af3_sequences:
        # A strict AF3 ectodomain may legitimately come from a signal-free
        # mature membrane record (P14437). Such a chain is valid for direct AF3
        # use, but its parent is not a complete precursor scaffold. Synthetic
        # full sequences must always inherit signal/TM/tail context from a
        # source record that passed the strict completeness rules.
        parent = full_by_id.get(candidate.source_full_record_id)
        if parent is None or not parent.is_complete:
            continue
        keys = [exact_name_key(candidate.allele), *(exact_name_key(alias) for alias in candidate.aliases)]
        target = next(
            (
                target_by_key.get(
                    (candidate.species, candidate.mhc_class, candidate.chain_role, key)
                )
                for key in keys
                if (
                    candidate.species,
                    candidate.mhc_class,
                    candidate.chain_role,
                    key,
                )
                in target_by_key
            ),
            None,
        )
        if target is None:
            continue
        segment = str(target["pseudo_segment"])
        if "X" in segment:
            continue
        group = (candidate.species, candidate.mhc_class, candidate.chain_role, candidate.locus)
        grouped[group].append((candidate, segment))
    return {
        key: tuple(sorted(value, key=lambda item: (item[0].allele, item[0].record_id)))
        for key, value in grouped.items()
    }


def _choose_template(
    target: Mapping[str, object],
    catalog: Mapping[Tuple[str, str, str, str], Tuple[Tuple[object, str], ...]],
) -> Tuple[object, str, int]:
    group = (
        str(target["species"]),
        str(target["mhc_class"]),
        str(target["chain_role"]),
        str(target["locus"]),
    )
    candidates = catalog.get(group, ())
    if not candidates:
        raise ValueError(f"no same-locus AF3 template for {target['allele']}")
    target_segment = str(target["pseudo_segment"])
    target_first = first_field(str(target["allele"]))
    ranked = []
    for candidate, segment in candidates:
        distance = hamming_distance(target_segment, segment)
        ranked.append(
            (
                0 if target_first and first_field(candidate.allele) == target_first else 1,
                distance,
                candidate.allele,
                candidate.record_id,
                candidate,
                segment,
            )
        )
    selected = min(ranked)
    return selected[4], selected[5], selected[1]


def _strict_pbr_for_af3(registry: MHCRegistry, af3_candidate: object) -> object:
    """Return the unique strict PBR cropped from one strict AF3 chain."""

    index = getattr(registry, "_strict_pbr_by_af3_record_id", None)
    if index is None:
        grouped: MutableMapping[str, List[object]] = defaultdict(list)
        for candidate in registry._pbr_sequences:
            grouped[candidate.source_af3_record_id].append(candidate)
        index = {key: tuple(value) for key, value in grouped.items()}
        setattr(registry, "_strict_pbr_by_af3_record_id", index)
    candidates = index.get(af3_candidate.record_id, ())
    if len(candidates) != 1:
        raise ValueError(
            f"AF3 record {af3_candidate.record_id} resolves to "
            f"{len(candidates)} strict PBR records"
        )
    return candidates[0]


def _contact_positions(
    registry: MHCRegistry,
    template: object,
    template_segment: str,
) -> Tuple[int, ...]:
    """Project one pseudo component onto its same-locus strict PBR.

    Pseudo contacts are PBR coordinates.  Both MHC classes therefore use the
    strict PBR crop here; using a complete MHC-II AF3 ectodomain in this path
    would make fallback compatibility use a different coordinate substrate
    from the final pseudo/PBR certificate.
    """

    target_pbr = _strict_pbr_for_af3(registry, template)
    if template.mhc_class == "I":
        positions = registry.metadata["pseudo_definitions"][PSEUDO_I_DEFINITION][
            "mature_reference_positions_1_based"
        ]
        reference_species = template.species
        reference_locus = template.locus
    else:
        locus_positions = registry.metadata.get("pseudo_derivation", {}).get(
            "mhc_ii_component_contact_positions_1based", {}
        )
        if template.locus in locus_positions:
            positions = locus_positions[template.locus]
            reference_species = template.species
            reference_locus = template.locus
        else:
            definition = registry.metadata["pseudo_definitions"][PSEUDO_II_DEFINITION]
            positions = definition[
                "alpha_mature_reference_positions_1_based"
                if template.chain_role == "alpha"
                else "beta_mature_reference_positions_1_based"
            ]
            # Mouse components currently use the published generic DR alpha or
            # beta coordinate frame, exactly as the final projection builder
            # does; no mouse-locus-specific contact list is defined.
            reference_species = "human"
            reference_locus = "DRA" if template.chain_role == "alpha" else "DRB1"

    reference_spec = next(
        (
            item
            for item in registry.metadata["pbr_definitions"]["reference_records"]
            if item["species"] == reference_species
            and item["locus"] == reference_locus
        ),
        None,
    )
    if reference_spec is None:
        raise ValueError(f"missing PBR reference for {template.allele}")
    reference_result = registry.allele_to_pbr(
        reference_spec["allele"],
        mhc_class=template.mhc_class,
        species=reference_species,
    )
    reference_candidates = tuple(
        candidate
        for candidate in reference_result.candidates
        if candidate.locus == reference_locus
        and candidate.chain_role == template.chain_role
    )
    if len(reference_candidates) != 1:
        raise ValueError(
            f"PBR reference {reference_spec['allele']} resolves to "
            f"{len(reference_candidates)} {reference_locus}/{template.chain_role} records"
        )
    try:
        return locate_contact_positions(
            reference_candidates[0].sequence,
            target_pbr.sequence,
            positions,
            template_segment,
        )
    except ValueError as error:
        raise ValueError(
            f"cannot map contacts for template {template.allele} "
            f"({template.record_id}, pseudo={template_segment}): {error}"
        ) from error


def _strict_candidates(registry: MHCRegistry, target: Mapping[str, object]) -> Tuple[object, ...]:
    key = exact_name_key(str(target["allele"]))
    return tuple(
        candidate
        for candidate in registry._af3_name_index.get(key, ())
        if candidate.species == target["species"]
        and candidate.mhc_class == target["mhc_class"]
        and candidate.chain_role == target["chain_role"]
    )


def _pseudo_compatible_strict_candidates(
    registry: MHCRegistry,
    target: Mapping[str, object],
) -> Tuple[object, ...]:
    """Return strict chains whose PBR contacts reproduce the target pseudo.

    A deposited chain remains a valid strict biological record even when an
    upstream DTU pseudosequence snapshot disagrees with the current deposited
    sequence.  Such a chain must not, however, satisfy the pseudo-to-PBR
    contract.  In that case the fallback builder emits a separately labelled
    pseudo-constrained scaffold instead of hiding the disagreement.
    """

    compatible = []
    segment = str(target["pseudo_segment"])
    for candidate in _strict_candidates(registry, target):
        try:
            positions = _contact_positions(
                registry,
                candidate,
                segment,
            )
        except ValueError:
            continue
        candidate_pbr = _strict_pbr_for_af3(registry, candidate)
        projected = "".join(
            candidate_pbr.sequence[position - 1] for position in positions
        )
        if all(
            expected == "X" or expected == observed
            for expected, observed in zip(segment, projected)
        ):
            compatible.append(candidate)
    return tuple(compatible)


def _deposited_candidates(registry: MHCRegistry, target: Mapping[str, object]) -> Tuple[object, ...]:
    key = exact_name_key(str(target["allele"]))
    unique = {}
    for candidate in registry._sequence_name_index.get(key, ()):
        if (
            candidate.species == target["species"]
            and candidate.mhc_class == target["mhc_class"]
            and candidate.chain_role == target["chain_role"]
        ):
            unique[candidate.record_id] = candidate
    return tuple(unique[key] for key in sorted(unique))


def build_fallback_files(data_dir: Path) -> Mapping[str, object]:
    registry = MHCRegistry(data_dir)
    targets = _component_targets(registry)
    catalog = _template_catalog(registry, targets)
    full_by_id = {candidate.record_id: candidate for candidate in registry._sequences}

    # A tiny number of secondary-mirror pseudo rows conflict with their exact
    # deposited AF3 chain. They are unsuitable as donors and are excluded here;
    # target rows remain constructible from another same-locus donor.
    contact_cache: Dict[str, Tuple[int, ...]] = {}
    valid_catalog: Dict[Tuple[str, str, str, str], Tuple[Tuple[object, str], ...]] = {}
    for group, entries in catalog.items():
        valid_entries = []
        for template, segment in entries:
            try:
                contact_cache[template.record_id] = _contact_positions(
                    registry,
                    template,
                    segment,
                )
            except ValueError:
                continue
            valid_entries.append((template, segment))
        valid_catalog[group] = tuple(valid_entries)
    catalog = valid_catalog

    # Thousands of allele names collapse to far fewer (first-field,
    # pseudosegment) template-selection problems.  Cache both the same-family
    # candidate pool and its deterministic nearest template; otherwise every
    # rare allele rescans an entire HLA locus and fallback rebuilding becomes
    # unnecessarily quadratic.
    template_pool_cache: Dict[
        Tuple[str, str, str, str, str], Tuple[Tuple[object, str], ...]
    ] = {}
    template_choice_cache: Dict[
        Tuple[str, str, str, str, str, str], Tuple[object, str, int]
    ] = {}

    def choose_template_cached(
        target: Mapping[str, object],
    ) -> Tuple[object, str, int]:
        group = (
            str(target["species"]),
            str(target["mhc_class"]),
            str(target["chain_role"]),
            str(target["locus"]),
        )
        target_first = first_field(str(target["allele"]))
        pool_key = (*group, target_first)
        if pool_key not in template_pool_cache:
            all_candidates = catalog.get(group, ())
            same_first = tuple(
                item
                for item in all_candidates
                if target_first and first_field(item[0].allele) == target_first
            )
            template_pool_cache[pool_key] = same_first or all_candidates
        target_segment = str(target["pseudo_segment"])
        choice_key = (*pool_key, target_segment)
        if choice_key not in template_choice_cache:
            candidates = template_pool_cache[pool_key]
            if not candidates:
                raise ValueError(f"no same-locus AF3 template for {target['allele']}")
            ranked = [
                (
                    hamming_distance(target_segment, segment),
                    candidate.allele,
                    candidate.record_id,
                    candidate,
                    segment,
                )
                for candidate, segment in candidates
            ]
            selected = min(ranked)
            template_choice_cache[choice_key] = (
                selected[3],
                selected[4],
                selected[0],
            )
        return template_choice_cache[choice_key]

    fallback_full: List[Dict[str, object]] = []
    fallback_af3: List[Dict[str, object]] = []
    curated_af3: Dict[str, Tuple[str, int]] = {"P14437": (EXACT_ECTODOMAIN_EVIDENCE, 195)}

    for target in targets:
        allele = str(target["allele"])
        aliases = set(target["aliases"])  # type: ignore[arg-type]
        pseudo_segment = str(target["pseudo_segment"])
        pseudo_offset = int(target["pseudo_offset"])
        pseudo_record_id = sorted(target["pseudo_record_ids"])[0]  # type: ignore[arg-type]
        deposited = _deposited_candidates(registry, target)
        strict_af3 = _pseudo_compatible_strict_candidates(
            registry,
            target,
        )

        complete_deposited = [candidate for candidate in deposited if candidate.is_complete]
        if not complete_deposited:
            for candidate in deposited:
                observed = tuple(
                    position
                    for position, residue in enumerate(candidate.sequence, start=1)
                    if residue != "X"
                )
                unknown = tuple(
                    pseudo_offset + index
                    for index, residue in enumerate(pseudo_segment, start=1)
                    if residue == "X"
                )
                fallback_full.append(
                    _row(
                        allele=allele,
                        aliases=aliases,
                        species=str(target["species"]),
                        mhc_class=str(target["mhc_class"]),
                        chain_role=str(target["chain_role"]),
                        locus=str(target["locus"]),
                        representation=PARTIAL_REPRESENTATION,
                        sequence=candidate.sequence,
                        evidence_level=PARTIAL_EVIDENCE,
                        observed_positions=observed,
                        imputed_positions=(),
                        source_record_ids=(candidate.record_id,),
                        template_record_id="",
                        template_allele="",
                        target_pseudo_record_id=pseudo_record_id,
                        pseudo_distance=0,
                        unknown_pseudo_positions=unknown,
                        quality_flags=(
                            *_quality_flags(allele, PARTIAL_EVIDENCE, unknown),
                            *(
                                ("source_contains_unknown_x",)
                                if "X" in candidate.sequence
                                else ()
                            ),
                        ),
                    )
                )

        curated = None
        for candidate in deposited:
            if not strict_af3 and candidate.source_accession in curated_af3:
                evidence, end = curated_af3[candidate.source_accession]
                sequence = candidate.sequence[:end]
                curated = _row(
                    allele=allele,
                    aliases=(*aliases, *candidate.aliases),
                    species=str(target["species"]),
                    mhc_class=str(target["mhc_class"]),
                    chain_role=str(target["chain_role"]),
                    locus=str(target["locus"]),
                    representation=AF3_DEFINITION,
                    sequence=sequence,
                    evidence_level=evidence,
                    observed_positions=range(1, len(sequence) + 1),
                    source_record_ids=(candidate.record_id,),
                    template_record_id="",
                    template_allele="",
                    target_pseudo_record_id=pseudo_record_id,
                    pseudo_distance=0,
                    unknown_pseudo_positions=(),
                    quality_flags=("exact_uniprot_extracellular_topology",),
                )
                fallback_af3.append(curated)
                break

        if strict_af3 or curated is not None:
            if complete_deposited:
                continue
            # A target may have an exact extracellular chain but no complete
            # precursor (notably signal-free mature mouse H2-Aa*s/P14437).
            # Preserve that exact ectodomain and add only the missing precursor
            # context from a full-capable same-locus template.
            template, _template_segment, pseudo_distance = choose_template_cached(target)
            template_full = full_by_id.get(template.source_full_record_id)
            if template_full is None or not template_full.is_complete:
                raise ValueError(
                    f"template {template.record_id} lacks a strict complete precursor parent"
                )
            exact_af3_sequence = (
                strict_af3[0].sequence
                if strict_af3
                else str(curated["sequence"])
            )
            if len(exact_af3_sequence) != len(template.sequence):
                raise ValueError(
                    f"exact AF3 chain/template length mismatch for {allele}: "
                    f"{len(exact_af3_sequence)} != {len(template.sequence)}"
                )
            start = template.precursor_start_1based
            end = start + len(exact_af3_sequence) - 1
            unknown = tuple(
                pseudo_offset + index
                for index, residue in enumerate(pseudo_segment, start=1)
                if residue == "X"
            )
            for source_candidate in deposited or (None,):
                source_sequence = (
                    source_candidate.sequence if source_candidate is not None else "X"
                )
                full_sequence, observed_full, alignment_flags = apply_observed_sequence(
                    template_full.sequence, source_sequence
                )
                full_residues = list(full_sequence)
                full_residues[start - 1 : end] = list(exact_af3_sequence)
                full_sequence = "".join(full_residues)
                source_ids = (
                    (source_candidate.record_id,)
                    if source_candidate is not None
                    else ()
                )
                evidence = (
                    IMPUTED_EVIDENCE if source_ids else PSEUDO_SCAFFOLD_EVIDENCE
                )
                flags = _quality_flags(
                    allele,
                    evidence,
                    unknown,
                    (
                        *alignment_flags,
                        "exact_ectodomain_grafted_into_complete_precursor_template",
                    ),
                )
                fallback_full.append(
                    _row(
                        allele=allele,
                        aliases=aliases,
                        species=str(target["species"]),
                        mhc_class=str(target["mhc_class"]),
                        chain_role=str(target["chain_role"]),
                        locus=str(target["locus"]),
                        representation=IMPUTED_PRECURSOR_DEFINITION,
                        sequence=full_sequence,
                        evidence_level=evidence,
                        observed_positions=observed_full,
                        source_record_ids=source_ids,
                        template_record_id=template_full.record_id,
                        template_allele=template_full.allele,
                        target_pseudo_record_id=pseudo_record_id,
                        pseudo_distance=pseudo_distance,
                        unknown_pseudo_positions=unknown,
                        quality_flags=flags,
                    )
                )
            continue

        template, template_segment, pseudo_distance = choose_template_cached(target)
        contact_positions = contact_cache[template.record_id]
        template_full = full_by_id.get(template.source_full_record_id)
        if template_full is None or not template_full.is_complete:
            raise ValueError(
                f"template {template.record_id} lacks a strict complete precursor parent"
            )

        source_candidates: Sequence[object | None] = deposited or (None,)
        for source_candidate in source_candidates:
            af3_sequence = template.sequence
            observed_af3: Tuple[int, ...] = ()
            alignment_flags: Tuple[str, ...] = ()
            source_ids: Tuple[str, ...] = ()
            source_grafted = False
            if source_candidate is not None:
                af3_sequence, observed_af3, alignment_flags = apply_observed_sequence(
                    af3_sequence, source_candidate.sequence
                )
                if "X" in source_candidate.sequence:
                    alignment_flags = (*alignment_flags, "source_contains_unknown_x")
                source_ids = (source_candidate.record_id,)
                source_grafted = True
            try:
                af3_sequence, _constrained_af3, unknown_local = apply_pseudosequence(
                    af3_sequence, contact_positions, pseudo_segment, observed_af3
                )
            except ValueError:
                # A legacy pseudo mirror can conflict with a deposited fragment.
                # Keep the pseudosequence contract and expose that the fragment
                # could not safely guide this construct.
                af3_sequence, _constrained_af3, unknown_local = apply_pseudosequence(
                    template.sequence, contact_positions, pseudo_segment, ()
                )
                observed_af3 = ()
                alignment_flags = (*alignment_flags, "partial_pseudo_conflict_not_grafted")
                source_grafted = False
            unknown = tuple(pseudo_offset + value for value in unknown_local)
            if source_grafted and len(observed_af3) == len(af3_sequence):
                af3_evidence = EXACT_ECTODOMAIN_EVIDENCE
                af3_representation = AF3_DEFINITION
                alignment_flags = (
                    *alignment_flags,
                    "deposited_fragment_covers_exact_ectodomain",
                )
            elif source_grafted:
                af3_evidence = IMPUTED_EVIDENCE
                af3_representation = IMPUTED_AF3_DEFINITION
            else:
                af3_evidence = PSEUDO_SCAFFOLD_EVIDENCE
                af3_representation = IMPUTED_AF3_DEFINITION
            flags = _quality_flags(
                allele, af3_evidence, unknown, alignment_flags
            )
            af3_row = _row(
                allele=allele,
                aliases=aliases,
                species=str(target["species"]),
                mhc_class=str(target["mhc_class"]),
                chain_role=str(target["chain_role"]),
                locus=str(target["locus"]),
                representation=af3_representation,
                sequence=af3_sequence,
                evidence_level=af3_evidence,
                observed_positions=observed_af3,
                source_record_ids=source_ids,
                template_record_id=template.record_id,
                template_allele=template.allele,
                target_pseudo_record_id=pseudo_record_id,
                pseudo_distance=pseudo_distance,
                unknown_pseudo_positions=unknown,
                quality_flags=flags,
            )
            fallback_af3.append(af3_row)

            if complete_deposited:
                continue
            full_evidence = (
                IMPUTED_EVIDENCE if source_grafted else PSEUDO_SCAFFOLD_EVIDENCE
            )
            full_scaffold, observed_full, full_alignment_flags = apply_observed_sequence(
                template_full.sequence,
                source_candidate.sequence if source_grafted else "X",
            )
            full_residues = list(full_scaffold)
            start = template.precursor_start_1based
            end = start + len(af3_sequence) - 1
            full_residues[start - 1 : end] = list(af3_sequence)
            full_sequence = "".join(full_residues)
            observed_full_set = {
                position
                for position in observed_full
                if position < start or position > end
            }
            observed_full_set.update(start + position - 1 for position in observed_af3)
            full_flags = _quality_flags(
                allele,
                full_evidence,
                unknown,
                (*alignment_flags, *full_alignment_flags),
            )
            fallback_full.append(
                _row(
                    allele=allele,
                    aliases=aliases,
                    species=str(target["species"]),
                    mhc_class=str(target["mhc_class"]),
                    chain_role=str(target["chain_role"]),
                    locus=str(target["locus"]),
                    representation=IMPUTED_PRECURSOR_DEFINITION,
                    sequence=full_sequence,
                    evidence_level=full_evidence,
                    observed_positions=observed_full_set,
                    source_record_ids=source_ids,
                    template_record_id=template_full.record_id,
                    template_allele=template_full.allele,
                    target_pseudo_record_id=pseudo_record_id,
                    pseudo_distance=pseudo_distance,
                    unknown_pseudo_positions=unknown,
                    quality_flags=full_flags,
                )
            )

    # Stable record-level de-duplication.
    full_unique = {str(row["record_id"]): row for row in fallback_full}
    af3_unique = {str(row["record_id"]): row for row in fallback_af3}
    _write_csv(data_dir / "fallback_full_sequences.csv", tuple(full_unique.values()))
    _write_csv(data_dir / "fallback_af3_sequences.csv", tuple(af3_unique.values()))
    # PBR fallback rows are a provenance-preserving crop of these finalized
    # AF3 rows, so they must be rebuilt only after fallback_af3_sequences.csv
    # has been written.
    from build_pbr_reference import (
        build_fallback_pbr_file,
        build_pseudo_pbr_projection_file,
    )

    fallback_pbr_count = build_fallback_pbr_file(data_dir)
    projection_counts = build_pseudo_pbr_projection_file(data_dir)
    return {
        "fallback_full": len(full_unique),
        "fallback_af3": len(af3_unique),
        "fallback_pbr": fallback_pbr_count,
        "component_targets": len(targets),
        **projection_counts,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=PROJECT_ROOT / "src" / "mhctransform" / "data",
    )
    args = parser.parse_args()
    counts = build_fallback_files(args.data_dir)
    update_fallback_metadata(args.data_dir, counts)
    for key, value in sorted(counts.items()):
        print(f"{key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
