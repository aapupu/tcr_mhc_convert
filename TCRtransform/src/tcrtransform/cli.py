"""Command-line interface for CSV TCR transformation workflows."""

from __future__ import annotations

import argparse
import csv
import os
import stat
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

from .errors import TCRTransformError
from .junctions import JunctionMetadata, load_junction_reference
from .references import (
    VGeneLoopReference,
    load_gene_sequences,
    normalize_species,
)
from .regions import CDR3Boundary, annotate_v_gene, define_cdr3_boundary
from .reverse import SpeciesReference, parse_chain
from .transform import reconstruct_tcr_chain


BUILD_ALIASES = {
    "trav": ("TRAV", "trav", "v_a_gene", "va_gene"),
    "traj": ("TRAJ", "traj", "j_a_gene", "ja_gene"),
    "cdr3a": ("CDR3A", "CDR3a", "cdr3a", "cdr3a_aa"),
    "trbv": ("TRBV", "trbv", "v_b_gene", "vb_gene"),
    "trbj": ("TRBJ", "trbj", "j_b_gene", "jb_gene"),
    "cdr3b": ("CDR3B", "CDR3b", "cdr3b", "cdr3b_aa"),
}
ALPHA_SEQUENCE_ALIASES = (
    "TCRa", "tcra", "TCRA", "A", "tcr_alpha", "alpha_sequence"
)
BETA_SEQUENCE_ALIASES = (
    "TCRb", "tcrb", "TCRB", "B", "tcr_beta", "beta_sequence"
)


def _boundary_columns(suffix: str) -> list[str]:
    base = f"CDR3{suffix}"
    return [
        f"{base}_v_supported",
        f"{base}_non_vj",
        f"{base}_j_supported",
        f"{base}_keep_prefix",
        f"{base}_keep_suffix",
        f"{base}_non_vj_start",
        f"{base}_non_vj_end",
        f"{base}_non_vj_length",
        f"{base}_fixed_mask",
        f"{base}_v_gene_start",
        f"{base}_v_gene_end",
        f"{base}_v_offset_from_end",
        f"{base}_j_gene_start",
        f"{base}_j_gene_end",
        f"{base}_j_offset_from_start",
        f"{base}_v_short_match_length",
        f"{base}_j_short_match_length",
        f"{base}_v_j_support_overlap",
        f"{base}_safe_for_optimization",
        f"{base}_boundary_status",
    ]


REGION_COLUMNS = [
    "tcr_species",
    "CDR1a", "CDR2a", "CDR2.5a", "CDR3a",
    *_boundary_columns("a"),
    "CDR1b", "CDR2b", "CDR2.5b", "CDR3b",
    *_boundary_columns("b"),
]
BUILD_COLUMNS = [
    *REGION_COLUMNS,
    "TCRa", "TCRb",
    "TCRa_reference_resolution", "TCRb_reference_resolution",
    "TCRa_reference_warnings", "TCRb_reference_warnings",
    "TCRa_productivity_status", "TCRb_productivity_status",
    "TCRa_input_cdr3", "TCRb_input_cdr3", "TCRa_sequence_source", "TCRb_sequence_source",
    "tcrtransform_status",
]
PARSE_COLUMNS = [
    "tcr_species",
    "TRAV", "TRAJ", "TRAC",
    "CDR1a", "CDR2a", "CDR2.5a", "CDR3a",
    *_boundary_columns("a"),
    "TRAV_identity", "TRAJ_identity", "alpha_parse_status",
    "TRBV", "TRBJ", "TRBC",
    "CDR1b", "CDR2b", "CDR2.5b", "CDR3b",
    *_boundary_columns("b"),
    "TRBV_identity", "TRBJ_identity", "beta_parse_status",
    "tcrtransform_status",
]


@dataclass
class RuntimeReference:
    genes: dict[str, str]
    junctions: dict[str, JunctionMetadata]
    reverse: SpeciesReference


class Runtime:
    def __init__(self) -> None:
        self.loops = VGeneLoopReference.from_tsv()
        self.by_species: dict[str, RuntimeReference] = {}

    def get(self, species_value: object) -> tuple[str, RuntimeReference]:
        species = normalize_species(species_value)
        if species not in self.by_species:
            genes = load_gene_sequences(species)
            junctions = load_junction_reference(species)
            self.by_species[species] = RuntimeReference(
                genes=genes,
                junctions=junctions,
                reverse=SpeciesReference(species, genes, junctions),
            )
        return species, self.by_species[species]


def _resolve_column(
    header: list[str],
    explicit: str | None,
    aliases: Iterable[str],
    label: str,
) -> str:
    lower_map: dict[str, list[str]] = {}
    for column in header:
        lower_map.setdefault(column.lower(), []).append(column)
    if explicit:
        matches = lower_map.get(explicit.strip().lower(), [])
    else:
        matches = []
        for alias in aliases:
            matches.extend(lower_map.get(alias.lower(), []))
        matches = list(dict.fromkeys(matches))
    if len(matches) != 1:
        raise TCRTransformError(
            f"could not identify one {label} column; candidates were {matches or 'none'}"
        )
    return matches[0]


def _species_for_row(
    row: dict[str, str],
    fixed_species: str,
    species_column: str | None,
) -> str:
    return normalize_species(row[species_column] if species_column else fixed_species)


def _boundary_output(boundary: CDR3Boundary, suffix: str) -> dict[str, object]:
    base = f"CDR3{suffix}"
    return {
        f"{base}_v_supported": boundary.v_supported_prefix,
        f"{base}_non_vj": boundary.non_vj_region,
        f"{base}_j_supported": boundary.j_supported_suffix,
        f"{base}_keep_prefix": boundary.keep_prefix,
        f"{base}_keep_suffix": boundary.keep_suffix,
        f"{base}_non_vj_start": boundary.non_vj_start,
        f"{base}_non_vj_end": boundary.non_vj_end,
        f"{base}_non_vj_length": boundary.non_vj_length,
        f"{base}_fixed_mask": boundary.fixed_mask,
        f"{base}_v_gene_start": boundary.v_gene_start,
        f"{base}_v_gene_end": boundary.v_gene_end,
        f"{base}_v_offset_from_end": boundary.v_offset_from_end,
        f"{base}_j_gene_start": boundary.j_gene_start,
        f"{base}_j_gene_end": boundary.j_gene_end,
        f"{base}_j_offset_from_start": boundary.j_offset_from_start,
        f"{base}_v_short_match_length": boundary.v_short_match_length,
        f"{base}_j_short_match_length": boundary.j_short_match_length,
        f"{base}_v_j_support_overlap": boundary.v_j_support_overlap,
        f"{base}_safe_for_optimization": boundary.safe_for_optimization,
        f"{base}_boundary_status": boundary.status,
    }


def _region_output(
    *,
    species: str,
    runtime: Runtime,
    reference: RuntimeReference,
    trav: object,
    traj: object,
    cdr3a: object,
    trbv: object,
    trbj: object,
    cdr3b: object,
    min_overlap: int,
    max_overlap: int,
    v_terminal_window: int,
    j_start_window: int,
) -> tuple[dict[str, object], CDR3Boundary, CDR3Boundary]:
    alpha_loops = annotate_v_gene(
        species, "A", trav, loop_reference=runtime.loops
    )
    beta_loops = annotate_v_gene(
        species, "B", trbv, loop_reference=runtime.loops
    )
    alpha_boundary = define_cdr3_boundary(
        cdr3=cdr3a,
        v_gene=trav,
        j_gene=traj,
        chain="A",
        species=species,
        gene_sequences=reference.genes,
        junctions=reference.junctions,
        min_overlap=min_overlap,
        max_overlap=max_overlap,
        v_terminal_window=v_terminal_window,
        j_start_window=j_start_window,
    )
    beta_boundary = define_cdr3_boundary(
        cdr3=cdr3b,
        v_gene=trbv,
        j_gene=trbj,
        chain="B",
        species=species,
        gene_sequences=reference.genes,
        junctions=reference.junctions,
        min_overlap=min_overlap,
        max_overlap=max_overlap,
        v_terminal_window=v_terminal_window,
        j_start_window=j_start_window,
    )
    values: dict[str, object] = {
        "tcr_species": species,
        "CDR1a": alpha_loops.cdr1,
        "CDR2a": alpha_loops.cdr2,
        "CDR2.5a": alpha_loops.cdr2_5,
        "CDR3a": alpha_boundary.cdr3,
        "CDR1b": beta_loops.cdr1,
        "CDR2b": beta_loops.cdr2,
        "CDR2.5b": beta_loops.cdr2_5,
        "CDR3b": beta_boundary.cdr3,
    }
    values.update(_boundary_output(alpha_boundary, "a"))
    values.update(_boundary_output(beta_boundary, "b"))
    values["tcrtransform_status"] = (
        "ok"
        if alpha_boundary.safe_for_optimization and beta_boundary.safe_for_optimization
        else (
            "warning: unresolved CDR3 boundary; "
            f"alpha={alpha_boundary.status}; beta={beta_boundary.status}"
        )
    )
    return values, alpha_boundary, beta_boundary


def _annotation_processor(
    args: argparse.Namespace,
    header: list[str],
    runtime: Runtime,
    *,
    include_full_sequences: bool,
) -> Callable[[dict[str, str]], dict[str, object]]:
    if include_full_sequences and args.prefer_observed and not (args.observed_alpha_column or args.observed_beta_column):
        raise TCRTransformError('--prefer-observed requires an observed alpha or beta column')
    columns = {
        key: _resolve_column(
            header,
            getattr(args, f"{key}_column", None),
            aliases,
            key.upper(),
        )
        for key, aliases in BUILD_ALIASES.items()
    }
    species_column = (
        _resolve_column(header, args.species_column, (args.species_column,), "species")
        if args.species_column
        else None
    )

    def process(row: dict[str, str]) -> dict[str, object]:
        species, reference = runtime.get(
            _species_for_row(row, args.species, species_column)
        )
        def annotate(cdr3a, cdr3b):
            return _region_output(
                species=species,
                runtime=runtime,
                reference=reference,
                trav=row[columns["trav"]],
                traj=row[columns["traj"]],
                cdr3a=cdr3a,
                trbv=row[columns["trbv"]],
                trbj=row[columns["trbj"]],
                cdr3b=cdr3b,
                min_overlap=args.min_overlap,
                max_overlap=args.max_overlap,
                v_terminal_window=args.v_terminal_window,
                j_start_window=args.j_start_window,
            )
        values = {} if include_full_sequences else annotate(row[columns['cdr3a']], row[columns['cdr3b']])[0]
        if include_full_sequences:
            def evidence(letter):
                column = getattr(args, f'observed_{letter}_column')
                if column and column not in row:
                    raise TCRTransformError(f'observed sequence column missing: {column}')
                return dict(observed_sequence=row.get(column) if column else None,
                    observed_kind=args.observed_kind,
                    prefer_observed=args.prefer_observed and bool(column),
                    complete_cdr3_from_observed=args.complete_cdr3_from_observed,
                    correct_terminal_to_f=args.correct_terminal_to_f)
            alpha = reconstruct_tcr_chain(
                **evidence('alpha'),
                cdr3=row[columns["cdr3a"]],
                v_gene=row[columns["trav"]],
                j_gene=row[columns["traj"]],
                chain="A",
                species=species,
                gene_sequences=reference.genes,
                junctions=reference.junctions,
                min_overlap=args.min_overlap,
                max_overlap=args.max_overlap,
                v_terminal_window=args.v_terminal_window,
                j_start_window=args.j_start_window,
            )
            beta = reconstruct_tcr_chain(
                **evidence('beta'),
                cdr3=row[columns["cdr3b"]],
                v_gene=row[columns["trbv"]],
                j_gene=row[columns["trbj"]],
                chain="B",
                species=species,
                gene_sequences=reference.genes,
                junctions=reference.junctions,
                min_overlap=args.min_overlap,
                max_overlap=args.max_overlap,
                v_terminal_window=args.v_terminal_window,
                j_start_window=args.j_start_window,
            )
            observed_used = alpha.sequence_source.startswith('observed') or beta.sequence_source.startswith('observed')
            if observed_used:
                # Do not label observed mutations as germline-supported loops.
                values = {column: '' for column in REGION_COLUMNS}
                values.update(tcr_species=species, CDR3a=alpha.cdr3, CDR3b=beta.cdr3)
            else:
                values = annotate(alpha.cdr3, beta.cdr3)[0]
            values.update(
                {
                    "TCRa": alpha.sequence,
                    "TCRb": beta.sequence,
                    "TCRa_reference_resolution": alpha.reference_resolution,
                    "TCRb_reference_resolution": beta.reference_resolution,
                    "TCRa_reference_warnings": alpha.reference_warnings,
                    "TCRb_reference_warnings": beta.reference_warnings,
                    "TCRa_productivity_status": alpha.productivity_status,
                    "TCRb_productivity_status": beta.productivity_status,
                    "TCRa_input_cdr3": alpha.input_cdr3,
                    "TCRb_input_cdr3": beta.input_cdr3,
                    "TCRa_sequence_source": alpha.sequence_source,
                    "TCRb_sequence_source": beta.sequence_source,
                    "tcrtransform_status": (
                        "warning: observed sequence preserved; germline region annotations omitted" if observed_used
                        else "warning: CDR3 corrected; inspect input and resolved sequences"
                        if alpha.input_cdr3 != alpha.cdr3 or beta.input_cdr3 != beta.cdr3 else "ok"
                    ),
                }
            )
        return values

    return process


def _parse_processor(
    args: argparse.Namespace,
    header: list[str],
    runtime: Runtime,
) -> Callable[[dict[str, str]], dict[str, object]]:
    alpha_column = _resolve_column(
        header, args.alpha_column, ALPHA_SEQUENCE_ALIASES, "TCR alpha"
    )
    beta_column = _resolve_column(
        header, args.beta_column, BETA_SEQUENCE_ALIASES, "TCR beta"
    )
    species_column = (
        _resolve_column(header, args.species_column, (args.species_column,), "species")
        if args.species_column
        else None
    )

    def process(row: dict[str, str]) -> dict[str, object]:
        species, reference = runtime.get(
            _species_for_row(row, args.species, species_column)
        )
        alpha = parse_chain(
            row[alpha_column],
            "A",
            reference.reverse,
            runtime.loops,
            max_leader_length=args.max_leader_length,
            min_v_identity=args.min_v_identity,
            min_j_identity=args.min_j_identity,
        )
        beta = parse_chain(
            row[beta_column],
            "B",
            reference.reverse,
            runtime.loops,
            max_leader_length=args.max_leader_length,
            min_v_identity=args.min_v_identity,
            min_j_identity=args.min_j_identity,
        )
        values: dict[str, object] = {
            "tcr_species": species,
            "TRAV": alpha.v_genes,
            "TRAJ": alpha.j_genes,
            "TRAC": alpha.constant_genes,
            "CDR1a": alpha.cdr1,
            "CDR2a": alpha.cdr2,
            "CDR2.5a": alpha.cdr2_5,
            "CDR3a": alpha.cdr3,
            "TRAV_identity": f"{alpha.v_identity:.6f}",
            "TRAJ_identity": f"{alpha.j_identity:.6f}",
            "alpha_parse_status": alpha.status,
            "TRBV": beta.v_genes,
            "TRBJ": beta.j_genes,
            "TRBC": beta.constant_genes,
            "CDR1b": beta.cdr1,
            "CDR2b": beta.cdr2,
            "CDR2.5b": beta.cdr2_5,
            "CDR3b": beta.cdr3,
            "TRBV_identity": f"{beta.v_identity:.6f}",
            "TRBJ_identity": f"{beta.j_identity:.6f}",
            "beta_parse_status": beta.status,
            "tcrtransform_status": "ok",
        }
        values.update(_boundary_output(alpha.cdr3_boundary, "a"))
        values.update(_boundary_output(beta.cdr3_boundary, "b"))
        return values

    return process


def _output_path(args: argparse.Namespace) -> Path:
    input_path = Path(args.input_csv).expanduser().resolve()
    if args.output:
        return Path(args.output).expanduser().resolve()
    return input_path.with_name(f"{input_path.stem}_{args.command}.csv")


def _run_csv(
    args: argparse.Namespace,
    added_columns: list[str],
    processor_factory: Callable[
        [argparse.Namespace, list[str], Runtime],
        Callable[[dict[str, str]], dict[str, object]],
    ],
) -> None:
    input_path = Path(args.input_csv).expanduser().resolve()
    output_path = _output_path(args)
    if output_path == input_path:
        raise TCRTransformError("output must differ from input")
    if output_path.exists() and not args.overwrite:
        raise TCRTransformError(f"{output_path} exists; use --overwrite")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    runtime = Runtime()

    temporary_handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=output_path.parent,
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        delete=False,
        buffering=1024 * 1024,
    )
    temporary_path = Path(temporary_handle.name)
    rows = 0
    errors = 0
    try:
        with input_path.open(
            "r", encoding="utf-8-sig", newline="", buffering=1024 * 1024
        ) as input_handle, temporary_handle:
            reader = csv.DictReader(input_handle)
            header = [str(column).strip() for column in (reader.fieldnames or [])]
            if not header or any(not column for column in header):
                raise TCRTransformError("input has no valid CSV header")
            if len(set(header)) != len(header):
                raise TCRTransformError("input has duplicate CSV column names")
            reader.fieldnames = header
            collisions = sorted(set(header) & set(added_columns))
            if collisions and not args.overwrite_columns:
                raise TCRTransformError(
                    f"output columns already exist: {collisions}; use --overwrite-columns"
                )
            output_header = header + [column for column in added_columns if column not in header]
            processor = processor_factory(args, header, runtime)
            writer = csv.DictWriter(
                temporary_handle,
                fieldnames=output_header,
                extrasaction="raise",
                lineterminator="\n",
            )
            writer.writeheader()
            for csv_row, row in enumerate(reader, start=2):
                if None in row:
                    raise TCRTransformError(
                        f"CSV row {csv_row} has more values than the header"
                    )
                rows += 1
                try:
                    values = processor(row)
                except (OSError, ValueError) as exc:
                    errors += 1
                    if args.on_error == "fail":
                        raise TCRTransformError(f"CSV row {csv_row}: {exc}") from exc
                    values = {column: "" for column in added_columns}
                    values["tcrtransform_status"] = f"error: {exc}"
                row.update(values)
                writer.writerow(row)
        if rows == 0:
            raise TCRTransformError("input contains no data rows")
        if input_path.exists():
            os.chmod(temporary_path, stat.S_IMODE(input_path.stat().st_mode))
        os.replace(temporary_path, output_path)
    except Exception:
        temporary_handle.close()
        temporary_path.unlink(missing_ok=True)
        raise
    print(f"Rows: {rows:,}")
    print(f"Errors retained: {errors:,}")
    print(f"Output: {output_path}")


def _add_common_csv_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input_csv")
    parser.add_argument("-o", "--output")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--overwrite-columns", action="store_true")
    parser.add_argument("--species", default="human")
    parser.add_argument("--species-column")
    parser.add_argument("--on-error", choices=("fail", "keep"), default="keep")


def _add_gene_input_options(parser: argparse.ArgumentParser) -> None:
    for key in BUILD_ALIASES:
        parser.add_argument(f"--{key}-column")
    parser.add_argument("--min-overlap", type=int, default=3)
    parser.add_argument("--max-overlap", type=int, default=30)
    parser.add_argument("--v-terminal-window", type=int, default=8)
    parser.add_argument("--j-start-window", type=int, default=10)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tcrtransform",
        description="Reconstruct and annotate human or mouse alpha/beta TCR proteins.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser(
        "build", help="Build full TCRa/TCRb and append all CDR annotations."
    )
    _add_common_csv_options(build)
    _add_gene_input_options(build)
    build.add_argument('--observed-alpha-column')
    build.add_argument('--observed-beta-column')
    build.add_argument('--observed-kind', choices=('full', 'variable'), default='full')
    build.add_argument('--prefer-observed', action='store_true')
    build.add_argument('--complete-cdr3-from-observed', action='store_true')
    build.add_argument('--correct-terminal-to-f', action='store_true')

    annotate = subparsers.add_parser(
        "annotate", help="Append CDR1/2/2.5 and CDR3 V/non-VJ/J regions."
    )
    _add_common_csv_options(annotate)
    _add_gene_input_options(annotate)

    parse = subparsers.add_parser(
        "parse", help="Infer genes and CDR regions from full TCRa/TCRb sequences."
    )
    _add_common_csv_options(parse)
    parse.add_argument("--alpha-column")
    parse.add_argument("--beta-column")
    parse.add_argument("--max-leader-length", type=int, default=60)
    parse.add_argument("--min-v-identity", type=float, default=0.95)
    parse.add_argument("--min-j-identity", type=float, default=0.90)
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "build":
            _run_csv(
                args,
                BUILD_COLUMNS,
                lambda namespace, header, runtime: _annotation_processor(
                    namespace, header, runtime, include_full_sequences=True
                ),
            )
        elif args.command == "annotate":
            columns = [*REGION_COLUMNS, "tcrtransform_status"]
            _run_csv(
                args,
                columns,
                lambda namespace, header, runtime: _annotation_processor(
                    namespace, header, runtime, include_full_sequences=False
                ),
            )
        else:
            _run_csv(args, PARSE_COLUMNS, _parse_processor)
    except (OSError, csv.Error, TCRTransformError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()

