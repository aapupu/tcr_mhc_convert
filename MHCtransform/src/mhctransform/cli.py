"""Command-line interface for exact MHC conversions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from .errors import MHCTransformError
from .registry import AF3_DEFINITION, MHCRegistry


def _sequence_argument(value: str) -> str:
    if value.startswith("@"):
        return Path(value[1:]).read_text(encoding="utf-8")
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mhctransform",
        description=(
            "Exact human/mouse MHC precursor, AF3 extracellular-chain and "
            "34-aa pseudosequence lookup"
        ),
    )
    parser.add_argument("--data-dir", help="use an unpacked registry directory")
    subparsers = parser.add_subparsers(dest="command", required=True)

    allele_full = subparsers.add_parser(
        "allele-to-full", help="chain allele/name -> precursor protein candidate(s)"
    )
    allele_full.add_argument("name")
    allele_full.add_argument("--class", dest="mhc_class", choices=("I", "II"))
    allele_full.add_argument("--species", choices=("human", "mouse"))
    allele_full.add_argument("--include-partial", action="store_true")
    allele_full.add_argument("--allow-two-field-projection", action="store_true")

    full_allele = subparsers.add_parser(
        "full-to-allele", help="precursor protein -> all exact chain allele/name candidates"
    )
    full_allele.add_argument("sequence", help="amino acids, or @path to a sequence file")
    full_allele.add_argument("--class", dest="mhc_class", choices=("I", "II"))
    full_allele.add_argument("--species", choices=("human", "mouse"))
    full_allele.add_argument("--include-partial", action="store_true")

    allele_pseudo = subparsers.add_parser(
        "allele-to-pseudo", help="allele/molecule name -> 34-aa pseudosequence candidate(s)"
    )
    allele_pseudo.add_argument("name")
    allele_pseudo.add_argument("--class", dest="mhc_class", choices=("I", "II"), required=True)
    allele_pseudo.add_argument("--definition", required=True)
    allele_pseudo.add_argument("--species", choices=("human", "mouse"))

    pseudo_allele = subparsers.add_parser(
        "pseudo-to-allele", help="34-aa pseudosequence -> all exact allele/molecule candidates"
    )
    pseudo_allele.add_argument("sequence", help="34 amino acids, or @path")
    pseudo_allele.add_argument("--class", dest="mhc_class", choices=("I", "II"), required=True)
    pseudo_allele.add_argument("--definition", required=True)
    pseudo_allele.add_argument("--species", choices=("human", "mouse"))

    molecule = subparsers.add_parser(
        "molecule-to-chains", help="MHC molecule -> two full precursor chain candidates"
    )
    molecule.add_argument("name")
    molecule.add_argument("--class", dest="mhc_class", choices=("I", "II"), required=True)
    molecule.add_argument("--include-partial", action="store_true")
    molecule.add_argument("--allow-two-field-projection", action="store_true")

    allele_af3 = subparsers.add_parser(
        "allele-to-af3",
        help="chain allele/name -> mature soluble extracellular AF3 chain",
    )
    allele_af3.add_argument("name")
    allele_af3.add_argument("--class", dest="mhc_class", choices=("I", "II"))
    allele_af3.add_argument("--species", choices=("human", "mouse"))
    allele_af3.add_argument("--definition", default=AF3_DEFINITION)
    allele_af3.add_argument("--allow-two-field-projection", action="store_true")

    af3_allele = subparsers.add_parser(
        "af3-to-allele", help="AF3 extracellular chain -> all exact allele/name candidates"
    )
    af3_allele.add_argument("sequence", help="amino acids, or @path to a sequence file")
    af3_allele.add_argument("--class", dest="mhc_class", choices=("I", "II"))
    af3_allele.add_argument("--species", choices=("human", "mouse"))
    af3_allele.add_argument("--definition", default=AF3_DEFINITION)

    molecule_af3 = subparsers.add_parser(
        "molecule-to-af3", help="MHC molecule -> two separate AF3 protein chains"
    )
    molecule_af3.add_argument("name")
    molecule_af3.add_argument(
        "--class", dest="mhc_class", choices=("I", "II"), required=True
    )
    molecule_af3.add_argument("--definition", default=AF3_DEFINITION)
    molecule_af3.add_argument("--allow-two-field-projection", action="store_true")

    subparsers.add_parser("stats", help="show bundled registry counts and provenance")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        registry = MHCRegistry(args.data_dir)
        if args.command == "allele-to-full":
            result = registry.allele_to_sequence(
                args.name,
                mhc_class=args.mhc_class,
                species=args.species,
                include_partial=args.include_partial,
                allow_two_field_projection=args.allow_two_field_projection,
            ).as_dict()
        elif args.command == "full-to-allele":
            result = registry.sequence_to_alleles(
                _sequence_argument(args.sequence),
                mhc_class=args.mhc_class,
                species=args.species,
                include_partial=args.include_partial,
            ).as_dict()
        elif args.command == "allele-to-pseudo":
            result = registry.allele_to_pseudo(
                args.name,
                mhc_class=args.mhc_class,
                definition_id=args.definition,
                species=args.species,
            ).as_dict()
        elif args.command == "pseudo-to-allele":
            result = registry.pseudo_to_alleles(
                _sequence_argument(args.sequence),
                mhc_class=args.mhc_class,
                definition_id=args.definition,
                species=args.species,
            ).as_dict()
        elif args.command == "molecule-to-chains":
            result = registry.molecule_to_chains(
                args.name,
                mhc_class=args.mhc_class,
                include_partial=args.include_partial,
                allow_two_field_projection=args.allow_two_field_projection,
            ).as_dict()
        elif args.command == "allele-to-af3":
            result = registry.allele_to_af3(
                args.name,
                mhc_class=args.mhc_class,
                species=args.species,
                definition_id=args.definition,
                allow_two_field_projection=args.allow_two_field_projection,
            ).as_dict()
        elif args.command == "af3-to-allele":
            result = registry.af3_to_alleles(
                _sequence_argument(args.sequence),
                mhc_class=args.mhc_class,
                species=args.species,
                definition_id=args.definition,
            ).as_dict()
        elif args.command == "molecule-to-af3":
            result = registry.molecule_to_af3_chains(
                args.name,
                mhc_class=args.mhc_class,
                definition_id=args.definition,
                allow_two_field_projection=args.allow_two_field_projection,
            ).as_dict()
        else:
            result = registry.stats()
    except (MHCTransformError, OSError) as error:
        parser.error(str(error))
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    return 0

