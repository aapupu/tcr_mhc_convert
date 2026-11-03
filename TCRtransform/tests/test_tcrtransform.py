from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from tcrtransform.cli import main as cli_main
from tcrtransform.junctions import (
    JunctionMetadata,
    load_junction_reference,
    match_j_terminal,
    match_v_terminal,
)
from tcrtransform.references import VGeneLoopReference, load_gene_sequences
from tcrtransform.regions import define_cdr3_boundary
from tcrtransform.reverse import SpeciesReference, parse_chain
from tcrtransform.transform import reconstruct_tcr_chain


class ReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.loops = VGeneLoopReference.from_tsv()

    def test_human_and_mouse_v_loops(self) -> None:
        human = self.loops.resolve("human", "A", "TRAV1-1")
        mouse = self.loops.resolve("mouse", "B", "TRBV12-1")
        self.assertEqual(
            (human.cdr1, human.cdr2, human.cdr2_5),
            ("TSGFYG", "NALDGL", "SRSDSY"),
        )
        self.assertEqual(
            (mouse.cdr1, mouse.cdr2, mouse.cdr2_5),
            ("SGHSN", "HYEKVE", "FDDYH"),
        )

    def test_unknown_explicit_allele_is_not_silently_replaced(self) -> None:
        with self.assertRaises(ValueError):
            self.loops.resolve("human", "A", "TRAV1-1*99")


class BoundaryTests(unittest.TestCase):
    def test_guarded_short_matches(self) -> None:
        v_meta = JunctionMetadata("human", "TRAVX", "TRAVX*01", "V", "C", 4)
        j_meta = JunctionMetadata("human", "TRAJX", "TRAJX*01", "J", "F", 3)
        self.assertIsNotNone(
            match_v_terminal(
                "AAAACPPE", "CQQF", min_overlap=1, max_overlap=2,
                terminal_window=4, junction=v_meta,
            )
        )
        self.assertIsNotNone(
            match_j_terminal(
                "CQQF", "AAAFGG", min_overlap=1, max_overlap=2,
                start_window=5, junction=j_meta,
            )
        )
        wrong_v = JunctionMetadata("human", "TRAVX", "TRAVX*01", "V", "C", 3)
        self.assertIsNone(
            match_v_terminal(
                "AAAACPPE", "CQQF", min_overlap=1, max_overlap=2,
                terminal_window=4, junction=wrong_v,
            )
        )

    def test_strict_threshold_reports_but_does_not_accept_short_hits(self) -> None:
        genes = {"TRAVX": "AAAACPPE", "TRAJX": "AAAFGG"}
        junctions = {
            "TRAVX": JunctionMetadata(
                "human", "TRAVX", "TRAVX*01", "V", "C", 4
            ),
            "TRAJX": JunctionMetadata(
                "human", "TRAJX", "TRAJX*01", "J", "F", 3
            ),
        }
        strict = define_cdr3_boundary(
            cdr3="CQQF", v_gene="TRAVX", j_gene="TRAJX", chain="A",
            gene_sequences=genes, junctions=junctions, min_overlap=3,
        )
        self.assertEqual(strict.keep_prefix, 0)
        self.assertEqual(strict.keep_suffix, 0)
        self.assertEqual(strict.v_short_match_length, 1)
        self.assertEqual(strict.j_short_match_length, 1)
        self.assertFalse(strict.safe_for_optimization)
        guarded = define_cdr3_boundary(
            cdr3="CQQF", v_gene="TRAVX", j_gene="TRAJX", chain="A",
            gene_sequences=genes, junctions=junctions, min_overlap=1,
        )
        self.assertEqual(guarded.v_supported_prefix, "C")
        self.assertEqual(guarded.j_supported_suffix, "F")
        self.assertEqual(guarded.non_vj_region, "QQ")
        self.assertTrue(guarded.safe_for_optimization)

    def test_fully_vj_supported_cdr3_has_zero_editable_residues_for_both_chains(self) -> None:
        cdr3 = "CASSYF"
        for chain, v_gene, j_gene in (
            ("A", "TRAVX", "TRAJX"),
            ("B", "TRBVX", "TRBJX"),
        ):
            genes = {v_gene: "AAAAACASQQ", j_gene: "SYFGGG"}
            junctions = {
                v_gene: JunctionMetadata(
                    "human", v_gene, f"{v_gene}*01", "V", "C", 5
                ),
                j_gene: JunctionMetadata(
                    "human", j_gene, f"{j_gene}*01", "J", "F", 2
                ),
            }
            boundary = define_cdr3_boundary(
                cdr3=cdr3,
                v_gene=v_gene,
                j_gene=j_gene,
                chain=chain,
                gene_sequences=genes,
                junctions=junctions,
                min_overlap=3,
            )
            self.assertEqual(boundary.keep_prefix, 3)
            self.assertEqual(boundary.keep_suffix, 3)
            self.assertEqual(boundary.non_vj_region, "")
            self.assertEqual(boundary.non_vj_length, 0)
            self.assertEqual(boundary.fixed_mask, "111111")
            self.assertEqual(boundary.status, "fully_vj_supported")


class RoundTripTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.species = "mouse"
        cls.genes = load_gene_sequences(cls.species)
        cls.junctions = load_junction_reference(cls.species)
        cls.reference = SpeciesReference(cls.species, cls.genes, cls.junctions)
        cls.loops = VGeneLoopReference.from_tsv()

    def test_mouse_forward_reverse_round_trip(self) -> None:
        alpha = reconstruct_tcr_chain(
            cdr3="CAASDNYQLIW",
            v_gene="TRAV14D-1",
            j_gene="TRAJ33",
            chain="A",
            species=self.species,
            gene_sequences=self.genes,
            junctions=self.junctions,
        )
        beta = reconstruct_tcr_chain(
            cdr3="CASSRANYEQYF",
            v_gene="TRBV12-1",
            j_gene="TRBJ2-7",
            chain="B",
            species=self.species,
            gene_sequences=self.genes,
            junctions=self.junctions,
        )
        self.assertEqual(alpha.sequence.count(alpha.cdr3), 1)
        self.assertEqual(beta.sequence.count(beta.cdr3), 1)
        self.assertIn("KLIIKPDIQNPE", alpha.sequence)

        parsed_a = parse_chain(
            alpha.sequence, "A", self.reference, self.loops
        )
        parsed_b = parse_chain(
            beta.sequence, "B", self.reference, self.loops
        )
        self.assertIn("TRAV14D-1", parsed_a.v_genes.split("|"))
        self.assertEqual(parsed_a.j_genes, "TRAJ33")
        self.assertEqual(parsed_a.cdr3, alpha.cdr3)
        self.assertIn("TRBV12-1", parsed_b.v_genes.split("|"))
        self.assertEqual(parsed_b.j_genes, "TRBJ2-7")
        self.assertEqual(parsed_b.cdr3, beta.cdr3)
        self.assertTrue(parsed_a.cdr3_boundary.safe_for_optimization)

    def test_build_cli_writes_full_sequences_and_regions(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "input.csv"
            output_path = Path(directory) / "output.csv"
            with input_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    ["species", "TRAV", "TRAJ", "CDR3A", "TRBV", "TRBJ", "CDR3B"]
                )
                writer.writerow(
                    [
                        "mouse", "TRAV14D-1*01", "TRAJ33*01", "CAASDNYQLIW",
                        "TRBV12-1*01", "TRBJ2-7*01", "CASSRANYEQYF",
                    ]
                )
            cli_main(
                [
                    "build", str(input_path), "-o", str(output_path),
                    "--species-column", "species",
                ]
            )
            with output_path.open("r", encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(row["tcrtransform_status"], "ok")
            self.assertIn(row["CDR3a"], row["TCRa"])
            self.assertIn(row["CDR3b"], row["TCRb"])
            self.assertIn("gene_level_fallback", row["TCRa_reference_resolution"])
            self.assertEqual(
                len(row["CDR3a_fixed_mask"]), len(row["CDR3a"])
            )


if __name__ == "__main__":
    unittest.main()

