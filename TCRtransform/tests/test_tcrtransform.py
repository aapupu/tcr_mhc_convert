from __future__ import annotations

import csv
import runpy
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
from tcrtransform.references import VGeneLoopReference, load_gene_sequences, get_gene_sequence
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

    def test_unknown_explicit_allele_falls_back_to_gene(self) -> None:
        self.assertEqual(self.loops.resolve("human", "A", "TRAV1-1*99"),
                         self.loops.resolve("human", "A", "TRAV1-1"))
        with self.assertRaises(ValueError):
            self.loops.resolve("human", "A", "TRAV999*01")

    def test_traj24_alleles_and_fallback(self) -> None:
        genes = load_gene_sequences("human")
        expected = {"01": "TTDSWGKFEFGAGTQVVVTP", "02": "TTDSWGKLQFGAGTQVVVTP",
                    "03": "TTDSWGKFQFGAGTQVVVTP"}
        for allele, sequence in expected.items():
            self.assertEqual(get_gene_sequence(genes, "TRAJ24*" + allele), sequence)
            cdr3 = "CAVR" + sequence[:sequence.rfind("FGAG") + 1]
            result = reconstruct_tcr_chain(cdr3=cdr3, v_gene="TRAV1-1*01",
                j_gene="TRAJ24*" + allele, chain="A", species="human")
            self.assertIn("J:exact(TRAJ24*" + allele + ")", result.reference_resolution)
            self.assertIn(sequence + "D", result.sequence)
        self.assertEqual(get_gene_sequence(genes, "TRAJ24"), expected["02"])
        self.assertEqual(get_gene_sequence(genes, "TRAJ24*99"), expected["02"])
        result = reconstruct_tcr_chain(cdr3="CAVRTTDSWGKLQF", v_gene="TRAV1-1*99",
            j_gene="TRAJ24*99", chain="A", species="human")
        self.assertIn("V:gene_level_fallback(TRAV1-1)", result.reference_resolution)
        self.assertIn("J:gene_level_fallback(TRAJ24)", result.reference_resolution)

    def test_mouse_corrections_aliases_and_nonfunctional_guards(self) -> None:
        genes = load_gene_sequences("mouse")
        self.assertTrue(genes["TRAV13D-1"].endswith("CAME"))
        self.assertIn("TVWEGGTTVLTCS", genes["TRAV14-1"])
        self.assertIn("TVWEGETAILNCS", genes["TRAV14D-1"])
        mouse_loops = self.loops.resolve("mouse", "A", "TRAV14-1*01")
        self.assertEqual(mouse_loops.cdr2, "ILSVSDK")
        self.assertEqual(genes["TRAJ49"], "NTGYQNFYFGKGTSLTVIP")
        self.assertEqual(get_gene_sequence(genes, "TRAV4-4-DV10*01"),
                         get_gene_sequence(genes, "TRAV4-4/DV10*01"))
        human = load_gene_sequences("human")
        for gene in ["TRAJ58", "TRAJ60", "TRBV12-1*01"]:
            self.assertIn('*', get_gene_sequence(human, gene))
        with self.assertRaises(ValueError):
            get_gene_sequence(genes, "TRBV12-2+TRBV13-2")
        # An incomplete functional allele is unavailable for full reconstruction.
        self.assertEqual(get_gene_sequence(genes, "TRAV8-1*02"), genes["TRAV8-1"])
        self.assertEqual(self.loops.resolve("mouse", "A", "TRAV8-1*02"),
                         self.loops.resolve("mouse", "A", "TRAV8-1"))

    def test_custom_reference_keeps_exact_then_gene_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "genes.csv"
            path.write_text("gene_name,amino_acid\nTRAVX,ACDE\nTRAVX*02,ACDF\n")
            genes = load_gene_sequences("human", path=path)
            self.assertEqual(get_gene_sequence(genes, "TRAVX*02"), "ACDF")
            self.assertEqual(get_gene_sequence(genes, "TRAVX*03"), "ACDE")
            self.assertEqual(get_gene_sequence(genes, "TRAVX"), "ACDE")
            self.assertNotIn("TRAV1-1", genes)

    def test_overlapping_j_motifs_use_the_conserved_downstream_anchor(self) -> None:
        junctions = load_junction_reference("human")
        self.assertEqual(junctions["TRAJ14*01"].cdr3_anchor_index, 6)  # IYSTFIFGSG...
        self.assertEqual(junctions["TRAJ21*01"].cdr3_anchor_index, 7)  # YNFNKFYFGSG...

    def test_rearrangement_removes_terminal_stops_before_validation(self) -> None:
        result = reconstruct_tcr_chain(species='mouse', chain='B', v_gene='TRBV4*01',
            j_gene='TRBJ1-5*02', cdr3='CASSQPLF')
        self.assertNotIn('*', result.sequence)
        self.assertIn('germline_stop_removed', result.reference_warnings)
        self.assertIn('V:exact(TRBV4*01)', result.reference_resolution)
        self.assertIn('J:exact(TRBJ1-5*02)', result.reference_resolution)
        self.assertEqual(result.productivity_status, 'not_assessed')
        result = reconstruct_tcr_chain(species='human', chain='A', v_gene='TRAV1-1*01',
            j_gene='TRAJ58*01', cdr3='CAVRETSGSRLTF')
        self.assertIn('germline_ORF', result.reference_warnings)
        self.assertNotIn('*', result.sequence)

    def test_retained_stops_fail_reconstruction_without_classifying_sample(self) -> None:
        for v, j, cdr3 in [('TRBV12-1*01', 'TRBJ2-7*01', 'CASSFNEQYF')]:
            with self.assertRaisesRegex(ValueError, 'retained_stop_in_reconstruction'):
                reconstruct_tcr_chain(species='human', chain='B' if v.startswith('TRB') else 'A',
                    v_gene=v, j_gene=j, cdr3=cdr3)
        with self.assertRaisesRegex(ValueError, 'retained_stop_in_reconstruction'):
            reconstruct_tcr_chain(species='human', chain='B', v_gene='TRBVX',
                j_gene='TRBJX', cdr3='CASSQYF',
                gene_sequences={'TRBVX':'AAAACASS', 'TRBJX':'QYFGPG*', 'TRBC1':'AAA'},
                junctions={'TRBVX':JunctionMetadata('human','TRBVX','TRBVX*01','V','C',4),
                           'TRBJX':JunctionMetadata('human','TRBJX','TRBJX*01','J','F',2)})

    def test_orf_label_alone_does_not_reject_reconstruction(self) -> None:
        result = reconstruct_tcr_chain(species='human', chain='A', v_gene='TRAV1-1*01',
            j_gene='TRAJ1*01', cdr3='CAVRYESITSQLQF')
        self.assertIn('germline_ORF', result.reference_warnings)
        self.assertEqual(result.productivity_status, 'not_assessed')

    def test_noncanonical_j_boundaries_reconstruct_and_parse(self) -> None:
        cases = [('mouse', 'TRAJ44', 'VTGSGGKLTLGAGTRLQVNL', 9, 'D'),
                 ('mouse', 'TRAJ7', 'DYSNNRLTLGKGTQVVVLP', 8, 'N'),
                 ('human', 'TRAJ61', 'YRVNRKLTFGANTRGIMKL', 8, 'N')]
        for species, gene, j_sequence, index, bridge in cases:
            genes = load_gene_sequences(species)
            junctions = load_junction_reference(species)
            v_gene = 'TRAV1*01' if species == 'mouse' else 'TRAV1-1*01'
            v_index = junctions[v_gene].cdr3_anchor_index
            self.assertEqual(junctions[gene].cdr3_anchor_index, index)
            self.assertEqual(junctions[gene + '*01'].cdr3_anchor, j_sequence[index])
            for suffix in [j_sequence[:index + 1], j_sequence[index - 2:index + 1]]:
                cdr3 = 'CAVR' + suffix
                with self.subTest(species=species, gene=gene, cdr3=cdr3):
                    result = reconstruct_tcr_chain(species=species, chain='A',
                        v_gene=v_gene, j_gene=gene + '*01', cdr3=cdr3)
                    expected = genes[v_gene][:v_index] + cdr3 + j_sequence[index + 1:] + bridge + genes['TRAC']
                    self.assertEqual(result.sequence, expected)
                    self.assertEqual(parse_chain(result.sequence, 'A', species=species).cdr3, cdr3)
                    self.assertIn('germline_ORF', result.reference_warnings)
            # Noncanonical L anchors also support guarded short matches.
            self.assertIsNotNone(match_j_terminal('CQQ' + j_sequence[index-1:index+1],
                j_sequence, min_overlap=2, max_overlap=2, junction=junctions[gene]))

    def test_j_anchor_generator_validates_annotated_exceptions(self) -> None:
        builder = Path(__file__).resolve().parents[1] / 'tools/build_imgt_reference.py'
        resolve = runpy.run_path(str(builder))['j_anchor']
        self.assertEqual(resolve('mouse', 'TRAJ44*01', 'VTGSGGKLTLGAGTRLQVNL'), 9)
        self.assertEqual(resolve('mouse', 'TRAJ7*01', 'DYSNNRLTLGKGTQVVVLP'), 8)
        self.assertEqual(resolve('human', 'TRAJ61*01', 'YRVNRKLTFGANTRGIMKL'), 8)
        self.assertEqual(resolve('human', 'UNKNOWN*01', 'YRVNRKLTFGANTRGIMKL'), -1)
        with self.assertRaisesRegex(ValueError, 'requires review'):
            resolve('human', 'TRAJ61*01', 'ARVNRKLTFGANTRGIMKL')


class BoundaryTests(unittest.TestCase):
    def test_allelic_v_anchor_and_single_f_after_anchor_filter(self) -> None:
        r = reconstruct_tcr_chain(species='human', chain='B', v_gene='TRBV12-3*01',
            j_gene='TRBJ2-7*02', cdr3='CASSRLGGTYEQYV')
        self.assertIn('J:exact(TRBJ2-7*02)', r.reference_resolution)
        self.assertEqual(r.j_match.sequence[-1], 'V')
        r = reconstruct_tcr_chain(species='human', chain='B', v_gene='TRBV12-3*01',
            j_gene='TRBJ2-1*01', cdr3='CASSQGSSGRRAF', min_overlap=1)
        self.assertEqual((r.j_match.length, r.j_match.gene_start), (1, 6))
        with self.assertRaises(ValueError):
            reconstruct_tcr_chain(species='human', chain='B', v_gene='TRBV12-3*01',
                j_gene='TRBJ2-1*01', cdr3='CASSQGSSGRRAF')

    def test_terminal_correction_is_explicit_and_context_limited(self) -> None:
        alpha = reconstruct_tcr_chain(species='human', chain='A', v_gene='TRAV1-1*01',
            j_gene='TRAJ6', cdr3='CASGGSYIPTY', min_overlap=1, correct_terminal_to_f=True)
        self.assertEqual(alpha.cdr3, 'CASGGSYIPTF')
        args = dict(species='human',chain='B',v_gene='TRBV12-3*01',j_gene='TRBJ2-7*01',cdr3='CASSLAGSKLYEQYC')
        with self.assertRaises(ValueError):
            reconstruct_tcr_chain(**args)
        r = reconstruct_tcr_chain(**args, correct_terminal_to_f=True)
        self.assertEqual(r.input_cdr3, args['cdr3'])
        self.assertEqual(r.cdr3, 'CASSLAGSKLYEQYF')
        with self.assertRaises(ValueError):
            reconstruct_tcr_chain(**dict(args,cdr3='CASSQQQQC'), correct_terminal_to_f=True)
        variable = r.sequence[:-len(load_gene_sequences('human')['TRBC2'])]
        pos = variable.index(r.cdr3) + len(r.cdr3) - 1
        observed = variable[:pos] + 'C' + variable[pos+1:]
        with self.assertRaisesRegex(ValueError, 'conflicts with observed'):
            reconstruct_tcr_chain(**args, correct_terminal_to_f=True,
                observed_sequence=observed, observed_kind='variable')
        kept = reconstruct_tcr_chain(**args, prefer_observed=True,
            observed_sequence=observed, observed_kind='variable')
        self.assertIn(args['cdr3'], kept.sequence)

    def test_missing_c_requires_observed_boundary_not_blind_prefix(self) -> None:
        genes = load_gene_sequences('human'); meta = load_junction_reference('human')
        v='TRAV1-1*01'; j='TRAJ17*01'; cdr3='CPAAGNKLTF'
        observed = genes[v][:meta[v].cdr3_anchor_index] + cdr3 + genes[j][meta[j].cdr3_anchor_index+1:]
        args=dict(species='human',chain='A',v_gene=v,j_gene=j,cdr3='PAAGNKLTF')
        with self.assertRaises(ValueError):
            reconstruct_tcr_chain(**args, complete_cdr3_from_observed=True)
        r = reconstruct_tcr_chain(**args, observed_sequence=observed,observed_kind='variable',
            complete_cdr3_from_observed=True,prefer_observed=True)
        self.assertEqual(r.cdr3,cdr3)
        self.assertTrue(r.sequence.startswith(observed))
        with self.assertRaises(ValueError):
            reconstruct_tcr_chain(**dict(args,cdr3='NTDKLIF'), observed_sequence=observed,
                observed_kind='variable',complete_cdr3_from_observed=True,prefer_observed=True)

    def test_observed_r_preserved_where_germline_has_stop(self) -> None:
        genes=load_gene_sequences('human'); meta=load_junction_reference('human')
        v='TRBV12-2*01'; j='TRBJ2-7*01'; cdr3='CASRLNEQYF'
        observed=genes[v][:meta[v].cdr3_anchor_index].replace('*','R')+cdr3+genes[j][meta[j].cdr3_anchor_index+1:]+genes['TRBC2']
        args=dict(species='human',chain='B',v_gene=v,j_gene=j,cdr3=cdr3)
        with self.assertRaisesRegex(ValueError,'retained_stop'):
            reconstruct_tcr_chain(**args)
        result=reconstruct_tcr_chain(**args,observed_sequence=observed,prefer_observed=True)
        self.assertEqual(result.sequence,observed)
        self.assertIn('52:*->R',result.reference_warnings)
        self.assertEqual(result.productivity_status,'not_assessed')
        self.assertIsNone(result.v_match)
        with self.assertRaises(ValueError):
            reconstruct_tcr_chain(**args,observed_sequence='A'*len(observed),prefer_observed=True)

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
    def test_cli_preserves_observed_stop_replacement_and_original_columns(self) -> None:
        genes=load_gene_sequences('human'); meta=load_junction_reference('human')
        v='TRBV12-2*01'; j='TRBJ2-7*01'; cdr3='CASRLNEQYF'
        observed=genes[v][:meta[v].cdr3_anchor_index].replace('*','R')+cdr3+genes[j][meta[j].cdr3_anchor_index+1:]+genes['TRBC2']
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'input.csv'; out=Path(directory)/'out.csv'
            with path.open('w',newline='') as handle:
                writer=csv.writer(handle)
                writer.writerow(['TRAV','TRAJ','CDR3A','TRBV','TRBJ','CDR3B','OTSbeta','productive'])
                for sequence in [observed, 'A'*len(observed)]:
                    writer.writerow(['TRAV1-1','TRAJ33','CAVRDSNYQLIW',v,j,cdr3,sequence,'True'])
            cli_main(['build',str(path),'-o',str(out),'--observed-beta-column','OTSbeta','--prefer-observed'])
            with out.open() as handle:
                rows=list(csv.DictReader(handle))
            self.assertEqual(rows[0]['TCRb'],observed)
            self.assertEqual(rows[0]['TCRb_sequence_source'],'observed_full')
            self.assertEqual(rows[0]['TCRb_input_cdr3'],cdr3)
            self.assertEqual(rows[0]['CDR2b'],'')
            self.assertIn('52:*->R',rows[0]['TCRb_reference_warnings'])
            self.assertEqual([r['productive'] for r in rows],['True','True'])
            self.assertEqual(len(rows),2)
            self.assertTrue(rows[1]['tcrtransform_status'].startswith('error:'))

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
                    ["species", "TRAV", "TRAJ", "CDR3A", "TRBV", "TRBJ", "CDR3B", "productive"]
                )
                writer.writerow(
                    [
                        "mouse", "TRAV14D-1*01", "TRAJ33*01", "CAASDNYQLIW",
                        "TRBV4*01", "TRBJ1-5*02", "CASSQPLF", "True",
                    ]
                )
                writer.writerow(['mouse', 'TRAV999*01', 'TRAJ33*01', 'CAASDNYQLIW',
                                 'TRBV4*01', 'TRBJ1-5*02', 'CASSQPLF', 'True'])
            cli_main(
                [
                    "build", str(input_path), "-o", str(output_path),
                    "--species-column", "species",
                ]
            )
            with output_path.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
                row = rows[0]
            self.assertEqual(len(rows), 2)
            self.assertEqual([r['productive'] for r in rows], ['True', 'True'])
            self.assertTrue(rows[1]['tcrtransform_status'].startswith('error:'))
            self.assertEqual(row['TCRb_productivity_status'], 'not_assessed')
            self.assertIn('germline_stop_removed', row['TCRb_reference_warnings'])
            self.assertEqual(row["tcrtransform_status"], "ok")
            self.assertIn(row["CDR3a"], row["TCRa"])
            self.assertIn(row["CDR3b"], row["TCRb"])
            self.assertIn("V:exact(TRAV14D-1*01)", row["TCRa_reference_resolution"])
            self.assertIn("J:exact(TRAJ33*01)", row["TCRa_reference_resolution"])
            self.assertEqual(
                len(row["CDR3a_fixed_mask"]), len(row["CDR3a"])
            )


if __name__ == "__main__":
    unittest.main()

