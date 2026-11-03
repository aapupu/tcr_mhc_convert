import unittest

from mhctransform import MHCRegistry, MHCTransformError
from mhctransform.registry import (
    AF3_DEFINITION,
    PSEUDO_I_DEFINITION,
    PSEUDO_II_DEFINITION,
)


class RegistryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry = MHCRegistry()

    def test_class_i_golden_pseudosequences(self):
        expected = {
            "HLA-A*02:01": "YFAMYGEKVAHTHVDTLYVRYHYYTWAVLAYTWY",
            "H-2-Kb": "YVEYYREKAGNSFVDTLYIVSQYYTWAELAYTWY",
        }
        for name, pseudo in expected.items():
            with self.subTest(name=name):
                result = self.registry.allele_to_pseudo(
                    name, mhc_class="I", definition_id=PSEUDO_I_DEFINITION
                )
                self.assertEqual(result.status, "exact_unique")
                self.assertEqual(result.candidates[0].pseudo_sequence, pseudo)

    def test_class_ii_golden_pseudosequences(self):
        expected = {
            "HLA-DRB1*15:01": "QEFFIASGAAVDAIMWPRFDYFDIQAATYHVVFT",
            "H-2-IAb": "YSYFLASGGQVVHVLYFGYTYHDIRTETVHGPHT",
        }
        for name, pseudo in expected.items():
            with self.subTest(name=name):
                result = self.registry.allele_to_pseudo(
                    name, mhc_class="II", definition_id=PSEUDO_II_DEFINITION
                )
                self.assertEqual(result.status, "exact_unique")
                self.assertEqual(result.candidates[0].pseudo_sequence, pseudo)

    def test_higher_field_hla_name_falls_back_without_selecting_a_full_sequence(self):
        result = self.registry.allele_to_pseudo(
            "HLA-A*02:01:01:01",
            mhc_class="I",
            definition_id=PSEUDO_I_DEFINITION,
        )
        self.assertEqual(result.status, "projected_unique")
        self.assertTrue(result.notes)

    def test_pseudosequence_reverse_is_explicitly_many_to_many(self):
        pseudo = "YFAMYGEKVAHTHVDTLYVRYHYYTWAVLAYTWY"
        result = self.registry.pseudo_to_alleles(
            pseudo, mhc_class="I", definition_id=PSEUDO_I_DEFINITION
        )
        self.assertEqual(result.status, "exact_ambiguous")
        self.assertGreater(result.candidate_count, 1)
        self.assertIn("HLA-A*02:01", {item.canonical_name for item in result.candidates})

    def test_class_ii_x_is_valid_but_length_is_strict(self):
        result = self.registry.allele_to_pseudo(
            "H-2-IAq", mhc_class="II", definition_id=PSEUDO_II_DEFINITION
        )
        self.assertEqual(result.candidates[0].pseudo_sequence[:2], "XX")
        reverse = self.registry.pseudo_to_alleles(
            result.candidates[0].pseudo_sequence,
            mhc_class="II",
            definition_id=PSEUDO_II_DEFINITION,
        )
        self.assertIn("H-2-IAq", {item.canonical_name for item in reverse.candidates})
        with self.assertRaises(MHCTransformError):
            self.registry.pseudo_to_alleles(
                "A" * 33, mhc_class="II", definition_id=PSEUDO_II_DEFINITION
            )

    def test_human_and_mouse_precursor_full_sequences(self):
        human = self.registry.allele_to_sequence("HLA-A*02:01:01:01", mhc_class="I")
        self.assertEqual(human.status, "exact_unique")
        self.assertEqual(len(human.candidates[0].sequence), 365)
        self.assertTrue(human.candidates[0].sequence.startswith("M"))

        mouse = self.registry.allele_to_sequence("H-2-Kb", mhc_class="I")
        self.assertEqual(mouse.status, "exact_unique")
        self.assertEqual(len(mouse.candidates[0].sequence), 369)

    def test_two_field_full_sequence_ambiguity_is_not_hidden(self):
        result = self.registry.allele_to_sequence("HLA-C*03:03", mhc_class="I")
        self.assertEqual(result.status, "exact_ambiguous")
        self.assertEqual({len(item.sequence) for item in result.candidates}, {366, 372})

    def test_pre_2010_pseudo_names_map_to_current_nomenclature(self):
        modern = self.registry.allele_to_pseudo(
            "HLA-A*02:101", mhc_class="I", definition_id=PSEUDO_I_DEFINITION
        )
        legacy = self.registry.allele_to_pseudo(
            "HLA-A9201", mhc_class="I", definition_id=PSEUDO_I_DEFINITION
        )
        self.assertEqual(modern.status, "exact_unique")
        self.assertEqual(legacy.status, "exact_unique")
        self.assertEqual(modern.candidates[0], legacy.candidates[0])
        renamed = self.registry.allele_to_pseudo(
            "HLA-B9530", mhc_class="I", definition_id=PSEUDO_I_DEFINITION
        )
        self.assertEqual(renamed.candidates[0].canonical_name, "HLA-B*15:27")
        self.assertEqual(
            self.registry.allele_to_pseudo(
                "HLA-B*95:30", mhc_class="I", definition_id=PSEUDO_I_DEFINITION
            ).status,
            "not_found",
        )
        cw_legacy = self.registry.allele_to_pseudo(
            "HLA-C0734", mhc_class="I", definition_id=PSEUDO_I_DEFINITION
        )
        self.assertEqual(cw_legacy.candidates[0].canonical_name, "HLA-C*07:27")

    def test_three_digit_dpb1_first_field_is_preserved(self):
        correct = self.registry.allele_to_pseudo(
            "HLA-DPA1*01:03/HLA-DPB1*100:01",
            mhc_class="II",
            definition_id=PSEUDO_II_DEFINITION,
        )
        wrong_split = self.registry.allele_to_pseudo(
            "HLA-DPA1*01:03/HLA-DPB1*10:001",
            mhc_class="II",
            definition_id=PSEUDO_II_DEFINITION,
        )
        self.assertEqual(correct.status, "exact_unique")
        self.assertEqual(wrong_split.status, "not_found")

    def test_mhcii_cross_family_pairs_are_rejected(self):
        for pair in (
            "HLA-DQA1*01:02/HLA-DPB1*04:01",
            "HLA-DPA1*01:03/HLA-DQB1*02:01",
        ):
            with self.subTest(pair=pair):
                self.assertEqual(
                    self.registry.molecule_to_chains(pair, mhc_class="II").status,
                    "not_found",
                )

    def test_colon_field_boundaries_do_not_collide_with_compact_names(self):
        partial = self.registry.allele_to_sequence(
            "HLA-A*02:11:02", mhc_class="I", include_partial=True
        )
        other = self.registry.allele_to_sequence("HLA-A*02:1102", mhc_class="I")
        self.assertEqual(partial.status, "exact_unique")
        self.assertEqual(partial.candidates[0].allele, "HLA-A*02:11")
        self.assertEqual(other.status, "exact_unique")
        self.assertEqual(other.candidates[0].allele, "HLA-A*02:1102")

    def test_source_conflict_is_quarantined_not_returned_as_ambiguity(self):
        result = self.registry.allele_to_pseudo(
            "HLA-DPA1*01:03/HLA-DPB1*06:01",
            mhc_class="II",
            definition_id=PSEUDO_II_DEFINITION,
        )
        self.assertEqual(result.status, "exact_unique")
        self.assertEqual(
            result.candidates[0].pseudo_sequence,
            "YAFFMFSGGAILNTLYLQFEYFDLEEVRMHLDVT",
        )

    def test_full_sequence_reverse_preserves_deposited_names(self):
        forward = self.registry.allele_to_sequence("HLA-A*02:01:01:01", mhc_class="I")
        reverse = self.registry.sequence_to_alleles(
            forward.candidates[0].sequence, mhc_class="I", species="human"
        )
        self.assertEqual(reverse.status, "exact_ambiguous")
        aliases = {alias for item in reverse.candidates for alias in item.aliases}
        self.assertIn("HLA-A*02:01:01:01", aliases)

    def test_class_i_molecule_has_species_matched_b2m(self):
        human = self.registry.molecule_to_chains(
            "HLA-A*02:01:01:01", mhc_class="I"
        )
        self.assertEqual(human.status, "exact_unique")
        self.assertEqual(
            [item.allele for item in human.candidates[0].chains],
            ["HLA-A*02:01", "B2M-human"],
        )
        self.assertEqual(len(human.candidates[0].chains[1].sequence), 119)

        mouse = self.registry.molecule_to_chains("H-2-Kb", mhc_class="I")
        self.assertEqual(
            [item.allele for item in mouse.candidates[0].chains],
            ["H-2-Kb", "B2M-mouse"],
        )

    def test_class_ii_human_and_mouse_pairs(self):
        human = self.registry.molecule_to_chains(
            "HLA-DRA*01:01/HLA-DRB1*15:01", mhc_class="II"
        )
        self.assertEqual(human.status, "exact_unique")
        self.assertEqual(
            [item.chain_role for item in human.candidates[0].chains], ["alpha", "beta"]
        )

        mouse = self.registry.molecule_to_chains("H-2-IEk", mhc_class="II")
        self.assertEqual(mouse.status, "exact_unique")
        self.assertEqual(
            [item.source_accession for item in mouse.candidates[0].chains],
            ["P04224", "Q31163"],
        )

    def test_partial_mouse_pair_requires_opt_in(self):
        default = self.registry.molecule_to_chains("H-2-IAq", mhc_class="II")
        self.assertEqual(default.status, "not_found")
        opted_in = self.registry.molecule_to_chains(
            "H-2-IAq", mhc_class="II", include_partial=True
        )
        self.assertEqual(opted_in.status, "exact_unique")
        self.assertFalse(opted_in.candidates[0].chains[0].is_complete)

    def test_exact_partial_hla_name_does_not_fall_back_to_another_allele(self):
        default = self.registry.allele_to_sequence(
            "HLA-A*02:01:17", mhc_class="I"
        )
        self.assertEqual(default.status, "not_found")
        self.assertTrue(any("include_partial" in note for note in default.notes))
        opted_in = self.registry.allele_to_sequence(
            "HLA-A*02:01:17", mhc_class="I", include_partial=True
        )
        self.assertEqual(opted_in.status, "exact_unique")
        self.assertEqual(len(opted_in.candidates[0].sequence), 321)

    def test_dr_beta_shorthand_discloses_alpha_assumption(self):
        result = self.registry.molecule_to_chains("HLA-DRB1*15:01", mhc_class="II")
        self.assertEqual(result.status, "exact_unique")
        self.assertTrue(any("DRA*01:01" in note for note in result.notes))

    def test_molecule_propagates_two_field_projection_notes(self):
        default = self.registry.molecule_to_chains(
            "HLA-A*02:01:999:999", mhc_class="I"
        )
        self.assertEqual(default.status, "not_found")
        result = self.registry.molecule_to_chains(
            "HLA-A*02:01:999:999",
            mhc_class="I",
            allow_two_field_projection=True,
        )
        self.assertEqual(result.status, "projected_unique")
        self.assertTrue(any("two-field" in note for note in result.notes))

    def test_nonexistent_high_field_pseudo_is_not_projected(self):
        result = self.registry.allele_to_pseudo(
            "HLA-A*02:01:999:999",
            mhc_class="I",
            definition_id=PSEUDO_I_DEFINITION,
        )
        self.assertEqual(result.status, "not_found")
        self.assertEqual(
            self.registry.molecule_to_chains(
                "HLA-DRB1*15:01:999:999", mhc_class="II"
            ).status,
            "not_found",
        )
        self.assertEqual(
            self.registry.molecule_to_af3_chains(
                "HLA-DRB1*15:01:999:999", mhc_class="II"
            ).status,
            "not_found",
        )

    def test_known_partial_name_is_never_projected_to_an_af3_chain(self):
        partial = self.registry.allele_to_sequence(
            "HLA-A*02:01:17", mhc_class="I", include_partial=True
        )
        self.assertEqual(partial.status, "exact_unique")
        self.assertFalse(partial.candidates[0].is_complete)
        af3 = self.registry.allele_to_af3(
            "HLA-A*02:01:17",
            mhc_class="I",
            allow_two_field_projection=True,
        )
        self.assertEqual(af3.status, "not_found")
        self.assertTrue(any("partial" in note for note in af3.notes))
        self.assertEqual(
            self.registry.molecule_to_af3_chains(
                "HLA-A*02:01:17",
                mhc_class="I",
                allow_two_field_projection=True,
            ).status,
            "not_found",
        )

        dr_partial = "HLA-DRB1*15:01:22"
        self.assertEqual(
            self.registry.allele_to_sequence(
                dr_partial, mhc_class="II", include_partial=True
            ).status,
            "exact_unique",
        )
        for name in (
            dr_partial,
            f"HLA-DRA*01:01/{dr_partial}",
        ):
            with self.subTest(name=name):
                self.assertEqual(
                    self.registry.molecule_to_af3_chains(
                        name,
                        mhc_class="II",
                        allow_two_field_projection=True,
                    ).status,
                    "not_found",
                )

    def test_dr_high_field_molecule_name_is_preserved(self):
        result = self.registry.molecule_to_af3_chains(
            "HLA-DRB1*15:01:01:01", mhc_class="II"
        )
        self.assertEqual(result.status, "exact_unique")
        self.assertEqual(
            result.candidates[0].canonical_name,
            "HLA-DRA*01:01/HLA-DRB1*15:01:01:01",
        )

    def test_expression_suffix_is_preserved_and_not_leaked_to_bare_name(self):
        self.assertEqual(
            self.registry.allele_to_sequence("HLA-A*03:03", mhc_class="I").status,
            "not_found",
        )
        null = self.registry.allele_to_sequence("HLA-A*03:03N", mhc_class="I")
        self.assertEqual(null.status, "exact_unique")
        self.assertEqual(null.candidates[0].allele, "HLA-A*03:03N")
        af3 = self.registry.allele_to_af3("HLA-A*03:03N", mhc_class="I")
        self.assertEqual(af3.status, "not_found")
        self.assertTrue(any("suffix N" in note for note in af3.notes))
        self.assertEqual(
            self.registry.allele_to_sequence("HLA-A*01:01L", mhc_class="I").status,
            "not_found",
        )
        low = self.registry.allele_to_sequence(
            "HLA-A*01:01:38L", mhc_class="I"
        )
        self.assertEqual(low.status, "exact_unique")
        self.assertEqual(low.candidates[0].allele, "HLA-A*01:01:38L")

    def test_af3_chain_golden_boundaries(self):
        expected = {
            "HLA-A*02:01": (284, 25, 308, "GSHSMRYFFTSVSRPGRGEP", "GLPKPLTLRWEPSSQPTIPI"),
            "B2M-human": (99, 21, 119, "IQRTPKIQVYSRHPAENGKS", "CRVNHVTLSQPKIVKWDRDM"),
            "HLA-DRA*01:01": (191, 26, 216, "IKEEHVIIQAEFYLNPDQSG", "EPLLKHWEFDAPSPLPETTE"),
            "HLA-DRB1*15:01": (198, 30, 227, "GDTRPRFLWQPKRECHFFNG", "SVTSPLTVEWRARSESAQSK"),
            "HLA-DQA1*01:02": (194, 24, 217, "EDIVADHVASCGVNLYQFYG", "QPLLKHWEPEIPAPMSELTE"),
            "HLA-DQB1*06:02": (198, 33, 230, "RDSPEDFVFQFKGMCYFTNG", "SLQSPITVEWRAQSESAQSK"),
            "HLA-DPA1*01:03": (194, 29, 222, "AGAIKADHVSTYAAFVQTHR", "QPLLKHWEAQEPIQMPETTE"),
            "HLA-DPB1*04:01": (196, 30, 225, "RATPENYLFQGRQECYAFNG", "SLDSPVTVEWKAQSDSARSK"),
            "H-2-Kb": (284, 22, 305, "GPHSLRYFVTAVSRPGLGEP", "GLPEPLTLRWEPPPSTVSNM"),
            "H2-Aa*b": (195, 24, 218, "EDDIEADHVGTYGISVYQSP", "EPVLKHWEPEIPAPMSELTE"),
            "H2-Ab1*b": (199, 28, 226, "GDSERHFVYQFMGECYFTNG", "SLKSPITVEWRAQSESAWSK"),
            "B2M-mouse": (99, 21, 119, "IQKTPQIQVYSRHPPENGKP", "CRVKHASMAEPKTVYWDRDM"),
        }
        for name, (length, start, end, prefix, suffix) in expected.items():
            with self.subTest(name=name):
                result = self.registry.allele_to_af3(name)
                self.assertEqual(result.status, "exact_unique")
                candidate = result.candidates[0]
                self.assertEqual(candidate.definition_id, AF3_DEFINITION)
                self.assertEqual(len(candidate.sequence), length)
                self.assertEqual(candidate.precursor_start_1based, start)
                self.assertEqual(candidate.precursor_end_1based, end)
                self.assertTrue(candidate.sequence.startswith(prefix))
                self.assertTrue(candidate.sequence.endswith(suffix))

    def test_every_af3_chain_is_an_exact_parent_subsequence(self):
        parents = {candidate.record_id: candidate for candidate in self.registry._sequences}
        for candidate in self.registry._af3_sequences:
            parent = parents[candidate.source_full_record_id]
            self.assertTrue(parent.is_complete)
            self.assertEqual(candidate.species, parent.species)
            self.assertEqual(candidate.mhc_class, parent.mhc_class)
            self.assertEqual(candidate.chain_role, parent.chain_role)
            self.assertNotIn("X", candidate.sequence)
            if candidate.species == "human":
                self.assertNotIn(candidate.allele[-1], "NLSQCA")
            observed = parent.sequence[
                candidate.precursor_start_1based - 1 : candidate.precursor_end_1based
            ]
            self.assertEqual(candidate.sequence, observed)

    def test_af3_reverse_is_non_lossy(self):
        forward = self.registry.allele_to_af3("HLA-A*02:01", mhc_class="I")
        reverse = self.registry.af3_to_alleles(
            forward.candidates[0].sequence, mhc_class="I", species="human"
        )
        self.assertEqual(reverse.status, "exact_ambiguous")
        names = {
            name
            for candidate in reverse.candidates
            for name in (candidate.allele, *candidate.aliases)
        }
        self.assertIn("HLA-A*02:01", names)

    def test_af3_molecule_returns_two_separate_chains(self):
        class_i = self.registry.molecule_to_af3_chains(
            "HLA-A*02:01", mhc_class="I"
        )
        self.assertEqual(class_i.status, "exact_unique")
        self.assertEqual(
            [(chain.chain_role, len(chain.sequence)) for chain in class_i.candidates[0].chains],
            [("heavy", 284), ("companion", 99)],
        )
        class_ii = self.registry.molecule_to_af3_chains(
            "HLA-DRA*01:01/HLA-DRB1*15:01", mhc_class="II"
        )
        self.assertEqual(class_ii.status, "exact_unique")
        self.assertEqual(
            [(chain.chain_role, len(chain.sequence)) for chain in class_ii.candidates[0].chains],
            [("alpha", 191), ("beta", 198)],
        )
        for pair in (
            "HLA-DQA1*01:02/HLA-DQB1*06:02",
            "HLA-DPA1*01:03/HLA-DPB1*04:01",
        ):
            with self.subTest(pair=pair):
                result = self.registry.molecule_to_af3_chains(pair, mhc_class="II")
                self.assertEqual(result.status, "exact_unique")
                self.assertEqual(
                    [chain.chain_role for chain in result.candidates[0].chains],
                    ["alpha", "beta"],
                )

    def test_registry_contains_both_species_and_classes(self):
        counts = self.registry.metadata["counts"]["pseudo_by_species_class"]
        for key in ("human/I", "human/II", "mouse/I", "mouse/II"):
            self.assertGreater(counts[key], 0)


if __name__ == "__main__":
    unittest.main()

