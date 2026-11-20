import unittest

from mhctransform import MHCRegistry, MHCTransformError
from mhctransform.fallback import expand_positions
from mhctransform.pbr import (
    PBR_I_DEFINITION,
    PBR_II_ALPHA_DEFINITION,
    PBR_II_BETA_DEFINITION,
)
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
            "HLA-C*03:46": "YYAGYREKYRQTDVSNLYIRYDYYTWAELAYLWY",
            "H-2-Kb": "YVEYYREKAGNSFVDTLYIVSQYYTWAELAYTWY",
        }
        for name, pseudo in expected.items():
            with self.subTest(name=name):
                result = self.registry.allele_to_pseudo(
                    name, mhc_class="I", definition_id=PSEUDO_I_DEFINITION
                )
                self.assertEqual(result.status, "exact_unique")
                self.assertEqual(result.candidates[0].pseudo_sequence, pseudo)

    def test_netmhc_historical_q_names_keep_current_ipd_canonical(self):
        full = self.registry.allele_to_sequence(
            "HLA-DQA1*01:07", mhc_class="II"
        )
        self.assertEqual(full.status, "exact_unique")
        self.assertEqual(full.candidates[0].allele, "HLA-DQA1*01:07Q")

        dq = self.registry.allele_to_pseudo(
            "HLA-DQA1*01:07/HLA-DQB1*02:01",
            mhc_class="II",
            definition_id=PSEUDO_II_DEFINITION,
        )
        dr = self.registry.allele_to_pseudo(
            "HLA-DRB1*15:13",
            mhc_class="II",
            definition_id=PSEUDO_II_DEFINITION,
        )
        self.assertEqual(dq.status, "exact_unique")
        self.assertEqual(dr.status, "exact_unique")
        self.assertEqual(dq.candidates[0].alpha_allele, "HLA-DQA1*01:07Q")
        self.assertEqual(dr.candidates[0].beta_allele, "HLA-DRB1*15:13Q")
        self.assertIn(
            "official_name_maps_to_current_q_suffix",
            dq.candidates[0].quality_flags,
        )

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

    def test_pseudo_provenance_separates_validation_from_derivation(self):
        validated = self.registry.allele_to_pseudo(
            "HLA-A*02:402",
            mhc_class="I",
            definition_id=PSEUDO_I_DEFINITION,
        )
        self.assertEqual(validated.status, "exact_unique")
        validated_candidate = validated.candidates[0]
        self.assertEqual(
            validated_candidate.evidence_level,
            "derived_contact_complete_partial",
        )
        self.assertTrue(validated_candidate.official_support)
        self.assertEqual(
            validated_candidate.validation_status, "official_training_exact"
        )
        self.assertIn("IPD-IMGT/HLA", validated_candidate.source_version)
        self.assertTrue(validated_candidate.source_record_ids)
        self.assertEqual(validated.evidence, "derived_contact_complete_partial")

        derived = self.registry.allele_to_pseudo(
            "HLA-A*01:19",
            mhc_class="I",
            definition_id=PSEUDO_I_DEFINITION,
        )
        self.assertEqual(derived.status, "exact_unique")
        derived_candidate = derived.candidates[0]
        self.assertEqual(
            derived_candidate.evidence_level, "derived_contact_complete_partial"
        )
        self.assertTrue(derived_candidate.official_support)
        self.assertEqual(
            derived_candidate.validation_status,
            "derived_contact_complete_partial",
        )
        self.assertTrue(derived_candidate.source_record_ids)
        self.assertNotIn("X", derived_candidate.pseudo_sequence)
        self.assertEqual(derived.evidence, "derived_contact_complete_partial")

        self.assertFalse(
            any(
                candidate.evidence_level == "official_reference"
                for candidate in self.registry._pseudos
            )
        )

    def test_expanded_netmhciipan_supported_pair_is_queryable_and_reversible(self):
        name = "HLA-DQA1*01:01/HLA-DQB1*03:105"
        forward = self.registry.allele_to_pseudo(
            name,
            mhc_class="II",
            definition_id=PSEUDO_II_DEFINITION,
        )
        self.assertEqual(forward.status, "exact_unique")
        candidate = forward.candidates[0]
        self.assertTrue(candidate.official_support)
        self.assertEqual(candidate.evidence_level, "derived_exact_pbr")
        self.assertEqual(
            candidate.validation_status,
            "official_supported_component_composition",
        )
        self.assertEqual(len(candidate.pseudo_sequence), 34)
        self.assertEqual(len(candidate.source_record_ids), 2)
        self.assertIn(
            "computational_pairing_not_haplotype_evidence",
            candidate.quality_flags,
        )

        reverse = self.registry.pseudo_to_alleles(
            candidate.pseudo_sequence,
            mhc_class="II",
            definition_id=PSEUDO_II_DEFINITION,
        )
        self.assertIn(name, {item.canonical_name for item in reverse.candidates})

    def test_contact_complete_partial_pseudos_have_no_imputed_contacts(self):
        candidates = [
            candidate
            for candidate in self.registry._pseudos
            if candidate.evidence_level == "derived_contact_complete_partial"
        ]
        self.assertGreater(len(candidates), 0)
        for candidate in candidates:
            with self.subTest(name=candidate.canonical_name):
                self.assertEqual(len(candidate.pseudo_sequence), 34)
                self.assertNotIn("X", candidate.pseudo_sequence)
                self.assertTrue(candidate.source_record_ids)
                self.assertFalse(
                    any("imput" in flag.lower() for flag in candidate.quality_flags)
                )

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
            if not parent.is_complete:
                self.assertEqual(candidate.source_accession, "P14437")
                self.assertEqual(candidate.allele, "H2-Aa*s")
                self.assertEqual(
                    candidate.boundary_method,
                    "reviewed UniProt exact extracellular topology on complete mature membrane sequence",
                )
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

    def test_partial_fallback_is_not_labelled_full(self):
        strict = self.registry.allele_to_sequence("HLA-B*08:03", mhc_class="I")
        self.assertEqual(strict.status, "not_found")
        fallback = self.registry.allele_to_sequence(
            "HLA-B*08:03", mhc_class="I", fallback="partial"
        )
        self.assertEqual(fallback.status, "exact_unique")
        self.assertEqual(fallback.representation, "source-partial-v1")
        self.assertEqual(fallback.evidence, "source_partial")
        self.assertEqual(len(fallback.candidates[0].sequence), 181)
        self.assertEqual(fallback.candidates[0].imputed_residue_count, 0)

    def test_partial_guided_af3_fallback_preserves_pseudo_contract(self):
        positions = self.registry.metadata["pseudo_definitions"][
            PSEUDO_I_DEFINITION
        ]["mature_reference_positions_1_based"]
        for name in (
            "HLA-B*08:03",
            "HLA-B*45:06",
            "HLA-A*02:50",
            "HLA-A*03:19",
            "H-2-Lq",
        ):
            with self.subTest(name=name):
                pseudo = self.registry.allele_to_pseudo(
                    name, mhc_class="I", definition_id=PSEUDO_I_DEFINITION
                ).candidates[0]
                result = self.registry.allele_to_af3(
                    name, mhc_class="I", fallback="pseudo-scaffold"
                )
                self.assertEqual(result.status, "exact_unique")
                candidate = result.candidates[0]
                self.assertEqual(candidate.evidence_level, "partial_guided_imputed")
                self.assertTrue(candidate.source_record_ids)
                self.assertTrue(candidate.template_allele)
                observed = "".join(candidate.sequence[position - 1] for position in positions)
                self.assertEqual(observed, pseudo.pseudo_sequence)

    def test_complete_ectodomain_fragment_is_not_labelled_imputed(self):
        result = self.registry.allele_to_af3(
            "HLA-A*02:134", mhc_class="I", fallback="partial"
        )
        self.assertEqual(result.status, "exact_unique")
        candidate = result.candidates[0]
        self.assertEqual(candidate.evidence_level, "source_exact_ectodomain")
        self.assertEqual(candidate.definition_id, AF3_DEFINITION)
        self.assertEqual(candidate.observed_residue_count, len(candidate.sequence))
        self.assertEqual(candidate.imputed_residue_count, 0)
        self.assertNotIn(
            "synthetic_construct_not_official_sequence", candidate.quality_flags
        )

    def test_pseudo_only_scaffold_is_explicit_lowest_evidence(self):
        name = "HLA-A*01:04"
        self.assertEqual(
            self.registry.allele_to_sequence(
                name, mhc_class="I", include_partial=True
            ).status,
            "not_found",
        )
        result = self.registry.allele_to_af3(
            name, mhc_class="I", fallback="pseudo-scaffold"
        )
        self.assertEqual(result.status, "exact_unique")
        candidate = result.candidates[0]
        self.assertEqual(candidate.evidence_level, "pseudo_only_scaffold")
        self.assertEqual(candidate.observed_residue_count, 0)
        self.assertEqual(candidate.imputed_residue_count, len(candidate.sequence))

    def test_ias_exact_ectodomain_and_iaq_iau_imputed_alpha(self):
        ias = self.registry.molecule_to_af3_chains("H-2-IAs", mhc_class="II")
        self.assertEqual(ias.status, "exact_unique")
        self.assertEqual(
            [(chain.allele, len(chain.sequence)) for chain in ias.candidates[0].chains],
            [("H2-Aa*s", 195), ("H2-Ab1*s", 197)],
        )
        self.assertEqual(
            self.registry.molecule_to_chains("H-2-IAs", mhc_class="II").status,
            "not_found",
        )
        ias_full = self.registry.molecule_to_chains(
            "H-2-IAs", mhc_class="II", fallback="imputed"
        )
        self.assertEqual(ias_full.status, "exact_unique")
        self.assertEqual(
            ias_full.candidates[0].chains[0].representation,
            "precursor-imputed-v1",
        )
        self.assertEqual(
            ias_full.candidates[0].chains[0].evidence_level,
            "partial_guided_imputed",
        )
        for name in ("H-2-IAq", "H-2-IAu"):
            with self.subTest(name=name):
                strict = self.registry.molecule_to_af3_chains(name, mhc_class="II")
                self.assertEqual(strict.status, "not_found")
                fallback = self.registry.molecule_to_af3_chains(
                    name, mhc_class="II", fallback="pseudo-scaffold"
                )
                self.assertEqual(fallback.status, "exact_unique")
                alpha, beta = fallback.candidates[0].chains
                self.assertEqual(alpha.evidence_level, "partial_guided_imputed")
                self.assertEqual(alpha.unknown_pseudo_positions, (1, 2))
                self.assertEqual(beta.allele, "H2-Ab1*" + name[-1].lower())

    def test_expression_qualified_fallback_requires_override(self):
        name = "HLA-B*07:44N"
        default = self.registry.allele_to_af3(
            name, mhc_class="I", fallback="pseudo-scaffold"
        )
        self.assertEqual(default.status, "not_found")
        explicit = self.registry.allele_to_af3(
            name,
            mhc_class="I",
            fallback="pseudo-scaffold",
            allow_expression_qualified=True,
        )
        self.assertEqual(explicit.status, "exact_unique")
        self.assertIn("expression_qualified_target", explicit.candidates[0].quality_flags)

    def test_fallback_reverse_is_opt_in_and_never_claims_biological_exactness(self):
        forward = self.registry.allele_to_af3(
            "HLA-B*08:03", mhc_class="I", fallback="imputed"
        )
        sequence = forward.candidates[0].sequence
        self.assertEqual(
            self.registry.af3_to_alleles(sequence, mhc_class="I").status,
            "not_found",
        )
        reverse = self.registry.af3_to_alleles(
            sequence, mhc_class="I", include_fallback="imputed"
        )
        self.assertTrue(reverse.status.startswith("compatible_"))
        self.assertTrue(any("engineered" in note for note in reverse.notes))

    def test_fallback_reverse_merges_strict_sequence_collisions(self):
        forward = self.registry.allele_to_af3(
            "HLA-A*01:04", mhc_class="I", fallback="pseudo-scaffold"
        )
        sequence = forward.candidates[0].sequence
        strict = self.registry.af3_to_alleles(sequence, mhc_class="I")
        strict_names = {
            name
            for candidate in strict.candidates
            for name in (candidate.allele, *candidate.aliases)
        }
        self.assertNotIn("HLA-A*01:04", strict_names)

        compatible = self.registry.af3_to_alleles(
            sequence, mhc_class="I", include_fallback="pseudo-scaffold"
        )
        compatible_names = {
            name
            for candidate in compatible.candidates
            for name in (candidate.allele, *candidate.aliases)
        }
        self.assertEqual(compatible.status, "compatible_ambiguous")
        self.assertIn("HLA-A*01:04", compatible_names)
        self.assertTrue(strict_names.issubset(compatible_names))

        full_forward = self.registry.allele_to_sequence(
            "HLA-A*01:04", mhc_class="I", fallback="pseudo-scaffold"
        )
        full_sequence = full_forward.candidates[0].sequence
        full_strict = self.registry.sequence_to_alleles(
            full_sequence, mhc_class="I"
        )
        full_compatible = self.registry.sequence_to_alleles(
            full_sequence,
            mhc_class="I",
            include_fallback="pseudo-scaffold",
        )
        full_strict_names = {
            name
            for candidate in full_strict.candidates
            for name in (candidate.allele, *candidate.aliases)
        }
        full_compatible_names = {
            name
            for candidate in full_compatible.candidates
            for name in (candidate.allele, *candidate.aliases)
        }
        self.assertNotIn("HLA-A*01:04", full_strict_names)
        self.assertIn("HLA-A*01:04", full_compatible_names)
        self.assertTrue(full_compatible.status.startswith("compatible_"))

    def test_partial_fallback_reverse_accepts_source_unknown_residue(self):
        forward = self.registry.allele_to_sequence(
            "HLA-B*07:44N", mhc_class="I", fallback="partial"
        )
        self.assertIn("X", forward.candidates[0].sequence)
        reverse = self.registry.sequence_to_alleles(
            forward.candidates[0].sequence,
            mhc_class="I",
            include_fallback="partial",
        )
        self.assertTrue(reverse.status.startswith("compatible_"))
        self.assertEqual(forward.candidates[0].imputed_residue_count, 0)
        self.assertIn(
            "source_contains_unknown_x", forward.candidates[0].quality_flags
        )

    def test_reverse_partial_modes_are_mutually_exclusive(self):
        forward = self.registry.allele_to_sequence(
            "HLA-B*08:03", mhc_class="I", fallback="partial"
        )
        with self.assertRaises(MHCTransformError):
            self.registry.sequence_to_alleles(
                forward.candidates[0].sequence,
                mhc_class="I",
                include_partial=True,
                include_fallback="partial",
            )

    def test_class_ii_coverage_honors_species_filter(self):
        coverage = self.registry.coverage(
            "HLA-DRA*01:01/HLA-DRB1*15:01",
            mhc_class="II",
            species="mouse",
        )
        for key, value in coverage.items():
            if key in {"query", "mhc_class"}:
                continue
            self.assertEqual(value["status"], "not_found", key)

    def test_pbr_golden_boundaries_and_parent_spans(self):
        expected = {
            "HLA-A*02:01": ("I", PBR_I_DEFINITION, 182),
            "HLA-DRA*01:01": ("II", PBR_II_ALPHA_DEFINITION, 84),
            "HLA-DRB1*15:01": ("II", PBR_II_BETA_DEFINITION, 95),
            "HLA-DQA1*01:01": ("II", PBR_II_ALPHA_DEFINITION, 96),
            "HLA-DPB1*01:01": ("II", PBR_II_BETA_DEFINITION, 92),
        }
        af3_by_id = {
            candidate.record_id: candidate for candidate in self.registry._af3_sequences
        }
        for name, (mhc_class, definition, length) in expected.items():
            with self.subTest(name=name):
                result = self.registry.allele_to_pbr(name, mhc_class=mhc_class)
                self.assertEqual(result.status, "exact_unique")
                candidate = result.candidates[0]
                self.assertEqual(candidate.definition_id, definition)
                self.assertEqual(len(candidate.sequence), length)
                parent = af3_by_id[candidate.source_af3_record_id]
                self.assertEqual(
                    candidate.sequence,
                    parent.sequence[
                        candidate.af3_start_1based - 1 : candidate.af3_end_1based
                    ],
                )
                self.assertEqual(
                    candidate.precursor_end_1based
                    - candidate.precursor_start_1based
                    + 1,
                    length,
                )

    def test_pbr_molecule_keeps_biological_chains_separate(self):
        class_i = self.registry.molecule_to_pbr_chains(
            "HLA-A*02:01", mhc_class="I"
        )
        self.assertEqual(class_i.status, "exact_unique")
        self.assertEqual(len(class_i.candidates[0].chains), 1)
        self.assertEqual(class_i.candidates[0].chains[0].chain_role, "heavy")

        class_ii = self.registry.molecule_to_pbr_chains(
            "HLA-DRA*01:01/HLA-DRB1*15:01", mhc_class="II"
        )
        self.assertEqual(class_ii.status, "exact_unique")
        self.assertEqual(
            [chain.chain_role for chain in class_ii.candidates[0].chains],
            ["alpha", "beta"],
        )
        self.assertEqual(
            [len(chain.sequence) for chain in class_ii.candidates[0].chains],
            [84, 95],
        )

    def test_reference_pseudos_are_exact_projections_of_the_pbr_chains(self):
        definitions = self.registry.metadata["pseudo_definitions"]

        class_i_pbr = self.registry.allele_to_pbr(
            "HLA-A*02:01", mhc_class="I"
        ).candidates[0]
        class_i_positions = definitions[PSEUDO_I_DEFINITION][
            "mature_reference_positions_1_based"
        ]
        class_i_projection = "".join(
            class_i_pbr.sequence[position - 1] for position in class_i_positions
        )
        class_i_pseudo = self.registry.allele_to_pseudo(
            "HLA-A*02:01",
            mhc_class="I",
            definition_id=PSEUDO_I_DEFINITION,
        ).candidates[0]
        self.assertEqual(class_i_projection, class_i_pseudo.pseudo_sequence)

        class_ii_pbr = self.registry.molecule_to_pbr_chains(
            "HLA-DRA*01:01/HLA-DRB1*15:01", mhc_class="II"
        ).candidates[0]
        alpha_positions = definitions[PSEUDO_II_DEFINITION][
            "alpha_mature_reference_positions_1_based"
        ]
        beta_positions = definitions[PSEUDO_II_DEFINITION][
            "beta_mature_reference_positions_1_based"
        ]
        class_ii_projection = "".join(
            class_ii_pbr.chains[0].sequence[position - 1]
            for position in alpha_positions
        ) + "".join(
            class_ii_pbr.chains[1].sequence[position - 1]
            for position in beta_positions
        )
        class_ii_pseudo = self.registry.allele_to_pseudo(
            "HLA-DRA*01:01/HLA-DRB1*15:01",
            mhc_class="II",
            definition_id=PSEUDO_II_DEFINITION,
        ).candidates[0]
        self.assertEqual(class_ii_projection, class_ii_pseudo.pseudo_sequence)

    def test_pbr_fallback_inherits_af3_evidence_and_remaps_spans(self):
        af3 = self.registry.allele_to_af3(
            "HLA-B*08:03", mhc_class="I", fallback="imputed"
        ).candidates[0]
        result = self.registry.allele_to_pbr(
            "HLA-B*08:03", mhc_class="I", fallback="imputed"
        )
        self.assertEqual(result.status, "exact_unique")
        candidate = result.candidates[0]
        self.assertEqual(candidate.source_af3_record_id, af3.record_id)
        self.assertEqual(candidate.sequence, af3.sequence[:182])
        self.assertEqual(candidate.evidence_level, af3.evidence_level)
        observed = set(expand_positions(candidate.observed_spans_1based))
        imputed = set(expand_positions(candidate.imputed_spans_1based))
        self.assertFalse(observed & imputed)
        self.assertEqual(observed | imputed, set(range(1, 183)))
        self.assertEqual(candidate.observed_residue_count, len(observed))
        self.assertEqual(candidate.imputed_residue_count, len(imputed))

    def test_pbr_reverse_fallback_is_explicitly_compatible(self):
        strict = self.registry.allele_to_pbr("HLA-A*02:01", mhc_class="I")
        reverse = self.registry.pbr_to_alleles(
            strict.candidates[0].sequence,
            mhc_class="I",
            definition_id=PBR_I_DEFINITION,
        )
        self.assertTrue(reverse.status.startswith("exact_"))

        fallback = self.registry.allele_to_pbr(
            "HLA-B*08:03", mhc_class="I", fallback="imputed"
        )
        hidden = self.registry.pbr_to_alleles(
            fallback.candidates[0].sequence, mhc_class="I"
        )
        compatible = self.registry.pbr_to_alleles(
            fallback.candidates[0].sequence,
            mhc_class="I",
            include_fallback="imputed",
        )
        self.assertNotIn(
            "HLA-B*08:03", {candidate.allele for candidate in hidden.candidates}
        )
        self.assertTrue(compatible.status.startswith("compatible_"))

    def test_every_strict_pbr_is_a_contiguous_af3_subsequence(self):
        af3_by_id = {
            candidate.record_id: candidate for candidate in self.registry._af3_sequences
        }
        for candidate in self.registry._pbr_sequences:
            parent = af3_by_id[candidate.source_af3_record_id]
            self.assertEqual(candidate.af3_start_1based, 1)
            self.assertEqual(
                candidate.sequence,
                parent.sequence[: candidate.af3_end_1based],
                candidate.allele,
            )

    def test_every_pseudo_name_has_an_explicit_af3_fallback(self):
        for pseudo in self.registry._pseudos:
            af3_result = self.registry.molecule_to_af3_chains(
                pseudo.canonical_name,
                mhc_class=pseudo.mhc_class,
                fallback="pseudo-scaffold",
                allow_expression_qualified=True,
            )
            if af3_result.status == "not_found":
                self.fail(f"missing AF3 fallback for {pseudo.canonical_name}")
            full_result = self.registry.molecule_to_chains(
                pseudo.canonical_name,
                mhc_class=pseudo.mhc_class,
                fallback="pseudo-scaffold",
            )
            if full_result.status == "not_found":
                self.fail(f"missing full fallback for {pseudo.canonical_name}")
            pbr_result = self.registry.molecule_to_pbr_chains(
                pseudo.canonical_name,
                mhc_class=pseudo.mhc_class,
                fallback="pseudo-scaffold",
                allow_expression_qualified=True,
            )
            if pbr_result.status == "not_found":
                self.fail(f"missing PBR fallback for {pseudo.canonical_name}")
            for molecule in full_result.candidates:
                for chain in molecule.chains:
                    self.assertNotEqual(
                        getattr(chain, "representation", "precursor_full"),
                        "source-partial-v1",
                        pseudo.canonical_name,
                    )

    def test_every_imputed_full_uses_a_complete_precursor_template(self):
        parents = {
            candidate.record_id: candidate for candidate in self.registry._sequences
        }
        for candidate in self.registry._fallback_sequences:
            if candidate.representation != "precursor-imputed-v1":
                continue
            self.assertIn(candidate.template_record_id, parents)
            self.assertTrue(
                parents[candidate.template_record_id].is_complete,
                candidate.allele,
            )


if __name__ == "__main__":
    unittest.main()

