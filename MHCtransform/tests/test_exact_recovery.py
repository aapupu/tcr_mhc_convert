import unittest
from importlib import resources

from mhctransform import MHCRegistry
from mhctransform.registry import PSEUDO_I_DEFINITION, PSEUDO_II_DEFINITION


class TargetedExactRecoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry = MHCRegistry()

    def test_mouse_class_i_strict_pbr_and_pseudo(self):
        expected = {
            "H-2-Kq": ("e3e3243da63d94dd3ad2", "YHSYYRNIADNSSVDTLYIRYEVYTWAARAYAWH", True),
            "H-2-Ds": ("6b3312606e5956f4fbaf", "YVEVYRNIAGNSFVNTAYVVSQHYTWGDYAYTWY", False),
            "H-2-Ks": ("615ee18ae1483fefb30b", "YVEYYRNNARSSFVNTAYIRYEFYTWGEFAYVWY", False),
            "H-2-Qa-1b": ("afec8a1c938bf54e26fb", "YTIMYREKADMNFVNTLYLWYCEYSSVEQAYPWY", False),
        }
        for name, (pbr_id, pseudo_sequence, official) in expected.items():
            with self.subTest(name=name):
                pbr = self.registry.allele_to_pbr(
                    name, mhc_class="I", species="mouse", require_pseudo_match=True
                )
                pseudo = self.registry.allele_to_pseudo(
                    name,
                    mhc_class="I",
                    species="mouse",
                    definition_id=PSEUDO_I_DEFINITION,
                )
                self.assertEqual(pbr.status, "exact_unique")
                self.assertEqual(pbr.candidates[0].record_id, pbr_id)
                self.assertNotIn("X", pbr.candidates[0].sequence)
                self.assertEqual(pseudo.status, "exact_unique")
                self.assertEqual(pseudo.candidates[0].pseudo_sequence, pseudo_sequence)
                self.assertEqual(pseudo.candidates[0].official_support, official)

    def test_mouse_class_ii_exact_components(self):
        expected = {
            "H-2-IAf": "YSYFWTSGGQTGHILSFGSDYYDYRTETVHGVPT",
            "H-2-IAg7": "YTYHLILGGQAEHILHFGLTYYSYRTETAHETPT",
            "H-2-IAr": "YTFHWTSGGQVVHILYFGLTYYDYRTETVHGPHT",
        }
        for name, pseudo_sequence in expected.items():
            with self.subTest(name=name):
                pbr = self.registry.molecule_to_pbr_chains(
                    name, mhc_class="II", species="mouse", require_pseudo_match=True
                )
                pseudo = self.registry.allele_to_pseudo(
                    name,
                    mhc_class="II",
                    species="mouse",
                    definition_id=PSEUDO_II_DEFINITION,
                )
                self.assertEqual(pbr.status, "exact_unique")
                self.assertEqual(len(pbr.candidates[0].chains), 2)
                self.assertTrue(all("X" not in chain.sequence for chain in pbr.candidates[0].chains))
                self.assertEqual(pseudo.candidates[0].pseudo_sequence, pseudo_sequence)
                self.assertFalse(pseudo.candidates[0].official_support)

    def test_iaf_af3_is_exact_ectodomain_fallback(self):
        molecule = self.registry.molecule_to_af3_chains(
            "H-2-IAf", mhc_class="II", species="mouse", fallback="partial"
        )
        self.assertEqual(molecule.status, "exact_unique")
        self.assertEqual(
            {chain.evidence_level for chain in molecule.candidates[0].chains},
            {"source_exact_ectodomain"},
        )
        self.assertTrue(all(chain.imputed_residue_count == 0 for chain in molecule.candidates[0].chains))

    def test_explicit_dra_0102_pairs_are_derived_not_official(self):
        betas = (
            "HLA-DRB1*03:01",
            "HLA-DRB1*13:01",
            "HLA-DRB3*01:01",
            "HLA-DRB3*02:02",
        )
        for beta in betas:
            name = f"HLA-DRA*01:02/{beta}"
            with self.subTest(name=name):
                pseudo = self.registry.allele_to_pseudo(
                    name,
                    mhc_class="II",
                    species="human",
                    definition_id=PSEUDO_II_DEFINITION,
                )
                pbr = self.registry.molecule_to_pbr_chains(
                    name, mhc_class="II", species="human", require_pseudo_match=True
                )
                self.assertEqual(pseudo.status, "exact_unique")
                self.assertFalse(pseudo.candidates[0].official_support)
                self.assertIn("pair_not_dtu_official", pseudo.candidates[0].quality_flags)
                self.assertEqual(
                    [chain.allele for chain in pbr.candidates[0].chains],
                    ["HLA-DRA*01:02", beta],
                )

    def test_deleted_b4401_resolves_only_to_documented_current_identity(self):
        full = self.registry.allele_to_sequence(
            "HLA-B*44:01", mhc_class="I", species="human"
        )
        pseudo = self.registry.allele_to_pseudo(
            "HLA-B*44:01",
            mhc_class="I",
            species="human",
            definition_id=PSEUDO_I_DEFINITION,
        )
        self.assertEqual([candidate.record_id for candidate in full.candidates], ["894f1d97499b82bf3658"])
        self.assertEqual(pseudo.candidates[0].canonical_name, "HLA-B*44:02")
        self.assertTrue(pseudo.candidates[0].official_support)
        self.assertTrue(any("retired" in note for note in pseudo.notes))
        self.assertTrue(any("current two-field allele" in note for note in pseudo.notes))

    def test_kq_af3_gap_and_eu_alpha_scope_are_explicit(self):
        self.assertEqual(
            self.registry.allele_to_af3(
                "H-2-Kq", mhc_class="I", species="mouse"
            ).status,
            "not_found",
        )
        partial = self.registry.allele_to_sequence(
            "H-2-Kq", mhc_class="I", species="mouse", include_partial=True
        )
        mf352193 = next(
            candidate for candidate in partial.candidates if candidate.record_id == "bc8a685111a592eeab3c"
        )
        self.assertIn("residues 275-284 absent", mf352193.source_version)

        for method in (
            self.registry.allele_to_sequence,
            self.registry.allele_to_af3,
            self.registry.allele_to_pbr,
        ):
            result = method("H2-Ea*u", mhc_class="II", species="mouse")
            self.assertEqual(result.status, "exact_unique")
            self.assertIn("P14439", result.candidates[0].source_accession)

    def test_ieu_has_exact_contacts_but_no_invented_beta_af3_chain(self):
        pseudo = self.registry.allele_to_pseudo(
            "H-2-IEu",
            mhc_class="II",
            species="mouse",
            definition_id=PSEUDO_II_DEFINITION,
        )
        self.assertEqual(pseudo.status, "exact_unique")
        self.assertEqual(
            pseudo.candidates[0].pseudo_sequence,
            "QEYFIASGAAVDAVMGSSFDYYEIQTATYHISFL",
        )
        self.assertEqual(
            pseudo.candidates[0].evidence_level,
            "derived_contact_complete_partial",
        )
        self.assertFalse(pseudo.candidates[0].official_support)

        pbr = self.registry.molecule_to_pbr_chains(
            "H-2-IEu",
            mhc_class="II",
            species="mouse",
            fallback="partial",
            require_pseudo_match=True,
        )
        self.assertEqual(pbr.status, "exact_unique")
        self.assertEqual(pbr.fallback_used, "partial")
        alpha, beta = pbr.candidates[0].chains
        self.assertEqual(alpha.allele, "H2-Ea*u")
        self.assertNotIn("X", alpha.sequence)
        self.assertEqual(beta.allele, "H2-Eb1*u")
        self.assertEqual(beta.sequence.count("X"), 1)
        self.assertEqual(beta.unknown_spans_1based, "1")
        self.assertEqual(beta.evidence_level, "source_partial_masked")

        beta_source = self.registry.allele_to_sequence(
            "H2-Eb1*u", mhc_class="II", species="mouse", include_partial=True
        )
        self.assertEqual(beta_source.status, "exact_unique")
        self.assertEqual(beta_source.candidates[0].source_accession, "M12826.1;AAA39645.1;Q31185")

        af3 = self.registry.molecule_to_af3_chains(
            "H-2-IEu", mhc_class="II", species="mouse", fallback="partial"
        )
        self.assertEqual(af3.status, "not_found")

    def test_category_four_names_are_absent_from_targeted_sidecars(self):
        data = resources.files("mhctransform").joinpath("data")
        text = "\n".join(
            data.joinpath(name).read_text(encoding="utf-8")
            for name in (
                "exact_recovery_molecule_pairs.csv",
                "exact_recovery_pseudo_sequences.csv",
            )
        )
        for excluded in (
            "HLA-DQA1*01:02/DRB1*15:01",
            "H-2-IEb",
            "H-2-IEg7",
        ):
            self.assertNotIn(excluded, text)

    def test_runtime_stats_include_targeted_recovery_overlays(self):
        expected = self.registry.metadata["targeted_recovery"]["runtime_counts"]
        stats = self.registry.stats()
        self.assertEqual(stats["sequence_records"], expected["full_sequences"])
        self.assertEqual(
            stats["complete_sequence_records"],
            expected["complete_full_sequences"],
        )
        self.assertEqual(stats["af3_sequence_records"], expected["af3_sequences"])
        self.assertEqual(stats["pbr_sequence_records"], expected["strict_pbr_sequences"])
        self.assertEqual(stats["fallback_af3_records"], expected["fallback_af3_sequences"])
        self.assertEqual(stats["fallback_pbr_records"], expected["fallback_pbr_sequences"])
        self.assertEqual(stats["pseudo_records"], expected["pseudo_sequences"])
        self.assertEqual(
            stats["pseudo_pbr_projection_records"],
            expected["pseudo_pbr_projection_components"],
        )
        self.assertEqual(stats["pair_records"], expected["molecule_pairs"])


if __name__ == "__main__":
    unittest.main()
