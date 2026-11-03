# MHCI curated sequences with species-matched B2M

Input: mhci_full_ipd_curated.csv
Rows: 22780

Files:
- mhci_full_ipd_curated_legacy_format_species_matched_b2m.csv: previous-style table with Unnamed: 0, MHC, full_aa, species, b2m_aa.
- mhci_full_ipd_curated_model_input_species_matched_b2m.csv: model-ready table with lengths and B2M metadata/status.
- mhci_full_ipd_curated_species_matched_b2m.fasta: paired FASTA, two records per allele: MHC_alpha and B2M.
- b2m_species_mapping.csv: one row per MHC prefix showing the B2M sequence used.

Important status notes:
- species_matched_reviewed / species_matched_structure_supported: use directly.
- species_matched_structure_supported_trimmed_expression_M: RCSB structure chain had an extra initial M; the output uses the mature chain without that M.
- species_unmatched_human_fallback__horse_sequence_not_embedded: Eqca rows use human B2M fallback and are NOT truly species-matched in this release. They are explicitly marked so you can replace them later with Equus caballus B2M if needed.

Counts by prefix/status:
mhc_prefix                                              b2m_match_status     n
      BoLA                           species_matched_structure_supported   147
       DLA      species_matched_structure_supported_trimmed_expression_M     2
      Eqca species_unmatched_human_fallback__horse_sequence_not_embedded    21
      Gogo                 species_matched_reviewed_same_mature_as_human    49
       HLA                                      species_matched_reviewed 20821
      Mamu                           species_matched_structure_supported  1335
      Patr                 species_matched_reviewed_same_mature_as_human   185
       SLA      species_matched_structure_supported_trimmed_expression_M   220
