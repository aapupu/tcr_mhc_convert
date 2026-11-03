# Reference data notice

This repository bundles reference material required for reproducible human and mouse TCR transformation.

## tcrdist3 snapshot

- Upstream project: <https://github.com/kmayerb/tcrdist3>
- Upstream file: `tcrdist/db/combo_xcr_2024-03-05.tsv`
- Git commit: `0e6efa39092e2c64c9c0b52a72c2f5961c161ebb`
- SHA-256: `1d7fcac37ad6d955b3a3f693526e36623186464aef2c528cf8e0322dc89706cd`
- Fields used here: human/mouse TRA/TRB V/J amino-acid references and V-loop annotations.

The underlying germline annotations originate from IMGT. Academic use should cite IMGT and follow its data-use terms. See the upstream [tcrdist3 input notice](https://github.com/kmayerb/tcrdist3/blob/master/docs/inputs.rst).

## Derived references

- `tcr_gene_aa_human_TRA_TRB_gene.csv`: gene-level human TRA/TRB V/J/C amino-acid reference inherited from the predecessor `tcr_design` project. SHA-256: `5f265101426753da8818d4b2a3f729f960476a7178a19703d8b57bf7741d502c`.
- `tcr_gene_aa_mouse_TRA_TRB_gene.csv`: gene-level mouse reference derived from the pinned tcrdist3 snapshot plus model-ready constant-region sequences. SHA-256: `4319446f30c28694c8e777e114640a8ee0905968acb1c79d59d65e0eb4730a05`. Provenance is recorded in its adjacent metadata JSON.
- `tcr_gene_junction_reference.csv`: zero-based V/J CDR3 anchors and TRAJ–TRAC split-codon bridge residues derived from the pinned tcrdist3 snapshot. SHA-256: `f2903e8d4098aa2d857a15bd1f913a3514940e3ee1b57958a583dc0c6f70c02e`. Provenance and unsupported records are recorded in the adjacent metadata JSON.

## Method references

- Heather JM et al. Stitchr: stitching coding TCR nucleotide sequences from V/J/CDR3 information. *Nucleic Acids Research*. 2022;50(12):e68. <https://doi.org/10.1093/nar/gkac190>
- IMGT FR-IMGT and CDR-IMGT definitions: <https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html>

No license for the new package code has been selected in this folder. Add one before public release. A code license does not override upstream database terms.

