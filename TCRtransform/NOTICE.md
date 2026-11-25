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

- `imgt_tcr_alleles.csv`: official IMGT human/mouse TRAV, TRAJ, TRBV and TRBJ reference directory snapshot, retrieved 2026-11-25. Accessions and original functionality are retained per row; source URLs and download SHA-256 values are in `imgt_tcr_alleles.metadata.json`. Source: [IMGT reference sequences](https://www.imgt.org/vquest/refseqh.html).
- `tcr_gene_aa_human_TRA_TRB_gene.csv` and `tcr_gene_aa_mouse_TRA_TRB_gene.csv`: inherited gene tables corrected against that snapshot, with explicit representative alleles and disabled nonfunctional records. Adjacent metadata files record current checksums and selections. Model-ready constant sequences are inherited unchanged.
- `tcr_gene_junction_reference.csv`: zero-based V/J CDR3 anchors and TRAJ–TRAC split-codon bridge residues derived from the official snapshot. Runtime also reads allele-specific junctions from the allele table. Adjacent metadata records the source.

## Method references

- Heather JM et al. Stitchr: stitching coding TCR nucleotide sequences from V/J/CDR3 information. *Nucleic Acids Research*. 2022;50(12):e68. <https://doi.org/10.1093/nar/gkac190>
- IMGT FR-IMGT and CDR-IMGT definitions: <https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html>

The package code is licensed under the MIT License; see `LICENSE`. The MIT License does not override or relicense upstream database content, and users should confirm the applicable terms before redistributing derived reference data.

