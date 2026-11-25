# IMGT TCR corrections — 2026-11-25

## Noncanonical J boundary fix

The `[FWC]XXG` heuristic missed mouse TRAJ44 (L, zero-based index 9), mouse
TRAJ7 (L, index 8), and human TRAJ61 (F, index 8; FGAN rather than FxxG).
Their pinned tcrdist CDR endpoints agree with the pre-update junction table and
map to identical ungapped IMGT sequences (M64239, M64239, M94081 respectively).
`j_anchor_overrides.json` now records those three sequence-validated annotations.
Both online and offline generation use the same resolver. Unknown sequences
are not given a permissive generic L/F anchor; changed annotated sequences fail
generation for review. The targeted refresh changes only three allele anchor
fields and adds their three canonical junction entries, retaining ORF warnings.

Regression tests cover two constructed CDR3 inputs per gene (six cases), full
sequence equality to independently assembled expected chains, reverse parsing,
guarded short J matches, and rejection of changed override sequences. These are
new regression examples, not the user's six private CSV examples.
All six were also executed against archived code/reference commit `e8d444f`;
the reconstructed complete chains match the corrected version character for
character. All 16 tests pass. A second targeted regeneration changes zero anchors.

## Follow-up: rearrangement-aware eligibility

The strict eligibility policy and counts below describe the initial correction
and are superseded by this section. IMGT explicitly notes that a V-terminal or
J-initial stop can disappear during rearrangement; germline F/ORF/P and rearranged
productive/unproductive are different annotations. See the
[official functionality definitions](https://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html).

Runtime now considers references with supported geometry regardless of germline
functionality, preserving literal stop positions until after junction matching.
Removed stops are reported; retained stops prevent model-chain output and require
review of the actual rearrangement. P/ORF warnings do not classify the sample.
Missing geometry remains unresolved rather than biologically unproductive.
The raw CSV `usable` column remains the historical strict-snapshot flag.

Examples tested: mouse `TRBV4*01` + `TRBJ1-5*02` + `CASSQPLF` removes both germline
end stops; human `TRAJ58*01` can lose its leading stop, and `TRAJ1*01` is not
rejected for ORF alone. Retained internal V or J stops still fail reconstruction.
All 14 tests pass, including preservation of failed CSV rows and their original
productive values. New productivity output is explicitly `not_assessed`.
No training/test records were deleted or regenerated in this correction.

This update fixes TCR V/J references and unifies allele resolution across chain
reconstruction, CDR loops and junction metadata. It does not modify MHC data or
rerun the historical Stitchr comparison.

## Sources and coverage

The eight FASTA files are the human and mouse TRAV/TRAJ/TRBV/TRBJ sets in the
[official IMGT reference directory](https://www.imgt.org/vquest/refseqh.html).
Download URLs/hashes are in `src/tcrtransform/data/imgt_tcr_alleles.metadata.json`;
individual accessions and functionality are retained in the CSV. `(F)`/`[F]`
are interpreted according to the [official functionality definitions](https://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html).

| Species | Raw V/J records | Enabled V alleles | Enabled J alleles | Gene representatives |
|---|---:|---:|---:|---:|
| Human | 433 | 263 | 76 | 158 |
| Mouse | 437 | 277 | 68 | 185 |

Enabled records have functional annotations, standard amino acids and the
expected CDR3 anchor. V records must retain the N-terminus and at least 80 aa.
These counts exclude constants. The complete raw snapshot includes P/ORF,
stop-bearing and partial records; their presence does not make them valid
functional-chain reconstruction candidates.

## Corrections

| Reference | Change |
|---|---|
| Human TRAJ24 | Removed the erroneous mixed 22-aa sequence. `*01` (X02887): `TTDSWGKFEFGAGTQVVVTP`; `*02` (M94081): `TTDSWGKLQFGAGTQVVVTP`; `*03` (IMGT000024): `TTDSWGKFQFGAGTQVVVTP`. Gene representative is `*02`. |
| Human TRAJ51, TRAJ55, TRAJ58, TRAJ60; TRBV12-1, TRBV26, TRBV26/OR9-2, TRBV3-2 | Restored deleted stop symbols in raw gene records and disabled their use for functional reconstruction. |
| Mouse TRAV14-1 / TRAV14D-1 | Corrected the swapped `*01` gene associations using IMGT accessions L77149 / L77150. Loops now follow the corrected association. |
| Mouse TRAV13D-1 | Corrected representative `*01` ending `CALE` to `CAME`. |
| Mouse TRAJ49 | Corrected J-region ending `LTVIPS` to `LTVIP`; the junction bridge comes from the split codon, not a copied terminal residue. |
| Mouse TRAJ3 | Corrected raw J-region ending `LRNPPY` to `LRNPP`; retained P status and excluded from functional reconstruction. |
| Mouse TRAJ4 and TRAV12D-3 | Changed gene representatives to functional `*02` records; the old `*01` representatives are nonfunctional. This is a representative-policy change, not a claim that the old allele sequence was mislabeled. |
| Mouse TRBV12-2+TRBV13-2 | Disabled the synthetic combined key as a standalone gene. |
| Shared TRA/DV names | Accept `-DV` spelling as an alias for official `/DV` spelling. |

## Resolution rule

1. Use a supported exact allele when supplied, regardless of whether its suffix
   is `*01`, `*02`, `*03` or another available value.
2. If no suffix is supplied, use the documented gene representative.
3. If the allele is unknown or its functional record is incomplete, use that
   same representative and report `gene_level_fallback`. No neighbouring gene
   substitution is performed. If the gene is absent too, report an error.
4. A known P/ORF or stop-bearing allele is rejected explicitly, so a nonfunctional
   reference cannot be silently converted into a functional chain.

`TRAJ24` and `TRAJ24*99` both resolve to the `*02` representative. C references
remain gene-level. Reverse parsing scores supported alleles but reports unique
gene names, retaining genuinely tied genes rather than counting a gene alias
and its representative allele as two different genes.

## Validation

Eleven focused regression tests cover the three TRAJ24 alleles and exact junctions,
unknown/partial allele fallback, mouse gene corrections, mouse TRA loop columns,
shared TRA/DV aliases, nonfunctional guards, custom references, CSV output and
mouse forward/reverse reconstruction. The existing boundary tests and overlapping
J motif anchors are included. Six human and six mouse enabled records were also
sampled to check consistency between sequence, loop and anchor outputs.
Run with `PYTHONPATH=src python -m unittest discover -s tests -v`.
