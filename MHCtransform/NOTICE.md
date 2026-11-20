# Data provenance and terms

The source code in this directory is licensed under the MIT License; see
`LICENSE`. Bundled biological reference records retain their original source
attribution and applicable terms. The MIT License does not relicense those
records.

- Human HLA deposited proteins and sequence-status metadata: pinned
  IPD-IMGT/HLA 3.66.0, IMGTHLA commit `b2055ea`. See
  <https://www.ebi.ac.uk/ipd/imgt/hla/> and
  <https://github.com/ANHIG/IMGTHLA>.
- Mouse MHC and human/mouse B2M proteins: UniProtKB accessions listed in
  `src/mhctransform/data/metadata.json`. See
  <https://www.uniprot.org/help/license>.
- AF3 extracellular-chain boundaries: reviewed UniProtKB topology features
  projected through the pinned IPD-IMGT/HLA 3.66.0 protein alignments. The
  project profile is an engineering convention, not an AlphaFold 3 official
  crop definition.
- Peptide-binding-region boundaries: reviewed UniProtKB alpha1/alpha2,
  alpha1, and beta1 REGION annotations projected through the same pinned
  alignments. Strict PBRs are verified contiguous subsequences of their
  deposited parent records.
- Official support-name snapshots: NetMHCpan 4.2f `allele.list` and
  NetMHCIIpan 4.3k `alleles_name.txt`, committed under `tools/reference/`
  with the source URLs and SHA-256 values recorded in `metadata.json`.
- 34-aa definitions: the original NetMHCpan / NetMHCIIpan publications cited
  in `AUDIT_REPORT.md`. Direct primary pseudosequence inputs are committed as
  `tools/reference/NetMHCpan-4.2e-training-pseudoseqs.dat` (SHA-256
  `bb8fe764cd187e55d7f83a465eef152419151ad9924825c411ba4e8ed84ef5a7`)
  and `tools/reference/NetMHCIIpan-4.3k-training-pseudosequence.2023.dat`
  (SHA-256
  `daf3b7d693e76d8e8102c374d130467ff1ef7f9311064a730d56a753a3d514b5`).
  The `4.2e` label identifies the supplied NetMHCpan pseudosequence artifact;
  the separate public support-name snapshot is NetMHCpan 4.2f. Derived and
  legacy entries retain their own evidence labels. The legacy PUFFIN mirror
  remains pinned at commit
  `eddf17d972285a89e80e2837adc1ee6308f5b96b`.
- The two DTU pseudosequence text snapshots above are bundled data inputs, but
  NetMHCpan and NetMHCIIpan executables, model files, complete training
  archives, and separately licensed software packages are not part of this
  repository. See the official service pages:
  <https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/> and
  <https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/>.
- `fallback_*_sequences.csv` contains project-derived partial annotations and
  synthetic same-locus scaffolds. Each row identifies its deposited source,
  complete template, pseudo record, and observed/imputed spans. These rows are
  not official IPD/UniProt allele protein sequences and remain subject to the
  upstream data terms above.
- `pseudo_pbr_projection.csv` is a project-generated audit table. It requires
  every pseudo contact index to map to a strict or explicitly labeled fallback
  PBR, requires every non-`X` residue to match exactly, and treats DTU `X`
  tokens only as documented alignment/unknown wildcards. It does not turn a
  synthetic fallback into an official allele protein sequence.

Source inputs are identified by SHA-256 digests in `metadata.json` and the
validation report; generated outputs are stably sorted and deterministic, but
the metadata does not claim to list a hash for every output file. Pseudo
validation coverage and quarantined conflicts are recorded in
`pseudo_validation.json`. Users redistributing data should review the current
terms of each upstream source.

