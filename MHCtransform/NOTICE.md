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
- 34-aa artifacts: the public PUFFIN data mirror pinned at commit
  `eddf17d972285a89e80e2837adc1ee6308f5b96b`. See
  <https://github.com/gifford-lab/PUFFIN> and the original NetMHCpan /
  NetMHCIIpan publications cited in `AUDIT_REPORT.md`.
- `fallback_*_sequences.csv` contains project-derived partial annotations and
  synthetic same-locus scaffolds. Each row identifies its deposited source,
  complete template, pseudo record, and observed/imputed spans. These rows are
  not official IPD/UniProt allele protein sequences and remain subject to the
  upstream data terms above.

Every generated source is identified by a SHA-256 digest in `metadata.json`.
Users redistributing data should review the current terms of each upstream
source.

