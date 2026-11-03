# tcr_mhc_convert

**English** | [简体中文](README.zh-CN.md)

`tcr_mhc_convert` is an auditable, offline-capable Python monorepo for converting human and mouse MHC and alpha/beta TCR annotations into protein sequences, with exact reverse lookup where the mapping permits it.

This repository was renamed from `aapupu/tcr_design`.

## Packages

| Package | Main capabilities | CLI |
|---|---|---|
| [`MHCtransform`](MHCtransform/) | Exact bidirectional conversion between human/mouse MHC-I or MHC-II allele names and precursor proteins, AF3-ready extracellular chains, or 34-aa pseudosequences; species-matched B2M for MHC-I and separate alpha/beta chains for MHC-II | `mhctransform` |
| [`TCRtransform`](TCRtransform/) | Human/mouse alpha/beta TCR extracellular-chain reconstruction from V + CDR3 + J + C, CDR1/2/2.5/3 annotation, CDR3 boundary analysis, and reverse parsing of complete extracellular chains | `tcrtransform` |

Both packages use only the Python standard library at runtime and support Python 3.9 or newer. They share one GitHub repository while retaining independent package names, versions, APIs, and command-line programs, so they can be installed separately or together.

## Installation

```bash
git clone https://github.com/aapupu/tcr_mhc_convert.git
cd tcr_mhc_convert

python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ./MHCtransform -e ./TCRtransform
```

## Quick start

### MHC allele to precursor, AF3 chain, or 34-aa pseudo

```bash
# HLA-A heavy-chain precursor
mhctransform allele-to-full 'HLA-A*02:01' --class I

# AF3-ready HLA-A extracellular heavy chain + human B2M
mhctransform molecule-to-af3 'HLA-A*02:01' --class I

# AF3-ready MHC-II alpha/beta chains
mhctransform molecule-to-af3 \
  'HLA-DRA*01:01/HLA-DRB1*15:01' --class II

# NetMHCpan-compatible 34-aa descriptor
mhctransform allele-to-pseudo 'HLA-A*02:01' \
  --class I --definition netmhcpan-34-v1
```

Use `full-to-allele`, `af3-to-allele`, and `pseudo-to-allele` for reverse lookup. Reverse queries retain every exact candidate and never present a many-to-one mapping as a unique allele. See the [`MHCtransform` documentation](MHCtransform/README.md) and [audit report](MHCtransform/AUDIT_REPORT.md) for definitions, coverage, limitations, and provenance.

### TCR V/J/CDR3 to AF3-ready alpha/beta extracellular chains

Prepare a CSV containing at least `TRAV, TRAJ, CDR3A, TRBV, TRBJ, CDR3B`:

```bash
tcrtransform build TCRtransform/examples/tcr_input.csv \
  --species-column species \
  -o tcr_built.csv
```

Use `tcrtransform annotate` to annotate CDR regions without reconstruction, or `tcrtransform parse` to infer V/J/C genes and CDRs from existing extracellular chains. See the [`TCRtransform` documentation](TCRtransform/README.md) and its [Stitchr/tcrdist3 comparative validation](TCRtransform/validation/REPORT.md).

## AlphaFold 3 pMHC-TCR five-chain input

Provide each molecule as an independent protein entity. Do not join chains with artificial linkers, and never use a non-contiguous 34-aa pseudosequence as a structural chain.

| AF3 entity | MHC-I complex | MHC-II complex |
|---|---|---|
| 1 | Peptide | Peptide |
| 2 | MHC heavy-chain extracellular region | MHC alpha-chain extracellular region |
| 3 | B2M | MHC beta-chain extracellular region |
| 4 | Mature TCR alpha V+C extracellular chain | Mature TCR alpha V+C extracellular chain |
| 5 | Mature TCR beta V+C extracellular chain | Mature TCR beta V+C extracellular chain |

The MHC profile `af3-soluble-extracellular-v1` removes the signal peptide, transmembrane helix, and cytoplasmic tail while retaining all extracellular folded domains and the membrane-proximal connector. `TCRtransform` produces mature TCR V+C extracellular structural chains. AlphaFold 3 specifies that these amino-acid chains should be supplied as separate entities but does not prescribe universal crop coordinates; these retained regions are the repository's audited engineering profile.

## Tests

```bash
PYTHONPATH=MHCtransform/src:TCRtransform/src \
  python -m unittest discover -s MHCtransform/tests -p 'test*.py' -v

PYTHONPATH=MHCtransform/src:TCRtransform/src \
  python -m unittest discover -s TCRtransform/tests -p 'test*.py' -v
```

The root [GitHub Actions workflow](.github/workflows/tests.yml) installs and tests both packages on Python 3.9 and 3.11.

## Repository layout

```text
tcr_mhc_convert/
|-- MHCtransform/          # mhctransform package, data, audit, and tests
|-- TCRtransform/          # tcrtransform package, validation, and tests
|-- legacy/tcr_design_v1/ # predecessor files retained for provenance
|-- .github/workflows/     # combined continuous integration
|-- LICENSE                # MIT code license
|-- README.zh-CN.md        # optional Chinese documentation
`-- README.md              # default English documentation
```

## License and data terms

Repository code is licensed under the [MIT License](LICENSE). Each package's `NOTICE.md` and data metadata record the pinned versions, checksums, and provenance of upstream resources including IPD-IMGT/HLA, UniProt, tcrdist3, and PUFFIN. The MIT License applies to repository code; it does not relicense bundled or derived biological reference data. Review the applicable upstream terms before redistributing those data.
