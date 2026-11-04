# MHCtransform

**English** | [简体中文](README.zh-CN.md)

`MHCtransform` is an independent Python package in the [`aapupu/tcr_mhc_convert`](https://github.com/aapupu/tcr_mhc_convert) monorepo.

`MHCtransform` is an auditable, offline-capable human/mouse MHC conversion package covering:

- MHC-I heavy chains with species-matched B2M;
- MHC-II alpha/beta chains;
- deposited IPD/UniProt proteins, with records confirmed by official `Full` status or explicit UniProt completeness checks labeled `precursor_full`;
- mature soluble extracellular chains for AF3 complex input under `af3-soluble-extracellular-v1`;
- two strictly separated 34-aa pseudosequence definitions;
- bidirectional lookup from allele/name to sequence and from sequence to **all exact candidates**.

It does not guess a unique allele from a pseudosequence, and it does not call a signal-peptide/TM-trimmed soluble chain full length.

## Representation definitions

| Representation/`definition_id` | Meaning | Primary use |
|---|---|---|
| `precursor_full` | Complete precursor in the source record: signal + extracellular + TM + cytoplasmic regions | Sequence verification and a starting point for expression-construct design; not the default AF3 chain |
| `af3-soluble-extracellular-v1` | Signal peptide, TM helix, and cytoplasmic tail removed; all extracellular folded domains and the membrane-proximal connector retained | AF3 pMHC-TCR multichain structural input |
| `netmhcpan-34-v1` | 34 peptide-contact positions on a mature MHC-I heavy-chain reference alignment; excludes B2M | Pan-allele model descriptor, not a structural chain |
| `netmhciipan-34-v1` | 15 MHC-II alpha-chain positions plus 19 beta-chain positions concatenated from reference alignments; `X` is allowed | Pan-allele model descriptor, not a structural chain |

Position definitions and source hashes are stored in [`metadata.json`](src/mhctransform/data/metadata.json). MHC-II pseudosequences cannot be cut with one universal set of precursor indices because signal-peptide lengths and alignment indels vary. The runtime therefore uses a versioned published pseudosequence artifact instead of slicing full chains on demand.

## Current local registry

| Content | human I | human II | mouse I | mouse II |
|---|---:|---:|---:|---:|
| Deposited protein records | 20,066 | 9,419 | 12 | 17 |
| Complete precursor records available by default | 11,283 | 2,525 | 9 | 14 |
| AF3-ready extracellular-chain records | 11,004 | 2,456 | 9 | 14 |
| 34-aa records | 2,953 | 5,620 | 7 | 8 |

Human HLA records come from pinned IPD-IMGT/HLA 3.66.0 at commit `b2055ea` and are evaluated together with the official `Full/Partial` status in the matching `Allele_status.txt`. Mouse MHC and human/mouse B2M records come from explicit UniProt accessions. AF3 boundaries are projected from reviewed UniProt topology onto the matching IPD protein alignments, and every output is verified as a contiguous subsequence of its parent precursor. The scientific 34-aa definitions come from the NetMHCpan/NetMHCIIpan literature; the distributed local artifact uses the public PUFFIN secondary mirror. Full source accessions, SHA-256 hashes, and generation dates are recorded in the metadata.

Deposited records include incomplete entries. By default, queries return human chains with official `Full` status and no unknown residues, or mouse chains that pass accession, fragment, and length checks. Use `--include-partial` to include fragments explicitly. The current mouse registry is a traceable set of commonly used H-2 alleles and does not claim exhaustive H-2 coverage.

## Usage

There are no third-party runtime dependencies:

```bash
cd MHCtransform
python3 -m pip install -e .
```

### Allele to full precursor

```bash
mhctransform allele-to-full 'HLA-A*02:01:01:01' --class I
mhctransform allele-to-full 'H-2-Kb' --class I
mhctransform allele-to-full 'HLA-DRB1*15:01' --class II
```

A two-field name can still match multiple complete proteins. For example, `HLA-C*03:03` has two official `Full` candidates of 366 and 372 aa in this release and is reported as `exact_ambiguous`; a higher-field name can disambiguate them. The older 364-aa `HLA-A*02:01` entry has official `Partial` status, so the default query returns only the 365-aa `Full` candidate.

A nonexistent high-field name is not silently reduced to two fields. For example, `HLA-A*02:01:999:999` returns `not_found`. Projection occurs only when the caller explicitly supplies `--allow-two-field-projection`, and the result is labeled `projected_unique` or `projected_ambiguous`, never exact. An official exact name for which only a Partial protein is known is not projected even when this option is enabled. HLA expression suffixes `N/L/S/Q/C/A` and their complete official fields are preserved: bare `HLA-A*03:03` does not match `HLA-A*03:03N`, and official `HLA-A*02:01:01:02L` is not rewritten into the nonexistent `HLA-A*02:01L`.

### Complete MHC molecule to two chains

```bash
# MHC-I: heavy chain + species-matched B2M
mhctransform molecule-to-chains 'H-2-Kb' --class I

# MHC-II: explicit alpha/beta pairing
mhctransform molecule-to-chains \
  'HLA-DRA*01:01/HLA-DRB1*15:01' --class II

# Mouse molecule alias
mhctransform molecule-to-chains 'H-2-IEk' --class II
```

When only a human DRB allele is supplied, the package can add `HLA-DRA*01:01` under the common shorthand convention, but this assumption is disclosed in `notes`. DQ and DP require an explicit alpha/beta pair. The package rejects unconstrained Cartesian products and cross-family mismatches such as DQ-alpha/DP-beta or DP-alpha/DQ-beta.

### AF3 extracellular chains and chain pairing

```bash
# One chain: allele -> AF3 chain
mhctransform allele-to-af3 'HLA-A*02:01' --class I

# Reverse: AF3 chain -> all exact candidates
mhctransform af3-to-allele @af3_chain.txt --class I --species human

# Complete MHC-I molecule: 284-aa heavy chain + 99-aa human B2M
mhctransform molecule-to-af3 'HLA-A*02:01' --class I

# Complete MHC-II molecule: separate alpha and beta chains
mhctransform molecule-to-af3 \
  'HLA-DRA*01:01/HLA-DRB1*15:01' --class II
```

`molecule-to-af3` never concatenates MHC chains or adds linkers/tags. MHC-I automatically receives species-matched B2M, while MHC-II requires a same-family alpha/beta pair. The default AF3 registry excludes partial records and human HLA records carrying `N/L/S/Q/C/A` expression suffixes. When the exclusion reason can be identified, it is returned in `notes`.

This crop is an auditable engineering profile, not a Google-prescribed pMHC-TCR coordinate set. The official AF3 input format specifies separate molecular entities with amino-acid sequences but does not prescribe universal crop coordinates. A 34-aa pseudo is a non-contiguous descriptor and must never be supplied as an AF3 protein chain.

### Bidirectional 34-aa conversion

```bash
mhctransform allele-to-pseudo 'HLA-A*02:01' \
  --class I --definition netmhcpan-34-v1

mhctransform pseudo-to-allele \
  YFAMYGEKVAHTHVDTLYVRYHYYTWAVLAYTWY \
  --class I --definition netmhcpan-34-v1

mhctransform allele-to-pseudo 'HLA-DRB1*15:01' \
  --class II --definition netmhciipan-34-v1

mhctransform pseudo-to-allele \
  YSYFLASGGQVVHVLYFGYTYHDIRTETVHGPHT \
  --class II --definition netmhciipan-34-v1
```

Reverse conversion is exact-only. A single 34mer is often shared by many alleles, so the result is commonly `exact_ambiguous`. An unknown 34mer returns `not_found`; the package never substitutes a nearest-neighbor prediction.

### Full-sequence reverse lookup

```bash
mhctransform full-to-allele @sequence.txt --class I --species human
```

One sequence candidate may collect multiple identical deposited HLA names in `aliases`; in that case the result remains `exact_ambiguous`.

### Python API

```python
from mhctransform import MHCRegistry
from mhctransform.registry import AF3_DEFINITION, PSEUDO_I_DEFINITION

registry = MHCRegistry()

forward = registry.allele_to_pseudo(
    "HLA-A*02:01",
    mhc_class="I",
    definition_id=PSEUDO_I_DEFINITION,
)

reverse = registry.pseudo_to_alleles(
    forward.candidates[0].pseudo_sequence,
    mhc_class="I",
    definition_id=PSEUDO_I_DEFINITION,
)
print(reverse.status, reverse.candidate_count)

af3_molecule = registry.molecule_to_af3_chains(
    "HLA-A*02:01",
    mhc_class="I",
    definition_id=AF3_DEFINITION,
)
for chain in af3_molecule.candidates[0].chains:
    print(chain.chain_role, len(chain.sequence), chain.sequence)
```

## Golden-reference results

| Query | Precursor | AF3 profile (1-based precursor coordinates) | 34 aa |
|---|---:|---|---|
| `HLA-A*02:01` | 365 aa | 25-308, 284 aa, `GSHSMRY...QPTIPI` | `YFAMYGEKVAHTHVDTLYVRYHYYTWAVLAYTWY` |
| `B2M-human` | 119 aa | 21-119, 99 aa, `IQRTPK...KWDRDM` | Not applicable |
| `H-2-Kb` | 369 aa | 22-305, 284 aa, `GPHSLR...STVSNM` | `YVEYYREKAGNSFVDTLYIVSQYYTWAELAYTWY` |
| `HLA-DRA*01:01` | 254 aa | 26-216, 191 aa | Combined with the next row |
| `HLA-DRB1*15:01` | 266 aa | 30-227, 198 aa | `QEFFIASGAAVDAIMWPRFDYFDIQAATYHVVFT` |
| `H-2-IAb` alpha/beta | 256/265 aa | 24-218 (195 aa) / 28-226 (199 aa) | `YSYFLASGGQVVHVLYFGYTYHDIRTETVHGPHT` |

The AF3 chain for HLA-A*02:01 ends with `...QPTIPI`; precursor position 309 is `V` and begins the UniProt-annotated TM segment at 309-332, so it must not be retained. The older 274-aa HLA-I crop instead omitted the membrane-proximal connector at positions 299-308 and is therefore not the default AF3 output.

## Regenerating the data

The builder reuses the local cache. Its first online run needs access to the pinned IPD-IMGT/HLA and PUFFIN sources and to UniProt:

```bash
python3 tools/build_reference.py
python3 tools/build_reference.py --offline
```

The builder does not inherit crops from the legacy `mhcii_model_input.csv`. The human HLA full registry is constructed by combining the matching-release FASTA and `Allele_status.txt`. AF3 chains use the pinned IPD protein alignments to project reviewed UniProt extracellular topology onto exact alleles, and every result must map uniquely back to one parent precursor. The pseudo artifact then normalizes legacy names through the official 2009-to-current nomenclature table and deleted/renamed tables. Mouse chains are retrieved by explicit accession. Mutable UniProt REST batches are pinned with expected SHA-256 values so changed content fails closed instead of silently changing a boundary. All generated files are stably sorted and record source hashes. Five conflicting records in the secondary pseudo mirror are not silently retained or resolved by choosing the first row; they are written to `quarantined_pseudosequences.csv` under an auditable rule.

## Tests

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

The 28 tests cover human/mouse, MHC-I/II, bidirectional full/pseudo/AF3 lookup, 12 golden AF3 boundaries, the parent-subsequence invariant across all 13,483 AF3 chains, MHC-II family pairing, B2M, valid `X` positions, official Partial fail-closed behavior, expression-suffix isolation, many-to-one ambiguity, pre-2010 HLA names, three-digit DPB1 first fields, compact-name collisions, and explicit projection status.

## Important boundaries

- `precursor_full` means that the source confirms coverage of the complete coding protein record, usually including signal peptide, extracellular domains, TM helix, and cytoplasmic tail. A true null/truncated allele may still be represented by a complete record. This is not the default AF3 input.
- AF3 structural input may use `af3-soluble-extracellular-v1` directly. Do not use a 34-aa pseudosequence as a foldable chain.
- P/G groups express identity over selected domains and must not be treated as proof of whole-precursor identity.
- The current pseudosequence artifact and latest full IPD snapshot have different coverage. A query without an exact pseudo record returns `not_found`.
- Source-data licensing and terms are not changed by this project's code license; see [`NOTICE.md`](NOTICE.md).

See [`AUDIT_REPORT.md`](AUDIT_REPORT.md) for the detailed audit, comparison with TCR tools, and AF3 five-chain recommendations.

## License and data terms

Package code is licensed under the [MIT License](LICENSE). Bundled and derived reference data retain their upstream attribution and applicable terms; the MIT License does not relicense data derived from IPD-IMGT/HLA, UniProt, or PUFFIN.
