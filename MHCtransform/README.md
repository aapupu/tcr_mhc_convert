# MHCtransform

**English** | [简体中文](README.zh-CN.md)

`MHCtransform` is an independent Python package in the [`aapupu/tcr_mhc_convert`](https://github.com/aapupu/tcr_mhc_convert) monorepo.

`MHCtransform` is an auditable, offline-capable human/mouse MHC conversion package covering:

- MHC-I heavy chains with species-matched B2M;
- MHC-II alpha/beta chains;
- deposited IPD/UniProt proteins, with records confirmed by official `Full` status or explicit UniProt completeness checks labeled `precursor_full`;
- mature soluble extracellular chains for AF3 complex input under `af3-soluble-extracellular-v1`;
- explicit partial, partial-guided imputed, and pseudo-only scaffold fallbacks without mixing them into the strict registry;
- two strictly separated 34-aa pseudosequence definitions;
- bidirectional lookup: strict sequences return **all exact candidates**, while explicitly included fallback constructs return `compatible_*` candidates.

It does not guess a unique allele from a pseudosequence, and it does not call a signal-peptide/TM-trimmed soluble chain full length.

## Representation definitions

| Representation/`definition_id` | Meaning | Primary use |
|---|---|---|
| `precursor_full` | Complete precursor in the source record: signal + extracellular + TM + cytoplasmic regions | Sequence verification and a starting point for expression-construct design; not the default AF3 chain |
| `af3-soluble-extracellular-v1` | Signal peptide, TM helix, and cytoplasmic tail removed; all extracellular folded domains and the membrane-proximal connector retained | AF3 pMHC-TCR multichain structural input |
| `source-partial-v1` | Exact deposited fragment; never relabeled as a full precursor | Evidence recovery and partial-sequence inspection |
| `precursor-imputed-v1` | Synthetic complete-length scaffold retaining aligned deposited residues and all non-`X` target pseudo residues; the remainder comes from a same-species, same-locus template | Explicit opt-in construct design, not an official allele sequence |
| `af3-soluble-extracellular-imputed-v1` | AF3-length construct under the same residue-level provenance rules | Explicit opt-in AF3 screening when no exact chain exists |
| `netmhcpan-34-v1` | 34 peptide-contact positions on a mature MHC-I heavy-chain reference alignment; excludes B2M | Pan-allele model descriptor, not a structural chain |
| `netmhciipan-34-v1` | 15 MHC-II alpha-chain positions plus 19 beta-chain positions concatenated from reference alignments; `X` is allowed | Pan-allele model descriptor, not a structural chain |

Position definitions and source hashes are stored in [`metadata.json`](src/mhctransform/data/metadata.json). MHC-II pseudosequences cannot be cut with one universal set of precursor indices because signal-peptide lengths and alignment indels vary. The runtime therefore uses a versioned published pseudosequence artifact instead of slicing full chains on demand.

Fallback is cumulative and disabled by default: `none` → `partial` → `imputed` → `pseudo-scaffold`. A 34mer is not invertible: 2,953 human MHC-I names collapse to 1,574 unique pseudos, while 5,620 human MHC-II names collapse to 1,160. Imputed outputs are therefore compatible structural scaffolds, not recovered biological full sequences.

## Current local registry

| Content | human I | human II | mouse I | mouse II |
|---|---:|---:|---:|---:|
| Deposited protein records | 20,066 | 9,419 | 12 | 17 |
| Complete precursor records available by default | 11,283 | 2,525 | 9 | 14 |
| Strict AF3-ready extracellular-chain records | 11,004 | 2,456 | 9 | 15 |
| 34-aa records | 2,953 | 5,620 | 7 | 8 |

Human HLA records come from pinned IPD-IMGT/HLA 3.66.0 at commit `b2055ea` and are evaluated together with the official `Full/Partial` status in the matching `Allele_status.txt`. Mouse MHC and human/mouse B2M records come from explicit UniProt accessions. AF3 boundaries are projected from reviewed UniProt topology onto the matching IPD protein alignments, and every output is verified as a contiguous subsequence of its parent precursor. The scientific 34-aa definitions come from the NetMHCpan/NetMHCIIpan literature; the distributed local artifact uses the public PUFFIN secondary mirror. Full source accessions, SHA-256 hashes, and generation dates are recorded in the metadata.

Deposited records include incomplete entries. By default, queries return human chains with official `Full` status and no unknown residues, or mouse chains that pass accession, fragment, and length checks. Use `--include-partial` to include fragments explicitly. The current mouse registry is a traceable set of commonly used H-2 alleles and does not claim exhaustive H-2 coverage.

Within the bundled 8,588-name pseudo registry, explicit fallback construction covers 100% of names: human I `2,953/2,953`, human II `5,620/5,620`, mouse I `7/7`, and mouse II `8/8`. Human expression-qualified `N/L/S/Q/C/A` targets remain blocked for AF3 unless `--allow-expression-qualified` is also supplied.

The isolated fallback assets contain 4,161 full/fragment rows and 2,097 AF3 rows. The latter comprise 67 exact deposited ectodomains, 1,990 partial-guided constructs, and 40 pseudo-only scaffolds.

Coverage by cumulative policy, counted over pseudo-registry names (`N/L/S/Q/C/A` AF3 targets explicitly enabled for this audit):

| Scope | Names | Strict full | `partial` full | `imputed` full | `pseudo-scaffold` full |
|---|---:|---:|---:|---:|---:|
| human I | 2,953 | 1,491 | 2,921 | 2,921 | 2,953 |
| human II | 5,620 | 2,454 | 5,619 | 5,619 | 5,620 |
| mouse I | 7 | 6 | 7 | 7 | 7 |
| mouse II | 8 | 5 | 8 | 8 | 8 |

| Scope | Names | Strict AF3 | `partial` AF3 | `imputed` AF3 | `pseudo-scaffold` AF3 |
|---|---:|---:|---:|---:|---:|
| human I | 2,953 | 1,491 | 1,551 | 2,921 | 2,953 |
| human II | 5,620 | 2,393 | 2,482 | 5,507 | 5,620 |
| mouse I | 7 | 6 | 6 | 7 | 7 |
| mouse II | 8 | 6 | 6 | 8 | 8 |

For AF3, `partial` accepts only a deposited sequence that covers the entire selected ectodomain; it never pads a short fragment. `imputed` requires usable deposited guidance and excludes pseudo-only constructs. Only the explicitly lowest-confidence `pseudo-scaffold` level reaches the full 34-aa registry.

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

# Exact deposited fragment only; representation=source-partial-v1
mhctransform allele-to-full 'HLA-B*08:03' --class I --fallback partial

# Complete-length synthetic construct guided by the deposited fragment
mhctransform allele-to-full 'HLA-B*08:03' --class I --fallback imputed
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

### Partial and imputed AF3 fallback

```bash
# Inspect every available representation without printing long sequences
mhctransform coverage 'HLA-B*08:03' --class I --species human

# Prefer a partial-guided construct; fails if only a pseudo-only scaffold exists
mhctransform molecule-to-af3 'HLA-B*08:03' \
  --class I --fallback imputed

# Permit the lowest-confidence pseudo-only scaffold when no deposited protein exists
mhctransform molecule-to-af3 'HLA-A*01:04' \
  --class I --fallback pseudo-scaffold

# MHC-II remains two separate chains
mhctransform molecule-to-af3 'H-2-IAq' \
  --class II --fallback pseudo-scaffold
```

Each fallback candidate reports `evidence_level`, deposited-observed/imputed spans and counts, source record IDs, template record/allele, pseudo distance, unresolved `X` positions, and quality flags. A deposited fragment that covers the complete extracellular target is promoted to `source_exact_ectodomain`; a shorter fragment is never passed directly to AF3 and instead requires `imputed` mode. Non-observed positions include template-derived and pseudo-constrained residues, with the pseudo source retained separately. `H-2-IAs` is not imputed: UniProt P14437 provides an exact 1–195 extracellular alpha chain, paired with the exact beta chain. `H-2-IAq` and `H-2-IAu` still require imputed alpha chains. Reverse lookup excludes fallback constructs unless `--include-fallback` is explicitly requested, and then reports construct compatibility rather than biological allele proof.

Representative results for the rare names in the audit set:

| Name | Deposited partial | AF3 output | Evidence/template |
|---|---:|---:|---|
| `HLA-B*08:03` | 181 aa | heavy 285 aa + B2M | partial-guided; `HLA-B*08:02` |
| `HLA-B*45:06` | 181 aa | heavy 285 aa + B2M | partial-guided; `HLA-B*45:01` |
| `HLA-A*02:50` | 181 aa | heavy 284 aa + B2M | partial-guided; `HLA-A*02:02` |
| `HLA-A*03:19` | 181 aa | heavy 284 aa + B2M | partial-guided; `HLA-A*03:01` |
| `H-2-Lq` | 206 aa | heavy 285 aa + mouse B2M | partial-guided; `H-2-Ld` |
| `H-2-IAs` | alpha 233 aa mature membrane record | alpha 195 aa + beta 197 aa | both strict exact; alpha from P14437 |
| `H-2-IAq` | alpha fragment | alpha 195 aa + beta 200 aa | alpha imputed from `H2-Aa*k`; beta exact |
| `H-2-IAu` | alpha fragment | alpha 195 aa + beta 197 aa | alpha imputed from `H2-Aa*k`; beta exact |

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
    "HLA-B*08:03",
    mhc_class="I",
    definition_id=AF3_DEFINITION,
    fallback="imputed",
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

# Rebuild only the deterministic fallback assets from existing strict assets
python3 tools/build_fallback_reference.py
```

The builder does not inherit crops from the legacy `mhcii_model_input.csv`. The human HLA full registry is constructed by combining the matching-release FASTA and `Allele_status.txt`. AF3 chains use the pinned IPD protein alignments to project reviewed UniProt extracellular topology onto exact alleles, and every result must map uniquely back to one parent precursor. The pseudo artifact then normalizes legacy names through the official 2009-to-current nomenclature table and deleted/renamed tables. Mouse chains are retrieved by explicit accession. Mutable UniProt REST batches are pinned with expected SHA-256 values so changed content fails closed instead of silently changing a boundary. All generated files are stably sorted and record source hashes. Five conflicting records in the secondary pseudo mirror are not silently retained or resolved by choosing the first row; they are written to `quarantined_pseudosequences.csv` under an auditable rule.

## Tests

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

The 41 tests cover the original strict behavior plus partial/imputed evidence separation, exact-ectodomain promotion, rare-HLA and H-2 fallbacks, P14437 exact I-A-s topology, expression-suffix and species gating, collision-safe opt-in fallback reverse lookup, and the invariants that all 8,588 pseudo names have AF3 and full fallbacks and every imputed full construct uses a complete precursor template. Strict AF3 chains still satisfy the parent-subsequence invariant; the one signal-free P14437 record is separately verified against its exact UniProt extracellular topology.

## Important boundaries

- `precursor_full` means that the source confirms coverage of the complete coding protein record, usually including signal peptide, extracellular domains, TM helix, and cytoplasmic tail. A true null/truncated allele may still be represented by a complete record. This is not the default AF3 input.
- AF3 structural input may use `af3-soluble-extracellular-v1` directly. Do not use a 34-aa pseudosequence as a foldable chain.
- `source-partial-v1` is a fragment, and `*-imputed-v1` is a synthetic construct. Neither is an official full allele sequence.
- A pseudo-only scaffold constrains 34 positions, not every peptide/TCR-facing polymorphism; use it for screening rather than allele-specific mechanistic interpretation.
- P/G groups express identity over selected domains and must not be treated as proof of whole-precursor identity.
- The current pseudosequence artifact and latest full IPD snapshot have different coverage. A query without an exact pseudo record returns `not_found`.
- Source-data licensing and terms are not changed by this project's code license; see [`NOTICE.md`](NOTICE.md).

See [`AUDIT_REPORT.md`](AUDIT_REPORT.md) for the detailed audit, comparison with TCR tools, and AF3 five-chain recommendations.

## License and data terms

Package code is licensed under the [MIT License](LICENSE). Bundled and derived reference data retain their upstream attribution and applicable terms; the MIT License does not relicense data derived from IPD-IMGT/HLA, UniProt, or PUFFIN.
