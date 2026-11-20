# MHCtransform

**English** | [简体中文](README.zh-CN.md)

`MHCtransform` is an independent Python package in the [`aapupu/tcr_mhc_convert`](https://github.com/aapupu/tcr_mhc_convert) monorepo.

`MHCtransform` is an auditable, offline-capable human/mouse MHC conversion package covering:

- MHC-I heavy chains with species-matched B2M;
- MHC-II alpha/beta chains;
- deposited IPD/UniProt proteins, with records confirmed by official `Full` status or explicit UniProt completeness checks labeled `precursor_full`;
- mature soluble extracellular chains for AF3 complex input under `af3-soluble-extracellular-v1`;
- contiguous peptide-binding regions: MHC-I alpha1+alpha2 as one chain, or MHC-II alpha1 and beta1 as two separate chains;
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
| `mhci-alpha1-alpha2-v1` | Contiguous MHC-I alpha1+alpha2 peptide-binding groove; normally 182 aa; excludes alpha3 and B2M | Domain comparison and auditable extraction/checking of the 34 contact positions |
| `mhcii-alpha1-v1` | Contiguous MHC-II alpha1 peptide-binding domain | First of two separate MHC-II PBR chains |
| `mhcii-beta1-v1` | Contiguous MHC-II beta1 peptide-binding domain | Second of two separate MHC-II PBR chains |
| `netmhcpan-34-v1` | 34 peptide-contact positions on a mature MHC-I heavy-chain reference alignment; excludes B2M | Pan-allele model descriptor, not a structural chain |
| `netmhciipan-34-v1` | 15 MHC-II alpha-chain positions plus 19 beta-chain positions concatenated from reference alignments; `X` is allowed | Pan-allele model descriptor, not a structural chain |

PBR and pseudo are complementary but not interchangeable. A PBR is a contiguous foldable structural domain and contains every pseudo contact site; a 34mer is a non-contiguous model descriptor. The build has a hard pseudo-to-PBR contract: all 34 pseudo indices must map to PBR coordinates, every non-`X` residue must match exactly, and DTU `X` tokens are recorded only as alignment/unknown wildcards—they are never inserted as amino acids into a PBR chain. The shorter PBR is useful for comparison and pseudo verification, but the complete extracellular chain remains the default AF3 pMHC-TCR input because it retains the membrane-proximal domains and their native interfaces.

Position definitions and source hashes are stored in [`metadata.json`](src/mhctransform/data/metadata.json). MHC-II pseudo positions are calibrated alignment columns rather than universal precursor indices. Fallback is cumulative and disabled by default: `none` → `partial` → `imputed` → `pseudo-scaffold`. A 34mer remains many-to-many and cannot uniquely determine residues outside its contact coordinates. Imputed outputs are therefore compatible structural scaffolds, not recovered biological full sequences.

## Current local registry

| Content | human I | human II | mouse I | mouse II |
|---|---:|---:|---:|---:|
| Deposited protein records | 20,066 | 9,419 | 12 | 17 |
| Complete precursor records available by default | 11,283 | 2,525 | 9 | 14 |
| Strict AF3-ready extracellular-chain records | 11,004 | 2,456 | 9 | 15 |
| Strict PBR chain records | 11,003 | 2,456 | 8 | 15 |
| 34-aa records | 15,478 | 8,968 | 11 | 8 |

Human HLA records come from pinned IPD-IMGT/HLA 3.66.0 at commit `b2055ea` and are evaluated with the matching official `Full/Partial` status. Mouse MHC and human/mouse B2M records come from explicit UniProt accessions. AF3 and PBR boundaries are projected from reviewed UniProt topology/domain annotations onto the matching IPD protein alignments. Every strict AF3 and PBR output is checked as a contiguous subsequence of its parent precursor.

The 34-aa registry uses two distinct DTU resources. NetMHCpan 4.2f and NetMHCIIpan 4.3k public support-name snapshots define the comparison namespace. The bundled `NetMHCpan-4.2e-training-pseudoseqs.dat` and `NetMHCIIpan-4.3k-training-pseudosequence.2023.dat` files are the direct primary source for training pseudosequences. The `4.2e` label is the version of the supplied MHC-I pseudosequence artifact; it is not a claim that the separate 4.2f support list is a 4.2e list. In the human/mouse scope, 8,744 of 10,395 supported MHC-I names and 8,976 of 11,344 supported MHC-II molecules have a local 34mer.

The 24,465 bundled rows (15,489 MHC-I and 8,976 MHC-II) have explicit evidence: 5,848 direct DTU pseudosequences, 3,339 rows composed from DTU-certified MHC-II chain components, 10,808 exact PBR derivations, 3,844 contact-complete deposited partial derivations, and 626 legacy secondary rows. `H-2-Kq` is now supplied directly by the MHC-I DTU snapshot. The raw text snapshots and their SHA-256 digests are retained for reproducibility; no DTU executable or model is bundled.

Current-name handling is explicit: historical `HLA-DQA1*01:07` and `HLA-DRB1*15:13` inputs resolve to current `HLA-DQA1*01:07Q` and `HLA-DRB1*15:13Q`; compact variable-width fields such as `DRB1_13100` resolve to `HLA-DRB1*13:100`, not `*131:00`. The DTU pseudo is authoritative for the pseudo representation, but an incompatible strict deposited PBR is never silently rewritten. The strict record remains available as biological evidence, while `--require-pseudo-match` selects a separately labelled synthetic fallback when necessary. The component registry contains 4,390 chain pseudos; 344 rejected/conflicting rows remain quarantined.

`pseudo_components.csv` is an audit table of resolved chain components. DTU component rows take precedence for the pseudo representation, including four DPB1 components that conflict with current deposited derivations. Their strict deposited PBR records are preserved unchanged; pseudo-compatible structural fallbacks are separately labelled. The resulting MHC-II pair rows are model-supported chain combinations, not evidence of haplotype or in-vivo stability.

Deposited records include incomplete entries. By default, queries return human chains with official `Full` status and no unknown residues, or mouse chains that pass accession, fragment, and length checks. Use explicit fallback modes to inspect or construct from fragments. `partial` never pads a short sequence; `imputed` requires deposited guidance; `pseudo-scaffold` is the lowest-confidence structural fallback. Human `N/L/S/Q/C/A` expression-qualified targets remain blocked for AF3 unless explicitly enabled. The mouse registry remains a traceable set of selected H-2 accessions rather than an exhaustive IPD-MHC catalogue.

The current build contains 13,219 full, 6,628 AF3, and 6,628 PBR fallback records. Its hard contract covers all 18,080 unique pseudo chain components: 11,559 select a strict PBR and 6,521 select an explicit fallback, with zero uncovered components. Only 17 components require a bounded one-residue coordinate recovery; 12 of those contain an alignment-gap contact, and the maximum shift is 1 aa. These counts are generated into [`metadata.json`](src/mhctransform/data/metadata.json); they describe constructs and source fragments, not additional verified biological full sequences.

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

### Peptide-binding-region chains

```bash
# MHC-I: one contiguous alpha1+alpha2 heavy-chain region (normally 182 aa)
mhctransform allele-to-pbr 'HLA-A*02:01' --class I
mhctransform molecule-to-pbr 'HLA-A*02:01' --class I

# MHC-II: two independent chains, alpha1 and beta1
mhctransform molecule-to-pbr \
  'HLA-DRA*01:01/HLA-DRB1*15:01' --class II

# Return only the strict/fallback PBR certified against the DTU pseudo
mhctransform molecule-to-pbr \
  'HLA-DRA*01:01/HLA-DRB1*13:10' --class II \
  --require-pseudo-match

# Reverse lookup returns every exact chain candidate
mhctransform pbr-to-allele @pbr_chain.txt \
  --class I --definition mhci-alpha1-alpha2-v1
```

For MHC-I, `molecule-to-pbr` returns the heavy-chain alpha1+alpha2 groove and deliberately omits alpha3 and B2M. For MHC-II it returns alpha1 and beta1 as two sequence objects, never a concatenated artificial chain. Typical reviewed domain lengths are 182 aa for MHC-I, 84–96 aa for MHC-II alpha1, and 90–95 aa for MHC-II beta1; exact boundaries are locus- and allele-aware. `--require-pseudo-match` uses [`pseudo_pbr_projection.csv`](src/mhctransform/data/pseudo_pbr_projection.csv) to choose the certified chain for each component. A compatible strict PBR is preferred; if the DTU pseudo conflicts with a deposited strict chain, the latter remains unchanged and the command returns a separately labelled synthetic fallback. MHC-II alpha1 and beta1 are certified independently. Fallback PBRs retain observed/imputed spans and evidence levels in the same way as AF3 fallbacks.

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

For MHC-II, a 34mer represents one alpha/beta molecule: 15 calibrated alpha positions followed by 19 calibrated beta positions. A DQ or DP combination in the supported-name space is a computational model input, not evidence that the pair forms a common biological haplotype.

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

| Query | Precursor | AF3 profile (1-based precursor coordinates) | PBR | 34 aa |
|---|---:|---|---|---|
| `HLA-A*02:01` | 365 aa | 25-308, 284 aa, `GSHSMRY...QPTIPI` | alpha1+alpha2, 182 aa | `YFAMYGEKVAHTHVDTLYVRYHYYTWAVLAYTWY` |
| `B2M-human` | 119 aa | 21-119, 99 aa, `IQRTPK...KWDRDM` | Not applicable | Not applicable |
| `H-2-Kb` | 369 aa | 22-305, 284 aa, `GPHSLR...STVSNM` | alpha1+alpha2, 182 aa | `YVEYYREKAGNSFVDTLYIVSQYYTWAELAYTWY` |
| `HLA-DRA*01:01` | 254 aa | 26-216, 191 aa | alpha1, 84 aa | Combined with the next row |
| `HLA-DRB1*15:01` | 266 aa | 30-227, 198 aa | beta1, 95 aa | `QEFFIASGAAVDAIMWPRFDYFDIQAATYHVVFT` |
| `H-2-IAb` alpha/beta | 256/265 aa | 24-218 (195 aa) / 28-226 (199 aa) | separate alpha1/beta1 | `YSYFLASGGQVVHVLYFGYTYHDIRTETVHGPHT` |

The AF3 chain for HLA-A*02:01 ends with `...QPTIPI`; precursor position 309 is `V` and begins the UniProt-annotated TM segment at 309-332, so it must not be retained. The older 274-aa HLA-I crop instead omitted the membrane-proximal connector at positions 299-308 and is therefore not the default AF3 output.

## Regenerating the data

The core builder reuses the local cache. Its first online run needs access to the pinned IPD-IMGT/HLA sources and UniProt. The two DTU pseudosequence snapshots required for the release registry are already pinned under `tools/reference/` and are used by default:

```bash
python3 tools/build_reference.py
python3 tools/build_reference.py --offline

# Rebuild only the deterministic fallback assets from existing strict assets
# (this also rebuilds fallback_pbr_sequences.csv)
python3 tools/build_fallback_reference.py

# Rebuild strict and fallback PBR assets directly when needed
python3 tools/build_pbr_reference.py
```

The support-name and pseudosequence snapshots are committed under `tools/reference/` with pinned SHA-256 values. The release uses these files directly:

```bash
python3 tools/build_netmhc_pseudo_reference.py
python3 tools/build_fallback_reference.py
```

Explicit `--official-pseudo-i` and `--official-pseudo-ii` paths can override the bundled defaults for a controlled rebuild. The fallback builder refreshes fallback PBRs and the hard-checked `pseudo_pbr_projection.csv`; it fails if any pseudo component lacks all contact coordinates or any non-`X` residue disagrees with its selected PBR. The builder does not inherit crops from the legacy `mhcii_model_input.csv`. Human full records combine the matching-release FASTA and status file. AF3/PBR chains project reviewed UniProt features through pinned IPD alignments and must map uniquely to a parent. MHC-I pseudo contact positions and the MHC-II alpha15+beta19 alignment mappings are versioned. Generated outputs and conflict reports are stably sorted and deterministic; `metadata.json` and `pseudo_validation.json` record hashes of source inputs, not hashes of every output file.

## Tests

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

Unit tests cover strict and fallback evidence separation, rare-HLA/H-2 handling, expression-suffix and species gating, collision-safe reverse lookup, complete-template requirements, PBR chain pairing/boundaries, DTU primary-source parsing, and the pseudo-to-PBR hard contract. Strict AF3 and PBR chains satisfy parent-subsequence invariants; the signal-free P14437 record is checked separately against its reviewed UniProt extracellular topology. Deterministic output was checked separately during release validation by repeating complete rebuilds and comparing their output hashes.

## Important boundaries

- `precursor_full` means that the source confirms coverage of the complete coding protein record, usually including signal peptide, extracellular domains, TM helix, and cytoplasmic tail. A true null/truncated allele may still be represented by a complete record. This is not the default AF3 input.
- AF3 structural input may use `af3-soluble-extracellular-v1` directly. Do not use a 34-aa pseudosequence as a foldable chain.
- PBR is a contiguous structural region, but it is intentionally shorter than the recommended full extracellular AF3 chain and omits stabilizing membrane-proximal domains (and B2M for MHC-I).
- `source-partial-v1` is a fragment, and `*-imputed-v1` is a synthetic construct. Neither is an official full allele sequence.
- A pseudo-only scaffold constrains 34 positions, not every peptide/TCR-facing polymorphism; use it for screening rather than allele-specific mechanistic interpretation.
- P/G groups express identity over selected domains and must not be treated as proof of whole-precursor identity.
- NetMHC support membership, availability of a validated 34mer, and availability of a full/PBR protein are different facts. A query without a safe pseudo source returns `not_found` even when the name appears in an official support list.
- Source-data licensing and terms are not changed by this project's code license; see [`NOTICE.md`](NOTICE.md).

See [`AUDIT_REPORT.md`](AUDIT_REPORT.md) for the detailed audit, comparison with TCR tools, and AF3 five-chain recommendations.

## License and data terms

Package code is licensed under the [MIT License](LICENSE). Bundled and derived reference data retain upstream attribution and applicable terms; MIT does not relicense data from IPD-IMGT/HLA, UniProt, DTU, or PUFFIN. The two DTU pseudosequence text snapshots are bundled as direct data inputs; NetMHCpan/NetMHCIIpan executables, model files, and complete training archives are not distributed by this repository.
