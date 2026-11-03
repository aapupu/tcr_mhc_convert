# TCRtransform

**English** | [Simplified Chinese](README.zh-CN.md)

`TCRtransform` is part of the [`aapupu/tcr_mhc_convert`](https://github.com/aapupu/tcr_mhc_convert) monorepo.

`TCRtransform` is a standalone Python package for human and mouse alpha/beta T-cell receptor sequence transformation. It provides:

- reconstruction of model-ready extracellular TCR chains from `V gene + CDR3 + J gene + constant region`, including sequences suitable for AlphaFold 3 input;
- V-gene-based annotation of CDR1, CDR2, and CDR2.5;
- definition of the V-supported prefix, J-supported suffix, and intervening non-V/J-explained region within CDR3;
- reverse inference of V/J/C genes and CDR1/2/2.5/3 from complete TCR alpha/beta protein chains;
- atomic and reproducible CSV batch processing.

The package uses only the Python standard library and supports Python 3.9 or newer.

## Biological definitions

### CDR1, CDR2, and CDR2.5

CDR1, CDR2, and CDR2.5 are obtained from the bundled tcrdist3-aligned reference. CDR2.5 corresponds to the pMHC-facing/HV4 loop used by tcrdist3; it is not an additional classical IMGT CDR.

### CDR3

The package accepts and reports inclusive-anchor CDR3 amino-acid sequences. A CDR3 begins at the conserved V-side `C` and ends at the annotated J-side anchor. The J anchor is usually `F` or `W`, while explicitly annotated noncanonical anchors are also supported.

CDR3 is represented as:

```text
[V-supported prefix][protein-level non-V/J interval][J-supported suffix]
```

`CDR3*_non_vj` is only the amino-acid interval not covered by the accepted V and J reference matches. It does not by itself prove nucleotide-level N/P insertion. Codon phase, convergent amino acids, and beta-chain D-gene contributions cannot be resolved from a protein sequence alone. Exact nucleotide-level attribution requires the rearranged V(D)J nucleotide sequence and a dedicated rearrangement annotation method.

Both alpha and beta chains may have `non_vj_length = 0`. The package does not force an editable position by removing germline-supported residues.

## Reconstruction logic

A TCR chain is not built by directly concatenating `V + CDR3 + J + C`. The package instead:

1. finds the longest trimming-aware match between the CDR3 prefix and the V terminal window;
2. finds the longest trimming-aware match between the CDR3 suffix and the J start window;
3. removes duplicated residues and retains the CDR3 exactly once;
4. maps `TRBJ1-*` to `TRBC1` and `TRBJ2-*` to `TRBC2` for beta chains;
5. restores the TRAJ-TRAC split-codon bridge residue from the junction metadata for alpha chains;
6. appends the model-ready extracellular constant domain.

The default `min_overlap=3` is a strict reconstruction QC threshold. When it is explicitly reduced to 1 or 2, a short match is accepted only if all of the following conditions hold:

- the V-side match starts from the conserved C;
- the J-side match ends at the annotated anchor for that J gene;
- the match occurs at the expected V or J boundary;
- the admissible position is unique within the configured boundary window.

Rejection of a 1-2 aa match by the default threshold means that reconstruction evidence is insufficient. It is not a claim that the TCR is biologically invalid.

## Installation

```bash
cd TCRtransform
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

After installation, use either `tcrtransform` or `python -m tcrtransform`.

## CSV commands

### 1. Reconstruct complete chains and annotate all CDR regions

The input must contain at least `TRAV, TRAJ, CDR3A, TRBV, TRBJ, CDR3B`:

```bash
tcrtransform build examples/tcr_input.csv \
  --species-column species \
  -o examples/tcr_built.csv
```

The output preserves the input columns and appends `TCRa`, `TCRb`, CDR1/2/2.5, CDR3 boundary evidence, fixed masks, reference-resolution information, and status fields.

### 2. Annotate CDR regions without reconstructing full chains

```bash
tcrtransform annotate examples/tcr_input.csv \
  --species-column species \
  -o examples/tcr_annotated.csv
```

### 3. Parse complete TCR chains

The input should contain `TCRa` and `TCRb`; common aliases such as `A` and `B` are also accepted:

```bash
tcrtransform parse full_tcr.csv \
  --species human \
  -o full_tcr_parsed.csv
```

Use `--species-column species` for mixed human/mouse tables. By default, any row-level error prevents the final output from being written. Use `--on-error keep` to retain failed rows and record the reason in `tcrtransform_status`.

Use `--overwrite` when the output file already exists. Existing annotation columns are protected by default; add `--overwrite-columns` only when they should be recalculated.

## Key output columns

The alpha-chain names below have corresponding beta-chain fields ending in `b`:

| Column | Meaning |
|---|---|
| `CDR1a`, `CDR2a`, `CDR2.5a`, `CDR3a` | CDR amino-acid sequences |
| `CDR3a_v_supported` | CDR3 prefix supported by the V reference boundary |
| `CDR3a_non_vj` | Contiguous interval not covered by accepted V/J matches |
| `CDR3a_j_supported` | CDR3 suffix supported by the J reference boundary |
| `CDR3a_keep_prefix`, `CDR3a_keep_suffix` | Terminal lengths recommended for fixation during optimization |
| `CDR3a_fixed_mask` | Mask with one character per CDR3 residue: `1` is V/J supported and `0` is uncovered |
| `CDR3a_non_vj_start/end` | Zero-based, half-open candidate editable interval |
| `CDR3a_safe_for_optimization` | Whether both V and J boundaries were resolved |
| `CDR3a_boundary_status` | `ok`, `fully_vj_supported`, below-threshold short match, or unresolved status |
| `TCRa_reference_resolution` | Exact versus gene-level fallback used for V/J/C reconstruction |

Downstream optimization should require at least `safe_for_optimization=True` and should use `fixed_mask` directly. Do not force `min_editable > 0` by overwriting V/J-supported residues.

## Python API

```python
from tcrtransform import (
    annotate_v_gene,
    build_paired_tcr,
    define_cdr3_boundary,
    parse_chain,
)

loops = annotate_v_gene("human", "B", "TRBV20-1")
print(loops.cdr1, loops.cdr2, loops.cdr2_5)

boundary = define_cdr3_boundary(
    species="human",
    chain="B",
    v_gene="TRBV30",
    j_gene="TRBJ1-1",
    cdr3="CAWGRGARAFIDTEAFF",
)
print(boundary.non_vj_region, boundary.fixed_mask)

pair = build_paired_tcr(
    species="mouse",
    trav="TRAV14D-1",
    traj="TRAJ33",
    cdr3a="CAASDNYQLIW",
    trbv="TRBV12-1",
    trbj="TRBJ2-7",
    cdr3b="CASSRANYEQYF",
)
print(pair.alpha.sequence)
print(pair.beta.sequence)

parsed = parse_chain(pair.alpha.sequence, "A", species="mouse")
print(parsed.v_genes, parsed.j_genes, parsed.cdr3)
```

## Reference data and scope

- Supported species: `Homo sapiens` and `Mus musculus`.
- The default reconstruction references are canonical, gene-level V/J/C amino-acid sequences and are suitable for common 10x gene-level V/J annotations. These inputs cannot recover unreported individual allele polymorphisms.
- If an explicit allele is supplied, CDR1/2/2.5 annotation uses that allele. Full-chain reconstruction falls back to the same gene-level reference when no allele-level sequence is bundled, and reports this in `TCR*_reference_resolution`. Faithful cloning should use measured full-length sequences or a custom allele-level reference.
- The bundled constant sequences are extracellular structural domains for model input. They do not include leaders, transmembrane segments, or cytoplasmic tails. Therefore, "complete chain" in this project means a model-ready extracellular chain, not a complete expression-vector CDS.
- Reverse parsing retains tied V/J genes separated by `|` when protein evidence cannot distinguish them. It does not invent allele resolution that is absent from the sequence.
- Gamma/delta TCRs are not supported.

Method and definition references:

- [Stitchr, Nucleic Acids Research 2022](https://doi.org/10.1093/nar/gkac190)
- [IMGT FR/CDR definition](https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html)
- [tcrdist3 database conventions](https://github.com/kmayerb/tcrdist3/blob/master/docs/inputs.rst)

Snapshot versions, commits, checksums, and redistribution notes are documented in `src/tcrtransform/data/README.md` and `NOTICE.md`.

`tools/build_mouse_reference.py` and `tools/build_junction_reference.py` regenerate the mouse gene-level reference and the human/mouse junction metadata from the pinned tcrdist snapshot. Both tools have been verified to reproduce the bundled files byte for byte. Install the package first or set `PYTHONPATH=src` before running them.

## Tests

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test*.py' -v
```

The tests cover human and mouse CDR references, guarded short matches, fully V/J-supported CDR3s with zero editable residues, TRAJ-TRAC bridging, CSV processing, and forward/reverse round trips.

## Comparative validation

`validation/` contains the random comparison experiment performed before the core implementation was extracted into this package: 100 complete paired human TCRs and 100 complete paired mouse TCRs were compared against tcrdist3 and Stitchr. See [`validation/REPORT.md`](validation/REPORT.md) for the summary. Per-chain results, mismatch classifications, sampled receptors, and the machine-readable summary are also included for auditing.

## Project layout

```text
TCRtransform/
|-- pyproject.toml
|-- README.md
|-- README.zh-CN.md
|-- NOTICE.md
|-- examples/
|-- validation/
|-- tools/
|-- src/tcrtransform/
|   |-- cli.py
|   |-- junctions.py
|   |-- references.py
|   |-- regions.py
|   |-- reverse.py
|   |-- transform.py
|   `-- data/
`-- tests/
```

The package code is licensed under the [MIT License](LICENSE). Bundled and derived reference data retain their upstream attribution and applicable terms; the MIT License does not relicense IMGT/tcrdist-derived data.

