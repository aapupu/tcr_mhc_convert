# TCR reference snapshots

The default V/J source is `imgt_tcr_alleles.csv`, downloaded from the official
[IMGT reference directory](https://www.imgt.org/vquest/refseqh.html) on 2026-11-25.
The adjacent metadata JSON contains all eight source URLs and their SHA-256
hashes. Each record retains its accession, functionality, raw translated
sequence, gapped sequence, and enabled/disabled status. There are 870 records;
684 are enabled for reconstruction. Missing/incomplete functional alleles fall
back to the gene representative; known P/ORF and stop-bearing records are
rejected. The gene tables identify their selected `reference_allele` explicitly.

CDR columns follow the tcrdist convention. Human TRA/TRB and mouse TRB use
1-based columns 27–38, 56–65, 81–86 and C anchor 104. Mouse TRA uses columns
28–39, 57–66, 82–88 and C anchor 106 because its alignment has additional
columns. These are alignment columns, not a redefinition of IMGT numbering.
CDR2.5 remains the tcrdist/HV4 convention.

Refresh with `tools/build_imgt_reference.py`; `--from-bundled` regenerates the
derived tables offline. The historical gene CSVs retain disabled records for
provenance; runtime also loads usable official alleles/genes not in those CSVs.

Historical tcrdist snapshot, retained for reproducibility and custom references:

```text
combo_xcr_2024-03-05.tsv
```

- Upstream: `kmayerb/tcrdist3`
- Git commit: `0e6efa39092e2c64c9c0b52a72c2f5961c161ebb`
- Upstream path: `tcrdist/db/combo_xcr_2024-03-05.tsv`
- SHA-256: `1d7fcac37ad6d955b3a3f693526e36623186464aef2c528cf8e0322dc89706cd`
- Included organisms used here: human and mouse
- Included chains used here: TRA (`A`) and TRB (`B`)

`alphabeta_gammadelta_db.tsv` is the older compatibility table referenced by
earlier tcrdist3 documentation. Default runtime lookups use the official IMGT
overlay; the snapshot above remains unchanged.

The underlying germline annotations originate from IMGT. Academic use should
cite IMGT and follow its data-use terms. See
https://github.com/kmayerb/tcrdist3/blob/master/docs/inputs.rst for the upstream
notice and database conventions.
