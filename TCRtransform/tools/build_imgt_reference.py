#!/usr/bin/env python3
"""Refresh the official IMGT V/J overlay; preserve raw sequence/functionality.

Run explicitly with network access. Generated tables are the offline runtime
reference. Only functional, complete-domain, standard-AA records are enabled.
"""
import concurrent.futures
import argparse
import csv
import hashlib
import json
import re
import runpy
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
CODONS = runpy.run_path(str(ROOT / 'tools/build_junction_reference.py'))['CODON_TABLE']
DATA = ROOT / 'src/tcrtransform/data'


def fetch(item):
    species, latin, group = item
    url = f'https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/{latin}/TR/{group}.fasta'
    raw = urllib.request.urlopen(url, timeout=90).read()
    records = []
    for chunk in raw.decode().split('>')[1:]:
        header, *lines = chunk.splitlines()
        fields = header.split('|')
        nt = ''.join(lines).upper()
        offset = int(fields[7]) - 1
        aligned = ''.join('.' if nt[i:i+3] == '...' else CODONS.get(nt[i:i+3], 'X') for i in range(offset, len(nt)-2, 3))
        aa = aligned.replace('.', '')
        region = group[-1]
        functionality = fields[3]
        anchor = -1
        loops = ['', '', '']
        bridge = codon = ''
        if region == 'V':
            # Mouse TRA alignments have additional columns; these are the
            # corresponding tcrdist CDR columns, not raw IMGT position numbers.
            mouse_alpha = species == 'mouse' and group == 'TRAV'
            c_column = 105 if mouse_alpha else 103
            spans = [(27,39),(56,66),(81,88)] if mouse_alpha else [(26,38),(55,65),(80,86)]
            if len(aligned) > c_column and aligned[c_column] == 'C':
                anchor = len(aligned[:c_column].replace('.', ''))
            loops = [aligned[a:b].replace('.', '') for a, b in spans]
        else:
            matches = list(re.finditer(r'(?=[FWC][A-Z]{2}G)', aa))
            if matches:
                anchor = matches[-1].start()
            coding = nt.replace('.', '')[offset:]
            remainder = coding[len(coding)//3*3:]
            if group == 'TRAJ' and len(remainder) == 1:
                codon = remainder + ('AT' if species == 'human' else 'AC')
                bridge = CODONS.get(codon, '')
        usable = functionality.strip('()[]') == 'F' and not (set(aa) - set('ACDEFGHIKLMNPQRSTVWY')) and anchor >= 0
        records.append(dict(species=species, allele=fields[1], gene=fields[1].split('*')[0], region=region,
            accession=fields[0], functionality=functionality, sequence=aa, aligned_sequence=aligned,
            usable=str(usable).lower(), cdr3_anchor=aa[anchor] if anchor >= 0 else '',
            cdr3_anchor_index=anchor, jc_bridge=bridge, jc_bridge_codon=codon,
            cdr1=loops[0], cdr2=loops[1], cdr2_5=loops[2]))
    return records, dict(url=url, sha256=hashlib.sha256(raw).hexdigest())


def update_gene_tables(rows):
    forced = {('human', 'TRAJ24'): 'TRAJ24*02'}
    for gene in ['TRAV13D-1', 'TRAV14-1', 'TRAV14D-1', 'TRAJ49']:
        forced['mouse', gene] = gene + '*01'
    for species in ['human', 'mouse']:
        path = DATA / f'tcr_gene_aa_{species}_TRA_TRB_gene.csv'
        with path.open() as handle:
            reader = csv.DictReader(handle)
            fields = list(reader.fieldnames)
            old = list(reader)
        for field in ['reference_allele', 'functionality', 'usable']:
            if field not in fields:
                fields.append(field)
        for row in old:
            gene = row['gene_name'].replace('-DV', '/DV')
            if not gene.startswith(('TRAV', 'TRBV', 'TRAJ', 'TRBJ')):
                row.update(reference_allele='', functionality='', usable='true')
                continue
            candidates = [r for r in rows if r['species'] == species and r['gene'].upper() == gene.upper()]
            if not candidates:
                row.update(reference_allele='', functionality='unrecognized', usable='false')
                continue
            eligible = [r for r in candidates if r['usable'] == 'true'] or candidates
            eligible.sort(key=lambda r: (r['allele'] != forced.get((species, gene), ''),
                r['sequence'] != row['amino_acid'], not r['allele'].endswith('*01'), r['allele']))
            chosen = eligible[0]
            row.update(amino_acid=chosen['sequence'], amino_acid_original=chosen['sequence'],
                af_trimmed_len=len(chosen['sequence']), af_trim_method='IMGT V/J region',
                reference_allele=chosen['allele'], functionality=chosen['functionality'], usable=chosen['usable'])
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator='\n')
            writer.writeheader(); writer.writerows(old)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--from-bundled', action='store_true', help='Rebuild derived tables from the existing downloaded snapshot')
    args = parser.parse_args()
    items = [(s,l,g) for s,l in [('human','Homo_sapiens'),('mouse','Mus_musculus')] for g in ['TRAV','TRAJ','TRBV','TRBJ']]
    rows, sources = [], []
    if args.from_bundled:
        with (DATA/'imgt_tcr_alleles.csv').open() as handle:
            rows = list(csv.DictReader(handle))
        sources = json.loads((DATA/'imgt_tcr_alleles.metadata.json').read_text())['sources']
    else:
        with concurrent.futures.ThreadPoolExecutor(8) as pool:
            for records, source in pool.map(fetch, items):
                rows.extend(records); sources.append(source)
    for row in rows:
        if row['region'] == 'J':
            # Overlapping motifs (e.g. FIFGSG) must select the downstream F.
            matches = list(re.finditer(r'(?=[FWC][A-Z]{2}G)', row['sequence']))
            anchor = matches[-1].start() if matches else -1
            row.update(cdr3_anchor_index=anchor, cdr3_anchor=row['sequence'][anchor] if anchor >= 0 else '')
        complete = row['region'] != 'V' or (not row['aligned_sequence'].startswith('.') and len(row['sequence']) >= 80)
        row['usable'] = str(row['usable'] == 'true' and complete).lower()
    with (DATA/'imgt_tcr_alleles.csv').open('w', newline='') as out:
        writer=csv.DictWriter(out, fieldnames=list(rows[0]), lineterminator='\n'); writer.writeheader(); writer.writerows(rows)
    update_gene_tables(rows)
    from tcrtransform.references import imgt_representatives
    junction_fields = ['species', 'gene_name', 'reference_allele', 'region', 'cdr3_anchor', 'cdr3_anchor_index', 'jc_bridge', 'jc_bridge_codon']
    with (DATA/'tcr_gene_junction_reference.csv').open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=junction_fields, lineterminator='\n')
        writer.writeheader()
        for species in ['human', 'mouse']:
            for gene, row in sorted(imgt_representatives(species).items()):
                record = {key: row[key] for key in junction_fields if key not in ['gene_name', 'reference_allele']}
                record.update(gene_name=gene, reference_allele=row['allele'])
                writer.writerow(record)
    retrieved_on = (json.loads((DATA/'imgt_tcr_alleles.metadata.json').read_text())['retrieved_on']
                    if args.from_bundled else datetime.now(timezone.utc).date().isoformat())
    metadata = dict(retrieved_on=retrieved_on, sources=sources,
        policy='Exact supported allele first; missing/incomplete alleles fall back to the gene representative. Runtime assesses geometry and retained stops after junction matching; P/ORF labels alone do not reject rearranged samples.',
        gene_policy='Preserve a functional exact existing representative where possible; otherwise prefer *01 then lowest allele. Human TRAJ24 uses *02; mouse TRAV13D-1, TRAV14-1, TRAV14D-1 and TRAJ49 use corrected *01.',
        usability='Historical strict subset only: F/(F)/[F], standard amino acids, conserved anchor; V N-terminus and at least 80 aa. Runtime also considers P/ORF and stop-bearing records with supported geometry.',
        citation='IMGT, the international ImMunoGeneTics information system; https://www.imgt.org',
        records=len(rows), usable_records=sum(r['usable']=='true' for r in rows))
    (DATA/'imgt_tcr_alleles.metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    (DATA/'tcr_gene_junction_reference.metadata.json').write_text(json.dumps(dict(
        source='imgt_tcr_alleles.csv', source_metadata='imgt_tcr_alleles.metadata.json',
        generated_by='tools/build_imgt_reference.py',
        note='Canonical gene rows; allele-specific junctions are loaded from the official allele table at runtime.'), indent=2)+'\n')
    for species in ['human', 'mouse']:
        table = DATA/f'tcr_gene_aa_{species}_TRA_TRB_gene.csv'
        provenance = dict(species=species, source_metadata='imgt_tcr_alleles.metadata.json',
            generated_by='tools/build_imgt_reference.py',
            table_sha256=hashlib.sha256(table.read_bytes()).hexdigest(),
            canonical_allele_rule=metadata['gene_policy'],
            selected_alleles={gene: row['allele'] for gene, row in sorted(imgt_representatives(species).items())},
            constant_region_sources='Unchanged model-ready constants inherited from tcr_design / pinned tcrdist build; mouse IMGT C-REGION display and UniProt P01851.',
            note='Runtime also includes official genes absent from the historical gene CSV; disabled raw rows in that CSV are not reconstruction candidates.')
        table.with_suffix('.metadata.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print(metadata['records'], metadata['usable_records'])


if __name__ == '__main__':
    main()
