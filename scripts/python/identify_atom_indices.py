#!/usr/bin/env python3
"""Verify and document the atom mapping of the 241-atom QMAP–DMSO system."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import numpy as np
from qmap_cpmd_common import first_xyz_frame, mic_distance, validate_symbols


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--geometry', required=True)
    p.add_argument('--mapping', default='config/atom_indices.json')
    p.add_argument('--outdir', default='data/processed/atom_mapping')
    args = p.parse_args()
    mapping = json.loads(Path(args.mapping).read_text(encoding='utf-8'))
    frame = first_xyz_frame(args.geometry)
    validate_symbols(frame.symbols, mapping)
    cell = mapping['system']['cell_angstrom']
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)

    qmap = mapping['qmap']
    rows=[]
    for label, idx in qmap['labels'].items():
        rows.append([label, idx, frame.symbols[idx-1], *frame.coordinates[idx-1]])
    with (out/'qmap_labeled_atoms.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['label','index_1based','element','x_A','y_A','z_A']); w.writerows(rows)

    dmso_rows=[]
    for mol in mapping['dmso_molecules']:
        s,o=mol['S'],mol['O']
        cs=mol['C']; hs=mol['H']
        dmso_rows.append([mol['id'],s,o,';'.join(map(str,cs)),';'.join(map(str,hs)),
                          mic_distance(frame.coordinates[s-1],frame.coordinates[o-1],cell),
                          mic_distance(frame.coordinates[s-1],frame.coordinates[cs[0]-1],cell),
                          mic_distance(frame.coordinates[s-1],frame.coordinates[cs[1]-1],cell)])
    with (out/'dmso_molecule_mapping.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['dmso_id','S','O','C_indices','H_indices','S_O_A','S_C1_A','S_C2_A']); w.writerows(dmso_rows)

    donor_o=qmap['labels']['O1']; donor_h=qmap['labels']['H1']
    report={
        'atom_count':len(frame.symbols),
        'qmap_atom_count':len(qmap['atom_indices']),
        'dmso_count':len(mapping['dmso_molecules']),
        'phenolic_O_H_distance_A':mic_distance(frame.coordinates[donor_o-1],frame.coordinates[donor_h-1],cell),
        'status':'PASS'
    }
    (out/'mapping_validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
