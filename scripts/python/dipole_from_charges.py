#!/usr/bin/env python3
"""Compute QMAP and finite-cluster dipoles from trajectory coordinates and charges.

The script supports one static charge per atom (or ``index charge`` pairs). For
frame-dependent charges, first convert them to a matrix with one frame per row
and one atom per column. Molecular coordinates are reconstructed with PBC-aware
connectivity before each neutral-molecule dipole is evaluated.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from qmap_cpmd_common import (read_xyz_frames,load_mapping,validate_symbols,frame_time_ps,
    unwrap_connected,DEBYE_PER_E_ANGSTROM)

def load_charges(path,natoms):
    lines=[x.strip() for x in Path(path).read_text(encoding='utf-8').splitlines() if x.strip() and not x.lstrip().startswith('#')]
    rows=[]
    for line in lines:
        try: rows.append([float(v) for v in line.replace(',',' ').split()])
        except ValueError: continue
    if len(rows)==natoms and all(len(r)==1 for r in rows):return np.asarray([r[0] for r in rows]),None
    if len(rows)==natoms and all(len(r)>=2 for r in rows):
        q=np.zeros(natoms);one=min(int(r[0]) for r in rows)==1
        for r in rows:q[int(r[0])-(1 if one else 0)]=r[1]
        return q,None
    arr=np.asarray(rows,float)
    if arr.ndim==2 and arr.shape[1]==natoms:return None,arr
    raise ValueError(f'Cannot interpret charges: expected {natoms} static rows or a frame x {natoms} matrix')

def mol_dipole(coords,q):
    origin=coords.mean(axis=0);res=float(q.sum());mu=np.sum((coords-origin)*q[:,None],axis=0)*DEBYE_PER_E_ANGSTROM
    return mu,res

def main():
    p=argparse.ArgumentParser();p.add_argument('--trajectory',required=True);p.add_argument('--charges',required=True);p.add_argument('--mapping',default='config/atom_indices.json');p.add_argument('--outdir',default='data/processed/figures03_05_dipoles');p.add_argument('--timestep-au',type=float,default=5.0);p.add_argument('--saved-every-steps',type=int,default=5);p.add_argument('--frame-stride',type=int,default=1);args=p.parse_args()
    m=load_mapping(args.mapping);nat=m['system']['atom_count'];static,matrix=load_charges(args.charges,nat);cell=m['system']['cell_angstrom'];out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True);rows=[]
    qids=m['qmap']['atom_indices'];qedges=m['qmap']['connectivity'];mols=[(qids,qedges)]
    for d in m['dmso_molecules']:
        ids=[d['S'],d['O'],*d['C'],*d['H']];edges=[[d['S'],d['O']]]+[[d['S'],c] for c in d['C']]
        for c,hs in zip(d['C'],[d['H'][:3],d['H'][3:]]):edges += [[c,h] for h in hs]
        mols.append((ids,edges))
    for fr in read_xyz_frames(args.trajectory):
        if fr.index%args.frame_stride:continue
        if fr.index==0:validate_symbols(fr.symbols,m)
        q=static if static is not None else matrix[fr.index]
        mus=[];res=[]
        for ids,edges in mols:
            uc=unwrap_connected(fr.coordinates,ids,edges,cell);qq=q[np.asarray(ids)-1];mu,r=mol_dipole(uc,qq);mus.append(mu);res.append(r)
        qm=mus[0];total=np.sum(mus,axis=0);t=frame_time_ps(fr.index,args.timestep_au,args.saved_every_steps)
        rows.append([fr.index,t,*qm,float(np.linalg.norm(qm)),*total,float(np.linalg.norm(total)),res[0],sum(res)])
    header=['frame','time_ps','qmap_mu_x_D','qmap_mu_y_D','qmap_mu_z_D','qmap_mu_D','cluster_mu_x_D','cluster_mu_y_D','cluster_mu_z_D','cluster_mu_D','qmap_charge_e','cluster_charge_e']
    with (out/'dipole_time_series.csv').open('w',newline='',encoding='utf-8') as f:w=csv.writer(f);w.writerow(header);w.writerows(rows)
    a=np.asarray(rows,float);summary={k:float(v) for k,v in zip(['qmap_mean_D','cluster_mean_D','qmap_std_D','cluster_std_D'],[a[:,5].mean(),a[:,9].mean(),a[:,5].std(ddof=1) if len(a)>1 else 0.0,a[:,9].std(ddof=1) if len(a)>1 else 0.0])}
    (out/'dipole_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8,4));plt.plot(a[:,1],a[:,5],label='QMAP');plt.plot(a[:,1],a[:,9],label='Full cluster');plt.xlabel('Time (ps)');plt.ylabel('|μ| (D)');plt.legend();plt.tight_layout();plt.savefig(out/'dipole_comparison.png',dpi=300);plt.close()
    except ImportError:pass
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
