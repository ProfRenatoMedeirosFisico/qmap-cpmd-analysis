#!/usr/bin/env python3
"""Trajectory-resolved QMAP internal geometry used for manuscript Figure 2."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import numpy as np
from qmap_cpmd_common import (read_xyz_frames, load_mapping, validate_symbols,
    frame_time_ps, mic_distance, angle_deg, dihedral_deg)


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument('--trajectory',required=True)
    p.add_argument('--mapping',default='config/atom_indices.json')
    p.add_argument('--outdir',default='data/processed/figure02_structural')
    p.add_argument('--timestep-au',type=float,default=5.0)
    p.add_argument('--saved-every-steps',type=int,default=5)
    p.add_argument('--frame-stride',type=int,default=1)
    p.add_argument('--max-time-ps',type=float)
    p.add_argument('--transition-ps',type=float,default=12.5)
    args=p.parse_args()
    m=load_mapping(args.mapping); cell=m['system']['cell_angstrom']; d=m['analysis']['structural']
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for fr in read_xyz_frames(args.trajectory):
        if fr.index % args.frame_stride: continue
        if fr.index==0: validate_symbols(fr.symbols,m)
        t=frame_time_ps(fr.index,args.timestep_au,args.saved_every_steps)
        if args.max_time_ps is not None and t>args.max_time_ps: break
        c=fr.coordinates
        i,j=[x-1 for x in d['distance_C1_N1']]
        a,b,k=[x-1 for x in d['angle_C2_C1_N1']]
        q=[x-1 for x in d['dihedral_C2_C1_C3_C4']]
        raw=float(np.linalg.norm(c[j]-c[i]))
        corrected=mic_distance(c[i],c[j],cell)
        ang=angle_deg(c[a],c[b],c[k],cell)
        dih=dihedral_deg(c[q[0]],c[q[1]],c[q[2]],c[q[3]],cell)
        rows.append([fr.index,t,raw,corrected,ang,dih])
    if not rows: raise RuntimeError('No frames were selected')
    header=['frame','time_ps','C1_N1_raw_A','C1_N1_minimum_image_A','C2_C1_N1_deg','C2_C1_C3_C4_deg']
    with (out/'structural_time_series.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)
    arr=np.asarray(rows,float)
    summary={}
    for col,name in enumerate(header[2:],2):
        values=arr[:,col]; summary[name]={'mean':float(np.nanmean(values)),'std':float(np.nanstd(values,ddof=1)) if len(values)>1 else 0.0}
        for tag,mask in [('pre',arr[:,1]<args.transition_ps),('post',arr[:,1]>=args.transition_ps)]:
            if np.any(mask): summary[name][f'{tag}_mean']=float(np.nanmean(values[mask])); summary[name][f'{tag}_std']=float(np.nanstd(values[mask],ddof=1)) if np.count_nonzero(mask)>1 else 0.0
    (out/'structural_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    try:
        import matplotlib.pyplot as plt
        labels=[('C1_N1_raw_A','C1–N1 raw distance (Å)'),('C2_C1_N1_deg','C2–C1–N1 angle (°)'),('C2_C1_C3_C4_deg','C2–C1–C3–C4 dihedral (°)')]
        for key,ylabel in labels:
            idx=header.index(key); plt.figure(figsize=(8,4));plt.plot(arr[:,1],arr[:,idx],lw=.8);plt.axhline(np.nanmean(arr[:,idx]),ls='--',lw=1);plt.xlabel('Time (ps)');plt.ylabel(ylabel);plt.tight_layout();plt.savefig(out/f'{key}.png',dpi=300);plt.close()
    except ImportError: pass
    print(f'Wrote {len(rows)} frames to {out}')
if __name__=='__main__': main()
