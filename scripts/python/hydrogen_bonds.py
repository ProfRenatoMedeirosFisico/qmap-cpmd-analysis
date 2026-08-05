#!/usr/bin/env python3
"""Directional phenolic O–H···O(DMSO) hydrogen-bond analysis."""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from qmap_cpmd_common import read_xyz_frames,load_mapping,validate_symbols,frame_time_ps,minimum_image,angle_deg

def main():
    p=argparse.ArgumentParser();p.add_argument('--trajectory',required=True);p.add_argument('--mapping',default='config/atom_indices.json');p.add_argument('--outdir',default='data/processed/figure08_hydrogen_bonds');p.add_argument('--timestep-au',type=float,default=5.0);p.add_argument('--saved-every-steps',type=int,default=5);p.add_argument('--frame-stride',type=int,default=1);p.add_argument('--max-time-ps',type=float,default=12.5);p.add_argument('--h-acceptor-cutoff-A',type=float,default=2.5);p.add_argument('--angle-cutoff-deg',type=float,default=150.0);args=p.parse_args()
    m=load_mapping(args.mapping);cell=m['system']['cell_angstrom'];cfg=m['analysis']['hydrogen_bonds'];D=cfg['donor_O']-1;H=cfg['donor_H']-1;acc=np.asarray(cfg['acceptor_O_DMSO'],int)-1
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True);rows=[];active={i:None for i in range(1,len(acc)+1)};events=[]
    for fr in read_xyz_frames(args.trajectory):
        if fr.index%args.frame_stride:continue
        if fr.index==0:validate_symbols(fr.symbols,m)
        t=frame_time_ps(fr.index,args.timestep_au,args.saved_every_steps)
        if t>args.max_time_ps:break
        c=fr.coordinates;dist=np.linalg.norm(minimum_image(c[acc]-c[H],cell),axis=1);ang=np.asarray([angle_deg(c[D],c[H],c[a],cell) for a in acc]);present=(dist<=args.h_acceptor_cutoff_A)&(ang>=args.angle_cutoff_deg)
        partners=np.where(present)[0]+1;rows.append([fr.index,t,int(present.sum()),';'.join(map(str,partners)),float(dist.min()),float(ang[dist.argmin()])])
        for k in active:
            ison=bool(present[k-1])
            if ison and active[k] is None:active[k]=(fr.index,t)
            elif not ison and active[k] is not None:
                sf,st=active[k];events.append([k,sf,fr.index-args.frame_stride,st,t-frame_time_ps(args.frame_stride,args.timestep_au,args.saved_every_steps),t-st]);active[k]=None
    if rows:
        lastfr,lastt=rows[-1][0],rows[-1][1]
        dt=frame_time_ps(args.frame_stride,args.timestep_au,args.saved_every_steps)
        for k,v in active.items():
            if v is not None:sf,st=v;events.append([k,sf,lastfr,st,lastt,lastt-st+dt])
    with (out/'hydrogen_bonds_by_frame.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['frame','time_ps','hbond_count','dmso_acceptor_ids','minimum_H_A_A','angle_for_minimum_H_A_deg']);w.writerows(rows)
    with (out/'hydrogen_bond_events.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['dmso_id','start_frame','end_frame','start_ps','end_ps','continuous_duration_ps']);w.writerows(events)
    counts=np.asarray([r[2] for r in rows],int);dt=frame_time_ps(args.frame_stride,args.timestep_au,args.saved_every_steps)
    summary={'frames':len(rows),'window_ps':args.max_time_ps,'criteria':{'H_A_max_A':args.h_acceptor_cutoff_A,'D_H_A_min_deg':args.angle_cutoff_deg},'fraction_frames_bonded':float(np.mean(counts>0)) if len(counts) else None,'mean_count':float(np.mean(counts)) if len(counts) else None,'maximum_count':int(counts.max()) if len(counts) else None,'continuous_events':len(events),'sample_interval_ps':dt}
    (out/'hydrogen_bond_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    if rows:
        try:
            import matplotlib.pyplot as plt
            a=np.asarray(rows,object);plt.figure(figsize=(8,3.5));plt.step(a[:,1].astype(float),a[:,2].astype(float),where='post',lw=.7);plt.xlabel('Time (ps)');plt.ylabel('H-bond count');plt.ylim(-.05,max(1.05,counts.max()+.05));plt.tight_layout();plt.savefig(out/'hydrogen_bond_count.png',dpi=300);plt.close()
        except ImportError:pass
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
