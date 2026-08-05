#!/usr/bin/env python3
"""Finite-cluster pair-distance distributions for O–O and N–S contacts."""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from qmap_cpmd_common import read_xyz_frames,load_mapping,validate_symbols,frame_time_ps,minimum_image


def hist_table(values,bins):
    count,edges=np.histogram(values,bins=bins); width=np.diff(edges); total=count.sum()
    pdf=count/(total*width) if total else np.zeros_like(width,float)
    prob=count/total if total else np.zeros_like(width,float)
    return (edges[:-1]+edges[1:])/2,count,pdf,prob

def main():
    p=argparse.ArgumentParser();p.add_argument('--trajectory',required=True);p.add_argument('--mapping',default='config/atom_indices.json');p.add_argument('--outdir',default='data/processed/figures06_07_pair_distances');p.add_argument('--timestep-au',type=float,default=5.0);p.add_argument('--saved-every-steps',type=int,default=5);p.add_argument('--frame-stride',type=int,default=1);p.add_argument('--bin-width-A',type=float,default=.05);p.add_argument('--rmax-A',type=float,default=10.0);p.add_argument('--write-raw',action='store_true');args=p.parse_args()
    m=load_mapping(args.mapping);cell=np.asarray([m['system']['cell_angstrom']]*3,float);a=m['analysis']['pair_distances'];out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    oq=a['O_QMAP']-1; od=np.asarray(a['O_DMSO'],int)-1; n1=a['N1_QMAP']-1;n2=a['N2_QMAP']-1;ss=np.asarray(a['S_DMSO'],int)-1
    data={'OQMAP_ODMSO':[],'N1QMAP_SDMSO':[],'N2QMAP_SDMSO':[],'NminQMAP_SDMSO':[]}
    raw=None
    if args.write_raw:
        raw=(out/'pair_distances_by_frame.csv').open('w',newline='',encoding='utf-8');wraw=csv.writer(raw);wraw.writerow(['frame','time_ps','pair','partner_id','distance_A'])
    nframes=0
    for fr in read_xyz_frames(args.trajectory):
        if fr.index%args.frame_stride:continue
        if fr.index==0:validate_symbols(fr.symbols,m)
        c=fr.coordinates;t=frame_time_ps(fr.index,args.timestep_au,args.saved_every_steps)
        do=np.linalg.norm(minimum_image(c[od]-c[oq],cell),axis=1)
        d1=np.linalg.norm(minimum_image(c[ss]-c[n1],cell),axis=1)
        d2=np.linalg.norm(minimum_image(c[ss]-c[n2],cell),axis=1)
        dm=np.minimum(d1,d2)
        for k,v in [('OQMAP_ODMSO',do),('N1QMAP_SDMSO',d1),('N2QMAP_SDMSO',d2),('NminQMAP_SDMSO',dm)]:data[k].extend(v.tolist())
        if raw:
            for pair,v in [('OQMAP_ODMSO',do),('N1QMAP_SDMSO',d1),('N2QMAP_SDMSO',d2),('NminQMAP_SDMSO',dm)]:
                for partner,x in enumerate(v,1):wraw.writerow([fr.index,t,pair,partner,x])
        nframes+=1
    if raw:raw.close()
    bins=np.arange(0,args.rmax_A+args.bin_width_A*1.001,args.bin_width_A)
    summary={'frames':nframes,'normalization':'probability density over sampled finite-cluster pair distances; not a bulk RDF'}
    for key,values in data.items():
        x=np.asarray(values);cent,count,pdf,prob=hist_table(x,bins)
        with (out/f'{key}_distribution.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.writer(f);w.writerow(['r_center_A','count','probability_density_A^-1','probability_per_bin']);w.writerows(zip(cent,count,pdf,prob))
        summary[key]={'samples':int(len(x)),'mean_A':float(x.mean()),'minimum_A':float(x.min()),'maximum_A':float(x.max())}
        try:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(6.5,4.5));plt.plot(cent,pdf,lw=1.2);plt.xlabel('r (Å)');plt.ylabel('Pair-distance probability density (Å⁻¹)');plt.tight_layout();plt.savefig(out/f'{key}_distribution.png',dpi=300);plt.close()
        except ImportError:pass
    (out/'pair_distance_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
