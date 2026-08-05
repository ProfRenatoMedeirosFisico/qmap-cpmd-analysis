#!/usr/bin/env python3
"""First- and second-rank orientational correlations of the unit cluster dipole."""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from qmap_cpmd_common import read_dipole_csv,autocorrelation_fft

def main():
 p=argparse.ArgumentParser();p.add_argument('--dipole-csv',required=True);p.add_argument('--prefix',default='cluster');p.add_argument('--outdir',default='data/processed/figure10_orientational');p.add_argument('--max-lag-ps',type=float);p.add_argument('--integration-cutoff-ps',type=float);args=p.parse_args()
 t,mu=read_dipole_csv(args.dipole_csv,args.prefix);norm=np.linalg.norm(mu,axis=1);valid=norm>1e-14
 if not np.all(valid):t=t[valid];mu=mu[valid];norm=norm[valid]
 u=mu/norm[:,None];c1_comp=autocorrelation_fft(u);c1=c1_comp.sum(axis=1)
 tensors=np.einsum('ni,nj->nij',u,u).reshape(len(u),9);c2raw=autocorrelation_fft(tensors).sum(axis=1);c2=1.5*c2raw-.5
 c1/=c1[0];c2/=c2[0];dt=float(np.median(np.diff(t)));lag=np.arange(len(c1))*dt
 mask=np.ones(len(lag),bool) if args.max_lag_ps is None else lag<=args.max_lag_ps
 out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
 with (out/'orientational_correlations.csv').open('w',newline='',encoding='utf-8') as f:w=csv.writer(f);w.writerow(['lag_ps','C1','C2']);w.writerows(zip(lag[mask],c1[mask],c2[mask]))
 cutoff=args.integration_cutoff_ps
 if cutoff is None:
  crossings=np.where(c1<=0)[0];cutoff=float(lag[crossings[0]]) if len(crossings) else float(min(lag[-1],5.0))
 im=lag<=cutoff;tau=float(np.trapezoid(c1[im],lag[im]));summary={'samples':len(t),'sampling_interval_ps':dt,'integration_cutoff_ps':cutoff,'apparent_tau1_ps':tau,'interpretation':'exploratory finite-cluster decorrelation; not a bulk Debye time'};(out/'orientational_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
 try:
  import matplotlib.pyplot as plt
  plt.figure(figsize=(7,4));plt.plot(lag[mask],c1[mask],label='C1(t)');plt.plot(lag[mask],c2[mask],label='C2(t)');plt.xlabel('Lag time (ps)');plt.ylabel('Correlation');plt.legend();plt.tight_layout();plt.savefig(out/'orientational_correlations.png',dpi=300);plt.close()
 except ImportError:pass
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
