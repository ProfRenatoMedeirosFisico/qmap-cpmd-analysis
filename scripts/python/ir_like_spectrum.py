#!/usr/bin/env python3
"""Qualitative IR-like spectrum from the total-dipole autocorrelation function."""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from qmap_cpmd_common import read_dipole_csv,autocorrelation_fft,LIGHT_SPEED_CM_S

def main():
 p=argparse.ArgumentParser();p.add_argument('--dipole-csv',required=True);p.add_argument('--prefix',default='cluster');p.add_argument('--outdir',default='data/processed/figure09_ir_like');p.add_argument('--max-wavenumber-cm1',type=float,default=4000);p.add_argument('--window',choices=['none','hann'],default='hann');p.add_argument('--remove-linear-trend',action='store_true');args=p.parse_args()
 t,mu=read_dipole_csv(args.dipole_csv,args.prefix);dt_ps=float(np.median(np.diff(t)));x=mu.copy()
 if args.remove_linear_trend:
  u=np.arange(len(x));A=np.column_stack([u,np.ones(len(u))]);x=x-A@np.linalg.lstsq(A,x,rcond=None)[0]
 else:x-=x.mean(axis=0)
 ac_components=autocorrelation_fft(x);ac=ac_components.sum(axis=1) if ac_components.ndim==2 else ac_components
 ac_norm=ac/ac[0];lag=np.arange(len(ac))*dt_ps
 out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
 with (out/'dipole_autocorrelation.csv').open('w',newline='',encoding='utf-8') as f:w=csv.writer(f);w.writerow(['lag_ps','C_mu_mu_D2','normalized_C']);w.writerows(zip(lag,ac,ac_norm))
 window=np.hanning(len(ac)) if args.window=='hann' else np.ones(len(ac));signal=ac*window;nfft=1<<(max(2,len(signal)*4)-1).bit_length();spec=np.fft.rfft(signal,n=nfft).real;freq=np.fft.rfftfreq(nfft,d=dt_ps*1e-12);wn=freq/LIGHT_SPEED_CM_S;intensity=np.maximum(spec,0)*wn**2;mask=(wn>0)&(wn<=args.max_wavenumber_cm1)
 if np.max(intensity[mask])>0:intensity=intensity/np.max(intensity[mask])
 with (out/'ir_like_spectrum.csv').open('w',newline='',encoding='utf-8') as f:w=csv.writer(f);w.writerow(['wavenumber_cm-1','relative_intensity']);w.writerows(zip(wn[mask],intensity[mask]))
 summary={'samples':len(t),'sampling_interval_ps':dt_ps,'trajectory_span_ps':float(t[-1]-t[0]),'window':args.window,'definition':'I(w) proportional to w^2 Re[FFT(<mu(0).mu(t)>)]','interpretation':'qualitative finite-cluster descriptor'};(out/'ir_like_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
 try:
  import matplotlib.pyplot as plt
  plt.figure(figsize=(7,4));plt.plot(wn[mask],intensity[mask],lw=.8);plt.xlabel('Wavenumber (cm⁻¹)');plt.ylabel('Relative IR-like intensity');plt.tight_layout();plt.savefig(out/'ir_like_spectrum.png',dpi=300);plt.close()
 except ImportError:pass
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
