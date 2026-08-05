#!/usr/bin/env python3
"""Extract coordinates and atomic charges from repeated ``XYZ= ... Q= ...`` log blocks.

This parser is provided for logs created by the interactive charge workflow used
in the project. It writes a frame-by-atom charge matrix accepted by
``dipole_from_charges.py`` and reports incomplete blocks instead of guessing.
"""
import argparse,re,csv
from pathlib import Path
pat=re.compile(r'(?:XYZ=\s*)?([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+Q=\s*([-+0-9.Ee]+)')
def main():
 p=argparse.ArgumentParser();p.add_argument('--log',required=True);p.add_argument('--atoms',type=int,default=241);p.add_argument('--out',default='charge_frames.csv');args=p.parse_args();values=[]
 for line in Path(args.log).open('r',encoding='utf-8',errors='ignore'):
  m=pat.search(line)
  if m:values.append(float(m.group(4)))
 complete=len(values)//args.atoms;rem=len(values)%args.atoms
 if complete==0:raise RuntimeError('No complete charge frame found')
 with Path(args.out).open('w',newline='',encoding='utf-8') as f:
  w=csv.writer(f);w.writerow([f'atom_{i}' for i in range(1,args.atoms+1)])
  for k in range(complete):w.writerow(values[k*args.atoms:(k+1)*args.atoms])
 print(f'frames={complete}, ignored_values={rem}, output={args.out}')
if __name__=='__main__':main()
