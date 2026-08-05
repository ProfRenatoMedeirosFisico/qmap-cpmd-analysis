#!/usr/bin/env python3
"""Run the available reconstruction workflow from a JSON configuration file."""
import argparse,json,subprocess,sys
from pathlib import Path

def run(script,args):
 cmd=[sys.executable,str(Path(__file__).with_name(script)),*args];print('+',' '.join(cmd));subprocess.run(cmd,check=True)
def main():
 p=argparse.ArgumentParser();p.add_argument('--config',default='config/reproduction_config.json');args=p.parse_args();cfg=json.loads(Path(args.config).read_text(encoding='utf-8'));traj=cfg.get('trajectory');mapping=cfg.get('mapping','config/atom_indices.json')
 if not traj:raise SystemExit('Set "trajectory" in the configuration file')
 common=['--trajectory',traj,'--mapping',mapping,'--timestep-au',str(cfg.get('timestep_au',5.0)),'--saved-every-steps',str(cfg.get('saved_every_steps',5))]
 run('structural_analysis.py',common);run('pair_distance_distributions.py',common);run('hydrogen_bonds.py',common)
 charges=cfg.get('charges')
 if charges:
  run('dipole_from_charges.py',[*common,'--charges',charges]);dip=cfg.get('dipole_csv','data/processed/figures03_05_dipoles/dipole_time_series.csv');run('ir_like_spectrum.py',['--dipole-csv',dip]);run('orientational_correlations.py',['--dipole-csv',dip])
 else:print('Dipole-dependent steps skipped: no charges file configured.')
if __name__=='__main__':main()
