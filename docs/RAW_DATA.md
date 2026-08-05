# Raw and processed data

## Files kept outside normal Git history

- `TRAJEC.xyz` — converted coordinate trajectory, about 700 MB;
- `TRAJECTORY` — original CPMD trajectory, about 2.36 GB;
- `RESTART.1` — CPMD restart, about 297 MB;
- `grafico08fev24.opj` — proprietary Origin project, about 67 MB.

The raw trajectories should be archived in a repository designed for large research data. The archive DOI and checksums should then be added here.

## Generating open processed data

```bash
python scripts/python/structural_analysis.py --trajectory /path/TRAJEC.xyz
python scripts/python/pair_distance_distributions.py --trajectory /path/TRAJEC.xyz
python scripts/python/hydrogen_bonds.py --trajectory /path/TRAJEC.xyz
```

For dipoles:

```bash
python scripts/python/dipole_from_charges.py --trajectory /path/TRAJEC.xyz --charges /path/charges.txt
python scripts/python/ir_like_spectrum.py --dipole-csv data/processed/figures03_05_dipoles/dipole_time_series.csv
python scripts/python/orientational_correlations.py --dipole-csv data/processed/figures03_05_dipoles/dipole_time_series.csv
```

CSV is the archival format; the Origin project is retained only as a historical working file.
