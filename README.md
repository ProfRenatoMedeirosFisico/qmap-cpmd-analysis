# QMAP–DMSO CPMD analysis

Reproducibility materials for **“Car-Parrinello Molecular Dynamics of QMAP-DMSO Microsolvation: First-Shell Structure and Solvent-Induced Polarization.”**

The production system is one QMAP plus 21 DMSO molecules (241 atoms) in a 16.0 Å cubic cell. The scripts were consolidated from analyses originally developed interactively and are explicitly documented as reimplementations.

## Included analyses

- QMAP internal coordinates used in Figure 2;
- finite-cluster O–O and N–S pair-distance distributions;
- directional O–H···O(DMSO) hydrogen bonds;
- QMAP and full-cluster dipoles from atomic charges;
- qualitative ω²-weighted dipole-autocorrelation spectrum;
- first- and second-rank orientational correlations.

No VMD/Tcl scripts were used; VMD was used only for visualization.

## Installation

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
```

## Atom mapping

The verified 1-based mapping is stored in `config/atom_indices.json`. It identifies QMAP as 31 atoms and reconstructs all 21 DMSO molecules from the production geometry.

## Coordinate-based analyses

```bash
python scripts/python/structural_analysis.py --trajectory /path/to/TRAJEC.xyz
python scripts/python/pair_distance_distributions.py --trajectory /path/to/TRAJEC.xyz
python scripts/python/hydrogen_bonds.py --trajectory /path/to/TRAJEC.xyz
```

The XYZ reader streams frames and ignores CPMD velocity columns, so the approximately 700 MB trajectory does not need to fit in memory.

## Dipole, IR-like and orientational analyses

```bash
python scripts/python/dipole_from_charges.py \
  --trajectory /path/to/TRAJEC.xyz \
  --charges /path/to/charges.txt

python scripts/python/ir_like_spectrum.py \
  --dipole-csv data/processed/figures03_05_dipoles/dipole_time_series.csv

python scripts/python/orientational_correlations.py \
  --dipole-csv data/processed/figures03_05_dipoles/dipole_time_series.csv
```

`charges.txt` may contain one charge per atom, `index charge` pairs, or a matrix with one frame per row and 241 atom columns.

## Important status

The scripts and atom mapping were validated on the real 241-atom reference geometry. Complete manuscript CSV tables require a local run on the full trajectory. Dipole-dependent figures additionally require the original charge series or Cartesian dipole series. The proprietary Origin project is not treated as the archival data format.

Large raw files are intentionally excluded from normal Git history. A permanent raw-data archive and DOI should be added before final submission.

## Citation and license

Citation metadata are provided in `CITATION.cff`. Code is released under the MIT License, subject to agreement of all authors.
