# Raw simulation data

The raw and processed data associated with this repository are intended to be archived separately in Zenodo.

## Dataset

**Title:** Raw and processed data for Car-Parrinello Molecular Dynamics of QMAP-DMSO Microsolvation

**DOI:** `PENDING — reserve in the Zenodo draft before submission`

## Main files

| File | Description |
|---|---|
| `TRAJEC.xyz.gz` | Converted XYZ trajectory used by the Python coordinate-analysis scripts |
| `TRAJECTORY.gz` | Native CPMD production trajectory |
| `RESTART.1.gz` | CPMD restart file |
| `ENERGIES.dat.gz` | Energy and dynamics monitoring data |
| `grafico08fev24.opj` | Legacy Origin project retained for provenance |
| `processed/*.csv` | Numerical data underlying manuscript figures |
| `atom_indices.json` | Verified atom mapping for one QMAP and 21 DMSO molecules |
| `SHA256SUMS.txt` | Integrity checksums for the archived files |

The production system contains one QMAP molecule and 21 DMSO molecules, totaling 241 atoms in a cubic cell with a side length of 16.0 Å.

## Files kept outside normal Git history

- `TRAJEC.xyz` — approximately 700 MB;
- `TRAJECTORY` — approximately 2.36 GB;
- `RESTART.1` — approximately 297 MB;
- `grafico08fev24.opj` — approximately 67 MB.

These files should not be committed to ordinary Git history. They should be uploaded to the separate Zenodo Dataset record described in `ZENODO_DATASET_METADATA.md`.

## Integrity verification

Generate SHA-256 checksums for all deposited files and include them in `SHA256SUMS.txt`.

## Generating open processed data

```bash
python scripts/python/structural_analysis.py --trajectory /path/TRAJEC.xyz
python scripts/python/pair_distance_distributions.py --trajectory /path/TRAJEC.xyz
python scripts/python/hydrogen_bonds.py --trajectory /path/TRAJEC.xyz
```

For dipole-dependent analyses:

```bash
python scripts/python/dipole_from_charges.py --trajectory /path/TRAJEC.xyz --charges /path/charges.txt
python scripts/python/ir_like_spectrum.py --dipole-csv data/processed/figures03_05_dipoles/dipole_time_series.csv
python scripts/python/orientational_correlations.py --dipole-csv data/processed/figures03_05_dipoles/dipole_time_series.csv
```

CSV is the archival format. The Origin project is retained only as a historical working file.

After the dataset DOI is reserved, replace the `PENDING` value above and update `docs/DATA_AVAILABILITY_TEXT.md`, `docs/REVIEWER_RESPONSE.md`, `README.md`, and the manuscript.
