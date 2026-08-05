# Zenodo dataset metadata

Use the fields below when creating the separate Zenodo **Dataset** record.

## Basic information

**Resource type:** Dataset

**Title:**

Raw and processed data for Car-Parrinello Molecular Dynamics of QMAP-DMSO Microsolvation

**Creators, in order:**

1. Medeiros, Renato
2. Osório, Francisco A. P.
3. Valverde, Clodoaldo
4. Camargo, Ademir J.

Add ORCID identifiers and affiliations in the Zenodo form when they are available and verified.

**Description:**

This dataset supports the manuscript “Car-Parrinello Molecular Dynamics of QMAP-DMSO Microsolvation: First-Shell Structure and Solvent-Induced Polarization”. It contains the converted XYZ production trajectory, the native CPMD trajectory and restart file, simulation energy data, the legacy Origin project used during data inspection, simulation inputs, representative geometries, atom-index mapping, checksums, and processed numerical tables associated with the manuscript analyses. The production system contains one QMAP molecule and 21 DMSO molecules, totaling 241 atoms in a cubic cell with a side length of 16.0 Å. The accompanying Python analysis software is maintained separately in the GitHub repository `ProfRenatoMedeirosFisico/qmap-cpmd-analysis` and will receive an independent software DOI through the Zenodo–GitHub integration.

**Keywords:**

- Car-Parrinello molecular dynamics
- microsolvation
- QMAP
- DMSO
- hydrogen bonding
- dipole moment
- orientational correlation
- molecular dynamics
- reproducible research

**Language:** English

**Access right:** Open access

**Recommended data license:** Creative Commons Attribution 4.0 International (CC BY 4.0), subject to agreement of all authors.

## DOI

In the Zenodo field asking whether the upload already has a DOI:

1. select **No**;
2. click **Get a DOI now!**;
3. copy the reserved DOI;
4. do not delete the draft, because deleting it also discards the reserved DOI.

**Reserved dataset DOI:** `PENDING`

## Files to upload

| File | Description |
|---|---|
| `TRAJEC.xyz.gz` | Converted XYZ trajectory used by the coordinate-analysis scripts |
| `TRAJECTORY.gz` | Native CPMD production trajectory |
| `RESTART.1.gz` | CPMD restart file |
| `ENERGIES.dat.gz` | Energy and dynamics monitoring data |
| `grafico08fev24.opj` | Legacy Origin project retained for provenance |
| `qmap241-wfnop.inp` | Wave-function optimization input |
| `qmap241-eq.inp` | Equilibration input |
| `qmap241-sim.inp` | Production input |
| `GEOMETRY.xyz` | Representative/final geometry |
| `atom_indices.json` | Verified 1-based atom mapping |
| `processed/*.csv` | Numerical data underlying manuscript figures |
| `README_DATASET.md` | Dataset description and usage notes |
| `SHA256SUMS.txt` | Integrity checksums |

## Related identifiers to add later

- Article DOI: relation **Is supplement to** — `PENDING`
- Software DOI: relation **Is supplemented by** — `PENDING`
- GitHub repository URL: relation **Is supplemented by**

## Draft status

Keep the record as a draft until:

- all files have finished uploading;
- SHA-256 checksums have been checked;
- the processed CSV files have been generated and inspected;
- author order, ORCIDs, affiliations, license, and description have been verified;
- the reserved DOI has been inserted into `docs/RAW_DATA.md`.
