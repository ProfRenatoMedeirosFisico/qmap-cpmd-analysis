# QMAP–DMSO CPMD analysis — manuscript release

**Tag:** `v1.0.0`

**Target branch:** `main`

## Release notes

This release provides the documented and consolidated Python analysis workflow associated with the manuscript:

**Car-Parrinello Molecular Dynamics of QMAP-DMSO Microsolvation: First-Shell Structure and Solvent-Induced Polarization**

### Included

- verified 1-based atom mapping for one QMAP and 21 DMSO molecules;
- structural analysis of QMAP internal coordinates;
- O(QMAP)–O(DMSO) and N(QMAP)–S(DMSO) pair-distance analyses;
- directional O–H···O(DMSO) hydrogen-bond analysis;
- QMAP and full-cluster dipole analysis from supplied charge data;
- qualitative dipole-autocorrelation/FFT spectrum;
- first- and second-rank orientational correlation functions;
- configuration files, reference geometry, validation data, tests, and figure-provenance documentation.

### Important provenance statement

The scripts were consolidated as standalone reimplementations of analyses originally developed interactively. They are not presented as byte-for-byte copies of the temporary interactive code.

### VMD

VMD was used for visual inspection and molecular representations only. No VMD/Tcl scripts were used in the quantitative analyses.

### Raw data

Large CPMD files are not attached to this software release. The production trajectory, native trajectory, restart file, energy data, legacy Origin project, and processed numerical tables are intended for a separate Zenodo Dataset record with its own DOI.

### License

Code: MIT License, subject to confirmation by all authors.
