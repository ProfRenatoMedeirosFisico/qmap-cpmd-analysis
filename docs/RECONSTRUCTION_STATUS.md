# Reconstruction and validation status

The original analyses were developed interactively in ChatGPT and were not saved as a complete script set. The scripts in this repository are therefore **documented reimplementations**, reconstructed from the production input, reference geometry, manuscript definitions and surviving analysis script.

## Validated

- system composition and 241-atom ordering;
- QMAP atom selection and manuscript label mapping;
- assignment of all 21 DMSO molecules by S–O, S–C and C–H connectivity;
- PBC-aware geometry, pair-distance and H-bond algorithms;
- syntax and single-frame execution on `structures/GEOMETRY_reference.xyz`;
- FFT definitions for dipole autocorrelation, qualitative IR-like intensity, C1 and C2.

## Requires the local raw files

The full `TRAJEC.xyz` is approximately 700 MB and cannot be transferred through the connected Drive interface used for this reconstruction. Run the scripts locally with that trajectory to generate the complete CSV tables. The raw trajectory should be deposited in Zenodo or an institutional repository rather than committed to normal Git history.

## Dipole-dependent analyses

Figures 3–5, 9 and 10 require atomic charges by frame or the already calculated Cartesian dipole series. The Origin project alone is not treated as an open numerical format. Once `dipole_time_series.csv` is generated or exported, `ir_like_spectrum.py` and `orientational_correlations.py` reproduce the corresponding analyses without Origin.

No VMD/Tcl scripts were used. VMD served only for visual inspection and molecular representations.
