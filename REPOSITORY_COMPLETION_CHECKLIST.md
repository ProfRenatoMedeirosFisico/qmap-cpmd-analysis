# Repository completion checklist

## Completed

- [x] Reimplement the structural analysis for Figure 2.
- [x] Reimplement QMAP and cluster dipoles from atomic charges.
- [x] Reimplement O–O and N–S pair-distance distributions.
- [x] Reimplement directional hydrogen-bond counting.
- [x] Reimplement the dipole-autocorrelation IR-like spectrum.
- [x] Reimplement first- and second-rank orientational correlations.
- [x] Verify one QMAP + 21 DMSO = 241 atoms.
- [x] Record exact atom ordering and selections.
- [x] Document the 5 a.u. time step and sampling every 5 MD steps.
- [x] Add PBC-aware molecular reconstruction and minimum-image distances.
- [x] State that no VMD/Tcl scripts were used.
- [x] Validate syntax and single-frame execution on the real reference geometry.

## Required before claiming full figure-level reproduction

- [ ] Run the coordinate-based scripts locally on the complete `TRAJEC.xyz`.
- [ ] Recover or regenerate the frame-dependent atomic charges or Cartesian dipole series.
- [ ] Generate and deposit the complete numerical CSV tables underlying Figures 2–10.
- [ ] Compare the regenerated figures and numerical means with the manuscript.
- [ ] Archive the large raw trajectory and restart files in a research-data repository.
- [ ] Insert the raw-data DOI in the manuscript and repository.
- [ ] Confirm the code license with all authors.
- [ ] Create a numbered GitHub release after final validation.
