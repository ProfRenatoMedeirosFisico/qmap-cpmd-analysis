# Repository completion checklist

Do not cite the repository as a complete reproducibility deposit until every required item below is resolved.

## Scientific files

- [x] Add `analisar_qmap_cpmd_correlacoes.py`.
- [ ] PENDING: add the exact structural-analysis script(s) used for manuscript Figure 2.
- [ ] PENDING: add the exact full-cluster and QMAP-only dipole script(s) used for Figures 3–5.
- [ ] PENDING: add the exact pair-distance distribution script(s) used for Figures 6–7.
- [ ] PENDING: add the exact hydrogen-bond script used for Figure 8.
- [ ] PENDING: add the exact IR-like autocorrelation/FFT script used for Figure 9.
- [ ] PENDING: add the exact first- and second-rank orientational-correlation script used for Figure 10.
- [ ] PENDING: add all VMD/Tcl scripts actually used.
- [ ] PENDING: add the processed numerical table underlying each figure.
- [ ] PENDING: add a figure-to-script/data map in `docs/FIGURE_PROVENANCE.md`.

## Verification

- [ ] Resolve the manuscript inconsistency concerning 10 versus 21 DMSO molecules.
- [ ] Verify the total atom count against the production trajectory.
- [ ] Record the exact atom order and indices in `docs/ATOM_SELECTIONS.md`.
- [ ] Record trajectory time step, saved-frame interval, and conversion from frame to ps.
- [ ] Verify whether periodic-boundary unwrapping was applied and document the procedure.
- [ ] Verify charge source, charge units, and whether charges vary by frame.
- [ ] Re-run all scripts in a clean environment.
- [ ] Compare regenerated figures/tables with the submitted manuscript.
- [ ] Remove confidential, proprietary, or unnecessary files.
- [ ] Confirm the code license with all authors.
- [ ] Create a numbered GitHub release.
- [ ] Archive the release and obtain a DOI, when appropriate.
- [ ] Replace the placeholder repository URL and DOI in the manuscript and response letter.
