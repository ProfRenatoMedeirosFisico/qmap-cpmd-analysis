# Figure provenance

The scripts below are documented reimplementations of analyses originally developed interactively.

| Figure | Quantity | Script | Main output | Current status |
|---|---|---|---|---|
| 2 | C1–N1 distance, C2–C1–N1 angle and C2–C1–C3–C4 torsion | `structural_analysis.py` | `structural_time_series.csv` | validated on reference geometry; full trajectory run pending |
| 3 | QMAP and cluster dipole magnitudes | `dipole_from_charges.py` | `dipole_time_series.csv` | requires original charge series |
| 4 | Cluster Cartesian dipole components | `dipole_from_charges.py` | `dipole_time_series.csv` | requires original charge series |
| 5 | QMAP dipole magnitude | `dipole_from_charges.py` | `dipole_time_series.csv` | requires original charge series |
| 6 | O(QMAP)–O(DMSO) pair-distance distribution | `pair_distance_distributions.py` | `OQMAP_ODMSO_distribution.csv` | algorithm validated; full trajectory run pending |
| 7 | N(QMAP)–S(DMSO) pair-distance distributions | `pair_distance_distributions.py` | N1, N2 and nearest-N CSV files | exact original N selection not preserved; all defensible variants exported |
| 8 | Phenolic O–H···O(DMSO) count and events | `hydrogen_bonds.py` | `hydrogen_bonds_by_frame.csv` and `hydrogen_bond_events.csv` | criteria validated; full trajectory run pending |
| 9 | Qualitative ω²-weighted dipole-autocorrelation spectrum | `ir_like_spectrum.py` | `dipole_autocorrelation.csv`, `ir_like_spectrum.csv` | code validated with synthetic vectors; original dipole series required |
| 10 | First- and second-rank dipole orientation correlations | `orientational_correlations.py` | `orientational_correlations.csv` | code validated with synthetic vectors; original dipole series required |
