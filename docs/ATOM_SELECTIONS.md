# Verified atom selections

The production system contains **one QMAP and 21 DMSO molecules (241 atoms)** in a cubic cell of side 16.0 Å. The CPMD XYZ ordering is grouped by element, not by molecule.

## Element blocks

| Element | 1-based indices | Count |
|---|---:|---:|
| H | 1–138 | 138 |
| C | 139–196 | 58 |
| N | 197–198 | 2 |
| O | 199–220 | 22 |
| S | 221–241 | 21 |

## QMAP selection

QMAP comprises H 1–12, C 139–154, N 197–198 and O 199. Key labels used in the manuscript are:

| Label | XYZ index | Role |
|---|---:|---|
| O1 | 199 | phenolic donor oxygen |
| H1 | 12 | phenolic hydrogen |
| N1 | 198 | quinoline nitrogen |
| N2 | 197 | imine nitrogen |
| C1 | 150 | carbon used in Figure 2 |
| C2 | 151 | carbon used in Figure 2 |
| C3 | 152 | carbon used in Figure 2 |
| C4 | 153 | carbon used in Figure 2 |

The full label map and all 21 DMSO molecule assignments are machine-readable in `config/atom_indices.json`.

## Analysis selections

- Figure 2 distance: C1–N1 = 150–198. Both raw and minimum-image values are exported; the manuscript plot used the raw coordinate difference and is PBC-sensitive.
- Figure 2 angle: C2–C1–N1 = 151–150–198.
- Figure 2 four-point torsion: C2–C1–C3–C4 = 151–150–152–153.
- Figure 6: O1(QMAP) = 199 against O(DMSO) = 200–220.
- Figure 7: N1 = 198 and N2 = 197 against S(DMSO) = 221–241. Separate N1, N2 and nearest-N curves are exported because the original interactive plotting record did not preserve the exact nitrogen selection.
- Figure 8 donor: O199–H12; acceptors: O200–O220; criteria H···A ≤ 2.5 Å and D–H···A ≥ 150°.

Validation uses `structures/GEOMETRY_reference.xyz` and `tests/test_reference_geometry.py`.
