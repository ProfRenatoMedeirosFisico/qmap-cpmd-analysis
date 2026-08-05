# QMAP–DMSO CPMD analysis

Reproducibility repository for the manuscript:

**Car-Parrinello Molecular Dynamics of QMAP-DMSO Microsolvation: First-Shell Structure and Solvent-Induced Polarization**

Authors: Renato Medeiros, Francisco A. P. Osório, Clodoaldo Valverde, and Ademir J. Camargo.

## Scope

This repository is intended to provide the analysis scripts and processed files required to reproduce the structural, hydrogen-bond, dipole-moment, qualitative IR-like, and orientational-correlation results reported in the manuscript.

> **Repository status:** preparation package. Before making the repository public and citing it in the manuscript, replace all `PENDING` items in `REPOSITORY_COMPLETION_CHECKLIST.md` with the actual scripts and processed outputs used in the study.

## Repository structure

```text
scripts/
  python/   Python analysis scripts
  vmd/      VMD/Tcl scripts
data/
  processed/
    figures/
    tables/
  example/  Small non-production examples, when permitted
docs/       Method notes, atom selections, commands, and provenance
```

## Available script

### `scripts/python/analisar_qmap_cpmd_correlacoes.py`

Reads a multi-frame XYZ or PDB trajectory and can calculate:

- QMAP molecular dipole from supplied atomic charges;
- minimum O(QMAP)–O(DMSO) distance per frame;
- a user-specified QMAP backbone dihedral;
- time-series and correlation plots;
- `qmap_correlacoes.csv` and an execution report.

The automatic QMAP/DMSO assignment must be validated visually. Publication-quality runs must use the verified atom indices documented in `docs/ATOM_SELECTIONS.md`.

## Installation

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
```

Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

Linux/macOS:

```bash
source .venv/bin/activate
pip install -r requirements.txt
```

## Basic execution

```bash
python scripts/python/analisar_qmap_cpmd_correlacoes.py \
  --traj path/to/trajectory.xyz \
  --charges path/to/charges.txt \
  --dihedral C2 C1 C3 C4 \
  --outdir results/correlations
```

The four dihedral indices are 1-based by default. Add `--zero-based` only when the supplied indices start at zero.

## Reproducibility requirements before publication

The public release should contain:

1. every Python script actually used to generate manuscript figures;
2. every VMD/Tcl script and the exact atom selections;
3. processed numerical tables underlying the figures;
4. a command log or workflow showing how each figure was generated;
5. software versions and dependencies;
6. verified trajectory composition, atom ordering, time step, sampling interval, and units;
7. a permanent archived release, preferably with a DOI.

Large trajectories should normally be archived in Zenodo or another data repository rather than committed directly to GitHub. Git LFS is an alternative when a GitHub-hosted trajectory is necessary.

## Citation

Use `CITATION.cff`. Add the article DOI and repository archive DOI when available.

## License

The preparation package proposes the MIT License for code. Confirm that all authors agree before the public release. Data and manuscript figures may require a separate data license and publisher-policy check.
