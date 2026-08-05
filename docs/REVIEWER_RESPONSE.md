# Suggested response to Reviewer 4

**Reviewer comment 4:** To improve reproducibility, the Python scripts and VMD scripts used in the analysis should be made available in a public repository, e.g., GitHub, or an institutional repository.

**Response:** We agree with the reviewer. The analysis routines that were originally developed interactively have been consolidated and reimplemented as standalone, documented Python scripts. The public repository contains routines for the QMAP internal-coordinate analysis, O(QMAP)–O(DMSO) and N(QMAP)–S(DMSO) pair-distance distributions, directional hydrogen-bond counting, molecular and cluster dipoles from atomic charges, the qualitative dipole-autocorrelation IR-like spectrum, and first- and second-rank orientational correlation functions. It also contains the verified atom-index mapping, a reference geometry, software dependencies, execution instructions, and explicit documentation of the remaining raw-data requirements. VMD was used only for visual inspection and molecular representations; no VMD/Tcl scripts were used in the quantitative analyses. The repository URL has been added to the Methods and Data availability sections.

Repository: https://github.com/ProfRenatoMedeirosFisico/qmap-cpmd-analysis

Permanent archived raw-data release/DOI: PENDING
