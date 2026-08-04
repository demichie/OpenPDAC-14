# Getting started

OpenPDAC-14 consists of the solver library, modified OpenFOAM libraries,
additional finite-volume models and function objects, and the `topoGrid`
utility. Tutorial cases are provided under `run/`.

The usual workflow is:

1. install and load OpenFOAM 14;
2. compile OpenPDAC and `topoGrid`;
3. optionally create the Python environment used by preprocessing and
   post-processing scripts;
4. select a tutorial and read its case-specific README;
5. execute the complete `Allrun` workflow or its numbered stages.

## Requirements

The solver requires a Linux installation of OpenFOAM 14 and its standard build
tools, including `wmake`. Parallel cases also require MPI.

Several tutorial workflows use Python for DEM processing, geometry generation,
sampling, analysis, and visualization. Their dependencies are defined in the
repository's `environment.yml` file.

[Install OpenPDAC](installation.md){ .md-button .md-button--primary }
[Run a tutorial](running.md){ .md-button }

