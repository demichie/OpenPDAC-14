# Installation

## Load OpenFOAM 14

With the standard Ubuntu package, load the OpenFOAM environment with:

```bash
source /opt/openfoam14/etc/bashrc
```

The path may differ on other systems. On HPC platforms, load the corresponding
OpenFOAM module instead.

## Clone the repository

```bash
git clone https://github.com/demichie/OpenPDAC-14.git
cd OpenPDAC-14
```

## Compile the solver and libraries

```bash
cd applications/OpenPDAC
./Allwmake
```

`Allwmake` compiles the modified `parcel`, `phaseSystem`, momentum-transport,
and thermophysical-transport libraries, followed by the OpenPDAC solver,
finite-volume models, and function objects.

## Compile `topoGrid`

```bash
cd ../utilities/topoGrid
wmake
```

The principal products are installed in the standard OpenFOAM user locations:

```text
$FOAM_USER_LIBBIN/libOpenPDACSolver.so
$FOAM_USER_APPBIN/topoGrid
```

Verify the installation with:

```bash
test -f "$FOAM_USER_LIBBIN/libOpenPDACSolver.so"
test -x "$FOAM_USER_APPBIN/topoGrid"
```

## Optional Python environment

From the repository root:

```bash
conda env create -f environment.yml
conda activate OpenPDACconda
```

The environment is not required to compile the solver, but is used by several
tutorial workflows.

## Run the solver

OpenPDAC is selected through a case's `controlDict` and run as a modular solver:

```bash
foamRun -solver OpenPDAC
```

For tutorial cases, use the supplied scripts instead of invoking the solver
directly.

