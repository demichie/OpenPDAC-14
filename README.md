# OpenPDAC-14

OpenPDAC is an OpenFOAM-based solver for compressible multiphase gas-particle
flows, with particular emphasis on pyroclastic density currents, volcanic
explosions, and related geophysical flows.

The solver is derived from the `multiphaseEuler` module distributed with
OpenFOAM. It extends the original formulation to represent multiple dispersed
solid phases using a modified kinetic theory for granular flows. An optional
Lagrangian particle model is also included for tracking particles transported
by the Eulerian gas-solid mixture.

OpenPDAC-14 is implemented as a modular solver and is run through `foamRun`.

## Main features

OpenPDAC includes:

- Eulerian modelling of compressible gas-solid multiphase flows;
- multiple dispersed solid phases with different physical properties;
- a modified kinetic-theory formulation for polydisperse granular mixtures;
- drag models for spherical and non-spherical particles;
- an adapted Lagrangian `parcel` library with one-way coupling to the
  Eulerian gas-solid mixture;
- distributions of Lagrangian particle size and density, and variable particle
  sphericity;
- hydrostatic initialization of stratified atmospheres over large vertical
  domains;
- boundary-condition support suitable for atmospheric inflow and outflow;
- adaptive mesh refinement and unrefinement during a simulation;
- mapping of Eulerian fields and Lagrangian particles when the mesh topology
  changes;
- dynamic mesh redistribution for load balancing in parallel simulations;
- the `topoGrid` utility for generating terrain-following meshes from
  topographic data;
- tutorial cases covering different geometries, flow regimes, and numerical
  workflows;
- tools for ensemble generation, log analysis, and post-processing.

## Eulerian and Lagrangian phases

The carrier flow and the main particulate components are represented using an
Eulerian multiphase formulation. The kinetic-theory models have been modified
to account for interactions within mixtures containing multiple dispersed
solid phases.

Lagrangian particles can be included in addition to the Eulerian phases. The
adapted `parcel` library uses the properties of the Eulerian gas-solid mixture
as carrier fields. The coupling is currently one-way: the Eulerian mixture
affects the Lagrangian particles, while the Lagrangian particles do not feed
momentum, mass, or energy back into the Eulerian phases.

## Adaptive mesh refinement

OpenPDAC supports adaptive mesh refinement (AMR) using the dynamic-mesh
infrastructure provided by OpenFOAM 14. The implementation handles the mapping
and reconstruction of the OpenPDAC fields when cells are refined or
unrefined. Lagrangian particle positions are preserved when the mesh changes or
is redistributed.

For parallel simulations, AMR can be combined with dynamic mesh
redistribution. This reduces the load imbalance produced when refinement is
concentrated in a limited part of the domain. Refinement criteria,
redistribution frequency, and imbalance thresholds are defined in the case
`dynamicMeshDict`.

The `run/Etna_refinement` tutorial provides an example combining AMR, dynamic
load balancing, terrain-following meshing, and a Lagrangian cloud.

## The `topoGrid` utility

`topoGrid` deforms an existing OpenFOAM mesh using topographic data provided as
an ESRI ASCII raster. It can be used to transform an initially flat-bottomed
mesh into a terrain-following mesh while retaining the original horizontal
mesh structure away from the lower boundary.

The utility supports decomposed meshes and synchronizes points shared by
different processor domains. Several tutorials use `topoGrid` after
`blockMesh`, domain decomposition, and cell-zone creation.

Case-specific options are read from the corresponding OpenFOAM dictionaries.
See `run/testVulcano` and `run/Etna_refinement` for complete examples.

## OpenFOAM reference sources

For reference and comparison, the repository contains copies of the original
OpenFOAM source code on which the OpenPDAC developments are based:

- `applications/multiphaseEuler`: original `multiphaseEuler` solver sources;
- `src/lagrangian/parcel`: original OpenFOAM `parcel` library sources.

The modified sources used by OpenPDAC are located under
`applications/OpenPDAC`. Keeping the upstream sources in the repository makes
it easier to identify and review the changes introduced in OpenPDAC.

## Repository structure

```text
applications/
├── OpenPDAC/                  OpenPDAC solver and modified libraries
├── multiphaseEuler/           reference OpenFOAM solver sources
└── utilities/
    └── topoGrid/              terrain-following mesh utility

src/
└── lagrangian/
    └── parcel/                reference OpenFOAM parcel sources

run/                           tutorials and ensemble workflow
packaging/easybuild/           EasyBuild template for LUMI
environment.yml                optional Python/Conda environment
CITATION.cff                   software citation metadata
```

## Requirements

The solver requires:

- a Linux environment with a working OpenFOAM 14 installation;
- the compiler and build tools required by OpenFOAM, including `wmake`;
- a shell environment in which the OpenFOAM variables have been loaded.

Some tutorials additionally require:

- MPI for parallel simulations;
- Python 3 and the packages listed in `environment.yml`;
- ParaView for visualization;
- external tools such as FFmpeg for selected post-processing workflows.

## OpenFOAM compatibility

This branch is developed for OpenFOAM 14 and has been tested with the following
Ubuntu package:

```text
openfoam14_20260724_amd64.deb
```

Compatibility with other OpenFOAM versions or package releases is not
guaranteed. Older OpenPDAC branches are available for earlier OpenFOAM
versions.

## Building OpenPDAC

First, load the OpenFOAM 14 environment. With the standard Ubuntu package this
can normally be done with:

```bash
source /opt/openfoam14/etc/bashrc
```

The path may differ for other installations or module-based HPC environments.

Clone the repository:

```bash
git clone https://github.com/demichie/OpenPDAC-14.git
cd OpenPDAC-14
```

Compile the OpenPDAC solver and its modified libraries:

```bash
cd applications/OpenPDAC
./Allwmake
```

The script compiles the modified `parcel`, `phaseSystem`, momentum-transport,
and thermophysical-transport libraries, followed by the OpenPDAC solver,
finite-volume models, and function objects.

Compile `topoGrid` separately:

```bash
cd ../utilities/topoGrid
wmake
```

The principal build products are installed in the standard OpenFOAM user
locations:

```text
$FOAM_USER_LIBBIN/libOpenPDACSolver.so
$FOAM_USER_APPBIN/topoGrid
```

To check that they are available:

```bash
test -f "$FOAM_USER_LIBBIN/libOpenPDACSolver.so"
test -x "$FOAM_USER_APPBIN/topoGrid"
```

OpenPDAC is a modular solver and is normally selected in a case through
`controlDict`, then executed with:

```bash
foamRun -solver OpenPDAC
```

The supplied `Allrun` scripts determine the application from each case's
`controlDict` and invoke it with the appropriate serial or parallel options.

## Optional Python environment

Several meshing, ensemble-generation, analysis, and post-processing scripts use
Python packages listed in `environment.yml`. A Conda environment can be created
with:

```bash
conda env create -f environment.yml
conda activate OpenPDACconda
```

This environment is not required to compile the solver, but it is required by
several of the workflows under `run`.

## Tutorials and workflows

The repository currently includes seven simulation tutorials and a separate
framework for ensemble generation and parametric studies.

| Directory | Description |
| --- | --- |
| `synthTopo2D` | Two-dimensional volcanic explosion from a subsurface conduit through a synthetic crater |
| `flow_wave` | Two-dimensional pyroclastic surge over synthetic topography |
| `fluidisedBed` | Polydisperse fluidized bed with two Eulerian solid phases |
| `Valentine2020_twoP` | Two-dimensional column collapse with coarse and fine particle phases, based on Valentine (2020) |
| `testVulcano` | Three-dimensional phreatic explosion over the real topography of Vulcano Island |
| `axisymmetricCrater` | Axisymmetric pressure-driven explosion using an Ubehebe Crater topographic transect |
| `Etna_refinement` | Three-dimensional pyroclastic density current over Mt Etna with AMR and dynamic load balancing |
| `ensembleWorkflow` | Latin-hypercube sampling and automated generation of simulation ensembles |

Each tutorial contains its own README and case-specific scripts. When an
`Allrun` script is provided, the complete workflow can normally be launched
with:

```bash
cd run/<case-name>
./Allrun
```

Many cases also provide:

- numbered scripts for running meshing, initialization, simulation, and
  post-processing separately;
- `Allclean` for restoring the initial case state;
- `Alltest` for checking basic completion criteria.

Always consult the README inside the selected case before running it. Some
three-dimensional examples require a DEM, a particular processor decomposition,
or additional Python and post-processing tools.

## Automated testing

GitHub Actions is used to compile OpenPDAC and `topoGrid` against OpenFOAM 14.
The current continuous-integration workflow runs a shortened `synthTopo2D`
regression case and invokes its validation script when available.

The repository also contains workflows for static analysis, C++ linting, and
coverage generation. The coverage workflow uses the `Valentine2020_twoP` case.

The test runner can be invoked locally with:

```bash
cd run
./run-all-tests.sh
```

At present, the list of cases selected by this script is intentionally limited;
it should not be interpreted as a complete validation of every tutorial.

## HPC installation

An EasyBuild template for installing OpenPDAC with OpenFOAM 14 on LUMI is
provided under:

```text
packaging/easybuild/OpenPDAC-LUMI.eb.in
```

The template installs the solver libraries and `topoGrid` in a dedicated
software prefix. It may be adapted to other module-based HPC systems, but the
toolchain and OpenFOAM dependency must be checked for the target platform.

## Citation

If you use OpenPDAC in scientific work, please cite the software using the
metadata provided in [`CITATION.cff`](CITATION.cff). GitHub also exposes these
metadata through the **Cite this repository** function.

When appropriate, please also cite the scientific publications describing the
physical and numerical formulations used by the selected OpenPDAC model and
tutorial.

## License

OpenPDAC is distributed under the GNU General Public License, version 3 or any
later version, as stated in the source files. It is provided without warranty;
see the source-code license notices for details.

The repository also contains source files derived from OpenFOAM and retains
their original copyright and license notices.

## Disclaimer

OpenPDAC is an independent development based on OpenFOAM. It is not approved,
supported, or endorsed by the OpenFOAM Foundation.

## Author and contact

OpenPDAC is developed by:

**Mattia de' Michieli Vitturi**  
Istituto Nazionale di Geofisica e Vulcanologia (INGV), Sezione di Pisa  

**Tomaso Esposti Ongaro**  
Istituto Nazionale di Geofisica e Vulcanologia (INGV), Sezione di Pisa  

**Federica Pardini**
Istituto Nazionale di Geofisica e Vulcanologia (INGV), Sezione di Pisa  

**Federico Brogi**
CINECA, Direzione HPC Supercalcolo, Sede di Bologna 

**Brandon Keim**
University at Buffalo, Department of Geology, Buffalo

**OpenPDAC Team**  
Istituto Nazionale di Geofisica e Vulcanologia (INGV), Sezione di Pisa  
[openpdac@ingv.it](mailto:openpdac@ingv.it)



Bug reports and reproducible examples can be submitted through the repository's
[GitHub issue tracker](https://github.com/demichie/OpenPDAC-14/issues).
