# Etna pyroclastic density current with AMR

This OpenPDAC case simulates the generation and propagation of a pyroclastic density current over the real topography of Mt Etna.

A hot mass of solid particles is introduced for 4 s inside a source region located above the terrain. After injection, the particle-rich mixture collapses under gravity, impacts the ground, and develops into a laterally spreading pyroclastic current. The collapse velocity and direction are not prescribed: they result from the multiphase dynamics resolved by OpenPDAC.

The case includes:

- a terrain-following mesh generated from a digital elevation model (DEM);
- a vertically graded background mesh, with finer resolution close to the ground;
- hydrostatic initialization of a stratified atmosphere;
- adaptive mesh refinement (AMR) around the particle-rich current;
- dynamic redistribution of the mesh during the parallel simulation;
- an automated three-stage workflow for meshing, initialization, and simulation.

## Numerical configuration

### Computational mesh

The precursor mesh is generated with `blockMesh`. It consists of horizontal layers with increased vertical resolution in the lowest part of the domain, where the current propagates and interacts with the terrain.

The mesh is then decomposed for parallel execution and processed in two steps:

1. `createZones` creates the zones required by the case, including the particle source region;
2. `topoGrid` deforms the mesh so that its lower boundary conforms to the topography read from the DEM.

Mesh quality is checked both before and after the terrain deformation.

### Atmospheric initialization

Before the physical simulation starts, OpenPDAC is run with the initialization dictionaries. In this stage, `hydrostaticInitialisation` is enabled and up to 100 hydrostatic correctors are used to establish a stratified atmosphere over the terrain-following mesh.

The initialization is performed at simulation time `0`, after which the main run starts from the latest available state.

### Particle source

The source is defined in `constant/fvModels` as a `massSource` acting on the `particles` phase inside the `particle_source` cell zone.

The mass-flow history is a square pulse:

- start time: `0 s`;
- duration: `4 s`;
- configured mass-flow-rate value: `1e8` in OpenFOAM SI units.

The thermal and material properties of the injected phase are specified by the run-stage phase and cloud dictionaries, including `constant/cloudProperties.run` and the corresponding field files.

### Adaptive mesh refinement

AMR is configured in `constant/dynamicMeshDict`. Refinement is driven by `magGradAlpha`, so additional resolution is concentrated around strong gradients in particle concentration, typically near the current margins and upper interface.

The current settings are:

| Parameter | Value |
|---|---:|
| Refinement interval | 40 time steps |
| Lower refinement threshold | `1e-4` |
| Upper refinement threshold | `1.0` |
| Buffer layers | 3 |
| Maximum refinement levels | 2 |
| Maximum number of cells | 650,000 |
| Write refinement level | enabled |
| Check mesh after refinement | enabled |

### Dynamic load balancing

The mesh is redistributed during the parallel run to reduce load imbalance caused by local refinement. Redistribution is checked every 50 time steps and is triggered when the estimated CPU imbalance exceeds 10%.

The `multiConstraint` option is enabled so that the Eulerian mesh load and any additional CPU-load fields can be handled as separate constraints.

### Main simulation

The main simulation uses `system/controlDict.run` and `system/fvSolution.run`.

The default run settings are:

| Parameter | Value |
|---|---:|
| Solver | `foamRun` with `OpenPDAC` |
| Start | latest initialized time |
| End time | `30 s` |
| Initial time step | `1e-4 s` |
| Time-step control | adaptive |
| Maximum Courant number | `0.25` |
| Output interval | `0.5 s` |
| Output format | binary |

## Requirements

The case requires:

- OpenFOAM 13 with a compatible OpenPDAC installation;
- the OpenFOAM utilities `blockMesh`, `decomposePar`, `checkMesh`, `createZones`, and `topoGrid`;
- a valid `system/decomposeParDict` for the desired number of MPI processes;
- the DEM and the terrain-processing dictionaries referenced by the case;
- a Conda environment named `OpenPDACconda`, as activated by the meshing script;
- Python 3 with `pandas` and `matplotlib` for the optional log-analysis utility.

The solver and OpenFOAM environment must be loaded before running the scripts. For example, `WM_PROJECT_DIR` must be defined and the OpenPDAC libraries must be available through the active OpenFOAM environment.

## Quick start

Run the complete workflow with:

```bash
./Allrun
```

`Allrun` stops immediately if one of the stages returns a non-zero exit status.

The three stages can also be executed separately:

```bash
./01_run_meshing.sh
./02_run_fieldInitialization.sh
./03_run_simulation.sh
```

This is useful when inspecting the mesh or initialized atmosphere before starting the main simulation.

## Workflow

### 1. `01_run_meshing.sh` — mesh generation

This script:

1. loads the standard OpenFOAM `RunFunctions`;
2. activates the `OpenPDACconda` Conda environment;
3. cleans previous results by running `Allclean`;
4. copies `system/controlDict.init` to `system/controlDict`;
5. generates the layered precursor mesh with `blockMesh`;
6. recreates the `0` directory from `org.0`;
7. decomposes the case with `decomposePar`;
8. checks the undeformed mesh in parallel and saves its log as `log.checkMesh0`;
9. creates the required zones with `createZones`;
10. deforms the mesh onto the DEM topography with `topoGrid`;
11. checks the final terrain-following mesh in parallel.

The mesh remains decomposed in the `processor*` directories.

### 2. `02_run_fieldInitialization.sh` — atmospheric initialization

This script copies the initialization dictionaries:

```text
system/controlDict.init       -> system/controlDict
system/fvSolution.init        -> system/fvSolution
constant/cloudProperties.init -> constant/cloudProperties
```

It then launches OpenPDAC in parallel. The initialization run establishes the hydrostatic stratification and writes its log to `log.foamRun0`.

### 3. `03_run_simulation.sh` — transient simulation

This script copies the run-stage dictionaries:

```text
system/controlDict.run       -> system/controlDict
system/fvSolution.run        -> system/fvSolution
constant/cloudProperties.run -> constant/cloudProperties
```

It then starts the transient OpenPDAC simulation in parallel. AMR and dynamic load balancing are active during this stage. The main solver output is written to `log.foamRun`.

## Validation and cleanup

To verify that the solver completed successfully, run:

```bash
./Alltest
```

The test checks that `log.foamRun` exists and contains the final `End` marker.

To remove generated time directories, logs, the decomposed case, and dictionaries copied by the workflow scripts, run:

```bash
./Allclean
```

The original templates such as `controlDict.init`, `controlDict.run`, `fvSolution.init`, and `fvSolution.run` are retained.

## Log diagnostics

The optional `analyzeLog.py` utility reads a large OpenFOAM log in streaming mode and generates numerical diagnostics without loading the complete file into memory.

Run it after the simulation with:

```bash
python3 analyzeLog.py log.foamRun
```

By default, outputs are written to `log_analysis/` and include:

- `log_diagnostics.csv`;
- execution-time diagnostics;
- time-step and Courant-number plots;
- pressure and PIMPLE convergence diagnostics;
- phase velocity diagnostics;
- temperature and energy-residual plots;
- volume-fraction, packing, and granular-temperature diagnostics.

Figures are saved in both PNG and PDF format.

A different output directory can be selected with:

```bash
python3 analyzeLog.py log.foamRun --output-dir my_log_analysis
```

## Particle VTK series

If the `cloudInfo1` function object is enabled and writes particle CSV files under
`postProcessing/cloudInfo1/<time>/output.csv`, the optional script
`VTKballistics_series.py` converts them into legacy VTK files and creates a ParaView time series:

```bash
python3 VTKballistics_series.py
```

The generated series is:

```text
postProcessing/cloudInfo1/particles.vtk.series
```

## Visualization

Because AMR changes the mesh topology during the run, the most direct approach is to keep the results decomposed and open the case in ParaView using its parallel OpenFOAM reader.

For example, create an empty case marker if one is not already present:

```bash
touch etnaPDC.foam
```

Then open `etnaPDC.foam` in ParaView and enable the decomposed/parallel case option as appropriate for the installed reader.

If reconstruction is required, dynamic-mesh reconstruction must begin from the earliest time at which the corresponding mesh exists. Starting reconstruction from a later time whose mesh was inherited from an earlier time can cause OpenFOAM to stop with a mesh-reconstruction error.

## Main files

```text
Allrun                         complete workflow
Allclean                       restore the case to its initial state
Alltest                        basic completion test
01_run_meshing.sh              layered mesh generation and DEM deformation
02_run_fieldInitialization.sh  hydrostatic atmospheric initialization
03_run_simulation.sh           main transient OpenPDAC run
constant/dynamicMeshDict       AMR and dynamic load-balancing settings
constant/fvModels              finite-volume particle mass source
system/controlDict.init        initialization-stage time controls
system/controlDict.run         main-run time and output controls
system/fvSolution.init         hydrostatic initialization settings
system/fvSolution.run          main PIMPLE and linear-solver settings
analyzeLog.py                  optional solver-log diagnostics
VTKballistics_series.py        optional particle CSV-to-VTK conversion
```

## Notes

- The MPI decomposition is controlled by `system/decomposeParDict`.
- The DEM path and terrain-deformation settings must be consistent with the dictionaries used by `topoGrid`.
- The source cell zone must be named `particle_source`, as expected by `constant/fvModels`.
- AMR can substantially increase the computational cost. Adjust `maxRefinement`, `maxCells`, the refinement thresholds, and the decomposition settings together.
- The refinement and redistribution intervals are expressed in solver time steps, not physical seconds.
- `03_run_simulation.sh` leaves the results in decomposed form; it does not call a reconstruction utility.
