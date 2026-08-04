# Running a tutorial

Each tutorial contains its own dictionaries, scripts, and README. A typical
case can be executed with:

```bash
cd run/<case-name>
./Allrun
```

Depending on the case, `Allrun` may perform mesh generation, domain
decomposition, hydrostatic initialization, the main simulation, and
post-processing.

Many cases also include:

- numbered scripts such as `01_run_meshing.sh` and
  `02_run_fieldInitialization.sh`;
- `Allclean`, which restores the case to its initial state;
- `Alltest`, which checks basic completion criteria.

!!! warning
    Read the README inside the selected tutorial before running it. Some cases
    require a DEM, a specific MPI decomposition, the `OpenPDACconda`
    environment, or additional post-processing software.

## Parallel cases

Parallel execution is controlled by the case's `decomposeParDict`. Several
three-dimensional examples keep their results decomposed, particularly when
AMR changes the mesh topology during the simulation.

## Hydrostatic initialization

Atmospheric cases normally use a preliminary initialization stage before the
physical transient simulation. Separate `.init` and `.run` dictionaries allow
the hydrostatic atmosphere to be established without activating the final
source or inlet conditions prematurely.

