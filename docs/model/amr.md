# Adaptive mesh refinement and load balancing

OpenPDAC supports adaptive mesh refinement and unrefinement through the
dynamic-mesh infrastructure of OpenFOAM 14.

When the topology changes, the solver maps the Eulerian fields, reconstructs
temporary internal fields where necessary, and preserves the global positions
of Lagrangian particles. A conservative time-step reduction can also be
applied around mesh updates.

## Dynamic load balancing

Local refinement can produce a severe imbalance between MPI ranks. OpenPDAC
cases can combine the mesh `refiner` with the OpenFOAM `distributor`, allowing
the mesh to be redistributed when the estimated imbalance exceeds a specified
threshold.

The main controls are defined in `dynamicMeshDict`:

- refinement field and thresholds;
- refinement interval;
- buffer layers and maximum refinement level;
- maximum cell count;
- redistribution interval;
- maximum permitted load imbalance.

## Etna example

The `Etna_refinement` tutorial refines the mesh around gradients in particle
concentration and periodically redistributes the decomposed mesh. It combines:

- a vertically graded background mesh;
- DEM-based terrain deformation;
- hydrostatic atmospheric initialization;
- a finite-duration particle source;
- AMR around the developing current;
- dynamic load balancing during the parallel run.

