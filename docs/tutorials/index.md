# Tutorials and workflows

The repository contains seven simulation tutorials and a separate framework
for ensemble generation.

| Directory | Description | Main capabilities |
| --- | --- | --- |
| `synthTopo2D` | Two-dimensional volcanic explosion through a synthetic crater | Local refinement, conduit initialization, hydrostatic atmosphere |
| `flow_wave` | Two-dimensional pyroclastic surge over synthetic topography | `snappyHexMesh`, `extrudeMesh`, two-stage initialization |
| `fluidisedBed` | Polydisperse fluidized bed | Two Eulerian solid phases, kinetic theory |
| `Valentine2020_twoP` | Collapsing column with coarse and fine particles | Two solid phases, dilute and concentrated currents |
| `testVulcano` | Three-dimensional phreatic explosion at Vulcano | DEM, `topoGrid`, Lagrangian particles, post-processing |
| `axisymmetricCrater` | Pressure-driven explosion at Ubehebe Crater | Axisymmetric geometry, multicomponent gas, real topography |
| `Etna_refinement` | Pyroclastic density current at Mt Etna | AMR, load balancing, terrain-following mesh |
| `ensembleWorkflow` | Automated parametric studies | Latin-hypercube sampling, templated ensemble generation |

## Choosing a starting point

- Use `fluidisedBed` to examine the polydisperse kinetic-theory formulation in
  a compact geometry.
- Use `synthTopo2D` for a complete but relatively contained atmospheric
  explosion workflow.
- Use `Valentine2020_twoP` to study the transition between dilute and
  concentrated currents.
- Use `testVulcano` to explore DEM processing, ballistic particles, and the
  full post-processing workflow.
- Use `Etna_refinement` as the reference for AMR and dynamic load balancing.

Consult the README inside each case for its current requirements and commands.

