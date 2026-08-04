# OpenPDAC-14

OpenPDAC is an OpenFOAM-based solver for the simulation of multiphase gas–particle flows, with particular emphasis on volcanic and geophysical applications.

The code is derived from the `multiphaseEuler` module distributed with OpenFOAM and extends its formulation to represent multiple dispersed solid phases.

## Main features

OpenPDAC includes:

* a modified kinetic-theory formulation for granular flows involving multiple dispersed solid phases;
* Eulerian modelling of the gas and solid phases;
* a Lagrangian particle-tracking library with one-way coupling to the Eulerian gas–solid mixture;
* adaptive mesh refinement (AMR), including dynamic load balancing for parallel simulations;
* initialization of hydrostatic atmospheric pressure profiles for simulations over large vertical domains;
* boundary-condition support suitable for atmospheric inflow and outflow;
* test cases covering different flow regimes, geometries, and numerical configurations.

## Test cases

The repository provides six test cases:

1. a three-dimensional explosion;
2. a two-dimensional explosion over flat terrain;
3. a two-dimensional dilute flow over a wavy surface;
4. a two-dimensional fluidized bed with two solid phases;
5. a two-dimensional impinging flow with two solid phases;
6. a three-dimensional flow over realistic topography using adaptive mesh refinement (AMR).

These cases are intended both as examples of the available modelling capabilities and as starting points for setting up new simulations.

## OpenFOAM compatibility

This version of OpenPDAC is developed for OpenFOAM 14 and has been tested with the following Ubuntu package:

```text
openfoam14_20260724_amd64.deb
```

Compatibility with other OpenFOAM versions or package releases is not guaranteed.

## Disclaimer

OpenPDAC is an independent development based on OpenFOAM. It is not approved or endorsed by the OpenFOAM Foundation or by ESI Group.
