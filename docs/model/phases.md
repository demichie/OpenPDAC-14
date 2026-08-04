# Eulerian and Lagrangian phases

## Eulerian formulation

The carrier gas and the main particulate components are represented with an
Eulerian multiphase formulation derived from OpenFOAM's `multiphaseEuler`
module. OpenPDAC modifies the kinetic-theory treatment to model mixtures
containing multiple dispersed solid phases.

Different solid phases can represent particle classes with distinct material
properties. The repository also provides drag models for spherical and
non-spherical particles.

## Lagrangian formulation

Lagrangian particles can be included in addition to the Eulerian phases. The
adapted `parcel` library uses mixture properties derived from the Eulerian
gas-solid phases as carrier fields.

The Lagrangian implementation supports particle-size and density distributions,
variable sphericity, injection models, and cloud function objects used for
analysis and post-processing.

## Coupling limitation

The coupling is currently one-way: the Eulerian mixture affects particle
trajectories, but Lagrangian particles do not feed mass, momentum, or energy
back into the Eulerian phases.

This limitation should be considered when the Lagrangian particle loading is
large enough to alter the carrier mixture significantly.

