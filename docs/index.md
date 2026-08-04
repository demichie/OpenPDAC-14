---
hide:
  - navigation
  - toc
---

# OpenPDAC-14

<div class="hero-copy" markdown>

**A multiphase solver for pyroclastic density currents, volcanic explosions,
and related gas-particle flows.**

OpenPDAC extends the OpenFOAM `multiphaseEuler` formulation to represent
multiple dispersed solid phases, adapted granular-flow models, and optional
Lagrangian particles transported by the Eulerian gas-solid mixture.

[Get started](getting-started/index.md){ .md-button .md-button--primary }
[View tutorials](tutorials/index.md){ .md-button }
[Source code](https://github.com/demichie/OpenPDAC-14){ .md-button }

</div>

<div class="feature-grid" markdown>

<div class="feature-card" markdown>

### Multiphase formulation

Compressible Eulerian gas-solid flows with multiple dispersed solid phases and
a modified kinetic theory for polydisperse granular mixtures.

</div>

<div class="feature-card" markdown>

### Lagrangian particles

An adapted `parcel` library with particle-size and density distributions,
variable sphericity, and one-way coupling to the Eulerian mixture.

</div>

<div class="feature-card" markdown>

### Adaptive meshes

Adaptive refinement and unrefinement, particle remapping, and dynamic load
balancing for parallel simulations.

</div>

<div class="feature-card" markdown>

### Real topography

The `topoGrid` utility transforms initially flat-bottomed meshes into
terrain-following meshes using ESRI ASCII elevation rasters.

</div>

</div>

## Applications

OpenPDAC is primarily developed for the simulation of volcanic gas-particle
flows. The supplied examples cover volcanic explosions, collapsing eruptive
columns, dilute and concentrated pyroclastic currents, ballistic particles,
and polydisperse fluidized beds.

The repository also includes workflows for hydrostatic atmospheric
initialization, DEM-based mesh generation, adaptive mesh refinement, ensemble
generation, log analysis, and scientific post-processing.

## OpenFOAM compatibility

This branch is developed for **OpenFOAM 14** and has been tested with:

```text
openfoam14_20260724_amd64.deb
```

Compatibility with other OpenFOAM versions or package releases is not
guaranteed.

!!! info "Modular solver"
    OpenPDAC-14 is loaded as a solver module and is normally executed with
    `foamRun -solver OpenPDAC`.

