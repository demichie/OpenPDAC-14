# `topoGrid`

`topoGrid` deforms an existing OpenFOAM mesh using elevations read from an ESRI
ASCII raster. Its main purpose is to convert an initially flat-bottomed mesh
into a terrain-following mesh.

## Typical workflow

```text
blockMesh
    ↓
decomposePar
    ↓
createZones
    ↓
topoGrid
    ↓
checkMesh
```

The exact sequence depends on the case. Both `testVulcano` and
`Etna_refinement` provide complete parallel examples.

## Parallel meshes

The utility supports decomposed meshes and synchronizes positions and
topological information associated with points shared by different processor
domains.

## Input topography

The elevation input is an ESRI ASCII raster, normally a digital elevation model
(DEM). The raster location and mesh-deformation settings must be consistent
with the dictionaries and scripts supplied by the case.

## Compilation

```bash
cd applications/utilities/topoGrid
wmake
```

The resulting executable is installed in `$FOAM_USER_APPBIN`.

