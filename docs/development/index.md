# Repository and build structure

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
environment.yml                optional Python environment
CITATION.cff                   software citation metadata
```

## Reference sources

The repository retains copies of the original OpenFOAM `multiphaseEuler` and
`parcel` sources. The modified implementations compiled for OpenPDAC are under
`applications/OpenPDAC`. This arrangement facilitates direct comparison and
review of the changes introduced by OpenPDAC.

## Build sequence

`applications/OpenPDAC/Allwmake` compiles:

1. the modified `parcel` library;
2. the modified phase system;
3. momentum-transport models;
4. thermophysical-transport models;
5. `libOpenPDACSolver.so`;
6. finite-volume models;
7. function objects.

`topoGrid` is compiled separately with `wmake`.

## HPC packaging

`packaging/easybuild/OpenPDAC-LUMI.eb.in` provides an EasyBuild template for
LUMI. It installs the solver libraries and utility in a dedicated prefix and
can serve as a starting point for other module-based HPC environments.

