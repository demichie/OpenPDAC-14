# Testing and continuous integration

The GitHub Actions workflows compile OpenPDAC and `topoGrid` against OpenFOAM
14. The current regression workflow runs a shortened `synthTopo2D` case and
invokes its `Alltest` validation when available.

The repository also contains workflows for:

- CodeQL static analysis;
- C++ linting;
- coverage instrumentation and report generation.

The coverage workflow uses the `Valentine2020_twoP` tutorial.

## Local test runner

```bash
cd run
./run-all-tests.sh
```

The script can shorten selected simulations when `CI=true`, run `Allrun`, call
`Alltest` when present, clean generated data, and restore temporarily modified
dictionaries.

!!! note
    The current case list in `run-all-tests.sh` is limited and does not
    represent validation of every tutorial in the repository.

