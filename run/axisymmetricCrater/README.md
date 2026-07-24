# Pressure-Driven Crater Explosion andPyroclastic Current Initiation Tutorial

## 1. Case Summary

This tutorial is for an axisymmetric, crater-based explosion using real-world topography. The DEM included is a 2D transect from Ubehebe Crater in Death Valley, California. 

The case models pressurized subsurface mixture of gas and particles, simulating a phreatomacmatic explosion in the moments after the conversion from magmatic thermal energy into kinetic energy, and after the initial fragmentation of the surrounding country rock.
The simulation features:

- **A multi-component gas phase** (H20 and air).

- **One solid phase** (`particles`)

- A **two-stage simulation process**, managed by swapping configuration files.

  - **Stage 1: Hydrostatic Initialization.** A preliminary run establishes a stable atmosphere. Using separate `.init` files is necessary because the inlet conditions would otherwise create an incorrect atmospheric profile.
  - **Stage 2: Main Simulation.** The main run uses different `.run` files to define the inlet conditions.

- The `Allrun` script executes two separate runs (`.init` and `.run`) for the two stages.

This tutorial is an excellent example of how to use OpenFOAM to replicate and explore the physics described in a scientific publication.

______________________________________________________________________


## 2. How to Run the Case

### Automated Execution

First, make the scripts executable:

```bash
chmod +x Allrun Allclean
```

To run the entire workflow from start to finish, execute:

./01_ 
./02_
decomposePar
./03_

The simulation will then be ready for a parallel job submission using a slurm script of your choice. 

------------------------------------------------------------------------
## 4. Cleaning the Case

To remove all generated data and reset the case to its original state, run the `Allclean` script:

```bash
./Allclean
```

______________________________________________________________________

## 5. Description of Key Files

- **`Allrun` / `Allclean`**: Master scripts for running the two simulation scenarios or cleaning the case.
- **`system/controlDict.init` & `.run`**: Control dictionaries for the two different simulation runs. They likely differ in end time, write intervals, or other run-time parameters.
- **`system/fvSolution.init` & `.run`**: Solution dictionaries, which might specify different numerical schemes or solver tolerances for each run.
- **`org.0/`**: A "template" directory for the initial conditions. It contains two complete sets of boundary conditions for each field, distinguished by `.init` and `.run` suffixes. The `Allrun` script renames the appropriate files before each simulation stage. This is a common and effective method for managing multiple case setups within a single directory.
