# Hydrostatic initialization

Atmospheric pressure varies substantially over the large vertical domains used
in volcanic-flow simulations. Directly starting from uniform pressure can
produce a strong artificial adjustment that contaminates the initial transient.

OpenPDAC therefore includes a hydrostatic initialization procedure for
stratified atmospheres. The pressure and thermodynamic fields are iteratively
adjusted before the physical simulation begins.

In the supplied atmospheric tutorials, initialization and simulation are
usually separated:

1. `.init` dictionaries establish the atmosphere;
2. the initialized state is written at the starting time;
3. `.run` dictionaries activate the physical source or inlet conditions;
4. the transient simulation starts from the initialized state.

This arrangement also permits atmospheric inflow/outflow boundary conditions
to operate around a consistent background state.

