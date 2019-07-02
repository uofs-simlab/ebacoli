eBACOLI tests
===

Problem definitions that have exact solutions.

Each test requires 2 files:
- Problem definition in `<problem.f>`,
- Driver routine in `test_<problem>.f95`.

Useful subroutines for comparing errors found in `utilities.f95`.

## Problems:

### Standard BACOLI

- [burg1](./burg1.f)
  - Viscous Burgers equation with single tanh wave front.
  - `npde` replicas of the basic system (`npde=1` in [test_burg1.f95](test_burg1.f95))

- [burg2](./burg2.f)
  - Viscous Burgers equation with two interacting wave fronts.
  - `npde` replicas of the basic system (`npde=4` in [test_burg2.f95](test_burg2.f95))

- [sincmads](./sincmads.f)
  - A single diffusion equation with sinusoidal spatial forcing

### parabolic + ODE

- [toy_par_plus_ode](./toy_par_plus_ode.f):
  - Toy model with NPDE=2, NU=1, NV=1

- [toy2_par_plus_ode.f](./toy2_par_plus_ode.f):
  - Toy model with NPDE=4, NU=2, NV=2.
  - Obtained by duplicating solution to toy_par_plus_ode, and interchanging the interaction terms of the system equations in a symmetric way

- [toy3_par_plus_ode](./toy3_par_plus_ode.f):
  - Toy model with NPDE=4, NU=2, NV=2.
  - Obtained by duplicating all of toy2_par_plus_ode, but replacing one of the boundary conditions to use spatial derivative of v
  - Has the same solution as toy2_par_plus_ode

### parabolic + elliptic


### parabolic + ODE + elliptic
