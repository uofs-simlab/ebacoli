# eBACOLI tests

Problem definitions that have exact solutions.

Each test requires 2 files:
- Problem definitions in `<problem.f>`,
- Driver routine in `test_<problem>.f95`.

Different systems are separated into their own subdirectories.

Useful subroutines for comparing errors found in `utilities.f95`.

- [parabolic problems](./par/)
- [parabolic + ODE problems](./par-ode/)
- [parabolic + elliptic problems](./par-ellip/)
- [parabolic + elliptic problems](./par-ellip/)
- [parabolic + ODE + elliptic problems](./par-ellip/)
