Parabolic + ODE problems
===

- [toy01-po](./toy01-po/toy01-po.f):
  - Toy model with NPDE=2, NU=1, NV=1

- [toy02-po.f](./toy02-po/toy02-po.f):
  - Toy model with NPDE=4, NU=2, NV=2.
  - Obtained by duplicating solution to toy_par_plus_ode, and interchanging the interaction terms of the system equations in a symmetric way

- [toy03-po](./toy03-po/toy03-po.f):
  - Toy model with NPDE=4, NU=2, NV=2.
  - Obtained by duplicating all of toy2_par_plus_ode, but replacing one of the boundary conditions to use spatial derivative of v
  - Has the same solution as toy2_par_plus_ode
