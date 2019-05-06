<b>E</b>xtended <b>B</b>-spline <b>A</b>daptive <b>COL</b>location with error control by <b>I</b>nterpolation.

Purpose
===

The purpose of eBACOLI is to solve m~u~ dimensional systems of second order parabolic partial differential equations (PDEs) in one space variable that are coupled to an additional m~v~ ordinary differential equations (ODEs) in the form:

u\_t(t,x) = f(t, x, u, u\_x, u\_xx, w, w_x, w_xx),

v\_t(t,x) = g(t, x, u, u_x, u_xx, v, w, w_x, w_xx),

0 = h(t, x, u, u_x, u_xx, v, w, w_x, w_xx)

where x\_a &lt; x &lt; x\_b and t &gt; t0, with initial conditions at time t = t0 are given by:

u(t0,x) = u\_0(x),

v(t0,x) = v\_0(x)

w(t0,x) = w\_0(x)

for a &lt;= x &lt;= b, subject to separated boundary conditions given by:

Ba(t, u(t,a), u\_x(t,a), v(t,a), v\_x(t,a), w(t,a), w\_x(t,a)) = 0,

Bb(t, u(t,b), u\_x(t,b), v(t,b), v\_x(t,b) w(t,b), w\_x(t,b)) = 0,

for t &gt; t0 and x = a, x = b, respectively.

Guide to the above notation:

u\_t - denotes the first partial derivative of u(t,x) with respect to the time variable t.

u\_x - denotes the first partial derivative of u(t,x) with respect to space variable x.

u\_xx - denotes the second partial derivative of u(t,x) with respect to space variable x.

Furthermore, the above functions f and u are m\_u dimensional vector functions, the functions g and v are m\_v dimensional vector functions, the functions h and w are m_\w dimensional vector functions, the boundary conditions Ba and Bb are m\_u+m\_w dimensional vector functions, and the total system size is m\_pde = m\_u + m\_v + m\_w.

How
===

The solution to this problem is adapted from BACOLI.

BACOLI is a method of lines algorithm which uses B-spline collocation to discretize the spatial domain \[a,b\]. BACOLI uses a secondary interpolant of one degree higher or lower to estimate the error of the current solution, and refine either the spatial or temporal grid (or both) accordingly.

The output is a vector of B-spline coefficients which can be used to calculate the approximate solution u(t,x) and its spatial derivatives at (tout,x) where a &lt;= x &lt;= b and t0 &lt; tout.

Building
===

To build:
```
mkdir build && cd build
cmake ..
make
```

`Debug` or `Release` builds can be selected by providing the `-DCMAKE_BUILD_TYPE=<Release|Debug>` flag to `cmake`.

To solve a specific problem, 2 additional things need to be linked to the library:

1.  A driver file (Fortran95) to act as the main program (ie., [driver-toy2-simple.f95](./examples/extended/driver-toy2-simple.f95))
2.  A system file (Fortran77) to specify the above systems of equations (ie., [toy2.f](./examples/extended/toy2.f))

Viewing solutions
===

A python library for handling the B-spline output of eBACOLI is found at [./lib/python/ebacoli.py](./lib/python/ebacoli.py).

If using bash, [./ebacoli.env](./ebacoli.env) sets the appropriate directory in your `PYTHONPATH` for its use. (Run with `source ebacoli.env`)

Example for how to structure the driver and post-processing to generate video frames is found in [./examples/movie](./examples/movie).

Visual results
---

Movies of various solutions found in [./media](./media).
