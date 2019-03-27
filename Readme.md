Extended B-spline Adaptive COLlocation with error control by Interpolation.

Purpose
=======

The purpose of eBACOLI is to solve m~u~ dimensional systems of second order parabolic partial differential equations (PDEs) in one space variable that are coupled to an additional m~v~ ordinary differential equations (ODEs) in the form:

u~t~(t,x) = f(t, x, v(t,x), u(t,x), u~x~(t,x), u~xx~(t,x)),

v~t~(t,x) = g(t, x, v(t,x), u(t,x)),

where x~a~ &lt; x &lt; x~b~ and t &gt; t0, with initial conditions at time t = t0 are given by:

u(t0,x) = u~0~(x),

v(t0,x) = v~0~(x)

for a &lt;= x &lt;= b, subject to separated boundary conditions given by:

Ba(t, u(t,a), u~x~(t,a), v(t,a), v~x~(t,a)) = 0,

Bb(t, u(t,b), u~x~(t,b), v(t,b), v~x~(t,b)) = 0,

for t &gt; t0 and x = a, x = b, respectively.

Guide to the above notation:

u~t~(t,x) - denotes the first partial derivative of u(t,x) with respect to the time variable t.

u~x~(t,x) - denotes the first partial derivative of u(t,x) with respect to space variable x.

u~xx~(t,x) - denotes the second partial derivative of u(t,x) with respect to space variable x.

Furthermore, the above functions f, u, Ba and Bb are m~u~ dimensional vector functions, while the functions g and v are m~v~ dimensional vector functions, and the total system size is m~pde~ = m~u~ + m~v~.

How
===

The solution to this problem is adapted from BACOLI.

BACOLI is a method of lines algorithm which uses B-spline collocation to discretize the spatial domain \[a,b\]. BACOLI uses a secondary interpolant of one degree higher or lower to estimate the error of the current solution, and refine either the spatial or temporal grid (or both) accordingly.

The output is a vector of B-spline coefficients which can be used to calculate the approximate solution u(t,x) and its spatial derivatives at (tout,x) where a &lt;= x &lt;= b and t0 &lt; tout.

Building
========

Library and example codes can be built by running \`make\` in this directory. The compiled library (\`libebacoli.so\`) and its \`f95\` module are stored in [./lib](./lib).

Make note of the \`debug\` and \`opt\` build targets.

To solve a specific problem, 2 things need to be supplied:

1.  A driver file (Fortran95) to act as the main program (ie., [driver-toy2-simple.f](./examples/extended/driver-toy2-simple.f95))
2.  A system file (Fortran77) to specify the above systems of equations (ie., [toy2.f](./examples/extended/toy2.f))

Viewing solutions
=================

A python library for handling the B-spline output of eBACOLI is found at [./lib/python/ebacoli.py](./lib/python/ebacoli.py). [./ebacoli.env](./ebacoli.env) sets the appropriate directory in your `PYTHONPATH` for its use.

Example for how to structure the driver and post-processing to generate video frames is found in [./examples/movie](./examples/movie).

Visual results
--------------

Movies of various solutions found in [./media](./media).
