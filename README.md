Extended B-spline Adaptive COLlocation with error control by Interpolation.

Purpose
=======

The purpose of eBACOLI is to solve $m_u$ dimensional systems of second order parabolic partial differential equations (PDEs) in one space variable that are coupled to an additional $m_v$ ordinary differential equations (ODEs) in the form:

$$
\begin{align}
u_t(t,x) & = f(t, x, v(t,x), u(t,x), u_x(t,x), u_xx(t,x)), \\
v_t(t,x) & = g(t, x, v(t,x), u(t,x))
\end{align}
$$

where $x_a \lt x \lt x_b$ and $t \gt t_0$, with initial conditions at time $t = t_0$ are given by:

$$
\begin{align}
u(t_0,x) & = u_0(x), \\
v(t_0,x) = v_0(x)
\end{align}
$$

for $a \le x \le b$, subject to separated boundary conditions given by:

$$
\begin{align}
B_a(t, u(t,a), u_x(t,a), v(t,a), v_x(t,a)) & = 0, \\
B_b(t, u(t,b), u_x(t,b), v(t,b), v_x(t,b)) & = 0,
\end{align}
$$

for $t \gt t_0$ and $x = a$, $x = b$, respectively.

Guide to the above notation:

$u_t(t,x)$ - denotes the first partial derivative of $u(t,x)$ with respect to the time variable $t$.

$u_x(t,x)$ - denotes the first partial derivative of $u(t,x)$ with respect to space variable $x$.

$u_{xx}(t,x)$ - denotes the second partial derivative of $u(t,x)$ with respect to space variable $x$.

Furthermore, the above functions $f$, $u$, $B_a$ and $B_b$ are $m_u$ dimensional vector functions, while the functions $g$ and $v$ are $m_v$ dimensional vector functions, and the total system size is $m_\text{pde} = m_u + m_v$.

How
===

The solution to this problem is adapted from BACOLI.

BACOLI is a method of lines algorithm which uses B-spline collocation to discretize the spatial domain $[a,b]$. BACOLI uses a secondary interpolant of one degree higher or lower to estimate the error of the current solution, and refine either the spatial or temporal grid (or both) accordingly.

The output is a vector of B-spline coefficients which can be used to calculate the approximate solution u(t,x) and its spatial derivatives at $(t_\text{out},x)$, where $a \le x \le b$ and $t_0 \lt t_\text{out}$.

Building
========

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
=================

A python library for handling the B-spline output of eBACOLI is found at [./lib/python/ebacoli.py](./lib/python/ebacoli.py).

If using bash, [./ebacoli.env](./ebacoli.env) sets the appropriate directory in your `PYTHONPATH` for its use. (Run with `source ebacoli.env`)

Example for how to structure the driver and post-processing to generate video frames is found in [./examples/movie](./examples/movie).

Visual results
--------------

Movies of various solutions found in [./media](./media).
