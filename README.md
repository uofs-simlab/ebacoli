# eBACOLI

Extended B-spline Adaptive COLlocation with error control by Interpolation.

## Purpose

The purpose of eBACOLI is to solve m_u dimensional systems of second order
parabolic partial differential equations (PDEs) in one space variable that are
coupled to an additional m_v ordinary differential equations (ODEs) in the
form:

![System definition](https://latex.codecogs.com/gif.latex?%5Cbegin%7Baligned%7D%20u_t%28t%2Cx%29%20%26%20%3D%20f%28t%2C%20x%2C%20v%28t%2Cx%29%2C%20u%28t%2Cx%29%2C%20u_x%28t%2Cx%29%2C%20u_%7Bxx%7D%28t%2Cx%29%29%2C%20%5C%5C%20v_t%28t%2Cx%29%20%26%20%3D%20g%28t%2C%20x%2C%20v%28t%2Cx%29%2C%20u%28t%2Cx%29%29%2C%20%5Cend%7Baligned%7D)

u_t(t,x) = f(t, x, v(t,x), u(t,x), u_x(t,x), u_xx(t,x)),
v_t(t,x) = g(t, x, v(t,x), u(t,x)),

where x_a < x < x_b and t > t0, with initial conditions at
time t = t0 are given by:

u(t0,x) = u_0(x),
v(t0,x) = v_0(x)

for a <= x <= b, subject to separated boundary conditions
given by:

Ba(t, u(t,a), u_x(t,a), v(t,a), v_x(t,a)) = 0,
Bb(t, u(t,b), u_x(t,b), v(t,b), v_x(t,b)) = 0,

for t > t0 and x = a, x = b, respectively.

Guide to the above notation:

u_t(t,x) - denotes the first partial derivative of u(t,x)
           with respect to the time variable t.

u_x(t,x) - denotes the first partial derivative of u(t,x)
           with respect to space variable x.

u_xx(t,x) - denotes the second partial derivative of u(t,x)
            with respect to space variable x.

Furthermore, the above functions f, u, Ba and Bb are m_u dimensional vector
functions, while the functions g and v are m_v dimensional vector functions,
and the total system size is m_pde = m_u + m_v.

## How

The solution to this problem is adapted from BACOLI.

BACOLI is a method of lines algorithm which uses B-spline collocation
to discretize the spatial domain [a,b]. BACOLI uses a secondary
interpolant of one degree higher or lower to estimate the error of the
current solution, and refine either the spatial or temporal grid (or
both) accordingly.

The output is a vector of B-spline coefficients which can be used to
calculate the approximate solution u(t,x) and its spatial derivatives
at (tout,x) where a <= x <= b and t0 < tout.

## Building

Library and example codes can be built by running `make` in this directory.

Make note of the `debug` and `opt` build targets.

To solve a specific problem, 2 things need to be supplied:
1. A driver file (Fortran95) to act as the main program (ie., [[./examples/extended/driver-toy2-simple.f95][driver-toy2-simple.f]])
2. A system file (Fortran77) to specify the above systems of equations (ie., [[./examples/extended/toy2.f][toy2.f]])

## TODO Viewing solutions

## TODO Visuals

gifs should be placed in-line here.
