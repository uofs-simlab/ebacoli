# Extended B-spline Adaptive COLlocation with error control by Interpolation.

## Purpose

The purpose of eBACOLI is to solve m_u dimensional systems of second-order parabolic partial differential equations (PDEs) in one space variable that are coupled to an additional m_v ordinary differential equations (ODEs) in the form:

u_t(t,x) = f(t, x, v(t,x), u(t,x), u_x(t,x), u_xx(t,x)),

v_t(t,x) = g(t, x, v(t,x), u(t,x)),

where x\_a &lt; x &lt; x\_b and t &gt; t0, with initial conditions at time t = t0 are given by:

u(t0,x) = u\_0(x),

v(t0,x) = v\_0(x)

for a &lt;= x &lt;= b, subject to separated boundary conditions given by:

Ba(t, u(t,a), u~x~(t,a), v(t,a), v\_x(t,a)) = 0,

Bb(t, u(t,b), u~x~(t,b), v(t,b), v\_x(t,b)) = 0,

for t &gt; t0 and x = a, x = b, respectively.

Guide to the above notation:

u\_t(t,x) - denotes the first partial derivative of u(t,x) with respect to the time variable t.

u\_x(t,x) - denotes the first partial derivative of u(t,x) with respect to space variable x.

u\_xx(t,x) - denotes the second partial derivative of u(t,x) with respect to space variable x.

Furthermore, the functions f, u, Ba, and Bb are m\_u-dimensional vector functions, the functions g and v are m\_v-dimensional vector functions, and the total system size is m\_pde = m\_u + m\_v.

## How

The solution to this problem is adapted from BACOLI.

BACOLI is a method-of-lines algorithm that uses B-spline collocation to discretize the spatial domain \[a,b\]. BACOLI uses a secondary interpolant of one degree higher or lower to estimate the error of the current solution and refines the spatial or temporal grid accordingly.

The output is a vector of B-spline coefficients that can be used to calculate the approximate solution u(t,x) and its spatial derivatives at (tout,x), where a &lt;= x &lt;= b and t0 &lt; tout.

## Building

To build:
```
mkdir build && cd build
cmake ..
make
```

`Debug` or `Release` builds can be selected by providing the `-DCMAKE_BUILD_TYPE=<Release|Debug>` flag to `cmake`.

To solve a specific problem, two additional things need to be linked to the library:

1.  A driver file (Fortran95) to act as the main program (e.g., [driver-toy2-simple.f95](./examples/extended/driver-toy2-simple.f95))
2.  A system file (Fortran77) to specify the above systems of equations (e.g., [toy2.f](./examples/extended/toy2.f))

## Linking with your own problems

### Without cmake

With the build instructions above (note there is currently no `install` target), the eBACOLI library should be found at `/path/to/ebacoli/build/src/libebacoli.{a,so}` and the required module file in `/path/to/ebacoli/build/modules/ebacoli95_mod.mod`. Therefore linking requires the command line options `-I/path/to/ebacoli/build/modules -I/path/to/ebacoli/build/src -lebacoli`.

### With cmake

Exploiting our setup for cmake, new examples can be added by following the conventions applied in the  `./examples/*/CMakeLists.txt` files.

If you wish to add a subdirectory for your code, you need to:
1. Inform the `CMakeLists.txt` file of the parent directory about your new directory using the `add_subdirectory(<dir name>)` function.
2. Create a `CMakeLists.txt` file in the new subdirectory.
3. Add the compile and link rules for your code. For simplicity, you can use our cmake macro `create_extended_executable`, which can be called as `create_extended_example(<driver file> <system file>)`.

Now, your code will be compiled when `make` is run from within the `build/` directory. The executable will be stored in the same relative path from `build/` as the source is from this directory. i.e., `build/examples/extended/...` and will have the same name as the driver file (without the `.f95` extension.


## Viewing solutions

A python library for handling the B-spline output of eBACOLI is found at [./lib/python/ebacoli.py](./lib/python/ebacoli.py).

If using bash, [./ebacoli.env](./ebacoli.env) sets the appropriate directory in your `PYTHONPATH` for its use. (Run with `source ebacoli.env`)

Example for how to structure the driver and post-processing to generate video frames is found in [./examples/movie](./examples/movie).

Visual results
--------------

Movies of various solutions found in [./media](./media).
