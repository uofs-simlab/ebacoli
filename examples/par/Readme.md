From Paul Muir's [website](http://cs.stmarys.ca/~muir/BACOLI-3_Webpage.htm).

Downloaded \[2017-05-31 Wed\]. -- Modified to run with the ebacoli solver.

**NOTE:** Files have been moved to capitals in file extension. This tells `gfortran` to run them through the C preprocessor before applying the fortran compiler. This is what allows us to compile executables for different sized systems with ease. System size is provided through the macros:
- `#define MAC_NU <SIZE>`,
- `#define MAC_NV 0`,
- `#define MAC_NW 0`.
Or by using our CMake macro: `create_macro_ebacoli_executable(<driver_file> <system_file> <size> 0 0)`.
