cmake_minimum_required(VERSION 2.6)

project(eBACOLI)
enable_language(Fortran)

# Be able to find the module file produced in building the library
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Option to disable compiling pieces beyond the libary
option(EBACOLI_BUILD_EXAMPLES "Build the example codes" ON)
option(EBACOLI_BUILD_TESTS "Build the test codes" ON)


############################################################################
# Useful macros and functions for creating and linking new executables

# Rule for creating and properly linking executables with eBACOLI
macro(create_ebacoli_executable driver_file system_file)
  add_executable(${driver_file}
    ${driver_file}.f95
    ${system_file}.f
    )
  target_link_libraries(${driver_file} ebacoli)
endmacro(create_ebacoli_executable)

# Rule for creating an executable that relies on macro definitions for numbers of fields
# - Because Fortran (note works with gfortran thanks to gnu C preprocessor, not sure about other compilers).
# - Intent is to have driver files that can be written independent of the system being solved.
#   - This does not allow for generic inclusion of system parameters (in common blocks for example).
# - Note the capitalization of the file extensions
macro(create_macro_ebacoli_executable driver_file system_file nu nv nw)
  add_executable(${driver_file}-${system_file}
    ${driver_file}.F95
    ${system_file}.F
    )
  target_compile_definitions(${driver_file}-${system_file} PUBLIC
    "-DMAC_NU=${nu}"
    "-DMAC_NV=${nv}"
    "-DMAC_NW=${nw}")
  target_link_libraries(${driver_file}-${system_file} ebacoli)
endmacro(create_macro_ebacoli_executable)

############################################################################

# subdirectories
add_subdirectory(src)
if(EBACOLI_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
if(EBACOLI_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
