cmake_minimum_required(VERSION 2.6)

project(eBACOLI)
enable_language(Fortran)

# Be able to find the module file produced in building the library
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

add_subdirectory(src)
add_subdirectory(examples)
