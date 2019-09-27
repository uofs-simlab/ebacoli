cmake_minimum_required(VERSION 2.6)

project(eBACOLI)
enable_language(Fortran)

# Be able to find the module file produced in building the library
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

enable_testing()

# Useful macros and functions for creating and linking new executables

# Rule for creating and properly linking the executables in this directory
macro(create_ebacoli_executable driver_file system_file)
  add_executable(${driver_file}
    ${driver_file}.f95
    ${system_file}.f
    )
  target_link_libraries(${driver_file} ebacoli)
endmacro(create_ebacoli_executable)


# subdirectories
add_subdirectory(src)
add_subdirectory(examples)

add_subdirectory(tests)
