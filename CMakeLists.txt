cmake_minimum_required(VERSION 3.0)

project(EBACOLI)
enable_language(Fortran)

# define library version
set(EBACOLI_VERSION_MAJOR 0 CACHE STRING "major version of ebacoli" FORCE)
set(EBACOLI_VERSION_MINOR 1 CACHE STRING "minor version of ebacoli" FORCE)
set(EBACOLI_VERSION ${EBACOLI_VERSION_MAJOR}.${EBACOLI_VERSION_MINOR} CACHE STRING "version of ebacoli" FORCE)

# set a debug postfix
set(CMAKE_DEBUG_POSTFIX "-dbg")

# Module location at compile time
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# installation destinations
if(UNIX OR VXWORKS)
    include(GNUInstallDirs)

    set(EBACOLI_INC_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}/ebacoli")
    set(EBACOLI_RUNTIME_INSTALL_DIR "${CMAKE_INSTALL_BINDIR}/ebacoli")
    set(EBACOLI_LIBRARY_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}")
    set(EBACOLI_ARCHIVE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}")
    set(EBACOLI_FRAMEWORK_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}")

    set(EBACOLI_CMAKE_CONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/ebacoli/cmake")
    set(EBACOLI_ADDITIONAL_FILES_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/ebacoli")
else()
    message(FATAL_ERROR "Could not set install folders for this platform!")
endif()

# Option to disable compiling pieces beyond the libary
option(EBACOLI_BUILD_EXAMPLES "Build the example codes" ON)
option(EBACOLI_BUILD_TESTS "Build the test codes" ON)

# Avoid verbose "Up-to-date" status information during installation:
#
MACRO(SET_IF_EMPTY _variable)
  IF("${${_variable}}" STREQUAL "")
    SET(${_variable} ${ARGN})
  ENDIF()
ENDMACRO()
SET_IF_EMPTY(CMAKE_INSTALL_MESSAGE "LAZY")

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

# create a dummy variable for source build
set(EBACOLI_LIB_COMPILE "")

# subdirectories
add_subdirectory(src)
if(EBACOLI_BUILD_EXAMPLES)
  add_subdirectory(examples)
  message(STATUS "ebacoli example codes to be installed at:
     ${CMAKE_INSTALL_PREFIX}/${EBACOLI_ADDITIONAL_FILES_INSTALL_DIR}")
  message(STATUS "ebacoli example binaries to be installed at:
     ${CMAKE_INSTALL_PREFIX}/${EBACOLI_RUNTIME_INSTALL_DIR}")
  install(DIRECTORY "${CMAKE_SOURCE_DIR}/examples" DESTINATION ${EBACOLI_ADDITIONAL_FILES_INSTALL_DIR})
endif()
if(EBACOLI_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE" "${CMAKE_CURRENT_SOURCE_DIR}/README.md" DESTINATION ${EBACOLI_ADDITIONAL_FILES_INSTALL_DIR})

install(EXPORT ebacoliTargets DESTINATION ${EBACOLI_CMAKE_CONFIG_INSTALL_DIR}
                                       FILE ebacoli-config.cmake)

# install the python tools
install(DIRECTORY ${CMAKE_SOURCE_DIR}/lib/python DESTINATION ${EBACOLI_LIBRARY_INSTALL_DIR})
