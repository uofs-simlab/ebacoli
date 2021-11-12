############################################################################
# Reaction-diffusion system examples

set(TARGET rcdsys)

set(TARGET_SRC
  driver95-trimesh.F95
  rcdsys.F
  )

set(TARGET_EXE rcdsys-trimesh)

set(TARGET_COMPILE_OPTIONS
"-DMAC_NU=4"
"-DMAC_NV=0"
"-DMAC_NW=0")

############################################################################
############################################################################
# This section is generally automatic with the above defined variables
cmake_minimum_required(VERSION 3.0)

# EBACOLI_LIB_COMPILE is only defined when building the lib
if(NOT DEFINED EBACOLI_LIB_COMPILE)
  # Find the installed ebacoli
  find_package(ebacoli 0.1 REQUIRED
    HINTS ${EBACOLI_DIR} ${ebacoli_DIR} ${eBACOLI_DIR} ../ ../../ $ENV{EBACOLI_DIR} $ENV{ebacoli_DIR} $ENV{eBACOLI_DIR}
    )
  # define these which would otherwise be defined in the library build
  PROJECT(${TARGET})
  enable_language(Fortran)
else()
  # specify an install target if built as part of the library
  install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_EXE}" DESTINATION ${EBACOLI_RUNTIME_INSTALL_DIR})
endif()

add_executable(${TARGET_EXE} ${TARGET_SRC})
target_compile_definitions(${TARGET_EXE} PUBLIC ${TARGET_COMPILE_OPTIONS})
target_link_libraries(${TARGET_EXE} ebacoli)
############################################################################
############################################################################
