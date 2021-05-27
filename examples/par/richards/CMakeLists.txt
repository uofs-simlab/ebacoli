##
# Richards equation example
##

set(TARGET richards)

set(TARGET_SRC
  driver-richards-frameoutput.f95
  richards.f
  )

set(TARGET_EXE driver-richards-frameoutput)

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
  PROJECT(${TARGET})
  enable_language(Fortran)
endif()

add_executable(${TARGET_EXE} ${TARGET_SRC})
target_link_libraries(${TARGET_EXE} ebacoli)
############################################################################
############################################################################
