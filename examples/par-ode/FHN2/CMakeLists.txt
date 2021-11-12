############################################################################
# Modified FHN Monodomain example

set(TARGET FHN2_monodomain)

set(TARGET_SRC
  driver-FHN2_monodomain-trimesh.f95
  FHN2_monodomain.f
  )

set(TARGET_EXE FHN2_monodomain)

############################################################################
############################################################################
# This section is generally automatic with the above defined variables
cmake_minimum_required(VERSION 3.0)

# EBACOLI_LIB_COMPILE is only defined when building the lib
if(NOT DEFINED EBACOLI_LIB_COMPILE)
  # define these which would otherwise be defined in the library build
  PROJECT(${TARGET})
  enable_language(Fortran)

  # Find the installed ebacoli and print status
  find_package(ebacoli 0.1 REQUIRED
    HINTS ${EBACOLI_DIR} ${ebacoli_DIR} ${eBACOLI_DIR} ../ ../../ $ENV{EBACOLI_DIR} $ENV{ebacoli_DIR} $ENV{eBACOLI_DIR}
    )
  get_target_property(libloc ebacoli IMPORTED_LOCATION_RELEASE)
  message(STATUS "Found libebacoli: ${libloc}")
else()
  # specify an install target if built as part of the library
  install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_EXE}" DESTINATION ${EBACOLI_RUNTIME_INSTALL_DIR})
endif()

add_executable(${TARGET_EXE} ${TARGET_SRC})
target_compile_definitions(${TARGET_EXE} PUBLIC ${TARGET_COMPILE_OPTIONS})
target_link_libraries(${TARGET_EXE} ebacoli)
############################################################################
############################################################################
