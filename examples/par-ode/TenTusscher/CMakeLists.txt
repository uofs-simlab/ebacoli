############################################################################
## ten Tusscher Monodomain example

set(TARGET TT-monodomain)

set(TARGET_SRC
  driver-TT_monodomain-trimesh.f95
  TenTusscher.f
  TT_ebacoli.f
  f_d_COL10.f
  f_d_COL12.f
  f_d_COL14.f
  f_d_COL16.f
  f_d_COL18.f
  f_d_COL2.f
  f_d_COL4.f
  f_d_COL6.f
  f_d_COL8.f
  f_d_COL1.f
  f_d_COL11.f
  f_d_COL13.f
  f_d_COL15.f
  f_d_COL17.f
  f_d_COL19.f
  f_d_COL3.f
  f_d_COL5.f
  f_d_COL7.f
  f_d_COL9.f
  )

set(TARGET_EXE TT-monodomain-trimesh)

############################################################################
############################################################################
# This section is generally automatic with the above defined variables
cmake_minimum_required(VERSION 3.0)

# EBACOLI_LIB_COMPILE is only defined when building the lib
if(NOT DEFINED EBACOLI_LIB_COMPILE)
  # Find the installed ebacoli and print status
  find_package(ebacoli 0.1 REQUIRED
    HINTS ${EBACOLI_DIR} ${ebacoli_DIR} ${eBACOLI_DIR} ../ ../../ $ENV{EBACOLI_DIR} $ENV{ebacoli_DIR} $ENV{eBACOLI_DIR}
    )
  get_target_property(libloc ebacoli IMPORTED_LOCATION_RELEASE)
  message(STATUS "Found libebacoli: ${libloc}")

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

#create_TenTusscher_executable(driver-TT_monodomain-simple)
#create_TenTusscher_executable(driver-TT_monodomain-trimesh)
