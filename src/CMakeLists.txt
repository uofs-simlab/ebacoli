
set(src
  ebacoli95.f95
  ebacoli.f
  ebacoli-aux.f
  d1mach_i1mach.f95
  abdtools.f
  DASUM.f
  utilities.f95
)

add_library(ebacoli ${src})
target_compile_definitions(ebacoli PUBLIC
                            EBACOLI_VERSION_MAJOR=${EBCOLI_VERSION_MAJOR}
                            EBACOLI_VERSION_MINOR=${EBACOLI_VERSION_MINOR})

set_target_properties(ebacoli PROPERTIES
                       OUTPUT_NAME "ebacoli-${EBACOLI_VERSION}"
                       POSITION_INDEPENDENT_CODE ON)

# Write/install version file
include(CMakePackageConfigHelpers)
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/cmake/ebacoli-config-version.cmake")
write_basic_package_version_file(${version_file}
                                 VERSION ${EBACOLI_VERSION}
                                 COMPATIBILITY AnyNewerVersion)

install(TARGETS ebacoli EXPORT ebacoliTargets
        RUNTIME       DESTINATION ${EBACOLI_RUNTIME_INSTALL_DIR}
        LIBRARY       DESTINATION ${EBACOLI_LIBRARY_INSTALL_DIR}
        ARCHIVE       DESTINATION ${EBACOLI_ARCHIVE_INSTALL_DIR}
        FRAMEWORK     DESTINATION ${EBACOLI_FRAMEWORK_INSTALL_DIR})

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${EBACOLI_INC_INSTALL_DIR})
