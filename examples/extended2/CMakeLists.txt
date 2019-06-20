
# Rule for creating and properly linking the executables in this directory
macro(create_extended_executable driver_file system_file)
  add_executable(${driver_file}
    ${driver_file}.f95
    ${system_file}.f
    )
  target_link_libraries(${driver_file} ebacoli)
endmacro(create_extended_executable)

############################################################################
# Toy problem examples

############################################################################
# FHN Bidomain examples

create_extended_executable(driver-FHN2_bidomain-simple FHN2_bidomain)
