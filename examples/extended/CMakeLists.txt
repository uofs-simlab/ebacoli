
# Rule for creating and properly linking the executables in this directory
macro(create_extended_executable driver_file system_file)
  add_executable(${driver_file}
    ${driver_file}.f95
    ${system_file}.f
    )
  target_link_libraries(${driver_file} ebacoli)
endmacro(create_extended_executable)

############################################################################
# FHN Monodomain examples

create_extended_executable(driver-FHN_monodomain-simple FHN_monodomain)
create_extended_executable(driver-FHN2_monodomain-simple FHN2_monodomain)

create_extended_executable(driver-FHN_monodomain-trimesh FHN_monodomain)
create_extended_executable(driver-FHN2_monodomain-trimesh FHN2_monodomain)

############################################################################
# Competition-diffusion examples

create_extended_executable(driver-CDM-simple CDM)

create_extended_executable(driver-CDM-trimesh CDM)

############################################################################
# Rabies example

create_extended_executable(driver-rabid-trimesh rabid)

############################################################################
# ten-Tusscher example

add_subdirectory(TenTusscher)
