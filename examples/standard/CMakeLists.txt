
# Rule for creating and properly linking the executables in this directory
macro(create_standard_executable driver_file system_file number_of_fields)
  add_executable(${driver_file}-${system_file}
    ${driver_file}.F95
    ${system_file}.F
    )
  target_compile_definitions(${driver_file}-${system_file} PUBLIC "-DMAC_NPDE=${number_of_fields}")
  target_link_libraries(${driver_file}-${system_file} ebacoli)
endmacro(create_standard_executable)

############################################################################
# Burgers equation examples

create_standard_executable(driver95-simple burg1 1)
create_standard_executable(driver95-simple burg2 1)

create_standard_executable(driver95-curve burg1 1)
create_standard_executable(driver95-curve burg2 1)

create_standard_executable(driver95-trimesh burg1 1)
create_standard_executable(driver95-trimesh burg2 1)

############################################################################
# "Steady" examples

create_standard_executable(driver95-simple steady 1)
create_standard_executable(driver95-curve steady 1)
create_standard_executable(driver95-trimesh steady 1)

############################################################################
# Cahn--Allen equation examples

create_standard_executable(driver95-simple cahn_allen 1)
create_standard_executable(driver95-curve cahn_allen 1)
create_standard_executable(driver95-trimesh cahn_allen 1)

############################################################################
# Reaction-diffusion system examples

create_standard_executable(driver95-simple rcdsys 4)
create_standard_executable(driver95-curve rcdsys 4)
create_standard_executable(driver95-trimesh rcdsys 4)

############################################################################
