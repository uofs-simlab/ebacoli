
# Source files defining the system
set(TenTusscher_SRC
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

# Rule for creating ten-Tusscher executables
macro(create_TenTusscher_executable driver_file)
  add_executable(${driver_file}
    ${driver_file}.f95
    ${TenTusscher_SRC}
    )
  target_link_libraries(${driver_file} ebacoli)
endmacro(create_TenTusscher_executable)

create_TenTusscher_executable(driver-TT_monodomain-simple)
create_TenTusscher_executable(driver-TT_monodomain-trimesh)
