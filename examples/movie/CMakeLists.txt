
add_executable(CDM-frameoutput
  driver-CDM-frameoutput.f95
  CDM.f
  )
target_link_libraries(CDM-frameoutput ebacoli)

# TODO find a clean way to run the movie generation commands
