
# Rule for creating and properly linking executables in these subdirectories
macro(create_test_executable problem_file)
  add_executable(test_${problem_file}
    test_${problem_file}.f95
    ${problem_file}.f
    utilities.f95
    )
  target_link_libraries(test_${problem_file} ebacoli)
  add_test(${problem_file} test_${problem_file})
endmacro(create_test_executable)

############################################################################
## Executables for all tests

add_subdirectory(par)
add_subdirectory(par-ode)
add_subdirectory(par-ellip)
add_subdirectory(par-ode-ellip)
