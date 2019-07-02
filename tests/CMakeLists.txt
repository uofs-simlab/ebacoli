
# Rule for creating and properly linking the executables in this directory
macro(create_test_executable problem_file)
  add_executable(test_${problem_file}
    test_${problem_file}.f95
    ${problem_file}.f
    utilities.f95
    )
  target_link_libraries(test_${problem_file} ebacoli)
endmacro(create_test_executable)

############################################################################
## Executables for all tests

## tests of standard BACOLI
create_test_executable(burg1)
create_test_executable(burg2)

create_test_executable(sincmads)

## tests of eBACOLI
# parabolic + ODE


# parabolic + elliptic


# parabolic + ODE + elliptic
