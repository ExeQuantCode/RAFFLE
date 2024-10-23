message(STATUS "Building tests")
foreach(execid
    misc
    misc_maths
    misc_linalg
    elements
    rw_geom
    geom_utils
    dist_calcs
    extended_geom
    atom_adder
    distribs_container
    evaluator_C
    evaluator_BTO
    generator
  )
  add_executable(test_${execid} test_${execid}.f90)
  # # Specify the include directories
  target_include_directories(test_${execid} PRIVATE "${CMAKE_BUILD_PREFIX}" "${CMAKE_BUILD_PREFIX}/mod")

  target_link_libraries(test_${execid} PRIVATE raffle)

  add_test(NAME test_${execid} COMMAND test_${execid})
endforeach()