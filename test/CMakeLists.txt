message(STATUS "Building tests")
foreach(execid
    misc
    misc_maths
    misc_linalg
    element_utils
    geom_rw
    geom_utils
    geom_extd
    dist_calcs
    place_methods
    viability
    distribs_container
    evaluator_C
    evaluator_BTO
    generator
    io_utils
    tools_infile
  )
  add_executable(test_${execid} test_${execid}.f90)
  # # Specify the include directories
  target_include_directories(test_${execid} PRIVATE "${CMAKE_BUILD_DIR}" "${CMAKE_BUILD_DIR}/mod")

  target_link_libraries(test_${execid} PRIVATE ${PROJECT_NAME})

  add_test(NAME test_${execid} COMMAND test_${execid})
endforeach()