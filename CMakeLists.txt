cmake_minimum_required(VERSION 3.17.5...3.27.7)

# define build environments
set( CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local/"
  CACHE STRING "Select where to install the library." )
set(CMAKE_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}
  CACHE STRING "Select where to build the library." )
set(MODULE_DIR ${CMAKE_BUILD_DIR}/mod)

if (DEFINED SKBUILD_PROJECT_NAME)
	set(SKBUILD_PROJECT_NAME ${SKBUILD_PROJECT_NAME})
else()
	set(SKBUILD_PROJECT_NAME "raffle")
endif()


# set compiler
set(CMAKE_Fortran_COMPILER gfortran
  CACHE STRING "Select Fortran compiler." )  # Change this to your desired compiler
set(CMAKE_C_COMPILER gcc
  CACHE STRING "Select C compiler." )  # Change this to your desired compiler
set(CMAKE_Fortran_STANDARD 2018)

# set the project name
project(raffle LANGUAGES C Fortran)

# set the library name
set( LIB_NAME ${PROJECT_NAME} )
set( PROJECT_DESCRIPTION
  "Fortran neural network" )
set( PROJECT_URL "https://github.com/nedtaylor/raffle" )
set( CMAKE_CONFIGURATION_TYPES "Release" "Parallel" "Serial" "Dev" "Debug" "Coverage"
  CACHE STRING "List of configurations types." )
set( CMAKE_BUILD_TYPE "Release"
  CACHE STRING "Select which configuration to build." )

# change name based on parallel
if (CMAKE_BUILD_TYPE MATCHES "Parallel*")
   project(raffle_omp NONE) # change project name to parallel
   message(FATAL_ERROR "Configuration stopped because Parallel is not yet set up")
endif()

# set options for building tests and examples
option(BUILD_TESTS "Build the unit tests" On)
option(BUILD_EXAMPLES "Build the examples" On)
option(BUILD_PYTHON "Build the python library" On)
option(BUILD_EXECUTABLE "Build the Fortran executable" On)
option(REMAKE_F90WRAP "Remake the f90wrap signature file" Off)

# set coverage compiler flags
if (CMAKE_BUILD_TYPE MATCHES "Coverage")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
  set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      include(CodeCoverage)
      setup_target_for_coverage_gcovr_html(
        NAME coverage
        EXECUTABLE ctest
        EXCLUDE "${CMAKE_SOURCE_DIR}/test/*")
  endif()
endif()

# enable testing
enable_testing()

# Define the sources
set(SRC_DIR src)
set(FORTRAN_SRC_DIR ${SRC_DIR}/fortran)
set(LIB_DIR ${FORTRAN_SRC_DIR}/lib)

set(LIB_FILES
	mod_io_utils.F90
	mod_constants.f90
	mod_misc.f90
    mod_tools_infile.f90
	mod_misc_maths.f90
	mod_misc_linalg.f90
	mod_dist_calcs.f90
	mod_geom_utils.f90
	mod_geom_extd.f90
	mod_element_utils.f90
	mod_evaluator.f90
	mod_place_methods.f90
	mod_viability.f90
	mod_distribs.f90
	mod_distribs_host.f90
)

set(SPECIAL_LIB_FILES
	mod_geom_rw.f90
	mod_distribs_container.f90
	mod_generator.f90
)


foreach(lib ${LIB_FILES})
    list(APPEND PREPENDED_LIB_FILES ${LIB_DIR}/${lib})
endforeach()
foreach(lib ${SPECIAL_LIB_FILES})
    list(APPEND PREPENDED_LIB_FILES ${LIB_DIR}/${lib})
endforeach()
message(STATUS "Modified LIB_FILES: ${PREPENDED_LIB_FILES}")



set(SRC_FILES
    raffle.f90
)
foreach(lib ${SPECIAL_LIB_FILES})
    list(APPEND F90WRAP_FORTRAN_SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/${LIB_DIR}/${lib})
endforeach()
foreach(src ${SRC_FILES})
    list(APPEND F90WRAP_FORTRAN_SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/${FORTRAN_SRC_DIR}/${src})
    list(APPEND PREPENDED_SRC_FILES ${FORTRAN_SRC_DIR}/${src})
endforeach()


set(EXECUTABLE_FILES
    mod_rw_vasprun.f90
    mod_read_structures.f90
    inputs.f90
    main.f90
)
set(APP_DIR app)
foreach(src ${EXECUTABLE_FILES})
    list(APPEND PREPENDED_EXECUTABLE_FILES ${APP_DIR}/${src})
endforeach()


# initialise flags
set(CPPFLAGS "")
set(CFLAGS "")
set(MODULEFLAGS "")
set(MPFLAGS "")
set(WARNFLAGS "")
set(DEVFLAGS "")
set(DEBUGFLAGS "")
set(MEMFLAGS "")
set(OPTIMFLAGS "")
set(FASTFLAGS "")

# set flags based on compiler
if (CMAKE_Fortran_COMPILER MATCHES ".*gfortran.*" OR CMAKE_Fortran_COMPILER MATCHES ".*gcc.*")
    message(STATUS "Using gfortran compiler")
    set(PPFLAGS -cpp)
    set(MPFLAGS -fopenmp)
    set(WARNFLAGS -Wall)
    set(DEVFLAGS -g -fbacktrace -fcheck=all -fbounds-check -Og)
    set(DEBUGFLAGS -fbounds-check)
    set(MEMFLAGS -mcmodel=large)
    set(OPTIMFLAGS -O3 -march=native)
    set(FASTFLAGS -Ofast -march=native)
    set(PYTHONFLAGS -c -O3 -fPIC)
elseif (CMAKE_Fortran_COMPILER MATCHES ".*nag.*")
    message(STATUS "Using nag compiler")
    set(PPFLAGS -f2018 -fpp)
    set(MPFLAGS -openmp)
    set(WARNFLAGS -Wall)
    set(DEVFLAGS -g -mtrace -C=all -colour -O0)
    set(DEBUGFLAGS -C=array)
    set(MEMFLAGS -mcmodel=large)
    set(OPTIMFLAGS -O3)
    set(FASTFLAGS -Ofast)
elseif (CMAKE_Fortran_COMPILER MATCHES ".*ifort.*" OR CMAKE_Fortran_COMPILER MATCHES ".*ifx.*")
    message(STATUS "Using intel compiler")
    set(PPFLAGS -fpp)
    set(MPFLAGS -qopenmp)
    set(WARNFLAGS -warn all)
    set(DEVFLAGS -check all -warn)
    set(DEBUGFLAGS -check all -fpe0 -warn -tracekback -debug extended)
    set(MEMFLAGS -mcmodel=large)
    set(OPTIMFLAGS -O3)
    set(FASTFLAGS -Ofast)
else()
    # Code for other Fortran compilers
    message(STATUS "Using a different Fortran compiler")
endif()



set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${PPFLAGS}")


# create the library
add_library(${PROJECT_NAME} STATIC ${PREPENDED_LIB_FILES} ${PREPENDED_SRC_FILES})
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC)

# replace ".f90" with ".mod"
string(REGEX REPLACE "\\.[^.]*$" ".mod" MODULE_FILES "${SRC_FILES}")

# installation
install(FILES ${MODULE_DIR}/${MODULE_FILES} DESTINATION ${SKBUILD_PROJECT_NAME}/include)
install(TARGETS ${PROJECT_NAME} DESTINATION ${SKBUILD_PROJECT_NAME}/lib)

# set compile options based on different build configurations
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Release>:${OPTIMFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Parallel>:${OPTIMFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Parallel>:${MPFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Dev>:${DEVFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Debug>:${DEBUGFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Debug>:${WARNFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Debug>:${PYTHONFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Release>:${PYTHONFLAGS}>")



if (BUILD_EXECUTABLE)
    add_executable(${PROJECT_NAME}_executable ${PREPENDED_EXECUTABLE_FILES})
    target_link_libraries(${PROJECT_NAME}_executable PRIVATE ${PROJECT_NAME})
    install(TARGETS ${PROJECT_NAME}_executable DESTINATION ${SKBUILD_PROJECT_NAME}/bin)
    set_target_properties(${PROJECT_NAME}_executable PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIR})
    target_compile_options(${PROJECT_NAME}_executable PUBLIC "$<$<CONFIG:Release>:${OPTIMFLAGS}>")
    target_compile_options(${PROJECT_NAME}_executable PUBLIC "$<$<CONFIG:Parallel>:${OPTIMFLAGS}>")
    target_compile_options(${PROJECT_NAME}_executable PUBLIC "$<$<CONFIG:Parallel>:${MPFLAGS}>")
    target_compile_options(${PROJECT_NAME}_executable PUBLIC "$<$<CONFIG:Dev>:${DEVFLAGS}>")
    target_compile_options(${PROJECT_NAME}_executable PUBLIC "$<$<CONFIG:Debug>:${DEBUGFLAGS}>")
    target_compile_options(${PROJECT_NAME}_executable PUBLIC "$<$<CONFIG:Release>:${PYTHONFLAGS}>")
endif()



if (BUILD_PYTHON)

    # # Get the directory where object files are generated
    get_target_property(OBJECTS ${PROJECT_NAME} EXTERNAL_OBJECT)
    # Print the object files directory
    set(OBJECTS_DIR ${CMAKE_BUILD_DIR}/CMakeFiles/${PROJECT_NAME}.dir)
    message(STATUS "Object files directory for ${PROJECT_NAME}: ${OBJECTS_DIR}")

    # Include f90wrap
    find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
    if(NOT DEFINED PYTHON_EXECUTABLE)
        set(PYTHON_EXECUTABLE ${Python_EXECUTABLE})
    endif()

    # Grab the variables from a local Python installation F2PY headers
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -c
                "import numpy.f2py; print(numpy.f2py.get_include())"
        OUTPUT_VARIABLE F2PY_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
    target_link_libraries(fortranobject PUBLIC Python::NumPy)
    target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
    set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

    set (F90WRAP_EXECUTABLE ${PYTHON_EXECUTABLE} -m f90wrap)
    set (F2PY_EXECUTABLE ${PYTHON_EXECUTABLE} -m f90wrap --f2py-f90wrap)

    # Run Python command to get the extension suffix
    execute_process(
        COMMAND ${Python_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE PYTHON_EXTENSION_MODULE_SUFFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Check if the suffix was retrieved successfully
    if (result EQUAL 0)
        message(STATUS "Python extension module suffix: ${PYTHON_EXTENSION_MODULE_SUFFIX}")
    else()
        message(FATAL_ERROR "Failed to retrieve Python extension module suffix")
    endif()
    set(F2PY_OUTPUT_FILE ${CMAKE_BUILD_DIR}/raffle/_${PROJECT_NAME}${PYTHON_EXTENSION_MODULE_SUFFIX})

    # Generate f90wrap signature file
    set(F90WRAP_FILE
        ${CMAKE_CURRENT_LIST_DIR}/src/wrapper/f90wrap_mod_distribs_container.f90
        ${CMAKE_CURRENT_LIST_DIR}/src/wrapper/f90wrap_mod_generator.f90
        ${CMAKE_CURRENT_LIST_DIR}/src/wrapper/f90wrap_mod_geom_rw.f90
        ${CMAKE_CURRENT_LIST_DIR}/src/wrapper/f90wrap_raffle.f90
    )
    if (REMAKE_F90WRAP)
        set(KIND_MAP ${CMAKE_SOURCE_DIR}/kind_map)
        set(F90WRAP_REMAKE_FILE ${CMAKE_BUILD_DIR}/f90wrap_${PROJECT_NAME}.f90)
        add_custom_command(
            OUTPUT ${F90WRAP_REMAKE_FILE}
            POST_BUILD
            COMMAND ${F90WRAP_EXECUTABLE}
            --default-to-inout
            -m ${PROJECT_NAME}
            -k ${KIND_MAP}
            ${F90WRAP_FORTRAN_SRC_FILES}
            --only raffle_generator_type stoichiometry_type basis_type distribs_container_type:
            DEPENDS ${F90WRAP_FORTRAN_SRC_FILES}
            WORKING_DIRECTORY ${CMAKE_BUILD_DIR}
            COMMENT "Generating f90wrap signature file"
            VERBATIM
        )
    endif()

    # Copy f90wrap edited files from edited_autogen_files to ${CMAKE_BUILD_DIR}
    add_custom_command(
        OUTPUT ${CMAKE_BUILD_DIR}/raffle/python_copied
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/${SRC_DIR}/raffle ${CMAKE_BUILD_DIR}/raffle
        COMMENT "Copying raffle python files"
    )

    # Create a Python module using f2py
    add_custom_command(
        OUTPUT ${F2PY_OUTPUT_FILE}
        POST_BUILD
        COMMAND CC="${CMAKE_C_COMPILER}" LDFLAGS="-L${CMAKE_BUILD_DIR}" LIBS="-lraffle" ${F2PY_EXECUTABLE}
        -c
        -m _${PROJECT_NAME}
        -I${MODULE_DIR}
        --f90flags="${PPFLAGS}"
        --backend meson
        ${F90WRAP_FILE}
        -L${CMAKE_BUILD_DIR}
        -lraffle
        # --build-dir ${CMAKE_BUILD_DIR}/Dmeson
        DEPENDS ${F90WRAP_FILE} ${CMAKE_BUILD_DIR}/raffle/python_copied ${PROJECT_NAME}
        WORKING_DIRECTORY ${CMAKE_BUILD_DIR}/raffle
        COMMENT "Creating Python module using f2py"
    )

    # Define output files
	set(PY_MODULE ${CMAKE_BUILD_DIR}/raffle/${PROJECT_NAME}.py ${CMAKE_BUILD_DIR}/raffle/__init__.py)

    # Create a custom target for the Python module
    add_custom_target(python_module ALL
        DEPENDS ${F90WRAP_REMAKE_FILE} ${F2PY_OUTPUT_FILE}
    )

    # Installation instructions
    install(FILES ${PY_MODULE} DESTINATION ${SKBUILD_PROJECT_NAME})
    install(FILES ${F2PY_OUTPUT_FILE} DESTINATION ${SKBUILD_PROJECT_NAME})

endif()


# include the build test directory
if(BUILD_TESTS)
  add_subdirectory(test)
endif()

# add coverage compiler flags
if ( ( CMAKE_Fortran_COMPILER MATCHES ".*gfortran.*" OR CMAKE_Fortran_COMPILER MATCHES ".*gcc.*" ) AND 
     ( CMAKE_BUILD_TYPE MATCHES "Coverage" ) )

    append_coverage_compiler_flags()

endif()

if(BUILD_TESTS)
    # Create a directory for test data
    add_custom_command(
        OUTPUT ${CMAKE_BUILD_DIR}/test/test_dir_created  # Temporary marker file to indicate the directory was created
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BUILD_DIR}/test/test
        COMMENT "Creating test directory"
    )

    # Copy test data files
    add_custom_command(
        OUTPUT ${CMAKE_BUILD_DIR}/test/test_data_copied  # Temporary marker file to indicate the data was copied
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/test/data ${CMAKE_BUILD_DIR}/test/test/data
        DEPENDS ${CMAKE_CURRENT_LIST_DIR}/test/data ${CMAKE_BUILD_DIR}/test/test_dir_created  # Ensure the source data exists and the directory was created
        COMMENT "Copying test data files"
    )

    # Create a custom target to ensure that the copy commands are executed
    add_custom_target(copy_test_data ALL
        DEPENDS ${CMAKE_BUILD_DIR}/test/test_data_copied  # Ensure it depends on the data copying command
    )

    # Add dependency to your main project target
    add_dependencies(${PROJECT_NAME} copy_test_data)

    # Ensure that ctest waits for copy_test_data
    add_custom_command(TARGET copy_test_data POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E echo "Test data files copied successfully"
    )

    # Get all registered tests
    get_property(all_tests GLOBAL PROPERTY TESTS)

    # Loop through all tests and add the dependency
    foreach(test_name IN LISTS all_tests)
        add_dependencies(${test_name} copy_test_data)
    endforeach()

endif()