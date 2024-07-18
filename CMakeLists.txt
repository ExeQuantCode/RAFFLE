cmake_minimum_required(VERSION 3.17.5)

# define build environments
set( CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local/raffle"
  CACHE STRING "Select where to install the library." )
execute_process(COMMAND pwd OUTPUT_VARIABLE CURRENT_WORKING_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
message("Current working directory: ${CURRENT_WORKING_DIR}")
set( CMAKE_BUILD_PREFIX ${CURRENT_WORKING_DIR}
  CACHE STRING "Select where to build the library." )
set(MODULE_DIR ${CMAKE_BUILD_PREFIX}/mod)

# set the project name
project(raffle NONE)

# set the library name
set( LIB_NAME ${PROJECT_NAME} )
set( PROJECT_DESCRIPTION
  "Fortran neural network" )
set( PROJECT_URL "https://github.com/nedtaylor/raffle" )
set( CMAKE_CONFIGURATION_TYPES "Release" "Parallel" "Serial" "Dev" "Debug" "Parallel_Dev"
  CACHE STRING "List of configurations types." )
set( CMAKE_BUILD_TYPE "Release"
  CACHE STRING "Select which configuration to build." )

# change name based on parallel
if (CMAKE_BUILD_TYPE MATCHES "Parallel*")
   project(raffle_omp NONE) # change project name to parallel
   message(FATAL_ERROR "Configuration stopped because Parallel is not yet set up")
endif()

# set compiler
set(CMAKE_Fortran_COMPILER gfortran
  CACHE STRING "Select fortran compiler." )  # Change this to your desired compiler
set(CMAKE_Fortran_STANDARD 2018)

# set language
enable_language(Fortran)

# get the user's home directory
set(HOME_DIR $ENV{HOME})

# set coverage compiler flags
if (CMAKE_BUILD_TYPE MATCHES "Debug*" OR CMAKE_BUILD_TYPE MATCHES "Dev*")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
  set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      include(CodeCoverage)
      setup_target_for_coverage_gcovr_html(
        NAME coverage
        EXECUTABLE ctest
        EXCLUDE "${PROJECT_SOURCE_DIR}/test/*")
  endif()
endif()

# enable testing
enable_testing()

# set options for building tests and examples
option(BUILD_TESTS "Build the unit tests" ON)
option(BUILD_EXAMPLES "Build the examples" ON)
option(ENABLE_ATHENA "Build energetic predictions with ATHENA" OFF)
option(BUILD_PYTHON "Build the python library" On)

# Define the sources
set(SRC_DIR src)
set(LIB_DIR ${SRC_DIR}/lib)

set(LIB_FILES
  mod_constants.f90
  mod_misc.f90
	mod_misc_maths.f90
	mod_misc_linalg.f90
	# mod_rw_geom.f90
	mod_rw_vasprun.f90
	mod_edit_geom.f90
	mod_elements.f90
	mod_evolver.f90
	mod_buildmap.f90
	mod_atom_adder.f90
	mod_read_structures.f90
)

set(SPECIAL_LIB_FILES
  mod_rw_geom.f90
  mod_generator.f90
  # mod_generator_sub.f90
)


# Add Athena library
if(ENABLE_ATHENA)
    # Specify the paths to the include and library directories
    set(ATHENA_ROOT ${HOME_DIR}/.local/athena)
    set(ATHENA_INCLUDE_DIR ${HOME_DIR}/.local/athena/include)
    set(ATHENA_LIBRARY ${HOME_DIR}/.local/athena/lib/libathena.a)  # or libfoo.a for static libraries

    find_library(ATHENA_LIBRARY NAMES athena)

    # Add the include directory
    include_directories(${ATHENA_INCLUDE_DIR})
    list(APPEND LIB_FILES mod_ml.f90)

    set(F2PY_ATHENA_LIBRARY_FLAGS
        -L${ATHENA_ROOT}/lib
        -I${ATHENA_ROOT}/include
        -lathena
    )
endif()

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
    list(APPEND F90WRAP_FORTRAN_SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/${SRC_DIR}/${src})
    list(APPEND PREPENDED_SRC_FILES ${SRC_DIR}/${src})
endforeach()

message(STATUS "Modified SRC_FILES: ${PREPENDED_SRC_FILES}")

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
install(FILES ${MODULE_DIR}/${MODULE_FILES} DESTINATION include)
install(TARGETS ${PROJECT_NAME} DESTINATION lib)

# set compile options based on different build configurations
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Release>:${OPTIMFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Parallel>:${OPTIMFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Parallel>:${MPFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Dev>:${DEVFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Debug>:${DEBUGFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Parallel_Dev>:${MPFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Parallel_Dev>:${DEVFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Release>:${PYTHONFLAGS}>")


# add coverage compiler flags
if (CMAKE_BUILD_TYPE MATCHES "Debug*" OR CMAKE_BUILD_TYPE MATCHES "Dev*")
  append_coverage_compiler_flags()
endif()





if (BUILD_PYTHON)

    # # Get the directory where object files are generated
    get_target_property(OBJECTS ${PROJECT_NAME} EXTERNAL_OBJECT)
    # Print the object files directory
    set(OBJECTS_DIR ${CMAKE_BUILD_PREFIX}/CMakeFiles/${PROJECT_NAME}.dir)
    message(STATUS "Object files directory for ${PROJECT_NAME}: ${OBJECTS_DIR}")


    # Include f90wrap
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
    find_program(F90WRAP_EXECUTABLE f90wrap)
    find_program(F2PY_EXECUTABLE f2py-f90wrap)

    if(NOT F90WRAP_EXECUTABLE)
        message(FATAL_ERROR "f90wrap not found. Please install f90wrap.")
    endif()

    # Generate f90wrap signature file
    set(F90WRAP_FILE ${CMAKE_BINARY_DIR}/f90wrap_*.f90)
    set(KIND_MAP ${CMAKE_SOURCE_DIR}/kind_map)
    add_custom_command(
        TARGET ${PROJECT_NAME}
        POST_BUILD
        COMMAND ${F90WRAP_EXECUTABLE}
        --default-to-inout
        -m ${PROJECT_NAME}
        -k ${KIND_MAP}
        ${F90WRAP_FORTRAN_SRC_FILES}
        # --only raffle_generator_type:
        DEPENDS ${F90WRAP_FORTRAN_SRC_FILES}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Generating f90wrap signature file"
        VERBATIM
    )

    # Copy f90wrap edited files from edited_autogen_files to ${CMAKE_BINARY_DIR}
    add_custom_command(
        TARGET ${PROJECT_NAME}
        POST_BUILD
        COMMAND cp ${CMAKE_CURRENT_LIST_DIR}/edited_autogen_files/* ${CMAKE_BINARY_DIR}
        COMMENT "Copying f90wrap edited files"
    )

    # Create a Python module using f2py
    add_custom_command(
        TARGET ${PROJECT_NAME}
        POST_BUILD
        COMMAND ${F2PY_EXECUTABLE}
        ${F2PY_ATHENA_LIBRARY_FLAGS}
        -I${MODULE_DIR}
        -c
        -m _${PROJECT_NAME}
        --f90flags="${PPFLAGS}"
        ${F90WRAP_FILE}
        ${OBJECTS_DIR}/src/*.o
        ${OBJECTS_DIR}/src/lib/*.o
        DEPENDS ${F90WRAP_FILE}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Creating Python module using f2py"
    )


    # Define output files
    set(PY_MODULE ${CMAKE_BINARY_DIR}/${PROJECT_NAME}.py)
    file(GLOB SO_MODULE "${CMAKE_BINARY_DIR}/_${PROJECT_NAME}*.so")

    # Create a custom target for the Python module
    add_custom_target(python_module ALL
        DEPENDS ${SO_MODULE} ${PY_MODULE}
    )

    # Installation instructions
    install(FILES ${PY_MODULE} DESTINATION lib)
    install(FILES ${SO_MODULE} DESTINATION lib)
endif()



# Print helpful messages
message(STATUS "Build configuration:")
message(STATUS "  Source directory: ${SRC_DIR}")
message(STATUS "  Output library: ${PROJECT_NAME}")
# message(STATUS "  Fortran modules directory: ${CMAKE_BINARY_DIR}/mod")
# message(STATUS "  Python module: ${PROJECT_NAME}.so")



target_link_libraries(raffle ${ATHENA_LIBRARY})