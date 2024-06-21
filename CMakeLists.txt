cmake_minimum_required(VERSION 3.17.5)

# define build environments
set( CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local/raffle"
  CACHE STRING "Select where to install the library." )
execute_process(COMMAND pwd OUTPUT_VARIABLE CURRENT_WORKING_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
message("Current working directory: ${CURRENT_WORKING_DIR}")
set( CMAKE_BUILD_PREFIX ${CURRENT_WORKING_DIR}
  CACHE STRING "Select where to build the library." )

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

# Define the sources
set(SRC_DIR src)
set(LIB_DIR ${SRC_DIR}/lib)

set(LIB_FILES
    mod_constants.f90
    mod_misc.f90
	mod_misc_maths.f90
	mod_misc_linalg.f90
	mod_rw_geom.f90
	mod_rw_vasprun.f90
	mod_edit_geom.f90
	mod_elements.f90
	mod_ml.f90
	mod_evolver.f90
	mod_buildmap.f90
	mod_atom_adder.f90
	mod_read_structures.f90
)
foreach(lib ${LIB_FILES})
    list(APPEND PREPENDED_LIB_FILES ${LIB_DIR}/${lib})
endforeach()
message(STATUS "Modified LIB_FILES: ${PREPENDED_LIB_FILES}")

if (CMAKE_BUILD_TYPE MATCHES "Parallel*")
    set(SRC_FILES raffle_omp.f90)
else()
    set(SRC_FILES inputs.f90 generator.f90 raffle.f90)
endif()
foreach(src ${SRC_FILES})
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


# Get the user's home directory
set(HOME_DIR $ENV{HOME})

# Specify the paths to the include and library directories
set(ATHENA_INCLUDE_DIR ${HOME_DIR}/.local/athena/include)
set(ATHENA_LIBRARY ${HOME_DIR}/.local/athena/lib/libathena.a)  # or libfoo.a for static libraries

find_library(ATHENA_LIBRARY NAMES athena)

# Add the include directory
include_directories(${ATHENA_INCLUDE_DIR})



set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${PPFLAGS}")


# create the library
add_library(${PROJECT_NAME} STATIC ${PREPENDED_LIB_FILES} ${PREPENDED_SRC_FILES})
set(MODULE_DIR ${CMAKE_BUILD_PREFIX}/modules)
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


# add coverage compiler flags
if (CMAKE_BUILD_TYPE MATCHES "Debug*" OR CMAKE_BUILD_TYPE MATCHES "Dev*")
  append_coverage_compiler_flags()
endif()

target_link_libraries(raffle ${ATHENA_LIBRARY})