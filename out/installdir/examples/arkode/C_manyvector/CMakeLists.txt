# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos and David J. Gardner @ LLNL
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------

# Set the minimum required cmake version
cmake_minimum_required(VERSION 3.26.3)

# Set cache variables for C compilers and flags
set(CMAKE_C_COMPILER
  "/Volumes/Workspace/spack/opt/spack/darwin-ventura-skylake/apple-clang-14.0.3/gcc-13.1.0-zexv4j6ljav6xmm3ajnobit6fdhvxwrn/bin/gcc"
  CACHE FILEPATH "C compiler")

set(CMAKE_C_FLAGS
  ""
  CACHE STRING "C compiler flags")

set(CMAKE_C_STANDARD
  "99"
  CACHE STRING "C standard")

# Set cache variables for C++ compilers and flags
set(CMAKE_CXX_COMPILER
  "/Volumes/Workspace/spack/opt/spack/darwin-ventura-skylake/apple-clang-14.0.3/gcc-13.1.0-zexv4j6ljav6xmm3ajnobit6fdhvxwrn/bin/c++"
  CACHE FILEPATH "C++ compiler")

set(CMAKE_CXX_FLAGS
  ""
  CACHE STRING "C++ compiler flags")

set(CMAKE_CXX_STANDARD
  "14"
  CACHE STRING "C++ standard")

# Specify project name and languages
project(examples C)

# Enable testing
include(CTest)


# ------------------------------------------------------------------------------

# Specify the path to SUNDIALSConfig.cmake
set(SUNDIALS_DIR
  "/Volumes/Workspace/sundials/out/installdir/lib/cmake/sundials"
  CACHE PATH "Location of SUNDIALSConfig.cmake")

# Find SUNDIALS
find_package(SUNDIALS
  COMPONENTS  nvecmanyvector nvecserial arkode
  REQUIRED NO_DEFAULT_PATH)

# Set the SUNDIALS targets
set(SUNDIALS_TARGETS
   SUNDIALS::nvecmanyvector SUNDIALS::nvecserial SUNDIALS::arkode)

# Set any additional libraries needed
set(EXTRA_LIBS
   -lm
  CACHE STRING "Additional libraries")

# Additional includes
include_directories(. )

# ------------------------------------------------------------------------------

# Set the names of the examples to be built and their dependencies
set(examples  ark_brusselator1D_manyvec.c)
set(examples_args " ")
set(examples_dependencies )

list(LENGTH examples n_elems) # Total number of the elements in the every list
math(EXPR len "${n_elems}-1") # The last index in the every list
 message("> ${examples_args}")

# Create targets for each example

# Create targets for each example
foreach(idx RANGE ${len})

  list(GET examples ${idx} example)
  list(GET examples_args ${idx} example_args)
  message(">>>> ${example_args}")

  # get filename without extension
  get_filename_component(example_target ${example} NAME_WE)

  # create target with example source files
  if(NOT TARGET ${example_target})
    add_executable(${example_target} ${example} ${examples_dependencies})
  endif()

  # libraries to link against
  target_link_libraries(${example_target} ${SUNDIALS_TARGETS} ${EXTRA_LIBS})

  if("${example_args}" STREQUAL "")
    set(test_name ${example})
  else()
    string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  endif()

  # add the example to ctest
  add_test(NAME ${test_name} COMMAND ${example_target} ${example_args})

endforeach()
