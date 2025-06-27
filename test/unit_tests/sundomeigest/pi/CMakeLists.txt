# ---------------------------------------------------------------
# Programmer(s): Mustafa Aggul @ SMU
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for sundomeigest PI examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases

# Examples using SUNDIALS PI dominant eigenvalue estimator
set(sundomeigest_pi_examples
    "test_sundomeigest_pi\;100 100 0 0\;"
    "test_sundomeigest_pi\;1000 100 0 0\;"
    "test_sundomeigest_pi\;10000 100 0 0\;"
    "test_sundomeigest_pi\;100000 100 0 0\;")

# Dependencies for dominant eigenvalue examples
set(sundomeigest_pi_dependencies test_sundomeigest)

# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against
set(DOMEIG_LIB sundials_sundomeigestpi)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS sundials_nvecserial ${DOMEIG_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${sundomeigest_pi_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add example
  # source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    sundials_add_executable(${example} ${example}.c)

    # link domeigest test utilities
    target_link_libraries(${example} PRIVATE test_sundomeigest_obj)

    # libraries to link against
    target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS})

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example})
  else()
    string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(
    ${test_name} ${example}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  # install example source files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.c ../test_sundomeigest.h ../test_sundomeigest.c
            DESTINATION ${EXAMPLES_INSTALL_PATH}/sundomeigest/pi)
  endif()

endforeach(example_tuple ${sundomeigest_pi_examples})

if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/sundomeigest/pi)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER_LIB "sundials_sundomeigestpi")

  examples2string(sundomeigest_pi_examples EXAMPLES)
  set(sundomeigest_pi_dependencies test_sundomeigest.c)
  examples2string(sundomeigest_pi_dependencies EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then be used
  # as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
    ${PROJECT_BINARY_DIR}/test/unit_tests/sundomeigest/pi/CMakeLists.txt @ONLY)

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/test/unit_tests/sundomeigest/pi/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/sundomeigest/pi)

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template for the
  # user's own programs.
  set(EXAMPLES_DEPENDENCIES)
  set(sundomeigest_pi_dependencies test_sundomeigest)
  examples2string(sundomeigest_pi_dependencies EXAMPLES_DEPENDENCIES)

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/test/unit_tests/sundomeigest/pi/Makefile_ex @ONLY)
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/test/unit_tests/sundomeigest/pi/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sundomeigest/pi
      RENAME Makefile)
  endif()

endif()
