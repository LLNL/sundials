# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Enable SUNDIALS Testing
# ---------------------------------------------------------------

# Enable testing with 'make test'
include(CTest)


# If development tests are enabled, Python is needed to use the test runner
if(SUNDIALS_TEST_DEVTESTS)

  find_package(PythonInterp)
  if(${PYTHON_VERSION_MAJOR} LESS 3)
    if(${PYTHON_VERSION_MINOR} LESS 7)
      print_warning("Python version must be 2.7.x or greater to run development tests"
        "Examples will build but 'make test' will fail.")
    endif()
  endif()

  # Default test output directory
  set(TEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/Testing/output)

  # Create test output directory
  if(NOT EXISTS ${TEST_OUTPUT_DIR})
    file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR})
  endif()

  # If a non-default output directory was provided make sure it exists
  if(SUNDIALS_TEST_OUTPUT_DIR)
    if(NOT EXISTS ${SUNDIALS_TEST_OUTPUT_DIR})
      file(MAKE_DIRECTORY ${SUNDIALS_TEST_OUTPUT_DIR})
    endif()
  endif()

  # If a non-default answer directory was provided make sure it exists
  if(SUNDIALS_TEST_ANSWER_DIR)
    if(NOT EXISTS ${SUNDIALS_TEST_ANSWER_DIR})
      print_error("SUNDIALS_TEST_ANSWER_DIR does not exist!")
    endif()
  endif()

  # look for the testRunner script in the test directory
  find_program(TESTRUNNER testRunner PATHS test)
  hide_variable(TESTRUNNER)

endif()


# If examples are installed, create post install smoke test targets
if(EXAMPLES_INSTALL)

  # Directories for installation testing
  set(TEST_INSTALL_DIR ${PROJECT_BINARY_DIR}/Testing_Install)
  set(TEST_INSTALL_ALL_DIR ${PROJECT_BINARY_DIR}/Testing_Install_All)

  # Create installation testing directories
  if(NOT EXISTS ${TEST_INSTALL_DIR})
    file(MAKE_DIRECTORY ${TEST_INSTALL_DIR})
  endif()

  if(NOT EXISTS ${TEST_INSTALL_ALL_DIR})
    file(MAKE_DIRECTORY ${TEST_INSTALL_ALL_DIR})
  endif()

  # Create test_install and test_install_all targets
  add_custom_target(test_install
    ${CMAKE_COMMAND} -E cmake_echo_color --cyan
    "All installation tests complete.")

  add_custom_target(test_install_all
    ${CMAKE_COMMAND} -E cmake_echo_color --cyan
    "All installation tests complete.")

endif()
