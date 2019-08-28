# ---------------------------------------------------------------
# Author: David J. Gardner @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
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
INCLUDE(CTest)


# If development tests are enabled, Python is needed to use the test runner
IF(SUNDIALS_DEVTESTS)

  find_package(PythonInterp)
  IF(${PYTHON_VERSION_MAJOR} LESS 3)
    IF(${PYTHON_VERSION_MINOR} LESS 7)
      PRINT_WARNING("Python version must be 2.7.x or greater to run development tests"
        "Examples will build but 'make test' will fail.")
    ENDIF()
  ENDIF()

  # Directory for test output
  SET(TEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/Testing/output)

  IF(NOT EXISTS ${TEST_OUTPUT_DIR})
    FILE(MAKE_DIRECTORY ${TEST_OUTPUT_DIR})
  ENDIF()

  # look for the testRunner script in the test directory
  FIND_PROGRAM(TESTRUNNER testRunner PATHS test)
  HIDE_VARIABLE(TESTRUNNER)

ENDIF()


# If examples are installed, create post install smoke test targets
IF(EXAMPLES_INSTALL)

  # Directories for installation testing
  SET(TEST_INSTALL_DIR ${PROJECT_BINARY_DIR}/Testing_Install)
  SET(TEST_INSTALL_ALL_DIR ${PROJECT_BINARY_DIR}/Testing_Install_All)

  # Create installation testing directories
  IF(NOT EXISTS ${TEST_INSTALL_DIR})
    FILE(MAKE_DIRECTORY ${TEST_INSTALL_DIR})
  ENDIF()

  IF(NOT EXISTS ${TEST_INSTALL_ALL_DIR})
    FILE(MAKE_DIRECTORY ${TEST_INSTALL_ALL_DIR})
  ENDIF()

  # Create test_install and test_install_all targets
  ADD_CUSTOM_TARGET(test_install
    ${CMAKE_COMMAND} -E cmake_echo_color --cyan
    "All installation tests complete.")

  ADD_CUSTOM_TARGET(test_install_all
    ${CMAKE_COMMAND} -E cmake_echo_color --cyan
    "All installation tests complete.")

ENDIF()
