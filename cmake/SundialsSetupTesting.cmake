# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
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

# Check if development tests are enabled
if (SUNDIALS_TEST_DEVTESTS OR BUILD_BENCHMARKS)
  # Python is needed to use the test runner
  find_package(PythonInterp REQUIRED)
  if(${PYTHON_VERSION_MAJOR} LESS 3)
    if(${PYTHON_VERSION_MINOR} LESS 7)
      print_warning("Python version must be 2.7.x or greater to run development tests"
        "Examples will build but 'make test' may fail.")
    endif()
  endif()

  # look for the testRunner script in the test directory
  find_program(TESTRUNNER testRunner PATHS test NO_DEFAULT_PATH)
  if(NOT TESTRUNNER)
    print_error("Could not locate testRunner. Set SUNDIALS_TEST_DEVTESTS=OFF or BUILD_BENCHMARKS=OFF to continue.")
  endif()
  message(STATUS "Found testRunner: ${TESTRUNNER}")
  set(TESTRUNNER ${TESTRUNNER} CACHE INTERNAL "")

endif()

# Check if development tests are enabled
if(SUNDIALS_TEST_DEVTESTS)

  message("SUNDIALS Development testing")

  # Create the default test output directory
  set(TEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/Testing/output)

  if(NOT EXISTS ${TEST_OUTPUT_DIR})
    file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR})
  endif()

  # If a non-default output directory was provided make sure it exists
  if(SUNDIALS_TEST_OUTPUT_DIR)
    message(STATUS "Using non-default test output directory: ${SUNDIALS_TEST_OUTPUT_DIR}")
    if(NOT EXISTS ${SUNDIALS_TEST_OUTPUT_DIR})
      file(MAKE_DIRECTORY ${SUNDIALS_TEST_OUTPUT_DIR})
    endif()
  endif()

  # If a non-default answer directory was provided make sure it exists
  if(SUNDIALS_TEST_ANSWER_DIR)
    message(STATUS "Using non-default test answer directory: ${SUNDIALS_TEST_ANSWER_DIR}")
    if(NOT EXISTS ${SUNDIALS_TEST_ANSWER_DIR})
      print_error("SUNDIALS_TEST_ANSWER_DIR does not exist!")
    endif()
  endif()

  # If a non-default caliper output directory was provided make sure it exists
  if(SUNDIALS_CALIPER_OUTPUT_DIR)
    message(STATUS "Using non-default caliper output directory: ${SUNDIALS_CALIPER_OUTPUT_DIR}")
    if(NOT EXISTS ${SUNDIALS_CALIPER_OUTPUT_DIR}/Example/${JOB_ID})
      file(MAKE_DIRECTORY ${SUNDIALS_CALIPER_OUTPUT_DIR}/Example/${JOB_ID})
    endif()
  endif()

  # Check if using non-default comparison precisions when testing
  if(SUNDIALS_TEST_FLOAT_PRECISION GREATER_EQUAL "0")
    message(STATUS "Using non-default float precision: ${SUNDIALS_TEST_FLOAT_PRECISION}")
  endif()

  if(SUNDIALS_TEST_INTEGER_PRECISION GREATER_EQUAL "0")
    message(STATUS "Using non-default integer precision: ${SUNDIALS_TEST_INTEGER_PRECISION}")
  endif()
  
  #
  # Target to run tests in CI containers
  # 
  if(NOT SUNDIALS_TEST_CONTAINER_EXE)
    find_program(container_exe docker)
    if(NOT container_exe)
      find_program(container_exe podman)
    endif()
    set(SUNDIALS_TEST_CONTAINER_EXE ${container_exe} CACHE PATH "Path to docker or podman" FORCE)
  endif()

  if(SUNDIALS_TEST_CONTAINER_EXE)
    add_custom_target(setup_local_ci 
      ${CMAKE_COMMAND} -E cmake_echo_color --cyan
      "Pulled SUNDIALS CI containers.")

    add_custom_target(test_local_ci
      ${CMAKE_COMMAND} -E cmake_echo_color --cyan
      "All testing with SUNDIALS CI containers complete.")

    macro(add_local_ci_target index_size precision tag)
      string(TOLOWER "${precision}" precision_)
      set(container sundials-ci-int${index_size}-${precision_})
      set(container_exe_args run ${SUNDIALS_TEST_CONTAINER_RUN_EXTRA_ARGS} -t -d --name ${container} --cap-add SYS_PTRACE 
          -v ${CMAKE_SOURCE_DIR}:${SUNDIALS_TEST_CONTAINER_MNT} ghcr.io/llnl/${container}:${tag})
      add_custom_target(setup_local_ci_${index_size}_${precision_}
        COMMENT "Pulling SUNDIALS CI container ghcr.io/llnl/${container}:${tag}"
        COMMAND ${SUNDIALS_TEST_CONTAINER_EXE} ${container_exe_args})
      add_dependencies(setup_local_ci setup_local_ci_${index_size}_${precision_})

      set(container_test_exe ./test_driver.sh)
      set(container_test_exe_args --testtype CUSTOM --env env/docker.sh --tpls --sunrealtype ${precision_} --indexsize ${index_size})
      set(container_exe_args exec -w ${SUNDIALS_TEST_CONTAINER_MNT}/test ${container} ${container_test_exe} ${container_test_exe_args})
      add_custom_target(test_local_ci_${index_size}_${precision_}
        COMMENT "Running tests in CI container ${container}:${tag}"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMAND ${SUNDIALS_TEST_CONTAINER_EXE} ${container_exe_args}
        VERBATIM)
      add_dependencies(test_local_ci test_local_ci_${index_size}_${precision_})

      unset(container)
      unset(container_exe_args)
      unset(container_test_exe)
      unset(container_test_exe_args)
    endmacro()

    add_local_ci_target(${SUNDIALS_INDEX_SIZE} ${SUNDIALS_PRECISION} latest)
  endif()

endif()

# Check if unit tests are enabled
if(SUNDIALS_TEST_UNITTESTS AND SUNDIALS_TEST_ENABLE_GTEST)
  find_package(GTest)
  if(NOT (TARGET GTest::gtest_main OR TARGET GTest::Main))
    include(FetchContent)
    FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
      GIT_TAG v1.14.0  
    )
    if(WIN32)
      # For Windows: Prevent overriding the parent project's compiler/linker settings
      set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    endif()
    FetchContent_MakeAvailable(googletest)
    include(GoogleTest)
  endif()
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

# If benchmarks are enabled, set up `make benchmark`
if(BUILD_BENCHMARKS)

  message("SUNDIALS Benchmarking")
  
  # Create benchmark targets
  add_custom_target(benchmark
    ${CMAKE_COMMAND} -E cmake_echo_color --cyan
    "All benchmarks complete."
  )
endif()
