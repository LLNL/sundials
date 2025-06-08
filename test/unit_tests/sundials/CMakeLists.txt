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

# List of test tuples of the form "name\;args"
set(unit_tests "test_sundials_datanode\;" "test_sundials_stlvector\;"
               "test_sundials_hashmap\;")

if(SUNDIALS_ENABLE_ERROR_CHECKS)
  list(APPEND unit_tests "test_sundials_errors\;")
endif()

# Add the build and install targets for each test
if(TARGET GTest::gtest_main AND TARGET GTest::gmock)
  foreach(test_tuple ${unit_tests})
    # parse the test tuple
    list(GET test_tuple 0 test)
    list(GET test_tuple 1 test_args)

    add_executable(
      ${test} ${test}.cpp ${PROJECT_SOURCE_DIR}/src/sundials/sundials_hashmap.c)
    target_link_libraries(
      ${test} PRIVATE sundials_core sundials_nvecserial sundials_sunmemsys_obj
                      GTest::gtest_main GTest::gmock)

    gtest_discover_tests(${test})
  endforeach()
endif()


# Tests for the domeig utility routine (estimates dominant eigenvalue; depends on LAPACK)
if(ENABLE_LAPACK)
  set(unit_tests_lapack
    "test_domeig\;100 10 10\;"
    "test_domeig\;1000 20 10\;"
    "test_domeig\;10000 50 10\;")

  # Add the build and install targets for each test
  foreach(test_tuple ${unit_tests_lapack})

    # parse the test tuple
    list(GET test_tuple 0 test)
    list(GET test_tuple 1 test_args)

    # check if this test has already been added, only need to add test source
    # files once for testing with different inputs
    if(NOT TARGET ${test})

      # test source files
      add_executable(${test} ${test}.c)

      set_target_properties(${test} PROPERTIES FOLDER "unit_tests")

      # include location of public and private header files
      target_include_directories(
        ${test} PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
                        ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src)

      # libraries to link against
      set(_domeig_libs sundials_core sundials_nvecserial sundials_sundomeigestpi sundials_sundomeigestarni LAPACK::LAPACK)
      if(CMAKE_VERSION VERSION_LESS 3.20)
        list(APPEND _domeig_libs BLAS::BLAS)
      endif()
      target_link_libraries(${test} ${_domeig_libs}
        ${EXE_EXTRA_LINK_LIBS})

    endif()

    # check if test args are provided and set the test name
    if("${test_args}" STREQUAL "")
      set(test_name ${test})
    else()
      string(REPLACE " " "_" test_name "${test}_${test_args}")
      string(REPLACE " " ";" test_args "${test_args}")
    endif()

    # add test to regression tests
    add_test(NAME ${test_name} COMMAND ${test} ${test_args})

  endforeach()

  message(STATUS "Added domeig units tests")
endif()
