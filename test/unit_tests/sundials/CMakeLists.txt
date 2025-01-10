# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2024, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------

if(SUNDIALS_ENABLE_ERROR_CHECKS)
  if(TARGET GTest::gtest_main AND TARGET GTest::gmock)
    add_executable(test_sundials_errors test_sundials_errors.cpp)
    target_link_libraries(
      test_sundials_errors PRIVATE sundials_core sundials_nvecserial
                                   GTest::gtest_main GTest::gmock)
    gtest_discover_tests(test_sundials_errors)
  endif()
endif()

# Tests for the arnoldi utility routine (estimates dominant eigenvalue)
if(ENABLE_LAPACK)
  set(unit_tests
    "test_arnoldi\;100 10\;"
    "test_arnoldi\;1000 20\;"
    "test_arnoldi\;10000 50\;")

  # Add the build and install targets for each test
  foreach(test_tuple ${unit_tests})

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
      set(_arnoldi_libs sundials_core sundials_nvecserial LAPACK::LAPACK)
      if(CMAKE_VERSION VERSION_LESS 3.20)
        list(APPEND _arnoldi_libs BLAS::BLAS)
      endif()
      target_link_libraries(${test} ${_arnoldi_libs}
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

  message(STATUS "Added arnoldi units tests")
endif()


add_subdirectory(reductions)
