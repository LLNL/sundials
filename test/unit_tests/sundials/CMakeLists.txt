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

# List of test tuples of the form "name\;args"
set(unit_tests
  "test_sundials_errors\;"
  "test_sundials_datanode\;"
  "test_sundials_arraylist\;"
  "test_sundials_hashmap\;")

# Add the build and install targets for each test
if(TARGET GTest::gtest_main AND TARGET GTest::gmock)
  foreach(test_tuple ${unit_tests})
    # parse the test tuple
    list(GET test_tuple 0 test)
    list(GET test_tuple 1 test_args)

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test}
      PRIVATE sundials_core sundials_nvecserial GTest::gtest_main GTest::gmock)

    gtest_discover_tests(${test})
  endforeach()
endif()

