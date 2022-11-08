# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------

if(ENABLE_CALIPER AND SUNDIALS_BUILD_WITH_PROFILING)
  set(_link_caliper_if_needed PRIVATE caliper)
endif()

add_executable(test_sundials_errors test_sundials_errors.cpp)
target_link_libraries(test_sundials_errors PRIVATE sundials_generic ${_link_caliper_if_needed} GTest::gtest_main)
gtest_discover_tests(test_sundials_errors)
