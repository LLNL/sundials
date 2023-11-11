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

add_executable(test_sundials_errors test_sundials_errors.cpp)
target_link_libraries(test_sundials_errors PRIVATE sundials_core GTest::gtest_main GTest::gmock)
gtest_discover_tests(test_sundials_errors)

add_subdirectory(reductions)
