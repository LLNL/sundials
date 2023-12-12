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

# include location of public and private header files

add_executable(test_arkode_error_handling test_arkode_error_handling.cpp)
target_include_directories(test_arkode_error_handling PRIVATE
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include> ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_arkode_error_handling
    PRIVATE  sundials_arkode sundials_nvecserial sundials_core GTest::gtest_main GTest::gmock)
gtest_discover_tests(test_arkode_error_handling)
