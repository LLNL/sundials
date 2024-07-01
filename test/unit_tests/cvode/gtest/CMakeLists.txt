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

# include location of public and private header files

add_executable(test_cvode_error_handling test_cvode_error_handling.cpp)
target_include_directories(
  test_cvode_error_handling
  PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
          ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src)

if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
  set(_fused_link_lib sundials_cvode_fused_stubs_obj)
endif()

# We explicitly choose which object libraries to link to and link in the cvode
# objects so that we have access to private functions w/o changing their
# visibility in the installed libraries.
target_link_libraries(
  test_cvode_error_handling
  PRIVATE $<TARGET_OBJECTS:sundials_cvode_obj>
          ${_fused_link_lib}
          sundials_sunmemsys_obj
          sundials_nvecserial_obj
          sundials_sunlinsolband_obj
          sundials_sunlinsoldense_obj
          sundials_sunnonlinsolnewton_obj
          ${EXE_EXTRA_LINK_LIBS})

# Tell CMake that we depend on the CVODE library since it does not pick that up
# from $<TARGET_OBJECTS:sundials_cvode_obj>.
add_dependencies(test_cvode_error_handling sundials_cvode_obj)

target_link_libraries(test_cvode_error_handling PRIVATE GTest::gtest_main
                                                        GTest::gmock)

gtest_discover_tests(test_cvode_error_handling)
