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

# include location of public and private header files

sundials_add_executable(test_arkode_error_handling
                        test_arkode_error_handling.cpp)
target_include_directories(
  test_arkode_error_handling
  PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
          ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src)

# Add variable DOMEIG_OBJECT_LIBRARIES
set(DOMEIG_OBJECT_LIBRARIES sundials_sundomeigestpi_obj)

if(ENABLE_LAPACK)
  list(APPEND DOMEIG_OBJECT_LIBRARIES sundials_sundomeigestarni_obj)
  set(LAPACK_LIBRARY LAPACK::LAPACK)
  if(CMAKE_VERSION VERSION_LESS 3.20)
    list(APPEND LAPACK_LIBRARY BLAS::BLAS)
  endif()
endif()

# We explicitly choose which object libraries to link to and link in the ARKODE
# objects so that we have access to private functions w/o changing their
# visibility in the installed libraries.
target_link_libraries(
  test_arkode_error_handling
  PRIVATE ${LAPACK_LIBRARY}
          $<TARGET_OBJECTS:sundials_arkode_obj>
          sundials_sunmemsys_obj
          sundials_nvecserial_obj
          sundials_nvecmanyvector_obj
          sundials_sunlinsolband_obj
          sundials_sunlinsoldense_obj
          sundials_sunnonlinsolnewton_obj
          sundials_sunadaptcontrollerimexgus_obj
          sundials_sunadaptcontrollersoderlind_obj
          ${DOMEIG_OBJECT_LIBRARIES}
          ${EXE_EXTRA_LINK_LIBS})

# Tell CMake that we depend on the ARKODE library since it does not pick that up
# from $<TARGET_OBJECTS:sundials_arkode_obj>.
add_dependencies(test_arkode_error_handling sundials_arkode_obj)

target_link_libraries(test_arkode_error_handling PRIVATE GTest::gtest_main
                                                         GTest::gmock)

gtest_discover_tests(test_arkode_error_handling)
