# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
# ARKode C++ serial unit_tests
# ---------------------------------------------------------------

# List of test tuples of the form "name\;args"
set(unit_tests
  "ark_test_analytic_sys_mri.cpp\;0"
  "ark_test_analytic_sys_mri.cpp\;1"
  "ark_test_dahlquist_ark.cpp\;0 0 0"
  "ark_test_dahlquist_ark.cpp\;0 0 1"
  "ark_test_dahlquist_ark.cpp\;0 1 0"
  "ark_test_dahlquist_ark.cpp\;0 1 1"
  "ark_test_dahlquist_ark.cpp\;1 0 0"
  "ark_test_dahlquist_ark.cpp\;1 0 1"
  "ark_test_dahlquist_ark.cpp\;1 1 0"
  "ark_test_dahlquist_ark.cpp\;1 1 1"
  "ark_test_dahlquist_ark.cpp\;2 0 0"
  "ark_test_dahlquist_ark.cpp\;2 0 1"
  "ark_test_dahlquist_ark.cpp\;2 1 0"
  "ark_test_dahlquist_ark.cpp\;2 1 1"
  "ark_test_dahlquist_erk.cpp\;0"
  "ark_test_dahlquist_erk.cpp\;1"
  "ark_test_dahlquist_mri.cpp\;"
  "ark_test_butcher.cpp\;"
  "ark_test_getjac.cpp\;"
  "ark_test_getjac_mri.cpp\;"
)

# Add the build and install targets for each test
foreach(test_tuple ${unit_tests})

  # Parse the test tuple
  list(GET test_tuple 0 test)
  list(GET test_tuple 1 test_args)

  # Extract the file name without extension
  get_filename_component(test_target ${test} NAME_WE)

  # Check if this test has already been added, only need to add test source
  # files once for testing with different inputs
  if(NOT TARGET ${test_target})

    # Test source files
    add_executable(${test_target} ${test})

    # Folder to organize targets in an IDE
    set_target_properties(${test_target} PROPERTIES FOLDER "unit_tests")

    # Include location of public and private header files
    target_include_directories(${test_target} PRIVATE
      $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
      ${CMAKE_SOURCE_DIR}/include
      ${CMAKE_SOURCE_DIR}/src)

    # We explicitly choose which object libraries to link to and link in the
    # arkode objects so that we have access to private functions w/o changing
    # their visibility in the installed libraries.
    target_link_libraries(${test_target}
      $<TARGET_OBJECTS:sundials_arkode_obj>
      sundials_sunmemsys_obj
      sundials_nvecserial_obj
      sundials_sunlinsolband_obj
      sundials_sunlinsoldense_obj
      sundials_sunnonlinsolnewton_obj
      sundials_sunnonlinsolfixedpoint_obj
      sundials_sunadaptcontrollerimexgus_obj
      sundials_sunadaptcontrollersoderlind_obj
      ${EXE_EXTRA_LINK_LIBS})

    # Tell CMake that we depend on the ARKODE library since it does not pick
    # that up from $<TARGET_OBJECTS:sundials_arkode_obj>.
    add_dependencies(${test_target} sundials_arkode_obj)

  endif()

  # Check if test args are provided and set the test name
  if("${test_args}" STREQUAL "")
    set(test_name ${test_target})
  else()
    string(REPLACE " " "_" test_name "${test_target}_${test_args}")
  endif()

  if(SUNDIALS_PRECISION MATCHES "DOUBLE")
    set(diff_output "")
  else()
    set(diff_output "NODIFF")
  endif()

  # add test to regression tests
  sundials_add_test(${test_name} ${test_target}
    TEST_ARGS ${test_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE "develop"
    ${diff_output}
  )

endforeach()

message(STATUS "Added ARKODE CXX serial unit tests")
