# ------------------------------------------------------------------------------
# Programmer(s): Cody Balos, Shahbaj Sohal @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------

macro(sundials_create_test_name example test_name)

  # macro options
  set(options )

  # macro keyword inputs followed by a single value
  set(oneValueArgs )

  # macro keyword inputs followed by multiple values
  # TEST_ARGS = command line arguments to pass to the test executable
  set(multiValueArgs "TEST_ARGS")

  # parse inputs and create variables SUNDIALS_ADD_TEST_<keyword>
  cmake_parse_arguments(sundials_create_test_name
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
  get_filename_component(file_wo_ext ${example} NAME_WE)
  if("${sundials_create_test_name_TEST_ARGS}" STREQUAL "<none>")
    set(${test_name} ${file_wo_ext})
  else()
    string(REGEX REPLACE " " "_" ${test_name} ${file_wo_ext}_${sundials_create_test_name_TEST_ARGS})
  endif()

endmacro()