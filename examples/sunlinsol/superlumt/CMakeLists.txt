# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for sunlinsol SuperLU_MT examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using the SuperLU_MT linear solver
# Note: threaded tests are excluded because of a potential bug in
# SuperLU_MT that causes random test failures
set(sunlinsol_superlumt_examples
  "test_sunlinsol_superlumt.c\;300 0 1 0\;passfail"
  "test_sunlinsol_superlumt.c\;300 1 1 0\;passfail"
  "test_sunlinsol_superlumt.c\;1000 0 3 0\;exclude"
  "test_sunlinsol_superlumt.c\;1000 1 3 0\;exclude"
)

# Dependencies for sunlinsol examples
set(sunlinsol_superlumt_dependencies
  test_sunlinsol.c
  test_sunlinsol.h
  )
# Adding extra files  
set(sunlinsol_superlumt_files
  ../test_sunlinsol.c
  ../test_sunlinsol.h
)

set(SUNLINSOL_SUPERLUMT_LIB sundials_sunlinsolsuperlumt sundials_sunmatrixdense sundials_sunmatrixsparse)
set(NVECS_LIB sundials_nvecserial)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${SUNLINSOL_SUPERLUMT_LIB} ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add source directory to include directories
include_directories(. ..)

# Add the build and install targets for each example
foreach(example_tuple ${sunlinsol_superlumt_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( sunlinsol_superlumt_examples 
    FILENAME_VAR example
    PARAMS_VAR example_args
    TYPE_VAR example_type
  )
  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example})

    # link vector test utilities 
    target_link_libraries(${example} PRIVATE test_sunlinsol_obj)

    # libraries to link against
    target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS})

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example})
  else()
    string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach(example_tuple ${sunlinsol_superlumt_examples})

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(sunlinsol sunlinsol_superlumt_examples
    CMAKE_TEMPLATE
      cmakelists_C_ex.in
    MAKE_TEMPLATE
      makefile_serial_C_ex.in
    SUNDIALS_TARGETS
      nvecserial
      sunlinsolsuperlumt
      sunmatrixdense
      sunmatrixsparse
    OTHER_TARGETS
      SUNDIALS::SUPERLUMT
    DESTINATION
      sunlinsol/superlumt
    EXAMPLES_DEPENDENCIES
      ${sunlinsol_superlumt_dependencies}
    EXTRA_FILES
      #README
      ${sunlinsol_superlumt_files}
    TEST_INSTALL
      superlumt
  )

endif()
