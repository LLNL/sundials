# ---------------------------------------------------------------
# Programmer(s): Daniel Reynolds @ SMU
#                David J. Gardner and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for dense sunmatrix examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS dense matrix
set(sunmatrix_dense_examples
  "test_sunmatrix_dense.c\; \"100 100 0\" \; \"\" "
  "test_sunmatrix_dense.c\; \"200 1000 0\" \; \"\" "
  "test_sunmatrix_dense.c\; \"2000 100 0\" \; \"\" "
  )

# If building F2003 tests
if (BUILD_FORTRAN_MODULE_INTERFACE)
  set(sunmatrix_dense_fortran_examples
    "test_fsunmatrix_dense_mod\; \"\" \; \"\" ")
endif()

# Dependencies for sunmatrix examples
set(sunmatrix_dense_dependencies
  test_sunmatrix.c
  test_sunmatrix.h
  )
# Adding extra files  
set(sunmatrix_dense_files
  ../test_sunmatrix.c
  ../test_sunmatrix.h
)

# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against
set(SUNMATRIX_DENSE_LIB sundials_sunmatrixdense)
set(NVECS_LIB sundials_nvecserial)
if (BUILD_FORTRAN_MODULE_INTERFACE)
  list(APPEND NVECS_LIB sundials_fnvecserial_mod)
  list(APPEND SUNMATRIX_DENSE_LIB sundials_fsunmatrixdense_mod)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${SUNMATRIX_DENSE_LIB} ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${sunmatrix_dense_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( sunmatrix_dense_examples 
    FILENAME_VAR example
    PARAMS_VAR example_args
    TYPE_VAR example_type
  )

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example})

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # link vector test utilities 
    target_link_libraries(${example} PRIVATE test_sunmatrix_obj)

    # libraries to link against
    target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS})    

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

endforeach(example_tuple ${sunmatrix_dense_examples})

# Add the build and install targets for each example
foreach(example_tuple ${sunmatrix_dense_fortran_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( sunmatrix_dense_fortran_examples 
    FILENAME_VAR example
    PARAMS_VAR example_args
    TYPE_VAR example_type
  )

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    add_executable(${example} ${example}.f90
      ${SUNDIALS_SOURCE_DIR}/examples/utilities/test_utilities.f90
      ${SUNDIALS_SOURCE_DIR}/examples/sunmatrix/test_sunmatrix.f90)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")
    set_target_properties(${example} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    
    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})
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

endforeach(example_tuple ${sunmatrix_dense_fortran_examples})

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(sunmatrix sunmatrix_dense_examples
    CMAKE_TEMPLATE
      cmakelists_C_ex.in
    MAKE_TEMPLATE
      makefile_serial_C_ex.in
    SUNDIALS_TARGETS
      nvecserial
      sunmatrixdense
    DESTINATION
      sunmatrix/dense
    EXAMPLES_DEPENDENCIES
      ${sunmatrix_dense_dependencies}
    EXTRA_FILES
      #README
      ${sunmatrix_dense_files}
    TEST_INSTALL
      dense
  )

endif()
