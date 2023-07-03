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
# CMakeLists.txt file for banded sunmatrix examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS banded matrix
set(sunmatrix_band_examples
  "test_sunmatrix_band.c\; \"10 2 3 0\" \; \"\" "
  "test_sunmatrix_band.c\; \"300 7 4 0\" \; \"\" "
  "test_sunmatrix_band.c\; \"1000 8 8 0\" \; \"\" "
  "test_sunmatrix_band.c\; \"5000 3 20 0\" \; \"\" "
)

# If building F2003 tests
if (BUILD_FORTRAN_MODULE_INTERFACE)
  set(sunmatrix_band_fortran_examples
    "test_fsunmatrix_band_mod\; \"\" \; \"\" ")
endif()

# Dependencies for sunmatrix examples
set(sunmatrix_band_dependencies
  test_sunmatrix.c
  test_sunmatrix.h
  )
# Adding extra files  
set(sunmatrix_band_files
  ../test_sunmatrix.c
  ../test_sunmatrix.h
)

# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against
set(SUNMATRIX_BAND_LIB sundials_sunmatrixband)
set(NVECS_LIB sundials_nvecserial)
if (BUILD_FORTRAN_MODULE_INTERFACE)
  list(APPEND NVECS_LIB sundials_fnvecserial_mod)
  list(APPEND SUNMATRIX_BAND_LIB sundials_fsunmatrixband_mod)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${SUNMATRIX_BAND_LIB} ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${sunmatrix_band_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( sunmatrix_band_examples 
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

endforeach(example_tuple ${sunmatrix_band_examples})

# Add the build and install targets for each example
foreach(example_tuple ${sunmatrix_band_fortran_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( sunmatrix_band_fortran_examples 
    FILENAME_VAR example
    PARAMS_VAR example_args
    TYPE_VAR example_type
  )

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.f90
      ${SUNDIALS_SOURCE_DIR}/examples/utilities/test_utilities.f90
      ${SUNDIALS_SOURCE_DIR}/examples/sunmatrix/test_sunmatrix.f90)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # set fortran module directory to avoid name collisions
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

endforeach(example_tuple ${sunmatrix_band_fortran_examples})

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(sunmatrix sunmatrix_band_examples
    CMAKE_TEMPLATE
      cmakelists_C_ex.in
    MAKE_TEMPLATE
      makefile_serial_C_ex.in
    SUNDIALS_TARGETS
      nvecserial
      sunmatrixband
    DESTINATION
      sunmatrix/band
    EXAMPLES_DEPENDENCIES
      ${sunmatrix_band_dependencies}
    EXTRA_FILES
      #README
      ${sunmatrix_band_files}
    TEST_INSTALL
      band
  )

endif()