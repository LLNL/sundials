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
# CMakeLists.txt file for sunlinsol LAPACK band examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using the LAPACK band linear solver
set(sunlinsol_lapackband_examples
  "test_sunlinsol_lapackband.c\;10 2 3 0 0\;passfail"
  "test_sunlinsol_lapackband.c\;300 7 4 0 0\;passfail"
  "test_sunlinsol_lapackband.c\;1000 8 8 0 0\;passfail"
  "test_sunlinsol_lapackband.c\;5000 3 100 0 0\;passfail"
  )

# Dependencies for nvector examples
set(sunlinsol_lapackband_dependencies
  test_sunlinsol.c
  test_sunlinsol.h
  )
# Adding extra files  
set(sunlinsol_lapackband_files
  ../test_sunlinsol.c
  ../test_sunlinsol.h
)

# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against
set(SUNLINSOL_BAND_LIB sundials_sunlinsollapackband sundials_sunmatrixband)
set(NVECS_LIB sundials_nvecserial)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${SUNLINSOL_BAND_LIB} ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each sunlinsol example
foreach(example_tuple ${sunlinsol_lapackband_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( sunlinsol_lapackband_examples 
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

    # libraries to link against
    target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS})

    # link vector test utilities 
    target_link_libraries(${example} PRIVATE test_sunlinsol_obj)

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

endforeach(example_tuple ${sunlinsol_lapackband_examples})

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(sunlinsol sunlinsol_lapackband_examples
    CMAKE_TEMPLATE
      cmakelists_C_ex.in
    MAKE_TEMPLATE
      makefile_serial_C_ex.in
    SUNDIALS_TARGETS
      nvecserial
      sunlinsollapackband
      sunmatrixband
    DESTINATION
      sunlinsol/lapackband
    EXAMPLES_DEPENDENCIES
      ${sunlinsol_lapackband_dependencies}
    EXTRA_FILES
      #README
      ${sunlinsol_lapackband_files}
    TEST_INSTALL
      lapackband
  )

endif()
