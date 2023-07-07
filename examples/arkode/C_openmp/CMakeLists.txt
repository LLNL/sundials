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
# CMakeLists.txt file for ARKODE OpenMP examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(ARKODE_examples
  "ark_brusselator1D_omp.c\; \"4\" \;exclude-single\;default\;default"
  "ark_heat1D_omp.c\; \"4\" \;develop\;default\;10"
  )

# Auxiliary files to install
set(ARKODE_extras_OMP
  plot_brusselator1D.py
  )

# Specify libraries to link against
set(ARKODE_LIB sundials_arkode)
set(NVECP_LIB sundials_nvecopenmp)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${ARKODE_LIB} ${NVECP_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each ARKODE example
foreach(example_tuple ${ARKODE_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( ARKODE_examples 
  FILENAME_VAR example
  PARAMS_VAR example_args
  TYPE_VAR example_type
  FLT_PRECISION_VAR example_float_precision
  INT_PRECISION_VAR example_int_precision
  ) 

  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})

    set_target_properties(${example} PROPERTIES FOLDER "Examples")
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example})
  else()
    string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  endif()

  # add to regression tests
  sundials_add_test(${example} ${example}
    TEST_ARGS ${example_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type}
    FLOAT_PRECISION ${example_float_precision}
    INTEGER_PRECISION ${example_int_precision})
  
  # find all .out files for this example
  file(GLOB example_out ${example}*.out)

endforeach(example_tuple ${ARKODE_examples})

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(arkode ARKODE_examples
    HAS_OPENMP
    CMAKE_TEMPLATE
      cmakelists_C_openmp_ex.in
    MAKE_TEMPLATE
      makefile_openmp_C_ex.in
    SUNDIALS_TARGETS
      nvecopenmp
      arkode
    DESTINATION
      arkode/C_openmp
    EXTRA_FILES
      README
      ${ARKODE_extras_OMP}
      # ${ARKODE_headers}
    TEST_INSTALL
      C_openmp
  )

endif()