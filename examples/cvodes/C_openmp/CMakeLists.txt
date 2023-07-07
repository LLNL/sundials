# ---------------------------------------------------------------
# Programmer(s): Ting Yan @ SMU
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
# CMakeLists.txt file for CVODES OpenMP examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(CVODES_omp_examples
  "cvsAdvDiff_bnd_omp.c\; \"4\" \;develop"
  )

# Specify libraries to link against
set(CVODES_LIB sundials_cvodes)
set(NVECOMP_LIB sundials_nvecopenmp)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${CVODES_LIB} ${NVECOMP_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${CVODES_omp_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( CVODES_omp_examples
    FILENAME_VAR example
    PARAMS_VAR example_args
    TYPE_VAR example_type
  ) 

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})

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
  sundials_add_test(${example} ${example}
    TEST_ARGS ${example_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})
  
  # find all .out files for this example
  file(GLOB example_out ${example}*.out)

endforeach(example_tuple ${CVODES_omp_examples})

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(cvodes CVODES_omp_examples
    HAS_OPENMP
    CMAKE_TEMPLATE
      cmakelists_C_openmp_ex.in
    MAKE_TEMPLATE
      makefile_openmp_C_ex.in
    SUNDIALS_TARGETS
      nvecopenmp
      cvodes
    DESTINATION
      cvodes/C_openmp
    EXTRA_FILES
      README
    TEST_INSTALL
      C_openmp
  )

endif()