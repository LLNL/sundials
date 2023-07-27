# ---------------------------------------------------------------
# Programmer:  Cody J. Balos @ LLNL
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
# CMakeLists.txt file for CVODE PETSc examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(CVODE_examples
  "cv_petsc_ex7.c\;<none>\;1\;1\;develop"
  "cvAdvDiff_petsc.c\;<none>\;1\;4\;develop"
  )

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Specify libraries to link against
set(CVODE_LIB sundials_cvode)
set(NVECP_LIB sundials_nvecpetsc)
set(SUNNLS_LIB sundials_sunnonlinsolpetscsnes)

if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
  list(APPEND CVODE_LIB sundials_cvode_fused_stubs)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${CVODE_LIB} ${NVECP_LIB} ${SUNNLS_LIB})

# Add the build and install targets for each example
foreach(example_tuple ${CVODE_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( CVODE_examples 
    FILENAME_VAR example
    NODES_VAR number_of_nodes
    NPROCS_VAR number_of_tasks
    TYPE_VAR example_type
  ) 

  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})

    if(NOT MPI_C_COMPILER)
      target_link_libraries(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
    endif()

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
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

endforeach(example_tuple ${CVODE_examples})

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
    set(_include_cvode_fused_stubs cvode_fused_stubs)
  endif()

  sundials_install_examples(cvode CVODE_examples
    HAS_MPI
    CMAKE_TEMPLATE
      cmakelists_C_MPI_ex.in
    MAKE_TEMPLATE
      makefile_parallel_C_ex.in
    SUNDIALS_TARGETS
      cvode
      nvecpetsc
      sunnonlinsolpetscsnes
      ${_include_cvode_fused_stubs}
    OTHER_TARGETS
      SUNDIALS::PETSC
    DESTINATION
      cvode/petsc
    EXTRA_FILES
      README
    TEST_INSTALL
      petsc
  )

endif()