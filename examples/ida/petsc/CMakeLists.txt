# ---------------------------------------------------------------
# Programmer(s): Slaven Peles, Eddy Banks and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for IDA PETSc examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args;\nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(IDA_examples
  "idaHeat2D_petsc_spgmr.c\;<none>\;1\;4\;develop"
  "idaHeat2D_petsc_snes.c\;<none>\;1\;4\;develop"
  "idaHeat2D_petsc_snes.c\;-pre\;1\;4\;develop"
  "idaHeat2D_petsc_snes.c\;-jac\;1\;4\;develop"
  "idaHeat2D_petsc_snes.c\;-jac -pre\;1\;4\;develop"
  )

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Specify libraries to link against
set(IDA_LIB sundials_ida)
set(NVECP_LIB sundials_nvecpetsc)
set(SUNNLS_LIB sundials_sunnonlinsolpetscsnes)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${IDA_LIB} ${NVECP_LIB} ${SUNNLS_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${IDA_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( IDA_examples 
    FILENAME_VAR example
    PARAMS_VAR example_args
    NODES_VAR number_of_nodes
    NPROCS_VAR number_of_tasks
    TYPE_VAR example_type
  ) 

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    add_executable(${example} ${example})

    # group examples in IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})

    if(NOT MPI_C_COMPILER)
      target_link_libraries(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
    endif()
  endif()

  # check if example args are provided and set the test name
  sundials_create_test_name(${example} test_name MPI_TEST_NAME
    TEST_ARGS ${example_args} MPI_NPROCS ${number_of_tasks})

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

  # find all .out files for this example
  file(GLOB example_out ${example}*.out)

endforeach(example_tuple ${IDA_examples})

if(EXAMPLES_INSTALL)

  sundials_install_examples(ida IDA_examples
    HAS_MPI
    CMAKE_TEMPLATE
      cmakelists_C_MPI_ex.in
    MAKE_TEMPLATE
      makefile_parallel_C_ex.in
    SUNDIALS_TARGETS
      ida
      nvecpetsc
      sunnonlinsolpetscsnes
    OTHER_TARGETS
      SUNDIALS::PETSC
    DESTINATION
      ida/petsc
    EXTRA_FILES
      README
    TEST_INSTALL
      petsc
  )

endif()
