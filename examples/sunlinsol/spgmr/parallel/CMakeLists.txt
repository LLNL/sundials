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
# CMakeLists.txt file for sunlinsol SPGMR parallel examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using the SUNDIALS SPGMR linear solver
set(sunlinsol_spgmr_examples
  "test_sunlinsol_spgmr_parallel.c\;100 1 1 50 1e-3 0\;1\;4\;passfail"
  "test_sunlinsol_spgmr_parallel.c\;100 1 2 50 1e-3 0\;1\;4\;passfail"
  "test_sunlinsol_spgmr_parallel.c\;100 2 1 50 1e-3 0\;1\;4\;passfail"
  "test_sunlinsol_spgmr_parallel.c\;100 2 2 50 1e-3 0\;1\;4\;passfail"
  )

# Dependencies for sunlinsol examples
set(sunlinsol_spgmr_dependencies
  test_sunlinsol.c
  test_sunlinsol.h
  )
# Adding extra files  
set(sunlinsol_spgmr_files
  ../../test_sunlinsol.c
  ../../test_sunlinsol.h
  )

# Add source directory to include directories
include_directories(. ../..)

# Specify libraries to link against
set(SUNLINSOL_SPGMR_LIB sundials_sunlinsolspgmr)
set(NVECS_LIB sundials_nvecparallel)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${SUNLINSOL_SPGMR_LIB} ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()


# Add the build and install targets for each example
foreach(example_tuple ${sunlinsol_spgmr_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( sunlinsol_spgmr_examples
    FILENAME_VAR example
    PARAMS_VAR example_args
    NODES_VAR number_of_nodes
    NPROCS_VAR number_of_tasks
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

    if(NOT MPI_C_COMPILER)
      target_link_libraries(${example} ${MPI_LIBRARIES})
    endif()
  endif()

  # check if example args are provided and set the test name
  sundials_create_test_name(${example} test_name MPI_TEST_NAME
    TEST_ARGS ${example_args} MPI_NPROCS ${number_of_tasks})

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_tasks}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach(example_tuple ${sunlinsol_spgmr_examples})

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(sunlinsol sunlinsol_spgmr_examples
    HAS_MPI
    CMAKE_TEMPLATE
      cmakelists_C_MPI_ex.in
    MAKE_TEMPLATE
      makefile_parallel_C_ex.in
    SUNDIALS_TARGETS
      nvecparallel
      sunlinsolspgmr
    DESTINATION
      sunlinsol/spgmr_parallel
    EXAMPLES_DEPENDENCIES
      ${sunlinsol_spgmr_dependencies}
    EXTRA_FILES
      #README
      ${sunlinsol_spgmr_files}
    TEST_INSTALL
      spgmr_parallel
  )

endif()
