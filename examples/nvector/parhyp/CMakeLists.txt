# ---------------------------------------------------------------
# Programmer(s): Steven Smith @ LLNL
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
# CMakeLists.txt file for hypre nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS parallel hypre nvector
set(nvector_parhyp_examples
  "test_nvector_parhyp.c\; \"1000 0\" \;1\;1\; \"\" " # run sequentially
  "test_nvector_parhyp.c\; \"1000 0\" \;1\;4\; \"\" " # run parallel on 4 procs
  )

# Dependencies for nvector parallel examples
set(nvector_examples_dependencies
  test_nvector.c
  test_nvector.h
  test_mpinvector.c
  )

# Adding extra files 
set(nvector_dependencies
  ../test_nvector.c
  ../test_nvector.h
  ../test_mpinvector.c
)

# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against
set(NVECS_LIB sundials_nvecparhyp)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Add the build and install targets for each example
foreach(example_tuple ${nvector_parhyp_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( nvector_parhyp_examples 
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
    add_executable(${example} ${example}  )

    # link vector test utilties
    target_link_libraries(${example} PRIVATE test_nvector_obj test_nvectormpi_obj)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS})

    if(NOT MPI_C_COMPILER)
      target_link_libraries(${example} PRIVATE ${MPI_LIBRARIES})
    endif()
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
    MPI_NPROCS ${number_of_tasks}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach(example_tuple ${nvector_parhyp_examples})

if(EXAMPLES_INSTALL)

  sundials_install_examples(nvector nvector_parhyp_examples
    HAS_MPI
    CMAKE_TEMPLATE
      cmakelists_C_MPI_ex.in
    MAKE_TEMPLATE
      makefile_parallel_C_ex.in
    SUNDIALS_TARGETS
      nvecparhyp
    OTHER_TARGETS
      SUNDIALS::HYPRE
    DESTINATION
      nvector/parhyp
    EXAMPLES_DEPENDENCIES
      ${nvector_examples_dependencies}
    EXTRA_FILES
      ${nvector_dependencies}
    TEST_INSTALL
      parhyp
  )

endif()