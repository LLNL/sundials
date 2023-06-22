# ---------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
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
# CMakeLists.txt file for IDAS parallel examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(IDAS_examples
  "idasBruss_ASAp_kry_bbd_p.c\;--\;1\;4\;exclude-single"
  "idasBruss_FSA_kry_bbd_p.c\;--\;1\;4\;exclude-single"
  "idasBruss_kry_bbd_p.c\;--\;1\;4\;exclude-single"
  "idasFoodWeb_kry_bbd_p.c\;--\;1\;4\;exclude-single"
  "idasFoodWeb_kry_p.c\;--\;1\;4\;exclude-single"
  "idasHeat2D_FSA_kry_bbd_p.c\;-sensi stg t\;1\;4\;exclude-single"
  "idasHeat2D_kry_bbd_p.c\;--\;1\;4\;exclude-single"
  "idasHeat2D_kry_p.c\;--\;1\;4\;exclude-single"
  )

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Specify libraries to link against
set(IDAS_LIB sundials_idas)
set(NVECP_LIB sundials_nvecparallel)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${IDAS_LIB} ${NVECP_LIB} ${EXE_EXTRA_LINK_LIBS})


# Add the build and install targets for each example
foreach(example_tuple ${IDAS_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 number_of_tasks)
  list(GET example_tuple 4 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example})

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})

    if(NOT MPI_C_COMPILER)
      target_link_libraries(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
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
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

  # find all .out files for this example
  file(GLOB example_out ${example}*.out)

endforeach(example_tuple ${IDAS_examples})

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(idas IDAS_examples
    CMAKE_TEMPLATE
      cmakelists_C_MPI_ex.in
    MAKE_TEMPLATE
      makefile_parallel_C_ex.in
    SUNDIALS_TARGETS
      nvecparallel
      idas
    DESTINATION
      idas/parallel
    EXTRA_FILES
      README
      ${IDAS_extras}
      ${IDAS_headers}
    TEST_INSTALL
      parallel
  )

endif()