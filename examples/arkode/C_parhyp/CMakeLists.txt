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
# CMakeLists.txt file for ARKODE parhyp C examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(ARKODE_examples
  "ark_diurnal_kry_ph.c\;<none>\;1\;4\;develop"
  )

if(MPI_C_COMPILER)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Specify libraries to link against
set(ARKODE_LIB sundials_arkode)
set(NVECPARHYP_LIB sundials_nvecparhyp)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${ARKODE_LIB} ${NVECPARHYP_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each ARKODE example
foreach(example_tuple ${ARKODE_examples})

 # parse the example tuple using sundials_parse_tuple macro
 sundials_parse_tuple( ARKODE_examples 
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

  # add to regression tests
  sundials_add_test(${test_name} ${example}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

endforeach(example_tuple ${ARKODE_examples})


if(EXAMPLES_INSTALL)

  sundials_install_examples(arkode ARKODE_examples
    HAS_MPI
    CMAKE_TEMPLATE
      cmakelists_C_MPI_ex.in
    MAKE_TEMPLATE
      makefile_parallel_C_ex.in
    SOLVER_LIBRARY
      sundials_arkode
    SUNDIALS_TARGETS
      arkode
      nvecparhyp
    OTHER_TARGETS
      SUNDIALS::HYPRE
    DESTINATION
      arkode/C_parhyp
    EXTRA_FILES
      README
    TEST_INSTALL
      C_parhyp
  )

endif()
