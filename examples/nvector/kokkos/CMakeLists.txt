# ---------------------------------------------------------------
# Programmer(s):  Cody Balos, and Daniel McGreer @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for Kokkos nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS kokkos nvector
set(nvector_kokkos_examples
  "test_nvector_kokkos\;1000 0\;\;\;"
  )

# Dependencies for nvector examples
set(nvector_examples_dependencies
  test_nvector
  sundials_nvector
  sundials_math
)

include_directories(..)

# Add the build and install targets for each example
foreach(example_tuple ${nvector_kokkos_examples})

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
    add_executable(${example} ${example}.cpp)

    # link vector test utilties
    target_link_libraries(${example} PRIVATE test_nvector_obj)
    target_link_libraries(${example} PRIVATE sundials_nveckokkos)

    # Organize examples for IDEs
    set_target_properties(${example} PROPERTIES FOLDER "Examples")
  endif()

  # check if example args are provided and set the test name
  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    if("${number_of_tasks}" STREQUAL "")
      set(test_name ${example})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks})
    endif()
  else()
    if("${number_of_tasks}" STREQUAL "")
      string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks}_${example_args})
    endif()
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_tasks}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  # install example source files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.cpp
      ../test_nvector.c
      ../test_nvector.h
      ../../../src/sundials/sundials_nvector.c
      ../../../src/sundials/sundials_math.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/kokkos)
  endif()

endforeach(example_tuple ${nvector_kokkos_examples})
