# ---------------------------------------------------------------
# Programmer(s):  Cody Balos, and Daniel McGreer @ LLNL
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
# CMakeLists.txt file for Kokkos nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the
# type is develop for examples excluded from 'make test' in releases
set(examples_list
  "test_nvector_kokkos.cpp\;1000 0\;"
)

include_directories(..)

# Add the build and install targets for each example
foreach(example_tuple ${examples_list})
  foreach(backend ${KOKKOS_EXAMPLES_BACKENDS})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_args)
    list(GET example_tuple 2 example_type)

    # extract the file name without extension
    get_filename_component(example_target ${example} NAME_WE)
    set(example_target "${example_target}.${backend}")

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    if(NOT TARGET ${example_target})
      # example source files
      add_executable(${example_target} ${example})

      # Organize examples for IDEs
      set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

      # which backend to use
      target_compile_definitions(${example_target} PRIVATE USE_${backend})

      # link vector test utilties
      target_link_libraries(${example_target} PRIVATE
        test_nvector_obj
        sundials_nveckokkos
      )
    endif()

    # check if example args are provided and set the test name
    if("${example_args}" STREQUAL "")
      set(test_name ${example_target})
    else()
      string(REGEX REPLACE " " "_" test_name ${example_target}_${example_args})
    endif()

    # add example to regression tests
    sundials_add_test(${test_name} ${example_target}
      TEST_ARGS ${example_args}
      EXAMPLE_TYPE ${example_type}
      NODIFF)

  endforeach()
endforeach()

if(EXAMPLES_INSTALL)

  sundials_install_examples(nvec_kokkos examples_list
    EXAMPLES_DEPENDENCIES
      test_nvector.c
    EXTRA_FILES
      ${SUNDIALS_SOURCE_DIR}/examples/nvector/test_nvector.c
      ${SUNDIALS_SOURCE_DIR}/examples/nvector/test_nvector.h
    CMAKE_TEMPLATE
      cmakelists_CXX_ex.in
    SUNDIALS_COMPONENTS
      nveckokkos
    SUNDIALS_TARGETS
      generic
    OTHER_TARGETS
      Kokkos::kokkos
    DESTINATION
      nvector/kokkos
    )

endif()
