# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

set(CVODE_examples
  "cv_kpr_ginkgo.cpp\;\;"
)

# if(SUNDIALS_GINKGO_BACKENDS MATCHES "CUDA")
#   set(vector nveccuda)
#   set(cuda_or_hip CUDA)
# elseif(SUNDIALS_GINKGO_BACKENDS MATCHES "HIP")
#   set(vector nvechip)
#   set(cuda_or_hip HIP)
# endif()

set(vector nvecserial)
set(cmake_template CXX)

# Add the build targets for each CVODE example
foreach(example_tuple ${CVODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  if (NOT TARGET ${example_target})
    # example source files
    add_executable(${example_target} ${example})

    # folder for IDEs
    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

    # directories to include
    target_include_directories(${example_target}
      PRIVATE
      ${PROJECT_SOURCE_DIR}/examples/utilities
    )

    # libraries to link against
    target_link_libraries(${example_target}
      PRIVATE
      sundials_cvode
      sundials_${vector}
      Ginkgo::ginkgo
      ${EXTRA_LINK_LIBS}
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
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

endforeach(example_tuple ${CVODE_examples})

# Install the targets
if(EXAMPLES_INSTALL AND CVODE_examples)

  sundials_install_examples(cvode CVODE_examples
    CMAKE_TEMPLATE
      cmakelists_${cmake_template}_ex.in
    SUNDIALS_TARGETS
      cvode
      ${vector}
      sunlinsolGINKGOdense
    DESTINATION
      cvode/ginkgo
  )

endif()
