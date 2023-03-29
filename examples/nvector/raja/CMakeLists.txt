# ---------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody Balos @ LLNL
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
# CMakeLists.txt file for raja nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the
# type is develop for examples excluded from 'make test' in releases
set(examples_list
  "test_nvector_raja.cpp\;1000 0\;"
  )

# Add source directory to include directories
include_directories(. ..
  ${SUNDIALS_SOURCE_DIR}/examples/utilities)

if(SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
  set(_lang CUDA)
else()
  set(_lang CXX)
endif()

if(SUNDIALS_RAJA_BACKENDS MATCHES "SYCL")
  set(_custom_memhelper custom_memory_helper_sycl.h)
else()
  set(_custom_memhelper custom_memory_helper_gpu.h)
endif()

if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  set(OTHER_LIBS OpenMP::OpenMP_CXX)
endif()

# Add the build and install targets for each example
foreach(example_tuple ${examples_list})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example_target})

    # add example source files
    add_executable(${example_target} ${example})

    set_source_files_properties(${example_target} PROPERTIES LANGUAGE ${_lang})

    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example_target}
      PRIVATE
      test_nvector_obj
      sundials_nvecraja
      ${EXE_EXTRA_LINK_LIBS}
      ${OTHER_LIBS})

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

if(EXAMPLES_INSTALL)

  if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
    set(EXAMPLES_FIND_PACKAGE "find_package(OpenMP REQUIRED)\n")
    set(_openmp_target OpenMP::OpenMP_CXX)
  endif()

  if(RAJA_NEEDS_THREADS)
    set(EXAMPLES_FIND_PACKAGE "${EXAMPLES_FIND_PACKAGE}find_package(Threads REQUIRED)\n")
  endif()

  sundials_install_examples(nvec_raja examples_list
    EXAMPLES_DEPENDENCIES
      test_nvector.c
    EXTRA_FILES
      ${SUNDIALS_SOURCE_DIR}/examples/nvector/test_nvector.c
      ${SUNDIALS_SOURCE_DIR}/examples/nvector/test_nvector.h
      ${SUNDIALS_SOURCE_DIR}/examples/utilities/${_custom_memhelper}
    CMAKE_TEMPLATE
      cmakelists_${_lang}_ex.in
    SUNDIALS_TARGETS
      generic
      nvecraja
    OTHER_TARGETS
      ${_openmp_target}
    DESTINATION
      nvector/raja
    )

endif()
