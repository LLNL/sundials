# ---------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for CVODE RAJA examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(CVODE_examples
  "cvAdvDiff_kry_raja\;develop"
  )

# Add source directory to include directories
include_directories(.)

if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
  set(_fused_stubs sundials_cvode_fused_stubs)
endif()

if(SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
  set(_lang CUDA)
else()
  set(_lang CXX)
endif()

# Add the build and install targets for each CVODE example
foreach(example_tuple ${CVODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_type)

  set_source_files_properties(${example}.cpp PROPERTIES LANGUAGE ${_lang})

  # add example source files
  add_executable(${example} ${example}.cpp)

  set_target_properties(${example} PROPERTIES FOLDER "Examples")

  # libraries to link against
  target_link_libraries(${example}
    PRIVATE
    sundials_cvode
    sundials_nvecraja
    ${_fused_stubs}
    ${EXTRA_LINK_LIBS})

  # add example to regression tests
  sundials_add_test(${example} ${example}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.cpp ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/raja)
  endif()

endforeach(example_tuple ${CVODE_examples})


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/raja)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "CVODE")
  set(SOLVER_LIB "sundials_cvode")
  set(NVECTOR_LIB "sundials_nvecraja")
  if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
    set(LIBS "-lsundials_cvode_fused_stubs ${LIBS}")
  endif()

  examples2string(CVODE_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_RAJA_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/raja/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/cvode/raja/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/raja
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_RAJA_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/raja/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/raja/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/raja
      RENAME Makefile
      )
  endif()

  # add test_install target
  sundials_add_test_install(cvode raja)

endif()
