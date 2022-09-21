# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# The macro:
#
#   SUNDIALS_INSTALL_EXAMPLES(<MODULE> <EXAMPLES_VAR>
#     DESTINATION path
#     CMAKE_TEMPLATE name
#     [MAKE_TEMPLATE name [SOLVER_LIBRARY target]]
#     [TEST_INSTALL target]
#     [SUNDIALS_TARGETS targets]
#     [EXTRA_FILES files]
#     [EXTRA_INCLUDES includes]
#   )
#
# adds an install target for examples in EXAMPLES_VAR that go with MODULE (e.g.
# arkode, nvecserial).
#
# The DESTINATION option is the path *within* EXAMPLES_INSTALL_PATH
# that the files should be installed.
#
# The CMAKE_TEMPLATE option is the name of the examples/templates
# CMake template to use (e.g. cmakelists_CXX_ex.in)
#
# The MAKE_TEMPLATE option is the name of the examples/templates Make
# template to use
#
# The SOLVER_LIBRARY option is used when a MAKE_TEMPLATE is
# provided. It should be the library name for SUNDIALS solver
# (e.g. arkode, cvode, ...)
#
# The TEST_INSTALL option adds a test_install target with the given
# target name for the MODULE.
#
# The SUNDIALS_TARGETS option is a list of CMake targets in the
# SUNDIALS:: namespace that the examples need to be linked to.
#
# The OTHER_TARGETS option is a list of CMake targets that the
# examples need to be linked to.
#
# The EXAMPLES_DEPENDENCIES option is a list of additional source
# files that the examples are dependent on.
#
# The EXTRA_FILES option is a list of files to install that are not
# example source code.
#
# The EXTRA_INCLUDES option is a list of additional includes to set
# with INCLUDE_DIRECTORIES.
# ------------------------------------------------------------------------------

# Add the build targets for each CVODE example
macro(sundials_add_examples_ginkgo EXAMPLES_VAR)

  set(options UNIT_TEST)
  set(oneValueArgs TARGETS)
  set(multiValueArgs BACKENDS)

  # Parse keyword arguments and options
  cmake_parse_arguments(arg
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  foreach(example_tuple ${${EXAMPLES_VAR}})
    foreach(backend ${arg_BACKENDS})

      # parse the example tuple
      list(GET example_tuple 0 example)
      list(GET example_tuple 1 example_args)
      list(GET example_tuple 2 example_type)

      if(NOT (SUNDIALS_GINKGO_BACKENDS MATCHES "${backend}"))
        continue()
      endif()

      if(backend MATCHES "CUDA")
        set_source_files_properties(${example} PROPERTIES LANGUAGE CUDA)
        set(vector nveccuda)
      elseif(backend MATCHES "HIP")
        set_source_files_properties(${example} PROPERTIES LANGUAGE CXX)
        set(vector nvechip)
      elseif(backend MATCHES "OMP")
        set(vector nvecserial)
      elseif(backend MATCHES "REF")
        set(vector nvecserial)
      endif()

      # extract the file name without extension
      get_filename_component(example_target ${example} NAME_WE)
      set(example_target "${example_target}.${backend}")

      if(NOT TARGET ${example_target})

        # create target
        add_executable(${example_target} ${example})

        # folder for IDEs
        set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

        # which backend to use
        target_compile_definitions(${example_target} PRIVATE USE_${backend})

        # directories to include
        target_include_directories(${example_target}
          PRIVATE
          "${PROJECT_SOURCE_DIR}/examples/utilities")

        # libraries to link against
        target_link_libraries(${example_target}
          PRIVATE
          ${arg_TARGETS}
          sundials_${vector}
          Ginkgo::ginkgo
          ${EXTRA_LINK_LIBS})

      endif()

      # check if example args are provided and set the test name
      if("${example_args}" STREQUAL "")
        set(test_name ${example_target})
      else()
        string(REGEX REPLACE " " "_" test_name ${example_target}_${example_args})
      endif()

      # add example to regression tests
      if(${arg_UNIT_TEST})
        sundials_add_test(${test_name} ${example_target}
          TEST_ARGS ${example_args}
          EXAMPLE_TYPE ${example_type}
          NODIFF)
      else()
        sundials_add_test(${test_name} ${example_target}
          TEST_ARGS ${example_args}
          ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
          ANSWER_FILE ${test_name}.out
          EXAMPLE_TYPE ${example_type})
      endif()

    endforeach()
  endforeach()

endmacro()
