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

# macro for building Ginkgo examples for various backends

# Add the build targets for each CVODE example
macro(sundials_add_examples_ginkgo EXAMPLES_VAR)

  set(options )
  set(oneValueArgs )
  set(multiValueArgs BACKENDS)

  # Parse keyword arguments and options
  cmake_parse_arguments(arg
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  foreach(example ${${EXAMPLES_VAR}})
    foreach(backend ${arg_BACKENDS})

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

    endforeach()
  endforeach()

endmacro()

# macro for installing Ginkgo examples

macro(sundials_install_examples_ginkgo MODULE)

  set(options )
  set(oneValueArgs SOLVER_LIBRARY DESTINATION)
  set(multiValueArgs CPU_EXAMPLES_VAR GPU_EXAMPLES_VAR CPU_GPU_EXAMPLES_VAR
    SUNDIALS_TARGETS OTHER_TARGETS)

  # Parse keyword arguments/options
  cmake_parse_arguments(arg
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Install the example source, header, and output file
  foreach(example_type CPU GPU CPU_GPU)
    foreach(example ${${arg_${example_type}_EXAMPLES_VAR}})

      # extract the file name without extension
      get_filename_component(example_noext ${example} NAME_WE)

      # get example header and output files
      file(GLOB example_header ${example_noext}.h*)
      file(GLOB example_out ${example_noext}*.out)

      # install files
      install(FILES ${example} ${example_header} ${example_out}
        "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
        DESTINATION ${EXAMPLES_INSTALL_PATH}/${arg_DESTINATION})

    endforeach()
  endforeach()

  # Prepare substitution variables for CMakeLists and/or Makefile templates
  string(TOUPPER "${MODULE}" SOLVER)

  set(EXAMPLES_DEPENDENCIES "${arg_EXAMPLES_DEPENDENCIES}")

  if(arg_CPU_EXAMPLES_VAR)
    examples2string(${arg_CPU_EXAMPLES_VAR} CPU_EXAMPLES)
  endif()
  if(arg_GPU_EXAMPLES_VAR)
    examples2string(${arg_GPU_EXAMPLES_VAR} GPU_EXAMPLES)
  endif()
  if(arg_CPU_GPU_EXAMPLES_VAR)
    examples2string(${arg_CPU_GPU_EXAMPLES_VAR} CPU_GPU_EXAMPLES)
  endif()

  # components for find_package
  list2string(arg_SUNDIALS_TARGETS
    EXAMPLES_CMAKE_COMPONENTS)

  set(target_list "")
  foreach(target ${arg_SUNDIALS_TARGETS})
    list(APPEND target_list SUNDIALS::${target})
  endforeach()
  foreach(target ${arg_OTHER_TARGETS})
    list(APPEND target_list ${target})
  endforeach()
  list2string(target_list EXAMPLES_CMAKE_TARGETS)

  # Generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_CXX_ginkgo_ex.in
    ${PROJECT_BINARY_DIR}/examples/${arg_DESTINATION}/CMakeLists.txt
    @ONLY
  )

  # Install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/${arg_DESTINATION}/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/${arg_DESTINATION}
  )

  # Add test_install target
  if(DEFINED arg_TEST_INSTALL)
    sundials_add_test_install(${MODULE} ${arg_TEST_INSTALL})
  endif()

endmacro()

# Examples that support CPU and GPU Ginkgo backends
set(cpu_gpu_examples
  cv_heat2D_ginkgo.cpp
)

sundials_add_examples_ginkgo(cpu_gpu_examples
  BACKENDS REF OMP CUDA HIP)

# Examples that only support CPU Ginkgo backends
set(cpu_examples
  cv_kpr_ginkgo.cpp
)

sundials_add_examples_ginkgo(cpu_examples
  BACKENDS REF OMP)

# Install the targets
if(EXAMPLES_INSTALL)

  if(SUNDIALS_GINKGO_BACKENDS MATCHES "CUDA")
    list(APPEND vectors nveccuda)
  endif()
  if(SUNDIALS_GINKGO_BACKENDS MATCHES "HIP")
    list(APPEND vectors nvechip)
  endif()
  if((SUNDIALS_GINKGO_BACKENDS MATCHES "OMP") OR
      (SUNDIALS_GINKGO_BACKENDS MATCHES "REF"))
    list(APPEND vectors nvecserial)
  endif()

  sundials_install_examples_ginkgo(cvode
    CPU_EXAMPLES_VAR
      cpu_examples
    CPU_GPU_EXAMPLES_VAR
      cpu_gpu_examples
    SUNDIALS_TARGETS
      cvode
      ${vectors}
    DESTINATION
      cvode/ginkgo
  )

endif()
