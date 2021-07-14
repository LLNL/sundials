# ---------------------------------------------------------------
# Programmer(s): Shelby Lockhart @ LLNL
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
# -----------------------------------------------------------------
# CMakeLists.txt file for KINSOL C++ parallel examples
# -----------------------------------------------------------------

# Example lists are tuples "name\;compile defs\;args\;nodes\;tasks\;type\;float precision"
# where the type is develop for examples excluded from 'make test' in releases.

# Examples to run.
set(serial_examples
  "kin_heat2D_nonlin_p\;\;--np 2 2\;1\;4\;develop\;default")

# Header files to install
set(KINSOL_headers "kin_heat2D_nonlin_p.hpp")

if(MPI_CXX_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  # disable C++ extensions (for known wrappers)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX -DLAM_BUILDING")
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# macro to loop through each set of examples and build/install them
macro(build_examples examples_to_build lang)

  # Add the build and install targets for each example
  foreach(example_tuple ${${examples_to_build}})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_defines)
    list(GET example_tuple 2 example_args)
    list(GET example_tuple 3 number_of_nodes)
    list(GET example_tuple 4 number_of_tasks)
    list(GET example_tuple 5 example_type)
    list(GET example_tuple 6 example_precision)

    if("${example_defines}" STREQUAL "")
      set(example_target ${example})
    else()
      string(REGEX REPLACE " " "." example_target ${example}.${example_defines})
    endif()

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    if(NOT TARGET ${example_target})
      set_source_files_properties(${example}.cpp PROPERTIES LANGUAGE ${lang})

      # example source files
      add_executable(${example_target} ${example}.cpp)

      # folder to organize targets in an IDE
      set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

      # ensure the linker language is reset to CXX
      set_target_properties(${example_target} PROPERTIES LINKER_LANGUAGE CXX)

      # libraries to link against
      if(NOT (${CMAKE_${lang}_COMPILER} STREQUAL ${MPI_CXX_COMPILER}))
        target_include_directories(${example_target} PRIVATE ${MPI_CXX_INCLUDE_DIRS})
        target_link_libraries(${example_target} PRIVATE ${OTHER_LIBS} ${SUNDIALS_LIBS} ${MPI_CXX_LIBRARIES})
      else()
        target_link_libraries(${example_target} PRIVATE ${OTHER_LIBS} ${SUNDIALS_LIBS})
      endif()

      # compile definitions
      if(NOT "${example_defines}" STREQUAL "")
        target_compile_definitions(${example_target} PUBLIC ${example_defines})
      endif()
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
      MPI_NPROCS ${number_of_tasks}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${test_name}.out
      EXAMPLE_TYPE ${example_type}
      FLOAT_PRECISION ${example_precision})

    if(EXAMPLES_INSTALL)
      # Find all .out files for this example
      file(GLOB example_out ${example_target}*.out)

      # Install example source and out files
      install(FILES ${example}.cpp ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/CXX_parallel)
    endif()
  endforeach()
endmacro()

# Specify libraries to link against for each configuration
# and build the examples.
set(SUNDIALS_LIBS sundials_kinsol sundials_nvecparallel)
set(OTHER_LIBS )
build_examples(serial_examples CXX)

if(EXAMPLES_INSTALL)
  # Install the header files
  foreach(extrafile ${KINSOL_headers})
    install(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/CXX_parallel)
  endforeach()

  # Install the extra files
  foreach(extrafile ${KINSOL_extras})
    install(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/CXX_parallel)
  endforeach()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "KINSOL")
  set(EXAMPLES_CMAKE_TARGETS "SUNDIALS::kinsol;SUNDIALS::nvecparallel")
  set(examples_to_install "${serial_examples}")
  examples2string(examples_to_install EXAMPLES)

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_MPI_ex.in
    ${PROJECT_BINARY_DIR}/examples/kinsol/CXX_parallel/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/kinsol/CXX_parallel/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/CXX_parallel
    )

  # add test_install target
  sundials_add_test_install(kinsol CXX_parallel)
endif()
