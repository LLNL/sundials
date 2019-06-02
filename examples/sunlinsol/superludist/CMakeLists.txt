# ------------------------------------------------------------------------------
# Programmer:  Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# CMakeLists.txt file for the superlu-dist SUNLinearSolver examples
# -----------------------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;type" where the
# 'develop' for examples excluded from 'make test' in releases

# Examples using the SuperLU_MT linear solver
SET(sunlinsol_superludist_examples
  "test_sunlinsol_superludist\;100 1 1 0\;1\;"
  "test_sunlinsol_superludist\;400 1 2 0\;2\;"
  "test_sunlinsol_superludist\;900 3 1 0\;3\;"
  "test_sunlinsol_superludist\;1000 2 2 0\;4\;"
)

# Dependencies for nvector examples
SET(sunlinsol_superludist_dependencies
  test_sunlinsol
  sundials_nvector
  sundials_matrix
  sundials_linearsolver
)

IF(MPI_C_COMPILER AND MPI_CXX_COMPILER)
  SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
ENDIF()

IF(SUPERLUDIST_OpenMP)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

# Add source directory to include directories
INCLUDE_DIRECTORIES(. ..)

# Add the build and install targets for each example
FOREACH(example_tuple ${sunlinsol_superludist_examples})

  # parse the example tuple
  LIST(GET example_tuple 0 example)
  LIST(GET example_tuple 1 example_args)
  LIST(GET example_tuple 2 number_of_nodes)
  LIST(GET example_tuple 3 example_type)

  # This is used to get around DLL linkage issue since we are
  # manually including sundials_nvector.c here, which is normally in
  # a library that is included.  If this is not set build system
  # thinks nvector is externally linked.
  IF(WIN32)
    ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)
  ENDIF(WIN32)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  IF(NOT TARGET ${example})
    # example source files
    ADD_EXECUTABLE(${example} ${example}.cpp
      ../test_sunlinsol.c
      ../../../src/sundials/sundials_matrix.c
      ../../../src/sundials/sundials_linearsolver.c
      ../../../src/sundials/sundials_nvector.c)

    # folder to organize targets in an IDE
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    IF (LINK_LIBRARY_TYPE MATCHES "static")
      TARGET_LINK_LIBRARIES(${example}
        PRIVATE
          sundials_sunlinsolsuperludist_static
          sundials_sunmatrixslunrloc_static
          sundials_sunmatrixdense_static
          sundials_nvecparallel_static
          sundials_nvecserial_static)
    ELSE()
      TARGET_LINK_LIBRARIES(${example}
        PRIVATE
          sundials_sunlinsolsuperludist_shared
          sundials_sunmatrixslunrloc_shared
          sundials_sunmatrixdense_shared
          sundials_nvecparallel_shared
          sundials_nvecserial_shared)
    ENDIF()

  ENDIF()

  # check if example args are provided and set the test name
  IF("${example_args}" STREQUAL "")
    IF("${number_of_nodes}" STREQUAL "")
      SET(test_name ${example})
    ELSE()
      STRING(REGEX REPLACE " " "_" test_name ${example}_${number_of_nodes})
    ENDIF()
  ELSE()
    IF("${number_of_nodes}" STREQUAL "")
      STRING(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    ELSE()
      STRING(REGEX REPLACE " " "_" test_name ${example}_${number_of_nodes}_${example_args})
    ENDIF()
  ENDIF()

  # add example to regression tests
  SUNDIALS_ADD_TEST(${test_name} ${example}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_nodes}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.cpp
      ../test_sunlinsol.h
      ../test_sunlinsol.c
      ../../../src/sundials/sundials_matrix.c
      ../../../src/sundials/sundials_linearsolver.c
      ../../../src/sundials/sundials_nvector.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/superludist)
  ENDIF()

ENDFOREACH(example_tuple ${sunlinsol_superludist_examples})


IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/superludist)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "SUNMatrix_SLUNRloc")
  STRING(CONCAT SUNDIALS_LIBS "-lsundials_cvode "
    "-lsundials_nvecparallel "
    "-lsundials_sunmatrixslunrloc "
    "-lsundials_sunlinsolsuperludist")

  EXAMPLES2STRING(sunlinsol_superludist_examples EXAMPLES)
  EXAMPLES2STRING(sunlinsol_superludist_dependencies EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_superludist_CXX_ex.in
    ${PROJECT_BINARY_DIR}/examples/sunlinsol/superludist/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/superludist/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/superludist
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_superludist_CXX_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunlinsol/superludist/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/superludist/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/superludist
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
