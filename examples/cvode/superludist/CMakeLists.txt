# ---------------------------------------------------------------
# Programmer:  Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for CVODE SuperLU_DIST examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;type" where the
# type is develop for examples excluded from 'make test' in releases

LIST(APPEND CVODE_examples "cvAdvDiff_sludist\;2 1\;2\;develop")

IF(MPI_CXX_COMPILER)
  SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  # disable C++ extensions (for known wrappers)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX -DLAM_BUILDING")
ELSE()
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

IF(SUPERLUDIST_OpenMP)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODE_LIB sundials_cvode_static)
  SET(NVECP_LIB sundials_nvecparallel_static)
  SET(SUNMAT_LIB sundials_sunmatrixslunrloc_static)
  SET(SUNLS_LIB sundials_sunlinsolsuperludist_static)
ELSE()
  SET(CVODE_LIB sundials_cvode_shared)
  SET(NVECP_LIB sundials_nvecparallel_shared)
  SET(SUNMAT_LIB sundials_sunmatrixslunrloc_shared)
  SET(SUNLS_LIB sundials_sunlinsolsuperludist_shared)
ENDIF()

# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${CVODE_LIB} ${NVECP_LIB} ${SUNMAT_LIB} ${SUNLS_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
FOREACH(example_tuple ${CVODE_examples})

  # parse the example tuple
  LIST(GET example_tuple 0 example)
  LIST(GET example_tuple 1 example_args)
  LIST(GET example_tuple 2 number_of_nodes)
  LIST(GET example_tuple 3 example_type)

  # example source files
  ADD_EXECUTABLE(${example} ${example}.cpp)

  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  SUNDIALS_ADD_TEST(${example} ${example}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_nodes}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}_${SUNDIALS_INDEX_SIZE}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})

  IF(NOT MPI_C_COMPILER)
    TARGET_LINK_LIBRARIES(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  ENDIF()

  # install example source and out files
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.cpp ${example}_${SUNDIALS_INDEX_SIZE}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/superludist)
  ENDIF()

ENDFOREACH(example_tuple ${CVODE_examples})


# create Makfile and CMakeLists.txt for examples
IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/superludist)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "CVODE")
  STRING(CONCAT SUNDIALS_LIBS "-lsundials_cvode "
    "-lsundials_nvecparallel "
    "-lsundials_sunmatrixslunrloc "
    "-lsundials_sunlinsolsuperludist")

  EXAMPLES2STRING(CVODE_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_superludist_CXX_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/superludist/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/cvode/superludist/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/superludist
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_superludist_CXX_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/superludist/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/superludist/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/superludist
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
