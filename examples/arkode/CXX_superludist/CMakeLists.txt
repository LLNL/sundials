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
# -----------------------------------------------------------------
# CMakeLists.txt file for ARKODE C++ SuperLU_DIST examples
# -----------------------------------------------------------------

# Example lists are tuples "name\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
LIST(APPEND ARKODE_examples "ark_brusselator1D_FEM_sludist\;1\;1\;develop")

# Auxiliary files to install
SET(ARKODE_extras
  )

IF(MPI_CXX_COMPILER)
  # use MPI wrapper as the compiler
  SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  # disable C++ extensions (for known wrappers)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX -DLAM_BUILDING")
ELSE()
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(ARKODE_LIB sundials_arkode_static)
  SET(SUNMAT_LIB sundials_sunmatrixslunrloc_static)
  SET(SUNLS_LIB sundials_sunlinsolsuperludist_static)
  SET(NVECP_LIB sundials_nvecopenmp_static)
ELSE()
  SET(ARKODE_LIB sundials_arkode_shared)
  SET(SUNMAT_LIB sundials_sunmatrixslunrloc_shared)
  SET(SUNLS_LIB sundials_sunlinsolsuperludist_shared)
  SET(NVECP_LIB sundials_nvecopenmp_shared)
ENDIF()

# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${ARKODE_LIB} ${NVECP_LIB} ${SUNMAT_LIB} ${SUNLS_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
FOREACH(example_tuple ${ARKODE_examples})

  # parse the example tuple
  LIST(GET example_tuple 0 example)
  LIST(GET example_tuple 1 number_of_nodes)
  LIST(GET example_tuple 2 number_of_tasks)
  LIST(GET example_tuple 3 example_type)

  # example source files
  ADD_EXECUTABLE(${example} ${example}.cpp)

  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  SUNDIALS_ADD_TEST(${example} ${example}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})

  IF(NOT MPI_CXX_COMPILER)
    TARGET_LINK_LIBRARIES(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  ENDIF()

  # install example source and out files
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.cpp ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_superludist)
  ENDIF()

ENDFOREACH(example_tuple ${ARKODE_examples})


# create Makfile and CMakeLists.txt for examples
IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_superludist)

  # Install the extra files
  FOREACH(extrafile ${ARKODE_extras})
    INSTALL(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_superludist)
  ENDFOREACH()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "ARKODE")
  STRING(CONCAT SUNDIALS_LIBS "-lsundials_arkode "
    "-lsundials_sunmatrixslunrloc "
    "-lsundials_sunlinsolsuperludist"
    "-lsundials_nvecopenmp")

  EXAMPLES2STRING(ARKODE_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_superludist_CXX_ex.in
    ${PROJECT_BINARY_DIR}/examples/arkode/CXX_superludist/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/arkode/CXX_superludist/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_superludist
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_superludist_CXX_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/CXX_superludist/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/arkode/CXX_superludist/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_superludist
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
