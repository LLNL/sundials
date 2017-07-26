# ---------------------------------------------------------------
# Programmer(s): Daniel Reynolds @ SMU
#                David Gardner, Carol Woodward, Slaven Peles @ LLNL
# -----------------------------------------------------------------
# LLNS/SMU Copyright Start
# Copyright (c) 2017, Southern Methodist University and 
# Lawrence Livermore National Security
#
# This work was performed under the auspices of the U.S. Department 
# of Energy by Southern Methodist University and Lawrence Livermore 
# National Laboratory under Contract DE-AC52-07NA27344.
# Produced at Southern Methodist University and the Lawrence 
# Livermore National Laboratory.
#
# All rights reserved.
# For details, see the LICENSE file.
# LLNS/SMU Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for sunmatrix examples

# Add variable sunmatrix_diagonal_examples with the names of the sunmatrix examples

SET(sunmatrix_diagonal_examples
  test_sunmatrix_diagonal
)

SET(sunmatrix_diagonal_dependencies
  test_sunmatrix
  sundials_nvector
  sundials_sunmatrix
  sundials_sunmatrixdiagonal
)

# Add source directory to include directories
INCLUDE_DIRECTORIES(. ..)

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(NVECS_LIB sundials_nvecserial_static)
  SET(NVECP_LIB sundials_nvecparallel_static)
  SET(SUNMATS_LIB sundials_sunmatrixdiagonal_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(NVECS_LIB sundials_nvecserial_shared)
  SET(NVECP_LIB sundials_nvecsparallel_shared)
  SET(SUNMATS_LIB sundials_sunmatrixdiagonal_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Set-up linker flags and link libraries

SET(SUNDIALS_LIBS ${NVECS_LIB} ${NVECP_LIB} ${SUNMATS_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each sunmatrix example

FOREACH(example ${sunmatrix_diagonal_examples})

  IF(WIN32)
    ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)
  ENDIF(WIN32)

  ADD_EXECUTABLE(${example} ${example}.c ../test_sunmatrix.c ../../../src/sundials/sundials_matrix.c ../../../src/sundials/sundials_nvector.c)
  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

  SUNDIALS_ADD_TEST( ${example}_ser_100_0  ${example} TEST_ARGS 100 0)
  SUNDIALS_ADD_TEST( ${example}_ser_1000_0 ${example} TEST_ARGS 1000 0)
  SUNDIALS_ADD_TEST( ${example}_par_100_0  ${example} MPI_NPROCS 4 TEST_ARGS 100 0)
  SUNDIALS_ADD_TEST( ${example}_par_1000_0 ${example} MPI_NPROCS 4 TEST_ARGS 1000 0)

  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})

  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ../test_sunmatrix.c ../test_sunmatrix.h ../../../src/sundials/sundials_matrix.c ../../../src/sundials/sundials_nvector.c DESTINATION ${EXAMPLES_INSTALL_PATH}/sunmatrix/diagonal)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${sunmatrix_diagonal_examples})

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/sunmatrix/diagonal)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER_LIB "sundials_sunmatrixdiagonal")
  LIST2STRING(sunmatrix_diagonal_examples EXAMPLES)
  LIST2STRING(sunmatrix_examples_dependencies EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunmatrix/diagonal/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/sunmatrix/diagonal/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/sunmatrix/diagonal
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunmatrix/diagonal/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/sunmatrix/diagonal/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunmatrix/diagonal
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
