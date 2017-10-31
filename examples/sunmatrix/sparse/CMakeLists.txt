# ---------------------------------------------------------------
# Programmer(s): Daniel Reynolds @ SMU
#                David Gardner @ LLNL
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

# Add variable sunmatrix_sparse_examples with the names of the sunmatrix examples

SET(sunmatrix_sparse_examples
  test_sunmatrix_sparse
)

SET(sunmatrix_sparse_dependencies
  test_sunmatrix
  sundials_nvector
  sundials_matrix
)

# Add source directory to include directories
INCLUDE_DIRECTORIES(. ..)

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(NVECS_LIB sundials_nvecserial_static)
  SET(SUNMATS_LIB sundials_sunmatrixdense_static)
  LIST(APPEND SUNMATS_LIB sundials_sunmatrixband_static)
  LIST(APPEND SUNMATS_LIB sundials_sunmatrixsparse_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(NVECS_LIB sundials_nvecserial_shared)
  SET(SUNMATS_LIB sundials_sunmatrixdense_shared)
  LIST(APPEND SUNMATS_LIB sundials_sunmatrixband_shared)
  LIST(APPEND SUNMATS_LIB sundials_sunmatrixsparse_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Set-up linker flags and link libraries

SET(SUNDIALS_LIBS ${NVECS_LIB} ${SUNMATS_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each sunmatrix example

FOREACH(example ${sunmatrix_sparse_examples})

  IF(WIN32)
    ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)
  ENDIF(WIN32)

  ADD_EXECUTABLE(${example} ${example}.c ../test_sunmatrix.c ../../../src/sundials/sundials_matrix.c ../../../src/sundials/sundials_nvector.c)
  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

  SUNDIALS_ADD_TEST( ${example}_400_400_0_0 ${example} TEST_ARGS 400 400 0 0)
  SUNDIALS_ADD_TEST( ${example}_450_450_1_0 ${example} TEST_ARGS 450 450 1 0)
  SUNDIALS_ADD_TEST( ${example}_200_1000_0_0 ${example} TEST_ARGS 200 1000 0 0)
  SUNDIALS_ADD_TEST( ${example}_6000_350_0_0 ${example} TEST_ARGS 6000 350 0 0)
  SUNDIALS_ADD_TEST( ${example}_500_5000_1_0 ${example} TEST_ARGS 500 5000 1 0)
  SUNDIALS_ADD_TEST( ${example}_4000_800_1_0 ${example} TEST_ARGS 4000 800 1 0)

  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})

  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ../test_sunmatrix.c ../test_sunmatrix.h ../../../src/sundials/sundials_matrix.c ../../../src/sundials/sundials_nvector.c DESTINATION ${EXAMPLES_INSTALL_PATH}/sunmatrix/sparse)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${sunmatrix_sparse_examples})

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/sunmatrix/sparse)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER_LIB "sundials_sunmatrixsparse")
  SET(LIBS "${LIBS} -lsundials_sunmatrixdense -lsundials_sunmatrixband")

  # Set the link directory for the dense and band sunmatrix libraries
  # The generated CMakeLists.txt does not use find_library() locate them
  SET(EXTRA_LIBS_DIR "${libdir}")

  LIST2STRING(sunmatrix_sparse_examples EXAMPLES)
  LIST2STRING(sunmatrix_sparse_dependencies EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunmatrix/sparse/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/sunmatrix/sparse/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/sunmatrix/sparse
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunmatrix/sparse/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/sunmatrix/sparse/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunmatrix/sparse
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
