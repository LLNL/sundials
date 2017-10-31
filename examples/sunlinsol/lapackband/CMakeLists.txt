# ---------------------------------------------------------------
# Programmer(s): Daniel Reynolds @ SMU
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
# CMakeLists.txt file for sunlinsol examples

# Add variable sunlinsol_lapackband_examples with the names of the sunlinsol examples

SET(sunlinsol_lapackband_examples
  test_sunlinsol_lapackband
)

SET(sunlinsol_lapackband_dependencies
  test_sunlinsol
  sundials_nvector
  sundials_matrix
  sundials_linearsolver
)

# Add source directory to include directories
INCLUDE_DIRECTORIES(. ..)

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(NVEC_LIB sundials_nvecserial_static)
  SET(SUNMAT_LIB sundials_sunmatrixband_static)
  SET(SUNLINSOL_LIB sundials_sunlinsollapackband_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(NVEC_LIB sundials_nvecserial_shared)
  SET(SUNMAT_LIB sundials_sunmatrixband_shared)
  SET(SUNLINSOL_LIB sundials_sunlinsollapackband_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Set-up linker flags and link libraries

SET(SUNDIALS_LIBS ${NVEC_LIB} ${SUNMAT_LIB} ${SUNLINSOL_LIB} ${EXTRA_LINK_LIBS})

LIST(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})

IF(BLAS_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${BLAS_LIBRARIES})
ENDIF(BLAS_FOUND)

# Add the build and install targets for each sunlinsol example

FOREACH(example ${sunlinsol_lapackband_examples})

  IF(WIN32)
    ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)
  ENDIF(WIN32)

  ADD_EXECUTABLE(${example} ${example}.c ../test_sunlinsol.c ../../../src/sundials/sundials_matrix.c ../../../src/sundials/sundials_linearsolver.c ../../../src/sundials/sundials_nvector.c)
  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

  SUNDIALS_ADD_TEST( ${example}_10_2_3_0 ${example} TEST_ARGS 10 2 3 0)
  SUNDIALS_ADD_TEST( ${example}_300_7_4_0 ${example} TEST_ARGS 300 7 4 0)
  SUNDIALS_ADD_TEST( ${example}_1000_8_8_0 ${example} TEST_ARGS 1000 8 8 0)
  SUNDIALS_ADD_TEST( ${example}_5000_3_100_0 ${example} TEST_ARGS 5000 3 100 0)

  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})

  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ../test_sunlinsol.h ../test_sunlinsol.c ../../../src/sundials/sundials_matrix.c ../../../src/sundials/sundials_linearsolver.c ../../../src/sundials/sundials_nvector.c DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/lapackband)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${sunlinsol_lapackband_examples})

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/lapackband)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER_LIB "sundials_sunlinsollapackband")
  SET(LIBS "${LIBS} -lsundials_sunmatrixband -lsundials_sunmatrixdense")

  # Set the link directory for the band and dense sunmatrix libraries
  # The generated CMakeLists.txt does not use find_library() locate them
  SET(EXTRA_LIBS_DIR "${libdir}")

  LIST2STRING(sunlinsol_lapackband_examples EXAMPLES_BL)
  LIST2STRING(sunlinsol_lapackband_dependencies EXAMPLES_DEPENDENCIES_BL)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunlinsol/lapackband/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/lapackband/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/lapackband
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunlinsol/lapackband/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/lapackband/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/lapackband
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
