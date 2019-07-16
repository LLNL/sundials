# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for the ARKode F2003 custom module examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
SET(FARKODE_examples
  "ark_brusselator1D_f2003\;develop"
  "ark_analytic_complex_f2003\;develop"
  )

SET(FARKODE_tests
  "test_fnvector_complex_mod\;develop"
  "test_fnvector_fortran_mod\;develop"
  "test_fsunmatrix_fortran_mod\;develop"
  "test_fsunlinsol_fortran_mod\;develop"
  )

# note the order matters when auto-generating the installed Makefile
SET(FARKODE_sharedsources
  fnvector_complex_mod.f90
  fnvector_fortran_mod.f90
  fsunmatrix_fortran_mod.f90
  fsunlinsol_fortran_mod.f90
  )

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(ARKODE_LIB sundials_arkode_static)
  SET(FARKODE_LIB sundials_farkode_mod_static)
ELSE()
  SET(ARKODE_LIB sundials_arkode_shared)
  SET(FARKODE_LIB sundials_farkode_mod_shared)
ENDIF()

# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${FARKODE_LIB} ${ARKODE_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
FOREACH(example_tuple ${FARKODE_examples})

  # parse the example tuple
  LIST(GET example_tuple 0 example)
  LIST(GET example_tuple 1 example_type)

  # Install fortran modules to a unique directory to avoid naming collisions
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)

  # example source files
  ADD_EXECUTABLE(${example} ${example}.f90 ${FARKODE_sharedsources})

  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  SUNDIALS_ADD_TEST(${example} ${example}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})

   # install example source and out files
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.f90 ${example}.out ${FARKODE_sharedsources}
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F2003_custom)
  ENDIF()

ENDFOREACH(example_tuple ${FARKODE_examples})

# Add the build and install targets for regression test
FOREACH(example_tuple ${FARKODE_tests})

  # parse the example tuple
  LIST(GET example_tuple 0 example)
  LIST(GET example_tuple 1 example_type)

  # Install fortran modules to a unique directory to avoid naming collisions
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)

  # example source files
  ADD_EXECUTABLE(${example} ${example}.f90 ${FARKODE_sharedsources})

  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  SUNDIALS_ADD_TEST(${example} ${example}
    TEST_ARGS ""
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  # libraries to link against
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})

  # install example source and out files
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.f90 ${FARKODE_sharedsources}
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F2003_custom)
  ENDIF()

ENDFOREACH(example_tuple ${FARKODE_tests})


# create Makfile and CMakeLists.txt for examples
IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F2003_custom)

  # Install the extra files
  FOREACH(extrafile ${ARKODE_extras})
    INSTALL(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F2003_custom)
  ENDFOREACH()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "ARKODE")
  SET(SOLVER_LIB "sundials_arkode")
  SET(SOLVER_FLIB "sundials_farkode_mod")

  EXAMPLES2STRING(FARKODE_examples EXAMPLES)
  EXAMPLES2STRING(FARKODE_tests TESTFILES)
  SET(EXAMPLES "${EXAMPLES} ${TESTFILES}")

  # shared sources with file extensions for generated CMakeLists.txt
  EXAMPLES2STRING(FARKODE_sharedsources EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_F2003_ex.in
    ${PROJECT_BINARY_DIR}/examples/arkode/F2003_custom/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/arkode/F2003_custom/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F2003_custom
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  IF(UNIX)
    # remove file extensions from dependencies for generated Makefile
    string(REPLACE ".f90" "" EXAMPLES_DEPENDENCIES "${EXAMPLES_DEPENDENCIES}")

    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_F2003_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/F2003_custom/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/arkode/F2003_custom/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F2003_custom
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
