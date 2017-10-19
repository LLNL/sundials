# ---------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
#                David J. Gardner @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department
# of Energy by Lawrence Livermore National Laboratory in part under
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ------------------------------------------------------------------
# CMakeLists.txt file for CVODES serial examples
# ------------------------------------------------------------------

# Example lists are lists of tuples where the first element is the
# example name and the optional second argument is test arguments

# Examples using SUNDIALS linear solvers
SET(CVODES_examples
  cvsAdvDiff_ASAi_bnd
  cvsAdvDiff_bnd
  "cvsAdvDiff_FSA_non\;-sensi stg t"
  "cvsAdvDiff_FSA_non\;-sensi sim t"
  cvsDirectDemo_ls
  "cvsDiurnal_FSA_kry\;-sensi stg t"
  "cvsDiurnal_FSA_kry\;-sensi sim t"
  cvsDiurnal_kry_bp
  cvsDiurnal_kry
  cvsFoodWeb_ASAi_kry
  cvsFoodWeb_ASAp_kry
  cvsHessian_ASA_FSA
  cvsKrylovDemo_ls
  cvsKrylovDemo_prec
  cvsRoberts_ASAi_dns
  cvsRoberts_dns
  cvsRoberts_dns_uw
  "cvsRoberts_FSA_dns\;-sensi stg1 t"
  )

# Examples using LAPACK linear solvers
SET(CVODES_examples_BL
  cvsAdvDiff_bndL
  cvsRoberts_dnsL
  )

# Examples using KLU linear solver
SET(CVODES_examples_KLU
  cvsRoberts_klu
  cvsRoberts_ASAi_klu
  "cvsRoberts_FSA_klu\;-sensi stg1 t"
  )

# Examples using SuperLU_MT linear solver
SET(CVODES_examples_SUPERLUMT
  cvsRoberts_sps
  cvsRoberts_ASAi_sps
  "cvsRoberts_FSA_sps\;-sensi stg1 t"
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODES_LIB sundials_cvodes_static)
  SET(NVECS_LIB sundials_nvecserial_static)
ELSE()
  SET(CVODES_LIB sundials_cvodes_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
ENDIF()

# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${CVODES_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})

# SGS is there a way to query CMAKE variable to get added executables
# rather than keeping our own list?
SET(ADDED_EXECUTABLES "")

# Add the build and install targets for each example
FOREACH(example_tuple ${CVODES_examples})

  # first item is example
  list(GET example_tuple 0 example)

  # check if example has already been added, only need to add an examples
  # once for testing with different inputs
  LIST(FIND ADDED_EXECUTABLES ${example} index)

  # if the example has not been added yet, add it to the list
  # and add the example source files
  IF(index EQUAL -1)
    LIST(APPEND ADDED_EXECUTABLES ${example})
    ADD_EXECUTABLE(${example} ${example}.c)
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
  ENDIF(index EQUAL -1)

  # Check if an optional 2nd argument is used and add test to regression tests
  LIST(LENGTH example_tuple n)
  IF(n EQUAL 2)
    LIST(GET example_tuple 1 test_args)
    STRING(REGEX REPLACE " " "_" test_name ${example}_${test_args})
    SUNDIALS_ADD_TEST(${test_name} ${example} TEST_ARGS ${test_args})
  ELSE()
    SUNDIALS_ADD_TEST(${example} ${example})
  ENDIF()

  # libraries to link against
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})

  # find all outfile files for this example
  FILE(GLOB example_out ${example}.out*)

  # install example source and out files
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ${example_out}
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
  ENDIF()

ENDFOREACH(example_tuple ${CVODES_examples})


# Add the build and install targets for each LAPACK example (if needed)
IF(LAPACK_FOUND)

  # Sundials LAPACK linear solver modules
  IF(LINK_LIBRARY_TYPE MATCHES "static")
    SET(SUNLINSOLLAPACK_LIBS
      sundials_sunlinsollapackband_static
      sundials_sunlinsollapackdense_static)
  ELSE()
    SET(SUNLINSOLLAPACK_LIBS
      sundials_sunlinsollapackband_shared
      sundials_sunlinsollapackdense_shared)
  ENDIF()

  # LAPACK libraries
  LIST(APPEND SUNLINSOLLAPACK_LIBS ${LAPACK_LIBRARIES})

  # BLAS libraries
  IF(BLAS_FOUND)
    LIST(APPEND SUNLINSOLLAPACK_LIBS ${BLAS_LIBRARIES})
  ENDIF(BLAS_FOUND)

  FOREACH(example ${CVODES_examples_BL})

    # example source files
    ADD_EXECUTABLE(${example} ${example}.c)

    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

    # add to regression tests
    SUNDIALS_ADD_TEST(${example} ${example})

    # libraries to link against
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS} ${SUNLINSOLLAPACK_LIBS})

    # find all outfile files for this example
    FILE(GLOB example_out ${example}.out*)

    # install example source and out files
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
    ENDIF()

  ENDFOREACH(example ${CVODES_examples_BL})

ENDIF(LAPACK_FOUND)


# Add the build and install targets for each KLU example (if needed)
IF(KLU_FOUND)

  # Sundials KLU linear solver module
  IF(LINK_LIBRARY_TYPE MATCHES "static")
    SET(SUNLINSOLKLU_LIBS sundials_sunlinsolklu_static)
  ELSE()
    SET(SUNLINSOLKLU_LIBS sundials_sunlinsolklu_shared)
  ENDIF()

  # KLU libraries
  LIST(APPEND SUNLINSOLKLU_LIBS ${KLU_LIBRARIES})

  SET(ADDED_EXECUTABLES "")

  FOREACH(example_tuple ${CVODES_examples_KLU})

    # first item is example
    LIST(GET example_tuple 0 example)
  
    # check if example has already been added, only need to add an examples
    # once for testing with different inputs
    LIST(FIND ADDED_EXECUTABLES ${example} index)

    # if the example has not been added yet, add it to the list
    # and add the example source files
    IF(index EQUAL -1)
      LIST(APPEND ADDED_EXECUTABLES ${example})
      ADD_EXECUTABLE(${example} ${example}.c)
      SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    ENDIF(index EQUAL -1)
  
    # Check if an optional 2nd argument is used and add test to regression tests
    LIST(LENGTH example_tuple n)
    IF(n EQUAL 2)
      LIST(GET example_tuple 1 test_args)
      STRING(REGEX REPLACE " " "_" test_name ${example}_${test_args})
      SUNDIALS_ADD_TEST(${test_name} ${example} TEST_ARGS ${test_args})
    ELSE()
      SUNDIALS_ADD_TEST(${example} ${example})
    ENDIF()
  
    # libraries to link against
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS} ${SUNLINSOLKLU_LIBS})

    # find all outfile files for this example
    FILE(GLOB example_out ${example}.out*)

    # install example source and out files
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
    ENDIF()

  ENDFOREACH(example_tuple ${CVODES_examples_KLU})

ENDIF(KLU_FOUND)


# Add the build and install targets for each SuperLU_MT example (if needed)
IF(SUPERLUMT_FOUND)

  # Sundials SuperLU_MT linear solver module
  IF(LINK_LIBRARY_TYPE MATCHES "static")
    SET(SUNLINSOLSLUMT_LIBS sundials_sunlinsolsuperlumt_static)
  ELSE()
    SET(SUNLINSOLSLUMT_LIBS sundials_sunlinsolsuperlumt_shared)
  ENDIF()

  # SuperLU_MT libraries
  LIST(APPEND SUNLINSOLSLUMT_LIBS ${SUPERLUMT_LIBRARIES})

  # BLAS libraries
  IF(BLAS_FOUND)
    LIST(APPEND SUNLINSOLSLUMT_LIBS ${BLAS_LIBRARIES})
  ENDIF(BLAS_FOUND)

  SET(ADDED_EXECUTABLES "")

  FOREACH(example_tuple ${CVODES_examples_SUPERLUMT})

    # first item is example
    LIST(GET example_tuple 0 example)
  
    # check if example has already been added, only need to add an examples
    # once for testing with different inputs
    LIST(FIND ADDED_EXECUTABLES ${example} index)

    # if the example has not been added yet, add it to the list
    # and add the example source files
    IF(index EQUAL -1)
      LIST(APPEND ADDED_EXECUTABLES ${example})
      ADD_EXECUTABLE(${example} ${example}.c)
      SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    ENDIF(index EQUAL -1)
  
    # Check if an optional 2nd argument is used and add test to regression tests
    LIST(LENGTH example_tuple n)
    IF(n EQUAL 2)
      LIST(GET example_tuple 1 test_args)
      STRING(REGEX REPLACE " " "_" test_name ${example}_${test_args})
      SUNDIALS_ADD_TEST(${test_name} ${example} TEST_ARGS ${test_args})
    ELSE()
      SUNDIALS_ADD_TEST(${example} ${example})
    ENDIF()
  
    # libraries to link against
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS} ${SUNLINSOLSLUMT_LIBS})

    # find all outfile files for this example
    FILE(GLOB example_out ${example}.out*)

    # install example source and out files
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
    ENDIF()

  ENDFOREACH(example_tuple ${CVODES_examples_SUPERLUMT})

ENDIF(SUPERLUMT_FOUND)


# create Makfile and CMakeLists.txt for examples
IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "CVODES")
  SET(SOLVER_LIB "sundials_cvodes")

  FOREACH(example_tuple ${CVODES_examples})
    list(GET example_tuple 0 example)
    list(APPEND TMP_EXAMPLES ${example})
  ENDFOREACH()
  list(REMOVE_DUPLICATES TMP_EXAMPLES)
  LIST2STRING(TMP_EXAMPLES EXAMPLES)

  IF(LAPACK_FOUND)
    LIST2STRING(CVODES_examples_BL EXAMPLES_BL)
  ELSE()
    SET(EXAMPLES_BL "")
  ENDIF()

  IF(KLU_FOUND)
    FOREACH(example_tuple ${CVODES_examples_KLU})
      list(GET example_tuple 0 example)
      list(APPEND TMP_EXAMPLES_KLU ${example})
    ENDFOREACH()
    list(REMOVE_DUPLICATES TMP_EXAMPLES_KLU)
  ELSE()
    SET(EXAMPLES_KLU "")
  ENDIF()
  LIST2STRING(TMP_EXAMPLES_KLU EXAMPLES_KLU)

  IF(SUPERLUMT_FOUND)
    FOREACH(example_tuple ${CVODES_examples_SUPERLUMT})
      list(GET example_tuple 0 example)
      list(APPEND TMP_EXAMPLES_SLUMT ${example})
    ENDFOREACH()
    list(REMOVE_DUPLICATES TMP_EXAMPLES_SLUMT)
  ELSE()
    SET(EXAMPLES_SLUMTU "")
  ENDIF()
  LIST2STRING(TMP_EXAMPLES_SLUMT EXAMPLES_SLUMT)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvodes/serial/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/cvodes/serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial 
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvodes/serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/cvodes/serial/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
