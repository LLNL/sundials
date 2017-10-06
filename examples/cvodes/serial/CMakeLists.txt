# ---------------------------------------------------------------
# $Revision: 4824 $
# $Date: 2016-07-24 16:39:04 -0700 (Sun, 24 Jul 2016) $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
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
# ---------------------------------------------------------------
# CMakeLists.txt file for CVODES serial examples

# Add variable CVODES_examples with the names of the serial CVODES examples
# This is a list of tuples, first element is example name.
# Optional second argument is test arguments
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

# Add variable CVODES_examples_BL with the names of the serial CVODES examples
# that use Lapack
SET(CVODES_examples_BL
  cvsAdvDiff_bndL
  cvsRoberts_dnsL
  )

SET(CVODES_examples_KLU
  cvsRoberts_klu
  cvsRoberts_ASAi_klu
  "cvsRoberts_FSA_klu\;-sensi stg1 t"
  )
 
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
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODES_LIB sundials_cvodes_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Set-up linker flags and link libraries

SET(SUNDIALS_LIBS ${CVODES_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
IF(LAPACK_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})
ENDIF(LAPACK_FOUND)

IF(KLU_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${KLU_LIBRARIES})
ENDIF(KLU_FOUND)

IF(SUPERLUMT_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${SUPERLUMT_LIBRARIES})
ENDIF(SUPERLUMT_FOUND)

# Add the build and install targets for each CVODES example

# SGS is there a way to query CMAKE variable to get added executables 
# rather than keeping our own list?
SET(ADDED_EXECUTABLES "")
FOREACH(example_tuple ${CVODES_examples})
  LIST(GET example_tuple 0 example)

  # Only need to add the executable once
  LIST(FIND ADDED_EXECUTABLES ${example} index)
  IF(index EQUAL -1)
    LIST(APPEND ADDED_EXECUTABLES ${example})
    ADD_EXECUTABLE(${example} ${example}.c)
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
  ENDIF(index EQUAL -1)

  LIST(LENGTH example_tuple n)
  IF(n EQUAL 2)
    LIST(GET example_tuple 1 test_args)
    STRING(REGEX REPLACE " " "_" test_name ${example}_${test_args})
    SUNDIALS_ADD_TEST(${test_name} ${example} TEST_ARGS ${test_args})
  ELSE()
    SUNDIALS_ADD_TEST(${example} ${example})
  ENDIF()

  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
  FILE(GLOB example_out ${example}.out*)
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ${example_out} DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example_tuple ${CVODES_examples})

# If Lapack support is enabled, add the build and install targets for
# the examples using Lapack

IF(LAPACK_FOUND)
  FOREACH(example ${CVODES_examples_BL})
    ADD_EXECUTABLE(${example} ${example}.c)
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    SUNDIALS_ADD_TEST(${example} ${example})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    FILE(GLOB example_out ${example}.out*)
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example_out} DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${CVODES_examples_BL})
ENDIF(LAPACK_FOUND)

# If KLU support is enabled, add the build and install targets for
# the examples using KLU
IF(KLU_FOUND)
  SET(ADDED_EXECUTABLES "")
  FOREACH(example_tuple ${CVODES_examples_KLU})
    LIST(GET example_tuple 0 example)
  
    # Only need to add the executable once
    LIST(FIND ADDED_EXECUTABLES ${example} index)
    IF(index EQUAL -1)
      LIST(APPEND ADDED_EXECUTABLES ${example})
      ADD_EXECUTABLE(${example} ${example}.c)
      SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    ENDIF(index EQUAL -1)
  
    LIST(LENGTH example_tuple n)
    IF(n EQUAL 2)
      LIST(GET example_tuple 1 test_args)
      STRING(REGEX REPLACE " " "_" test_name ${example}_${test_args})
      SUNDIALS_ADD_TEST(${test_name} ${example} TEST_ARGS ${test_args})
    ELSE()
      SUNDIALS_ADD_TEST(${example} ${example})
    ENDIF()
  
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    FILE(GLOB example_out ${example}.out*)
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example_out} DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example_tuple ${CVODES_examples_KLU})
ENDIF(KLU_FOUND)

# If SUPERLUMT support is enabled, add the build and install targets for
# the examples using SUPERLUMT
IF(SUPERLUMT_FOUND)
  SET(ADDED_EXECUTABLES "")
  FOREACH(example_tuple ${CVODES_examples_SUPERLUMT})
    LIST(GET example_tuple 0 example)
  
    # Only need to add the executable once
    LIST(FIND ADDED_EXECUTABLES ${example} index)
    IF(index EQUAL -1)
      LIST(APPEND ADDED_EXECUTABLES ${example})
      ADD_EXECUTABLE(${example} ${example}.c)
      SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    ENDIF(index EQUAL -1)
  
    LIST(LENGTH example_tuple n)
    IF(n EQUAL 2)
      LIST(GET example_tuple 1 test_args)
      STRING(REGEX REPLACE " " "_" test_name ${example}_${test_args})
      SUNDIALS_ADD_TEST(${test_name} ${example} TEST_ARGS ${test_args})
    ELSE()
      SUNDIALS_ADD_TEST(${example} ${example})
    ENDIF()
  
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    FILE(GLOB example_out ${example}.out*)
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example_out} DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example_tuple ${CVODES_examples_SUPERLUMT})
ENDIF(SUPERLUMT_FOUND)

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "CVODES")
  SET(SOLVER_LIB "sundials_cvodes")
  FOREACH(example_tuple ${CVODES_examples})
    list(GET example_tuple 0 example)
    LIST2STRING(example EXAMPLES)
  ENDFOREACH(example_tuple ${CVODES_examples})

  IF(LAPACK_FOUND)
    LIST2STRING(CVODES_examples_BL EXAMPLES_BL)
  ELSE(LAPACK_FOUND)
    SET(EXAMPLES_BL "")
  ENDIF(LAPACK_FOUND)

  IF(KLU_FOUND)
    FOREACH(example_tuple ${CVODES_examples_KLU})
      list(GET example_tuple 0 example)
      LIST2STRING(example EXAMPLES)
    ENDFOREACH(example_tuple ${CVODES_examples_KLU})
  ENDIF(KLU_FOUND)

  IF(SUPERLUMT_FOUND)
    FOREACH(example_tuple ${CVODES_examples_SUPERLUMT})
      list(GET example_tuple 0 example)
      LIST2STRING(example EXAMPLES)
    ENDFOREACH(example_tuple ${CVODES_examples_SUPERLUMT})
  ENDIF(SUPERLUMT_FOUND)

  STRING (REPLACE ";" " " TMP_STR ${EXAMPLES})
  SET(EXAMPLES ${TMP_STR})
  list(REMOVE_DUPLICATES EXAMPLES)

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
