# ---------------------------------------------------------------
# $Revision: 4957 $
# $Date: 2016-09-23 12:21:47 -0700 (Fri, 23 Sep 2016) $
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
# CMakeLists.txt file for the FCVODE serial examples


# Add variable FCVODE_examples with the names of the serial FCVODE examples

SET(FCVODE_examples
  fcvAdvDiff_bnd
  fcvDiurnal_kry_bp
  fcvDiurnal_kry
  fcvRoberts_dns
  )

SET(FCVODE_examples_KLU
  fcvRoberts_klu
  )

SET(FCVODE_examples_SUPERLUMT
  fcvRoberts_sps
  )

# Add variable fcvode_bl_examples with the names of the serial FCVODE examples
# that use Lapack

SET(FCVODE_examples_BL
  fcvRoberts_dnsL
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODE_LIB sundials_cvode_static)
  SET(NVECS_LIB sundials_nvecserial_static)
  SET(FNVECS_LIB sundials_fnvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODE_LIB sundials_cvode_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
  SET(FNVECS_LIB sundials_fnvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Only static FCMIX libraries are available

SET(FCVODE_LIB sundials_fcvode_static)

# Set-up linker flags and link libraries

SET(SUNDIALS_LIBS ${FCVODE_LIB} ${CVODE_LIB} ${FNVECS_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
IF(LAPACK_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})
ENDIF(LAPACK_FOUND)

IF(KLU_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${KLU_LIBRARIES})
ENDIF(KLU_FOUND)

IF(SUPERLUMT_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${SUPERLUMT_LIBRARIES})
ENDIF(SUPERLUMT_FOUND)

IF(HYPRE_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${HYPRE_LIBRARIES})
ENDIF(HYPRE_FOUND)

# Add the build and install targets for each FCVODE example

FOREACH(example ${FCVODE_examples})
  ADD_EXECUTABLE(${example} ${example}.f)
  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
  SUNDIALS_ADD_TEST(${example} ${example})
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.f ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_serial)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${FCVODE_examples})

# Add the build and install targets for each Lapack FCVODE example (if needed)

IF(LAPACK_FOUND)
  FOREACH(example ${FCVODE_examples_BL})
    ADD_EXECUTABLE(${example} ${example}.f)
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    SUNDIALS_ADD_TEST(${example} ${example})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.f ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${FCVODE_examples_BL})
ENDIF(LAPACK_FOUND)

# If KLU support is enabled, add the build and install targets for
# the examples using KLU
IF(KLU_FOUND)
  FOREACH(example ${FCVODE_examples_KLU})
    ADD_EXECUTABLE(${example} ${example}.f)
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    SUNDIALS_ADD_TEST(${example} ${example})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.f ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${FCVODE_examples_KLU})
ENDIF(KLU_FOUND)

# If SUPERLUMT support is enabled, add the build and install targets for
# the examples using SUPERLUMT
IF(SUPERLUMT_FOUND)
  FOREACH(example ${FCVODE_examples_SUPERLUMT})
    ADD_EXECUTABLE(${example} ${example}.f)
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    SUNDIALS_ADD_TEST(${example} ${example})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.f ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${FCVODE_examples_SUPERLUMT})
ENDIF(SUPERLUMT_FOUND)

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "CVODE")
  SET(SOLVER_LIB "sundials_cvode")
  SET(SOLVER_FLIB "sundials_fcvode")
  LIST2STRING(FCVODE_examples EXAMPLES)
  IF(LAPACK_FOUND)
    LIST2STRING(FCVODE_examples_BL EXAMPLES_BL)
  ELSE(LAPACK_FOUND)
    SET(EXAMPLES_BL "")
  ENDIF(LAPACK_FOUND)
  
  IF(KLU_FOUND)
    LIST2STRING(FCVODE_examples_KLU EXAMPLES)
  ENDIF(KLU_FOUND)
  
  IF(SUPERLUMT_FOUND)
    LIST2STRING(FCVODE_examples_SUPERLUMT EXAMPLES)
  ENDIF(SUPERLUMT_FOUND)
  
  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_F77_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_serial/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_serial 
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_F77_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_serial/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_serial 
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
