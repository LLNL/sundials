# ---------------------------------------------------------------
# $Revision: 1.3 $
# $Date: 2007-10-25 20:03:35 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the FIDA serial examples


# Add variable ida_examples with the names of the serial FIDA examples

SET(FIDA_examples
  fidaRoberts_dns
  )

# Add variable FIDA_examples_BL with the names of the serial FIDA examples
# that use Lapack

SET(FIDA_examples_BL
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDA_LIB sundials_ida_static)
  SET(IDA_BL_LIB sundials_ida_bl_static)
  SET(NVECS_LIB sundials_nvecserial_static)
  SET(FNVECS_LIB sundials_fnvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDA_LIB sundials_ida_shared)
  SET(IDA_BL_LIB sundials_ida_bl_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
  SET(FNVECS_LIB sundials_fnvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Only static FCMIX libraries are available

SET(FIDA_LIB sundials_fida_static)
SET(FIDA_BL_LIB sundials_fida_bl_static)

# Add the build and install targets for each FIDA example

FOREACH(example ${FIDA_examples})
  # add build target for example
  ADD_EXECUTABLE(${example} ${example}.f)
  # link to the sundials libraries
  TARGET_LINK_LIBRARIES(${example} 
    ${FIDA_LIB} ${IDA_LIB}
    ${FNVECS_LIB} ${NVECS_LIB}
    ${EXTRA_LINK_LIBS})
  IF(EXAMPLES_INSTALL)
    # install example source and output files
    INSTALL(FILES ${example}.f ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/fcmix_serial)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${FIDA_examples})

# If Lapack support is enabled, add the build and install targets for
# the examples using Lapack

IF(LAPACK_OK)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  FOREACH(example ${FIDA_examples_BL})
    ADD_EXECUTABLE(${example} ${example}.f)
    TARGET_LINK_LIBRARIES(${example} 
      ${FIDA_BL_LIB} ${IDA_BL_LIB}
      ${FIDA_LIB} ${IDA_LIB}
      ${FNVECS_LIB} ${NVECS_LIB}
      ${EXTRA_LINK_LIBS})
    TARGET_LINK_LIBRARIES(${example} ${LAPACK_LIBRARIES})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.f ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/fcmix_serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${FIDA_examples_BL})
ENDIF(LAPACK_OK)

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/fcmix_serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "IDA")
  SET(SOLVER_LIB "sundials_ida")
  SET(SOLVER_FLIB "sundials_fida")
  LIST2STRING(FIDA_examples EXAMPLES)
  IF(LAPACK_OK)
    SET(SOLVER_BL_LIB "sundials_ida_bl")
    SET(SOLVER_BL_FLIB "sundials_fida_bl")
    LIST2STRING(FIDA_examples_BL EXAMPLES_BL)
  ELSE(LAPACK_OK)
    SET(EXAMPLES_BL "")
    SET(SOLVER_BL_LIB "")
    SET(SOLVER_BL_FLIB "")
  ENDIF(LAPACK_OK)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_F77_ex.in
      ${PROJECT_BINARY_DIR}/examples/ida/fcmix_serial/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/ida/fcmix_serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/fcmix_serial 
    )

  # On UNIX-type platforms, we have the option of generating and instaling
  # a makefile for building the examples. This makefile can then be used
  # as a template for the user's own programs.

  IF(EXAMPLES_GENERATE_MAKEFILES)

    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_F77_ex.in
      ${PROJECT_BINARY_DIR}/examples/ida/fcmix_serial/Makefile_ex
      @ONLY
      )

    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/ida/fcmix_serial/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/fcmix_serial 
      RENAME Makefile
      )

  ENDIF(EXAMPLES_GENERATE_MAKEFILES)

ENDIF(EXAMPLES_INSTALL)
