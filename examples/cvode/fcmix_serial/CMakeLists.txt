# ---------------------------------------------------------------
# $Revision: 1.1 $
# $Date: 2007-01-29 17:40:00 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the FCVODE serial examples


# Add variable fcvode_examples with the names of the serial FCVODE examples

SET(fcvode_examples
  fcvbanx
  fcvdenx
  fcvkryx
  fcvkryx_bp
  )

# Add variable fcvode_bl_examples with the names of the serial FCVODE examples
# that use Lapack

SET(fcvode_bl_examples
  fcvdenx_lap
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODE_LIB sundials_cvode_static)
  SET(CVODE_BL_LIB sundials_cvode_bl_static)
  SET(NVECS_LIB sundials_nvecserial_static)
  SET(FNVECS_LIB sundials_fnvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODE_LIB sundials_cvode_shared)
  SET(CVODE_BL_LIB sundials_cvode_bl_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
  SET(FNVECS_LIB sundials_fnvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Only static FCMIX libraries are available

SET(FCVODE_LIB sundials_fcvode_static)
SET(FCVODE_BL_LIB sundials_fcvode_bl_static)

# Add the build and install targets for each FCVODE example

FOREACH(example ${fcvode_examples})
  # add build target for example
  ADD_EXECUTABLE(${example} ${example}.f)
  # link to the sundials libraries
  TARGET_LINK_LIBRARIES(${example} 
    ${FCVODE_LIB} ${CVODE_LIB}
    ${FNVECS_LIB} ${NVECS_LIB}
    ${EXTRA_LINK_LIBS})
  # install example source and output files
  INSTALL(FILES ${example}.f ${example}.out DESTINATION examples/cvode/fcmix_serial)
ENDFOREACH(example ${fcvode_examples})

# If Lapack support is enabled, add the build and install targets for
# the examples using Lapack

IF(LAPACK_OK)
  FOREACH(example ${fcvode_bl_examples})
    ADD_EXECUTABLE(${example} ${example}.f)
    TARGET_LINK_LIBRARIES(${example} 
      ${FCVODE_BL_LIB} ${CVODE_BL_LIB}
      ${FCVODE_LIB} ${CVODE_LIB}
      ${FNVECS_LIB} ${NVECS_LIB}
      ${EXTRA_LINK_LIBS})
    TARGET_LINK_LIBRARIES(${example} ${LAPACK_LIBRARIES})
    INSTALL(FILES ${example}.f ${example}.out DESTINATION examples/cvode/fcmix_serial)
  ENDFOREACH(example ${fcvode_examples})
ENDIF(LAPACK_OK)

# Regardless of the platform we're on, we will generate and install 
# CMakeLists.txt file for building the examples. This file  can then 
# be used as a template for the user's own programs.

# generate CMakelists.txt in the binary directory
IF(LAPACK_OK)

  FILE(WRITE ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_serial/CMakeLists.txt
    "PROJECT(FCVODE_serial_examples Fortran)\n"
    "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
    "SET(fcvode_examples ${fcvode_examples})\n"
    "SET(fcvode_bl_examples ${fcvode_bl_examples})\n"
    "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
    "FIND_LIBRARY(SUNDIALS_FCVODE_BL_LIB \n"
    "  sundials_fcvode_bl \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"FCVODE Lapack wrapper library\")\n"
    "FIND_LIBRARY(SUNDIALS_CVODE_BL_LIB\n"
    "  sundials_cvode_bl \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"CVODE Lapack wrapper library\")\n"
    "FIND_LIBRARY(SUNDIALS_FCVODE_LIB\n"
    "  sundials_fcvode \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"FCVODE wrapper library\")\n"
    "FIND_LIBRARY(SUNDIALS_CVODE_LIB\n"
    "  sundials_cvode \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"CVODE library\")\n"
    "FIND_LIBRARY(SUNDIALS_FNVECS_LIB \n"
    "  sundials_fnvecserial \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"FNVECSERIAL wrapper library\")\n"
    "FIND_LIBRARY(SUNDIALS_NVECS_LIB \n"
    "  sundials_nvecserial \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"NVECSERIAL library\")\n"
    "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
    "SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE STRING \"Lapack libraries\")\n"
    "FOREACH(example \${fcvode_examples})\n"
    "  ADD_EXECUTABLE(\${example} \${example}.f)\n"
    "  TARGET_LINK_LIBRARIES(\${example}\n"
    "    \${SUNDIALS_FCVODE_LIB} \${SUNDIALS_CVODE_LIB}\n"
    "    \${SUNDIALS_FNVECS_LIB} \${SUNDIALS_NVECS_LIB}\n"
    "    \${SUNDIALS_EXTRA_LIB})\n"
    "ENDFOREACH(example \${fcvode_examples})\n"
    "FOREACH(example \${fcvode_bl_examples})\n"
    "  ADD_EXECUTABLE(\${example} \${example}.f)\n"
    "  TARGET_LINK_LIBRARIES(\${example}\n"
    "    \${SUNDIALS_FCVODE_BL_LIB} \${SUNDIALS_CVODE_BL_LIB}\n"
    "    \${SUNDIALS_FCVODE_LIB} \${SUNDIALS_CVODE_LIB}\n"
    "    \${SUNDIALS_FNVECS_LIB} \${SUNDIALS_NVECS_LIB}\n"
    "    \${SUNDIALS_EXTRA_LIB})\n"
    "  TARGET_LINK_LIBRARIES(\${example} \${LAPACK_LIBRARIES})\n"
    "ENDFOREACH(example \${fcvode_bl_examples})\n"
    )

ELSE(LAPACK_OK)

  FILE(WRITE ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_serial/CMakeLists.txt
    "PROJECT(FCVODE_serial_examples Fortran)\n"
    "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
    "SET(fcvode_examples ${fcvode_examples})\n"
    "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
    "FIND_LIBRARY(SUNDIALS_FCVODE_LIB\n"
    "  sundials_fcvode \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"FCVODE wrapper library\")\n"
    "FIND_LIBRARY(SUNDIALS_CVODE_LIB\n"
    "  sundials_cvode \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"CVODE library\")\n"
    "FIND_LIBRARY(SUNDIALS_FNVECS_LIB\n"
    "  sundials_fnvecserial \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"FNVECSERIAL wrapper library\")\n"
    "FIND_LIBRARY(SUNDIALS_NVECS_LIB\n"
    "  sundials_nvecserial \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"NVECSERIAL library\")\n"
    "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
    "FOREACH(example \${fcvode_examples})\n"
    "  ADD_EXECUTABLE(\${example} \${example}.f)\n"
    "  TARGET_LINK_LIBRARIES(\${example}\n"
    "    \${SUNDIALS_FCVODE_LIB} \${SUNDIALS_CVODE_LIB}\n"
    "    \${SUNDIALS_FNVECS_LIB} \${SUNDIALS_NVECS_LIB}\n"
    "    \${SUNDIALS_EXTRA_LIB})\n"
    "ENDFOREACH(example \${fcvode_examples})\n"
    )

ENDIF(LAPACK_OK)

# install CMakelists.txt
INSTALL(
  FILES ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_serial/CMakeLists.txt
  DESTINATION examples/cvode/fcmix_serial 
  )

# On UNIX-type platforms, we have the option of generating and instaling
# a makefile for building the examples. This makefile can then be used
# as a template for the user's own programs.

IF(GENERATE_EX_MAKEFILES)

  # generate Makefile and place it in the binary dir
  IF(LAPACK_OK)
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/cvode/fcmix_serial/Makefile_ex_bl.in
      ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_serial/Makefile_ex
      @ONLY
      )
  ELSE(LAPACK_OK)
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/cvode/fcmix_serial/Makefile_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_serial/Makefile_ex
      @ONLY
      )
  ENDIF(LAPACK_OK)

  # install the configured Makefile_ex as Makefile
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_serial/Makefile_ex 
    DESTINATION examples/cvode/fcmix_serial 
    RENAME Makefile
    )

ENDIF(GENERATE_EX_MAKEFILES)
