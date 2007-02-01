# ---------------------------------------------------------------
# $Revision: 1.2 $
# $Date: 2007-02-01 16:42:08 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the FKINSOL serial examples


# Add variable kinsol_examples with the names of the serial FKINSOL examples

SET(fkinsol_examples
  fkinkryx
  )

# Add variable fkinsol_bl_examples with the names of the serial FKINSOL examples
# that use Lapack

SET(fkinsol_bl_examples
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(KINSOL_LIB sundials_kinsol_static)
  SET(KINSOL_BL_LIB sundials_kinsol_bl_static)
  SET(NVECS_LIB sundials_nvecserial_static)
  SET(FNVECS_LIB sundials_fnvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(KINSOL_LIB sundials_kinsol_shared)
  SET(KINSOL_BL_LIB sundials_kinsol_bl_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
  SET(FNVECS_LIB sundials_fnvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Only static FCMIX libraries are available

SET(FKINSOL_LIB sundials_fkinsol_static)
SET(FKINSOL_BL_LIB sundials_fkinsol_bl_static)

# Add the build and install targets for each FKINSOL example

FOREACH(example ${fkinsol_examples})
  # add build target for example
  ADD_EXECUTABLE(${example} ${example}.f)
  # link to the sundials libraries
  TARGET_LINK_LIBRARIES(${example} 
    ${FKINSOL_LIB} ${KINSOL_LIB}
    ${FNVECS_LIB} ${NVECS_LIB}
    ${EXTRA_LINK_LIBS})
  IF(EXAMPLES_INSTALL)
    # install example source and output files
    INSTALL(FILES ${example}.f ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/fcmix_serial)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${fkinsol_examples})

# If Lapack support is enabled, add the build and install targets for
# the examples using Lapack

IF(LAPACK_OK)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  FOREACH(example ${fkinsol_bl_examples})
    ADD_EXECUTABLE(${example} ${example}.f)
    TARGET_LINK_LIBRARIES(${example} 
      ${FKINSOL_BL_LIB} ${KINSOL_BL_LIB}
      ${FKINSOL_LIB} ${KINSOL_LIB}
      ${FNVECS_LIB} ${NVECS_LIB}
      ${EXTRA_LINK_LIBS})
    TARGET_LINK_LIBRARIES(${example} ${LAPACK_LIBRARIES})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.f ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/fcmix_serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${fkinsol_examples})
ENDIF(LAPACK_OK)

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/fcmix_serial)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  IF(LAPACK_OK)

    FILE(WRITE ${PROJECT_BINARY_DIR}/examples/kinsol/fcmix_serial/CMakeLists.txt
      "PROJECT(FKINSOL_serial_examples Fortran)\n"
      "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
      "SET(fkinsol_examples ${fkinsol_examples})\n"
      "SET(fkinsol_bl_examples ${fkinsol_bl_examples})\n"
      "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
      "FIND_LIBRARY(SUNDIALS_FKINSOL_BL_LIB \n"
      "  sundials_fkinsol_bl \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"FKINSOL Lapack wrapper library\")\n"
      "FIND_LIBRARY(SUNDIALS_KINSOL_BL_LIB\n"
      "  sundials_kinsol_bl \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"KINSOL Lapack wrapper library\")\n"
      "FIND_LIBRARY(SUNDIALS_FKINSOL_LIB\n"
      "  sundials_fkinsol \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"FKINSOL wrapper library\")\n"
      "FIND_LIBRARY(SUNDIALS_KINSOL_LIB\n"
      "  sundials_kinsol \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"KINSOL library\")\n"
      "FIND_LIBRARY(SUNDIALS_FNVECS_LIB \n"
      "  sundials_fnvecserial \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"FNVECSERIAL wrapper library\")\n"
      "FIND_LIBRARY(SUNDIALS_NVECS_LIB \n"
      "  sundials_nvecserial \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"NVECSERIAL library\")\n"
      "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
      "SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE STRING \"Lapack libraries\")\n"
      "FOREACH(example \${fkinsol_examples})\n"
      "  ADD_EXECUTABLE(\${example} \${example}.f)\n"
      "  TARGET_LINK_LIBRARIES(\${example}\n"
      "    \${SUNDIALS_FKINSOL_LIB} \${SUNDIALS_KINSOL_LIB}\n"
      "    \${SUNDIALS_FNVECS_LIB} \${SUNDIALS_NVECS_LIB}\n"
      "    \${SUNDIALS_EXTRA_LIB})\n"
      "ENDFOREACH(example \${fkinsol_examples})\n"
      "SET(CMAKE_EXE_LINKER_FLAGS \"\${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}\")\n"
      "FOREACH(example \${fkinsol_bl_examples})\n"
      "  ADD_EXECUTABLE(\${example} \${example}.f)\n"
      "  TARGET_LINK_LIBRARIES(\${example}\n"
      "    \${SUNDIALS_FKINSOL_BL_LIB} \${SUNDIALS_KINSOL_BL_LIB}\n"
      "    \${SUNDIALS_FKINSOL_LIB} \${SUNDIALS_KINSOL_LIB}\n"
      "    \${SUNDIALS_FNVECS_LIB} \${SUNDIALS_NVECS_LIB}\n"
      "    \${SUNDIALS_EXTRA_LIB})\n"
      "  TARGET_LINK_LIBRARIES(\${example} \${LAPACK_LIBRARIES})\n"
      "ENDFOREACH(example \${fkinsol_bl_examples})\n"
      )

  ELSE(LAPACK_OK)

    FILE(WRITE ${PROJECT_BINARY_DIR}/examples/kinsol/fcmix_serial/CMakeLists.txt
      "PROJECT(FKINSOL_serial_examples Fortran)\n"
      "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
      "SET(fkinsol_examples ${fkinsol_examples})\n"
      "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
      "FIND_LIBRARY(SUNDIALS_FKINSOL_LIB\n"
      "  sundials_fkinsol \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"FKINSOL wrapper library\")\n"
      "FIND_LIBRARY(SUNDIALS_KINSOL_LIB\n"
      "  sundials_kinsol \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"KINSOL library\")\n"
      "FIND_LIBRARY(SUNDIALS_FNVECS_LIB\n"
      "  sundials_fnvecserial \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"FNVECSERIAL wrapper library\")\n"
      "FIND_LIBRARY(SUNDIALS_NVECS_LIB\n"
      "  sundials_nvecserial \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"NVECSERIAL library\")\n"
      "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
      "FOREACH(example \${fkinsol_examples})\n"
      "  ADD_EXECUTABLE(\${example} \${example}.f)\n"
      "  TARGET_LINK_LIBRARIES(\${example}\n"
      "    \${SUNDIALS_FKINSOL_LIB} \${SUNDIALS_KINSOL_LIB}\n"
      "    \${SUNDIALS_FNVECS_LIB} \${SUNDIALS_NVECS_LIB}\n"
      "    \${SUNDIALS_EXTRA_LIB})\n"
      "ENDFOREACH(example \${fkinsol_examples})\n"
      )

  ENDIF(LAPACK_OK)

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/kinsol/fcmix_serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/fcmix_serial 
    )

  # On UNIX-type platforms, we have the option of generating and instaling
  # a makefile for building the examples. This makefile can then be used
  # as a template for the user's own programs.

  IF(EXAMPLES_GENERATE_MAKEFILES)

    # generate Makefile and place it in the binary dir
    IF(LAPACK_OK)
      CONFIGURE_FILE(
        ${PROJECT_SOURCE_DIR}/examples/kinsol/fcmix_serial/Makefile_ex_bl.in
        ${PROJECT_BINARY_DIR}/examples/kinsol/fcmix_serial/Makefile_ex
        @ONLY
        )
    ELSE(LAPACK_OK)
      CONFIGURE_FILE(
        ${PROJECT_SOURCE_DIR}/examples/kinsol/fcmix_serial/Makefile_ex.in
        ${PROJECT_BINARY_DIR}/examples/kinsol/fcmix_serial/Makefile_ex
        @ONLY
        )
    ENDIF(LAPACK_OK)

    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/kinsol/fcmix_serial/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/fcmix_serial 
      RENAME Makefile
      )

  ENDIF(EXAMPLES_GENERATE_MAKEFILES)

ENDIF(EXAMPLES_INSTALL)
