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
# CMakeLists.txt file for CVODES serial examples


# Add variable cvodes_examples with the names of the serial CVODES examples

SET(cvodes_examples
  cvsadjbanx
  cvsadjdenx
  cvsadjkryx_int
  cvsadjkryx_pnt
  cvsbanx
  cvsdenx
  cvsdenx_uw
  cvsdirectdem
  cvsfwddenx
  cvsfwdkryx
  cvsfwdnonx
  cvskrydem_lin
  cvskrydem_pre
  cvskryx_bp
  cvskryx
  )

# Add variable cvodes_bl_examples with the names of the serial CVODES examples
# that use Lapack

SET(cvodes_bl_examples
  cvsbanx_lap
  cvsdenx_lap
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODES_LIB sundials_cvodes_static)
  SET(CVODES_BL_LIB sundials_cvodes_bl_static)
  SET(NVECS_LIB sundials_nvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODES_LIB sundials_cvodes_shared)
  SET(CVODES_BL_LIB sundials_cvodes_bl_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Add the build and install targets for each CVODES example

FOREACH(example ${cvodes_examples})
  # add build target for example
  ADD_EXECUTABLE(${example} ${example}.c)
  # link to the sundials libraries
  TARGET_LINK_LIBRARIES(${example} ${CVODES_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
  # find all sample output files for this example
  FILE(GLOB example_out ${example}.out*)
  IF(EXAMPLES_INSTALL)
    # install example files
    INSTALL(FILES ${example}.c ${example_out} DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${cvodes_examples})

# If Lapack support is enabled, add the build and install targets for
# the examples using Lapack

IF(LAPACK_OK)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  FOREACH(example ${cvodes_bl_examples})
    ADD_EXECUTABLE(${example} ${example}.c)
    TARGET_LINK_LIBRARIES(${example} ${CVODES_BL_LIB} ${CVODES_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
    TARGET_LINK_LIBRARIES(${example} ${LAPACK_LIBRARIES})
    FILE(GLOB example_out ${example}.out*)
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example_out} DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${cvodes_bl_examples})
ENDIF(LAPACK_OK)

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  IF(LAPACK_OK)

    FILE(WRITE ${PROJECT_BINARY_DIR}/examples/cvodes/serial/CMakeLists.txt
      "PROJECT(CVODES_serial_examples C)\n"
      "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
      "SET(cvodes_examples ${cvodes_examples})\n"
      "SET(cvodes_bl_examples ${cvodes_bl_examples})\n"
      "SET(SUNDIALS_INC_DIR\n"
      "  ${CMAKE_INSTALL_PREFIX}/include\n"
      "  CACHE STRING\n"
      "  \"Location of SUNDIALS header files\")\n"
      "INCLUDE_DIRECTORIES(\${SUNDIALS_INC_DIR})\n"
      "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
      "FIND_LIBRARY(SUNDIALS_CVODES_BL_LIB\n"
      "  sundials_cvodes_bl \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"CVODES Lapack wrapper library\")\n"
      "FIND_LIBRARY(SUNDIALS_CVODES_LIB\n"
      "  sundials_cvodes \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"CVODES library\")\n"
      "FIND_LIBRARY(SUNDIALS_NVECS_LIB\n"
      "  sundials_nvecserial \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"NVECSERIAL library\")\n"
      "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
      "SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE STRING \"Lapack libraries\")\n"
      "FOREACH(example \${cvodes_examples})\n"
      "  ADD_EXECUTABLE(\${example} \${example}.c)\n"
      "  TARGET_LINK_LIBRARIES(\${example} \${SUNDIALS_CVODES_LIB} \${SUNDIALS_NVECS_LIB} \${SUNDIALS_EXTRA_LIB})\n"
      "ENDFOREACH(example \${cvodes_examples})\n"
      "SET(CMAKE_EXE_LINKER_FLAGS \"\${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}\")\n"
      "FOREACH(example \${cvodes_bl_examples})\n"
      "  ADD_EXECUTABLE(\${example} \${example}.c)\n"
      "  TARGET_LINK_LIBRARIES(\${example} \${SUNDIALS_CVODES_BL_LIB} \${SUNDIALS_CVODES_LIB} \${SUNDIALS_NVECS_LIB} \${SUNDIALS_EXTRA_LIB})\n"
      "  TARGET_LINK_LIBRARIES(\${example} \${LAPACK_LIBRARIES})\n"
      "ENDFOREACH(example \${cvodes_bl_examples})\n"
      )

  ELSE(LAPACK_OK)

    FILE(WRITE ${PROJECT_BINARY_DIR}/examples/cvodes/serial/CMakeLists.txt
      "PROJECT(CVODES_serial_examples C)\n"
      "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
      "SET(cvodes_examples ${cvodes_examples})\n"
      "SET(SUNDIALS_INC_DIR\n"
      "  ${CMAKE_INSTALL_PREFIX}/include\n"
      "  CACHE STRING\n"
      "  \"Location of SUNDIALS header files\")\n"
      "INCLUDE_DIRECTORIES(\${SUNDIALS_INC_DIR})\n"
      "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
      "FIND_LIBRARY(SUNDIALS_CVODES_LIB\n"
      "  sundials_cvodes \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"CVODES library\")\n"
      "FIND_LIBRARY(SUNDIALS_NVECS_LIB\n"
      "  sundials_nvecserial \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"NVECSERIAL library\")\n"
      "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
      "FOREACH(example \${cvodes_examples})\n"
      "  ADD_EXECUTABLE(\${example} \${example}.c)\n"
      "  TARGET_LINK_LIBRARIES(\${example} \${SUNDIALS_CVODES_LIB} \${SUNDIALS_NVECS_LIB} \${SUNDIALS_EXTRA_LIB})\n"
      "ENDFOREACH(example \${cvodes_examples})\n"
      )
    
  ENDIF(LAPACK_OK)

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/cvodes/serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial 
    )

  # On UNIX-type platforms, we have the option of generating and instaling
  # a makefile for building the examples. This makefile can then be used
  # as a template for the user's own programs.

  IF(EXAMPLES_GENERATE_MAKEFILES)

    # generate Makefile and place it in the binary dir
    IF(LAPACK_OK)
      CONFIGURE_FILE(
        ${PROJECT_SOURCE_DIR}/examples/cvodes/serial/Makefile_ex_bl.in
        ${PROJECT_BINARY_DIR}/examples/cvodes/serial/Makefile_ex
        @ONLY
        )
    ELSE(LAPACK_OK)
      CONFIGURE_FILE(
        ${PROJECT_SOURCE_DIR}/examples/cvodes/serial/Makefile_ex.in
        ${PROJECT_BINARY_DIR}/examples/cvodes/serial/Makefile_ex
        @ONLY
        )
    ENDIF(LAPACK_OK)

    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/cvodes/serial/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial 
      RENAME Makefile
      )

  ENDIF(EXAMPLES_GENERATE_MAKEFILES)

ENDIF(EXAMPLES_INSTALL)
