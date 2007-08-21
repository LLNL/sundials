# ---------------------------------------------------------------
# $Revision: 1.6 $
# $Date: 2007-08-21 16:55:07 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for IDAS serial examples


# Add variable idas_examples with the names of the serial IDAS examples

SET(idas_examples
  idasHeat2D_b
  idasFoodWeb_b
  idasRobertson
  idasKrylovDemo
  idasHeat2D_k
  idasSliderCrank
  idasSliderCrankFSA
  idasRobertsonFSA
  idasRobertsonASA
  idasRobertsonHessian
  idasAkzoNobelASA
  )

# Add variable idas_bl_examples with the names of the serial IDAS examples
# that use Lapack

SET(idas_bl_examples
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDAS_LIB sundials_idas_static)
  SET(IDAS_BL_LIB sundials_idas_bl_static)
  SET(NVECS_LIB sundials_nvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDAS_LIB sundials_idas_shared)
  SET(IDAS_BL_LIB sundials_idas_bl_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Add the build and install targets for each IDAS example

FOREACH(example ${idas_examples})
  # add build target for example
  ADD_EXECUTABLE(${example} ${example}.c)
  # link to the sundials libraries
  TARGET_LINK_LIBRARIES(${example} ${IDAS_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
  IF(EXAMPLES_INSTALL)
    # install example files
    INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${idas_examples})

# If Lapack support is enabled, add the build and install targets for
# the examples using Lapack

IF(LAPACK_OK)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  FOREACH(example ${idas_bl_examples})
    ADD_EXECUTABLE(${example} ${example}.c)
    TARGET_LINK_LIBRARIES(${example} ${IDAS_BL_LIB} ${IDAS_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
    TARGET_LINK_LIBRARIES(${example} ${LAPACK_LIBRARIES})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${idas_bl_examples})
ENDIF(LAPACK_OK)

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  IF(LAPACK_OK)

    FILE(WRITE ${PROJECT_BINARY_DIR}/examples/idas/serial/CMakeLists.txt
      "PROJECT(IDAS_serial_examples C)\n"
      "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
      "SET(idas_examples ${idas_examples})\n"
      "SET(idas_bl_examples ${idas_bl_examples})\n"
      "SET(SUNDIALS_INC_DIR\n"
      "  ${CMAKE_INSTALL_PREFIX}/include\n"
      "  CACHE STRING\n"
      "  \"Location of SUNDIALS header files\")\n"
      "INCLUDE_DIRECTORIES(\${SUNDIALS_INC_DIR})\n"
      "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
      "FIND_LIBRARY(SUNDIALS_IDAS_BL_LIB\n"
      "  sundials_idas_bl \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"IDAS Lapack wrapper library\")\n"
      "FIND_LIBRARY(SUNDIALS_IDAS_LIB\n"
      "  sundials_idas \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"IDAS library\")\n"
      "FIND_LIBRARY(SUNDIALS_NVECS_LIB\n"
      "  sundials_nvecserial \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"NVECSERIAL library\")\n"
      "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
      "SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE STRING \"Lapack libraries\")\n"
      "FOREACH(example \${idas_examples})\n"
      "  ADD_EXECUTABLE(\${example} \${example}.c)\n"
      "  TARGET_LINK_LIBRARIES(\${example} \${SUNDIALS_IDAS_LIB} \${SUNDIALS_NVECS_LIB} \${SUNDIALS_EXTRA_LIB})\n"
      "ENDFOREACH(example \${idas_examples})\n"
      "SET(CMAKE_EXE_LINKER_FLAGS \"\${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}\")\n"
      "FOREACH(example \${idas_bl_examples})\n"
      "  ADD_EXECUTABLE(\${example} \${example}.c)\n"
      "  TARGET_LINK_LIBRARIES(\${example} \${SUNDIALS_IDAS_BL_LIB} \${SUNDIALS_IDAS_LIB} \${SUNDIALS_NVECS_LIB} \${SUNDIALS_EXTRA_LIB})\n"
      "  TARGET_LINK_LIBRARIES(\${example} \${LAPACK_LIBRARIES})\n"
      "ENDFOREACH(example \${idas_bl_examples})\n"
      )

  ELSE(LAPACK_OK)

    FILE(WRITE ${PROJECT_BINARY_DIR}/examples/idas/serial/CMakeLists.txt
      "PROJECT(IDAS_serial_examples C)\n"
      "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
      "SET(idas_examples ${idas_examples})\n"
      "SET(SUNDIALS_INC_DIR\n"
      "  ${CMAKE_INSTALL_PREFIX}/include\n"
      "  CACHE STRING\n"
      "  \"Location of SUNDIALS header files\")\n"
      "INCLUDE_DIRECTORIES(\${SUNDIALS_INC_DIR})\n"
      "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
      "FIND_LIBRARY(SUNDIALS_IDAS_LIB\n"
      "  sundials_idas \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"IDAS library\")\n"
      "FIND_LIBRARY(SUNDIALS_NVECS_LIB\n"
      "  sundials_nvecserial \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"NVECSERIAL library\")\n"
      "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
      "FOREACH(example \${idas_examples})\n"
      "  ADD_EXECUTABLE(\${example} \${example}.c)\n"
      "  TARGET_LINK_LIBRARIES(\${example} \${SUNDIALS_IDAS_LIB} \${SUNDIALS_NVECS_LIB} \${SUNDIALS_EXTRA_LIB})\n"
      "ENDFOREACH(example \${idas_examples})\n"
      )
    
  ENDIF(LAPACK_OK)

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/idas/serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial 
    )

  # On UNIX-type platforms, we have the option of generating and instaling
  # a makefile for building the examples. This makefile can then be used
  # as a template for the user's own programs.

  IF(EXAMPLES_GENERATE_MAKEFILES)

    # generate Makefile and place it in the binary dir
    IF(LAPACK_OK)
      CONFIGURE_FILE(
        ${PROJECT_SOURCE_DIR}/examples/idas/serial/Makefile_ex_bl.in
        ${PROJECT_BINARY_DIR}/examples/idas/serial/Makefile_ex
        @ONLY
        )
    ELSE(LAPACK_OK)
      CONFIGURE_FILE(
        ${PROJECT_SOURCE_DIR}/examples/idas/serial/Makefile_ex.in
        ${PROJECT_BINARY_DIR}/examples/idas/serial/Makefile_ex
        @ONLY
        )
    ENDIF(LAPACK_OK)

    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/idas/serial/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial 
      RENAME Makefile
      )

  ENDIF(EXAMPLES_GENERATE_MAKEFILES)

ENDIF(EXAMPLES_INSTALL)
