# ---------------------------------------------------------------
# $Revision: 1.1 $
# $Date: 2007-01-29 17:40:01 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for IDA serial examples


# Add variable ida_examples with the names of the serial IDA examples

SET(ida_examples
  idabanx1
  idabanx2
  idadenx
  idakrydem_lin
  idakryx
  )

# Add variable ida_bl_examples with the names of the serial IDA examples
# that use Lapack

SET(ida_bl_examples
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDA_LIB sundials_ida_static)
  SET(IDA_BL_LIB sundials_ida_bl_static)
  SET(NVECS_LIB sundials_nvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDA_LIB sundials_ida_shared)
  SET(IDA_BL_LIB sundials_ida_bl_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Add the build and install targets for each IDA example

FOREACH(example ${ida_examples})
  # add build target for example
  ADD_EXECUTABLE(${example} ${example}.c)
  # link to the sundials libraries
  TARGET_LINK_LIBRARIES(${example} ${IDA_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
  # install example files
  INSTALL(FILES ${example}.c ${example}.out DESTINATION examples/ida/serial)
ENDFOREACH(example ${ida_examples})

# If Lapack support is enabled, add the build and install targets for
# the examples using Lapack

IF(LAPACK_OK)
  FOREACH(example ${ida_bl_examples})
    ADD_EXECUTABLE(${example} ${example}.c)
    TARGET_LINK_LIBRARIES(${example} ${IDA_BL_LIB} ${IDA_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
    TARGET_LINK_LIBRARIES(${example} ${LAPACK_LIBRARIES})
    INSTALL(FILES ${example}.c ${example}.out DESTINATION examples/ida/serial)
  ENDFOREACH(example ${ida_bl_examples})
ENDIF(LAPACK_OK)

# Regardless of the platform we're on, we will generate and install 
# CMakeLists.txt file for building the examples. This file  can then 
# be used as a template for the user's own programs.

# generate CMakelists.txt in the binary directory
IF(LAPACK_OK)

  FILE(WRITE ${PROJECT_BINARY_DIR}/examples/ida/serial/CMakeLists.txt
    "PROJECT(IDA_serial_examples C)\n"
    "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
    "SET(ida_examples ${ida_examples})\n"
    "SET(ida_bl_examples ${ida_bl_examples})\n"
    "SET(SUNDIALS_INC_DIR\n"
    "  ${CMAKE_INSTALL_PREFIX}/include\n"
    "  CACHE STRING\n"
    "  \"Location of SUNDIALS header files\")\n"
    "INCLUDE_DIRECTORIES(\${SUNDIALS_INC_DIR})\n"
    "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
    "FIND_LIBRARY(SUNDIALS_IDA_BL_LIB\n"
    "  sundials_ida_bl \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"IDA Lapack wrapper library\")\n"
    "FIND_LIBRARY(SUNDIALS_IDA_LIB\n"
    "  sundials_ida \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"IDA library\")\n"
    "FIND_LIBRARY(SUNDIALS_NVECS_LIB\n"
    "  sundials_nvecserial \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"NVECSERIAL library\")\n"
    "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
    "SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE STRING \"Lapack libraries\")\n"
    "FOREACH(example \${ida_examples})\n"
    "  ADD_EXECUTABLE(\${example} \${example}.c)\n"
    "  TARGET_LINK_LIBRARIES(\${example} \${SUNDIALS_IDA_LIB} \${SUNDIALS_NVECS_LIB} \${SUNDIALS_EXTRA_LIB})\n"
    "ENDFOREACH(example \${ida_examples})\n"
    "FOREACH(example \${ida_bl_examples})\n"
    "  ADD_EXECUTABLE(\${example} \${example}.c)\n"
    "  TARGET_LINK_LIBRARIES(\${example} \${SUNDIALS_IDA_BL_LIB} \${SUNDIALS_IDA_LIB} \${SUNDIALS_NVECS_LIB} \${SUNDIALS_EXTRA_LIB})\n"
    "  TARGET_LINK_LIBRARIES(\${example} \${LAPACK_LIBRARIES})\n"
    "ENDFOREACH(example \${ida_bl_examples})\n"
    )

ELSE(LAPACK_OK)

  FILE(WRITE ${PROJECT_BINARY_DIR}/examples/ida/serial/CMakeLists.txt
    "PROJECT(IDA_serial_examples C)\n"
    "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
    "SET(ida_examples ${ida_examples})\n"
    "SET(SUNDIALS_INC_DIR\n"
    "  ${CMAKE_INSTALL_PREFIX}/include\n"
    "  CACHE STRING\n"
    "  \"Location of SUNDIALS header files\")\n"
    "INCLUDE_DIRECTORIES(\${SUNDIALS_INC_DIR})\n"
    "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
    "FIND_LIBRARY(SUNDIALS_IDA_LIB\n"
    "  sundials_ida \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"IDA library\")\n"
    "FIND_LIBRARY(SUNDIALS_NVECS_LIB\n"
    "  sundials_nvecserial \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"NVECSERIAL library\")\n"
    "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
    "FOREACH(example \${ida_examples})\n"
    "  ADD_EXECUTABLE(\${example} \${example}.c)\n"
    "  TARGET_LINK_LIBRARIES(\${example} \${SUNDIALS_IDA_LIB} \${SUNDIALS_NVECS_LIB} \${SUNDIALS_EXTRA_LIB})\n"
    "ENDFOREACH(example \${ida_examples})\n"
    )
 
ENDIF(LAPACK_OK)

# install CMakelists.txt
INSTALL(
  FILES ${PROJECT_BINARY_DIR}/examples/ida/serial/CMakeLists.txt
  DESTINATION examples/ida/serial 
  )

# On UNIX-type platforms, we have the option of generating and instaling
# a makefile for building the examples. This makefile can then be used
# as a template for the user's own programs.

IF(GENERATE_EX_MAKEFILES)

  # generate Makefile and place it in the binary dir
  IF(LAPACK_OK)
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/ida/serial/Makefile_ex_bl.in
      ${PROJECT_BINARY_DIR}/examples/ida/serial/Makefile_ex
      @ONLY
      )
  ELSE(LAPACK_OK)
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/ida/serial/Makefile_ex.in
      ${PROJECT_BINARY_DIR}/examples/ida/serial/Makefile_ex
      @ONLY
      )
  ENDIF(LAPACK_OK)

  # install the configured Makefile_ex as Makefile
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/ida/serial/Makefile_ex 
    DESTINATION examples/ida/serial 
    RENAME Makefile
    )

ENDIF(GENERATE_EX_MAKEFILES)
