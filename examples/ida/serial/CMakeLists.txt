# ---------------------------------------------------------------
# $Revision: 1.8 $
# $Date: 2012-09-18 23:01:43 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for IDA serial examples


# Add variable IDA_examples with the names of the serial IDA examples

SET(IDA_examples
  idaFoodWeb_bnd
  idaHeat2D_bnd
  idaHeat2D_kry
  idaKrylovDemo_ls
  idaRoberts_dns
  idaSlCrank_dns
  )

# Add variable IDA_examples_BL with the names of the serial IDA examples
# that use Lapack

SET(IDA_examples_BL
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDA_LIB sundials_ida_static)
  SET(NVECS_LIB sundials_nvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDA_LIB sundials_ida_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Set-up linker flags and link libraries

SET(SUNDIALS_LIBS ${IDA_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
IF(LAPACK_FOUND)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  SET(SUNDIALS_LIBS "${SUNDIALS_LIBS} ${LAPACK_LIBRARIES}")
ENDIF(LAPACK_FOUND)

# Add the build and install targets for each IDA example

FOREACH(example ${IDA_examples})
  ADD_EXECUTABLE(${example} ${example}.c)
  ADD_TEST(${example} ${TESTRUNNER} ${example})
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/serial)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${IDA_examples})

# If Lapack support is enabled, add the build and install targets for
# the examples using Lapack

IF(LAPACK_FOUND)
  FOREACH(example ${IDA_examples_BL})
    ADD_EXECUTABLE(${example} ${example}.c)
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${IDA_examples_BL})
ENDIF(LAPACK_FOUND)

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "IDA")
  SET(SOLVER_LIB "sundials_ida")
  LIST2STRING(IDA_examples EXAMPLES)
  IF(LAPACK_FOUND)
    LIST2STRING(IDA_examples_BL EXAMPLES_BL)
  ELSE(LAPACK_FOUND)
    SET(EXAMPLES_BL "")
  ENDIF(LAPACK_FOUND)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/ida/serial/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/ida/serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/serial 
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.
  
  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/ida/serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/ida/serial/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/serial 
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
