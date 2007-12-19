# ---------------------------------------------------------------
# $Revision: 1.4 $
# $Date: 2007-12-19 20:34:00 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for KINSOL serial examples


# Add variable KINSOL_examples with the names of the serial KINSOL examples

SET(KINSOL_examples
  kinFerTron_dns
  kinFoodWeb_kry
  kinKrylovDemo_ls
  kinLaplace_bnd
  kinRoboKin_dns
  )

# Add variable KINSOL_examples_BL with the names of the serial KINSOL examples
# that use Lapack

SET(KINSOL_examples_BL
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(KINSOL_LIB sundials_kinsol_static)
  SET(NVECS_LIB sundials_nvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(KINSOL_LIB sundials_kinsol_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Add the build and install targets for each KINSOL example

FOREACH(example ${KINSOL_examples})
  # add build target for example
  ADD_EXECUTABLE(${example} ${example}.c)
  # link to the sundials libraries
  TARGET_LINK_LIBRARIES(${example} ${KINSOL_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
  IF(EXAMPLES_INSTALL)
    # install example files
    INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${KINSOL_examples})

# If Lapack support is enabled, add the build and install targets for
# the examples using Lapack

IF(LAPACK_OK)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  FOREACH(example ${KINSOL_examples_BL})
    ADD_EXECUTABLE(${example} ${example}.c)    
    TARGET_LINK_LIBRARIES(${example} ${KINSOL_LIB} ${KINSOL_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
    TARGET_LINK_LIBRARIES(${example} ${LAPACK_LIBRARIES})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${KINSOL_examples_BL})
ENDIF(LAPACK_OK)

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "KINSOL")
  SET(SOLVER_LIB "sundials_kinsol")
  LIST2STRING(KINSOL_examples EXAMPLES)
  IF(LAPACK_OK)
    LIST2STRING(KINSOL_examples_BL EXAMPLES_BL)
  ELSE(LAPACK_OK)
    SET(EXAMPLES_BL "")
  ENDIF(LAPACK_OK)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/kinsol/serial/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/kinsol/serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial 
    )

  # On UNIX-type platforms, we have the option of generating and instaling
  # a makefile for building the examples. This makefile can then be used
  # as a template for the user's own programs.

  IF(EXAMPLES_GENERATE_MAKEFILES)

    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/kinsol/serial/Makefile_ex
      @ONLY
      )

    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/kinsol/serial/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial 
      RENAME Makefile
      )

  ENDIF(EXAMPLES_GENERATE_MAKEFILES)

ENDIF(EXAMPLES_INSTALL)
