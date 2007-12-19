# ---------------------------------------------------------------
# $Revision: 1.3 $
# $Date: 2007-12-19 20:26:43 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for CPODES serial examples


# Add variable CPODES_examples with the names of the serial CPODES examples

SET(CPODES_examples
  cpsAdvDiff_bnd
  cpsAdvDiff_non
  cpsNewtCrd_dns
  cpsPend_dns
  cpsRoberts_dns
  cpsSlCrank_dns
  cpsVanDPol_non
  )

# Add variable CPODES_examples_BL with the names of the serial CPODES examples
# that use Lapack

SET(CPODES_examples_BL
  cpsAdvDiff_bndL
  cpsAkzoNob_dnsL
  cpsPend_dnsL
  cpsRoberts_dnsL
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CPODES_LIB sundials_cpodes_static)
  SET(NVECS_LIB sundials_nvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CPODES_LIB sundials_cpodes_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Add the build and install targets for each CPODES example

FOREACH(example ${CPODES_examples})
  # add build target for example
  ADD_EXECUTABLE(${example} ${example}.c)
  # link to the sundials libraries
  TARGET_LINK_LIBRARIES(${example} ${CPODES_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
  # find all sample output files for this example
  FILE(GLOB example_out ${example}.out*)
  IF(EXAMPLES_INSTALL)
    # install example files
    INSTALL(FILES ${example}.c ${example_out} DESTINATION ${EXAMPLES_INSTALL_PATH}/cpodes/serial)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${CPODES_examples})

# If Lapack support is enabled, add the build and install targets for
# the examples using Lapack

IF(LAPACK_OK)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  FOREACH(example ${CPODES_examples_BL})
    ADD_EXECUTABLE(${example} ${example}.c)
    TARGET_LINK_LIBRARIES(${example} ${CPODES_LIB} ${CPODES_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
    TARGET_LINK_LIBRARIES(${example} ${LAPACK_LIBRARIES})
    FILE(GLOB example_out ${example}.out*)
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example_out} DESTINATION ${EXAMPLES_INSTALL_PATH}/cpodes/serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${CPODES_examples_BL})
ENDIF(LAPACK_OK)

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cpodes/serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "CPODES")
  SET(SOLVER_LIB "sundials_cpodes")
  LIST2STRING(CPODES_examples EXAMPLES)
  IF(LAPACK_OK)
    LIST2STRING(CPODES_examples_BL EXAMPLES_BL)
  ELSE(LAPACK_OK)
    SET(EXAMPLES_BL "")
  ENDIF(LAPACK_OK)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cpodes/serial/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/cpodes/serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cpodes/serial 
    )

  # On UNIX-type platforms, we have the option of generating and instaling
  # a makefile for building the examples. This makefile can then be used
  # as a template for the user's own programs.

  IF(EXAMPLES_GENERATE_MAKEFILES)

    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cpodes/serial/Makefile_ex
      @ONLY
      )

    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/cpodes/serial/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cpodes/serial 
      RENAME Makefile
      )

  ENDIF(EXAMPLES_GENERATE_MAKEFILES)

ENDIF(EXAMPLES_INSTALL)
