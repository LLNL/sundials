# ---------------------------------------------------------------
# $Revision: 1.8 $
# $Date: 2012-08-28 23:23:42 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for CVODES serial examples


# Add variable CVODES_examples with the names of the serial CVODES examples

SET(CVODES_examples
  cvsAdvDiff_ASAi_bnd
  cvsAdvDiff_bnd
  cvsAdvDiff_FSA_non
  cvsDirectDemo_ls
  cvsDiurnal_FSA_kry
  cvsDiurnal_kry_bp
  cvsDiurnal_kry
  cvsFoodWeb_ASAi_kry
  cvsFoodWeb_ASAp_kry
  cvsHessian_ASA_FSA
  cvsKrylovDemo_ls
  cvsKrylovDemo_prec
  cvsRoberts_ASAi_dns
  cvsRoberts_dns
  cvsRoberts_dns_uw
  cvsRoberts_FSA_dns
  )

# Add variable CVODES_examples_BL with the names of the serial CVODES examples
# that use Lapack

SET(CVODES_examples_BL
  cvsAdvDiff_bndL
  cvsRoberts_dnsL
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODES_LIB sundials_cvodes_static)
  SET(NVECS_LIB sundials_nvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODES_LIB sundials_cvodes_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Set-up linker flags and link libraries

SET(SUNDIALS_LIBS ${CVODES_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
IF(LAPACK_FOUND)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  SET(SUNDIALS_LIBS "${SUNDIALS_LIBS} ${LAPACK_LIBRARIES}")
ENDIF(LAPACK_FOUND)

# Add the build and install targets for each CVODES example

FOREACH(example ${CVODES_examples})
  ADD_EXECUTABLE(${example} ${example}.c)
  ADD_TEST(${example} ${example})
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
  FILE(GLOB example_out ${example}.out*)
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ${example_out} DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${CVODES_examples})

# If Lapack support is enabled, add the build and install targets for
# the examples using Lapack

IF(LAPACK_FOUND)
  FOREACH(example ${CVODES_examples_BL})
    ADD_EXECUTABLE(${example} ${example}.c)
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    FILE(GLOB example_out ${example}.out*)
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example_out} DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${CVODES_examples_BL})
ENDIF(LAPACK_FOUND)

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "CVODES")
  SET(SOLVER_LIB "sundials_cvodes")
  LIST2STRING(CVODES_examples EXAMPLES)
  IF(LAPACK_FOUND)
    LIST2STRING(CVODES_examples_BL EXAMPLES_BL)
  ELSE(LAPACK_FOUND)
    SET(EXAMPLES_BL "")
  ENDIF(LAPACK_FOUND)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvodes/serial/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/cvodes/serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial 
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvodes/serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/cvodes/serial/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial 
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
