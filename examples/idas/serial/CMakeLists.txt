# ---------------------------------------------------------------
# $Revision: 1.12 $
# $Date: 2012-09-18 23:01:45 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for IDAS serial examples
#
# LEBanks 3/25/2013
#    Added SuperLU logic.


# Add variable IDAS_examples with the names of the serial IDAS examples

SET(IDAS_examples
  idasAkzoNob_ASAi_dns
  idasAkzoNob_dns
  idasFoodWeb_bnd
  idasHeat2D_bnd
  idasHeat2D_kry
  idasHessian_ASA_FSA
  idasKrylovDemo_ls
  idasRoberts_ASAi_dns
  idasRoberts_dns
  idasRoberts_FSA_dns
  idasSlCrank_dns
  idasSlCrank_FSA_dns
  )

# Add variable IDAS_examples_BL with the names of the serial IDAS examples
# that use Lapack

SET(IDAS_examples_BL
  )
  
SET(IDAS_examples_SUPERLU
  idasRoberts_sps
  )

SET(IDAS_examples_KLU
  idasRoberts_klu
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDAS_LIB sundials_idas_static)
  SET(NVECS_LIB sundials_nvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDAS_LIB sundials_idas_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Set-up linker flags and link libraries

SET(SUNDIALS_LIBS ${IDAS_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
IF(LAPACK_FOUND)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  SET(SUNDIALS_LIBS "${SUNDIALS_LIBS} ${LAPACK_LIBRARIES}")
ENDIF(LAPACK_FOUND)

IF(SUPERLU_FOUND)
  SET(SUNDIALS_LIBS "${SUNDIALS_LIBS} ${SUPERLU_LIBRARIES}")
ENDIF(SUPERLU_FOUND)

# Add the build and install targets for each IDAS example

FOREACH(example ${IDAS_examples})
  ADD_EXECUTABLE(${example} ${example}.c)
  ADD_TEST(${example} ${TESTRUNNER} ${example})
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${IDAS_examples})

# If Lapack support is enabled, add the build and install targets for
# the examples using Lapack
IF(LAPACK_FOUND)
  FOREACH(example ${IDAS_examples_BL})
    ADD_EXECUTABLE(${example} ${example}.c)
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${IDAS_examples_BL})
ENDIF(LAPACK_FOUND)

# If SuperLU support is enabled, add the build and install targets for
# the examples using SuperLU
IF(SUPERLU_FOUND)
  FOREACH(example ${IDAS_examples_SUPERLU})
    ADD_EXECUTABLE(${example} ${example}.c)
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${IDAS_examples_SUPERLU})
ENDIF(SUPERLU_FOUND)

# If KLU support is enabled, add the build and install targets for
# the examples using KLU
IF(KLU_FOUND)
  SET(SUNDIALS_LIBS "${SUNDIALS_LIBS} ${KLU_LIBRARIES}")

  FOREACH(example ${IDAS_examples_KLU})
    ADD_EXECUTABLE(${example} ${example}.c)
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${IDAS_examples_KLU})
ENDIF(KLU_FOUND)

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "IDAS")
  SET(SOLVER_LIB "sundials_idas")
  LIST2STRING(IDAS_examples EXAMPLES)

  IF(LAPACK_FOUND)
    LIST2STRING(IDAS_examples_BL EXAMPLES_BL)
  ELSE(LAPACK_FOUND)
    SET(EXAMPLES_BL "")
  ENDIF(LAPACK_FOUND)

  # SGS TODO what are these used for?  Don't seem to be reference anywhere
  IF(SUPERLU_FOUND)
    LIST2STRING(IDAS_examples_SUPERLU EXAMPLES_SUPERLU)
  ELSE(SUPERLU_FOUND)
    SET(EXAMPLES_SUPERLU "")
  ENDIF(SUPERLU_FOUND)

  IF(KLU_FOUND)
    LIST2STRING(IDAS_examples_KLU EXAMPLES_KLU)
  ELSE(KLU_FOUND)
    SET(EXAMPLES_KLU "")
  ENDIF(KLU_FOUND)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/idas/serial/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/idas/serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial 
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/idas/serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/idas/serial/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial 
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
