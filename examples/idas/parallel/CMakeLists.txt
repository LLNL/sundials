# ---------------------------------------------------------------
# $Revision: 1.7 $
# $Date: 2012-09-18 23:01:44 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for IDAS parallel examples


# Add variable IDAS_examples with the names of the parallel IDAS examples

SET(IDAS_examples
  idasBruss_ASAp_kry_bbd_p
  idasBruss_FSA_kry_bbd_p
  idasBruss_kry_bbd_p
  idasFoodWeb_kry_bbd_p
  idasFoodWeb_kry_p
  idasHeat2D_FSA_kry_bbd_p
  idasHeat2D_kry_bbd_p
  idasHeat2D_kry_p
  )

# Check whether we use MPI compiler scripts.
# If yes, then change the C compiler to the MPICC script.
# If not, then add the MPI include directory for MPI headers.

IF(MPI_MPICC)
  # use MPI_MPICC as the compiler
  SET(CMAKE_C_COMPILER ${MPI_MPICC})
ELSE(MPI_MPICC)
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(MPI_MPICC)

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDAS_LIB sundials_idas_static)
  SET(NVECP_LIB sundials_nvecparallel_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDAS_LIB sundials_idas_shared)
  SET(NVECP_LIB sundials_nvecparallel_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Set-up linker flags and link libraries

SET(SUNDIALS_LIBS ${IDAS_LIB} ${NVECP_LIB} ${EXTRA_LINK_LIBS})
IF(LAPACK_FOUND)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  SET(SUNDIALS_LIBS "${SUNDIALS_LIBS} ${LAPACK_LIBRARIES}")
ENDIF(LAPACK_FOUND)

# Add the build and install targets for each IDAS example

FOREACH(example ${IDAS_examples})
  ADD_EXECUTABLE(${example} ${example}.c)
  ADD_TEST(${example} ${TESTRUNNER} ${example})
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
  IF(NOT MPI_MPICC)
    TARGET_LINK_LIBRARIES(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  ENDIF(NOT MPI_MPICC)
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/parallel)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${IDAS_examples})

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/parallel)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "IDAS")
  SET(SOLVER_LIB "sundials_idas")
  LIST2STRING(IDAS_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_parallel_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/idas/parallel/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/idas/parallel/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/parallel 
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_parallel_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/idas/parallel/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/idas/parallel/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/parallel
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
