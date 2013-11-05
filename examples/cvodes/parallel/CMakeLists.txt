# ---------------------------------------------------------------
# $Revision: 1.6 $
# $Date: 2012-09-18 23:01:38 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the CVODES parallel examples


# Add variable CVODES_examples with the names of the parallel CVODES examples
# "name\;nodes\;tasks"
SET(CVODES_examples
  "cvsAdvDiff_ASAp_non_p\;2\;2"
  "cvsAdvDiff_FSA_non_p\;2\;2"
  "cvsAdvDiff_non_p\;2\;2"
  "cvsAtmDisp_ASAi_kry_bbd_p\;2\;4"
  "cvsDiurnal_FSA_kry_p\;2\;4"
  "cvsDiurnal_kry_bbd_p\;2\;4"
  "cvsDiurnal_kry_p\;2\;4"
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
  SET(CVODES_LIB sundials_cvodes_static)
  SET(NVECP_LIB sundials_nvecparallel_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODES_LIB sundials_cvodes_shared)
  SET(NVECP_LIB sundials_nvecparallel_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Set-up linker flags and link libraries

SET(SUNDIALS_LIBS ${CVODES_LIB} ${NVECP_LIB} ${EXTRA_LINK_LIBS})
IF(LAPACK_FOUND)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  SET(SUNDIALS_LIBS "${SUNDIALS_LIBS} ${LAPACK_LIBRARIES}")
ENDIF(LAPACK_FOUND)

# Add the build and install targets for each CVODES example

FOREACH(example_tuple ${CVODES_examples})
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_nodes)
  list(GET example_tuple 2 number_of_tasks)
  # first item is example
  ADD_EXECUTABLE(${example} ${example}.c)
  IF(MPI_RUN_COMMAND MATCHES "srun")
      ADD_TEST(${example} ${TESTRUNNER} ${example} "srun -N${number_of_nodes} -n${number_of_tasks} -ppdebug")
  ELSE(MPI_RUN_COMMAND MATCHES "srun")
      ADD_TEST(${example} ${TESTRUNNER} ${example} "mpirun -n ${number_of_tasks}")
  ENDIF(MPI_RUN_COMMAND MATCHES "srun")
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
  IF(NOT MPI_MPICC)
    TARGET_LINK_LIBRARIES(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  ENDIF(NOT MPI_MPICC)
  FILE(GLOB example_out ${example}.out*)
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ${example_out} DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/parallel)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example_tuple ${CVODES_examples})

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/parallel)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "CVODES")
  SET(SOLVER_LIB "sundials_cvodes")
  LIST2STRING(CVODES_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_parallel_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvodes/parallel/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/cvodes/parallel/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/parallel 
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_parallel_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvodes/parallel/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/cvodes/parallel/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/parallel
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
