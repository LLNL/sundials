# ---------------------------------------------------------------
# $Revision: 4942 $
# $Date: 2016-09-21 17:50:22 -0700 (Wed, 21 Sep 2016) $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the CVODES parallel examples


# Add variable CVODES_examples with the names of the parallel CVODES examples
# "name\;nodes\;tasks"
SET(CVODES_examples
  "cvsAdvDiff_ASAp_non_p\;2\;2"
  "cvsAdvDiff_FSA_non_p\;2\;2\;-sensi stg t"
  "cvsAdvDiff_FSA_non_p\;2\;2\;-sensi sim t"
  "cvsAdvDiff_non_p\;2\;2"
  "cvsAtmDisp_ASAi_kry_bbd_p\;2\;8"
  "cvsDiurnal_FSA_kry_p\;2\;4\;-sensi stg t"
  "cvsDiurnal_FSA_kry_p\;2\;4\;-sensi sim t"
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
  LIST(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})
ENDIF(LAPACK_FOUND)

IF(KLU_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${KLU_LIBRARIES})
ENDIF(KLU_FOUND)

IF(SUPERLUMT_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${SUPERLUMT_LIBRARIES})
ENDIF(SUPERLUMT_FOUND)

# Add the build and install targets for each CVODES example

# SGS is there a way to query CMAKE variable to get added executables 
# rather than keeping our own list?
SET(ADDED_EXECUTABLES "")
FOREACH(example_tuple ${CVODES_examples})
  # first item is example
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_nodes)
  list(GET example_tuple 2 number_of_tasks)

  # Only need to add the executable once
  LIST(FIND ADDED_EXECUTABLES ${example} index)
  IF(index EQUAL -1)
    LIST(APPEND ADDED_EXECUTABLES ${example})
    ADD_EXECUTABLE(${example} ${example}.c)
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
  ENDIF(index EQUAL -1)

  # Optional 4th element is test arguments
  LIST(LENGTH example_tuple n)
  IF(n EQUAL 4)
    LIST(GET example_tuple 3 test_args)
    STRING(REGEX REPLACE " " "_" test_name ${example}_${test_args})
    SUNDIALS_ADD_TEST(${test_name} ${example} MPI_NPROCS ${number_of_tasks} TEST_ARGS ${test_args})
  ELSE()
    SUNDIALS_ADD_TEST(${example} ${example} MPI_NPROCS ${number_of_tasks})
  ENDIF()

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
  FOREACH(example_tuple ${CVODES_examples})
  	list(GET example_tuple 0 example)
  	LIST2STRING(example EXAMPLES)
  ENDFOREACH(example_tuple ${CVODES_examples})

  STRING (REPLACE ";" " " TMP_STR ${EXAMPLES})
  SET(EXAMPLES ${TMP_STR})
  list(REMOVE_DUPLICATES EXAMPLES)

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
