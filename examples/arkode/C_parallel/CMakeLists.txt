# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# LLNS/SMU Copyright Start
# Copyright (c) 2015, Southern Methodist University and 
# Lawrence Livermore National Security
#
# This work was performed under the auspices of the U.S. Department 
# of Energy by Southern Methodist University and Lawrence Livermore 
# National Laboratory under Contract DE-AC52-07NA27344.
# Produced at Southern Methodist University and the Lawrence 
# Livermore National Laboratory.
#
# All rights reserved.
# For details, see the LICENSE file.
# LLNS/SMU Copyright End
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
# -----------------------------------------------------------------
# CMakeLists.txt file for ARKODE parallel C examples


# Add variable ARKODE_examples with the names of the parallel ARKODE examples
# "name\;nodes\;tasks"
SET(ARKODE_examples
  "ark_diurnal_kry_p\;1\;4"
  "ark_diurnal_kry_bbd_p\;1\;4"
  )

# Add variable ARKODE_extras with the names of auxiliary files to install
SET(ARKODE_extras
  # none yet
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
  SET(ARKODE_LIB sundials_arkode_static)
  SET(NVECP_LIB sundials_nvecparallel_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(ARKODE_LIB sundials_arkode_shared)
  SET(NVECP_LIB sundials_nvecparallel_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${ARKODE_LIB} ${NVECP_LIB} ${EXTRA_LINK_LIBS})
IF(LAPACK_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})
ENDIF(LAPACK_FOUND)

IF(KLU_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${KLU_LIBRARIES})
ENDIF(KLU_FOUND)

IF(SUPERLUMT_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${SUPERLUMT_LIBRARIES})
ENDIF(SUPERLUMT_FOUND)

# Add the build and install targets for each ARKODE example
FOREACH(example_tuple ${ARKODE_examples})
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_nodes)
  list(GET example_tuple 2 number_of_tasks)
  # first item is example
  ADD_EXECUTABLE(${example} ${example}.c)
  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
  SUNDIALS_ADD_TEST(${example} ${example} MPI_NPROCS ${number_of_tasks})
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
  IF(NOT MPI_MPICC)
    TARGET_LINK_LIBRARIES(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  ENDIF(NOT MPI_MPICC)
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_parallel)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example_tuple ${ARKODE_examples})

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_parallel)

  # Install the extra files
  FOREACH(extrafile ${ARKODE_extras})
    INSTALL(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_parallel)
  ENDFOREACH(extrafile ${ARKODE_extras})

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "ARKODE")
  SET(SOLVER_LIB "sundials_arkode")
  FOREACH(example_tuple ${ARKODE_examples})
  	list(GET example_tuple 0 example)
  	LIST2STRING(example EXAMPLES)
  ENDFOREACH(example_tuple ${ARKODE_examples})


  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_parallel_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/C_parallel/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/arkode/C_parallel/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_parallel 
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_parallel_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/C_parallel/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/arkode/C_parallel/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_parallel
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
