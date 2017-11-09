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
# CMakeLists.txt file for the FIDA parallel examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
SET(FIDA_examples
  "fidaHeat2D_kry_bbd_p\;1\;4\;develop"
  )

# Check whether we use MPI compiler scripts.
# If yes, then change the Fortran compiler to the MPIF77 script.
# If not, then add the MPI include directory for MPI headers.
IF(MPI_MPIF77)
  # use MPI_MPIF77 as the compiler
  SET(CMAKE_Fortran_COMPILER ${MPI_MPIF77})
ELSE()
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDA_LIB sundials_ida_static)
  SET(NVECP_LIB sundials_nvecparallel_static)
  SET(FNVECP_LIB sundials_fnvecparallel_static)
ELSE()
  SET(IDA_LIB sundials_ida_shared)
  SET(NVECP_LIB sundials_nvecparallel_shared)
  SET(FNVECP_LIB sundials_fnvecparallel_shared)
ENDIF()

# Only static FCMIX libraries are available
SET(FIDA_LIB sundials_fida_static)

# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${FIDA_LIB} ${IDA_LIB} ${FNVECP_LIB} ${NVECP_LIB} ${EXTRA_LINK_LIBS})


# Add the build and install targets for each example
FOREACH(example_tuple ${FIDA_examples})

  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_nodes)
  list(GET example_tuple 2 number_of_tasks)
  list(GET example_tuple 3 example_type)

  # example source files
  ADD_EXECUTABLE(${example} ${example}.f)

  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  SUNDIALS_ADD_TEST(${example} ${example}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})

  IF(NOT MPI_MPIF77)
    TARGET_LINK_LIBRARIES(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  ENDIF()

  # install example source and out files
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.f ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/fcmix_parallel)
  ENDIF()

ENDFOREACH(example_tuple ${FIDA_examples})


# create Makfile and CMakeLists.txt for examples
IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/fcmix_parallel)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "IDA")
  SET(SOLVER_LIB "sundials_ida")
  SET(SOLVER_FLIB "sundials_fida")

  EXAMPLES2STRING(FIDA_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_parallel_F77_ex.in
    ${PROJECT_BINARY_DIR}/examples/ida/fcmix_parallel/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/ida/fcmix_parallel/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/fcmix_parallel
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_parallel_F77_ex.in
      ${PROJECT_BINARY_DIR}/examples/ida/fcmix_parallel/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/ida/fcmix_parallel/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/fcmix_parallel
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
