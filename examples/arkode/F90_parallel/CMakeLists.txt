# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
# -----------------------------------------------------------------
# CMakeLists.txt file for the ARKODE parallel F90 examples


# Add variable FARKODE_examples with the names of the parallel ARKODE examples

SET(FARKODE_examples
  fark_heat2D
  )

# Add variable ARKODE_extras with the names of auxiliary files to install

SET(ARKODE_extras
  plot_heat2D.py
  )

# Check whether we use MPI compiler scripts.
# If yes, then change the Fortran compiler to the MPIF90 script.
# If not, then add the MPI include directory for MPI headers.

IF(MPI_MPIF90 )
  # use MPI_MPIF90 as the compiler
  SET(CMAKE_Fortran_COMPILER ${MPI_MPIF90})
ELSE(MPI_MPIF90)
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(MPI_MPIF90)

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(ARKODE_LIB sundials_arkode_static)
  SET(NVECP_LIB sundials_nvecparallel_static)
  SET(FNVECP_LIB sundials_fnvecparallel_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(ARKODE_LIB sundials_arkode_shared)
  SET(NVECP_LIB sundials_nvecparallel_shared)
  SET(FNVECP_LIB sundials_fnvecparallel_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Only static FCMIX libraries are available

SET(FARKODE_LIB sundials_farkode_static)

# Set-up linker flags and link libraries

SET(SUNDIALS_LIBS ${FARKODE_LIB} ${ARKODE_LIB} ${FNVECP_LIB} ${NVECP_LIB} ${EXTRA_LINK_LIBS})
IF(LAPACK_FOUND)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  SET(SUNDIALS_LIBS "${SUNDIALS_LIBS} ${LAPACK_LIBRARIES}")
ENDIF(LAPACK_FOUND)

# Add the build and install targets for each ARKODE example

FOREACH(example ${FARKODE_examples})
  ADD_EXECUTABLE(${example} ${example}.f90)
  ADD_TEST(${example} ${TESTRUNNER} ${example})
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
  IF(NOT MPI_MPIF90)
    TARGET_LINK_LIBRARIES(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  ENDIF(NOT MPI_MPIF90)
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.f90 ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_parallel)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${FARKODE_examples})

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_parallel)

  # Install the extra files
  FOREACH(extrafile ${ARKODE_extras})
    INSTALL(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_parallel)
  ENDFOREACH(extrafile ${ARKODE_extras})

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "ARKODE")
  SET(SOLVER_LIB "sundials_arkode")
  SET(SOLVER_FLIB "sundials_farkode")
  LIST2STRING(FARKODE_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_parallel_F90_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/F90_parallel/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/arkode/F90_parallel/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_parallel 
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_parallel_F90_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/F90_parallel/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/arkode/F90_parallel/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_parallel 
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
