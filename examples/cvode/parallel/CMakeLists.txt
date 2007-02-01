# ---------------------------------------------------------------
# $Revision: 1.2 $
# $Date: 2007-02-01 16:42:08 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for CVODE parallel examples


# Add variable cvode_examples with the names of the parallel CVODE examples

SET(cvode_examples
  cvnonx_p
  cvkryx_p
  cvkryx_bbd_p
  )

# Check whether we use MPI compiler scripts.
# If yes, then change the C compiler to the MPICC script.
# If not, then add the MPI include directory for MPI headers.

IF(MPI_USE_MPISCRIPTS)
  # use MPI_MPICC as the compiler
  SET(CMAKE_C_COMPILER ${MPI_MPICC})
ELSE(MPI_USE_MPISCRIPTS)
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(MPI_USE_MPISCRIPTS)

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODE_LIB sundials_cvode_static)
  SET(NVECP_LIB sundials_nvecparallel_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODE_LIB sundials_cvode_shared)
  SET(NVECP_LIB sundials_nvecparallel_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Add the build and install targets for each CVODE example

FOREACH(example ${cvode_examples})
  # add build target for example
  ADD_EXECUTABLE(${example} ${example}.c)
  # link to the sundials libraries
  TARGET_LINK_LIBRARIES(${example} ${CVODE_LIB} ${NVECP_LIB} ${EXTRA_LINK_LIBS})
  # if needed, link to the MPI libraries
  IF(NOT MPI_USE_MPISCRIPTS)
    TARGET_LINK_LIBRARIES(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  ENDIF(NOT MPI_USE_MPISCRIPTS)
  IF(EXAMPLES_INSTALL)
    # install example files
    INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/parallel)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${cvode_examples})

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/parallel)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  IF(MPI_USE_MPISCRIPTS)
    FILE(WRITE ${PROJECT_BINARY_DIR}/examples/cvode/parallel/CMakeLists.txt
      "PROJECT(CVODE_parallel_examples C)\n"
      "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
      "SET(cvode_examples ${cvode_examples})\n"
      "SET(MPI_MPICC ${MPI_MPICC} CACHE STRING \"MPICC compiler script\")\n"
      "SET(CMAKE_C_COMPILER \${MPI_MPICC})\n"
      "SET(SUNDIALS_INC_DIR\n"
      "  ${CMAKE_INSTALL_PREFIX}/include\n"
      "  CACHE STRING\n"
      "  \"Location of SUNDIALS header files\")\n"
      "INCLUDE_DIRECTORIES(\${SUNDIALS_INC_DIR})\n"
      "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
      "FIND_LIBRARY(SUNDIALS_CVODE_LIB\n"
      "  sundials_cvode \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"CVODE library\")\n"
      "FIND_LIBRARY(SUNDIALS_NVECP_LIB\n"
      "  sundials_nvecparallel \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"NVECPARALLEL library\")\n"
      "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
      "FOREACH(example \${cvode_examples})\n"
      "  ADD_EXECUTABLE(\${example} \${example}.c)\n"
      "  TARGET_LINK_LIBRARIES(\${example} \${SUNDIALS_CVODE_LIB} \${SUNDIALS_NVECP_LIB} \${SUNDIALS_EXTRA_LIB})\n"
      "ENDFOREACH(example \${cvode_examples})\n"
      )
  ELSE(MPI_USE_MPISCRIPTS)
    FILE(WRITE ${PROJECT_BINARY_DIR}/examples/cvode/parallel/CMakeLists.txt
      "PROJECT(CVODE_parallel_examples C)\n"
      "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
      "SET(cvode_examples ${cvode_examples})\n"
      "SET(SUNDIALS_INC_DIR\n"
      "  ${CMAKE_INSTALL_PREFIX}/include\n"
      "  CACHE STRING\n"
      "  \"Location of SUNDIALS header files\")\n"
      "INCLUDE_DIRECTORIES(\${SUNDIALS_INC_DIR})\n"
      "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
      "FIND_LIBRARY(SUNDIALS_CVODE_LIB\n"
      "  sundials_cvode \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"CVODE library\")\n"
      "FIND_LIBRARY(SUNDIALS_NVECP_LIB\n"
      "  sundials_nvecparallel \${SUNDIALS_LIB_DIR}\n"
      "  DOC \"NVECPARALLEL library\")\n"
      "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
      "SET(MPI_INC_DIR ${MPI_INCLUDE_PATH} CACHE STRING \"Location of MPI header files\")\n"
      "INCLUDE_DIRECTORIES(\${MPI_INC_DIR})\n"
      "SET(MPI_LIB ${MPI_LIBRARY} CACHE STRING \"Main MPI library\")\n"
      "SET(MPI_EXTRA_LIB ${MPI_EXTRA_LIBRARIES} CACHE STRING \"Additional MPI libraries\")\n"
      "FOREACH(example \${cvode_examples})\n"
      "  ADD_EXECUTABLE(\${example} \${example}.c)\n"
      "  TARGET_LINK_LIBRARIES(\${example} \${SUNDIALS_CVODE_LIB} \${SUNDIALS_NVECP_LIB} \${SUNDIALS_EXTRA_LIB} \${MPI_LIB} \${MPI_EXTRA_LIB})\n"
      "ENDFOREACH(example \${cvode_examples})\n"
      )
  ENDIF(MPI_USE_MPISCRIPTS)

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/cvode/parallel/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/parallel 
    )

  # On UNIX-type platforms, we have the option of generating and instaling
  # a makefile for building the examples. This makefile can then be used
  # as a template for the user's own programs.

  IF(EXAMPLES_GENERATE_MAKEFILES)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/cvode/parallel/Makefile_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/parallel/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/parallel/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/parallel
      RENAME Makefile
      )
  ENDIF(EXAMPLES_GENERATE_MAKEFILES)

ENDIF(EXAMPLES_INSTALL)
