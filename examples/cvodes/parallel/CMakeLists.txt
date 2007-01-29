# ---------------------------------------------------------------
# $Revision: 1.1 $
# $Date: 2007-01-29 17:40:01 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the CVODES parallel examples


# Add variable cvodes_examples with the names of the parallel CVODES examples

SET(cvodes_examples
  cvsadjkryx_p
  cvsadjnonx_p
  cvsfwdkryx_p
  cvsfwdnonx_p
  cvskryx_bbd_p
  cvskryx_p
  cvsnonx_p
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
  SET(CVODES_LIB sundials_cvodes_static)
  SET(NVECP_LIB sundials_nvecparallel_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODES_LIB sundials_cvodes_shared)
  SET(NVECP_LIB sundials_nvecparallel_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Add the build and install targets for each CVODES example

FOREACH(example ${cvodes_examples})
  # add build target for example
  ADD_EXECUTABLE(${example} ${example}.c)
  # link to the sundials libraries
  TARGET_LINK_LIBRARIES(${example} ${CVODES_LIB} ${NVECP_LIB} ${EXTRA_LINK_LIBS})
  # if needed, link to the MPI libraries
  IF(NOT MPI_USE_MPISCRIPTS)
    TARGET_LINK_LIBRARIES(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  ENDIF(NOT MPI_USE_MPISCRIPTS)
  # find all sample output files for this example
  FILE(GLOB example_out ${example}.out*)
  # install example files
  INSTALL(FILES ${example}.c ${example_out} DESTINATION examples/cvodes/parallel)
ENDFOREACH(example ${cvodes_examples})

# Regardless of the platform we're on, we will generate and install 
# CMakeLists.txt file for building the examples. This file  can then 
# be used as a template for the user's own programs.

# generate CMakelists.txt in the binary directory
IF(MPI_USE_MPISCRIPTS)
  FILE(WRITE ${PROJECT_BINARY_DIR}/examples/cvodes/parallel/CMakeLists.txt
    "PROJECT(CVODES_parallel_examples C)\n"
    "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
    "SET(cvodes_examples ${cvodes_examples})\n"
    "SET(MPI_MPICC ${MPI_MPICC} CACHE STRING \"MPICC compiler script\")\n"
    "SET(CMAKE_C_COMPILER \${MPI_MPICC})\n"
    "SET(SUNDIALS_INC_DIR\n"
    "  ${CMAKE_INSTALL_PREFIX}/include\n"
    "  CACHE STRING\n"
    "  \"Location of SUNDIALS header files\")\n"
    "INCLUDE_DIRECTORIES(\${SUNDIALS_INC_DIR})\n"
    "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
    "FIND_LIBRARY(SUNDIALS_CVODES_LIB\n"
    "  sundials_cvodes \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"CVODES library\")\n"
    "FIND_LIBRARY(SUNDIALS_NVECP_LIB\n"
    "  sundials_nvecparallel \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"NVECPARALLEL library\")\n"
    "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
    "FOREACH(example \${cvodes_examples})\n"
    "  ADD_EXECUTABLE(\${example} \${example}.c)\n"
    "  TARGET_LINK_LIBRARIES(\${example} \${SUNDIALS_CVODES_LIB} \${SUNDIALS_NVECP_LIB} \${SUNDIALS_EXTRA_LIB})\n"
    "ENDFOREACH(example \${cvodes_examples})\n"
    )
ELSE(MPI_USE_MPISCRIPTS)
  FILE(WRITE ${PROJECT_BINARY_DIR}/examples/cvodes/parallel/CMakeLists.txt
    "PROJECT(CVODES_parallel_examples C)\n"
    "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
    "SET(cvodes_examples ${cvodes_examples})\n"
    "SET(SUNDIALS_INC_DIR\n"
    "  ${CMAKE_INSTALL_PREFIX}/include\n"
    "  CACHE STRING\n"
    "  \"Location of SUNDIALS header files\")\n"
    "INCLUDE_DIRECTORIES(\${SUNDIALS_INC_DIR})\n"
    "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
    "FIND_LIBRARY(SUNDIALS_CVODES_LIB\n"
    "  sundials_cvodes \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"CVODES library\")\n"
    "FIND_LIBRARY(SUNDIALS_NVECP_LIB\n"
    "  sundials_nvecparallel \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"NVECPARALLEL library\")\n"
    "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
    "SET(MPI_INC_DIR ${MPI_INCLUDE_PATH} CACHE STRING \"Location of MPI header files\")\n"
    "INCLUDE_DIRECTORIES(\${MPI_INC_DIR})\n"
    "SET(MPI_LIB ${MPI_LIBRARY} CACHE STRING \"Main MPI library\")\n"
    "SET(MPI_EXTRA_LIB ${MPI_EXTRA_LIBRARIES} CACHE STRING \"Additional MPI libraries\")\n"
    "FOREACH(example \${cvodes_examples})\n"
    "  ADD_EXECUTABLE(\${example} \${example}.c)\n"
    "  TARGET_LINK_LIBRARIES(\${example} \${SUNDIALS_CVODES_LIB} \${SUNDIALS_NVECP_LIB} \${SUNDIALS_EXTRA_LIB} \${MPI_LIB} \${MPI_EXTRA_LIB})\n"
    "ENDFOREACH(example \${cvodes_examples})\n"
    )
ENDIF(MPI_USE_MPISCRIPTS)

# install CMakelists.txt
INSTALL(
  FILES ${PROJECT_BINARY_DIR}/examples/cvodes/parallel/CMakeLists.txt
  DESTINATION examples/cvodes/parallel 
  )

# On UNIX-type platforms, we have the option of generating and instaling
# a makefile for building the examples. This makefile can then be used
# as a template for the user's own programs.

IF(GENERATE_EX_MAKEFILES)
  # generate Makefile and place it in the binary dir
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/cvodes/parallel/Makefile_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvodes/parallel/Makefile_ex
    @ONLY
    )
  # install the configured Makefile_ex as Makefile
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/cvodes/parallel/Makefile_ex
    DESTINATION examples/cvodes/parallel
    RENAME Makefile
    )
ENDIF(GENERATE_EX_MAKEFILES)
