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
# CMakeLists.txt file for the FKINSOL parallel examples


# Add variable kinsol_examples with the names of the parallel KINSOL examples

SET(fkinsol_examples
  fkinkryx_p
  )

# Check whether we use MPI compiler scripts.
# If yes, then change the Fortran compiler to the MPIF77 script.
# If not, then add the MPI include directory for MPI headers.

IF(MPI_USE_MPISCRIPTS )
  # use MPI_MPIF77 as the compiler
  SET(CMAKE_Fortran_COMPILER ${MPI_MPIF77})
ELSE(MPI_USE_MPISCRIPTS)
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(MPI_USE_MPISCRIPTS)

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(KINSOL_LIB sundials_kinsol_static)
  SET(NVECP_LIB sundials_nvecparallel_static)
  SET(FNVECP_LIB sundials_fnvecparallel_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(KINSOL_LIB sundials_kinsol_shared)
  SET(NVECP_LIB sundials_nvecparallel_shared)
  SET(FNVECP_LIB sundials_fnvecparallel_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Only static FCMIX libraries are available

SET(FKINSOL_LIB sundials_fkinsol_static)

# Add the build and install targets for each KINSOL example

FOREACH(example ${fkinsol_examples})
  # add build target for example
  ADD_EXECUTABLE(${example} ${example}.f)
  # link to the sundials libraries
  TARGET_LINK_LIBRARIES(${example} 
    ${FKINSOL_LIB} ${KINSOL_LIB}
    ${FNVECP_LIB} ${NVECP_LIB}
    ${EXTRA_LINK_LIBS})
  # if needed, link to the MPI libraries
  IF(NOT MPI_USE_MPISCRIPTS)
    TARGET_LINK_LIBRARIES(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  ENDIF(NOT MPI_USE_MPISCRIPTS)
  # install example files
  INSTALL(FILES ${example}.f ${example}.out DESTINATION examples/kinsol/fcmix_parallel)
ENDFOREACH(example ${fkinsol_examples})

# Regardless of the platform we're on, we will generate and install 
# CMakeLists.txt file for building the examples. This file  can then 
# be used as a template for the user's own programs.

# generate CMakelists.txt in the binary directory
IF(MPI_USE_MPISCRIPTS)
  FILE(WRITE ${PROJECT_BINARY_DIR}/examples/kinsol/fcmix_parallel/CMakeLists.txt
    "PROJECT(FKINSOL_parallel_examples Fortran)\n"
    "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
    "SET(fkinsol_examples ${fkinsol_examples})\n"
    "SET(MPI_MPIF77 ${MPI_MPIF77} CACHE STRING \"MPIF77 compiler script\")\n"
    "SET(CMAKE_Fortran_COMPILER \${MPI_MPIF77})\n"
    "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
    "FIND_LIBRARY(SUNDIALS_FKINSOL_LIB \n"
    "  sundials_fkinsol \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"FKINSOL wrapper library\")\n"
    "FIND_LIBRARY(SUNDIALS_KINSOL_LIB \n"
    "  sundials_kinsol \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"KINSOL library\")\n"
    "FIND_LIBRARY(SUNDIALS_FNVECP_LIB \n"
    "  sundials_fnvecparallel \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"FNVECPARALLEL wrapper library\")\n"
    "FIND_LIBRARY(SUNDIALS_NVECP_LIB \n"
    "  sundials_nvecparallel \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"NVECPARALLEL library\")\n"
    "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
    "FOREACH(example \${fkinsol_examples})\n"
    "  ADD_EXECUTABLE(\${example} \${example}.f)\n"
    "  TARGET_LINK_LIBRARIES(\${example} \n"
    "    \${SUNDIALS_FKINSOL_LIB} \${SUNDIALS_KINSOL_LIB} \n"
    "    \${SUNDIALS_FNVECP_LIB} \${SUNDIALS_NVECP_LIB} \n"
    "    \${SUNDIALS_EXTRA_LIB})\n"
    "ENDFOREACH(example \${fkinsol_examples})\n"
    )
ELSE(MPI_USE_MPISCRIPTS)
  FILE(WRITE ${PROJECT_BINARY_DIR}/examples/kinsol/fcmix_parallel/CMakeLists.txt
    "PROJECT(FKINSOL_parallel_examples Fortran)\n"
    "MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)\n"
    "SET(fkinsol_examples ${fkinsol_examples})\n"
    "SET(SUNDIALS_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)\n"
    "FIND_LIBRARY(SUNDIALS_FKINSOL_LIB \n"
    "  sundials_fkinsol \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"FKINSOL wrapper library\")\n"
    "FIND_LIBRARY(SUNDIALS_KINSOL_LIB \n"
    "  sundials_kinsol \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"KINSOL library\")\n"
    "FIND_LIBRARY(SUNDIALS_FNVECP_LIB \n"
    "  sundials_fnvecparallel \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"FNVECPARALLEL wrapper library\")\n"
    "FIND_LIBRARY(SUNDIALS_NVECP_LIB \n"
    "  sundials_nvecparallel \${SUNDIALS_LIB_DIR}\n"
    "  DOC \"NVECPARALLEL library\")\n"
    "SET(SUNDIALS_EXTRA_LIB ${EXTRA_LINK_LIBS} CACHE STRING \"Additional libraries\")\n"
    "SET(MPI_INC_DIR ${MPI_INCLUDE_PATH} CACHE STRING \"Location of MPI header files\")\n"
    "INCLUDE_DIRECTORIES(\${MPI_INC_DIR})\n"
    "SET(MPI_LIB ${MPI_LIBRARY} CACHE STRING \"Main MPI library\")\n"
    "SET(MPI_EXTRA_LIB ${MPI_EXTRA_LIBRARIES} CACHE STRING \"Additional MPI libraries\")\n"
    "FOREACH(example \${fkinsol_examples})\n"
    "  ADD_EXECUTABLE(\${example} \${example}.f)\n"
    "  TARGET_LINK_LIBRARIES(\${example} \n"
    "    \${SUNDIALS_FKINSOL_LIB} \${SUNDIALS_KINSOL_LIB} \n"
    "    \${SUNDIALS_FNVECP_LIB} \${SUNDIALS_NVECP_LIB} \n"
    "    \${SUNDIALS_EXTRA_LIB} \n"
    "    \${MPI_LIB} \${MPI_EXTRA_LIB})\n"
    "ENDFOREACH(example \${fkinsol_examples})\n"
    )
ENDIF(MPI_USE_MPISCRIPTS)

# install CMakelists.txt
INSTALL(
  FILES ${PROJECT_BINARY_DIR}/examples/kinsol/fcmix_parallel/CMakeLists.txt
  DESTINATION examples/kinsol/fcmix_parallel 
  )

# On UNIX-type platforms, we have the option of generating and instaling
# a makefile for building the examples. This makefile can then be used
# as a template for the user's own programs.

IF(GENERATE_EX_MAKEFILES)
  # generate Makefile and place it in the binary dir
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/kinsol/fcmix_parallel/Makefile_ex.in
    ${PROJECT_BINARY_DIR}/examples/kinsol/fcmix_parallel/Makefile_ex
    @ONLY
    )
  # install the configured Makefile_ex as Makefile
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/kinsol/fcmix_parallel/Makefile_ex
    DESTINATION examples/kinsol/fcmix_parallel 
    RENAME Makefile
    )
ENDIF(GENERATE_EX_MAKEFILES)
