# ---------------------------------------------------------------
# $Revision: 1.1 $
# $Date: 2007-10-26 16:36:33 $
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
  pvkt
  pvfkt
  pvakt
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

# Add the build target for each CVODES example

FOREACH(example ${cvodes_examples})
  # add build target for example
  ADD_EXECUTABLE(${example} ${example}.c)
  # link to the sundials libraries
  TARGET_LINK_LIBRARIES(${example} ${CVODES_LIB} ${NVECP_LIB} ${EXTRA_LINK_LIBS})
  # if needed, link to the MPI libraries
  IF(NOT MPI_USE_MPISCRIPTS)
    TARGET_LINK_LIBRARIES(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  ENDIF(NOT MPI_USE_MPISCRIPTS)
ENDFOREACH(example ${cvodes_examples})

