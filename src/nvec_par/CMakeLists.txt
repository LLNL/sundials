# ---------------------------------------------------------------
# $Revision: 1.1 $
# $Date: 2007-01-29 17:40:00 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the parallel NVECTOR library

IF(MPI_USE_MPISCRIPTS)
  # use MPI_MPICC as the compiler
  SET(CMAKE_C_COMPILER ${MPI_MPICC})
ELSE(MPI_USE_MPISCRIPTS)
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(MPI_USE_MPISCRIPTS)

# Add variable nvecparallel_SOURCES with the sources for the NVECPARALLEL lib
SET(nvecparallel_SOURCES nvector_parallel.c)

# Add variable nvecparallel_HEADERS with the exported NVECPARALLEL header files
SET(nvecparallel_HEADERS nvector_parallel.h)

# Add prefix with complete path to the NVECPARALLEL header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nvecparallel_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
IF(BUILD_STATIC_LIBS)

  # Add the build target for the NVECPARALLEL library
  ADD_LIBRARY(sundials_nvecparallel_static STATIC ${nvecparallel_SOURCES})
  
  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_nvecparallel_static
    PROPERTIES OUTPUT_NAME sundials_nvecparallel CLEAN_DIRECT_OUTPUT 1)

  # Install the NVECPARALLEL library
  INSTALL(TARGETS sundials_nvecparallel_static DESTINATION lib)

ENDIF(BUILD_STATIC_LIBS)

# Build the shared library
IF(BUILD_SHARED_LIBS)

  # Add the build target for the NVECPARALLEL library
  ADD_LIBRARY(sundials_nvecparallel_shared SHARED ${nvecparallel_SOURCES})
  
  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_nvecparallel_shared
    PROPERTIES OUTPUT_NAME sundials_nvecparallel CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  SET_TARGET_PROPERTIES(sundials_nvecparallel_shared 
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  
  # Install the NVECPARALLEL library
  INSTALL(TARGETS sundials_nvecparallel_shared DESTINATION lib)

ENDIF(BUILD_SHARED_LIBS)

# Install the NVECPARALLEL header files
INSTALL(FILES ${nvecparallel_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled and MPI-F77 works, build and install the FNVECPARALLEL library
IF(MPIF_OK)
  SET(fnvecparallel_SOURCES fnvector_parallel.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fnvecparallel_static STATIC ${fnvecparallel_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecparallel_static
      PROPERTIES OUTPUT_NAME sundials_fnvecparallel CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fnvecparallel_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fnvecparallel_shared SHARED ${fnvecparallel_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecparallel_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecparallel CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fnvecparallel_shared
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    INSTALL(TARGETS sundials_fnvecparallel_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(MPIF_OK)
