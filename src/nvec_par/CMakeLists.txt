# ---------------------------------------------------------------
# $Revision: 1.3 $
# $Date: 2009-02-17 02:58:48 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the parallel NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_PARALLEL\n\")")

IF(MPI_MPICC)
  # use MPI_MPICC as the compiler
  SET(CMAKE_C_COMPILER ${MPI_MPICC})
ELSE(MPI_MPICC)
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(MPI_MPICC)

# Add variable nvecparallel_SOURCES with the sources for the NVECPARALLEL lib
SET(nvecparallel_SOURCES nvector_parallel.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECPARALLEL library
SET(shared_SOURCES sundials_math.c)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nvecparallel_HEADERS with the exported NVECPARALLEL header files
SET(nvecparallel_HEADERS nvector_parallel.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nvecparallel_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECPARALLEL library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECSERIAL library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_nvecparallel_static STATIC ${nvecparallel_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecparallel_static
    PROPERTIES OUTPUT_NAME sundials_nvecparallel CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_nvecparallel_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECPARALLEL library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECSERIAL library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_nvecparallel_shared SHARED ${nvecparallel_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecparallel_shared
    PROPERTIES OUTPUT_NAME sundials_nvecparallel CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_nvecparallel_shared 
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nvecparallel_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECPARALLEL header files
INSTALL(FILES ${nvecparallel_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled and MPI-F77 works, build and install the FNVECPARALLEL library
IF(FCMIX_ENABLE AND MPIF_FOUND)
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
ENDIF(FCMIX_ENABLE AND MPIF_FOUND)
#
MESSAGE(STATUS "Added NVECTOR_PARALLEL module")
