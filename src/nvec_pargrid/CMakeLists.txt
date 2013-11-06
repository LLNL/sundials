# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2013, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the parallel grid NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_PARALLEL_GRID\n\")")

IF(MPI_MPICC)
  # use MPI_MPICC as the compiler
  SET(CMAKE_C_COMPILER ${MPI_MPICC})
ELSE(MPI_MPICC)
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(MPI_MPICC)

# Add variable nvecparallelgrid_SOURCES with the sources for the NVECPARALLELGRID lib
SET(nvecparallelgrid_SOURCES nvector_parallel_grid.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECPARALLELGRID library
SET(shared_SOURCES sundials_math.c)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nvecparallelgrid_HEADERS with the exported NVECPARALLELGRID header files
SET(nvecparallelgrid_HEADERS nvector_parallel_grid.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nvecparallelgrid_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECPARALLELGRID library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECPARALLELGRID library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_nvecparallelgrid_static STATIC ${nvecparallelgrid_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecparallelgrid_static
    PROPERTIES OUTPUT_NAME sundials_nvecparallelgrid CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_nvecparallelgrid_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECPARALLELGRID library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECPARALLELGRID library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_nvecparallelgrid_shared SHARED ${nvecparallelgrid_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecparallelgrid_shared
    PROPERTIES OUTPUT_NAME sundials_nvecparallelgrid CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_nvecparallelgrid_shared 
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nvecparallelgrid_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECPARALLELGRID header files
INSTALL(FILES ${nvecparallelgrid_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled and MPI-F77 works, build and install the FNVECPARALLELGRID library
IF(FCMIX_ENABLE AND MPIF_FOUND)
  SET(fnvecparallelgrid_SOURCES fnvector_parallel_grid.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fnvecparallelgrid_static STATIC ${fnvecparallelgrid_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecparallelgrid_static
      PROPERTIES OUTPUT_NAME sundials_fnvecparallelgrid CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fnvecparallelgrid_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fnvecparallelgrid_shared SHARED ${fnvecparallelgrid_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecparallelgrid_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecparallelgrid CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fnvecparallelgrid_shared
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    INSTALL(TARGETS sundials_fnvecparallelgrid_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(FCMIX_ENABLE AND MPIF_FOUND)
#
MESSAGE(STATUS "Added NVECTOR_PARALLEL_GRID module")
