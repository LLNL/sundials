# ---------------------------------------------------------------
# $Revision$
# $Date$
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
# CMakeLists.txt file for the CPODES library

INSTALL(CODE "MESSAGE(\"\nInstall CPODES\n\")")

# Add variable cpodes_SOURCES with the sources for the CPODES library
SET(cpodes_SOURCES
  cpodes_band.c
  cpodes_bandpre.c
  cpodes_bbdpre.c
  cpodes.c
  cpodes_dense.c
  cpodes_direct.c
  cpodes_ic.c
  cpodes_io.c
  cpodes_nls.c
  cpodes_proj.c
  cpodes_root.c
  cpodes_spbcgs.c
  cpodes_spgmr.c
  cpodes_spils.c
  cpodes_sptfqmr.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the CPODES library
SET(shared_SOURCES
  sundials_nvector.c
  sundials_math.c
  sundials_direct.c
  sundials_band.c
  sundials_dense.c
  sundials_iterative.c
  sundials_spbcgs.c
  sundials_spgmr.c
  sundials_sptfqmr.c
  )

# Add prefix with complete path to the common SUNDIALS sources
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable cpodes_HEADERS with the exported CPODES header files
SET(cpodes_HEADERS
  cpodes_band.h
  cpodes_bandpre.h
  cpodes_bbdpre.h
  cpodes_dense.h
  cpodes_direct.h
  cpodes.h
  cpodes_spbcgs.h
  cpodes_spgmr.h
  cpodes_spils.h
  cpodes_sptfqmr.h
  )

# Add prefix with complete path to the CPODES header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/cpodes/ cpodes_HEADERS)

# If Blas/Lapack support was enabled, set-up additional file lists
IF(LAPACK_FOUND)
  SET(cpodes_BL_SOURCES cpodes_lapack.c)
  SET(cpodes_BL_HEADERS cpodes_lapack.h)
  ADD_PREFIX(${sundials_SOURCE_DIR}/include/cpodes/ cpodes_BL_HEADERS)
ELSE(LAPACK_FOUND)
  SET(cpodes_BL_SOURCES "")
  SET(cpodes_BL_HEADERS "")
ENDIF(LAPACK_FOUND)

# Add source directories to include directories for access to
# implementation only header files.
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
IF(BUILD_STATIC_LIBS)

  # Add the build target for the static CPODES library
  ADD_LIBRARY(sundials_cpodes_static STATIC
    ${cpodes_SOURCES} ${cpodes_BL_SOURCES} ${shared_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_cpodes_static
    PROPERTIES OUTPUT_NAME sundials_cpodes CLEAN_DIRECT_OUTPUT 1)

  # Install the CPODES library
  INSTALL(TARGETS sundials_cpodes_static DESTINATION lib)

ENDIF(BUILD_STATIC_LIBS)

# Build the shared library
IF(BUILD_SHARED_LIBS)

  # Add the build target for the CPODES library
  ADD_LIBRARY(sundials_cpodes_shared SHARED
    ${cpodes_SOURCES} ${cpodes_BL_SOURCES} ${shared_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_cpodes_shared
    PROPERTIES OUTPUT_NAME sundials_cpodes CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  SET_TARGET_PROPERTIES(sundials_cpodes_shared
    PROPERTIES VERSION ${cpodeslib_VERSION} SOVERSION ${cpodeslib_SOVERSION})

  # Install the CPODES library
  INSTALL(TARGETS sundials_cpodes_shared DESTINATION lib)

ENDIF(BUILD_SHARED_LIBS)

# Install the CPODES header files
INSTALL(FILES ${cpodes_HEADERS} ${cpodes_BL_HEADERS} DESTINATION include/cpodes)

# Install the CPODES implementation header file
INSTALL(FILES cpodes_impl.h DESTINATION include/cpodes)

#
MESSAGE(STATUS "Added CPODES module")
