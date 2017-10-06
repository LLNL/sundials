# ---------------------------------------------------------------
# $Revision: 4759 $
# $Date: 2016-05-18 16:25:17 -0700 (Wed, 18 May 2016) $
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
# CMakeLists.txt file for the IDAS library

INSTALL(CODE "MESSAGE(\"\nInstall IDAS\n\")")

# Add variable idas_SOURCES with the sources for the IDAS library
SET(idas_SOURCES
  idas.c
  idaa.c
  idas_io.c
  idas_ic.c
  idaa_io.c
  idas_direct.c
  idas_band.c
  idas_dense.c
  idas_sparse.c
  idas_spils.c
  idas_spbcgs.c
  idas_spgmr.c
  idas_sptfqmr.c
  idas_bbdpre.c
  )

IF(KLU_FOUND)
    LIST(APPEND idas_SOURCES idas_klu.c)
ENDIF()

IF(SUPERLUMT_FOUND)
    LIST(APPEND idas_SOURCES idas_superlumt.c)
ENDIF()

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the IDAS library
SET(shared_SOURCES
  sundials_nvector.c
  sundials_math.c
  sundials_direct.c
  sundials_band.c
  sundials_dense.c
  sundials_iterative.c
  sundials_sparse.c
  sundials_spbcgs.c
  sundials_spgmr.c
  sundials_sptfqmr.c
  )

# Add prefix with complete path to the common SUNDIALS sources
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable idas_HEADERS with the exported IDAS header files
SET(idas_HEADERS
  idas_band.h
  idas_bbdpre.h
  idas_dense.h
  idas_direct.h
  idas.h
  idas_sparse.h
  idas_spbcgs.h
  idas_spgmr.h
  idas_spils.h
  idas_sptfqmr.h
  )

IF(KLU_FOUND)
    LIST(APPEND idas_HEADERS idas_klu.h)
ENDIF()

IF(SUPERLUMT_FOUND)
    LIST(APPEND idas_HEADERS idas_superlumt.h)
ENDIF()

# Add prefix with complete path to the IDAS header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/idas/ idas_HEADERS)

# If Blas/Lapack support was enabled, set-up additional file lists
IF(LAPACK_FOUND)
  SET(idas_BL_SOURCES idas_lapack.c)
  SET(idas_BL_HEADERS idas_lapack.h)
  ADD_PREFIX(${sundials_SOURCE_DIR}/include/idas/ idas_BL_HEADERS)
ELSE(LAPACK_FOUND)
  SET(idas_BL_SOURCES "")
  SET(idas_BL_HEADERS "")
ENDIF(LAPACK_FOUND)

IF(SUPERLUMT_FOUND)
  SET(idas_SLU_SOURCES idas_superlumt.c)
  SET(idas_SLU_HEADERS idas_superlumt.h)
  ADD_PREFIX(${sundials_SOURCE_DIR}/include/idas/ idas_SLU_HEADERS)
ELSE(SUPERLUMT_FOUND)
  SET(idas_SLU_SOURCES "")
  SET(idas_SLU_HEADERS "")
ENDIF(SUPERLUMT_FOUND)

# Add source directories to include directories for access to
# implementation only header files.
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
IF(BUILD_STATIC_LIBS)

  # Add the build target for the static IDAS library
  ADD_LIBRARY(sundials_idas_static STATIC 
    ${idas_SOURCES} ${idas_BL_SOURCES} ${idas_SLU_SOURCES}  ${shared_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_idas_static
    PROPERTIES OUTPUT_NAME sundials_idas CLEAN_DIRECT_OUTPUT 1)

  # Install the IDA library
  INSTALL(TARGETS sundials_idas_static DESTINATION lib)

ENDIF(BUILD_STATIC_LIBS)

# Build the shared library
IF(BUILD_SHARED_LIBS)

  # Add the build target for the IDAS library
  ADD_LIBRARY(sundials_idas_shared SHARED 
    ${idas_SOURCES}  ${idas_BL_SOURCES} ${idas_SLU_SOURCES} ${shared_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_idas_shared
    PROPERTIES OUTPUT_NAME sundials_idas CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  SET_TARGET_PROPERTIES(sundials_idas_shared
    PROPERTIES VERSION ${idaslib_VERSION} SOVERSION ${idaslib_SOVERSION})

  # Install the IDAS library
  INSTALL(TARGETS sundials_idas_shared DESTINATION lib)

ENDIF(BUILD_SHARED_LIBS)

# Install the IDAS header files
INSTALL(FILES ${idas_HEADERS} ${idas_BL_HEADERS} ${idas_SLU_HEADERS} DESTINATION include/idas)

# Install the IDAS implementation header file
INSTALL(FILES idas_impl.h DESTINATION include/idas)

#
MESSAGE(STATUS "Added IDAS module")
