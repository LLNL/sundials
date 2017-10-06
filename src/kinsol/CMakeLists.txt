# ---------------------------------------------------------------
# $Revision: 4917 $
# $Date: 2016-09-18 23:24:29 -0700 (Sun, 18 Sep 2016) $
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
# CMakeLists.txt file for the KINSOL library

INSTALL(CODE "MESSAGE(\"\nInstall KINSOL\n\")")

# Add variable kinsol_SOURCES with the sources for the KINSOL library
SET(kinsol_SOURCES
  kinsol.c
  kinsol_io.c
  kinsol_direct.c
  kinsol_band.c
  kinsol_dense.c
  kinsol_sparse.c
  kinsol_spils.c
  kinsol_spbcgs.c
  kinsol_spfgmr.c
  kinsol_spgmr.c
  kinsol_sptfqmr.c
  kinsol_bbdpre.c
  )

IF(KLU_FOUND)
    LIST(APPEND kinsol_SOURCES kinsol_klu.c)
ENDIF()

IF(SUPERLUMT_FOUND)
    LIST(APPEND kinsol_SOURCES kinsol_superlumt.c)
ENDIF()

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the KINSOL library
SET(shared_SOURCES
  sundials_nvector.c
  sundials_math.c
  sundials_direct.c
  sundials_band.c
  sundials_dense.c
  sundials_iterative.c
  sundials_sparse.c
  sundials_spbcgs.c
  sundials_spfgmr.c
  sundials_spgmr.c
  sundials_sptfqmr.c
  )

# Add prefix with complete path to the common SUNDIALS sources
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable kinsol_HEADERS with the exported KINSOL header files
SET(kinsol_HEADERS
  kinsol_band.h
  kinsol_bbdpre.h
  kinsol_dense.h
  kinsol_direct.h
  kinsol.h
  kinsol_spbcgs.h
  kinsol_spfgmr.h
  kinsol_spgmr.h
  kinsol_spils.h
  kinsol_sptfqmr.h
  )

IF(KLU_FOUND)
    LIST(APPEND kinsol_HEADERS kinsol_klu.h kinsol_sparse.h)
ENDIF()

IF(SUPERLUMT_FOUND)
    LIST(APPEND kinsol_HEADERS kinsol_superlumt.h kinsol_sparse.h)
ENDIF()

# Add prefix with complete path to the KINSOL header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/kinsol/ kinsol_HEADERS)

# If Blas/Lapack support was enabled, set-up additional file lists
IF(LAPACK_FOUND)
  SET(kinsol_BL_SOURCES kinsol_lapack.c)
  SET(kinsol_BL_HEADERS kinsol_lapack.h)
  ADD_PREFIX(${sundials_SOURCE_DIR}/include/kinsol/ kinsol_BL_HEADERS)
ELSE(LAPACK_FOUND)
  SET(kinsol_BL_SOURCES "")
  SET(kinsol_BL_HEADERS "")
ENDIF(LAPACK_FOUND)

IF(SUPERLUMT_FOUND)
  SET(kinsol_SLU_SOURCES kinsol_superlumt.c)
  SET(kinsol_SLU_HEADERS kinsol_superlumt.h)
  ADD_PREFIX(${sundials_SOURCE_DIR}/include/kinsol/ kinsol_SLU_HEADERS)
ELSE(SUPERLUMT_FOUND)
  SET(kinsol_SLU_SOURCES "")
  SET(kinsol_SLU_HEADERS "")
ENDIF(SUPERLUMT_FOUND)

IF(KLU_FOUND)
  SET(kinsol_KLU_SOURCES kinsol_klu.c)
  SET(kinsol_KLU_HEADERS kinsol_klu.h)
  ADD_PREFIX(${sundials_SOURCE_DIR}/include/kinsol/ kinsol_KLU_HEADERS)
ELSE(KLU_FOUND)
  SET(kinsol_KLU_SOURCES "")
  SET(kinsol_KLU_HEADERS "")
ENDIF(KLU_FOUND)

# Add source directories to include directories for access to
# implementation only header files.
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
IF(BUILD_STATIC_LIBS)

  # Add the build target for the static KINSOL library
  ADD_LIBRARY(sundials_kinsol_static STATIC
    ${kinsol_SOURCES} ${kinsol_BL_SOURCES} ${kinsol_SLU_SOURCES} ${kinsol_KLU_SOURCES} ${shared_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_kinsol_static
    PROPERTIES OUTPUT_NAME sundials_kinsol CLEAN_DIRECT_OUTPUT 1)

  # Install the KINSOL library
  INSTALL(TARGETS sundials_kinsol_static DESTINATION lib)

ENDIF(BUILD_STATIC_LIBS)

# Build the shared library
IF(BUILD_SHARED_LIBS)

  # Add the build target for the KINSOL library
  ADD_LIBRARY(sundials_kinsol_shared SHARED
    ${kinsol_SOURCES} ${kinsol_BL_SOURCES} ${kinsol_SLU_SOURCES} ${kinsol_KLU_SOURCES} ${shared_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_kinsol_shared
    PROPERTIES OUTPUT_NAME sundials_kinsol CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  SET_TARGET_PROPERTIES(sundials_kinsol_shared
    PROPERTIES VERSION ${kinsollib_VERSION} SOVERSION ${kinsollib_SOVERSION})

  # Install the KINSOL library
  INSTALL(TARGETS sundials_kinsol_shared DESTINATION lib)

ENDIF(BUILD_SHARED_LIBS)

# Install the KINSOL header files
INSTALL(FILES ${kinsol_HEADERS} ${kinsol_BL_HEADERS} ${kinsol_SLU_HEADERS} ${kinsol_KLU_HEADERS} DESTINATION include/kinsol)

# Install the KINSOL implementation header file
INSTALL(FILES kinsol_impl.h DESTINATION include/kinsol)

#
MESSAGE(STATUS "Added KINSOL module")
