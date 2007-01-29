# ---------------------------------------------------------------
# $Revision: 1.1 $
# $Date: 2007-01-29 17:39:59 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the CVODES library


# Add variable cvodes_SOURCES with the sources for the CVODES library
SET(cvodes_SOURCES
  cvodes.c
  cvodea.c
  cvodes_io.c
  cvodea_io.c
  cvodes_direct.c
  cvodes_band.c
  cvodes_dense.c
  cvodes_diag.c
  cvodes_spils.c
  cvodes_spbcgs.c
  cvodes_spgmr.c
  cvodes_sptfqmr.c
  cvodes_bandpre.c
  cvodes_bbdpre.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the CVODES library
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

# Add variable cvodes_HEADERS with the exported CVODES header files
SET(cvodes_HEADERS
  cvodes_band.h
  cvodes_bandpre.h
  cvodes_bbdpre.h
  cvodes_dense.h
  cvodes_diag.h
  cvodes_direct.h
  cvodes.h
  cvodes_spbcgs.h
  cvodes_spgmr.h
  cvodes_spils.h
  cvodes_sptfqmr.h
  )

# Add prefix with complete path to the CVODES header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/cvodes/ cvodes_HEADERS)

# Add source directories to include directories for access to
# implementation only header files.
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
IF(BUILD_STATIC_LIBS)

  # Add the build target for the static CVODES library
  ADD_LIBRARY(sundials_cvodes_static STATIC ${cvodes_SOURCES} ${shared_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_cvodes_static
    PROPERTIES OUTPUT_NAME sundials_cvodes CLEAN_DIRECT_OUTPUT 1)

  # Install the CVODES library
  INSTALL(TARGETS sundials_cvodes_static DESTINATION lib)

ENDIF(BUILD_STATIC_LIBS)

# Build the shared library
IF(BUILD_SHARED_LIBS)

  # Add the build target for the CVODES library
  ADD_LIBRARY(sundials_cvodes_shared SHARED ${cvodes_SOURCES} ${shared_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_cvodes_shared
    PROPERTIES OUTPUT_NAME sundials_cvodes CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  SET_TARGET_PROPERTIES(sundials_cvodes_shared
    PROPERTIES VERSION ${cvodeslib_VERSION} SOVERSION ${cvodeslib_SOVERSION})

  # Install the CVODES library
  INSTALL(TARGETS sundials_cvodes_shared DESTINATION lib)

ENDIF(BUILD_SHARED_LIBS)

# Install the CVODES header files
INSTALL(FILES ${cvodes_HEADERS} DESTINATION include/cvodes)

# Install the CVODES implementation header file
INSTALL(FILES cvodes_impl.h DESTINATION include/cvodes)

# If Blas/Lapack support was enabled, create and install the wrapper library
# and install the interface headers
IF(LAPACK_OK)
  SET(cvodes_BL_SOURCES cvodes_lapack.c)
  SET(shared_BL_SOURCES sundials_lapack.c)
  ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_BL_SOURCES)
  SET(cvodes_BL_HEADERS cvodes_lapack.h)
  ADD_PREFIX(${sundials_SOURCE_DIR}/include/cvodes/ cvodes_BL_HEADERS)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_cvodes_bl_static STATIC ${cvodes_BL_SOURCES} ${shared_BL_SOURCES})
    SET_TARGET_PROPERTIES(sundials_cvodes_bl_static
      PROPERTIES OUTPUT_NAME sundials_cvodes_bl CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_cvodes_bl_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_cvodes_bl_shared SHARED ${cvodes_BL_SOURCES} ${shared_BL_SOURCES})
    SET_TARGET_PROPERTIES(sundials_cvodes_bl_shared
      PROPERTIES OUTPUT_NAME sundials_cvodes_bl CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_cvodes_bl_shared
      PROPERTIES VERSION ${cvodeslib_VERSION} SOVERSION ${cvodeslib_SOVERSION})
    INSTALL(TARGETS sundials_cvodes_bl_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
  INSTALL(FILES ${cvodes_BL_HEADERS} DESTINATION include/cvodes)
ENDIF(LAPACK_OK)
