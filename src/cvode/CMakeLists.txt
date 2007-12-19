# ---------------------------------------------------------------
# $Revision: 1.3 $
# $Date: 2007-12-19 20:34:00 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the CVODE library

INSTALL(CODE "MESSAGE(\"\nInstall CVODE\n\")")

# Add variable cvode_SOURCES with the sources for the CVODE library
SET(cvode_SOURCES
  cvode.c
  cvode_io.c
  cvode_direct.c        
  cvode_band.c
  cvode_dense.c
  cvode_diag.c
  cvode_spils.c
  cvode_spbcgs.c
  cvode_spgmr.c
  cvode_sptfqmr.c
  cvode_bandpre.c
  cvode_bbdpre.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the CVODE library
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

# Add variable cvode_HEADERS with the exported CVODE header files
SET(cvode_HEADERS
  cvode_band.h
  cvode_bandpre.h
  cvode_bbdpre.h
  cvode_dense.h
  cvode_diag.h
  cvode_direct.h
  cvode.h
  cvode_spbcgs.h
  cvode_spgmr.h
  cvode_spils.h
  cvode_sptfqmr.h
  )

# Add prefix with complete path to the CVODE header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/cvode/ cvode_HEADERS)

# If Blas/Lapack support was enabled, set-up additional file lists
IF(LAPACK_OK)
  SET(cvode_BL_SOURCES cvode_lapack.c)
  SET(shared_BL_SOURCES sundials_lapack.c)
  ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_BL_SOURCES)
  SET(cvode_BL_HEADERS cvode_lapack.h)
  ADD_PREFIX(${sundials_SOURCE_DIR}/include/cvode/ cvode_BL_HEADERS)
ELSE(LAPACK_OK)
  SET(cvode_BL_SOURCES "")
  SET(shared_BL_SOURCES "")
  SET(cvode_BL_HEADERS "")
ENDIF(LAPACK_OK)

# Add source directories to include directories for access to
# implementation only header files.
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
IF(BUILD_STATIC_LIBS)

  # Add the build target for the static CVODE library
  ADD_LIBRARY(sundials_cvode_static STATIC 
    ${cvode_SOURCES}  ${cvode_BL_SOURCES} 
    ${shared_SOURCES} ${shared_BL_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_cvode_static
    PROPERTIES OUTPUT_NAME sundials_cvode CLEAN_DIRECT_OUTPUT 1)

  # Install the CVODE library
  INSTALL(TARGETS sundials_cvode_static DESTINATION lib)

ENDIF(BUILD_STATIC_LIBS)

# Build the shared library
IF(BUILD_SHARED_LIBS)

  # Add the build target for the CVODE library
  ADD_LIBRARY(sundials_cvode_shared SHARED 
    ${cvode_SOURCES}  ${cvode_BL_SOURCES} 
    ${shared_SOURCES} ${shared_BL_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_cvode_shared
    PROPERTIES OUTPUT_NAME sundials_cvode CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  SET_TARGET_PROPERTIES(sundials_cvode_shared
    PROPERTIES VERSION ${cvodelib_VERSION} SOVERSION ${cvodelib_SOVERSION})

  # Install the CVODE library
  INSTALL(TARGETS sundials_cvode_shared DESTINATION lib)

ENDIF(BUILD_SHARED_LIBS)

# Install the CVODE header files
INSTALL(FILES ${cvode_HEADERS} ${cvode_BL_HEADERS} DESTINATION include/cvode)

# Install the CVODE implementation header file
INSTALL(FILES cvode_impl.h DESTINATION include/cvode)

