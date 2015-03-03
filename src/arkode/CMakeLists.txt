# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# LLNS/SMU Copyright Start
# Copyright (c) 2015, Southern Methodist University and 
# Lawrence Livermore National Security
#
# This work was performed under the auspices of the U.S. Department 
# of Energy by Southern Methodist University and Lawrence Livermore 
# National Laboratory under Contract DE-AC52-07NA27344.
# Produced at Southern Methodist University and the Lawrence 
# Livermore National Laboratory.
#
# All rights reserved.
# For details, see the LICENSE file.
# LLNS/SMU Copyright End
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the ARKODE library

INSTALL(CODE "MESSAGE(\"\nInstall ARKODE\n\")")

# Add variable arkode_SOURCES with the sources for the ARKODE library
SET(arkode_SOURCES
  arkode.c
  arkode_band.c
  arkode_bandpre.c
  arkode_bbdpre.c
  arkode_butcher.c
  arkode_dense.c
  arkode_direct.c        
  arkode_io.c
  arkode_pcg.c
  arkode_spbcgs.c
  arkode_spfgmr.c
  arkode_spgmr.c
  arkode_spils.c
  arkode_sptfqmr.c
  arkode_sparse.c
  )

IF(KLU_FOUND)
    LIST(APPEND arkode_SOURCES arkode_klu.c)
ENDIF()

IF(SUPERLUMT_FOUND)
    LIST(APPEND arkode_SOURCES arkode_superlumt.c)
ENDIF()

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the ARKODE library
SET(shared_SOURCES
  sundials_nvector.c
  sundials_math.c
  sundials_direct.c
  sundials_band.c
  sundials_dense.c
  sundials_iterative.c
  sundials_pcg.c
  sundials_spbcgs.c
  sundials_spfgmr.c
  sundials_spgmr.c
  sundials_sptfqmr.c
  sundials_sparse.c
  )

# Add prefix with complete path to the common SUNDIALS sources
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable arkode_HEADERS with the exported ARKODE header files
SET(arkode_HEADERS
  arkode.h
  arkode_band.h
  arkode_bandpre.h
  arkode_bbdpre.h
  arkode_dense.h
  arkode_direct.h
  arkode_pcg.h
  arkode_spbcgs.h
  arkode_spfgmr.h
  arkode_spgmr.h
  arkode_spils.h
  arkode_sptfqmr.h
  arkode_sparse.h
  )

IF(KLU_FOUND)
    LIST(APPEND arkode_HEADERS arkode_klu.h)
ENDIF()

IF(SUPERLUMT_FOUND)
    LIST(APPEND arkode_HEADERS arkode_superlumt.h)
ENDIF()

# Add prefix with complete path to the ARKODE header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/arkode/ arkode_HEADERS)

# If Blas/Lapack support was enabled, set-up additional file lists
IF(LAPACK_FOUND)
  SET(arkode_BL_SOURCES arkode_lapack.c)
  SET(arkode_BL_HEADERS arkode_lapack.h)
  ADD_PREFIX(${sundials_SOURCE_DIR}/include/arkode/ arkode_BL_HEADERS)
ELSE(LAPACK_FOUND)
  SET(arkode_BL_SOURCES "")
  SET(arkode_BL_HEADERS "")
ENDIF(LAPACK_FOUND)

# Add source directories to include directories for access to
# implementation only header files.
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
IF(BUILD_STATIC_LIBS)

  # Add the build target for the static ARKODE library
  ADD_LIBRARY(sundials_arkode_static STATIC 
    ${arkode_SOURCES}  ${arkode_BL_SOURCES} ${shared_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_arkode_static
    PROPERTIES OUTPUT_NAME sundials_arkode CLEAN_DIRECT_OUTPUT 1)

  # Install the ARKODE library
  INSTALL(TARGETS sundials_arkode_static DESTINATION lib)

ENDIF(BUILD_STATIC_LIBS)

# Build the shared library
IF(BUILD_SHARED_LIBS)

  # Add the build target for the ARKODE library
  ADD_LIBRARY(sundials_arkode_shared SHARED 
    ${arkode_SOURCES}  ${arkode_BL_SOURCES}  ${shared_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_arkode_shared
    PROPERTIES OUTPUT_NAME sundials_arkode CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  SET_TARGET_PROPERTIES(sundials_arkode_shared
    PROPERTIES VERSION ${arkodelib_VERSION} SOVERSION ${arkodelib_SOVERSION})

  # Install the ARKODE library
  INSTALL(TARGETS sundials_arkode_shared DESTINATION lib)

ENDIF(BUILD_SHARED_LIBS)

# Install the ARKODE header files
INSTALL(FILES ${arkode_HEADERS} ${arkode_BL_HEADERS} DESTINATION include/arkode)

# Install the ARKODE implementation header file
INSTALL(FILES arkode_impl.h DESTINATION include/arkode)

#
MESSAGE(STATUS "Added ARKODE module")
