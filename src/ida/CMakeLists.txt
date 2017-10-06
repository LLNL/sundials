# ---------------------------------------------------------------
# $Revision: 4961 $
# $Date: 2016-09-26 07:58:48 -0700 (Mon, 26 Sep 2016) $
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
# CMakeLists.txt file for the IDA library

INSTALL(CODE "MESSAGE(\"\nInstall IDA\n\")")

# Add variable ida_SOURCES with the sources for the IDA library
SET(ida_SOURCES
  ida.c
  ida_io.c
  ida_ic.c
  ida_direct.c
  ida_band.c
  ida_dense.c
  ida_sparse.c
  ida_spils.c
  ida_spbcgs.c
  ida_spgmr.c
  ida_sptfqmr.c
  ida_bbdpre.c
  )

IF(KLU_FOUND)
    LIST(APPEND ida_SOURCES ida_klu.c)
ENDIF()

IF(SUPERLUMT_FOUND)
    LIST(APPEND ida_SOURCES ida_superlumt.c)
ENDIF()

# IF(PETSC_FOUND)
#    # use MPI_MPICC as the compiler
#    SET(CMAKE_C_COMPILER ${MPI_MPICC})
#    LIST(APPEND ida_SOURCES ida_petsc.c)
# ENDIF()

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the IDA library
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

# Add variable ida_HEADERS with the exported IDA header files
SET(ida_HEADERS
  ida_band.h
  ida_bbdpre.h
  ida_dense.h
  ida_direct.h
  ida.h
  ida_spbcgs.h
  ida_spgmr.h
  ida_spils.h
  ida_sptfqmr.h
  )

IF(KLU_FOUND)
    LIST(APPEND ida_HEADERS ida_klu.h ida_sparse.h)
ENDIF()

IF(SUPERLUMT_FOUND)
    LIST(APPEND ida_HEADERS ida_superlumt.h ida_sparse.h)
ENDIF()

# IF(PETSC_FOUND)
#     LIST(APPEND ida_HEADERS ida_petsc.h)
# ENDIF()

# Add prefix with complete path to the IDA header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/ida/ ida_HEADERS)

# If Blas/Lapack support was enabled, set-up additional file lists
IF(LAPACK_FOUND)
  SET(ida_BL_SOURCES ida_lapack.c)
  SET(ida_BL_HEADERS ida_lapack.h)
  ADD_PREFIX(${sundials_SOURCE_DIR}/include/ida/ ida_BL_HEADERS)
ELSE(LAPACK_FOUND)
  SET(ida_BL_SOURCES "")
  SET(ida_BL_HEADERS "")
ENDIF(LAPACK_FOUND)

# Add source directories to include directories for access to
# implementation only header files.
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
IF(BUILD_STATIC_LIBS)

  # Add the build target for the static IDA library
  ADD_LIBRARY(sundials_ida_static STATIC 
    ${ida_SOURCES}  ${ida_BL_SOURCES}  ${shared_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_ida_static
    PROPERTIES OUTPUT_NAME sundials_ida CLEAN_DIRECT_OUTPUT 1)

  # Install the IDA library
  INSTALL(TARGETS sundials_ida_static DESTINATION lib)

ENDIF(BUILD_STATIC_LIBS)

# Build the shared library
IF(BUILD_SHARED_LIBS)

  # Add the build target for the IDA library
  ADD_LIBRARY(sundials_ida_shared SHARED
    ${ida_SOURCES}  ${ida_BL_SOURCES} ${shared_SOURCES})

  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_ida_shared
    PROPERTIES OUTPUT_NAME sundials_ida CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  SET_TARGET_PROPERTIES(sundials_ida_shared
    PROPERTIES VERSION ${idalib_VERSION} SOVERSION ${idalib_SOVERSION})

  # Install the IDA library
  INSTALL(TARGETS sundials_ida_shared DESTINATION lib)

ENDIF(BUILD_SHARED_LIBS)

# Install the IDA header files
INSTALL(FILES ${ida_HEADERS} ${ida_BL_HEADERS} DESTINATION include/ida)

# Install the IDA implementation header file
INSTALL(FILES ida_impl.h DESTINATION include/ida)

#
MESSAGE(STATUS "Added IDA module")
