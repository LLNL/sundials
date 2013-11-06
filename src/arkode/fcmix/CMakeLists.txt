# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# begincopyright(llns)
# Copyright (c) 2013, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# endcopyright(llns)
# ---------------------------------------------------------------
# CMakeLists.txt file for the FARKODE library

# Add variable farkode_SOURCES with the sources for the FARKODE library
SET(farkode_SOURCES
  farkadapt.c
  farkband.c
  farkbbd.c
  farkbp.c
  farkdense.c
  farkewt.c
  farkexpstab.c
  farkjtimes.c
  farkode.c
  farkpreco.c
  farkroot.c
  )

IF(LAPACK_FOUND)
  SET(farkode_BL_SOURCES farklapack.c)
ELSE(LAPACK_FOUND)
  SET(farkode_BL_SOURCES "")
ENDIF(LAPACK_FOUND)

# Add source directories to include directories for access to
# implementation only header files (both for farkode and arkode)
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(..)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Only build STATIC libraries (we cannot build shared libraries 
# for the FCMIX interfaces due to unresolved symbol errors 
# coming from inexistent user-provided functions)

# Add the build target for the FARKODE library
ADD_LIBRARY(sundials_farkode_static STATIC ${farkode_SOURCES} ${farkode_BL_SOURCES})

# Set the library name and make sure it is not deleted
SET_TARGET_PROPERTIES(sundials_farkode_static
  PROPERTIES OUTPUT_NAME sundials_farkode CLEAN_DIRECT_OUTPUT 1)

# Install the FARKODE library
INSTALL(TARGETS sundials_farkode_static DESTINATION lib)

#
MESSAGE(STATUS "Added ARKODE FCMIX module")
