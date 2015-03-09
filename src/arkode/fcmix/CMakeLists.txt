# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the FARKODE library

# Add variable farkode_SOURCES with the sources for the FARKODE library
SET(farkode_SOURCES
  farkadapt.c
  farkband.c
  farkbandmass.c
  farkbbd.c
  farkbp.c
  farkdense.c
  farkdensemass.c
  farkewt.c
  farkexpstab.c
  farkjtimes.c
  farkmtimes.c
  farkode.c
  farkpreco.c
  farkmasspreco.c
  farkroot.c
  )

IF(LAPACK_FOUND)
  SET(farkode_BL_SOURCES farklapack.c farklapdense.c farklapdensemass.c farklapband.c farklapbandmass.c)
ELSE(LAPACK_FOUND)
  SET(farkode_BL_SOURCES "")
ENDIF(LAPACK_FOUND)

IF(KLU_FOUND OR SUPERLUMT_FOUND)
  SET(farkode_SPARSE_SOURCES farksparse.c farksparsemass.c)
ELSE(SPARSE_FOUND)
  SET(farkode_SPARSE_SOURCES "")
ENDIF(KLU_FOUND OR SUPERLUMT_FOUND)

IF(KLU_FOUND)
  SET(farkode_KLU_SOURCES farkklu.c)
ELSE(KLU_FOUND)
  SET(farkode_KLU_SOURCES "")
ENDIF(KLU_FOUND)

IF(SUPERLUMT_FOUND)
  SET(farkode_SUPERLUMT_SOURCES farksuperlumt.c)
ELSE(SUPERLUMT_FOUND)
  SET(farkode_SUPERLUMT_SOURCES "")
ENDIF(SUPERLUMT_FOUND)

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
ADD_LIBRARY(sundials_farkode_static STATIC ${farkode_SOURCES} ${farkode_BL_SOURCES} ${farkode_SPARSE_SOURCES} ${farkode_KLU_SOURCES} ${farkode_SUPERLUMT_SOURCES})

# Set the library name and make sure it is not deleted
SET_TARGET_PROPERTIES(sundials_farkode_static
  PROPERTIES OUTPUT_NAME sundials_farkode CLEAN_DIRECT_OUTPUT 1)

# Install the FARKODE library
INSTALL(TARGETS sundials_farkode_static DESTINATION lib)

#
MESSAGE(STATUS "Added ARKODE FCMIX module")
