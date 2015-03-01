# CMakeLists.txt file for the FIDA library

# Add variable fida_SOURCES with the sources for the FIDA library
SET(fida_SOURCES
  fidaband.c
  fidabbd.c
  fida.c
  fidadense.c
  fidaewt.c
  fidajtimes.c
  fidapreco.c
  fidaroot.c
  )

IF(LAPACK_FOUND)
  SET(fida_BL_SOURCES fidalapack.c fidalapdense.c fidalapband.c)
ELSE(LAPACK_FOUND)
  SET(fida_BL_SOURCES "")
ENDIF(LAPACK_FOUND)

IF(KLU_FOUND OR SUPERLUMT_FOUND)
  SET(fida_SPARSE_SOURCES fidasparse.c)
ELSE(SPARSE_FOUND)
  SET(fida_SPARSE_SOURCES "")
ENDIF(KLU_FOUND OR SUPERLUMT_FOUND)
 
IF(KLU_FOUND)
  SET(fida_KLU_SOURCES fidaklu.c)
ELSE(KLU_FOUND)
  SET(fida_KLU_SOURCES "")
ENDIF(KLU_FOUND)
 
IF(SUPERLUMT_FOUND)
  SET(fida_SUPERLUMT_SOURCES fidasuperlumt.c)
ELSE(SUPERLUMT_FOUND)
  SET(fida_SUPERLUMT_SOURCES "")
ENDIF(SUPERLUMT_FOUND)

# Add source directories to include directories for access to
# implementation only header files (both for fida and ida)
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(..)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Only build STATIC libraries (we cannot build shared libraries 
# for the FCMIX interfaces due to unresolved symbol errors 
# coming from inexistent user-provided functions)

# Add the build target for the FIDA library
ADD_LIBRARY(sundials_fida_static STATIC ${fida_SOURCES} ${fida_BL_SOURCES} ${fida_SPARSE_SOURCES} ${fida_KLU_SOURCES} ${fida_SUPERLUMT_SOURCES})

# Set the library name and make sure it is not deleted
SET_TARGET_PROPERTIES(sundials_fida_static
  PROPERTIES OUTPUT_NAME sundials_fida CLEAN_DIRECT_OUTPUT 1)

# Install the FIDA library
INSTALL(TARGETS sundials_fida_static DESTINATION lib)
#
MESSAGE(STATUS "Added IDA FCMIX module")
