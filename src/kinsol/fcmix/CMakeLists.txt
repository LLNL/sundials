# CMakeLists.txt file for the FKINSOL library

# Add variable fcvode_SOURCES with the sources for the FCVODE library
SET(fkinsol_SOURCES
  fkinband.c
  fkinbbd.c
  fkindense.c
  fkinjtimes.c
  fkinpreco.c
  fkinsol.c
  )

IF(LAPACK_FOUND)
  SET(fkinsol_BL_SOURCES fkinlapack.c fkinlapdense.c fkinlapband.c)
ELSE(LAPACK_FOUND)
  SET(fkinsol_BL_SOURCES "")
ENDIF(LAPACK_FOUND)

IF(KLU_FOUND OR SUPERLUMT_FOUND)
  SET(fkinsol_SPARSE_SOURCES fkinsparse.c)
ELSE(SPARSE_FOUND)
  SET(fkinsol_SPARSE_SOURCES "")
ENDIF(KLU_FOUND OR SUPERLUMT_FOUND)
 
IF(KLU_FOUND)
  SET(fkinsol_KLU_SOURCES fkinklu.c)
ELSE(KLU_FOUND)
  SET(fkinsol_KLU_SOURCES "")
ENDIF(KLU_FOUND)
 
IF(SUPERLUMT_FOUND)
  SET(fkinsol_SUPERLUMT_SOURCES fkinsuperlumt.c)
ELSE(SUPERLUMT_FOUND)
  SET(fkinsol_SUPERLUMT_SOURCES "")
ENDIF(SUPERLUMT_FOUND)

# Add source directories to include directories for access to
# implementation only header files (both for fkinsol and kinsol)
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(..)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Only build STATIC libraries (we cannot build shared libraries 
# for the FCMIX interfaces due to unresolved symbol errors 
# coming from inexistent user-provided functions)

# Add the build target for the FKINSOL library
ADD_LIBRARY(sundials_fkinsol_static STATIC ${fkinsol_SOURCES} ${fkinsol_BL_SOURCES} ${fkinsol_SPARSE_SOURCES} ${fkinsol_KLU_SOURCES} ${fkinsol_SUPERLUMT_SOURCES})

# Set the library name and make sure it is not deleted
SET_TARGET_PROPERTIES(sundials_fkinsol_static
  PROPERTIES OUTPUT_NAME sundials_fkinsol CLEAN_DIRECT_OUTPUT 1)

# Install the FKINSOL library
INSTALL(TARGETS sundials_fkinsol_static DESTINATION lib)

#
MESSAGE(STATUS "Added KINSOL FCMIX module")
