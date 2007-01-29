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
ADD_LIBRARY(sundials_fida_static STATIC ${fida_SOURCES})

# Set the library name and make sure it is not deleted
SET_TARGET_PROPERTIES(sundials_fida_static
  PROPERTIES OUTPUT_NAME sundials_fida CLEAN_DIRECT_OUTPUT 1)

# Install the FIDA library
INSTALL(TARGETS sundials_fida_static DESTINATION lib)

# If Blas/Lapack support was enabled, create and install the wrapper library
# and install the interface headers
IF(LAPACK_OK)
  SET(fida_BL_SOURCES fidalapack.c fidalapdense.c fidalapband.c)
  ADD_LIBRARY(sundials_fida_bl_static STATIC ${fida_BL_SOURCES})
  SET_TARGET_PROPERTIES(sundials_fida_bl_static
    PROPERTIES OUTPUT_NAME sundials_fida_bl CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_fida_bl_static DESTINATION lib)
ENDIF(LAPACK_OK)
