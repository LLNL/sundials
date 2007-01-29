# ---------------------------------------------------------------
# $Revision: 1.1 $
# $Date: 2007-01-29 17:40:00 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the serial NVECTOR library


# Add variable nvecserial_SOURCES with the sources for the NVECSERIAL lib
SET(nvecserial_SOURCES nvector_serial.c)

# Add variable nvecserial_HEADERS with the exported NVECSERIAL header files
SET(nvecserial_HEADERS nvector_serial.h)

# Add prefix with complete path to the NVECSERIAL header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nvecserial_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
IF(BUILD_STATIC_LIBS)

  # Add the build target for the NVECSERIAL library
  ADD_LIBRARY(sundials_nvecserial_static STATIC ${nvecserial_SOURCES})
  
  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_nvecserial_static
    PROPERTIES OUTPUT_NAME sundials_nvecserial CLEAN_DIRECT_OUTPUT 1)

  # Install the NVECSERIAL library
  INSTALL(TARGETS sundials_nvecserial_static DESTINATION lib)

ENDIF(BUILD_STATIC_LIBS)

# Build the shared library
IF(BUILD_SHARED_LIBS)

  # Add the build target for the NVECSERIAL library
  ADD_LIBRARY(sundials_nvecserial_shared SHARED ${nvecserial_SOURCES})
  
  # Set the library name and make sure it is not deleted
  SET_TARGET_PROPERTIES(sundials_nvecserial_shared
    PROPERTIES OUTPUT_NAME sundials_nvecserial CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  SET_TARGET_PROPERTIES(sundials_nvecserial_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  
  # Install the NVECSERIAL library
  INSTALL(TARGETS sundials_nvecserial_shared DESTINATION lib)

ENDIF(BUILD_SHARED_LIBS)


# Install the NVECSERIAL header files
INSTALL(FILES ${nvecserial_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled, build and install the FNVECSERIAL library
IF(FCMIX_OK)
  SET(fnvecserial_SOURCES fnvector_serial.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fnvecserial_static STATIC ${fnvecserial_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecserial_static
      PROPERTIES OUTPUT_NAME sundials_fnvecserial CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fnvecserial_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fnvecserial_shared ${fnvecserial_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecserial_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecserial CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fnvecserial_shared 
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    INSTALL(TARGETS sundials_fnvecserial_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(FCMIX_OK)
