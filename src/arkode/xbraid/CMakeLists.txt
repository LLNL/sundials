# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the ARKODE + XBraid interface library
# ---------------------------------------------------------------

# Set install message
install(CODE "MESSAGE(\"\nInstall ARKODE XBraid interface\n\")")

# Use MPI compiler wrapper if available
if(MPI_C_COMPILER)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Source files for the ARKODE + XBraid interface library
set(arkode_SOURCES
  arkode_xbraid.c
  )

# Shared source files for the ARKODE + XBraid interface library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_xbraid.c
  )

# Add variable arkode_HEADERS with the exported ARKODE header files
set(arkode_HEADERS
  arkode_xbraid.h
  )

# Add prefix with complete path to the ARKODE header files
add_prefix(${sundials_SOURCE_DIR}/include/arkode/ arkode_HEADERS)

# Define C preprocessor flag
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library
if(BUILD_STATIC_LIBS)
  add_library(sundials_arkode_xbraid_static STATIC ${arkode_SOURCES} ${shared_SOURCES})
  target_include_directories(sundials_arkode_xbraid_static PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_include_directories(sundials_arkode_xbraid_static PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(sundials_arkode_xbraid_static m sundials_arkode_static XBRAID::XBRAID)
  set_target_properties(sundials_arkode_xbraid_static
    PROPERTIES
    OUTPUT_NAME sundials_arkode_xbraid
    CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_arkode_xbraid_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if(BUILD_SHARED_LIBS)
  add_library(sundials_arkode_xbraid_shared SHARED ${arkode_SOURCES} ${shared_SOURCES})
  target_include_directories(sundials_arkode_xbraid_shared PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_include_directories(sundials_arkode_xbraid_shared PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(sundials_arkode_xbraid_shared m sundials_arkode_shared  XBRAID::XBRAID)
  set_target_properties(sundials_arkode_xbraid_shared
    PROPERTIES
    OUTPUT_NAME sundials_arkode_xbraid
    CLEAN_DIRECT_OUTPUT 1
    VERSION ${arkodelib_VERSION}
    SOVERSION ${arkodelib_SOVERSION})
  install(TARGETS sundials_arkode_xbraid_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# Install the ARKODE header files
install(FILES ${arkode_HEADERS} DESTINATION include/arkode)

# Finished
message(STATUS "Added ARKODE Xbraid interface")
