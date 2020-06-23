# ---------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
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
# CMakeLists.txt file for the CPODES library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall CPODES\n\")")

# Add variable cpodes_SOURCES with the sources for the CPODES library
set(cpodes_SOURCES
  cpodes_band.c
  cpodes_bandpre.c
  cpodes_bbdpre.c
  cpodes.c
  cpodes_dense.c
  cpodes_direct.c
  cpodes_ic.c
  cpodes_io.c
  cpodes_nls.c
  cpodes_proj.c
  cpodes_root.c
  cpodes_spbcgs.c
  cpodes_spgmr.c
  cpodes_spils.c
  cpodes_sptfqmr.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the CPODES library
set(shared_SOURCES
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
add_prefix(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable cpodes_HEADERS with the exported CPODES header files
set(cpodes_HEADERS
  cpodes_band.h
  cpodes_bandpre.h
  cpodes_bbdpre.h
  cpodes_dense.h
  cpodes_direct.h
  cpodes.h
  cpodes_spbcgs.h
  cpodes_spgmr.h
  cpodes_spils.h
  cpodes_sptfqmr.h
  )

# Add prefix with complete path to the CPODES header files
add_prefix(${sundials_SOURCE_DIR}/include/cpodes/ cpodes_HEADERS)

# If Blas/Lapack support was enabled, set-up additional file lists
if(LAPACK_FOUND)
  set(cpodes_BL_SOURCES cpodes_lapack.c)
  set(cpodes_BL_HEADERS cpodes_lapack.h)
  add_prefix(${sundials_SOURCE_DIR}/include/cpodes/ cpodes_BL_HEADERS)
else(LAPACK_FOUND)
  set(cpodes_BL_SOURCES "")
  set(cpodes_BL_HEADERS "")
endif(LAPACK_FOUND)

# Add source directories to include directories for access to
# implementation only header files.
include_directories(.)
include_directories(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
if(BUILD_STATIC_LIBS)

  # Add the build target for the static CPODES library
  add_library(sundials_cpodes_static STATIC
    ${cpodes_SOURCES} ${cpodes_BL_SOURCES} ${shared_SOURCES})

  # Set the library name and make sure it is not deleted
  set_target_properties(sundials_cpodes_static
    PROPERTIES OUTPUT_NAME sundials_cpodes CLEAN_DIRECT_OUTPUT 1)

  # Install the CPODES library
  install(TARGETS sundials_cpodes_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_STATIC_LIBS)

# Build the shared library
if(BUILD_SHARED_LIBS)

  # Add the build target for the CPODES library
  add_library(sundials_cpodes_shared SHARED
    ${cpodes_SOURCES} ${cpodes_BL_SOURCES} ${shared_SOURCES})

  # Set the library name and make sure it is not deleted
  set_target_properties(sundials_cpodes_shared
    PROPERTIES OUTPUT_NAME sundials_cpodes CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  set_target_properties(sundials_cpodes_shared
    PROPERTIES VERSION ${cpodeslib_VERSION} SOVERSION ${cpodeslib_SOVERSION})

  # Install the CPODES library
  install(TARGETS sundials_cpodes_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_SHARED_LIBS)

# Install the CPODES header files
install(FILES ${cpodes_HEADERS} ${cpodes_BL_HEADERS} DESTINATION include/cpodes)

# Install the CPODES implementation header file
install(FILES cpodes_impl.h DESTINATION include/cpodes)

#
message(STATUS "Added CPODES module")
