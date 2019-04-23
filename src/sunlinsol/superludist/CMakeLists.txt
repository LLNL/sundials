# ------------------------------------------------------------------------------
# Programmer:  Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# CMakeLists.txt file for the superlu-dist SUNLinearSolver library
# -----------------------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_SUPERLUDIST\n\")")

# Add variable sunlinsolsuperludist_SOURCES with the sources for the
# SUNLINSOLSUPERLUDIST lib
set(sunlinsolsuperludist_SOURCES sunlinsol_superludist.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLSUPERLUDIST library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  )

# Add variable sunlinsolsuperludist_HEADERS with the exported
# SUNLINSOLSUPERLUDIST header files
set(sunlinsolsuperludist_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_superludist.h
  )

# Use MPI
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

if(BUILD_STATIC_LIBS)
  add_library(sundials_sunlinsolsuperludist_static STATIC
    ${sunlinsolsuperludist_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsolsuperludist_static PROPERTIES
    OUTPUT_NAME sundials_sunlinsolsuperludist CLEAN_DIRECT_OUTPUT 1
  )

  target_link_libraries(sundials_sunlinsolsuperludist_static
    PUBLIC
      SuperLU_DIST::SuperLU_DIST
      sundials_sunmatrixslunrloc_static
    PRIVATE
      m
  )

  if(SUPERLUDIST_OpenMP)
    target_compile_options(sundials_sunlinsolsuperludist_static
      PUBLIC ${OpenMP_C_FLAGS})
  endif()

  target_compile_definitions(sundials_sunlinsolsuperludist_static
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsolsuperludist_static
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

if(BUILD_SHARED_LIBS)
  add_library(sundials_sunlinsolsuperludist_shared SHARED
    ${sunlinsolsuperludist_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsolsuperludist_shared PROPERTIES
    OUTPUT_NAME sundials_sunlinsolsuperludist CLEAN_DIRECT_OUTPUT 1
    VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION}
  )

  target_link_libraries(sundials_sunlinsolsuperludist_shared
    PUBLIC
      SuperLU_DIST::SuperLU_DIST
      sundials_sunmatrixslunrloc_shared
    PRIVATE
      m
  )

  if(SUPERLUDIST_OpenMP AND OPENMP_ENABLE AND OPENMP_FOUND)
    target_compile_options(sundials_sunlinsolsuperludist_shared
      PUBLIC ${OpenMP_C_FLAGS})
  endif()

  target_compile_definitions(sundials_sunlinsolsuperludist_shared
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsolsuperludist_shared
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the SUNLINSOLSUPERLUDIST header files
install(FILES ${sunlinsolsuperludist_HEADERS} DESTINATION include/sunlinsol)

message(STATUS "Added SUNLINSOL_SUPERLUDIST module")
