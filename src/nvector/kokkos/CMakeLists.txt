# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos, and Daniel McGreer @ LLNL
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
# CMakeLists.txt file for the Kokkos NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_KOKKOS\n\")")

# Add variable nveckokkos_SOURCES with the sources for the NVECKOKKOS lib
set(nveckokkos_SOURCES
  nvector_kokkos.cpp
)

# Create the library
sundials_add_library(sundials_nveckokkos
  SOURCES
    nvector_kokkos.cpp
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_kokkos.h
  INCLUDE_SUBDIR
    nvector
  OBJECT_LIBRARIES
    sundials_generic_obj
    # ${_sunmemlib}
  LINK_LIBRARIES
    PUBLIC Kokkos::kokkos
  COMPILE_DEFINITIONS
    PRIVATE ${_compile_defs}
  COMPILE_FEATURES
    PUBLIC cxx_std_14
  OUTPUT_NAME
    ${_lib_output_name}
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

# Install the NVECKOKKOS header files
install(FILES ${nveckokkos_HEADERS} DESTINATION include/nvector)

message(STATUS "Added NVECTOR_KOKKOS module")
